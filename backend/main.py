from fastapi import FastAPI, Depends, HTTPException, status, Body, File, UploadFile, Header
from typing import List, Optional, Dict, Any, Tuple
import json
import io
import re
from bson import ObjectId
import asyncio
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, EmailStr
from passlib.context import CryptContext
from jose import jwt, JWTError
from datetime import datetime, timedelta
from motor.motor_asyncio import AsyncIOMotorClient
import os
from dotenv import load_dotenv
import pandas as pd

load_dotenv()

app = FastAPI(title="Multiomic Data Orchestrator API", version="1.0.0")

# CORS Configuration
# Get CORS origins from environment variable
# Default to localhost for local development (allows credentials)
cors_origins_str = os.getenv("CORS_ORIGINS", "http://localhost:5173,http://localhost:3000")
# Parse and clean origins - remove trailing slashes and whitespace
origins = []
for origin in cors_origins_str.split(","):
    cleaned = origin.strip().rstrip("/")
    if cleaned:
        origins.append(cleaned)

# Log CORS configuration (will be printed on startup)
print(f"[STARTUP] CORS_ORIGINS env: {cors_origins_str}")
print(f"[STARTUP] CORS allowed origins: {origins}")

# CORS middleware configuration
# Note: When allow_credentials=True, you CANNOT use allow_origins=["*"]
# This is a CORS specification requirement. We always use specific origins.
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["GET", "POST", "PUT", "DELETE", "OPTIONS", "PATCH", "HEAD"],
    allow_headers=["*"],
    expose_headers=["*"],
)

# Configuration
MONGODB_URI = os.getenv("MONGODB_URI")
if not MONGODB_URI:
    raise ValueError("MONGODB_URI environment variable is required")

SECRET_KEY = os.getenv("JWT_SECRET")
if not SECRET_KEY:
    raise ValueError("JWT_SECRET environment variable is required")

ALGORITHM = "HS256"
ACCESS_TOKEN_EXPIRE_MINUTES = int(os.getenv("JWT_EXPIRES_IN", "60"))

# Database
DATABASE_NAME = os.getenv("DATABASE_NAME", "mdo")
client = AsyncIOMotorClient(MONGODB_URI)
db = client.get_database(DATABASE_NAME)

# Password hashing - Using Argon2 (no password length limit, more secure than bcrypt)
pwd_context = CryptContext(schemes=["argon2"], deprecated="auto")

# ============== Schema Templates with Validation Rules ==============

SCHEMA_TEMPLATES = {
    "illumina-ngs-run-v1.2": {
        "name": "Illumina NGS Run",
        "fields": [
            {"name": "Run_ID", "type": "string", "required": True},
            {"name": "Sample_ID", "type": "string", "required": True},
            {"name": "Library_ID", "type": "string", "required": True},
            {"name": "Lane", "type": "integer", "required": True, "min": 1, "max": 8},
            {"name": "Index_I1", "type": "string", "required": True, "pattern": r"^[ATGCN]+$"},
            {"name": "Index_I2", "type": "string", "required": False, "pattern": r"^[ATGCN]*$"},
        ]
    },
    "10x-single-cell-v2.0": {
        "name": "10x Single-Cell",
        "fields": [
            {"name": "Library_ID", "type": "string", "required": True},
            {"name": "Sample_ID", "type": "string", "required": True},
            {"name": "Chemistry", "type": "string", "required": True},
            {"name": "Expected_Cells", "type": "integer", "required": False, "min": 1},
        ]
    },
    "spatial-visium-v1.0": {
        "name": "Spatial - Visium",
        "fields": [
            {"name": "Slide_ID", "type": "string", "required": True},
            {"name": "Capture_Area", "type": "string", "required": True, "pattern": r"^[A-D][1-4]$"},
            {"name": "Library_ID", "type": "string", "required": True},
            {"name": "Block_ID", "type": "string", "required": True},
        ]
    },
}

# ============== Validation Engine ==============

def validate_csv_data(file_content: bytes, file_name: str, template_id: str, column_mapping: Dict[str, str]) -> List[Dict[str, Any]]:
    """
    Validate CSV data against schema template rules.
    Returns a list of validation issues.
    """
    issues = []
    issue_id = 0
    
    template = SCHEMA_TEMPLATES.get(template_id)
    if not template:
        issues.append({
            "id": str(issue_id := issue_id + 1),
            "severity": "Blocker",
            "fileName": file_name,
            "rowIndex": None,
            "columnName": None,
            "description": f"Unknown template: {template_id}"
        })
        return issues
    
    # Parse CSV
    try:
        df = pd.read_csv(io.BytesIO(file_content))
    except Exception as e:
        issues.append({
            "id": str(issue_id := issue_id + 1),
            "severity": "Blocker",
            "fileName": file_name,
            "rowIndex": None,
            "columnName": None,
            "description": f"Failed to parse CSV: {str(e)}"
        })
        return issues
    
    # Track values for uniqueness checks
    seen_values = {}
    
    for field in template["fields"]:
        canonical_name = field["name"]
        csv_column = column_mapping.get(canonical_name)
        
        # Check if required field is mapped
        if not csv_column:
            if field["required"]:
                issues.append({
                    "id": str(issue_id := issue_id + 1),
                    "severity": "Blocker",
                    "fileName": file_name,
                    "rowIndex": None,
                    "columnName": canonical_name,
                    "description": f"Required field '{canonical_name}' is not mapped to any CSV column"
                })
            continue
        
        # Check if mapped column exists in CSV
        if csv_column not in df.columns:
            issues.append({
                "id": str(issue_id := issue_id + 1),
                "severity": "Blocker",
                "fileName": file_name,
                "rowIndex": None,
                "columnName": canonical_name,
                "description": f"Mapped column '{csv_column}' not found in CSV file"
            })
            continue
        
        # Initialize uniqueness tracking for this field
        seen_values[canonical_name] = {}
        
        # Validate each row
        for idx, value in df[csv_column].items():
            row_num = idx + 2  # +2 for header row and 0-indexing
            
            # Required field check
            if field["required"] and (pd.isna(value) or str(value).strip() == ""):
                issues.append({
                    "id": str(issue_id := issue_id + 1),
                    "severity": "Blocker",
                    "fileName": file_name,
                    "rowIndex": row_num,
                    "columnName": canonical_name,
                    "description": f"Required field '{canonical_name}' is empty"
                })
                continue
            
            # Skip further validation if value is empty (for non-required fields)
            if pd.isna(value) or str(value).strip() == "":
                continue
            
            str_value = str(value).strip()
            
            # Type validation
            if field["type"] == "integer":
                try:
                    int_val = int(float(value))
                    # Range validation
                    if "min" in field and int_val < field["min"]:
                        issues.append({
                            "id": str(issue_id := issue_id + 1),
                            "severity": "Blocker",
                            "fileName": file_name,
                            "rowIndex": row_num,
                            "columnName": canonical_name,
                            "description": f"Value {int_val} is below minimum {field['min']}"
                        })
                    if "max" in field and int_val > field["max"]:
                        issues.append({
                            "id": str(issue_id := issue_id + 1),
                            "severity": "Blocker",
                            "fileName": file_name,
                            "rowIndex": row_num,
                            "columnName": canonical_name,
                            "description": f"Value {int_val} exceeds maximum {field['max']}"
                        })
                except (ValueError, TypeError):
                    issues.append({
                        "id": str(issue_id := issue_id + 1),
                        "severity": "Blocker",
                        "fileName": file_name,
                        "rowIndex": row_num,
                        "columnName": canonical_name,
                        "description": f"Expected integer, got '{value}'"
                    })
            
            # Pattern validation (e.g., DNA sequences, capture area format)
            if "pattern" in field:
                if not re.match(field["pattern"], str_value):
                    issues.append({
                        "id": str(issue_id := issue_id + 1),
                        "severity": "Warning",
                        "fileName": file_name,
                        "rowIndex": row_num,
                        "columnName": canonical_name,
                        "description": f"Value '{str_value}' doesn't match expected format"
                    })
            
            # Uniqueness check for ID fields
            if canonical_name.endswith("_ID"):
                if str_value in seen_values[canonical_name]:
                    first_row = seen_values[canonical_name][str_value]
                    issues.append({
                        "id": str(issue_id := issue_id + 1),
                        "severity": "Warning",
                        "fileName": file_name,
                        "rowIndex": row_num,
                        "columnName": canonical_name,
                        "description": f"Duplicate value '{str_value}' (first seen in row {first_row})"
                    })
                else:
                    seen_values[canonical_name][str_value] = row_num
    
    # If no issues found, add an info message
    if len(issues) == 0:
        issues.append({
            "id": str(issue_id := issue_id + 1),
            "severity": "Info",
            "fileName": file_name,
            "rowIndex": None,
            "columnName": None,
            "description": "All validations passed successfully!"
        })
    
    return issues

# ============== Pydantic Models ==============

class UserInDB(BaseModel):
    email: EmailStr
    name: str
    password_hash: str

class UserResponse(BaseModel):
    email: EmailStr
    name: str

class Token(BaseModel):
    access_token: str
    token_type: str

class LoginRequest(BaseModel):
    email: EmailStr
    password: str

class SignUpRequest(BaseModel):
    email: EmailStr
    name: str
    password: str

class MappingCreate(BaseModel):
    name: str
    mappings: list

class MappingResponse(BaseModel):
    id: str
    name: str
    mappings: list
    createdAt: str

# ============== Helper Functions ==============

def verify_password(plain_password: str, hashed_password: str) -> bool:
    return pwd_context.verify(plain_password, hashed_password)

def hash_password(password: str) -> str:
    """
    Hash a password using Argon2.
    Argon2 has no password length limit and is more secure than bcrypt.
    """
    return pwd_context.hash(password)

def validate_password(password: str) -> Tuple[bool, Optional[str]]:
    """
    Validate password strength.
    Returns (is_valid, error_message)
    Note: Argon2 has no password length limit, but we enforce reasonable limits for security.
    """
    if len(password) < 8:
        return False, "Password must be at least 8 characters long"
    
    # Reasonable upper limit (not a technical limitation, just good practice)
    if len(password) > 128:
        return False, "Password is too long (maximum 128 characters). Please use a shorter password."
    
    # Check for at least one letter and one number
    has_letter = any(c.isalpha() for c in password)
    has_number = any(c.isdigit() for c in password)
    if not (has_letter and has_number):
        return False, "Password must contain at least one letter and one number"
    return True, None

def create_access_token(data: dict, expires_delta: Optional[timedelta] = None) -> str:
    to_encode = data.copy()
    if expires_delta:
        expire = datetime.utcnow() + expires_delta
    else:
        expire = datetime.utcnow() + timedelta(minutes=15)
    to_encode.update({"exp": expire})
    encoded_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)
    return encoded_jwt

async def get_user_by_email(email: str) -> Optional[dict]:
    user = await db.users.find_one({"email": email})
    return user

# ============== JWT Authentication Dependency ==============

async def get_current_user(authorization: Optional[str] = Header(None)) -> dict:
    """
    Dependency to extract and validate JWT token from Authorization header.
    Returns the user document from the database.
    """
    credentials_exception = HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Could not validate credentials",
        headers={"WWW-Authenticate": "Bearer"},
    )
    
    if authorization is None:
        raise credentials_exception
    
    # Extract token from "Bearer <token>"
    parts = authorization.split()
    if len(parts) != 2 or parts[0].lower() != "bearer":
        raise credentials_exception
    
    token = parts[1]
    
    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
        email: str = payload.get("sub")
        if email is None:
            raise credentials_exception
    except JWTError:
        raise credentials_exception
    
    user = await get_user_by_email(email)
    if user is None:
        raise credentials_exception
    
    return user

# ============== Startup Event - Seed User ==============

@app.on_event("startup")
async def startup_event():
    """Create seed user if configured and doesn't exist"""
    seed_email = os.getenv("SEED_USER_EMAIL")
    seed_password = os.getenv("SEED_USER_PASSWORD")
    seed_name = os.getenv("SEED_USER_NAME", "Admin")
    
    # Only create seed user if email and password are provided
    if seed_email and seed_password:
        existing_user = await db.users.find_one({"email": seed_email})
        if not existing_user:
            user_doc = {
                "email": seed_email,
                "name": seed_name,
                "password_hash": hash_password(seed_password),
                "created_at": datetime.utcnow()
            }
            await db.users.insert_one(user_doc)
            print(f"[STARTUP] Created seed user: {seed_email}")
        else:
            print(f"[STARTUP] Seed user already exists: {seed_email}")
    else:
        print(f"[STARTUP] Seed user not configured (SEED_USER_EMAIL and SEED_USER_PASSWORD not set)")

# ============== Health Check ==============

@app.get("/healthz")
async def health_check():
    """
    Health check endpoint that verifies database connectivity.
    Returns status and database connection status.
    """
    try:
        # Ping MongoDB to verify connection
        await client.admin.command('ping')
        db_status = "connected"
    except Exception as e:
        db_status = f"disconnected: {str(e)}"
    
    return {
        "status": "ok" if db_status == "connected" else "degraded",
        "db_status": db_status,
        "timestamp": datetime.utcnow().isoformat()
    }

# ============== Authentication Endpoints ==============

@app.post("/api/v1/auth/signup", response_model=Token)
async def signup(request: SignUpRequest = Body(...)):
    """
    Register a new user account and return JWT token.
    """
    # Validate password strength
    is_valid, error_msg = validate_password(request.password)
    if not is_valid:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=error_msg,
        )
    
    # Check if user already exists
    existing_user = await get_user_by_email(request.email)
    if existing_user:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Email already registered",
        )
    
    # Create new user
    password_hash = hash_password(request.password)
    
    user_doc = {
        "email": request.email,
        "name": request.name,
        "password_hash": password_hash,
        "created_at": datetime.utcnow()
    }
    
    result = await db.users.insert_one(user_doc)
    
    # Log the signup event
    await db.audit_logs.insert_one({
        "action": "USER_SIGNUP",
        "user_email": request.email,
        "timestamp": datetime.utcnow(),
        "details": {}
    })
    
    # Auto-login: generate and return JWT token
    access_token_expires = timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES)
    access_token = create_access_token(
        data={"sub": request.email}, expires_delta=access_token_expires
    )
    
    return {"access_token": access_token, "token_type": "bearer"}

@app.post("/api/v1/auth/login", response_model=Token)
async def login(request: LoginRequest = Body(...)):
    """
    Authenticate user and return JWT token.
    """
    user = await get_user_by_email(request.email)
    
    if not user or not verify_password(request.password, user["password_hash"]):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect email or password",
        )
    
    access_token_expires = timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES)
    access_token = create_access_token(
        data={"sub": user["email"]}, expires_delta=access_token_expires
    )
    
    # Log the login event
    await db.audit_logs.insert_one({
        "action": "USER_LOGIN",
        "user_email": user["email"],
        "timestamp": datetime.utcnow(),
        "details": {}
    })
    
    return {"access_token": access_token, "token_type": "bearer"}

@app.get("/api/v1/auth/me", response_model=UserResponse)
async def get_current_user_info(current_user: dict = Depends(get_current_user)):
    """
    Get current logged-in user's information.
    Requires valid JWT token.
    """
    return {
        "email": current_user["email"],
        "name": current_user["name"]
    }

@app.post("/api/v1/auth/logout")
async def logout(current_user: dict = Depends(get_current_user)):
    """
    Logout the current user.
    Logs the logout event for audit purposes.
    Note: Since JWTs are stateless, we can't invalidate the token server-side,
    but we log the logout event for audit trails.
    """
    # Log the logout event
    await db.audit_logs.insert_one({
        "action": "USER_LOGOUT",
        "user_email": current_user["email"],
        "timestamp": datetime.utcnow(),
        "details": {}
    })
    
    return {"status": "success", "message": "Logged out successfully"}

# ============== Mapping Configuration Endpoints ==============

@app.get("/api/v1/mappings")
async def get_mappings(current_user: dict = Depends(get_current_user)):
    """
    Get all saved mapping configurations for the current user.
    """
    user_id = current_user["_id"]
    mappings = await db.mapping_configurations.find({"user_id": user_id}).to_list(100)
    
    result = []
    for m in mappings:
        result.append({
            "id": str(m["_id"]),
            "name": m["name"],
            "mappings": m["mappings"],
            "createdAt": m["created_at"].isoformat() if isinstance(m["created_at"], datetime) else m["created_at"]
        })
    
    return result

@app.post("/api/v1/mappings")
async def save_mapping(mapping: MappingCreate, current_user: dict = Depends(get_current_user)):
    """
    Save a new mapping configuration for the current user.
    """
    user_id = current_user["_id"]
    
    mapping_doc = {
        "user_id": user_id,
        "name": mapping.name,
        "mappings": mapping.mappings,
        "created_at": datetime.utcnow()
    }
    
    result = await db.mapping_configurations.insert_one(mapping_doc)
    
    # Log the action
    await db.audit_logs.insert_one({
        "action": "MAPPING_SAVED",
        "user_email": current_user["email"],
        "timestamp": datetime.utcnow(),
        "details": {"mapping_id": str(result.inserted_id), "mapping_name": mapping.name}
    })
    
    return {
        "id": str(result.inserted_id),
        "name": mapping.name,
        "mappings": mapping.mappings,
        "createdAt": mapping_doc["created_at"].isoformat()
    }

@app.delete("/api/v1/mappings/{mapping_id}")
async def delete_mapping(mapping_id: str, current_user: dict = Depends(get_current_user)):
    """
    Delete a saved mapping configuration.
    Only the owner can delete their mappings.
    """
    user_id = current_user["_id"]
    
    # Verify ownership
    mapping = await db.mapping_configurations.find_one({
        "_id": ObjectId(mapping_id),
        "user_id": user_id
    })
    
    if not mapping:
        raise HTTPException(status_code=404, detail="Mapping not found")
    
    await db.mapping_configurations.delete_one({"_id": ObjectId(mapping_id)})
    
    # Log the action
    await db.audit_logs.insert_one({
        "action": "MAPPING_DELETED",
        "user_email": current_user["email"],
        "timestamp": datetime.utcnow(),
        "details": {"mapping_id": mapping_id}
    })
    
    return {"status": "success"}

# ============== Harmonization Run Endpoints ==============

@app.post("/api/v1/runs")
async def start_run(
    files: List[UploadFile] = File(...),
    mapping: str = Body(...),
    current_user: dict = Depends(get_current_user)
):
    """
    Initiate a new harmonization and validation run.
    Accepts uploaded CSV files and mapping configuration.
    Now performs REAL validation against schema templates.
    """
    mapping_data = json.loads(mapping)
    user_id = current_user["_id"]
    
    # Read file contents before they're consumed
    file_contents = {}
    for file in files:
        content = await file.read()
        file_contents[file.filename] = content
        await file.seek(0)  # Reset file pointer
    
    run_doc = {
        "user_id": user_id,
        "status": "pending",
        "mapping": mapping_data,
        "files": [file.filename for file in files],
        "created_at": datetime.utcnow(),
        "validation_issues": []
    }
    
    result = await db.harmonization_runs.insert_one(run_doc)
    run_id = result.inserted_id
    
    # Log the action
    await db.audit_logs.insert_one({
        "action": "RUN_STARTED",
        "user_email": current_user["email"],
        "timestamp": datetime.utcnow(),
        "details": {"run_id": str(run_id), "files": run_doc["files"]}
    })

    async def complete_run():
        """Background task to perform REAL validation"""
        await asyncio.sleep(1)  # Small delay for UI feedback
        
        all_issues = []
        
        # Process each file with its mapping
        for file_mapping in mapping_data:
            file_name = file_mapping.get("fileName")
            template_id = file_mapping.get("templateId")
            column_mapping = file_mapping.get("mapping", {})
            
            if file_name and file_name in file_contents:
                # Run real validation
                issues = validate_csv_data(
                    file_content=file_contents[file_name],
                    file_name=file_name,
                    template_id=template_id,
                    column_mapping=column_mapping
                )
                all_issues.extend(issues)
        
        # If no files were processed, add an info message
        if len(all_issues) == 0:
            all_issues.append({
                "id": "1",
                "severity": "Info",
                "fileName": "N/A",
                "rowIndex": None,
                "columnName": None,
                "description": "No files were processed for validation"
            })
        
        await db.harmonization_runs.update_one(
            {"_id": run_id},
            {"$set": {
                "status": "complete",
                "validation_issues": all_issues,
                "completed_at": datetime.utcnow()
            }}
        )
        
        # Log completion
        await db.audit_logs.insert_one({
            "action": "RUN_COMPLETED",
            "user_email": current_user["email"],
            "timestamp": datetime.utcnow(),
            "details": {
                "run_id": str(run_id),
                "blockers": len([i for i in all_issues if i["severity"] == "Blocker"]),
                "warnings": len([i for i in all_issues if i["severity"] == "Warning"]),
                "infos": len([i for i in all_issues if i["severity"] == "Info"])
            }
        })

    asyncio.create_task(complete_run())
    
    return {"status": "run started", "run_id": str(run_id)}

@app.get("/api/v1/runs")
async def list_runs(
    current_user: dict = Depends(get_current_user),
    limit: int = 20
):
    """
    Get all harmonization runs for the current user.
    Returns a list of runs sorted by most recent first.
    """
    user_id = current_user["_id"]
    
    # Query runs for this user, sorted by most recent first
    cursor = db.harmonization_runs.find({"user_id": user_id}).sort("created_at", -1).limit(limit)
    runs = await cursor.to_list(length=limit)
    
    result = []
    for run in runs:
        # Calculate validation summary
        validation_issues = run.get("validation_issues", [])
        blockers = len([i for i in validation_issues if i.get("severity") == "Blocker"])
        warnings = len([i for i in validation_issues if i.get("severity") == "Warning"])
        infos = len([i for i in validation_issues if i.get("severity") == "Info"])
        
        run_data = {
            "id": str(run["_id"]),
            "status": run.get("status", "pending"),
            "files": run.get("files", []),
            "created_at": run["created_at"].isoformat() if isinstance(run.get("created_at"), datetime) else run.get("created_at"),
            "completed_at": run["completed_at"].isoformat() if isinstance(run.get("completed_at"), datetime) else run.get("completed_at"),
            "validation_summary": {
                "blockers": blockers,
                "warnings": warnings,
                "infos": infos,
                "total": len(validation_issues)
            }
        }
        result.append(run_data)
    
    return result

@app.get("/api/v1/runs/{run_id}")
async def get_run(run_id: str, current_user: dict = Depends(get_current_user)):
    """
    Get the status and results of a harmonization run.
    Only the owner can view their runs.
    """
    user_id = current_user["_id"]
    
    run = await db.harmonization_runs.find_one({
        "_id": ObjectId(run_id),
        "user_id": user_id
    })
    
    if not run:
        raise HTTPException(status_code=404, detail="Run not found")
    
    run["_id"] = str(run["_id"])
    run["user_id"] = str(run["user_id"])
    
    # Convert datetime objects to ISO strings
    if "created_at" in run and isinstance(run["created_at"], datetime):
        run["created_at"] = run["created_at"].isoformat()
    if "completed_at" in run and isinstance(run["completed_at"], datetime):
        run["completed_at"] = run["completed_at"].isoformat()
    
    return run

# ============== Run the application ==============

if __name__ == "__main__":
    import uvicorn
    port = int(os.getenv("PORT", "8000"))
    uvicorn.run(app, host="0.0.0.0", port=port)
