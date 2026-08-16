# Multiomic Data Orchestrator (MDO)

Deterministic metadata harmonization for FFPE-based multiomic studies.

Metadata in multiomic labs arrives fragmented across platforms — spatial transcriptomics, sequencing assays, library preps — each with its own spreadsheet conventions. A mismatch that slips through costs a failed instrument run in reagents and time. MDO enforces a canonical data model, validates referential integrity across it, and refuses to export until the blockers are cleared.

## The problem it solves

A study's entities form a chain: **Block → Slide → ROI → Library → Run**. MDO validates that every link in that chain resolves — that no library points at a missing ROI, no run at a missing library — and that identifiers and formats are normalised the same way every time. The result is an auditable Chain of Identity suitable for downstream pipelines and regulatory submission.

## The workflow

1. **Ingest** — upload CSVs and map columns onto versioned schema templates (Illumina, 10x Genomics, Spatial)
2. **Harmonize** — canonical identifiers enforced, formats normalised
3. **Validate** — a versioned rules engine checks referential integrity and reports every violation
4. **Remediate** — fix errors offline in the source files, re-upload, re-validate
5. **Export** — once no Blockers remain, emit a bundle: canonical tables, join index, and manifest

Export is *gated*. Blocker-class errors cannot be waived — that is the point of the tool.

## Who it's for

| Persona | Use |
|---|---|
| Platform Operators | QC metadata before an instrument run |
| Bioinformaticians | Consume harmonized metadata to drive analysis pipelines |
| Principal Investigators | Confidence in multiomic data integrity |
| QA / Regulatory | Auditable, reproducible data lineage |

## Stack

| Layer | Technology |
|---|---|
| Backend | FastAPI (Python) |
| Frontend | React, TypeScript, Vite, Tailwind CSS, shadcn/ui |
| Design | Deterministic and explicit rules engine — no heuristics, no external API dependencies |

## Running it locally

**Backend**

```bash
cd backend
cp env.example .env
pip install -r requirements.txt
python start.py
```

**Frontend**

```bash
cd frontend
npm install
npm run dev
```

## Documentation

The full product specification — personas, rulesets, entity model, MVP success criteria — is in [`frontend/PRD.md`](frontend/PRD.md). Backend design notes are in [`Backend-dev-plan.md`](Backend-dev-plan.md).
