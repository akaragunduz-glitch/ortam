# ORTAM — Architecture

Version: 0.1 · Date: 2026-11-07

This is the in-repository summary that Claude refers to. The full
architecture, technology stack and roadmap are in the project plan document.

## Operating modes

| Mode | Example | Required components |
|------|---------|---------------------|
| Standalone module (CLI) | `ortam run groundwater.three-point input.yaml` | Python + the module's engine only |
| Standalone module (UI) | Groundwater model screen only | UI + API (+ job queue for heavy models) |
| Integrated workflow | Sampling → modelling → risk → remediation | All layers |

## Layers

```
┌──────────────────────────────────────────────────────────────┐
│ User interface (React + TypeScript)                           │
│  2D map (MapLibre) · 3D (CesiumJS) · charts/boring logs (D3)  │
└───────────────────────────┬──────────────────────────────────┘
                            │ HTTPS / JSON
┌───────────────────────────▼──────────────────────────────────┐
│ API server (FastAPI) — validation, auth, fast computations    │
└──────┬────────────────────┬─────────────────────┬────────────┘
       │                    │                     │
┌──────▼──────┐   ┌─────────▼────────┐   ┌────────▼───────┐
│ PostgreSQL  │   │ Job queue        │   │ File store     │
│ + PostGIS   │   │ (Celery + Redis) │   │ (MinIO)        │
└──────▲──────┘   └─────────┬────────┘   └────────▲───────┘
       │ results            │ jobs                │ files
┌──────┴────────────────────▼─────────────────────┴───────────┐
│ Compute engines (Docker containers)                           │
│  Python modules: RAGS, J&E, HHRAP, ProUCL, kriging, vadose    │
│  Fortran: MODFLOW 6, MT3D-USGS, AERMOD, VLEACH, VS2DT, MINTEQA2│
│  Libraries: IPhreeqc, GemPy, PyVista                          │
└───────────────────────────────────────────────────────────────┘
```

## Design principles

0. **Every module runs standalone:** every module implements
   `ortam_core.model.Module` (`validate` → `run` → `summary`). Input can be a
   YAML/JSON file or a form; the database is optional. See
   `docs/interfaces/model-contract.md`.
1. **Science first, UI second:** computations in `modules/` are independent
   of the UI and tested on their own.
2. **Contract first:** each module's input/output schema is defined before
   coding.
3. **Modules never call each other; workflows connect them:** passing one
   module's output to another is the job of `workflows/` (e.g.
   `workflows/risk_chain/`: groundwater → risk → remediation).
4. **Fortran codes are called through file-based wrappers;** small numerical
   kernels (Richards solver) are bound with f2py.
5. **Reproducibility:** inputs, code version and parameters of every model
   run are stored.

## Technology stack (summary)

| Layer | Tools |
|-------|-------|
| Languages | Python 3.12, TypeScript |
| Scientific | NumPy, SciPy, pandas, statsmodels, PyKrige, GSTools, pint |
| Fortran bridge | subprocess, f2py, FloPy |
| API | FastAPI, Pydantic, SQLAlchemy, GeoAlchemy2, Alembic |
| Database | PostgreSQL + PostGIS |
| Job queue | Celery + Redis |
| UI | React, MapLibre GL JS, CesiumJS, Plotly, D3.js, i18n (react-i18next) |
| Packaging | Docker, Docker Compose |
| Testing | pytest, GitHub Actions |
