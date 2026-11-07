# CLAUDE.md — Project working rules

This is the first file every Claude session reads in this repository.
Keep it short, current and binding.

## Project

**ORTAM** (Open Risk, Transport And Modelling): a modular environmental
modelling platform with modules for groundwater, vadose zone, vapor
intrusion, air dispersion, surface runoff, food chain, geochemistry,
statistics, site characterization, sampling design, risk assessment and
remediation.

**Core principle: every module runs standalone.** A user can work only with
the groundwater model; risk assessment is just one of the workflows that
connect modules (`workflows/risk_chain/`).

- Scope: `docs/00-scope.md`
- Architecture: `docs/01-architecture.md`
- Module contract: `docs/interfaces/model-contract.md`
- Daily journal: `docs/journal/YYYY-MM-DD.md`
- Reference documents (EPA, USGS PDFs): `docs/references/`

## Language

- **The software is in English:** code, identifiers, comments, docstrings,
  CLI messages, UI text, repository documentation and commit messages.
- **Talk to the project owner in Turkish**, in academic language; write
  units as symbols (m3, m3/h, %, mg/kg, µg/L, µg/m3).
- Design the UI for internationalization (i18n) from the start; English is
  the default locale.

## About the project owner

The project owner is experienced in scientific numerical methods and has
written Fortran in the past; software development experience is limited.
- Explain code using Fortran analogies.
- Explain in one sentence what every terminal command does.

## Binding rules

1. **No equations from memory.** Every equation comes from a document in
   `docs/references/`; the code cites source, page and equation number in a
   comment, e.g. `# RAGS Part B (1991), p. 23, Eq. 4`.
2. **Tests first.** Before a new computation, write a `pytest` test based on
   reference values (EPA RSL/VISL calculator output, guidance example,
   analytical solution).
3. **Explicit units.** Inputs and outputs document their units; use `pint`
   where practical. Variable names may carry units (`c_soil_mg_kg`).
4. **Non-detects** are never stored as zero or DL/2; they keep a separate
   qualifier field.
5. **Modules never call each other.** Every module implements
   `ortam_core.model.Module` (`validate` → `run` → `summary`) and must run
   from an input file alone, without the database. Connections between
   modules live only in `workflows/`.
6. **Three ways to run every module:** CLI (`ortam run <module> <input.yaml>`),
   Python (`module.run(...)`) and the web UI. Every module has at least one
   input file in `examples/` and reference tests in `tests/`.
7. **Fortran codes** live in `engines/<model>/` and are compiled inside
   Docker with a pinned compiler version; any change to the source is
   documented as a diff.
8. **Plan first.** For any task that touches more than one file, present a
   plan and wait for approval.
9. **End of every session:** write a journal entry in `docs/journal/`
   (done, remaining, open questions).

## Folder structure

```
ortam_core/        Shared core: module contract, units, CLI
modules/           Standalone modules (each runs on its own)
  groundwater/     Three-point, analytical solutions, MODFLOW 6 flow + GWT, MODPATH
  vadose/          SSL leaching, VLEACH, Richards solver, UZF
  vapor/           Johnson & Ettinger vapor intrusion
  air/             AERMOD, emission factors
  surface/         RUSLE, surface runoff
  foodchain/       HHRAP plant and food chain uptake
  geochem/         PHREEQC, MINTEQA2
  stats/           ProUCL methods
  characterization/ Boring logs, potentiometric surface, kriging
  sampling/        Sampling design (VSP methods)
  risk/            RAGS exposure and risk
  remediation/     Risk-based cleanup goals and strategy
  registry.py      Module registry (read by the CLI and the UI)
  <module>/examples/  Example input files (YAML)
  <module>/tests/     Reference tests
workflows/         Workflows that connect modules
  risk_chain/      Sampling -> modelling -> risk -> remediation
engines/           Fortran and external codes (MODFLOW, AERMOD, VLEACH, ...)
api/               FastAPI server (Phase 2)
worker/            Celery background workers (Phase 5)
web/               React + TypeScript UI (Phase 3)
data/reference/    Reference datasets for verification
docs/              Documentation
```

## Commands

- Run tests: `uv run pytest`
- List modules: `uv run ortam list`
- Run a module standalone:
  `uv run ortam run groundwater.three-point modules/groundwater/examples/three_point.yaml`
- Add a dependency: `uv add <package>`

## Current status

- Phase 0 started (2026-11-07). Folder structure, rules, module contract and
  CLI in place. First standalone tool: `groundwater.three-point` (tested).
- Next task: D0.2 — development environment and CI.
