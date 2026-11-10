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
10. **Every completed deliverable updates the roadmap**
    `docs/ORTAM_Yol_Haritasi.xlsx` (in Turkish, for the project owner):
    - `Teslimatlar`: edit only the input cells (yellow, blue text) —
      Durum, İlerleme (0-1), Tamamlanma tarihi — plus Gerçekleşen sonuç and
      Commit; add a row for any new sub-deliverable (copy the formulas of
      columns F and M from the row above);
    - `Doğrulama`: add one row per verification (reference, compared,
      matched) and keep the total formulas covering all rows;
    - `Açık sorular`: add new questions, mark resolved ones `Çözüldü`;
    - never overwrite formula cells; recalculate (LibreOffice) and check
      that there are no formula errors before saving.

## Folder structure

```
ortam_core/        Shared core: module contract, units, CLI, chemical database
modules/           Standalone modules (each runs on its own)
  groundwater/     Three-point, analytical solutions, MODFLOW 6 flow + GWT, MODPATH
  vadose/          SSL leaching, VLEACH, Richards solver, UZF
  vapor/           Vapor intrusion: VISL screening levels, Johnson & Ettinger
  air/             AERMOD, emission factors
  surface/         RUSLE, surface runoff
  foodchain/       HHRAP plant and food chain uptake
  napl/            LNAPL: Raoult equilibrium, BIOSCREEN, HSSM
  surfacewater/    Watershed (HHRAP B-4), rivers and lakes (Mills et al. 1985)
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
  soil_volatiles/  air.ssl-vf -> risk soil tools (volatile chemicals)
  soil_to_groundwater/ risk.resident-tapwater -> vadose.ssl-leaching
  napl_plume/      napl.equilibrium -> napl.source-plume -> risk.resident-tapwater
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
- D0.1 complete: repository on GitHub (`akaragunduz-glitch/ortam`).
- D0.2 complete: local tests and CI green; Docker verified.
- D1.1 complete (2026-11-07): `risk.resident-soil`; RSL Nov 2024 data in
  `data/reference/rsl/2024-11/`; 533 chemicals match the RSL Summary Table.
- D1.1b complete (2026-11-07): mutagens, TCE, VC; `risk.resident-tapwater`
  and `risk.resident-air`. Shared equations in `modules/risk/equations/`.
- D1.1c complete (2026-11-07): `risk.worker-soil`, `risk.worker-air`,
  `risk.construction-soil`, `risk.recreator-soil`,
  `risk.recreator-surface-water`, `risk.fish`.
- D1.2a complete (2026-11-07): chemical database moved to
  `ortam_core/chemicals/`; `air.ssl-vf` (VF, VFsc, Csat), `air.ssl-pef`
  (PEF, PEFsc, PEF'sc); workflow `workflows/soil_volatiles` (volatile soil
  SLs match RSL 308/308 resident and industrial).
- D1.2b complete (2026-11-07): `vadose.ssl-leaching` (Method 1/2, DAF);
  workflow `workflows/soil_to_groundwater`; MCL-based SSL 91/91, risk-based
  SSL 741/747 (6 documented deviations).
- D1.2c complete (2026-11-08): `vapor.visl` (AF_ss 0.03, AF_gw 0.001, OSWER
  2015) and `vapor.johnson-ettinger` (EPA J&E v6.0, 2017; data in
  `data/reference/je/v6-2017/`); workflow `workflows/vapor_intrusion`
  (`visl_screening`, `visl_risk_from_subsurface`, `je_risk`). VISL calculator
  4/4; J&E spreadsheet 16/16 cases to 1e-9.
- D1.1d complete (2026-11-08): doses (CDI, DAD, EC) by route and age segment
  in all nine risk tools (`modules/risk/equations/dose.py`); RAGS E Exhibit B-3
  209/209 within rounding; risk/HQ from doses = tool risk/HQ for all chemicals.
- D1.3a complete (2026-11-08): `vadose.vleach` (VLEACH 2.2a, SI units; .INP
  converter); scheme identified from EPA's executable under Wine; 6 reference
  cases (EPA sample + 5 generated) within 5e-4; workflow
  `soil_to_groundwater.vleach_to_groundwater`.
- D1.3b part 1 (2026-11-08): HHRAP database import; `foodchain.media` (App. B),
  `risk.foodchain-ingestion` (App. C), workflow `workflows/foodchain`; App. A-2
  correlations reproduce the database (2139 of 2141 values).
- D1.3b complete (2026-11-08): `risk.farmer` (RAIS §4.8, food and soil
  scenarios) matches the RAIS calculator output for 8 chemicals; three RAIS
  internal deviations documented in `docs/methods/foodchain.md`.
- Milestone M1 deliverables complete (D1.1-D1.3).
- Scope extended 2026-11-08 (owner): milestone M1b, LNAPL (service stations) and
  surface water (rivers and lakes) before Phase 2. Order: D1.5a, D1.5b, D1.4a,
  D1.4b, D1.5c; then D2.1 (database and API).
- D1.5a complete (2026-11-08): `napl.equilibrium` (Raoult's law Ko HSSM Eq. 44-45,
  effective solubility, NAPL presence/composition/saturation in soil, PVI LNAPL
  indicators); HSSM RAOULT 9/9, BIOSCREEN 15/16 (1 source misprint), RSL Csat
  206/206.
- D1.5b complete (2026-11-10): `napl.source-plume` (BIOSCREEN 1.4: Domenico model, layer-cake
  source, source decay, no degradation / first order / instantaneous reaction, mass balance, 1.4
  mass flux, receptors, numerical plume mass); engine `modules/napl/bioscreen.py` in ft/yr, SI at
  the module boundary; `bioscreen_compatible` (default) reproduces the workbook incl. three
  inconsistencies Q1-Q3 (`docs/methods/napl.md` §4.6), false corrects them. BIOSCRN4.xls
  recalculated in LibreOffice for 14 cases: 163/163 to 1e-9 (data in
  `data/reference/bioscreen/v1.4/`); Examples 1 (Keesler) and 2 (Hill) match the printed output.
  Workflow `workflows/napl_plume`. Instantaneous model: default `instantaneous_source:
  mass-balance` (BC constant, source mass and concentration vanish together; §4.8),
  `bioscreen` keeps BIOSCREEN's depletion. Q1-Q3 are not reported to the developers (owner).
  Next: D1.4a `surfacewater.watershed` (HHRAP B-4).
