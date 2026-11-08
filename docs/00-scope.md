# ORTAM — Purpose and Scope

Version: 0.1 · Date: 2026-11-07

## 1. Purpose

ORTAM (Open Risk, Transport And Modelling) is a **modular environmental
modelling platform** for soil, groundwater, soil gas, air and the food chain.

The platform is used in two ways:

1. **Standalone modules:** every model can run on its own. A user may work
   only with the groundwater flow and transport model, only with AERMOD air
   dispersion, or only with geochemical speciation. No risk assessment,
   database or other module is required.
2. **Integrated workflows:** modules are connected through the Conceptual
   Site Model (source → transport pathway → exposure pathway → receptor).
   The risk chain, from sampling design to risk-based cleanup goals, is one
   such workflow; others (e.g. groundwater monitoring and trend analysis)
   can be added.

## 2. Scope

### 2.1 In scope (every module runs standalone)

| No | Module | Content |
|----|--------|---------|
| M1 | Project and site | Site boundary, coordinate reference systems, basemaps |
| M2 | Sampling design | Number and location of samples (VSP methods) |
| M3 | Field and laboratory data | Shallow and deep soil, groundwater, soil gas; GPS coordinates; QA/QC |
| M4 | Site characterization | Boring logs, potentiometric surface, flow direction, kriging, 3D geological model |
| M5 | Statistics | UCL95, background comparison, non-detects, trend analysis |
| M6 | Groundwater flow | MODFLOW 6, MODPATH 7 |
| M7 | Contaminant transport | MODFLOW 6 GWT, MT3D-USGS, Domenico analytical solution |
| M8 | Vadose zone | SSL leaching equation, VLEACH, Richards solver, MODFLOW 6 UZF |
| M9 | Vapor intrusion | Johnson & Ettinger (1991), EPA VISL |
| M10 | Air dispersion | AERMOD, AERMET, AERMAP; SSL VF and PEF |
| M11 | Surface runoff and erosion | RUSLE, HHRAP equations |
| M12 | Plant and food chain | EPA HHRAP (2005) |
| M13 | Geochemical speciation | PHREEQC (minteq.v4), MINTEQA2 |
| M14 | Exposure and risk | RAGS A/B/E/F; HQ, HI, ILCR; Monte Carlo |
| M15 | Remediation goals and strategy | Risk-based back-calculation, technology screening, remediation volume (m3) |
| M16 | Reporting | PDF and Word reports |
| M17 | NAPL | Raoult equilibrium and effective solubility, BIOSCREEN source-plume, HSSM (KOPT, OILENS, TSGPLUME) |
| M18 | Surface water | Watershed loading and water body (HHRAP B-4), rivers and lakes (Mills et al. 1985) |

### 2.2 Out of scope (version 1.0)

- Ecological risk assessment (to be considered in later versions)
- Radionuclides
- Multiphase flow (numerical, finite-difference/finite-element modelling of NAPL
  migration); NAPL equilibrium (Raoult) and the semi-analytical HSSM screening model
  are in scope (module `napl`, D1.5)
- Real-time sensor data streams

## 3. Target users

- Environmental engineers and hydrogeologists (characterization, modelling)
- Risk assessors
- Regulatory reviewers (read-only, report review)

## 4. Success criteria (version 1.0)

1. Every module can run standalone from the CLI, from Python and from the
   web UI, using only its own input file.
2. For an integrated workflow, risk results and cleanup goals can be produced
   in the same software with full traceability to their sources.
3. Every computational module is verified against the relevant EPA/USGS
   reference outputs.
4. Every result is reproducible from its samples, model run, parameters and
   toxicity table version.

## 5. Regulatory and standards framework

- Turkish Regulation on the Control of Soil Pollution and Point Source
  Contaminated Sites, and its technical guidance
- US EPA RAGS, RSL, VISL, SSL, HHRAP, ProUCL
- ASTM E1739, ASTM E2081

## 6. Assumptions and constraints

- Developed with Claude, module by module.
- No time constraint; quality and verifiability come first.
- All components are free/open source or public domain.
- The software, UI and documentation are in English; the UI is designed for
  internationalization.
- Web-based architecture; local (Docker) and desktop (Tauri) operation
  supported.
