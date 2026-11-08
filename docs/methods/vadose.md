# Vadose zone module: methods

## `vadose.ssl-leaching` — soil screening level protective of groundwater (task D1.2b)

| Quantity | Equation | Source |
|----------|----------|--------|
| Method 1 (partitioning, infinite source) | SSL = Cw·[Kd + (θw + θa·H')/ρb] | RSL UG §4.8.3, PDF p. 190; SSG 2002 Eq. 4-10 |
| Method 2 (mass-limit) | SSL = Cw·I·ED / (ρb·ds) | RSL UG §4.8.4, p. 193; SSG Eq. 4-14 |
| Target leachate concentration | Cw = SL(µg/L)/1000·DAF | RSL UG p. 190 |
| Dilution attenuation factor | DAF = 1 + K·i·d / (I·L) | RSL UG §4.8.5, p. 194; SSG Eq. 4-11 |
| Mixing zone depth | d = (0.0112·L²)^0.5 + da·[1 − exp(−L·I/(K·i·da))] | RSL UG p. 194; SSG Eq. 4-12 |
| Porosities | n = 1 − ρb/ρs; θa = n − θw | RSL UG p. 190 |

Defaults (`data/reference/rsl/2024-11/soil_to_groundwater.yaml`, RSL
workbook sheet "Soil to Groundwater"): DAF 1, θw 0.3, ρb 1.5 kg/L, ρs 2.65
kg/L, foc 0.002 (subsurface soil), I 0.18 m/yr, ED 70 yr. A DAF of 20 may be
used for a source up to 0.5 acre (RSL UG §4.8.5).

**Kd and H'.** Inorganics: table Kd (SSG 1996 Exhibit C-4 at pH 6.8, or Baes
1984; listed in RSL UG §4.8, pp. 185-188). Organics: Kd = Koc·foc.
Inorganic = no log Kow in the RSL table. H' = 0 for inorganics except
mercury (SSG Eq. 4-10, note).

**Target concentration.** A user value (µg/L) or the MCL from the RSL table
(`target: mcl`). Some MCL cells carry a note (e.g. "6.0E+01(G)"); the numeric
part is used and the note is reported. The module does not call the risk
module (rule 5): the risk-based SSL is produced by the workflow
`workflows/soil_to_groundwater` (risk.resident-tapwater → vadose.ssl-leaching),
which uses the unrounded tap water SL.

**Mixing zone depth.** SSG 2002 gives no cap of d at the aquifer thickness
(the 1996 SSG User's Guide, not in `docs/references/`, does); ORTAM warns
when d > da and does not change d.

No Csat substitution (RSL UG §4.8, p. 184).

### Verification (RSL Summary Table, November 2024; DAF = 1)

| Column | Chemicals | Agreement (2 s.f.) |
|--------|-----------|--------------------|
| MCL-based SSL | 91 | 91 |
| Risk-based SSL (workflow) | 747 | 741 (+ 6 documented deviations, all < 15 %) |

Deviations: the five tap water deviations of `risk.resident-tapwater`
(auramine, chlordane alpha, 4,4'-methylene-bis aniline, Aroclor 1254,
tetrachlorotoluene) and 3,3'-dichlorobenzidine (0.00085 vs 0.00082 mg/kg;
open question: check Koc and the tap water dermal parameters against the RSL
calculator output).

## `vadose.vleach` — VLEACH 2.2a one-dimensional leaching (task D1.3a)

Source: VLEACH UG = Ravi and Johnson, *VLEACH, A One-Dimensional Finite
Difference Vadose Zone Leaching Model, Version 2.2* (EPA CSMOS), printed pages.

| Quantity | Equation | Source |
|----------|----------|--------|
| Three-phase partitioning | C_T = θ·C_l + (n − θ)·C_g + ρ_b·C_s; C_g = H·C_l; C_s = K_d·C_l | VLEACH UG Eq. 4-7, pp. 9-10 |
| Distribution coefficient | K_d = K_oc·f_oc | Eq. 8, p. 10 |
| Gas diffusion coefficient | D = D_air·(n − θ)^(10/3)/n² | Eq. 9, p. 10 |
| Liquid advection | ∂C_l/∂t = −(q/θ)·∂C_l/∂z | Eq. 2, p. 9 |
| Gas diffusion | ∂C_g/∂t = D·∂²C_g/∂z² | Eq. 3, p. 9 |

**Numerical scheme.** The guide says only "implicit finite difference"
(p. 10). The scheme was identified by numerical experiments with EPA's own
executable (VleachForGUI.exe, extracted from the VLEACH 2.2a setup and run under
Wine) and then reproduced:
1. upwind Crank-Nicolson advection of the liquid phase only,
   w = q·Δt/(θ·Δz), inflow concentration = recharge water;
2. fully implicit diffusion of the gas phase only (storage n − θ); the boundary
   cells exchange with the atmosphere gas concentration, or with H·C_gw at the
   water table, over a distance Δz; a negative input value (null in ORTAM) makes
   the boundary impermeable;
3. both mass changes are computed from the same old state and added (additive
   operator splitting); each cell is then re-equilibrated (Eq. 5-7).

Findings:
- The 2.2b release notes say the Millington exponent was changed to 7/3 in
  version 2.2; the 2.2a executable uses 10/3, as in the guide (Eq. 9).
- VLEACH converts D_air from m2/day to ft2/yr with 3.2809 ft/m (0.7 m2/day →
  2750.3 ft2/yr in the .PRM file). In metres this is a factor 1.0000366 on D,
  kept so that ORTAM reproduces the executable.
- The solubility is only checked (warning), not enforced, as in VLEACH.
- Groundwater impact: rate in the last time step before each printout;
  cumulative mass = sum over all time steps (GWIMP.DAT holds every step).
- The EPA sample problem uses Δt = 100 yr with Δz = 1 ft (Courant number 100);
  ORTAM warns when q·Δt/Δz > 1.

Input in SI units (m, m/yr, g/cm3, ug/kg, mg/L, m2/day); `modules/vadose/vleach_inp.py`
converts an original VLEACH .INP file (feet).

**Verification** (`modules/vadose/tests/test_vleach.py`): six cases in
`data/reference/vleach/2.2a/` — EPA's sample problem from the setup and five
cases generated with the executable (`generate_vleach_reference.py`): recharge
concentration and both gas boundaries, impermeable boundaries with fine cells,
three polygons, diffusion only, advection only. Phase masses, boundary fluxes,
liquid/gas/sorbed profiles and the groundwater rate of every time step agree
within 5e-4 relative (≤ 1e-4 except the diffusion-only case, 3e-4; the
executable computes in single precision).

**Workflow.** `workflows/soil_to_groundwater.vleach_to_groundwater`: leachate
concentration C_L = rate/(area·q) (mg/L) and groundwater concentration
C_L/DAF (RSL UG 4.8.5; SSG 2002 Eq. 4-11), peak, first exceedance of a target.
