# Vapor intrusion module: methods

## `vapor.visl` — vapor intrusion screening levels (task D1.2c)

| Quantity | Equation | Source |
|----------|----------|--------|
| Sub-slab / exterior soil gas target | C_sg = C_ia / AF_ss | VISL UG (2014) p. 8 |
| Groundwater target | C_gw = C_ia / (H' · AF_gw · 1000) | VISL UG p. 8 |
| H' at 25 °C | H' = HLC·1000 / (0.082057·298.15) | RSL UG 4.9.9, p. 228 |
| H' at groundwater temperature | Clausius-Clapeyron with ΔHv,S | RSL UG 4.9.9, p. 228 (the VISL UG print omits the /(R·T) step) |
| Pure-component vapor concentration (NVT check) | C_vp = VP/760 · 1e9 · MW / (R·T) | VISL UG p. 7 |
| Vapor at solubility (NVT check) | C_hc = S · 1000 · H' · 1000 | VISL UG p. 7 |

Defaults: AF_ss = 0.03, AF_gw = 0.001 (OSWER VI Technical Guide 2015,
Table 6-1, PDF p. 132; 0.0005 for laterally extensive fine-grained layers).
The VISL UG (2014) still lists 0.1 for soil gas; the 2015 guide supersedes it.
Target indoor air: risk.resident-air or risk.worker-air (indoor), connected
in `workflows/vapor_intrusion` (`visl_screening`, `visl_risk_from_subsurface`).

Verification: VISL calculator output (resident, THQ 0.1, saved 2026-11-08,
`docs/references/visl_calculator_resident_2026-11-08.xlsx`) for benzene,
PCE, TCE and vinyl chloride: IA, sub-slab, groundwater, C_vp, C_hc and the
forward risk/HQ from sub-slab concentrations, 4/4 to three significant figures.

## `vapor.johnson-ettinger` — J&E model, EPA implementation v6.0 (task D1.2c)

Source: J&E v6 = *Documentation for EPA's Implementation of the Johnson and
Ettinger Model ..., Version 6.0*, September 2017 (page numbers as printed).
Where the documentation leaves a detail open, the cell formula of EPA's
spreadsheet (`je_spreadsheet_v6_2017.xlsm`, worksheet MODEL) is followed.

| Quantity | Equation | Source |
|----------|----------|--------|
| Source vapor concentration | C_S = H'_S · C_gw · 1000 (groundwater); C_S = C_sg (soil gas) | Eq. 1, p. 12; MODEL!F14 |
| H' at T_S | H'_S = exp[−(ΔHv,S/R_c)(1/T_S − 1/T_R)] · H_c/(R·T_S) | Eq. 2, p. 12; MODEL!F27 |
| Enthalpy at T_S | ΔHv,S = ΔHv,b · [(1 − T_S/T_C)/(1 − T_B/T_C)]^n | Eq. 3, p. 13; n from Table 1, p. 14 |
| Capillary zone water-filled porosity | θw,cz = θr + (n − θr)/2^M | Eq. 4b, p. 15 |
| Capillary zone height | h_cz = 0.15/(0.2·D) cm | Eq. 7b, p. 19 |
| Effective diffusivity of a layer | D_eff = D_a·(n − n_w)^3.33/n² + (D_w/H'_S)·n_w^3.33/n² | Eq. 5b, p. 17; Eq. 8, p. 23 |
| Total effective diffusivity | D_T = L_T / Σ(L_i/D_i) | Eq. 9, p. 23; MODEL!G108 |
| A, B, C parameters | A = D_T(A_bf + 4L_b√A_bf)·0.36/(Q_b(L_s − L_b)); B = (Q_soil/Q_b)·Q_b·L_f/(D_A·η·(A_bf + 4L_b√A_bf)·0.36); C = Q_soil/Q_b | Eq. 13-15, p. 26 |
| Attenuation coefficient | dirt floor α = A/(1 + A); slab α = A/(1 + A·e^−B + (A/C)(1 − e^−B)) | Eq. 16a, 16b, p. 28 |
| Sub-slab source | α = Q_soil/Q_b | p. 60; MODEL!G93 |
| Q_soil/Q_b adjusted for ach | 0.003 · 0.45/ach | Eq. 17, p. 29 |
| Building ventilation rate | Q_b = A_bf · H_b · ach | Eq. 18, p. 29 |
| Indoor air; sub-slab | C_ia = α · C_S; C_ss = C_ia · Q_b/Q_soil | p. 60; MODEL!G96, G100 |

**Spreadsheet constants.** T = T_S + 273 (not 273.15), T_R = 298 K,
R = 0.000082057 atm·m3/(mol·K), R_c = 1.9872 cal/(mol·K). They are kept for exact
agreement with the spreadsheet; the VISL module uses 273.15 and 298.15 (RSL UG).
The difference in H' is below 0.1 %.

**Path geometry.** Strata A, B, C from the surface down; their thicknesses
must add up to L_s. The diffusion path starts at the base of the foundation,
so stratum A counts with h_A − L_b. For groundwater, the capillary zone (soil
of the deepest stratum) is taken out of the deepest stratum. Checks marked
"ERROR" in the spreadsheet (MODEL!J11-J85) are rejected; range checks give
warnings. A groundwater concentration above the water solubility gives an
invalid result (α is still reported, C_ia is not).

**Data.** `data/reference/je/v6-2017/` (imported by `modules/vapor/je_import.py`):
283 chemicals (CHEM_DATA, physical-chemical properties), 12 SCS soils
(SOIL_DATA), 5 residential buildings (BLDG_DATA). The J&E chemical table, not
the RSL table, is used for J&E, so that the model reproduces EPA's spreadsheet;
every property can be overridden per chemical.

**Verification.** 16 cases run in EPA's spreadsheet (LibreOffice recalculation,
`modules/vapor/tests/generate_je_reference_cases.py`): three source types,
five foundation types, 1-3 strata, 10-25 °C, six chemicals. H'_S, C_S, D_eff of
every layer, h_cz, D_T, A, B, α, C_ia and C_ss agree to a relative 1e-9. Case 16
is the case saved in the downloaded spreadsheet; its result equals the value
cached by Excel (α = 9.93074e-4).

**Risk (workflow `je_risk`).** C_ia is passed to risk.resident-air or
risk.worker-air (indoor). Against the spreadsheet risk output (MODEL!G153-G154,
16 cases): equal to 1e-6 except (i) vinyl chloride HQ, because the RSL Nov 2024
RfC is 0.0511 mg/m3 (J&E 2017: 0.1), and (ii) TCE risk, 1.2e-4 relative, because
the spreadsheet rounds the TCE IUR parts to 1.0e-6 and 3.1e-6 while the RSL uses
4.1e-6 × CAF 0.244 / MAF 0.756 (RSL UG 4.1.3.5, p. 101).
