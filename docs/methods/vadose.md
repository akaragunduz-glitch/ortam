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
