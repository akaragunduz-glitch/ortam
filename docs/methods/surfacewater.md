# Surface water methods

Module `surfacewater` (milestone M1b). D1.4a: `surfacewater.watershed` (HHRAP 2005 Appendix B,
Tables B-4-1 to B-4-25). Planned: D1.4b river and lake models (Mills et al. 1985).
Workflows: `workflows/surface_water` (`plume_to_water_body`, `water_body_fish_risk`).

Source: U.S. EPA (2005), Human Health Risk Assessment Protocol for Hazardous Waste Combustion
Facilities, Final, EPA530-R-05-006 (HHRAP): Appendix B (`hhrap_2005_appendix_b.pdf`), Chapter 5
(`hhrap_2005_ch5.pdf`). Pages: printed page [PDF page]; in Appendix B, PDF = printed + 3.
Chemical parameters (Kds, Kdsw, Kdbs, H, Dw, Da, ksg, Fv): HHRAP Companion Database
(`data/reference/hhrap/2005/parameters.csv`), overridable. Defaults with their table and page:
`data/reference/hhrap/2005/defaults.yaml` (section `surface_water`).

## 1. Watershed soil (Tables B-4-1 to B-4-6)

The watershed soil concentration is either measured (mg/kg) or calculated from air deposition.

| Quantity | Equation | Source |
|----------|----------|--------|
| Deposition term | Ds = 100·Q/(Zs·BD)·[Fv·Dytwv + Dytwp·(1 − Fv)] | B-4-1, p. B-169 [172] |
| Noncarcinogens | CstD = Ds·[1 − exp(−ks·tD)]/ks | B-4-1, p. B-169 [172] |
| Carcinogens, T1 < tD < T2 | Cs = [(Ds·tD − CstD)/ks + (CstD/ks)(1 − exp(−ks(T2 − tD)))]/(T2 − T1) | B-4-1, p. B-168 [171] |
| Carcinogens, T2 ≤ tD | Cs = Ds/(ks(tD − T1))·([tD + exp(−ks·tD)/ks] − [T1 + exp(−ks·T1)/ks]) | B-4-1, p. B-168 [171] |
| Soil loss constant | ks = ksg + kse + ksr + ksl + ksv | B-4-2, p. B-176 [179] |
| Erosion | kse = 0.1·Xe·SD·ER/(BD·Zs)·Kds·BD/(θsw + Kds·BD); default kse = 0 | B-4-3, p. B-180 [183] |
| Runoff | ksr = RO/(θsw·Zs)·1/(1 + Kds·BD/θsw) | B-4-4, p. B-185 [188] |
| Leaching | ksl = (P + I − RO − Ev)/(θsw·Zs·[1.0 + BD·Kds/θsw]) | B-4-5, p. B-189 [192] |
| Volatilization | ksv = [3.1536e7·H/(Zs·Kds·R·Ta·BD)]·(Da/Zs)·[1 − BD/ρsoil − θsw] | B-4-6, p. B-194 [197] |

With deposition-derived soil both bases are reported: Cs (carcinogens, averaged over the exposure
duration) and CstD (noncarcinogens, highest annual average) (Table B-4-1 description, p. B-168 [171]).
Watershed (Dytwv, Dytwp) and water body (Dytwv, Dytwp, Cywv) unitized air model results are separate
inputs: Ds and LRI use the watershed values, LDEP and Ldif the water body values.

## 2. Loads to the water body (Tables B-4-7 to B-4-14)

| Quantity | Equation | Source |
|----------|----------|--------|
| Total load | LT = LDEP + Ldif + LRI + LR + LE (+ external loads) | B-4-7, p. B-199 [202] |
| Direct deposition | LDEP = Q·[Fv·Dytwv + (1 − Fv)·Dytwp]·Aw | B-4-8, p. B-201 [204] |
| Impervious runoff | LRI = Q·[Fv·Dytwv + (1.0 − Fv)·Dytwp]·AI | B-4-9, p. B-204 [207] |
| Pervious runoff | LR = RO·(AL − AI)·Cs·BD/(θsw + Kds·BD)·0.01 | B-4-10, p. B-207 [210] |
| Erosion | LE = Xe·(AL − AI)·SD·ER·Cs·Kds·BD/(θsw + Kds·BD)·0.001 | B-4-11, p. B-211 [214] |
| Vapor diffusion | Ldif = Kv·Q·Fv·Cywv·Aw·1e-6/(H/(R·Twk)) | B-4-12, p. B-215 [218] |
| USLE | Xe = RF·K·LS·C·PF·907.18/4047 | B-4-13, p. B-219 [222] |
| Sediment delivery ratio | SD = a·AL^(−b); a from the area (≤ 0.1, 1, 10, 100, > 100 sq mi: 2.1, 1.9, 1.4, 1.2, 0.6), b = 0.125 | B-4-14, pp. B-223 to B-224 [226-227] |

**External loads (ORTAM extension).** HHRAP has no groundwater pathway to surface water. The
module adds user loads (g/yr) to LT, e.g. the mass flux of a groundwater plume discharging to a
gaining stream (`napl.source-plume`; workflow `plume_to_water_body`) or a point discharge. As
Table B-4-7 is a sum of independent loads and B-4-15 is linear in LT, the addition is consistent
with the HHRAP water body model; it is not itself in HHRAP.

## 3. Water body (Tables B-4-15 to B-4-25)

| Quantity | Equation | Source |
|----------|----------|--------|
| Total water body | Cwtot = LT/(Vfx·fwc + kwt·Aw·(dwc + dbs)) | B-4-15, p. B-226 [229] |
| Fractions | fwc = (1 + Kdsw·TSS·1e-6)·dwc/dz / [(1 + Kdsw·TSS·1e-6)·dwc/dz + (θbs + Kdbs·CBS)·dbs/dz]; fbs = 1 − fwc; dz = dwc + dbs | B-4-16, pp. B-230, B-232 [233, 235] |
| Dissipation | kwt = fwc·kv + fbs·kb | B-4-17, p. B-235 [238] |
| Volatilization | kv = Kv/(dz·(1 + Kdsw·TSS·1e-6)) | B-4-18, p. B-237 [240] |
| Overall transfer | Kv = [KL⁻¹ + (KG·H/(R·Twk))⁻¹]⁻¹·θ^(Twk − 293) | B-4-19, p. B-241 [244] |
| Liquid phase, river | KL = (1e-4·Dw·u/dz)^½·3.1536e7 | B-4-20, p. B-244 [247] |
| Liquid phase, lake | KL = (Cd^0.5·W)(ρa/ρw)^0.5(k^0.33/λz)(μw/(ρw·Dw))^−0.67·3.1536e7 | B-4-20, p. B-244 [247] |
| Gas phase | river KG = 36500 m/yr; lake KG = (Cd^0.5·W)(k^0.33/λz)(μa/(ρa·Da))^−0.67·3.1536e7 | B-4-21, p. B-249 [252] |
| Benthic burial | kb = ((Xe·AL·SD·1e3 − Vfx·TSS)/(Aw·TSS))·(TSS·1e-6/(CBS·dbs)); kb < 0 → 0 | B-4-22, p. B-253 [256]; Ch. 5, p. 5-81 [82] |
| Water column | Cwctot = fwc·Cwtot·(dwc + dbs)/dwc | B-4-23, p. B-257 [260] |
| Dissolved | Cdw = Cwctot/(1 + Kdsw·TSS·1e-6) | B-4-24, p. B-260 [263] |
| Bed sediment | Csb = fbs·Cwtot·(Kdbs/(θbs + Kdbs·CBS))·(dwc + dbs)/dbs | B-4-25, p. B-263 [266] |

kb > 1: warning to re-evaluate Xe (Chapter 5, p. 5-81 [82]). If H = 0, Kv = 0 and Ldif = 0.

## 4. Mercury

For total mercury (CAS 7439-97-6) the module follows the mercury notes of the tables: the emission
rate is 0.48·Q (loss to the global cycle) with Fv = 0.85 in Ds, LDEP, LRI and Ldif (B-4-1, -8, -9,
-12); Ds is split into 98 % Hg²⁺ and 2 % methyl mercury, each with its own soil parameters
(mercuric chloride, methyl mercury) and soil loss constant (B-4-1, p. B-169 [172]); LR and LE are
the sums of the two species (B-4-11, p. B-211 [214]); the water body uses the mercuric chloride
parameters and Cdw is split into 85 % Hg²⁺ and 15 % MHg (B-4-24, p. B-260 [263]). Measured total
soil mercury is split 98/2 % like the deposition (ORTAM assumption; noted in the output).
The workflow `water_body_fish_risk` passes both species to `foodchain.media`.

## 5. Verification

HHRAP gives no worked numerical example for these tables. The reference is an independent hand
calculation: `data/reference/hhrap/2005/b4_hand_calculation.ods`, built by
`modules/surfacewater/tests/build_b4_hand_calculation.py`, in which every equation is typed as a
worksheet formula from the printed table (one column per case, table and page in column C) and
recalculated by LibreOffice; the values are stored in `b4_hand_calculation.json`.

| Reference | Compared | Agreement |
|-----------|----------|-----------|
| Hand calculation, 6 cases: river/lake; measured soil (organic, inorganic); deposition with T2 ≤ tD, T2 > tD, CstD basis; external load; calculated kse, irrigation, negative kb | 177 values (all intermediate and final quantities) | relative 1e-11 |
| Hand calculation, mercury (two soil species and the water body) | 16 values | relative 1e-11 |
| Analytical: continuity of the two Cs branches at T2 = tD; ks → 0 limits (CstD → Ds·tD, Cs → Ds·tD/2); Kv limits (liquid- and gas-film control); no-loss dilution; USLE units; SD table | 12 | 1e-4 to exact |
| Workflow: plume flux = external load; fish risk = foodchain workflow | 6 | relative 1e-12 |

Observations: Table B-4-13 gives K = 0.39 ton/acre as the default cited from NC DEHNR (1997) and
U.S. EPA (1994), while its reference list (p. B-222 [225]) says U.S. EPA (1994) recommends 0.36;
0.39 is used. In B-4-15 the flow term is Vfx·fwc·Cwtot, i.e. Vfx·Cwctot·dwc/dz: the outflow uses
the water-column mass per total depth (a difference of dbs/dz, < 3 % for dwc ≥ 1 m).
