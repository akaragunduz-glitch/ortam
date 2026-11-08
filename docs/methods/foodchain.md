# Food chain: methods (task D1.3b)

Sources: HHRAP 2005 = *Human Health Risk Assessment Protocol for Hazardous Waste
Combustion Facilities, Final*, EPA530-R-05-006 (Appendices A, B, C; printed page,
PDF page in brackets); HHRAP Companion Database (`docs/references/hhrap_2005_chemical_database.mdb`).

## Data

`ortam_core/chemicals/hhrap_import.py` copies the Access database (read with the
pure-Python package access-parser) into `data/reference/hhrap/2005/`:
`parameters.csv` (204 chemicals, 7344 values, with the database's explanation code
and reference for every value), `explanations.csv`, `references.csv`, `metadata.yaml`
(SHA-256). Default parameters with their table and page: `defaults.yaml`.
The database's 2005 toxicity values are not used; risk uses the RSL tables.

## `foodchain.media` — concentrations in food (HHRAP Appendix B)

| Food | Equation | Source |
|------|----------|--------|
| Aboveground produce, root uptake | Pr_ag = Cs·Br_ag | Table B-2-9, p. B-80 [83] |
| Exposed aboveground produce | Pd + Pv + Pr_ag (Pd, Pv optional input from an air model) | Tables B-2-7, B-2-8 |
| Belowground produce | Pr_bg = Cs·Br_rootveg·VG_rootveg; VG = 0.01 if log Kow > 4, else 1 | Table B-2-10, pp. B-83, B-85 [86, 88] |
| Forage, silage / grain | Pr = Cs·Br_forage / Cs·Br_grain | Table B-3-9, pp. B-134, B-135 [137-138] |
| Beef, milk, pork | A = (Σ F_i·Qp_i·P_i + Qs·Cs·Bs)·Ba·MF | Tables B-3-10 to B-3-12, pp. B-137 to B-154 |
| Eggs, chicken | A = (Σ F_i·Qp_i·P_i + Qs·Cs·Bs)·Ba | Tables B-3-13, B-3-14, pp. B-158 to B-165 |
| Fish | C_dw·BCF; C_dw·BAF; C_sb·f_lipid·BSAF/OC_sed | Tables B-4-26 to B-4-28, pp. B-267 to B-273 |

Defaults (`defaults.yaml`): Qp beef 8.8/2.5/0.47, milk 13.2/4.1/3.0 (forage/silage/grain),
pork 1.4/3.3 (silage/grain), eggs and chicken 0.2 (grain) kg DW/day; Qs 0.5, 0.4, 0.37,
0.022 kg/day; F = Bs = 1; MF = 0.01 for bis(2-ethylhexyl)phthalate, 1 otherwise (beef,
milk, pork only); f_lipid 0.07, OC_sed 0.04. Fish factor: BSAF (PCDD/F, PCBs) when a
sediment concentration is given, BAF when the database has one (log Kow ≥ 4), otherwise BCF
(App. A-2-2.13.4, p. A-2-28).

**Transfer factors.** From the database; for organic chemicals that are not in it, from the
Appendix A-2 correlations with the RSL log Kow and Koc (`modules/foodchain/btf.py`):
RCF Eq. A-2-14/15 (DW with 87 % moisture), Br_rootveg = RCF/Kds Eq. A-2-16, Br Eq. A-2-17/18
(Travis and Arms), Bv Eq. A-2-19/20 (÷100, database explanation 53), Ba Eq. A-2-21 to A-2-23
(RTI 2005; fat fractions milk 0.04, beef 0.19, pork 0.23, chicken 0.14, eggs 0.08); every
factor can be overridden. Metals, ionizing organics (pKa needed) and fish BCF (Meylan
correction factors) have database values only.

## `risk.foodchain-ingestion` — intake and risk (HHRAP Appendix C)

| Quantity | Equation | Source |
|----------|----------|--------|
| Produce | I_ag = [(Pd+Pv+Pr)·CR_ag + Pr·CR_pp + Pr_bg·CR_bg]·F_ag | Table C-1-2, p. C-6 [7] |
| Animal products | I_j = A_j·CR_j·F_j | Table C-1-3, p. C-10 [11] |
| Fish | I_fish = C_fish·CR_fish·F_fish | Table C-1-4, p. C-14 [15] |
| Cancer risk | I·ED·EF·CSF/(AT·365), AT = 70 yr | Table C-1-7, p. C-23 [24] |
| Hazard quotient | I·ED·EF/(RfD·AT·365), AT = ED | Table C-1-8, p. C-26 [27] |

Receptors: farmer, farmer child, resident, resident child, fisher, fisher child;
consumption rates kg/kg-day (Tables C-1-2 to C-1-4); ED 40/6/30/6/30/6 yr, EF 350 d/yr.
Doses are reported as CDI rows per food (D1.1d structure). HHRAP 2005 applies no ADAF;
mutagens are flagged.

## Workflow

`workflows/foodchain.foodchain_risk`: medium concentrations → `foodchain.media` →
`risk.foodchain-ingestion`; optionally plus direct soil contact from
`risk.resident-soil` (RSL equations).

## Verification

- `modules/foodchain/tests/test_btf.py`: the Appendix A-2 correlations recomputed from the
  database's own log Kow, H and Koc reproduce the database values that carry the matching
  explanation code: Br_ag, Br_forage, Br_grain 187 each; RCF 187; Br_rootveg 185; Bv 169
  (above-ground and forage); Ba beef/milk/pork/chicken/eggs 174 each; 2141 values in all (within 1 %; Br_rootveg
  5 %, a quotient of two rounded values). Two documented exceptions: OCDF RCF
  (database value corresponds to log Kow 8.2, its log Kow is 8.0) and 1,1-dichloroethane
  Br_rootveg (Kds stored 0.06, the value implies 0.0635).
- Appendix B and C equations by hand (benzene, arsenic, BEHP; all foods and receptors).
- RAIS farmer comparison: pending (D1.3b part 2).
