# Risk module: methods

## `risk.resident-soil` — resident soil (task D1.1)

Screening levels (SL), cancer risk and hazard quotients for a resident
exposed to soil by incidental ingestion, dermal contact and inhalation of
particulates (PEF) and, when a volatilization factor is supplied, volatiles
(VF).

### Sources

| Item | Source |
|------|--------|
| Noncancer SL, child (ingestion, dermal, inhalation, total) | RSL User's Guide (2026-03), §4.1.1.1, PDF pp. 36-37 |
| Noncancer SL, adult | RSL User's Guide, §4.1.1.2, PDF pp. 38-39 |
| Cancer SL, age-adjusted; IFSres-adj and DFSres-adj | RSL User's Guide, §4.1.1.3, PDF pp. 40-41 |
| Risk = C·TR/SL, HQ = C·THQ/SL | RSL User's Guide, §2.6.1, PDF p. 28 |
| One-hit rule, Risk = 1 − exp(−CDI·SF) above 0.01 | RSL User's Guide, §2.6.2, PDF p. 28; RAGS Part A (1989), p. 8-11 |
| Total risk, Risk_T = Σ Risk_i | RAGS Part A (1989), p. 8-12 |
| Hazard index, HI = Σ HQ_i | RAGS Part A (1989), p. 8-13 |
| Ceiling limit 1E+05 mg/kg (flagged, not applied) | RSL User's Guide, §5.13 |
| Default exposure parameters | RSL Default Exposure Parameters, Nov 2024 (`data/reference/rsl/2024-11/exposure_resident.yaml`) |
| Toxicity values, GIABS, ABSd | RSL Summary Table, Nov 2024 (`data/reference/rsl/2024-11/toxicity.csv`) |

The RSL User's Guide does not number its equations; code comments cite the
section and PDF page.

### Equations

Child noncancer (AT = 365 · EDres-c):

- SL_ing = THQ · AT · BW / [(RBA/RfD) · 10⁻⁶ · EF · ED · IRS]
- SL_derm = THQ · AT · BW / [(1/(RfD · GIABS)) · 10⁻⁶ · EF · ED · SA · AF · ABSd]
- SL_inh = THQ · AT / [(1/RfC) · EF · ED · ET · (1/24) · (1/VF + 1/PEF)]
- SL_total = 1 / (1/SL_ing + 1/SL_derm + 1/SL_inh)

Adult noncancer: the same with adult parameters and ED = EDres (26 yr),
AT = 365 · EDres.

Cancer (AT = 365 · LT):

- SL_ing = TR · AT / (SFO · 10⁻⁶ · RBA · IFSres-adj)
- SL_derm = TR · AT / [(SFO/GIABS) · 10⁻⁶ · DFSres-adj · ABSd]
- SL_inh = TR · AT / [IUR · 1000 · EFres · EDres · ETres · (1/24) · (1/VF + 1/PEF)]
- IFSres-adj = EFc·EDc·IRSc/BWc + EFa·EDa·IRSa/BWa = 36 750 mg/kg
- DFSres-adj = EFc·EDc·SAc·AFc/BWc + EFa·EDa·SAa·AFa/BWa = 103 390 mg/kg

Reported screening level: the lower of the cancer SL and the child
noncancer SL, as in the RSL Generic Tables (verified by the regression test).

### Decisions

- **THQ = 1 by default.** The Default Exposure Parameters workbook lists
  THQ = 0.1 as a user selection, but the RSL Summary Table used for
  verification is computed with THQ = 1.
- **RBA:** 0.6 for inorganic arsenic (CAS 7440-38-2), 1 otherwise
  (Default Exposure Parameters, Miscellaneous Variables, row 5).
- **Overrides:** any toxicity value, GIABS, ABSd, RBA, VF or exposure
  parameter may be supplied by the user; every override is listed in the
  output and in the warnings.
- **Standalone use:** `toxicity_table: none` runs without the RSL table;
  all values then come from the input file.
- **Errors rather than wrong numbers:** mutagens (task D1.1b) and volatile
  chemicals without a user-supplied VF (VF calculation: task D1.2) are
  rejected.

### Verification

`modules/risk/tests/test_resident_soil.py`:

- IFSres-adj and DFSres-adj against the RSL User's Guide values;
- Acephate route-by-route hand calculation;
- ten selected chemicals covering every route combination;
- **all 533 non-volatile, non-mutagenic chemicals with a resident soil RSL
  in the November 2024 Summary Table agree to two significant figures**;
- forward risk/HQ, one-hit rule, HI, overrides, error cases, CLI example.

`modules/risk/tests/test_rsl_import.py` checks that the committed data files
are reproduced exactly by the importer and that values and exposure defaults
match the EPA workbooks.

### Planned extensions

| Task | Content |
|------|---------|
| D1.1b | Mutagens (ADAF), TCE and vinyl chloride equations; tap water and air |
| D1.1c | Composite and construction worker, recreator (sediment, surface water), fish |
| D1.2 | VF and PEF calculation (SSL), Johnson & Ettinger |
| D1.3 | `foodchain` module (RAIS farmer, HHRAP) and food ingestion in `risk` |
