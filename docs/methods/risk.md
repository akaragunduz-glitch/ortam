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
- **Errors rather than wrong numbers:** volatile chemicals without a
  user-supplied VF (VF calculation: task D1.2) are rejected.

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

## Cancer equation sets (task D1.1b)

Selected per chemical (`cancer_method` in the output):

| Method | Applies to | Soil | Tap water | Air |
|--------|-----------|------|-----------|-----|
| standard | all other carcinogens | §4.1.1.3, pp. 40-41 | §4.1.2.3, pp. 70-75 | §4.1.3.2, p. 98 |
| mutagen | RSL mutagen flag "M" (ADAF 10/3/3/1 over ages 0-2, 2-6, 6-16, 16-26) | §4.1.1.4, pp. 42-45 | §4.1.2.4, pp. 75-81 | §4.1.3.3, p. 99 |
| vinyl chloride | CAS 75-01-4 (early-life term without averaging time) | §4.1.1.5, pp. 46-49 | §4.1.2.5, pp. 81-84 | §4.1.3.4, p. 100 |
| trichloroethylene | CAS 79-01-6 (CAF/MAF fractions: oral 0.804/0.202, inhalation 0.756/0.244) | §4.1.1.6, pp. 50-56 | §4.1.2.6, pp. 84-90 | §4.1.3.5, p. 101 |

Age-adjusted factors (`modules/risk/equations/factors.py`), checked against
the values printed in the guide: IFS 36 750 mg/kg, DFS 103 390 mg/kg,
IFSM 166 833.3 mg/kg, DFSM 428 260 mg/kg, IFW 327.95 L/kg, IFWM 1019.9 L/kg,
DFW 2 610 650 cm2-event/kg, DFWM 8 191 633 cm2-event/kg,
ET_event-adj = ET_event-madj = 0.67077 h/event.

## `risk.resident-tapwater` — resident tap water (task D1.1b)

Routes: ingestion; dermal contact (§4.9.8); inhalation of volatiles during
household use, volatile chemicals only, Andelman K = 0.5 L/m3 (§4.1.2, p. 62).

Dermal route (`modules/risk/equations/dermal_water.py`):

- **Inorganic** (no log Kow in the RSL Chemical Parameters table):
  SL = DA_event · 1000 / (Kp · ET_event).
- **Organic**: calculated only inside the effective predictive domain (EPD)
  and when FA > 0. If ET_event ≤ t*: SL = DA_event · 1000 /
  (2 · FA · Kp · √(6 · τ · ET_event / π)); otherwise SL = DA_event · 1000 /
  (FA · Kp · [ET_event/(1+B) + 2τ(1+3B+3B²)/(1+B)²]) (§4.1.2.1, p. 64).
- **FA**: regression in log B and log τ_event (§4.9.8, p. 221), clipped to
  [0, 1] and rounded to the nearest tenth (without rounding only 772 of 849
  chemicals match).
- B, τ_event and t* are taken from the RSL Chemical Parameters table.

**EPD — correction of the RSL User's Guide.** §4.9.8 (PDF p. 217) prints
both boundaries as −0.06831 ≤ −5.103·10⁻⁴·MW + 0.05616·log Kow ≤ 0.5577 / 0.1758.
The primary source, RAGS Part E (2004) p. 3-7, gives:

- Eq. 3.9: −0.06831 ≤ **+**5.103·10⁻⁴·MW + 0.05616·log Kow ≤ 0.5577
- Eq. 3.10: **−0.3010** ≤ −5.103·10⁻⁴·MW + 0.05616·log Kow ≤ 0.1758

With the RSL UG form, 766 (upper 0.5577) or 805 (upper 0.1758) of 849
chemicals match the RSL tap water table; with RAGS Part E, 840 match.
ORTAM uses RAGS Part E.

**Known deviations** (9 of 849, all within 30 %; listed in
`test_resident_tapwater.py`):

| Chemical | RSL ug/L | ORTAM ug/L | Reason |
|----------|----------|------------|--------|
| Aroclor 1242, 1248, 1254, 1260 | 0.0078 | 0.00785 | PCB special case (RSL UG §5.8) |
| Mercuric chloride (and other mercury salts) | 5.7 | 4.17 | special case (§5.19); RSL appears to treat it as inorganic |
| Auramine | 0.078 | 0.0766 | dermal parameters; open |
| Chlordane (alpha) | 3.6 | 3.30 | dermal parameters; open |
| 4,4'-Methylene-bis(N,N-dimethyl)aniline | 0.7 | 0.745 | dermal parameters; open |
| p-alpha,alpha,alpha-Tetrachlorotoluene | 0.0017 | 0.00155 | dermal parameters; open |

## `risk.resident-air` — resident air (task D1.1b)

Inhalation of ambient air; no assumption about the source of the air
concentration (§4.1.3, p. 96). Noncancer uses EDres (26 yr) and
AT = 365 · EDres; cancer uses AT = 365 · LT.

## Verification summary (RSL Summary Table, November 2024)

| Tool | Chemicals compared | Agreement (2 s.f.) |
|------|--------------------|--------------------|
| `risk.resident-soil` (non-volatile, incl. 21 mutagens) | 554 | 554 |
| `risk.resident-tapwater` | 849 | 840 (+ 9 documented deviations) |
| `risk.resident-air` | 365 | 365 |

Excluded from the comparisons: rows without toxicity values (lead:
regulatory levels, RSL UG §5.2) and refractory ceramic fibres (fibre units,
§5.23). Volatile chemicals in soil (including VC and TCE) need a VF (task
D1.2); their ingestion, dermal and inhalation equations are verified by hand
calculation.

### Planned extensions

| Task | Content |
|------|---------|
| D1.1b | Done: mutagens (ADAF), TCE, vinyl chloride; tap water and air |
| D1.1d | Doses (CDI, DAD, EC) as outputs, verified against RAGS A/E/F |
| D1.1c | Composite and construction worker, recreator (sediment, surface water), fish |
| D1.2 | VF and PEF calculation (SSL), Johnson & Ettinger |
| D1.3 | `foodchain` module (RAIS farmer, HHRAP) and food ingestion in `risk` |
