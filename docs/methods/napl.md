# NAPL methods

Module `napl` (light non-aqueous phase liquids, mainly petroleum fuels at service
stations). D1.5a: `napl.equilibrium`; D1.5b: `napl.source-plume` (BIOSCREEN 1.4).
Planned: `napl.hssm` (HSSM KOPT/OILENS/TSGPLUME, D1.5c).

Sources (printed page [PDF page]):

| Short name | Document |
|------------|----------|
| HSSM1 | Weaver, J.W., Charbeneau, R.J., Tauxe, J.D., Lien, B.K., Provost, J.B. (1994), The Hydrocarbon Spill Screening Model (HSSM) Volume 1: User's Guide, EPA/600/R-94/039a (`hssm_vol1_users_guide.pdf`) |
| HSSM2 | Charbeneau, R.J., Weaver, J.W., Lien, B.K. (1994), HSSM Volume 2: Theoretical Background and Source Codes, EPA/600/R-94/039b (`hssm_vol2_theory_source.pdf`) |
| BSCR | Newell, C.J., McLeod, R.K., Gonzales, J.R. (1996), BIOSCREEN Natural Attenuation Decision Support System, User's Manual Version 1.3, EPA/600/R-96/087, with the Version 1.4 revisions (July 1997) (`bioscreen_users_manual_1996.pdf`) |
| BSCR14 | BIOSCREEN Version 1.4 Revisions, July 1997 (pages 64-85 of `bioscreen_users_manual_1996.pdf`; also `BSCRN14.PDF` in `bioscreen_1.4.zip`) |
| XLS | BIOSCRN4.xls (`bioscreen_1.4.zip`): worksheet formulas, cited as sheet!cell; VBA "Macro Module" |
| PVI | EPA (2015), Technical Guide for Addressing Petroleum Vapor Intrusion at Leaking Underground Storage Tank Sites, EPA 510-R-15-001 (`epa_pvi_guide_2015.pdf`) |
| SSG / RSL UG | Soil Screening Guidance (2002), Eq. 4-9; RSL User's Guide §4.9.4, §5.12 |

## 1. NAPL/water partitioning (Raoult's law)

| Quantity | Equation | Source |
|----------|----------|--------|
| Linear partitioning | co = Ko · cw | HSSM1 Eq. 43, p. 172 [187] |
| Ko, composition basis | Ko = Tk · Σj(coj/Tj) / (sk · γk) | HSSM1 Eq. 44, p. 172 [187] |
| Ko, average molecular weight basis | Ko = 10^6 · (ρo/To) / (Sk/Tk) | HSSM1 Eq. 45, p. 173 [189] |
| NAPL density, average MW | ρo = Σj coj; To = Σj Tj·coj / Σj coj (mass-weighted) | HSSM2 RAOULT source, p. 256 [269] |
| Unit factors | ×1000 mg/g × 1000 cm3/L (coj in g/cm3, sk in mg/L) | HSSM2 p. 256 [269] |
| Effective solubility | S_eff,k = xk · Sk (mole fraction × pure-phase solubility) | BSCR p. 3 [68] |

With γk the activity coefficient of Eq. 44, co/Ko from Eq. 43-44 equals xk·γk·Sk, so the
effective solubility and Ko are two forms of the same relation (checked in the tests).
To in Eq. 45 is the mass-weighted mean molecular weight, as computed by the RAOULT program
(not the number-average); with it the RAOULT outputs are reproduced to 4 decimals.

Components without a CAS number (classes such as "other paraffins") need their own
molecular weight and solubility; solubility 0 marks an insoluble component (HSSM1
Table 106, RAOULT.DAT). Components that are solid at soil temperature (MP ≥ 20 °C, RSL
UG 5.12) are used with the solubility of the solid and flagged: the subcooled-liquid
solubility that Raoult's law would require is not given in the reference documents (open
question).

## 2. NAPL in soil: presence, composition and saturation

Unknown: whether a NAPL phase exists for measured total soil concentrations Ct,k (mg/kg)
of the NAPL components, and if so its composition xk and volume fraction θo.

Per litre of soil, the mass of component k is distributed over water, sorbed, gas and NAPL
phases. HSSM2 Eq. 25-26 (p. 16 [30]) write the bulk partition coefficient as
φB = θw + θo·ko + ρb·kd (air neglected); the SSG Csat equation (SSG Eq. 4-9, RSL UG 5.12)
uses the same three-phase denominator with the air term. ORTAM keeps the air term:

    B'k = θw + ρb·Kd,k + θa·H'k                     (SSG Eq. 4-9 denominator)
    Ct,k·ρb = cw,k·B'k + θo·co,k                    (HSSM2 Eq. 25-26 with the air term)

With cw,k = xk·γk·Sk (Eq. 43-44, BSCR p. 3) and θo·co,k = xk·Mk·No (No: mmol of NAPL per
litre of soil), in molar units (Nk = Ct,k·ρb/Mk, mmol/L):

    Nk = xk·ak + xk·No,       ak = γk·Sk·B'k/Mk
    xk = Nk/(ak + No),        Σk xk = 1

- At No → 0, Σk Nk/ak = Σk Ct,k/(γk·Csat,k), where Csat,k = Sk·B'k/ρb is the
  single-chemical Csat of SSG Eq. 4-9. **A NAPL phase exists if Σk Ct,k/Csat,k > 1**
  (ideal solution). For one component this is Ct > Csat (RSL UG 5.12); for a mixture,
  NAPL forms even when every component is below its own Csat.
- Otherwise No is the single root of Σk Nk/(ak + No) = 1 in (0, Σk Nk] (the left side
  decreases monotonically); solved with Brent's method (scipy).
- An insoluble component (S = 0) is entirely in the NAPL, so any amount of it forms a NAPL.
- NAPL mass per litre of soil: No·Σk xk·Mk; θo = that mass / ρo; So = θo/n.
- co,k = xk·Mk·No/θo; Ko,k then follows from Eq. 44 with the solved composition.
- Without NAPL: cw,k = Ct,k·ρb/B'k.

Assumptions: equilibrium; ideal solution unless γ is given; the NAPL volume does not
reduce θw and θa (a warning is given when θo > 0.1·θa); Kd = Koc·foc for organics (RSL UG
4.9.4); n = 1 − ρb/ρs, θa = n − θw (RSL UG 4.9.4). The NAPL density is a user input (for
gasoline HSSM1 p. 173 quotes 0.74 g/mL from Cline et al. 1991). If the residual saturation
Sor is given (site data), So > Sor is reported as "potentially mobile".

The derivation combines the cited equations; it is not itself printed in a reference
document. It is verified by (a) the single-component limit, which must equal the RSL Csat,
and (b) a forward calculation of a known NAPL in soil with HSSM2 Eq. 25-26 (tests).

## 3. LNAPL indicators (EPA 2015 PVI guide)

PVI p. 58 [64] (`data/reference/napl/lnapl_indicators.yaml`):

| Indicator | Threshold |
|-----------|-----------|
| Benzene in groundwater | 5 mg/L |
| TPH in groundwater | 30 mg/L |
| Benzene in soil | 10 mg/kg |
| TPH in soil, unweathered gasoline | 100 mg/kg |
| TPH in soil, diesel or weathered gasoline | 250 mg/kg |
| Groundwater / effective solubility (Bruce et al. 1991, footnote 29) | 0.2 |

The guide notes there is no precise threshold and recommends several lines of evidence;
the module reports each indicator separately and draws no combined conclusion.

## 4. Dissolved plume from the source zone: `napl.source-plume` (BIOSCREEN 1.4)

Engine `modules/napl/bioscreen.py` (BIOSCREEN units ft, yr, mg/L, kg, L); module
`modules/napl/source_plume.py` (SI input and output; 1 ft = 0.3048 m, 1 yr = 365 d, exact).
Workflow `workflows/napl_plume.napl_plume_risk`: napl.equilibrium (effective solubility as the
central source-zone concentration, optional) → napl.source-plume → risk.resident-tapwater.

### 4.1 Input relations

| Quantity | Equation | Source |
|----------|----------|--------|
| Seepage velocity | Vs = K·i/n | BSCR p. 16 [22]; XLS Input!E5 |
| Longitudinal dispersivity | αx = 3.28·0.83·(log10(Lp/3.28))^2.414 (ft) | Xu and Eckstein (1995); BSCR p. 19 [25], p. 53 [59]; Input!E12 |
| Transverse, vertical | αy = 0.1·αx; αz = 0 (replaced by 1e-99 ft) | BSCR p. 19 [25]; Input!E13, E73 |
| Retardation | R = 1 + ρb·Koc·foc/n | BSCR p. 21 [27]; Input!E19 |
| Decay coefficient | λ = 0.693/t½ | BSCR p. 23 [29]; Input!E26 |
| Biodegradation capacity | BC = Σ C(ea)n/UFn; UF: O2 3.14, NO3 4.9, Fe2+ 21.8, SO4 4.7, CH4 0.78 | BSCR14 p. 8 [72]; BSCR p. 25-27 [31-33]; Input!F47-F52 |

### 4.2 Source

Up to three zones ("layer cake", BSCR p. 47 [53], p. 30 [36]): the input lists them from the centre
outward (BIOSCREEN zones 3, 2/4, 1/5). The superposed layers have full widths
Y1 = 2Yo1 + 2Yo2 + Yo3, Y2 = 2Yo2 + Yo3, Y3 = Yo3 and concentration increments C1 = Co1,
C2 = Co2 − Co1, C3 = Co3 − Co2 (Input!L45-L47, Q45-Q47).

| Quantity | Equation | Source |
|----------|----------|--------|
| Flow through the source | Q = Vs·n·Z·W (W = largest layer width) | BSCR p. 41 [47]; Input!F59 |
| Average source concentration | Cavg = (2·Yo1·Co1 + 2·Yo2·Co2 + Yo3·Co3)/W | Input!F60 |
| Source decay | ks = Q·Cavg/M0; C(t) = C0·e^(−ks·t); t½ = 0.693/ks | BSCR App. A.3 p. 50-51 [56-57], p. 33 [39]; Input!F62-F63 |
| Instantaneous model | C0 → C0 + BC; ks,inst = Q·(Cavg + BC)/M0 | BSCR p. 35 [41], p. 47 [53]; Input!F65-F68, T45-T47 |
| Source time | concentration at (x, t) uses the source at t − x/v, v = Vs/R, not below 0 | BSCR p. 35 [41] |

### 4.3 Domenico model (z = 0)

BSCR14 p. 8 [72] (replaces BSCR p. 46 [52]); vertical dispersion downward only (BSCR14 p. 7 [71]):

    C = Σk Ck·e^(−ks(t − x/v)) / 8 · Fx · Fy,k · Fz
    Fx (no decay)    = erfc[(x − v·t) / (2·(αx·v·t)^½)]
    Fx (first order) = exp[x/(2αx)·(1 − (1 + 4λαx/v)^½)] · erfc[(x − v·t·(1 + 4λαx/v)^½) / (2·(αx·v·t)^½)]
    Fy,k = erf[(y + Yk/2)/(2·(αy·x)^½)] − erf[(y − Yk/2)/(2·(αy·x)^½)]
    Fz   = erf[Z/(2·(αz·x)^½)] − erf[−Z/(2·(αz·x)^½)]

Instantaneous reaction (superposition, BSCR App. A.2 p. 48-49 [54-55]): the no-decay model with
source C0 + BC, minus BC, not below zero ('Plume Output'!C119).

Workbook implementation details reproduced: the source column is evaluated at x = 1e-4 ft and
there C = Σ Ck·e^(−ks·t) over the layers that contain y ('Plume Output'!C63); the first-order
model is computed as C(no decay)·Fx(λ)/Fx(0) ('Centerline Output'!C60); the error function is
the workbook's own `gerf` (a Chebyshev approximation of erfc with fractional error < 1.2e-7),
not Excel's.

### 4.4 Mass balance and mass flux (11 × 5 array)

Grid x = 0, L/10, …, L and y = 0, ±W/4, ±W/2 (BSCR p. 38 [44]).

| Quantity | Equation | Source |
|----------|----------|--------|
| Mass that has left the source | M0·(1 − e^(−ks·t)), or Q·Cavg·t for an infinite source | BSCR p. 39 [45]; 'Plume Output'!P88, P111 |
| Plume mass on the grid | (ΣC − ½·ΣC(x=0))·Δx·Δy·Z·n·R | BSCR p. 39 [45]; P91, P104, P116, P121 |
| Fraction biodegraded | (M_nodecay − M_model)/M_nodecay on the grid | P107, P125 |
| Actual plume mass | mass left source × (1 − fraction) | BSCR p. 39 [45]; M19 |
| "Can't Calc." | grid mass outside 50-150 % of the mass left in the source | BSCR p. 40 [46]; M17 |
| Electron acceptor change | removed mass × C(ea)k/BC, − for O2, NO3, SO4; + for Fe2+, CH4 | BSCR p. 40 [46]; I25-M25 |
| Plume volume | number of cells with C > 0.005 mg/L × (Δx − 1e-4 ft)·Δy·Z·n | BSCR p. 41 [47]; M30 |
| Mass flux (1.4) | Vs·n·Z·Δy·(C(0) + 2C(Δy) + 2C(2Δy)); x = 0: Q·Cavg; not with αz > 0 | BSCR14 p. 2 [66]; B11-L11 |

In addition the module integrates C(x, y) numerically on a fine grid over the modelled area
(trapezoidal rule, default 201 × 201 points): plume mass n·Z·R·∬C dx dy and mass flux
Vs·n·Z·∫C dy at the array x positions (αz = 0 only). BIOSCREEN itself warns that its 5-point
mass flux is accurate to about ±50 % (BSCR14 p. 2 [66]).

### 4.5 Receptors

Concentrations of the three models at user points (x, y) and, optionally, as a time series.
The Domenico solution is evaluated at the water table (z = 0) as in BIOSCREEN.

### 4.6 Workbook inconsistencies (bioscreen_compatible)

`bioscreen_compatible: true` (default) reproduces BIOSCREEN 1.4 exactly, including its constants
(28.3 L/ft3, 0.693, 3.28 ft/m, `gerf`) and three inconsistencies of the workbook found while
transcribing its formulas. `false` uses exact constants and erf and corrects them:

- **Q1** Instantaneous model with a finite source: the middle layer (C2) decays with ks instead of
  ks,inst. In 'Centerline Output'!B80 and 'Plume Output'!C113-C115 the zone-2 term is
  bioc2·e^(−ks,inst·τ)/(C2·e^(−ks,inst·τ))·C(zone 2), and C(zone 2) contains C2·e^(−ks·τ);
  zones 1 and 3 use e^(−ks·τ) in the denominator and the factor cancels correctly.
- **Q2** First-order plume mass without the retardation factor ('Plume Output'!P104), while the
  no-decay and instantaneous masses include R (P91, P116, P121): for R > 1 the first-order
  fraction biodegraded is overstated.
- **Q3** BC is added to a layer whenever the next outer increment is zero (Input!T46-T47), so with
  equal concentrations in two zones it is added twice. The SI input requires positive zone
  concentrations; Q3 can still occur when two zones have the same concentration.

Further observations: the manual text says 1837 kg remain in the source of Example 1, the
output screen (Fig. 4) and the workbook 1873.7 kg; the Example 2 data table gives
K = 8.05E-03 cm/s and R = 1.3, the input screen (Fig. 6) 8.1E-03 cm/s (Vs 1609.1 ft/yr) and
R = 1.2; the plume volume counts cells above 0.005 mg/L (M30), the manual says "> 0"; the mass
flux at x = 0 is the source flux Q·Cavg for all three models.

### 4.7 Verification

| Reference | Compared | Agreement |
|-----------|----------|-----------|
| BIOSCRN4.xls recalculated in LibreOffice (14 cases: Examples 1-2, one to three zones, finite/infinite source, R up to 57, αz > 0, BC = 0, t = 1-100 yr; `data/reference/bioscreen/v1.4/reference_cases.json`): source terms, centreline, 3 arrays, grid masses, mass balance screens and mass flux of the three models | 14 cases, 163 tests | relative 1e-9 |
| Excel values cached in the distributed BIOSCRN4.xls vs the LibreOffice recalculation | 39 values | relative 1e-12 |
| Example 1, Keesler AFB (BSCR14 Fig. 3, 4, 6 [82-85]): centreline (33 values), masses, electron acceptors, source mass, volume, flow, mass flux (W = 200 ft and 50 ft), array W = 50 ft | all printed values | printed precision |
| Example 2, Hill AFB (BSCR Fig. 7, 8 [98-99], Version 1.3 output) | 33 centreline values, array row | ≤ 0.001 mg/L |
| Steady state, infinite source: closed-form limit of the Domenico solution | 3 | relative 1e-12 |
| No decay, constant source: integrated plume mass vs Q·C0·t; mass flux at L/10 vs Q·C0 | 2 | 1 %, 2 % |

The LibreOffice recalculation needs the workbook's VBA functions `gerf`/`gerfC`, which LibreOffice
imports but does not resolve from worksheet formulas; the generator script copies the two
functions verbatim from the workbook's "Macro Module" into the document's Standard library
(`modules/napl/tests/generate_bioscreen_reference_cases.py`).

## 5. Verification of napl.equilibrium

| Reference | Compared | Agreement |
|-----------|----------|-----------|
| HSSM1 Table 105 (p. 185 [200]): density 0.7168, average MW 104.0458, benzene Ko 311.6757 and 307.0647 | 4 | 4 decimals |
| HSSM1 Table 107 (p. 187 [202]): benzene S = 1780 mg/L, Ko 306.4227 and 301.8895 | 2 | 4 decimals |
| HSSM1 Table 100 (p. 174 [189]): calculated Ko benzene (312), toluene (1202), o-xylene (4440) | 3 | 3 significant figures |
| BSCR pp. 4-5 [68-70]: effective solubility of BTEX in fresh, two weathered gasolines and JP-4 | 16 | 15 (see below) |
| RSL Summary Table Nov 2024: Csat of 206 organics | 206 | 3 significant figures |
| Forward calculation, HSSM2 Eq. 25-26, θo = 1e-4, 0.01, 0.05 | 3 | relative 1e-9 |

BIOSCREEN inconsistency: for weathered gasoline #2 the table prints a xylene mole fraction
of 0.084 with S = 198 mg/L but an effective solubility of 15 mg/L (0.084 × 198 = 16.6;
15 mg/L corresponds to x = 0.076). The 15 mg/L is consistent with the printed BTEX total
(44 mg/L); the mole fraction is probably a misprint. The test keeps this case as a known
discrepancy.
