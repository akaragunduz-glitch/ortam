# NAPL methods

Module `napl` (light non-aqueous phase liquids, mainly petroleum fuels at service
stations). D1.5a: `napl.equilibrium`. Planned: `napl.source-plume` (BIOSCREEN 1.4,
D1.5b) and `napl.hssm` (HSSM KOPT/OILENS/TSGPLUME, D1.5c).

Sources (printed page [PDF page]):

| Short name | Document |
|------------|----------|
| HSSM1 | Weaver, J.W., Charbeneau, R.J., Tauxe, J.D., Lien, B.K., Provost, J.B. (1994), The Hydrocarbon Spill Screening Model (HSSM) Volume 1: User's Guide, EPA/600/R-94/039a (`hssm_vol1_users_guide.pdf`) |
| HSSM2 | Charbeneau, R.J., Weaver, J.W., Lien, B.K. (1994), HSSM Volume 2: Theoretical Background and Source Codes, EPA/600/R-94/039b (`hssm_vol2_theory_source.pdf`) |
| BSCR | Newell, C.J., McLeod, R.K., Gonzales, J.R. (1996), BIOSCREEN Natural Attenuation Decision Support System, User's Manual Version 1.3, EPA/600/R-96/087, with the Version 1.4 revisions (July 1997) (`bioscreen_users_manual_1996.pdf`) |
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

## 4. Verification

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
