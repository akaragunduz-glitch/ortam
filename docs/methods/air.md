# Air module: methods

## `air.ssl-vf` — volatilization factor and Csat (task D1.2a)

| Quantity | Equation | Source |
|----------|----------|--------|
| Apparent diffusivity DA | [(θa^(10/3)·Dia·H' + θw^(10/3)·Diw)/n²] / (ρb·Kd + θw + θa·H') | RSL UG §4.9.4, PDF p. 204; SSG 2002 Eq. 4-8 |
| Porosities | n = 1 − ρb/ρs; θa = n − θw | RSL UG §4.9.4, p. 205 |
| Kd | Koc·foc (organics); table Kd (inorganics) | RSL UG §4.9.4, p. 205 |
| VFulim (chronic) | Q/C·(3.14·DA·T)^½·10⁻⁴ / (2·ρb·DA) | RSL UG §4.9.4, p. 204; SSG Eq. 4-8 |
| VFmlim (chronic) | Q/C·T(yr)·3.15·10⁷ / (ρb·ds·10⁶) | RSL UG §4.9.5, p. 206; SSG Eq. 4-13 |
| VFulim-sc (construction) | Q/Csa·(1/FD)·(3.14·DA·Tc)^½·10⁻⁴ / (2·ρb·DA) | RSL UG §4.9.6, p. 213; SSG Eq. 5-14 |
| VFmlim-sc (construction) | Q/Csa·(1/FD)·Tc / (ρb·ds·10⁶) | RSL UG §4.9.7, p. 215; SSG Eq. 5-17 |
| Q/C | A·exp[(ln As − B)²/C] | RSL UG §4.9.1/4.9.4/4.9.6; SSG Eq. 5-6, 5-15 |
| FD | 0.1852 + 5.3537/tc − 9.6318/tc², tc = ED·EW·7·24 h | RSL UG §4.9.6, p. 214 |
| Csat | S/ρb·(Kd·ρb + θw + H'·θa) | SSG Eq. 4-9; RSL UG §5.12 |

Defaults (`data/reference/rsl/2024-11/air_emission_factors.yaml`, from the RSL
parameter workbook): ρb 1.5 g/cm3, ρs 2.65 g/cm3, θw 0.15, foc 0.006,
Q/Cvol 68.18 (0.5 acre), T 8.2·10⁸ s (26 yr); construction A 2.4538,
B 17.566, C 189.0426, ED 1 yr, EW 50 wk/yr.

**Csat.** Liquid = melting point < 20 °C (RSL UG §5.12). Following the RSL
tables, SLs above Csat are flagged, never replaced.

**Verification** (`workflows/soil_volatiles`): VF from `air.ssl-vf` passed to
`risk.resident-soil` and `risk.worker-soil` reproduces the RSL Summary Table
(Nov 2024) for all volatile chemicals: Resident Soil 308/308, Industrial
Soil 308/308 (including VC and TCE). One industrial value (n-nonane,
72.501 vs 72) lies on a rounding boundary and is accepted by a documented
tolerance. Where VF cannot be calculated (H', Dia, Diw or Kd missing; 17
chemicals, e.g. HCl, phosphine), the RSL evaluates inhalation with the PEF
only; the workflow does the same and reports it.

## `air.ssl-pef` — particulate emission factors (task D1.2a)

| Method | Equation | Source |
|--------|----------|--------|
| wind (chronic) | PEF = Q/Cwind·3600 / [0.036·(1−V)·(Um/Ut)³·F(x)] | RSL UG §4.9.1, p. 196; SSG Eq. 4-5 |
| construction-traffic | PEFsc = Q/Csr·(1/FD)·T·AR / [556·(W/3)^0.4·((365−p)/365)·ΣVKT] | SSG Eq. 5-5, p. 5-11 (default form) |
| construction-traffic, `form: rsl-ug` | 556 replaced by 2.6·(s/12)^0.8/(Mdry/0.2)^0.3·281.9 | RSL UG §4.9.2, p. 199 |
| construction-other | PEF'sc = Q/Csa·(1/FD)·(1/⟨J'T⟩), ⟨J'T⟩ = ΣM/(Asurf·Tt); M_wind, M_excav, M_doz, M_grade, M_till | RSL UG §4.9.3, pp. 201-203 |

**F(x).** The polynomial of RSL UG §4.9.4 gives F(x) = 0.1931 at the default
wind speeds (Um 4.69, Ut 11.32 m/s); the parameter workbook lists 0.194.
Only 0.194 reproduces the RSL default PEF of 1.36·10⁹ m3/kg, so the workbook
value is used with the default wind speeds and the polynomial is used when
Um or Ut is site-specific.

**PEFsc forms.** SSG 2002 Eq. 5-5 (constant 556) is the source cited by the
RSL parameter workbook and is the default. The RSL User's Guide prints a
newer AP-42 form with silt content s and moisture Mdry but gives no default
values; it is available as `form: rsl-ug` and requires both.

**Missing source material.** The SSG 2002 file in `docs/references/` ends at
the References; its appendices are separate EPA documents. Appendix D
(Q/C constants A, B, C for 29 cities) is therefore not available: the tools
accept Q/C directly or user-supplied A, B, C. Appendix E:
https://semspub.epa.gov/src/document/hq/174068.pdf
