"""Johnson & Ettinger (1991) vapor intrusion model, EPA implementation v6.0 (2017).

Source: "Documentation for EPA's Implementation of the Johnson and Ettinger Model
to Evaluate Site Specific Vapor Intrusion into Buildings, Version 6.0", U.S. EPA,
September 2017 (docs/references/je_documentation_v6_2017.pdf; cited "J&E v6",
page numbers as printed). Where the documentation leaves a detail open, the
formula of EPA's spreadsheet (docs/references/je_spreadsheet_v6_2017.xlsm,
worksheet MODEL) is followed and the cell is cited, because ORTAM is verified
against that spreadsheet.

Spreadsheet constants (not the textbook ones), kept for exact agreement:
  absolute temperature T = Ts + 273 (MODEL!F27), reference TR = 298 K,
  R = 0.000082057 atm-m3/(mol K), Rc = 1.9872 cal/(mol K).
The VISL module (modules/vapor/equations.py) uses 273.15 and 298.15 as in the
RSL User's Guide; the difference in H' is below 0.1 %.

Units: lengths in m unless stated, diffusion coefficients in cm2/s, flows in m3/h,
concentrations in ug/m3 (soil gas, indoor air) and ug/L (groundwater).
"""

from __future__ import annotations

import math

R_ATM_M3 = 0.000082057   # atm-m3/(mol K), MODEL!F26-F27
RC_CAL = 1.9872          # cal/(mol K), J&E v6 Eq. 2, p. 12
T_REF_K = 298.0          # K, MODEL!F26-F27
T_OFFSET = 273.0         # C -> K as in MODEL!F27
CM2S_TO_M2H = 0.36       # cm2/s -> m2/h (1e-4 * 3600), J&E v6 p. 27


def henry_reference(hc_atm_m3_mol: float) -> float:
    """H'R at TR = 298 K (MODEL!F26): Hc / (R * TR). Dimensionless."""
    return hc_atm_m3_mol / (R_ATM_M3 * T_REF_K)


def enthalpy_at_temperature(dhvb_cal_mol: float, ts_c: float, tb_k: float, tc_k: float) -> float:
    """dHv,S at the system temperature (J&E v6 Eq. 3, p. 13), cal/mol.

    dHv,S = dHv,b * [(1 - TS/TC) / (1 - TB/TC)]^n, with n by TB/TC:
    < 0.57 -> 0.3; 0.57-0.71 -> 0.74 TB/TC - 0.116; > 0.71 -> 0.41 (J&E v6 Table 1, p. 14).
    """
    ratio = tb_k / tc_k
    n = 0.3 if ratio < 0.57 else (0.41 if ratio > 0.71 else 0.74 * ratio - 0.116)
    ts_k = ts_c + T_OFFSET
    return dhvb_cal_mol * ((1.0 - ts_k / tc_k) / (1.0 - ratio)) ** n


def henry_at_temperature(hc_atm_m3_mol: float, ts_c: float, tb_k: float, tc_k: float,
                         dhvb_cal_mol: float) -> float:
    """H'S, dimensionless Henry's law constant at TS (J&E v6 Eq. 2, p. 12; MODEL!F27).

    H'S = exp[-(dHv,S / Rc) (1/TS - 1/TR)] * Hc / (R * TS).
    """
    ts_k = ts_c + T_OFFSET
    dhvs = enthalpy_at_temperature(dhvb_cal_mol, ts_c, tb_k, tc_k)
    return math.exp(-(dhvs / RC_CAL) * (1.0 / ts_k - 1.0 / T_REF_K)) * hc_atm_m3_mol / (R_ATM_M3 * ts_k)


def source_vapor_concentration(c_medium: float, hs: float, groundwater: bool) -> float:
    """CS, vapor concentration at the source, ug/m3 (J&E v6 Eq. 1, p. 12; MODEL!F14).

    Groundwater: CS = H'S * Cgw[ug/L] * 1000 L/m3. Soil gas: CS = the measured value.
    """
    return c_medium * 1000.0 * hs if groundwater else c_medium


def capillary_water_filled_porosity(n: float, theta_r: float, m_vg: float) -> float:
    """theta_w,cz = theta_r + (n - theta_r) / 2^M (J&E v6 Eq. 4b, p. 15; SOIL_DATA!G)."""
    return theta_r + (n - theta_r) / 2.0 ** m_vg


def capillary_zone_height_m(mean_grain_diameter_cm: float) -> float:
    """hcz = 0.15 / (0.2 * D) cm (J&E v6 Eq. 7b, p. 19; SOIL_DATA!H); returned in m."""
    return 0.15 / (0.2 * mean_grain_diameter_cm) / 100.0


def effective_diffusivity(dair: float, dwater: float, n: float, nw: float, hs: float) -> float:
    """Millington-Quirk effective diffusion coefficient, cm2/s (J&E v6 Eq. 5b, p. 17 and Eq. 8, p. 23).

    Deff = Da * (n - nw)^3.33 / n^2 + (Dw / H'S) * nw^3.33 / n^2 (MODEL!G104-G107).
    """
    return (dair * (n - nw) ** 3.33 + dwater * nw ** 3.33 / hs) / n ** 2


def total_effective_diffusivity(deff: list[float], thickness_m: list[float], lb_m: float,
                                deff_cz: float | None = None, hcz_m: float | None = None) -> float:
    """DTeff, cm2/s, for layers in series (J&E v6 Eq. 9, p. 23; MODEL!G108).

    DTeff = LT / sum(Li / Dieff). The path starts at the base of the foundation, so the
    top stratum counts with thickness hA - Lb. For a groundwater source the capillary
    zone (hcz, deff_cz) is the lowest layer and is taken out of the deepest stratum.
    """
    lengths = list(thickness_m)
    lengths[0] -= lb_m
    total = sum(lengths)
    if deff_cz is not None:
        lengths[-1] -= hcz_m
        resistance = sum(L / d for L, d in zip(lengths, deff)) + hcz_m / deff_cz
    else:
        resistance = sum(L / d for L, d in zip(lengths, deff))
    return total / resistance


def ventilation_rate(abf_m2: float, hb_m: float, ach_per_h: float) -> float:
    """Qb = Abf * Hb * ach, m3/h (J&E v6 Eq. 18, p. 29)."""
    return abf_m2 * hb_m * ach_per_h


def qsoil_qb_adjusted(ach_per_h: float, default_ratio: float = 0.003, default_ach: float = 0.45) -> float:
    """[Qsoil/Qb]adj = 0.003 * 0.45 / ach (J&E v6 Eq. 17, p. 29)."""
    return default_ratio * default_ach / ach_per_h


def exposed_area_m2(abf_m2: float, lb_m: float) -> float:
    """Floor plus wall area below grade, Abf + 4 Lb sqrt(Abf), m2 (J&E v6 p. 27)."""
    return abf_m2 + 4.0 * lb_m * math.sqrt(abf_m2)


def a_param(dteff: float, abf_m2: float, lb_m: float, qb_m3_h: float, ls_m: float) -> float:
    """A = DTeff (Abf + 4 Lb sqrt(Abf)) 0.36 / (Qb (Ls - Lb)) (J&E v6 Eq. 13, p. 26; MODEL!G111)."""
    return dteff * exposed_area_m2(abf_m2, lb_m) * CM2S_TO_M2H / (qb_m3_h * (ls_m - lb_m))


def b_param(qsoil_qb: float, qb_m3_h: float, lf_m: float, deff_a: float, eta: float, abf_m2: float,
            lb_m: float) -> float:
    """B = Qsoil/Qb Qb Lf / (DAeff eta (Abf + 4 Lb sqrt(Abf)) 0.36) (J&E v6 Eq. 14, p. 26; MODEL!G112).

    Peclet number of transport through the foundation cracks.
    """
    return qsoil_qb * qb_m3_h * lf_m / (deff_a * eta * exposed_area_m2(abf_m2, lb_m) * CM2S_TO_M2H)


def attenuation_factor(a: float, b: float, c: float, dirt_floor: bool) -> float:
    """Steady-state attenuation coefficient alpha (J&E v6 Eq. 16a and 16b, p. 28; MODEL!G93).

    Dirt floor: alpha = A / (1 + A).
    Slab:       alpha = A / (1 + A e^-B + (A / C)(1 - e^-B)).
    For very large B, e^-B -> 0 and alpha tends to Eq. 12 (p. 26), as in the spreadsheet.
    """
    if dirt_floor:
        return a / (1.0 + a)
    e = math.exp(-b)
    return a / (1.0 + a * e + (a / c) * (1.0 - e))


def saturated_vapor_concentration(vp_mmhg: float, mw_g_mol: float) -> float:
    """Saturated vapor concentration, ug/m3 (CHEM_DATA!G): VP/760 * 1e9 * MW / 24.46."""
    return vp_mmhg / 760.0 * 1.0e9 * mw_g_mol / 24.46
