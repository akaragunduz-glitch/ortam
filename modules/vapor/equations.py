"""Vapor intrusion screening equations (VISL).

Sources:
  VISL UG = docs/references/visl_users_guide_2014.pdf (equations on PDF pp. 6-9)
  OSWER   = docs/references/oswer_vi_technical_guide_2015.pdf (Table 6-1, PDF p. 132)
  RSL UG  = docs/references/rsl_users_guide_2026-03.pdf (4.9.9, PDF p. 228)
Pure functions (Fortran: PURE FUNCTION).
"""

from __future__ import annotations

import math

R_L_ATM = 0.082057         # L-atm/mol-K (VISL UG p. 7)
R_M3_ATM = 8.2057e-5       # atm-m3/mol-K (RSL UG 4.9.9)
RC_CAL = 1.9872            # cal/mol-K (RSL UG 4.9.9)
T_REF_K = 298.15           # K


def subslab_target(c_ia: float, af_ss: float) -> float:
    """C_soil-gas = C_target,ia / AF_ss, ug/m3. VISL UG p. 8."""
    return c_ia / af_ss


def groundwater_target(c_ia: float, h: float, af_gw: float) -> float:
    """C_gw = C_target,ia / (H' * AF_gw * 1000 L/m3), ug/L. VISL UG p. 9."""
    return c_ia / (h * af_gw * 1000.0)


def henry_dimensionless_25c(hlc_atm_m3_mol: float) -> float:
    """H'25 = Hc25 * 1000 L/m3 / (R * 298.15 K). VISL UG p. 9."""
    return hlc_atm_m3_mol * 1000.0 / (R_L_ATM * T_REF_K)


def henry_at_temperature(hlc_atm_m3_mol: float, t_c: float, tb_k: float, tc_k: float, dhvb_cal_mol: float) -> float:
    """H' at groundwater temperature, dimensionless. RSL UG 4.9.9, PDF p. 228.

    dHv,gw = dHv,b * [(1 - Tgw/Tc)/(1 - Tb/Tc)]^n, n from Tb/Tc;
    H'(Tgw) = exp[-dHv,gw/Rc * (1/Tgw - 1/298.15)] * HLC / (R * Tgw)
    (The VISL UG 2014 p. 9 print omits the division by R*Tgw; the RSL UG form is used.)
    """
    tgw = t_c + 273.15
    ratio = tb_k / tc_k
    if ratio < 0.57:
        n = 0.3
    elif ratio > 0.71:
        n = 0.41
    else:
        n = 0.74 * ratio - 0.116
    dhv = dhvb_cal_mol * ((1.0 - tgw / tc_k) / (1.0 - tb_k / tc_k)) ** n
    return math.exp(-dhv / RC_CAL * (1.0 / tgw - 1.0 / T_REF_K)) * hlc_atm_m3_mol / (R_M3_ATM * tgw)


def pure_vapor_concentration(vp_mmhg: float, mw_g_mol: float) -> float:
    """Maximum (pure component) vapor concentration at 25 C, ug/m3.

    From the ppbv conversion of VISL UG p. 7, C[ppbv] = C[ug/m3]*1e9*1e-3*R*T/(MW*1e6),
    with the pure-component partial pressure VP/760 atm = 1e9 ppbv * (VP/760):
    C = (VP/760)*1e9 * MW / (R*T).
    """
    return vp_mmhg / 760.0 * 1.0e9 * mw_g_mol / (R_L_ATM * T_REF_K)


def vapor_at_solubility(s_mg_l: float, h: float) -> float:
    """Vapor concentration in equilibrium with water at the solubility limit, ug/m3.

    Inverse of the groundwater equation (VISL UG p. 9) with AF = 1: C_vap = C_w(ug/L) * H' * 1000 L/m3.
    """
    return s_mg_l * 1000.0 * h * 1000.0
