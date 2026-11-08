"""Biotransfer and bioconcentration factors for organic chemicals from log Kow (HHRAP 2005 Appendix A-2).

Source: HHRAP 2005 Appendix A-2 (docs/references/hhrap_2005_appendix_a.pdf), printed page
[PDF page]. Used when a chemical is not in the HHRAP database (data/reference/hhrap/2005);
the database values were computed with these same equations (explanation codes in
explanations.csv), which is how the functions are verified (tests/test_btf.py).

  RCF       Eq. A-2-14 (log Kow >= 2), A-2-15 (< 2), p. A-2-16 [45]; log Kow capped to
            [-0.57, 8.2]; fresh weight converted to dry weight with 87 % moisture (p. A-2-16).
  Br_rootveg = RCF / Kds                     Eq. A-2-16, p. A-2-17 [46]
  Br_ag, Br_forage (also silage, grain)      Eq. A-2-17, A-2-18, p. A-2-18 [47]; log Kow capped
            to [1.15, 9.35]
  Bv        Eq. A-2-19, A-2-20, p. A-2-20 [49]; rho_air 1.19 g/L, rho_forage 770 g/L,
            f_water 0.85, T = 298.1 K; reduced by a factor of 100 for organics other than
            PCDDs/PCDFs (database explanation 53, Lorber 1997)
  Ba_fat    Eq. A-2-21, p. A-2-23 [52]; log Kow capped to [-0.67, 8.2]
  Ba = Ba_fat x fat fraction: milk 0.04, beef 0.19 (Eq. A-2-22, A-2-23, p. A-2-23 [52]),
            pork 0.23 (p. A-2-26 [55]), chicken 0.14, eggs 0.08 (p. A-2-27 [56])
Not covered: metals (Baes et al. 1984 / EPA 1992 uptake slopes: database values only),
ionizing organics (Ba with Eq. A-2-24/25 needs pKa), fish BCF (Meylan correction factors).
"""

from __future__ import annotations

import math

R_ATM_M3 = 8.205e-5           # atm-m3/(mol K)
T_BV_K = 298.1                # K (p. A-2-20)
RHO_AIR_G_L = 1.19
RHO_FORAGE_G_L = 770.0
F_WATER_FORAGE = 0.85
BV_REDUCTION = 100.0
ROOT_MOISTURE = 0.87
FAT = {"milk": 0.04, "beef": 0.19, "pork": 0.23, "chicken": 0.14, "eggs": 0.08}


def _cap(x: float, lo: float, hi: float) -> float:
    return min(max(x, lo), hi)


def rcf_dw(log_kow: float) -> float:
    """Root concentration factor, (ug/g DW)/(ug/mL)."""
    lk = _cap(log_kow, -0.57, 8.2)
    fw = 10 ** (0.77 * lk - 1.52) if lk >= 2.0 else 10 ** (0.77 * lk - 1.52) + 0.82
    return fw / (1.0 - ROOT_MOISTURE)


def br_rootveg(rcf: float, kds_ml_g: float) -> float:
    return rcf / kds_ml_g


def br_aboveground(log_kow: float) -> float:
    """Br_ag = Br_forage (also silage and grain): plant DW / soil, unitless."""
    return 10 ** (1.588 - 0.578 * _cap(log_kow, 1.15, 9.35))


def bv(log_kow: float, henry_atm_m3_mol: float) -> float:
    """Air-to-plant biotransfer factor (DW), already divided by 100."""
    log_bvol = 1.065 * log_kow - math.log10(henry_atm_m3_mol / (R_ATM_M3 * T_BV_K)) - 1.654
    return RHO_AIR_G_L * 10 ** log_bvol / ((1.0 - F_WATER_FORAGE) * RHO_FORAGE_G_L) / BV_REDUCTION


def ba(log_kow: float, product: str) -> float:
    """Biotransfer factor, day/kg FW, for milk, beef, pork, chicken or eggs."""
    lk = _cap(log_kow, -0.67, 8.2)
    return 10 ** (-0.099 * lk * lk + 1.07 * lk - 3.56) * FAT[product]
