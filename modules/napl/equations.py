"""NAPL/water equilibrium: Raoult's law partitioning, effective solubility, soil saturation.

Every function is a pure function (Fortran: PURE FUNCTION). Sources, printed page [PDF page]:
  HSSM1 = docs/references/hssm_vol1_users_guide.pdf (Weaver, Charbeneau, Tauxe, Lien, Provost, 1994)
  HSSM2 = docs/references/hssm_vol2_theory_source.pdf (Charbeneau, Weaver, Lien, 1995)
  BSCR  = docs/references/bioscreen_users_manual_1996.pdf (BIOSCREEN 1.4 revisions, July 1997)
  RSL UG = docs/references/rsl_users_guide_2026-03.pdf; SSG = Soil Screening Guidance (2002)

  HSSM1 Eq. 43, p. 172 [187]   co = Ko * cw                    (linear NAPL/water partitioning)
  HSSM1 Eq. 44, p. 172 [187]   Ko = Tk * sum_j(coj/Tj) / (sk * gamma_k)   (composition basis)
  HSSM1 Eq. 45, p. 173 [189]   Ko = 1e6 * (rho_o/T_o) / (Sk/Tk)           (average molecular weight basis)
  HSSM2 RAOULT source, pp. 256-257 [269-270]: rho_o = sum_j coj; T_o = sum_j(Tj*coj)/sum_j coj
       (mass-weighted mean); unit factors 1000 mg/g and 1000 cm3/L.
  BSCR p. 3 [68]               effective solubility = mole fraction x pure-phase solubility
  HSSM2 Eq. 25-26, p. 16 [30]  phi*B = theta_w + theta_o*ko + rho_b*kd    (bulk partition coefficient)
  SSG Eq. 4-9; RSL UG 5.12, p. [242]   Csat = S/rho_b * (Kd*rho_b + theta_w + H'*theta_a)
  RSL UG 4.9.4, p. [205]       n = 1 - rho_b/rho_s; theta_a = n - theta_w; Kd = Koc*foc (organics)
"""

from __future__ import annotations

from scipy.optimize import brentq

MG_PER_G = 1000.0     # HSSM2 RAOULT source p. 256 [269], ucmgpg
CM3_PER_L = 1000.0    # HSSM2 RAOULT source p. 256 [269], ucccpl


# ------------------------------------------------------------ NAPL composition
def napl_density(co_g_cm3: list[float]) -> float:
    """Hydrocarbon density = sum_j coj, g/cm3 (HSSM2 RAOULT source, p. 256 [269])."""
    return sum(co_g_cm3)


def average_mw(co_g_cm3: list[float], mw_g_mol: list[float]) -> float:
    """Average molecular weight = sum_j(Tj*coj)/sum_j(coj), g/mol (HSSM2 RAOULT source, p. 256 [269])."""
    return sum(c * t for c, t in zip(co_g_cm3, mw_g_mol)) / sum(co_g_cm3)


def molar_concentration(co_g_cm3: list[float], mw_g_mol: list[float]) -> float:
    """sum_j(coj/Tj), mol/cm3 of NAPL: the sum in HSSM1 Eq. 44."""
    return sum(c / t for c, t in zip(co_g_cm3, mw_g_mol))


def mole_fractions(co_g_cm3: list[float], mw_g_mol: list[float]) -> list[float]:
    """x_k = (cok/Tk) / sum_j(coj/Tj), -. Mole fraction as used in BSCR p. 3 [68]."""
    n = molar_concentration(co_g_cm3, mw_g_mol)
    return [c / t / n for c, t in zip(co_g_cm3, mw_g_mol)]


def ko_composition(mw_k: float, s_k_mg_l: float, gamma_k: float,
                   co_g_cm3: list[float], mw_g_mol: list[float]) -> float:
    """NAPL/water partition coefficient, composition basis, -. HSSM1 Eq. 44, p. 172 [187].

    Ko = Tk * sum_j(coj/Tj) / (sk*gamma_k); with coj in g/cm3 and sk in mg/L the RAOULT
    program multiplies by 1000 mg/g * 1000 cm3/L (HSSM2 p. 256 [269]).
    """
    return mw_k * molar_concentration(co_g_cm3, mw_g_mol) / s_k_mg_l / gamma_k * MG_PER_G * CM3_PER_L


def ko_average_mw(rho_o_g_cm3: float, mw_o_g_mol: float, s_k_mg_l: float, mw_k: float) -> float:
    """NAPL/water partition coefficient, average molecular weight basis, -. HSSM1 Eq. 45, p. 173 [189]."""
    return 1.0e6 * (rho_o_g_cm3 / mw_o_g_mol) / (s_k_mg_l / mw_k)


def effective_solubility(x_k: float, s_k_mg_l: float, gamma_k: float = 1.0) -> float:
    """Effective solubility = mole fraction x pure-phase solubility, mg/L. BSCR p. 3 [68].

    gamma_k is the activity coefficient of HSSM1 Eq. 44 (1 for an ideal solution); with it
    co/Ko from Eq. 43-44 equals x_k * gamma_k * S_k.
    """
    return x_k * gamma_k * s_k_mg_l


# ------------------------------------------------------------ soil
def porosities(rho_b: float, rho_s: float, theta_w: float) -> tuple[float, float]:
    """n = 1 - rho_b/rho_s; theta_a = n - theta_w. RSL UG 4.9.4, p. [205]."""
    n = 1.0 - rho_b / rho_s
    return n, n - theta_w


def bulk_partition(theta_w: float, rho_b: float, kd_l_kg: float, theta_a: float, h: float) -> float:
    """theta_w + rho_b*Kd + theta_a*H', L water-equivalent per L soil.

    The denominator of SSG Eq. 4-9 (RSL UG 5.12); equals phi*B of HSSM2 Eq. 25 (p. 16 [30])
    without the NAPL term and with the air term that HSSM neglects.
    """
    return theta_w + rho_b * kd_l_kg + theta_a * h


def csat(s_mg_l: float, rho_b: float, bp: float) -> float:
    """Csat = S/rho_b * (Kd*rho_b + theta_w + H'*theta_a), mg/kg. SSG Eq. 4-9; RSL UG 5.12, p. [242]."""
    return s_mg_l / rho_b * bp


def solve_napl_moles(n_k: list[float], a_k: list[float]) -> tuple[float, list[float]]:
    """Moles of NAPL per litre of soil and its mole fractions for given component totals.

    Mass balance of component k per litre of soil (derivation in docs/methods/napl.md):
        N_k = x_k * a_k + x_k * N_o,   a_k = gamma_k*S_k*B'_k/M_k
    where N_k = Ct_k*rho_b/M_k (mmol/L), B'_k from bulk_partition, N_o = mmol NAPL per L soil.
    It combines co = Ko*cw (HSSM1 Eq. 43), Ko from Eq. 44 and the bulk coefficient of
    HSSM2 Eq. 25-26. The NAPL exists if sum_k N_k/a_k > 1 (= sum_k Ct_k/Csat_k > 1); then
    N_o solves sum_k N_k/(a_k + N_o) = 1 (monotone, one root in (0, sum N_k]).
    Returns (N_o, x_k); N_o = 0 and x_k = 0 when no NAPL forms.
    """
    insoluble = any(a == 0.0 and n > 0.0 for n, a in zip(n_k, a_k))
    if not insoluble and sum(n / a for n, a in zip(n_k, a_k) if n > 0.0) <= 1.0:
        return 0.0, [0.0] * len(n_k)

    def f(no: float) -> float:
        return sum(n / (a + no) for n, a in zip(n_k, a_k) if n > 0.0) - 1.0

    total = sum(n_k)
    # Lower bracket: the insoluble components alone have x summing to 1 at N_o = their moles.
    lo = sum(n for n, a in zip(n_k, a_k) if a == 0.0)
    if f(total) == 0.0:
        no = total
    else:
        no = brentq(f, lo, total, xtol=total * 1e-15, rtol=1e-14, maxiter=500)
    return no, [n / (a + no) if n > 0.0 else 0.0 for n, a in zip(n_k, a_k)]
