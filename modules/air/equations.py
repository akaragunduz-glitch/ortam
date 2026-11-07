"""Soil-to-air emission factors: volatilization (VF), particulate emission (PEF), Csat.

Every function is a pure function (Fortran: PURE FUNCTION). Sources:
  RSL UG = docs/references/rsl_users_guide_2026-03.pdf (section, PDF page)
  SSG    = docs/references/ssg_supplemental_2002.pdf (EPA OSWER 9355.4-24, Dec 2002;
           equation number; the appendices are separate EPA documents and are not
           included in this file)
"""

from __future__ import annotations

import math

ACRE_M2 = 4047.0          # RSL UG 4.9.3, p. 203 (4047 m2/acre)
SECONDS_PER_YEAR = 3.15e7  # RSL UG 4.9.5, p. 206


def q_over_c(a: float, b: float, c: float, area_acres: float) -> float:
    """Dispersion factor Q/C = A*exp((ln(As) - B)^2 / C), g/m2-s per kg/m3.

    RSL UG 4.9.1 (p. 196), 4.9.4 (p. 204), 4.9.6 (p. 213); SSG Eq. 5-6, 5-15.
    """
    return a * math.exp((math.log(area_acres) - b) ** 2 / c)


def f_of_x(x: float) -> float:
    """F(x) for wind erosion, x = 0.886*Ut/Um. RSL UG 4.9.3/4.9.4, PDF p. 204."""
    if x < 2.0:
        return 1.91207 - 0.0278085 * x + 0.48113 * x ** 2 - 1.09871 * x ** 3 + 0.335341 * x ** 4
    return 0.18 * (8.0 * x ** 3 + 12.0 * x) * math.exp(-x ** 2)


def dispersion_correction(tc_h: float) -> float:
    """FD = 0.1852 + 5.3537/tc + (-9.6318)/tc^2, tc in h. RSL UG 4.9.2 (p. 200), 4.9.6 (p. 214)."""
    return 0.1852 + 5.3537 / tc_h + (-9.6318) / tc_h ** 2


# ------------------------------------------------------------ soil properties
def porosities(rho_b: float, rho_s: float, theta_w: float) -> tuple[float, float]:
    """n = 1 - rho_b/rho_s; theta_a = n - theta_w. RSL UG 4.9.4, p. 205."""
    n = 1.0 - rho_b / rho_s
    return n, n - theta_w


def apparent_diffusivity(h: float, dia: float, diw: float, kd: float,
                         rho_b: float, rho_s: float, theta_w: float) -> float:
    """DA, cm2/s. RSL UG 4.9.4, p. 204; SSG Eq. 4-8.

    DA = [(theta_a^(10/3)*Dia*H' + theta_w^(10/3)*Diw) / n^2] / (rho_b*Kd + theta_w + theta_a*H')
    """
    n, ta = porosities(rho_b, rho_s, theta_w)
    return ((ta ** (10.0 / 3.0) * dia * h + theta_w ** (10.0 / 3.0) * diw) / n ** 2) / (
        rho_b * kd + theta_w + ta * h)


def csat(s: float, kd: float, h: float, rho_b: float, rho_s: float, theta_w: float) -> float:
    """Csat = S/rho_b * (Kd*rho_b + theta_w + H'*theta_a), mg/kg. SSG Eq. 4-9; RSL UG 5.12."""
    _, ta = porosities(rho_b, rho_s, theta_w)
    return s / rho_b * (kd * rho_b + theta_w + h * ta)


# -------------------------------------------------------- volatilization (VF)
def vf_unlimited(qc: float, da: float, t_s: float, rho_b: float) -> float:
    """VFulim = Q/C * (3.14*DA*T)^0.5 * 1e-4 / (2*rho_b*DA), m3/kg. RSL UG 4.9.4, p. 204; SSG Eq. 4-8."""
    return qc * (3.14 * da * t_s) ** 0.5 * 1e-4 / (2.0 * rho_b * da)


def vf_mass_limit(qc: float, t_yr: float, rho_b: float, ds_m: float) -> float:
    """VFmlim = Q/C * T(yr)*3.15e7 / (rho_b(Mg/m3)*ds*1e6), m3/kg. RSL UG 4.9.5, p. 206; SSG Eq. 4-13."""
    return qc * t_yr * SECONDS_PER_YEAR / (rho_b * ds_m * 1.0e6)


def vf_subchronic_unlimited(qc: float, fd: float, da: float, tc_s: float, rho_b: float) -> float:
    """VFulim-sc = Q/Csa * (1/FD) * (3.14*DA*Tc)^0.5 * 1e-4 / (2*rho_b*DA). RSL UG 4.9.6, p. 213; SSG Eq. 5-14."""
    return qc / fd * (3.14 * da * tc_s) ** 0.5 * 1e-4 / (2.0 * rho_b * da)


def vf_subchronic_mass_limit(qc: float, fd: float, tc_s: float, rho_b: float, ds_m: float) -> float:
    """VFmlim-sc = Q/Csa * (1/FD) * Tc / (rho_b*ds*1e6). RSL UG 4.9.7, p. 215; SSG Eq. 5-17."""
    return qc / fd * tc_s / (rho_b * ds_m * 1.0e6)


# --------------------------------------------------- particulate emission (PEF)
def pef_wind(qc: float, v: float, um: float, ut: float, fx: float) -> float:
    """PEF = Q/Cwind * 3600 / (0.036*(1-V)*(Um/Ut)^3*F(x)), m3/kg. RSL UG 4.9.1, p. 196; SSG Eq. 4-5."""
    return qc * 3600.0 / (0.036 * (1.0 - v) * (um / ut) ** 3 * fx)


def pef_traffic_ssg(qc: float, fd: float, t_s: float, ar_m2: float, w_t: float, p_d: float,
                    vkt_km: float) -> float:
    """PEFsc = Q/Csr * (1/FD) * T*AR / (556*(W/3)^0.4*((365-p)/365)*sum VKT). SSG Eq. 5-5, p. 5-11."""
    return qc / fd * t_s * ar_m2 / (556.0 * (w_t / 3.0) ** 0.4 * (365.0 - p_d) / 365.0 * vkt_km)


def pef_traffic_rsl(qc: float, fd: float, t_s: float, ar_m2: float, w_t: float, p_d: float, vkt_km: float,
                    silt_pct: float, m_dry_pct: float) -> float:
    """PEFsc, RSL UG 4.9.2 form (PDF p. 199):
    Q/Csr * (1/FD) * T*AR / ([2.6*(s/12)^0.8*(W/3)^0.4 / (Mdry/0.2)^0.3] * ((365-p)/365) * 281.9 * sum VKT)
    """
    ef = 2.6 * (silt_pct / 12.0) ** 0.8 * (w_t / 3.0) ** 0.4 / (m_dry_pct / 0.2) ** 0.3
    return qc / fd * t_s * ar_m2 / (ef * (365.0 - p_d) / 365.0 * 281.9 * vkt_km)


def m_wind(v: float, um: float, ut: float, fx: float, asurf_m2: float, ed_yr: float) -> float:
    """M_wind = 0.036*(1-V)*(Um/Ut)^3*F(x)*Asurf*ED*8760, g. RSL UG 4.9.3, p. 202."""
    return 0.036 * (1.0 - v) * (um / ut) ** 3 * fx * asurf_m2 * ed_yr * 8760.0


def m_excavation(um: float, moisture_pct: float, rho_soil_mg_m3: float, area_m2: float, depth_m: float,
                 n_dump: float) -> float:
    """M_excav = 0.35*0.0016*(Um/2.2)^1.3/(M/2)^1.4 * rho_soil*A*d*NA*1000, g. RSL UG 4.9.3, p. 202."""
    return (0.35 * 0.0016 * (um / 2.2) ** 1.3 / (moisture_pct / 2.0) ** 1.4
            * rho_soil_mg_m3 * area_m2 * depth_m * n_dump * 1000.0)


def vkt_blade(area_acres: float, blade_m: float, passes: float) -> float:
    """sum VKT = A(acre)*4047/B(m)/1000*NA, km. RSL UG 4.9.3, p. 202."""
    return area_acres * ACRE_M2 / blade_m / 1000.0 * passes


def m_dozing(silt_pct: float, moisture_pct: float, vkt_km: float, speed_kph: float) -> float:
    """M_doz = 0.75*0.45*s^1.5/M^1.4 * VKT/S * 1000, g. RSL UG 4.9.3, p. 202."""
    return 0.75 * 0.45 * silt_pct ** 1.5 / moisture_pct ** 1.4 * vkt_km / speed_kph * 1000.0


def m_grading(speed_kph: float, vkt_km: float) -> float:
    """M_grade = 0.60*0.0056*S^2*VKT*1000, g. RSL UG 4.9.3, p. 202."""
    return 0.60 * 0.0056 * speed_kph ** 2 * vkt_km * 1000.0


def m_tilling(silt_pct: float, area_acres: float, passes: float) -> float:
    """M_till = 1.1*s^0.6*A(acre)*4047*1e-4*1000*NA, g. RSL UG 4.9.3, p. 202."""
    return 1.1 * silt_pct ** 0.6 * area_acres * ACRE_M2 * 1.0e-4 * 1000.0 * passes


def pef_other(qc: float, fd: float, jt_g_m2_s: float) -> float:
    """PEF'sc = Q/Csa * (1/FD) * (1/J'T), m3/kg. RSL UG 4.9.3, p. 202."""
    return qc / fd / jt_g_m2_s
