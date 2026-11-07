"""Soil screening levels protective of groundwater (EPA SSL).

Sources: RSL UG = docs/references/rsl_users_guide_2026-03.pdf; SSG = ssg_supplemental_2002.pdf.
Pure functions (Fortran: PURE FUNCTION).
"""

from __future__ import annotations

import math


def ssl_partitioning(cw_mg_l: float, kd_l_kg: float, theta_w: float, theta_a: float, h: float,
                     rho_b_kg_l: float) -> float:
    """Method 1: SSL = Cw*[Kd + (theta_w + theta_a*H')/rho_b], mg/kg. RSL UG 4.8.3, p. 190; SSG Eq. 4-10."""
    return cw_mg_l * (kd_l_kg + (theta_w + theta_a * h) / rho_b_kg_l)


def ssl_mass_limit(cw_mg_l: float, infiltration_m_yr: float, ed_yr: float, rho_b_kg_l: float,
                   ds_m: float) -> float:
    """Method 2: SSL = Cw*I*ED/(rho_b*ds), mg/kg. RSL UG 4.8.4, p. 193; SSG Eq. 4-14."""
    return cw_mg_l * infiltration_m_yr * ed_yr / (rho_b_kg_l * ds_m)


def mixing_zone_depth(length_m: float, infiltration_m_yr: float, k_m_yr: float, gradient: float,
                      aquifer_m: float) -> float:
    """d = (0.0112*L^2)^0.5 + da*[1 - exp(-L*I/(K*i*da))], m. RSL UG 4.8.5, p. 194; SSG Eq. 4-12."""
    return (0.0112 * length_m ** 2) ** 0.5 + aquifer_m * (
        1.0 - math.exp(-length_m * infiltration_m_yr / (k_m_yr * gradient * aquifer_m)))


def dilution_attenuation_factor(k_m_yr: float, gradient: float, d_m: float, infiltration_m_yr: float,
                                length_m: float) -> float:
    """DAF = 1 + K*i*d/(I*L), -. RSL UG 4.8.5, p. 194; SSG Eq. 4-11."""
    return 1.0 + k_m_yr * gradient * d_m / (infiltration_m_yr * length_m)
