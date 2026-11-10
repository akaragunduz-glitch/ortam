"""HHRAP 2005 Appendix B, Tables B-4-1 to B-4-25: watershed soil, loads to the water body and
water body concentrations (U.S. EPA 2005, EPA530-R-05-006, hhrap_2005_appendix_b.pdf).

Units are those of the tables (printed page [PDF page] in each docstring). Plain functions of
floats, one per printed equation, like a Fortran module of pure functions.
"""
from __future__ import annotations

import math

S_PER_YR = 3.1536e7       # Tables B-4-6, B-4-20, B-4-21 (s/yr)


# --- Watershed soil (Tables B-4-1 to B-4-6) ------------------------------------------------------
def deposition_term(q_g_s, zs_cm, bd_g_cm3, fv, dytwv, dytwp):
    """Ds = 100 Q / (Zs BD) [Fv Dytwv + Dytwp (1 - Fv)], mg/kg-yr. Table B-4-1, p. B-169 [172]."""
    return 100.0 * q_g_s / (zs_cm * bd_g_cm3) * (fv * dytwv + dytwp * (1.0 - fv))


def soil_conc_at_td(ds, ks, td):
    """CstD = Ds [1 - exp(-ks tD)] / ks, mg/kg (noncarcinogens). Table B-4-1, p. B-169 [172]."""
    return ds * (1.0 - math.exp(-ks * td)) / ks


def soil_conc_average(ds, ks, td, t1, t2):
    """Cs averaged over the exposure duration (carcinogens), mg/kg. Table B-4-1, p. B-168 [171]:
    T1 < tD < T2: Cs = [(Ds tD - CstD)/ks + (CstD/ks)(1 - exp(-ks (T2 - tD)))] / (T2 - T1)
    T2 <= tD:     Cs = Ds / (ks (tD - T1)) ([tD + exp(-ks tD)/ks] - [T1 + exp(-ks T1)/ks])"""
    if t2 <= td:
        return ds / (ks * (td - t1)) * ((td + math.exp(-ks * td) / ks) - (t1 + math.exp(-ks * t1) / ks))
    cstd = soil_conc_at_td(ds, ks, td)
    return ((ds * td - cstd) / ks + cstd / ks * (1.0 - math.exp(-ks * (t2 - td)))) / (t2 - t1)


def k_erosion(xe, sd, er, bd, zs, kds, tsw):
    """kse = 0.1 Xe SD ER / (BD Zs) * Kds BD / (theta_sw + Kds BD), 1/yr. Table B-4-3, p. B-180 [183]
    (recommended default: 0)."""
    return 0.1 * xe * sd * er / (bd * zs) * (kds * bd / (tsw + kds * bd))


def k_runoff(ro, tsw, zs, kds, bd):
    """ksr = RO / (theta_sw Zs) * 1 / (1 + Kds BD / theta_sw), 1/yr. Table B-4-4, p. B-185 [188]."""
    return ro / (tsw * zs) * (1.0 / (1.0 + kds * bd / tsw))


def k_leaching(p, i, ro, ev, tsw, zs, bd, kds):
    """ksl = (P + I - RO - Ev) / (theta_sw Zs [1.0 + BD Kds / theta_sw]), 1/yr. Table B-4-5, p. B-189 [192]."""
    return (p + i - ro - ev) / (tsw * zs * (1.0 + bd * kds / tsw))


def k_volatilization(h, zs, kds, r, ta, bd, da, rho_soil, tsw):
    """ksv = [3.1536e7 H / (Zs Kds R Ta BD)] (Da/Zs) [1 - BD/rho_soil - theta_sw], 1/yr.
    Table B-4-6, p. B-194 [197]."""
    return (S_PER_YR * h / (zs * kds * r * ta * bd)) * (da / zs) * (1.0 - bd / rho_soil - tsw)


def k_soil(ksg, kse, ksr, ksl, ksv):
    """ks = ksg + kse + ksr + ksl + ksv, 1/yr. Table B-4-2, p. B-176 [179]."""
    return ksg + kse + ksr + ksl + ksv


# --- Erosion (Tables B-4-13, B-4-14) ---------------------------------------------------------------
def unit_soil_loss(rf, k, ls, c, pf):
    """Xe = RF K LS C PF 907.18 / 4047, kg/m2-yr (USLE). Table B-4-13, p. B-219 [222]."""
    return rf * k * ls * c * pf * 907.18 / 4047.0


def sd_intercept(al_m2, table, m2_per_sq_mi):
    """Empirical intercept a from the watershed area (sq. miles). Table B-4-14, p. B-223 [226]."""
    area = al_m2 / m2_per_sq_mi
    for upper, a in table:
        if area <= upper:
            return a
    return table[-1][1]


def sediment_delivery_ratio(a, al_m2, b):
    """SD = a AL^(-b). Table B-4-14, p. B-223 [226]."""
    return a * al_m2 ** (-b)


# --- Loads (Tables B-4-7 to B-4-12) ----------------------------------------------------------------
def load_deposition(q, fv, dytwv, dytwp, aw):
    """LDEP = Q [Fv Dytwv + (1 - Fv) Dytwp] Aw, g/yr. Table B-4-8, p. B-201 [204]."""
    return q * (fv * dytwv + (1.0 - fv) * dytwp) * aw


def load_impervious(q, fv, dytwv, dytwp, ai):
    """LRI = Q [Fv Dytwv + (1.0 - Fv) Dytwp] AI, g/yr. Table B-4-9, p. B-204 [207]."""
    return q * (fv * dytwv + (1.0 - fv) * dytwp) * ai


def load_runoff(ro, al, ai, cs, bd, tsw, kds):
    """LR = RO (AL - AI) Cs BD / (theta_sw + Kds BD) 0.01, g/yr. Table B-4-10, p. B-207 [210]."""
    return ro * (al - ai) * cs * bd / (tsw + kds * bd) * 0.01


def load_erosion(xe, al, ai, sd, er, cs, kds, bd, tsw):
    """LE = Xe (AL - AI) SD ER Cs Kds BD / (theta_sw + Kds BD) 0.001, g/yr. Table B-4-11, p. B-211 [214]."""
    return xe * (al - ai) * sd * er * cs * kds * bd / (tsw + kds * bd) * 0.001


def load_diffusion(kv_m_yr, q, fv, cywv, aw, h, r, twk):
    """Ldif = Kv Q Fv Cywv Aw 1e-6 / (H / (R Twk)), g/yr. Table B-4-12, p. B-215 [218]."""
    return kv_m_yr * q * fv * cywv * aw * 1e-6 / (h / (r * twk))


# --- Water body (Tables B-4-15 to B-4-25) ----------------------------------------------------------
def fraction_water_column(kdsw, tss, dwc, dbs, tbs, kdbs, cbs):
    """fwc = (1 + Kdsw TSS 1e-6) dwc/dz / [(1 + Kdsw TSS 1e-6) dwc/dz + (theta_bs + Kdbs CBS) dbs/dz];
    fbs = 1 - fwc; dz = dwc + dbs. Table B-4-16, p. B-230 [233], dz p. B-232 [235]."""
    dz = dwc + dbs
    w = (1.0 + kdsw * tss * 1e-6) * dwc / dz
    return w / (w + (tbs + kdbs * cbs) * dbs / dz)


def kl_river(dw, u, dz):
    """KL = sqrt(1e-4 Dw u / dz) 3.1536e7, m/yr (flowing streams or rivers). Table B-4-20, p. B-244 [247]."""
    return math.sqrt(1e-4 * dw * u / dz) * S_PER_YR


def kl_lake(cd, w, rho_a, rho_w, k, lz, mu_w, dw):
    """KL = (Cd^0.5 W)(rho_a/rho_w)^0.5 (k^0.33/lambda_z)(mu_w/(rho_w Dw))^-0.67 3.1536e7, m/yr
    (quiescent lakes or ponds). Table B-4-20, p. B-244 [247]."""
    return (cd ** 0.5 * w) * (rho_a / rho_w) ** 0.5 * (k ** 0.33 / lz) * (mu_w / (rho_w * dw)) ** -0.67 * S_PER_YR


def kg_lake(cd, w, k, lz, mu_a, rho_a, da):
    """KG = (Cd^0.5 W)(k^0.33/lambda_z)(mu_a/(rho_a Da))^-0.67 3.1536e7, m/yr (quiescent lakes or ponds);
    flowing streams or rivers: KG = 36500 m/yr. Table B-4-21, p. B-249 [252]."""
    return (cd ** 0.5 * w) * (k ** 0.33 / lz) * (mu_a / (rho_a * da)) ** -0.67 * S_PER_YR


def overall_transfer(kl, kg, h, r, twk, theta):
    """Kv = [KL^-1 + (KG H/(R Twk))^-1]^-1 theta^(Twk - 293), m/yr. Table B-4-19, p. B-241 [244]."""
    return (1.0 / kl + 1.0 / (kg * h / (r * twk))) ** -1 * theta ** (twk - 293.0)


def k_volatilization_water(kv_overall, dz, kdsw, tss):
    """kv = Kv / (dz (1 + Kdsw TSS 1e-6)), 1/yr. Table B-4-18, p. B-237 [240]."""
    return kv_overall / (dz * (1.0 + kdsw * tss * 1e-6))


def k_burial(xe, al, sd, vfx, tss, aw, cbs, dbs):
    """kb = ((Xe AL SD 1e3 - Vfx TSS) / (Aw TSS)) (TSS 1e-6 / (CBS dbs)), 1/yr. Table B-4-22, p. B-253 [256].
    A negative value is replaced by 0 by the caller (HHRAP Chapter 5, p. 5-81 [82])."""
    return ((xe * al * sd * 1e3 - vfx * tss) / (aw * tss)) * (tss * 1e-6 / (cbs * dbs))


def k_water_total(fwc, kv, kb):
    """kwt = fwc kv + fbs kb, fbs = 1 - fwc, 1/yr. Table B-4-17, p. B-235 [238]."""
    return fwc * kv + (1.0 - fwc) * kb


def total_water_body(lt, vfx, fwc, kwt, aw, dwc, dbs):
    """Cwtot = LT / (Vfx fwc + kwt Aw (dwc + dbs)), g/m3 = mg/L. Table B-4-15, p. B-226 [229]."""
    return lt / (vfx * fwc + kwt * aw * (dwc + dbs))


def total_water_column(fwc, cwtot, dwc, dbs):
    """Cwctot = fwc Cwtot (dwc + dbs)/dwc, mg/L. Table B-4-23, p. B-257 [260]."""
    return fwc * cwtot * (dwc + dbs) / dwc


def dissolved(cwctot, kdsw, tss):
    """Cdw = Cwctot / (1 + Kdsw TSS 1e-6), mg/L. Table B-4-24, p. B-260 [263]."""
    return cwctot / (1.0 + kdsw * tss * 1e-6)


def bed_sediment(fwc, cwtot, kdbs, tbs, cbs, dwc, dbs):
    """Csb = fbs Cwtot (Kdbs / (theta_bs + Kdbs CBS)) (dwc + dbs)/dbs, mg/kg. Table B-4-25, p. B-263 [266]."""
    return (1.0 - fwc) * cwtot * (kdbs / (tbs + kdbs * cbs)) * (dwc + dbs) / dbs
