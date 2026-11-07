"""Soil-to-air volatilization factors and Csat.

References:
  VFulim  RSL UG 4.9.4 (PDF pp. 204-206); SSG 2002 Eq. 4-8 (p. 4-15)
  VFmlim  RSL UG 4.9.5 (p. 206); SSG 2002 Eq. 4-13
  VFsc    RSL UG 4.9.6 (pp. 213-214); SSG 2002 Eq. 5-14, 5-15
  VFmlim-sc RSL UG 4.9.7 (pp. 215-216); SSG 2002 Eq. 5-17
  Csat    RSL UG 5.12; SSG 2002 Eq. 4-9
Volatile soil SLs built from these VFs are verified against the RSL tables in
workflows/soil_volatiles/tests.
"""

import math

import pytest

from modules.air.ssl_vf import SslVfModule

# Benzene (RSL Chemical Parameters, Nov 2024)
H, DIA, DIW, KOC, S = 0.227, 0.0895, 1.03e-5, 145.8, 1790.0


def run(chemicals, **kw):
    m = SslVfModule()
    return m.run(m.validate({"chemicals": chemicals, **kw}))


def benzene_da(foc=0.006, rb=1.5, tw=0.15):
    n = 1 - rb / 2.65
    ta = n - tw
    kd = KOC * foc
    return ((ta ** (10 / 3) * DIA * H + tw ** (10 / 3) * DIW) / n ** 2) / (rb * kd + tw + ta * H)


BENZENE = {"name": "Benzene", "henry_dimensionless": H, "dia_cm2_s": DIA, "diw_cm2_s": DIW,
           "koc_l_kg": KOC, "solubility_mg_l": S, "mp_c": 5.5}


def test_vf_unlimited_by_hand():
    da = benzene_da()
    expected = 68.18 * (3.14 * da * 8.2e8) ** 0.5 * 1e-4 / (2 * 1.5 * da)
    c = run([BENZENE], toxicity_table="none").chemicals[0]
    assert c.da_cm2_s == pytest.approx(da)
    assert c.vf_m3_kg == pytest.approx(expected)


def test_csat_by_hand_and_liquid_flag():
    kd = KOC * 0.006
    expected = S / 1.5 * (kd * 1.5 + 0.15 + H * (1 - 1.5 / 2.65 - 0.15))
    c = run([BENZENE], toxicity_table="none").chemicals[0]
    assert c.csat_mg_kg == pytest.approx(expected)
    assert c.liquid is True   # MP < 20 C (RSL UG 5.12)


def test_vf_mass_limit_by_hand():
    # VFmlim = Q/C * T(yr)*3.15e7 / (rho_b(Mg/m3) * ds(m) * 1e6)
    expected = 68.18 * 26 * 3.15e7 / (1.5 * 2.0 * 1e6)
    c = run([BENZENE], toxicity_table="none", source_model="mass-limit", source_depth_m=2.0).chemicals[0]
    assert c.vf_m3_kg == pytest.approx(expected)


def test_subchronic_constants():
    # RSL UG 4.9.6-4.9.7: Tc = 30,240,000 s, tc = 8,400 h, FD = 0.18584; Q/Csa(0.5 acre) = 14.31
    from modules.air.equations import dispersion_correction, q_over_c
    assert dispersion_correction(8400) == pytest.approx(0.18584, abs=5e-6)
    assert q_over_c(2.4538, 17.566, 189.0426, 0.5) == pytest.approx(14.31, abs=0.005)


def test_vf_subchronic_by_hand():
    da = benzene_da()
    qc = 2.4538 * math.exp((math.log(0.5) - 17.566) ** 2 / 189.0426)
    tc_h = 1 * 50 * 7 * 24
    fd = 0.1852 + 5.3537 / tc_h - 9.6318 / tc_h ** 2
    expected = qc / fd * (3.14 * da * tc_h * 3600) ** 0.5 * 1e-4 / (2 * 1.5 * da)
    c = run([BENZENE], toxicity_table="none", scenario="subchronic").chemicals[0]
    assert c.vf_m3_kg == pytest.approx(expected)


def test_table_lookup_and_inorganic_kd():
    c = run([{"cas": "7439-97-6"}]).chemicals[0]     # mercury (elemental): inorganic, Kd from table
    assert c.kd_l_kg is not None and c.vf_m3_kg > 0


def test_nonvolatile_flagged():
    o = run([{"cas": "7440-50-8"}])                   # copper
    assert o.chemicals[0].vf_m3_kg is None
    assert any("volatile" in w for w in o.warnings)


def test_site_specific_q_over_c_from_constants():
    o = run([BENZENE], toxicity_table="none", dispersion={"a": 14.0, "b": 18.4385, "c": 209.7845,
                                                          "area_acres": 2.0})
    qc = 14.0 * math.exp((math.log(2.0) - 18.4385) ** 2 / 209.7845)
    assert o.q_over_c == pytest.approx(qc)
