"""VISL conversions: target indoor air -> sub-slab/soil gas and groundwater screening levels.

References:
  VISL UG (2014), pp. 7-9: C_sg = C_ia/AF_ss; C_gw = C_ia/(H' * AF_gw * 1000 L/m3);
      ppbv conversion with R = 0.082057 L-atm/mol-K, T = 298.15 K; NVT rules.
  OSWER VI Technical Guide (2015), Table 6-1, PDF p. 132: AF_ss = 0.03, AF_gw = 0.001
      (0.0005 for laterally extensive fine-grained layers).
  RSL UG 4.9.9, PDF p. 228: H' at groundwater temperature.
"""

import math

import pytest

from modules.vapor.visl import ViSlModule

R, T = 0.082057, 298.15


def run(chemicals, **kw):
    m = ViSlModule()
    return m.run(m.validate({"chemicals": chemicals, **kw}))


def test_defaults_are_oswer_2015_attenuation_factors():
    o = run([{"cas": "71-43-2", "target_ia_ug_m3": 0.36}])
    assert o.af_subslab == 0.03 and o.af_groundwater == 0.001


def test_benzene_by_hand():
    # Benzene: H' 0.22690106295994 at 25 C (RSL Nov 2024)
    h = 0.22690106295994
    c = run([{"cas": "71-43-2", "target_ia_ug_m3": 0.36}]).chemicals[0]
    assert c.target_subslab_ug_m3 == pytest.approx(0.36 / 0.03)
    assert c.target_groundwater_ug_l == pytest.approx(0.36 / (h * 0.001 * 1000))


def test_pure_vapor_concentration_by_hand():
    # VISL UG p. 7 ppbv relation inverted for the pure component: C = (VP/760)*1e9 * MW / (R*T), ug/m3
    vp, mw = 94.8, 78.115           # benzene, RSL Nov 2024
    c = run([{"cas": "71-43-2", "target_ia_ug_m3": 0.36}]).chemicals[0]
    assert c.pure_vapor_ug_m3 == pytest.approx(vp / 760 * 1e9 * mw / (R * T), rel=1e-6)


def test_groundwater_temperature_correction_by_hand():
    # RSL UG 4.9.9, p. 228, with illustrative properties (test of the equation, not of data)
    hlc, tb, tc, dhvb, tgw = 0.00555, 353.2, 562.2, 7342.0, 15.0 + 273.15
    ratio = tb / tc
    n = 0.3 if ratio < 0.57 else (0.41 if ratio > 0.71 else 0.74 * ratio - 0.116)
    dhv = dhvb * ((1 - tgw / tc) / (1 - tb / tc)) ** n
    h_t = math.exp(-dhv / 1.9872 * (1 / tgw - 1 / 298.15)) * hlc / (8.2057e-5 * tgw)
    c = run([{"name": "Test", "henry_atm_m3_mol": hlc, "boiling_point_k": tb, "critical_temp_k": tc,
              "enthalpy_vap_cal_mol": dhvb, "target_ia_ug_m3": 1.0, "volatile": True}],
            toxicity_table="none", groundwater_temp_c=15.0).chemicals[0]
    assert c.henry_used == pytest.approx(h_t)
    assert c.target_groundwater_ug_l == pytest.approx(1.0 / (h_t * 0.001 * 1000))


def test_temperature_without_properties_falls_back_to_25c():
    o = run([{"cas": "71-43-2", "target_ia_ug_m3": 0.36}], groundwater_temp_c=12.0)
    assert o.chemicals[0].henry_temp_c == 25.0
    assert any("25" in w for w in o.warnings)


def test_fine_grained_attenuation_factor():
    a = run([{"cas": "71-43-2", "target_ia_ug_m3": 0.36}]).chemicals[0].target_groundwater_ug_l
    b = run([{"cas": "71-43-2", "target_ia_ug_m3": 0.36}], af_groundwater=0.0005).chemicals[0].target_groundwater_ug_l
    assert b == pytest.approx(2 * a)


def test_nvt_groundwater_when_above_solubility():
    # A chemical whose groundwater target exceeds its solubility is 'not sufficiently volatile/toxic'
    c = run([{"name": "Low solubility", "henry_dimensionless": 0.01, "solubility_mg_l": 0.001,
              "vp_mmhg": 10.0, "mw_g_mol": 100.0, "target_ia_ug_m3": 50.0, "volatile": True}],
            toxicity_table="none").chemicals[0]
    assert c.nvt_groundwater is True and c.target_groundwater_ug_l is None


def test_nvt_indoor_air_when_above_pure_vapor():
    c = run([{"name": "Nonvolatile", "henry_dimensionless": 1e-6, "solubility_mg_l": 1.0,
              "vp_mmhg": 1e-12, "mw_g_mol": 300.0, "target_ia_ug_m3": 1.0, "volatile": True}],
            toxicity_table="none").chemicals[0]
    assert c.nvt_indoor_air is True
    assert c.target_subslab_ug_m3 is None and c.target_groundwater_ug_l is None


def test_nonvolatile_gets_air_target_only():
    # VISL UG (2024 note in 'What's New'): sub-slab and groundwater targets only for volatile chemicals
    c = run([{"cas": "7440-38-2", "target_ia_ug_m3": 0.00065}]).chemicals[0]
    assert c.target_subslab_ug_m3 is None and c.target_groundwater_ug_l is None
