"""Particulate emission factors.

References:
  PEF (wind)  RSL UG 4.9.1 (PDF p. 196); SSG 2002 Eq. 4-5
  PEFsc       SSG 2002 Eq. 5-5, 5-6 (pp. 5-11 to 5-13); RSL UG 4.9.2 (pp. 198-200)
  PEF'sc      RSL UG 4.9.3 (pp. 201-203); SSG 2002 Appendix E, Eq. E-15 to E-26
"""

import math

import pytest

from modules.air.ssl_pef import SslPefModule


def run(**kw):
    m = SslPefModule()
    return m.run(m.validate(kw))


def test_default_pef_is_rsl_value():
    # Q/Cwind 93.77, V 0.5, Um 4.69, Ut 11.32 -> PEF = 1.36E+09 m3/kg (RSL Default Exposure Parameters)
    o = run(method="wind")
    assert float(f"{o.pef_m3_kg:.3g}") == 1.36e9


def test_fx_polynomial_close_to_workbook():
    # RSL UG 4.9.4 polynomial gives 0.1931 at Um 4.69, Ut 11.32; the workbook lists 0.194 (0.5 % apart).
    from modules.air.equations import f_of_x
    assert f_of_x(0.886 * 11.32 / 4.69) == pytest.approx(0.194, rel=0.01)


def test_site_specific_wind_uses_polynomial():
    from modules.air.equations import f_of_x
    o = run(method="wind", wind={"um_m_s": 5.0, "ut_m_s": 11.32})
    assert o.details["f_x"] == pytest.approx(f_of_x(0.886 * 11.32 / 5.0))


def test_traffic_pefsc_ssg_example():
    # SSG 2002 p. 5-12 example: W = 8 t, sum VKT = 175.5 km; road 45 m across a 0.5 acre site
    lr_ft = 45 / 0.3048
    o = run(method="construction-traffic", traffic={"road_length_ft": lr_ft, "mean_vehicle_weight_t": 8.0,
                                                    "rain_days_per_yr": 130, "vkt_km": 175.5})
    qc = 12.9351 * math.exp((math.log(0.5) - 5.7383) ** 2 / 71.7711)
    tc_h = 1 * 50 * 7 * 24
    fd = 0.1852 + 5.3537 / tc_h - 9.6318 / tc_h ** 2
    t_s = 1 * 250 * 8 * 3600
    ar = lr_ft * 20 * 0.092903
    expected = qc / fd * t_s * ar / (556 * (8 / 3) ** 0.4 * (365 - 130) / 365 * 175.5)
    assert qc == pytest.approx(23.02, abs=0.01)
    assert o.pef_m3_kg == pytest.approx(expected)


def test_traffic_rsl_ug_form_needs_silt_and_moisture():
    with pytest.raises(Exception):
        run(method="construction-traffic", traffic={"road_length_ft": 150, "mean_vehicle_weight_t": 8,
                                                    "rain_days_per_yr": 130, "vkt_km": 175.5, "form": "rsl-ug"})


def test_other_activities_by_hand():
    act = {"area_acres": 0.5, "excavation_area_m2": 1000, "excavation_depth_m": 1.0,
           "dozing_area_acres": 0.5, "dozer_blade_m": 2.5, "dozing_passes": 1,
           "grading_area_acres": 0.5, "grader_blade_m": 3.0, "grading_passes": 1, "tilling_area_acres": 0.25}
    o = run(method="construction-other", other=act)
    asurf = 0.5 * 4047
    fx = 0.194  # workbook value (PEF Mechanical Other, row 15)
    m_wind = 0.036 * (1 - 0) * (4.69 / 11.32) ** 3 * fx * asurf * 1 * 8760
    m_excav = 0.35 * 0.0016 * (4.69 / 2.2) ** 1.3 / (12 / 2) ** 1.4 * 1.68 * 1000 * 1.0 * 2 * 1000
    vkt_doz = 0.5 * 4047 / 2.5 / 1000 * 1
    m_doz = 0.75 * 0.45 * 6.9 ** 1.5 / 7.9 ** 1.4 * vkt_doz / 11.4 * 1000
    vkt_grade = 0.5 * 4047 / 3.0 / 1000 * 1
    m_grade = 0.60 * 0.0056 * 11.4 ** 2 * vkt_grade * 1000
    m_till = 1.1 * 18 ** 0.6 * 0.25 * 4047 * 1e-4 * 1000 * 2
    t_s = 1 * 250 * 8 * 3600
    jt = (m_wind + m_excav + m_doz + m_grade + m_till) / (asurf * t_s)
    qc = 2.4538 * math.exp((math.log(0.5) - 17.566) ** 2 / 189.0426)
    fd = 0.1852 + 5.3537 / 8400 - 9.6318 / 8400 ** 2
    assert o.pef_m3_kg == pytest.approx(qc / fd / jt, rel=1e-9)
