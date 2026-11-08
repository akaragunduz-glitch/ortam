"""workflows.soil_to_groundwater.vleach_to_groundwater: VLEACH flux -> leachate -> groundwater."""

from pathlib import Path

import pytest

from modules.vadose.vleach_inp import read_inp
from workflows.soil_to_groundwater import vleach_to_groundwater

REF = Path(__file__).resolve().parents[3] / "data" / "reference" / "vleach" / "2.2a"


def test_leachate_is_flux_over_recharge_and_dilution_by_daf():
    inp = read_inp(REF / "case1_boundaries" / "input.inp")
    r = vleach_to_groundwater(inp, daf=20.0, target_mg_l=0.005)
    p, v = r["polygons"][0], r["vleach"].polygons[0]
    q, area = inp["polygons"][0]["recharge_m_yr"], v.area_m2
    t, rate = v.groundwater_series[10]
    assert p["leachate_mg_l"][10] == (t, pytest.approx(rate / (area * q)))
    assert p["groundwater_mg_l"][10][1] == pytest.approx(rate / (area * q) / 20.0)
    assert p["peak_groundwater_mg_l"] == pytest.approx(p["peak_leachate_mg_l"] / 20.0)
    assert p["exceeds_target"] is (p["peak_groundwater_mg_l"] > 0.005)


def test_no_recharge_is_reported():
    inp = read_inp(REF / "case4_diffusion_only" / "input.inp")
    r = vleach_to_groundwater(inp)
    assert r["polygons"] == [] and "no recharge" in r["notes"][0]
