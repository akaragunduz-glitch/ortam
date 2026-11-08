"""Tests for vapor.johnson-ettinger against EPA's J&E spreadsheet v6.0 (2017).

Reference values: data/reference/je/v6-2017/reference_cases.csv. The 16 cases were
computed by EPA's own spreadsheet (docs/references/je_spreadsheet_v6_2017.xlsm),
recalculated in LibreOffice by modules/vapor/tests/generate_je_reference_cases.py.
Case 16 (acetaldehyde) is the case saved in the downloaded spreadsheet; its
LibreOffice result equals the value cached by Excel (alpha = 9.93074221054135e-4).
The cases cover the three source types, the five residential foundation types,
1-3 soil strata and temperatures from 10 to 25 C.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from modules.vapor import je_equations as je
from modules.vapor.johnson_ettinger import JohnsonEttingerModule

ROOT = Path(__file__).resolve().parents[3]
CASES = list(csv.DictReader((ROOT / "data/reference/je/v6-2017/reference_cases.csv").open(encoding="utf-8")))
SOURCE = {"Groundwater": "groundwater", "Exterior Soil Gas": "exterior-soil-gas",
          "Sub-slab Soil Gas": "subslab-soil-gas"}
REL = 1e-9


def _num(s: str) -> float | None:
    return None if s in ("", "NA") else float(s)


def _inputs(case: dict) -> dict:
    strata = [{"soil": case[f"soil_{k}"], "thickness_m": float(case[f"h_{k}_m"])}
              for k in "abc" if case[f"soil_{k}"] != "Not Present"]
    return {
        "source": SOURCE[case["source"]],
        "depth_to_source_m": float(case["ls_m"]),
        "temperature_c": float(case["ts_c"]),
        "foundation": case["foundation"],
        "strata": strata,
        "chemicals": [{"name": case["chemical"], "concentration": float(case["c_medium"])}],
    }


@pytest.fixture(scope="module")
def results():
    m = JohnsonEttingerModule()
    return {c["case"]: m.run(m.validate(_inputs(c))) for c in CASES}


def test_sixteen_reference_cases_present():
    assert len(CASES) == 16
    assert all(c["flag"] == "" for c in CASES)


@pytest.mark.parametrize("case", CASES, ids=lambda c: f"case{c['case']}-{c['chemical']}")
def test_case_matches_spreadsheet(case, results):
    out = results[case["case"]]
    r = out.chemicals[0]
    assert r.hs == pytest.approx(float(case["Hs"]), rel=REL)
    assert r.cs_ug_m3 == pytest.approx(float(case["Cs_ug_m3"]), rel=REL)
    assert out.qb_m3_h == pytest.approx(float(case["Qb_m3_h"]), rel=REL)
    for i, k in enumerate("ABC"):
        ref = _num(case[f"Deff{k}_cm2_s"])
        if ref is not None:
            assert r.deff_strata_cm2_s[i] == pytest.approx(ref, rel=REL)
    if case["source"] == "Groundwater":
        assert out.capillary_zone_height_m == pytest.approx(float(case["hcz_m"]), rel=REL)
        assert out.capillary_zone_total_porosity == pytest.approx(float(case["ncz"]), rel=REL)
        assert out.capillary_zone_water_filled_porosity == pytest.approx(float(case["nwcz"]), rel=REL)
        assert r.deff_cz_cm2_s == pytest.approx(float(case["DeffCZ_cm2_s"]), rel=REL)
    assert r.deff_t_cm2_s == pytest.approx(float(case["DeffT_cm2_s"]), rel=REL)
    assert r.a_param == pytest.approx(float(case["A_param"]), rel=REL)
    assert r.b_param == pytest.approx(float(case["B_param"]), rel=REL, abs=1e-15)
    assert r.alpha == pytest.approx(float(case["alpha"]), rel=REL)
    assert r.cia_ug_m3 == pytest.approx(float(case["Cia_ug_m3"]), rel=REL)
    css = _num(case["Css_ug_m3"])
    if css is None:
        assert r.css_ug_m3 is None          # dirt floor: no sub-slab concentration
    else:
        assert r.css_ug_m3 == pytest.approx(css, rel=REL)


def test_capillary_zone_matches_documentation_table_13():
    # J&E v6 Table 13, p. 70: sand n,cz = 0.253, h,cz = 17.05 cm; silty clay 0.424, 192.31 cm
    assert je.capillary_water_filled_porosity(0.375, 0.053, 0.6852) == pytest.approx(0.253, abs=5e-4)
    assert je.capillary_zone_height_m(0.044) * 100 == pytest.approx(17.05, abs=5e-3)
    assert je.capillary_water_filled_porosity(0.481, 0.111, 0.243) == pytest.approx(0.424, abs=5e-4)
    assert je.capillary_zone_height_m(0.0039) * 100 == pytest.approx(192.31, abs=5e-3)


def test_dirt_floor_limit_and_large_peclet_limit():
    # Eq. 16a for dirt floors; for B -> infinity Eq. 16b tends to Eq. 12: A / (1 + A / C)
    assert je.attenuation_factor(0.01, 0.0, 0.003, dirt_floor=True) == pytest.approx(0.01 / 1.01)
    assert je.attenuation_factor(0.01, 800.0, 0.003, dirt_floor=False) == pytest.approx(0.01 / (1 + 0.01 / 0.003))


def test_air_exchange_adjusts_qsoil_qb():
    # J&E v6 Eq. 17, p. 29: [Qsoil/Qb]adj = 0.003 * 0.45 / ach
    m = JohnsonEttingerModule()
    inp = _inputs(CASES[0])
    inp["building"] = {"ach_per_h": 0.9}
    out = m.run(m.validate(inp))
    assert out.qsoil_qb == pytest.approx(0.0015)
    assert any("Eq. 17" in w for w in out.warnings)


def test_geometry_errors_are_rejected():
    m = JohnsonEttingerModule()
    bad = _inputs(CASES[0])
    bad["strata"][0]["thickness_m"] = 2.5          # strata do not reach the source (Ls = 3 m)
    with pytest.raises(ValueError, match="thickness"):
        m.run(m.validate(bad))
    bad = _inputs(CASES[0])
    bad["depth_to_source_m"] = bad["strata"][0]["thickness_m"] = 0.2   # Ls < Lb + hcz
    with pytest.raises(ValueError, match="capillary"):
        m.run(m.validate(bad))


def test_groundwater_above_solubility_is_flagged():
    m = JohnsonEttingerModule()
    inp = _inputs(CASES[0])
    inp["chemicals"][0]["concentration"] = 2.0e6     # benzene S = 1790 mg/L
    r = m.run(m.validate(inp)).chemicals[0]
    assert r.valid is False and r.cia_ug_m3 is None
    assert r.alpha == pytest.approx(float(CASES[0]["alpha"]), rel=REL)   # alpha does not depend on C


def test_non_detect_cannot_be_entered_as_zero():
    m = JohnsonEttingerModule()
    inp = _inputs(CASES[0])
    inp["chemicals"][0]["concentration"] = 0.0
    with pytest.raises(Exception):
        m.validate(inp)


def test_chemical_lookup_by_cas_and_overrides():
    m = JohnsonEttingerModule()
    inp = _inputs(CASES[0])
    inp["chemicals"] = [{"cas": "71-43-2", "concentration": 100.0},
                        {"cas": "71-43-2", "concentration": 100.0, "henry_atm_m3_mol": 0.0111}]
    out = m.run(m.validate(inp))
    a, b = out.chemicals
    assert a.name == "Benzene" and a.alpha == pytest.approx(float(CASES[0]["alpha"]), rel=REL)
    assert b.hs == pytest.approx(2 * a.hs, rel=1e-12)
