"""Recreator soil/sediment and surface water.

The recreator equations (RSL UG 4.6, PDF pp. 124-182) have the same form as
the resident equations with recreator parameters, derived as ED-weighted means
of age-segment values (4.6.1.7, 4.6.2.7). Equivalence tests: with the resident
parameters, the recreator tools must reproduce the verified resident results.
"""

import csv

import pytest

from modules.risk.recreator_soil import RecreatorSoilModule
from modules.risk.recreator_surface_water import RecreatorSurfaceWaterModule
from modules.risk.resident_soil import ResidentSoilModule
from modules.risk.resident_tapwater import ResidentTapwaterModule
from ortam_core.chemicals.rsl_import import DATA

RES_SOIL = {"ef_d_yr": [350] * 4, "et_h_d": [24] * 4, "irs_mg_d": [200, 200, 100, 100],
            "sa_soil_cm2_d": [2373, 2373, 6032, 6032], "af_mg_cm2": [0.2, 0.2, 0.07, 0.07]}
# Resident tap water: IRW 0.78 / 2.5 L/day, ET_event 0.54 / 0.71 h, EV 1 -> IRW (L/h) = IRW_day / ET_event
RES_WATER = {"ef_d_yr": [350] * 4, "ev_per_d": [1] * 4, "et_event_h": [0.54, 0.54, 0.71, 0.71],
             "irw_l_h": [0.78 / 0.54, 0.78 / 0.54, 2.5 / 0.71, 2.5 / 0.71]}


def table():
    with (DATA / "2024-11" / "toxicity.csv").open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def tox(r):
    return any(r[k] for k in ("sfo_per_mg_kg_day", "iur_per_ug_m3", "rfd_mg_kg_day", "rfc_mg_m3"))


def run(module, chemicals, **kw):
    return module.run(module.validate({"chemicals": chemicals, **kw}))


def test_soil_equivalent_to_resident_for_all_nonvolatile_chemicals():
    rec, res = RecreatorSoilModule(), ResidentSoilModule()
    rows = [r for r in table() if r["cas"] and r["volatile"] != "V" and tox(r)]
    assert len(rows) > 550
    for r in rows:
        ch = [{"cas": r["cas"], "name": r["name"]}]
        a = run(rec, ch, exposure=RES_SOIL).chemicals[0]
        b = run(res, ch).chemicals[0]
        for f in ("sl_nc_child_mg_kg", "sl_nc_adult_mg_kg", "sl_ca_mg_kg"):
            for route in ("ingestion", "dermal", "inhalation", "total"):
                x, y = getattr(getattr(a, f), route), getattr(getattr(b, f), route)
                assert (x is None and y is None) or x == pytest.approx(y, rel=1e-9), (r["name"], f, route)


@pytest.mark.parametrize("cas", ["75-01-4", "79-01-6"])
def test_soil_vc_tce_equivalent_to_resident(cas):
    a = run(RecreatorSoilModule(), [{"cas": cas, "vf_m3_kg": 3000.0}], exposure=RES_SOIL).chemicals[0]
    b = run(ResidentSoilModule(), [{"cas": cas, "vf_m3_kg": 3000.0}]).chemicals[0]
    assert a.sl_ca_mg_kg.total == pytest.approx(b.sl_ca_mg_kg.total, rel=1e-9)


def test_water_equivalent_to_resident_tapwater_ingestion_and_dermal():
    rec, res = RecreatorSurfaceWaterModule(), ResidentTapwaterModule()
    rows = [r for r in table() if r["cas"] and tox(r)]
    for r in rows[::5]:   # every fifth chemical keeps the test fast
        ch = [{"cas": r["cas"], "name": r["name"]}]
        a = run(rec, ch, exposure=RES_WATER).chemicals[0]
        b = run(res, ch).chemicals[0]
        for f in ("sl_nc_child_ug_l", "sl_nc_adult_ug_l", "sl_ca_ug_l"):
            for route in ("ingestion", "dermal"):
                x, y = getattr(getattr(a, f), route), getattr(getattr(b, f), route)
                assert (x is None and y is None) or x == pytest.approx(y, rel=1e-9), (r["name"], f, route)
            assert getattr(a, f).inhalation is None   # no inhalation for surface water (4.6.2, p. 153)


def test_water_child_ingestion_by_hand():
    # RSL UG 4.6.2.1 (p. 153): copper RfD 0.04; EF 52 d/yr, EV 1, ET_event 2 h, IRW 0.12 L/h
    ex = {"ef_d_yr": [52] * 4, "ev_per_d": [1] * 4, "et_event_h": [2] * 4}
    expected = 1.0 * 365 * 6 * 15 / ((1 / 0.04) / 1000 * 52 * 6 * 1 * 2 * 0.12)
    c = run(RecreatorSurfaceWaterModule(), [{"cas": "7440-50-8"}], exposure=ex).chemicals[0]
    assert c.sl_nc_child_ug_l.ingestion == pytest.approx(expected)


def test_site_specific_inputs_required():
    with pytest.raises(Exception):
        RecreatorSoilModule().validate({"chemicals": [{"cas": "7440-50-8"}]})
