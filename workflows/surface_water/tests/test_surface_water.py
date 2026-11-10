"""Workflows plume_to_water_body and water_body_fish_risk."""
import copy
from pathlib import Path

import pytest
import yaml

from modules.napl.source_plume import SourcePlumeModule
from modules.surfacewater.watershed import WatershedModule
from workflows.foodchain import foodchain_risk
from workflows.surface_water import plume_to_water_body, water_body_fish_risk

ROOT = Path(__file__).resolve().parents[3]


def load(p):
    return yaml.safe_load((ROOT / p).read_text(encoding="utf-8"))


PLUME = load("modules/napl/examples/gas_station_plume.yaml")
RIVER = load("modules/surfacewater/examples/river_measured_soil.yaml")
HG = load("modules/surfacewater/examples/lake_deposition_mercury.yaml")


def river_benzene():
    w = copy.deepcopy(RIVER)
    w["chemical"] = {"cas": "71-43-2"}
    w["soil"] = {"concentration_mg_kg": 0.0}      # only the groundwater load
    return w


@pytest.mark.parametrize("flux", ["integrated", "bioscreen"])
def test_plume_flux_becomes_external_load(flux):
    res = plume_to_water_body(PLUME, river_benzene(), discharge_distance_m=60, model="first_order", flux=flux)
    p = copy.deepcopy(PLUME)
    p["general"]["model_length_m"] = 60
    pm = SourcePlumeModule()
    plume = pm.run(pm.validate(p))
    if flux == "integrated":
        expect = plume.integrated[1].mass_flux_mg_day[-1]
    else:
        expect = plume.mass_balance[1].mass_flux_mg_day[-1]
    assert res["mass_flux_mg_day"] == pytest.approx(expect, rel=1e-12)
    assert res["load_g_yr"] == pytest.approx(expect * 365 / 1000, rel=1e-12)
    w = river_benzene()
    w["external_loads"] = [{"name": "gw", "load_g_yr": res["load_g_yr"]}]
    wm = WatershedModule()
    direct = wm.run(wm.validate(w)).results[0]
    r = res["water_body"].results[0]
    assert r.loads.total_g_yr == pytest.approx(res["load_g_yr"], rel=1e-12)     # soil load is zero
    assert r.dissolved_mg_l == pytest.approx(direct.dissolved_mg_l, rel=1e-12) and r.dissolved_mg_l > 0


def test_fish_risk_measured_soil_matches_foodchain_workflow():
    res = water_body_fish_risk(RIVER)
    r = res["water_body"].results[0]
    direct = foodchain_risk([{"cas": "91-20-3", "surface_water_dissolved_mg_l": r.dissolved_mg_l,
                              "sediment_mg_kg": r.bed_sediment_mg_kg}], receptor="fisher")
    assert res["hazard_index"] == pytest.approx(direct["hazard_index"], rel=1e-12)
    assert res["total_risk"] == pytest.approx(direct["total_risk"], rel=1e-12)


def test_fish_risk_mercury_bases_and_species():
    res = water_body_fish_risk(HG)
    assert set(res["risk"]) == {"carcinogens", "noncarcinogens"}
    nc = res["risk"]["noncarcinogens"]
    assert [c.cas for c in nc["media"].chemicals] == ["7487-94-7", "22967-92-6"]
    assert res["hazard_index"] == nc["hazard_index"] > 0
    assert res["total_risk"] == res["risk"]["carcinogens"]["total_risk"]


def test_errors():
    with pytest.raises(ValueError, match="model"):
        plume_to_water_body(PLUME, river_benzene(), 60, model="x")
    with pytest.raises(ValueError, match="flux"):
        plume_to_water_body(PLUME, river_benzene(), 60, flux="x")
