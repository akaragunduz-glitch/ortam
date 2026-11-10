"""Workflow napl.equilibrium -> napl.source-plume -> risk.resident-tapwater."""
from pathlib import Path

import pytest
import yaml

from modules.napl.equilibrium import NaplEquilibriumModule
from modules.napl.source_plume import SourcePlumeModule
from modules.risk.resident_tapwater import ResidentTapwaterModule
from workflows.napl_plume import napl_plume_risk

ROOT = Path(__file__).resolve().parents[3]
EX = ROOT / "modules/napl/examples"


def load(name):
    return yaml.safe_load((EX / name).read_text(encoding="utf-8"))


def test_effective_solubility_becomes_source_and_risk_matches_modules():
    plume, napl = load("gas_station_plume.yaml"), load("gasoline_composition.yaml")
    res = napl_plume_risk(plume, model="first_order", napl_inputs=napl, component="benzene")
    # 1. effective solubility of benzene from napl.equilibrium
    eqm = NaplEquilibriumModule()
    seff = next(c for c in eqm.run(eqm.validate(napl)).components if c.name == "benzene").effective_solubility_mg_l
    assert res["plume"].derived.layer_increments_mg_l[2] == pytest.approx(seff - 0.8)   # central layer increment
    # 2. the plume run equals a direct module run with that concentration
    plume["source"]["zones"][0]["concentration_mg_l"] = seff
    pm = SourcePlumeModule()
    direct = pm.run(pm.validate(plume))
    for a, b in zip(res["receptors"], direct.receptors):
        assert a["concentration_mg_l"] == pytest.approx(b.first_order_mg_l, rel=1e-12)
    # 3. risk equals risk.resident-tapwater at that concentration
    tw = ResidentTapwaterModule()
    r0 = res["receptors"][0]
    o = tw.run(tw.validate({"chemicals": [{"cas": "71-43-2", "concentration_ug_l": r0["concentration_mg_l"] * 1000}]}))
    assert r0["risk"] == pytest.approx(o.total_risk, rel=1e-12)
    assert r0["exceeds_sl"] == (r0["concentration_mg_l"] * 1000 > o.chemicals[0].sl_ug_l)


def test_measured_source_without_equilibrium():
    res = napl_plume_risk(load("gas_station_plume.yaml"), model="instantaneous")
    assert res["equilibrium"] is None and len(res["receptors"]) == 2
    assert res["receptors"][0]["concentration_mg_l"] == pytest.approx(res["plume"].receptors[0].instantaneous_mg_l)


def test_errors():
    with pytest.raises(ValueError, match="model"):
        napl_plume_risk(load("gas_station_plume.yaml"), model="x")
    with pytest.raises(ValueError, match="component"):
        napl_plume_risk(load("gas_station_plume.yaml"), napl_inputs=load("gasoline_composition.yaml"))
