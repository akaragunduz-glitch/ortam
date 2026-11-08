"""workflows.foodchain: soil -> food concentrations -> ingestion risk."""

import pytest

from workflows.foodchain import foodchain_risk


def test_chain_is_consistent():
    out = foodchain_risk([{"cas": "71-43-2", "soil_mg_kg": 1.0}], receptor="farmer")
    media, ing = out["media"].chemicals[0], out["ingestion"].chemicals[0]
    assert ing.intake_mg_kg_d["beef"] == pytest.approx(media.beef_mg_kg_fw * 0.00122)
    assert ing.intake_mg_kg_d["bg"] == pytest.approx(media.produce_belowground_mg_kg_dw * 0.00017)
    assert out["total_risk"] == pytest.approx(out["ingestion"].total_risk)


def test_direct_soil_added():
    a = foodchain_risk([{"cas": "7440-38-2", "soil_mg_kg": 10.0}], receptor="farmer_child")
    b = foodchain_risk([{"cas": "7440-38-2", "soil_mg_kg": 10.0}], receptor="farmer_child", direct_soil=True)
    assert b["total_risk"] == pytest.approx(a["total_risk"] + b["soil"].total_risk)
    assert b["hazard_index"] == pytest.approx(a["hazard_index"] + b["soil"].hazard_index_child)
