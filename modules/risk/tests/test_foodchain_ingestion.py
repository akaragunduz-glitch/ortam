"""risk.foodchain-ingestion: HHRAP 2005 Appendix C equations by hand."""

from __future__ import annotations

import pytest

from modules.risk.foodchain_ingestion import FoodchainIngestionModule


def run(chemicals, **kw):
    m = FoodchainIngestionModule()
    return m.run(m.validate({"chemicals": chemicals, **kw}))


FOODS = {"produce_exposed_aboveground_mg_kg_dw": 2.0, "produce_protected_aboveground_mg_kg_dw": 1.5,
         "produce_belowground_mg_kg_dw": 3.0, "beef_mg_kg_fw": 0.1, "milk_mg_kg_fw": 0.02,
         "pork_mg_kg_fw": 0.05, "chicken_mg_kg_fw": 0.03, "eggs_mg_kg_fw": 0.04, "fish_mg_kg_fw": 0.5}


def test_farmer_benzene_by_hand():
    o = run([{"cas": "71-43-2", "foods": FOODS}], receptor="farmer")
    c = o.chemicals[0]
    # Table C-1-2: I_ag = (2.0*0.00047 + 1.5*0.00064 + 3.0*0.00017) * 1
    assert c.intake_mg_kg_d["ag"] + c.intake_mg_kg_d["pp"] + c.intake_mg_kg_d["bg"] == pytest.approx(
        2.0 * 0.00047 + 1.5 * 0.00064 + 3.0 * 0.00017)
    # Table C-1-3: I_beef = 0.1 * 0.00122
    assert c.intake_mg_kg_d["beef"] == pytest.approx(0.1 * 0.00122)
    assert "fish" not in c.intake_mg_kg_d                          # farmer has no fish rate
    I = sum(c.intake_mg_kg_d.values())
    # Table C-1-7 with RSL SFO 0.055: risk = I*ED*EF*CSF/(AT*365), ED 40, EF 350, AT 70
    assert c.risk == pytest.approx(I * 40 * 350 * 0.055 / (70 * 365))
    # Table C-1-8 with RSL RfD 0.004: HQ = I*ED*EF/(RfD*ED*365)
    assert c.hq == pytest.approx(I * 350 / (0.004 * 365))
    assert o.total_risk == pytest.approx(c.risk) and o.hazard_index == pytest.approx(c.hq)


def test_fisher_child_and_doses():
    c = run([{"cas": "71-43-2", "foods": FOODS}], receptor="fisher_child").chemicals[0]
    assert c.intake_mg_kg_d["fish"] == pytest.approx(0.5 * 0.00088)        # Table C-1-4
    assert "beef" not in c.intake_mg_kg_d
    cdi = {d.segment: d.value for d in c.doses if d.basis == "cancer"}
    assert cdi["fish"] == pytest.approx(0.5 * 0.00088 * 6 * 350 / (70 * 365))
    assert sum(cdi.values()) * 0.055 == pytest.approx(c.risk)


def test_override_and_mutagen_note():
    o = run([{"cas": "50-32-8", "foods": {"beef_mg_kg_fw": 0.01}}], receptor="farmer",
            exposure={"consumption_kg_kg_d": {"beef": 0.002}, "ed_yr": 30})
    c = o.chemicals[0]
    assert c.intake_mg_kg_d["beef"] == pytest.approx(0.01 * 0.002)
    assert any("ADAF" in n for n in c.notes) and o.warnings
