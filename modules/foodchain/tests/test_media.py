"""foodchain.media: HHRAP 2005 Appendix B equations, by hand with the database transfer factors."""

from __future__ import annotations

import pytest

from modules.foodchain.media import FoodchainMediaModule


def run(chemicals, **kw):
    m = FoodchainMediaModule()
    return m.run(m.validate({"chemicals": chemicals, **kw}))


def test_benzene_soil_chain_by_hand():
    # HHRAP database: Br_ag = Br_forage = Br_grain = 2.37, Br_rootveg = 80.1, log Kow 2.1 (VG = 1),
    # Ba_beef = 3.3819e-3, Ba_milk = 7.1198e-4 day/kg FW
    cs = 10.0
    r = run([{"cas": "71-43-2", "soil_mg_kg": cs}]).chemicals[0]
    assert r.produce_protected_aboveground_mg_kg_dw == pytest.approx(cs * 2.37)              # B-2-9
    assert r.produce_belowground_mg_kg_dw == pytest.approx(cs * 80.1 * 1.0)                 # B-2-10
    p = cs * 2.37                                                                              # B-3-9
    beef_intake = 8.8 * p + 2.5 * p + 0.47 * p + 0.5 * cs * 1.0                              # B-3-10
    assert r.beef_mg_kg_fw == pytest.approx(beef_intake * 0.0033819221094889646)
    milk_intake = 13.2 * p + 4.1 * p + 3.0 * p + 0.4 * cs                                     # B-3-11
    assert r.milk_mg_kg_fw == pytest.approx(milk_intake * 0.0007119836019976768)
    pork_intake = 1.4 * p + 3.3 * p + 0.37 * cs                                               # B-3-12
    assert r.pork_mg_kg_fw == pytest.approx(pork_intake * 0.004093905711486641)
    egg_intake = 0.2 * p + 0.022 * cs                                                         # B-3-13, B-3-14
    assert r.eggs_mg_kg_fw == pytest.approx(egg_intake * 0.0014239672039953535)
    assert r.chicken_mg_kg_fw == pytest.approx(egg_intake * 0.002491942606991869)


def test_air_contribution_added_to_exposed_produce_and_forage_only():
    r0 = run([{"cas": "71-43-2", "soil_mg_kg": 1.0}]).chemicals[0]
    r1 = run([{"cas": "71-43-2", "soil_mg_kg": 1.0,
               "air": {"exposed_produce_mg_kg": 0.5, "forage_mg_kg": 0.2}}]).chemicals[0]
    assert r1.produce_exposed_aboveground_mg_kg_dw == pytest.approx(r0.produce_exposed_aboveground_mg_kg_dw + 0.5)
    assert r1.produce_protected_aboveground_mg_kg_dw == r0.produce_protected_aboveground_mg_kg_dw
    assert r1.beef_mg_kg_fw == pytest.approx(r0.beef_mg_kg_fw + 8.8 * 0.2 * 0.0033819221094889646)


def test_behp_metabolism_factor_and_vg():
    # bis(2-ethylhexyl)phthalate: MF = 0.01 for beef, milk, pork (Table B-3-10, p. B-140); log Kow 5.1 > 4 -> VG 0.01
    r = run([{"cas": "117-81-7", "soil_mg_kg": 1.0}]).chemicals[0]
    p = 0.0437
    beef = (11.77 * p + 0.5) * 0.0398807760670795 * 0.01
    assert r.beef_mg_kg_fw == pytest.approx(beef)
    assert r.chicken_mg_kg_fw == pytest.approx((0.2 * p + 0.022) * 0.029385834996795424)   # no MF
    assert r.produce_belowground_mg_kg_dw == pytest.approx(0.000854 * 0.01)


def test_fish_method_selection():
    # benzene: BCF 8.26 (log Kow < 4); BEHP: BAF 194.29 (log Kow >= 4); BSAF needs sediment
    rb = run([{"cas": "71-43-2", "surface_water_dissolved_mg_l": 0.01}]).chemicals[0]
    assert rb.fish_mg_kg_fw == pytest.approx(0.01 * 8.26) and rb.fish_method.startswith("BCF")
    rp = run([{"cas": "117-81-7", "surface_water_dissolved_mg_l": 0.01}]).chemicals[0]
    assert rp.fish_mg_kg_fw == pytest.approx(0.01 * 194.29390242614014) and rp.fish_method.startswith("BAF")
    rs = run([{"cas": "71-43-2", "sediment_mg_kg": 2.0, "bsaf_fish": 0.5}]).chemicals[0]
    assert rs.fish_mg_kg_fw == pytest.approx(2.0 * 0.07 * 0.5 / 0.04)                        # B-4-28


def test_arsenic_uses_database_metal_factors():
    r = run([{"cas": "7440-38-2", "soil_mg_kg": 20.0}]).chemicals[0]
    assert r.produce_protected_aboveground_mg_kg_dw == pytest.approx(20 * 0.00633)
    assert r.forage_mg_kg_dw == pytest.approx(20 * 0.036) and r.grain_mg_kg_dw == pytest.approx(20 * 0.004)
    assert r.pork_mg_kg_fw is None and any("pork" in n for n in r.notes)


def test_chemical_not_in_hhrap_uses_correlations():
    # metolachlor is not in the HHRAP database: Appendix A-2 correlations with the RSL log Kow 3.13, Koc 488.5
    r = run([{"cas": "51218-45-2", "soil_mg_kg": 1.0}]).chemicals[0]
    assert r.factors["br_ag"] == pytest.approx(10 ** (1.588 - 0.578 * 3.13))
    assert r.factors["br_rootveg"] == pytest.approx(10 ** (0.77 * 3.13 - 1.52) / 0.13 / (488.5 * 0.01))
    assert r.factors["ba_beef"] == pytest.approx(10 ** (-0.099 * 3.13 ** 2 + 1.07 * 3.13 - 3.56) * 0.19)
    assert "correlation" in r.factor_source["br_ag"] and r.notes


def test_user_override():
    r = run([{"cas": "71-43-2", "soil_mg_kg": 1.0, "br_ag": 1.0}]).chemicals[0]
    assert r.produce_protected_aboveground_mg_kg_dw == 1.0 and r.factor_source["br_ag"] == "user"
