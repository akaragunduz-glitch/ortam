"""surfacewater.watershed against an independent hand calculation and analytical checks.

Reference: data/reference/hhrap/2005/b4_hand_calculation.ods / .json: the equations of HHRAP 2005
Appendix B, Tables B-4-1 to B-4-25, transcribed as worksheet formulas (one column per case) and
recalculated by LibreOffice (build_b4_hand_calculation.py). The module must reproduce every row.
"""
import json
import math
from pathlib import Path

import pytest

from modules.surfacewater import equations as eq
from modules.surfacewater.watershed import WatershedModule, defaults

ROOT = Path(__file__).resolve().parents[3]
REF = json.loads((ROOT / "data/reference/hhrap/2005/b4_hand_calculation.json").read_text(encoding="utf-8"))["cases"]
M = WatershedModule()
REL = 1e-11


def close(a, b, rel=REL, abs_=1e-300):
    return abs(a - b) <= rel * max(abs(a), abs(b)) + abs_


def module_input(case: str, mercury=False) -> dict:
    i = REF[case]["inputs"]
    chem = {"name": case, "organic": i["er"] == 3, "kds_ml_g": i["kds"], "kdsw_l_kg": i["kdsw"], "kdbs_l_kg": i["kdbs"],
            "henry_atm_m3_mol": i["h"], "dw_cm2_s": i["dw"], "da_cm2_s": i["da"], "ksg_per_yr": i["ksg"], "fv": i["fv"]}
    if mercury:
        chem = {"cas": "7439-97-6"}
    data = {
        "chemical": chem, "chemical_table": "hhrap-2005" if mercury else "none",
        "watershed": {"area_m2": i["al"], "impervious_area_m2": i["ai"], "soil_mixing_depth_cm": i["zs"],
                      "bulk_density_g_cm3": i["bd"], "soil_water_content": i["tsw"], "particle_density_g_cm3": i["rhos"],
                      "air_temperature_k": i["ta"], "runoff_cm_yr": i["ro"], "precipitation_cm_yr": i["p"],
                      "irrigation_cm_yr": i["irr"], "evapotranspiration_cm_yr": i["ev"],
                      "usle": {"rainfall_rf_per_yr": i["rf"], "erodibility_k_ton_acre": i["k_usle"],
                               "length_slope_ls": i["ls"], "cover_management_c": i["c_usle"],
                               "supporting_practice_pf": i["pf"]},
                      "sd_b": i["b_sd"], "enrichment_ratio": i["er"],
                      "erosion_loss_constant": "calculated" if i["kse_mode"] == 1 else "zero"},
        "soil": {"concentration_mg_kg": i["cs_in"], "deposition_period_yr": i["td"], "t1_yr": i["t1"],
                 "exposure_duration_yr": i["t2"]},
        "water_body": {"type": "river" if i["river"] == 1 else "lake", "area_m2": i["aw"], "flow_m3_yr": i["vfx"],
                       "water_column_depth_m": i["dwc"], "benthic_depth_m": i["dbs"], "tss_mg_l": i["tss"],
                       "bed_sediment_concentration_g_cm3": i["cbs"], "bed_sediment_porosity": i["tbs"],
                       "temperature_k": i["twk"], "theta": i["theta"], "wind_speed_m_s": i["wind"],
                       "drag_coefficient": i["cd"], "von_karman": i["vk"], "viscous_sublayer": i["lz"],
                       "air_density_g_cm3": i["rho_a"], "water_density_g_cm3": i["rho_w"],
                       "air_viscosity_g_cm_s": i["mu_a"], "water_viscosity_g_cm_s": i["mu_w"]},
        "external_loads": [{"name": "external", "load_g_yr": i["l_ext"]}] if i["l_ext"] else [],
    }
    if i["river"] == 1:
        data["water_body"]["current_velocity_m_s"] = i["u"]
    if i["q"]:
        data["air"] = {"emission_rate_g_s": i["q"], "watershed_dytwv_s_m2_yr": i["dytwv_ws"],
                       "watershed_dytwp_s_m2_yr": i["dytwp_ws"], "water_body_dytwv_s_m2_yr": i["dytwv_wb"],
                       "water_body_dytwp_s_m2_yr": i["dytwp_wb"], "water_body_cywv_ug_s_g_m3": i["cywv"]}
    return data


CASES = [c for c in REF if not c.startswith("hg_")]


def basis_of(case):
    i = REF[case]["inputs"]
    if i["cs_in"] is not None:
        return "measured"
    return "noncarcinogens" if i["basis"] == 1 else "carcinogens"


@pytest.mark.parametrize("case", CASES)
def test_hand_calculation(case):
    v = REF[case]["values"]
    out = M.run(M.validate(module_input(case)))
    s, t = out.soil[0], out.transfer
    r = next(x for x in out.results if x.basis == basis_of(case))
    pairs = {
        "xe": t.unit_soil_loss_kg_m2_yr, "a_sd": t.sd_intercept_a, "sd": t.sediment_delivery_ratio,
        "kse": s.kse_per_yr, "ksr": s.ksr_per_yr, "ksl": s.ksl_per_yr, "ksv": s.ksv_per_yr, "ks": s.ks_per_yr,
        "kl": t.kl_m_yr, "kg": t.kg_m_yr, "kv_overall": t.kv_overall_m_yr, "fwc": t.fwc, "fbs": t.fbs,
        "kv": t.kv_per_yr, "kb_raw": t.kb_calculated_per_yr, "kb": t.kb_per_yr, "kwt": t.kwt_per_yr,
        "cs": r.soil_mg_kg["COPC"], "ldep": r.loads.deposition_g_yr, "lri": r.loads.impervious_runoff_g_yr,
        "lr": r.loads.pervious_runoff_g_yr, "le": r.loads.erosion_g_yr, "ldif": r.loads.diffusion_g_yr,
        "lt": r.loads.total_g_yr, "cwtot": r.total_water_body_mg_l, "cwctot": r.total_water_column_mg_l,
        "cdw": r.dissolved_mg_l, "csb": r.bed_sediment_mg_kg,
    }
    if s.ds_mg_kg_yr is not None:
        pairs.update({"ds": s.ds_mg_kg_yr, "cstd": s.cs_td_mg_kg, "cs_dep": s.cs_average_mg_kg})
    for k, ours in pairs.items():
        assert close(ours, v[k]), f"{case} {k}: {ours!r} vs hand calculation {v[k]!r}"


def test_hand_calculation_mercury():
    out = M.run(M.validate(module_input("hg_water_body", mercury=True)))
    hg2, mhg, w = REF["hg_soil_hg2"]["values"], REF["hg_soil_mhg"]["values"], REF["hg_water_body"]["values"]
    soil = {s.species: s for s in out.soil}
    assert close(soil["Hg2+"].ks_per_yr, hg2["ks"]) and close(soil["MHg"].ks_per_yr, mhg["ks"])
    assert close(soil["Hg2+"].ds_mg_kg_yr, hg2["ds"]) and close(soil["MHg"].ds_mg_kg_yr, mhg["ds"])
    r = next(x for x in out.results if x.basis == "carcinogens")
    assert close(r.soil_mg_kg["Hg2+"], hg2["cs"]) and close(r.soil_mg_kg["MHg"], mhg["cs"])
    for k, ours in {"ldep": r.loads.deposition_g_yr, "lri": r.loads.impervious_runoff_g_yr, "ldif": r.loads.diffusion_g_yr,
                    "lr_sum": r.loads.pervious_runoff_g_yr, "le_sum": r.loads.erosion_g_yr, "lt": r.loads.total_g_yr,
                    "fwc": out.transfer.fwc, "cwtot": r.total_water_body_mg_l, "cdw": r.dissolved_mg_l,
                    "csb": r.bed_sediment_mg_kg, "cdw_hg2": r.dissolved_species_mg_l["Hg2+"],
                    "cdw_mhg": r.dissolved_species_mg_l["MHg"]}.items():
        assert close(ours, w[k]), f"mercury {k}: {ours!r} vs {w[k]!r}"
    assert out.mercury and any("Mercury" in n for n in out.notes)


# --- Analytical checks ---------------------------------------------------------------------------
def test_soil_average_branches_continuous_at_t2_equal_td():
    ds, ks, td = 2.0, 0.3, 30.0
    a = eq.soil_conc_average(ds, ks, td, 0.0, td)               # branch T2 <= tD
    b = eq.soil_conc_average(ds, ks, td, 0.0, td * (1 + 1e-12))  # branch T1 < tD < T2
    assert a == pytest.approx(b, rel=1e-9)


def test_soil_small_ks_limits():
    ds, td = 1.5, 30.0
    assert eq.soil_conc_at_td(ds, 1e-9, td) == pytest.approx(ds * td, rel=1e-6)        # CstD -> Ds tD
    assert eq.soil_conc_average(ds, 1e-6, td, 0.0, td) == pytest.approx(ds * td / 2, rel=1e-4)   # mean of a ramp


def test_kv_limits():
    kl, kg, r, t = 1.0, 100.0, 8.205e-5, 293.0
    assert eq.overall_transfer(kl, kg, 1e3, r, t, 1.026) == pytest.approx(kl, rel=1e-4)     # liquid-film control
    h = 1e-9
    assert eq.overall_transfer(kl, kg, h, r, t, 1.026) == pytest.approx(kg * h / (r * t), rel=1e-4)


def test_no_losses_water_body_is_flow_dilution():
    assert eq.total_water_body(100.0, 1e6, 1.0, 0.0, 1e4, 2.0, 0.03) == pytest.approx(1e-4)
    assert eq.dissolved(0.5, 0.0, 10.0) == 0.5
    assert eq.fraction_water_column(0.0, 10, 2.0, 0.0, 0.6, 0.0, 1.0) == 1.0


def test_usle_units_and_sd_table():
    assert eq.unit_soil_loss(1, 1, 1, 1, 1) == pytest.approx(907.18 / 4047)
    sdd = defaults()["sediment_delivery"]
    for area_sq_mi, a in ((0.05, 2.1), (0.1, 2.1), (0.5, 1.9), (5, 1.4), (50, 1.2), (500, 0.6)):
        assert eq.sd_intercept(area_sq_mi * 2.59e6, sdd["a_by_area_sq_mi"], sdd["m2_per_sq_mi"]) == a


def test_defaults_from_tables():
    d = defaults()
    assert d["gas_transfer_river_m_yr"]["value"] == 36500 and d["tss_mg_l"]["value"] == 10
    assert d["benthic_depth_m"]["value"] == 0.03 and d["bed_sediment_porosity"]["value"] == 0.6
    assert d["usle"]["erodibility_k_ton_acre"]["value"] == 0.39


def test_negative_kb_set_to_zero_with_warning():
    out = M.run(M.validate(module_input("river_negative_kb_kse_irr")))
    assert out.transfer.kb_calculated_per_yr < 0 and out.transfer.kb_per_yr == 0
    assert any("kb" in w for w in out.warnings)


def test_external_load_is_additive():
    data = module_input("river_measured_organic")
    base = M.run(M.validate(data)).results[0]
    data["external_loads"] = [{"name": "groundwater", "load_g_yr": 1000.0}]
    ext = M.run(M.validate(data)).results[0]
    assert ext.loads.total_g_yr == pytest.approx(base.loads.total_g_yr + 1000.0)
    assert ext.total_water_body_mg_l / base.total_water_body_mg_l == pytest.approx(
        ext.loads.total_g_yr / base.loads.total_g_yr)


def test_database_parameters_and_enrichment_ratio():
    data = module_input("river_measured_organic")
    data["chemical"] = {"cas": "71-43-2"}
    data["chemical_table"] = "hhrap-2005"
    data["watershed"].pop("enrichment_ratio")
    out = M.run(M.validate(data))
    assert out.chemical == "Benzene" and out.parameters["COPC"]["kds_ml_g"] == 0.12
    v = REF["river_measured_organic"]["values"]
    assert close(out.results[0].dissolved_mg_l, v["cdw"])      # database values = hand-calculation inputs, ER 3


def test_validation():
    data = module_input("river_measured_organic")
    data["soil"] = {}
    with pytest.raises(ValueError, match="deposition"):
        M.validate(data)
    data = module_input("river_measured_organic")
    data["water_body"].pop("current_velocity_m_s")
    with pytest.raises(ValueError, match="current_velocity"):
        M.validate(data)


@pytest.mark.parametrize("example", sorted((ROOT / "modules/surfacewater/examples").glob("*.yaml")), ids=lambda p: p.name)
def test_examples_run(example):
    import yaml
    from modules.registry import get_module
    m = get_module("surfacewater.watershed")
    out = m.run(m.validate(yaml.safe_load(example.read_text(encoding="utf-8"))))
    assert out.results and m.summary(out)
