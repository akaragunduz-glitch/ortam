"""napl.source-plume (SI module) against BIOSCREEN 1.4.

1. Every workbook reference case, given in SI, reproduces the workbook (1e-9): checks the unit
   conversion at the module boundary.
2. The printed results of the two manual examples, at their printed precision:
   Example 1, Keesler AFB (BIOSCREEN 1.4 Revisions, Fig. 3, 4, 6; PDF p. 82-85);
   Example 2, Hill AFB Site 870 (User's Manual 1.3, Fig. 7, 8; PDF p. 98-99).
3. Analytical checks of the numerical integration and of the consistent mode.
"""
import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest
import yaml

from modules.napl.source_plume import FT, FT3, SourcePlumeModule

ROOT = Path(__file__).resolve().parents[3]
REF = json.loads((ROOT / "data/reference/bioscreen/v1.4/reference_cases.json").read_text(encoding="utf-8"))
CASES = {c["id"]: c for c in REF["cases"]}
EX1 = yaml.safe_load((ROOT / "modules/napl/examples/bioscreen_keesler_example1.yaml").read_text(encoding="utf-8"))
ACRE_FT_M3 = 43560 * FT3
M = SourcePlumeModule()


def run(data):
    return M.run(M.validate(data))


def si_input(rec) -> dict:
    """Workbook inputs (after recalculation) as SI module input."""
    ic = rec["input_cells"]
    zones = []
    for y, cc in (("y3_ft", "c3_mg_l"), ("y2_ft", "c2_mg_l"), ("y1_ft", "c1_mg_l")):
        if ic[y]:
            zones.append({"width_m": ic[y] * FT, "concentration_mg_l": ic[cc]})
    return {
        "chemical_table": "none", "instantaneous_source": "bioscreen",
        "hydrogeology": {"seepage_velocity_m_yr": ic["vs_ft_yr"] * FT, "effective_porosity": ic["porosity"]},
        "dispersion": {"alpha_x_m": ic["alpha_x_ft"] * FT, "alpha_y_m": ic["alpha_y_ft"] * FT,
                       "alpha_z_m": (ic["alpha_z_ft"] or 0) * FT},
        "adsorption": {"retardation_factor": ic["retardation"]},
        "biodegradation": {"decay_rate_per_yr": ic["lambda_per_yr"], "delta_oxygen_mg_l": ic["do_mg_l"] or 0,
                           "delta_nitrate_mg_l": ic["no3_mg_l"] or 0, "ferrous_iron_mg_l": ic["fe_mg_l"] or 0,
                           "delta_sulfate_mg_l": ic["so4_mg_l"] or 0, "methane_mg_l": ic["ch4_mg_l"] or 0},
        "general": {"model_length_m": ic["length_ft"] * FT, "model_width_m": ic["width_ft"] * FT,
                    "simulation_time_yr": ic["time_yr"]},
        "source": {"thickness_m": ic["thickness_ft"] * FT, "zones": zones,
                   "soluble_mass_kg": None if ic["mass_kg"] == "Infinite" else ic["mass_kg"]},
    }


def close(a, b, rel=1e-9, abs_=1e-12):
    return abs(a - b) <= rel * max(abs(a), abs(b)) + abs_


# Cases whose zone layout the SI input can express (outer zones need inner zones; no zero increments)
SI_CASES = [k for k in CASES if k != "equal_zone_concentrations"]


@pytest.mark.parametrize("cid", SI_CASES)
def test_si_module_reproduces_workbook(cid):
    rec = CASES[cid]
    out = run(si_input(rec))
    o = rec["outputs"]
    for p, a, b, cc in zip(out.centerline, o["Centerline Output!B51:L51"], o["Centerline Output!B60:L60"],
                           o["Centerline Output!B83:L83"]):
        assert close(p.no_decay_mg_l, a) and close(p.first_order_mg_l, b) and close(p.instantaneous_mg_l, cc)
    for key, cells in (("no_decay_mg_l", "C88:M92"), ("first_order_mg_l", "C103:M107"),
                       ("instantaneous_mg_l", "C119:M123")):
        ours = np.array(getattr(out.array, key)).ravel()
        ref = np.array(o[f"Plume Output!{cells}"]).ravel()
        assert all(close(a, b) for a, b in zip(ours, ref)), key
    d = out.derived
    assert close(d.source_flux_kg_yr, o["Input!F61"]) and close(d.biodegradation_capacity_mg_l, o["Input!F52"])
    assert close(d.source_flow_m3_yr, o["Input!F59"] / 28.3 * FT3)
    for mb, key in zip(out.mass_balance, ("t", "f", "i")):
        ref = rec["models"][key]
        if ref["M17"] != "Can't Calc.":
            assert close(mb.actual_plume_mass_kg, ref["M19"])
            assert close(mb.plume_volume_m3 / ACRE_FT_M3, ref["M30"])
            assert close(mb.source_flow_m3_yr / ACRE_FT_M3, ref["M31"])
        if mb.mass_flux_mg_day is not None:
            assert all(close(a, b) for a, b in zip(mb.mass_flux_mg_day, ref["B11:L11"]))


# --- Example 1: Keesler AFB (printed, 3 decimals) --------------------------------------------------
FIG3 = {  # BSCR14 Fig. 3, p. A-7 [82]
    "no_decay": [13.544, 6.575, 5.280, 4.581, 4.107, 3.754, 3.474, 3.241, 3.040, 2.861, 2.697],
    "first_order": [13.544, 3.117, 1.186, 0.488, 0.208, 0.090, 0.040, 0.018, 0.008, 0.004, 0.002],
    "instantaneous": [12.021, 5.463, 4.248, 3.500, 2.860, 2.257, 1.678, 1.114, 0.559, 0.004, 0.000],
}
FIG4_FLUX = [9.1e3, 7.2e3, 5.6e3, 4.6e3, 3.8e3, 3.0e3, 2.2e3, 1.5e3, 7.4e2, 5.2e0, 0.0]   # mg/day, Fig. 4
FIG6_FLUX = [9.1e3, 5.9e3, 5.1e3, 4.1e3, 3.2e3, 2.3e3, 1.5e3, 8.6e2, 3.2e2, 1.3e0, 0.0]   # Fig. 6 (W = 50 ft)
FIG6_ROWS = {  # Fig. 6, instantaneous model, y = +-12.5 ft and +-25 ft
    "dy": [1.554, 4.231, 3.603, 2.984, 2.391, 1.821, 1.272, 0.737, 0.210, 0.000, 0.000],
    "2dy": [1.554, 2.024, 1.962, 1.526, 1.036, 0.556, 0.096, 0.000, 0.000, 0.000, 0.000],
}


def sig2(v):
    return 0.0 if v == 0 else float(f"{v:.1e}")


@pytest.fixture(scope="module")
def ex1():
    return run(EX1)


def test_example1_centerline(ex1):
    for k, printed in FIG3.items():
        ours = [getattr(p, f"{k}_mg_l") for p in ex1.centerline]
        assert [round(v, 3) for v in ours] == pytest.approx(printed, abs=1e-12), k


def test_example1_inputs_and_source(ex1):
    d = ex1.derived
    assert d.seepage_velocity_m_yr / FT == pytest.approx(113.8, abs=0.05)          # Fig. 2: 113.8 ft/yr
    assert round(d.source_half_life_inst_yr, -1) == 60 and round(d.source_half_life_yr, -2) == 400   # Fig. 2
    assert d.retardation_factor == pytest.approx(1.0, abs=0.05)


def test_example1_mass_balance(ex1):
    i = next(m for m in ex1.mass_balance if m.model == "Instantaneous reaction")
    assert round(i.mass_left_source_kg, 1) == 126.3            # Fig. 4 (and BSCR14 p. A-4)
    assert round(i.actual_plume_mass_kg, 1) == 7.8
    assert round(i.mass_removed_kg, 1) == 118.6
    assert round(100 * i.fraction_removed) == 94
    ea = {k: round(v, 1) for k, v in i.electron_acceptor_change_kg.items()}
    assert ea == {"oxygen": -13.3, "nitrate": -5.7, "ferrous_iron": 134.3, "sulfate": -181.2, "methane": 53.4}
    assert round(i.source_mass_now_kg, 1) == 1873.7            # Fig. 4 (the text on p. A-4 misprints 1837)
    assert round(i.plume_volume_m3 / ACRE_FT_M3, 1) == 1.0
    assert round(i.source_flow_m3_yr / ACRE_FT_M3, 3) == 1.019
    assert [sig2(v) for v in i.mass_flux_mg_day] == FIG4_FLUX


def test_example1_width_50ft_array_and_flux():
    data = copy.deepcopy(EX1)
    data["general"]["model_width_m"] = 50 * FT
    out = run(data)
    i = next(m for m in out.mass_balance if m.model == "Instantaneous reaction")
    assert i.cant_calc and i.actual_plume_mass_kg is None       # Fig. 6: "Can't Calc."
    assert [sig2(v) for v in i.mass_flux_mg_day] == FIG6_FLUX
    a = out.array.instantaneous_mg_l
    for row, key in ((1, "dy"), (3, "dy"), (0, "2dy"), (4, "2dy")):
        assert [round(v, 3) for v in a[row]] == pytest.approx(FIG6_ROWS[key], abs=1e-12)
    assert any("Can't Calc" in w for w in out.warnings)


# --- Example 2: Hill AFB (printed; Version 1.3 output, alpha_z = 0) ----------------------------------
EX2 = {
    "chemical_table": "none",
    "hydrogeology": {"hydraulic_conductivity_m_s": 8.1e-5, "hydraulic_gradient": 0.048, "effective_porosity": 0.25},
    "dispersion": {"alpha_x_m": 28.5 * FT, "alpha_y_m": 2.85 * FT},
    "adsorption": {"bulk_density_kg_l": 1.7, "koc_l_kg": 38, "foc": 8e-4},
    "biodegradation": {"solute_half_life_yr": 0.1, "delta_oxygen_mg_l": 5.78, "delta_nitrate_mg_l": 17,
                       "ferrous_iron_mg_l": 11.3, "delta_sulfate_mg_l": 100, "methane_mg_l": 0.414},
    "general": {"model_length_m": 1450 * FT, "model_width_m": 320 * FT, "simulation_time_yr": 5},
    "source": {"thickness_m": 10 * FT, "zones": [{"width_m": 100 * FT, "concentration_mg_l": 9},
                                                {"width_m": 25 * FT, "concentration_mg_l": 2.8},
                                                {"width_m": 50 * FT, "concentration_mg_l": 0.07}]},
}
FIG7 = {  # BSCR Fig. 7, p. 63 [98]
    "no_decay": [9.000, 8.467, 7.466, 6.684, 6.089, 5.624, 5.250, 4.940, 4.679, 4.455, 4.260],
    "first_order": [9.000, 4.348, 1.969, 0.905, 0.424, 0.201, 0.096, 0.047, 0.023, 0.011, 0.005],
    "instantaneous": [9.000, 8.466, 7.407, 6.350, 5.268, 4.192, 3.152, 2.168, 1.245, 0.385, 0.000],
}


def test_example2_centerline_and_array():
    out = run(EX2)
    assert out.derived.seepage_velocity_m_yr / FT == pytest.approx(1609.1, abs=0.05)   # Fig. 6
    assert round(out.derived.retardation_factor, 1) == 1.2
    for k, printed in FIG7.items():
        ours = [getattr(p, f"{k}_mg_l") for p in out.centerline]
        # printed Version 1.3 values: largest difference 0.00074 mg/L
        assert ours == pytest.approx(printed, abs=1e-3), k
    a = out.array.instantaneous_mg_l       # Fig. 8: y = +-80 ft row
    assert [round(v, 3) for v in a[1][:3]] == [0.070, 0.538, 0.000]
    assert all(v == 0 for v in a[0])
    i = out.mass_balance[2]
    assert i.source_mass_now_kg is None and i.cant_calc                # "Infinite", "Can't Calc."


# --- Receptors, integration, consistent mode, validation -------------------------------------------
def test_receptor_equals_centerline_point(ex1):
    data = copy.deepcopy(EX1)
    data["receptors"] = [{"name": "x4", "x_m": ex1.centerline[4].x_m}]
    data["receptor_times_yr"] = [6]
    out = run(data)
    r, p = out.receptors[0], ex1.centerline[4]
    assert r.instantaneous_mg_l == pytest.approx(p.instantaneous_mg_l, rel=1e-12)
    assert r.time_series[0]["first_order_mg_l"] == pytest.approx(p.first_order_mg_l, rel=1e-12)


def test_integrated_mass_conservation():
    """No decay, constant source, R = 1, alpha_z = 0, plume inside the area: the integrated plume mass
    equals the mass that has left the source, Q C0 t (to the accuracy of the Domenico approximation)."""
    data = {
        "chemical_table": "none", "bioscreen_compatible": False,
        "hydrogeology": {"seepage_velocity_m_yr": 30, "effective_porosity": 0.3},
        "dispersion": {"alpha_x_m": 3, "alpha_y_m": 0.3},
        "general": {"model_length_m": 400, "model_width_m": 80, "simulation_time_yr": 5},
        "source": {"thickness_m": 3, "zones": [{"width_m": 10, "concentration_mg_l": 10}]},
        "integration_points": 801,
    }
    out = run(data)
    q_c_t = 30 * 0.3 * 3 * 10 * 10 * 5 / 1000      # m3/yr * mg/L * yr -> kg (1000 L/m3, 1e-6 kg/mg)
    nd = out.integrated[0]
    assert nd.plume_mass_kg == pytest.approx(q_c_t, rel=0.01)
    # mass flux across the first cross-section ~ source flux Q C0 (mg/day)
    assert nd.mass_flux_mg_day[1] == pytest.approx(30 * 0.3 * 3 * 10 * 10 * 1000 / 365, rel=0.02)


def test_consistent_mode_close_to_compatible_without_quirks():
    data = copy.deepcopy(EX1)
    data["source"]["soluble_mass_kg"] = None          # infinite source: Q1 has no effect
    a = run(data)
    data["bioscreen_compatible"] = False
    b = run(data)
    for p, q in zip(a.centerline, b.centerline):
        assert q.no_decay_mg_l == pytest.approx(p.no_decay_mg_l, rel=1e-6, abs=1e-6)
        assert q.instantaneous_mg_l == pytest.approx(p.instantaneous_mg_l, rel=1e-5, abs=1e-5)
    assert b.derived.source_flux_kg_yr / a.derived.source_flux_kg_yr == pytest.approx(28.316846592 / 28.3)


def test_koc_from_rsl_table():
    data = copy.deepcopy(EX1)
    data.pop("chemical"), data.pop("chemical_table")
    data["chemical"] = {"cas": "71-43-2"}
    data["adsorption"] = {"bulk_density_kg_l": 1.7, "foc": 0.001}
    out = run(data)
    assert any("Koc" in n and "Benzene" in n for n in out.notes)
    assert out.derived.retardation_factor > 1


def test_xu_eckstein_dispersivity_si():
    data = copy.deepcopy(EX1)
    data["dispersion"] = {"plume_length_m": 280 * FT}      # Example 1, original: alpha_x 13.3 ft, alpha_y 1.3 ft
    d = run(data).derived
    assert round(d.alpha_x_m / FT, 1) == 13.3 and round(d.alpha_y_m / FT, 1) == 1.3


@pytest.mark.parametrize("patch,msg", [
    ({"hydrogeology": {"effective_porosity": 0.3}}, "seepage_velocity"),
    ({"dispersion": {"alpha_z_m": 0}}, "alpha_x_m"),
    ({"biodegradation": {"decay_rate_per_yr": 1, "solute_half_life_yr": 1}}, "not both"),
])
def test_validation(patch, msg):
    data = copy.deepcopy(EX1)
    data.update(patch)
    with pytest.raises(ValueError, match=msg):
        M.validate(data)


def test_vertical_dispersion_disables_flux():
    data = copy.deepcopy(EX1)
    data["dispersion"]["alpha_z_m"] = 0.5 * FT
    out = run(data)
    assert out.integrated is None and all(m.mass_flux_mg_day is None for m in out.mass_balance)


def test_examples_run_from_cli_files():
    for f in ("bioscreen_keesler_example1.yaml", "gas_station_plume.yaml"):
        data = yaml.safe_load((ROOT / "modules/napl/examples" / f).read_text(encoding="utf-8"))
        out = run(data)
        assert out.centerline[0].no_decay_mg_l > 0 and M.summary(out)


# --- Mass-balance source of the instantaneous model (default) -------------------------------------
def test_default_instantaneous_source_is_mass_balance():
    data = copy.deepcopy(EX1)
    data.pop("instantaneous_source")
    out = run(data)
    d = out.derived
    assert d.instantaneous_source == "mass-balance"
    assert d.source_half_life_inst_yr is None and d.source_decay_rate_inst_per_yr is None
    b = d.biodegradation_capacity_mg_l / d.source_concentration_avg_mg_l
    assert d.source_depletion_time_inst_yr == pytest.approx(math.log((1 + b) / b) / d.source_decay_rate_per_yr)
    i = out.mass_balance[2]
    f = (1 + b) * math.exp(-d.source_decay_rate_per_yr * 6) - b
    assert i.source_mass_now_kg == pytest.approx(2000 * f) and i.mass_left_source_kg == pytest.approx(2000 * (1 - f))
    # no-degradation and first-order results are unchanged; instantaneous plume is longer
    ex1 = run(EX1)
    for p, q in zip(ex1.centerline, out.centerline):
        assert q.no_decay_mg_l == p.no_decay_mg_l and q.first_order_mg_l == p.first_order_mg_l
        assert q.instantaneous_mg_l >= p.instantaneous_mg_l - 1e-12
    assert any("mass-balance source" in n for n in out.notes)


def test_mass_balance_receptor_time_series_reaches_zero_after_depletion():
    data = copy.deepcopy(EX1)
    data["instantaneous_source"] = "mass-balance"
    te = run(data).derived.source_depletion_time_inst_yr
    data["receptors"] = [{"name": "source", "x_m": 0.0}]
    data["receptor_times_yr"] = [0.5 * te, 0.99 * te, 1.01 * te]
    ts = run(data).receptors[0].time_series
    assert ts[0]["instantaneous_mg_l"] > ts[1]["instantaneous_mg_l"] > 0 and ts[2]["instantaneous_mg_l"] == 0
