"""BIOSCREEN 1.4 engine (modules/napl/bioscreen.py) against the workbook BIOSCRN4.xls.

Reference values: data/reference/bioscreen/v1.4/reference_cases.json, produced by
generate_bioscreen_reference_cases.py (LibreOffice recalculation of the workbook for 14 cases;
for the distributed data set the recalculation equals the values cached by Excel in the .xls).
Compatible mode must reproduce every compared cell to 1e-9 (relative).
"""
import json
import math
from pathlib import Path

import numpy as np
import pytest

from modules.napl import bioscreen as bs

ROOT = Path(__file__).resolve().parents[3]
REF = json.loads((ROOT / "data/reference/bioscreen/v1.4/reference_cases.json").read_text(encoding="utf-8"))
CASES = {c["id"]: c for c in REF["cases"]}
COMPAT = bs.Settings(compatible=True)
EA_KEYS = {"do_mg_l": "oxygen", "no3_mg_l": "nitrate", "fe_mg_l": "ferrous_iron", "so4_mg_l": "sulfate",
           "ch4_mg_l": "methane"}


def close(a, b, rel=1e-9, abs_=1e-12):
    return abs(a - b) <= rel * max(abs(a), abs(b)) + abs_


def assert_close(ours, ref, rel=1e-9, abs_=1e-12, what=""):
    ours = np.asarray(ours, dtype=float).ravel()
    ref = np.asarray(ref, dtype=float).ravel()
    assert ours.shape == ref.shape, what
    for k, (a, b) in enumerate(zip(ours, ref)):
        assert close(a, b, rel, abs_), f"{what}[{k}]: {a!r} vs workbook {b!r}"


def case_from_workbook(rec) -> bs.Case:
    ic = rec["input_cells"]
    mass = ic["mass_kg"]
    return bs.Case(
        vs_ft_yr=ic["vs_ft_yr"], porosity=ic["porosity"], alpha_x_ft=ic["alpha_x_ft"], alpha_y_ft=ic["alpha_y_ft"],
        alpha_z_ft=ic["alpha_z_ft"], retardation=ic["retardation"], lambda_per_yr=ic["lambda_per_yr"],
        length_ft=ic["length_ft"], width_ft=ic["width_ft"], time_yr=ic["time_yr"], thickness_ft=ic["thickness_ft"],
        zones=[(ic["y1_ft"] or 0.0, ic["c1_mg_l"] or 0.0), (ic["y2_ft"] or 0.0, ic["c2_mg_l"] or 0.0),
               (ic["y3_ft"] or 0.0, ic["c3_mg_l"] or 0.0)],
        mass_kg=None if mass == "Infinite" else mass,
        electron_acceptors={v: ic[k] or 0.0 for k, v in EA_KEYS.items()})


@pytest.fixture(scope="module", params=list(CASES))
def pair(request):
    rec = CASES[request.param]
    c = case_from_workbook(rec)
    return rec, c, bs.run(c, COMPAT)


def test_reference_file_integrity():
    assert len(CASES) == 14
    assert REF["sha256_xls"] and REF["sha256_zip"]


def test_erf_bioscreen_accuracy():
    # Chebyshev approximation of erfc, fractional error < 1.2e-7 (XLS 'Macro Module', gerf)
    from scipy.special import erf, erfc
    x = np.linspace(-6, 6, 2001)
    assert np.max(np.abs(bs.erf_bioscreen(x) - erf(x))) < 1.2e-7
    xp = x[(x > 0) & (x <= 2.5)]   # beyond, 1 - erf loses the digits of erfc
    assert np.max(np.abs((1 - bs.erf_bioscreen(xp)) / erfc(xp) - 1)) < 1.2e-7


def test_excel_cached_values_equal_recalculation():
    """The LibreOffice recalculation equals Excel's cached values for the distributed data set."""
    cached = REF["excel_cached_values"]
    rec = CASES["keesler_as_distributed"]
    for key in ("Centerline Output!B51:L51", "Centerline Output!B60:L60", "Centerline Output!B83:L83"):
        assert_close(rec["outputs"][key], cached[key], rel=1e-12, what=key)
    m = rec["models"][cached["Plume Output!P4"]]
    for cell in ("M17", "M19", "M21", "M28", "M30", "M31"):
        assert close(m[cell], cached[f"Plume Output!{cell}"], 1e-12), cell


def test_input_relations(pair):
    rec, c, _ = pair
    ic, given = rec["input_cells"], rec["inputs"]
    if "vs_ft_yr" not in given:      # E5 formula
        assert close(bs.seepage_velocity_ft_yr(ic["k_cm_s"], ic["gradient"], ic["porosity"]), ic["vs_ft_yr"])
    if "alpha_x_ft" not in given:    # E12 formula
        assert close(bs.alpha_x_xu_eckstein_ft(ic["plume_length_ft"], COMPAT), ic["alpha_x_ft"])
    if "alpha_y_ft" not in given:    # E13 formula
        assert close(bs.alpha_y_default(ic["alpha_x_ft"]), ic["alpha_y_ft"])
    if "retardation" not in given:   # E19 formula
        assert close(bs.retardation_factor(ic["rho_kg_l"], ic["koc_l_kg"], ic["foc"], ic["porosity"]), ic["retardation"])
    if "lambda_per_yr" not in given:  # E26 formula
        assert close(bs.decay_rate_per_yr(ic["half_life_yr"], COMPAT), ic["lambda_per_yr"])


def test_source_terms(pair):
    rec, c, r = pair
    o, s = rec["outputs"], r.source
    assert_close(s.widths_ft, [o["Input!L45"], o["Input!L46"], o["Input!L47"]], what="Y")
    assert_close(s.increments_mg_l, [o["Input!Q45"], o["Input!Q46"], o["Input!Q47"]], what="C")
    assert_close(s.bioc_mg_l, [o["Input!T45"], o["Input!T46"], o["Input!T47"]], what="bioc")
    for ours, cell in ((s.biodegradation_capacity_mg_l, "F52"), (s.flow_l_yr, "F59"), (s.c_avg_mg_l, "F60"),
                       (s.flux_kg_yr, "F61"), (s.ks_per_yr, "F62"), (s.c_avg_inst_mg_l, "F65"),
                       (s.flux_inst_kg_yr, "F66"), (s.ks_inst_per_yr, "F67")):
        assert close(ours, o[f"Input!{cell}"]), cell
    for ours, cell in ((s.half_life_yr, "F63"), (s.half_life_inst_yr, "F68")):
        ref = o[f"Input!{cell}"]
        assert (math.isinf(ours) and ref == "Infinite") or close(ours, ref), cell


def test_centerline(pair):
    rec, c, r = pair
    o = rec["outputs"]
    assert_close(r.x_ft, o["Centerline Output!B35:L35"], what="x")
    assert_close(r.centerline["no_decay"], o["Centerline Output!B51:L51"], what="no decay")
    assert_close(r.centerline["first_order"], o["Centerline Output!B60:L60"], what="first order")
    assert_close(r.centerline["instantaneous"], o["Centerline Output!B83:L83"], what="instantaneous")


def test_arrays(pair):
    rec, c, r = pair
    o = rec["outputs"]
    assert_close(r.x_ft, o["Plume Output!C57:M57"], what="x")
    assert_close(np.abs(r.y_ft), o["Plume Output!B88:B92"], what="y")
    assert_close(r.arrays["no_decay"], o["Plume Output!C88:M92"], what="no decay")
    assert_close(r.arrays["first_order"], o["Plume Output!C103:M107"], what="first order")
    assert_close(r.inst_before_array, o["Plume Output!C111:M115"], what="inst. before BC")
    assert_close(r.arrays["instantaneous"], o["Plume Output!C119:M123"], what="instantaneous")
    for name, key in (("no_decay", "t"), ("first_order", "f"), ("instantaneous", "i")):
        assert_close(r.arrays[name], rec["models"][key]["B5:L9"], what=f"display {key}")


def test_grid_masses(pair):
    rec, c, r = pair
    o = rec["outputs"]
    for ours, cell in ((r.mass_left_source_kg, "P88"), (r.mass_left_source_kg, "P101"),
                       (r.grid_masses_kg["no_decay"], "P91"), (r.grid_masses_kg["first_order"], "P104"),
                       (r.fraction_biodegraded["first_order"], "P107"), (r.mass_left_source_inst_kg, "P111"),
                       (r.grid_masses_kg["inst_before"], "P116"), (r.grid_masses_kg["instantaneous"], "P121"),
                       (r.fraction_biodegraded["instantaneous"], "P125")):
        assert close(ours, o[f"Plume Output!{cell}"]), cell


@pytest.mark.parametrize("model,key", [("no_decay", "t"), ("first_order", "f"), ("instantaneous", "i")])
def test_mass_balance_screen(pair, model, key):
    rec, c, r = pair
    ref, m = rec["models"][key], r.models[model]
    if ref["M17"] == "Can't Calc.":
        assert m.cant_calc and m.actual_plume_mass_kg is None
        assert ref["M30"] == "Can't Calc." and m.plume_volume_ft3 is None
    else:
        assert not m.cant_calc
        assert close(m.plume_mass_if_no_biodegradation_kg, ref["M17"])
        assert close(m.actual_plume_mass_kg, ref["M19"])
        assert close(m.mass_removed_kg, ref["M21"], abs_=1e-9)
        assert ref["N22"] == "0 %" if key == "t" else close(-m.fraction_removed, ref["N22"])
        assert close(m.plume_volume_ft3 / bs.FT2_PER_ACRE, ref["M30"])
        assert close(m.source_flow_ft3_yr / bs.FT2_PER_ACRE, ref["M31"])
        if key == "i":
            if m.electron_acceptor_change_kg is None:      # BC = 0: workbook divides by zero
                assert all(isinstance(v, str) for v in ref["I25:M25"])
            else:
                assert_close(list(m.electron_acceptor_change_kg.values()), ref["I25:M25"], what="EA")
        else:
            assert ref["I25:M25"] == ["na"] * 5
    if ref["M28"] == "Infinite":
        assert m.source_mass_now_kg is None
    else:
        assert close(m.source_mass_now_kg, ref["M28"])


@pytest.mark.parametrize("model,key", [("no_decay", "t"), ("first_order", "f"), ("instantaneous", "i")])
def test_mass_flux(pair, model, key):
    rec, c, r = pair
    ref, m = rec["models"][key]["B11:L11"], r.models[model]
    if c.alpha_z_ft != 0:
        assert m.mass_flux_mg_day is None and ref == ["-"] * 11
    else:
        assert_close(m.mass_flux_mg_day, ref, what="flux")


# --- Quirks: compatible mode reproduces them, consistent mode corrects them -----------------------
def test_q1_middle_layer_source_decay_in_instantaneous_model():
    """Q1: with a finite source and a middle layer, BIOSCREEN decays the middle layer with ks."""
    c = case_from_workbook(CASES["keesler_example1"])
    src_c = bs.source_terms(c, COMPAT)
    cons = bs.Settings(compatible=False)
    # Compatible: at the source column (centreline) the instantaneous concentration equals
    # bioc1 e^(-ksi t) + bioc2 e^(-ks t) + bioc3 e^(-ksi t) - BC
    t = c.time_yr - bs.X_SOURCE_FT * c.retardation / c.vs_ft_yr    # source time at x = 1e-4 ft
    b1, b2, b3 = src_c.bioc_mg_l
    expect = b1 * math.exp(-src_c.ks_inst_per_yr * t) + b2 * math.exp(-src_c.ks_per_yr * t) \
        + b3 * math.exp(-src_c.ks_inst_per_yr * t) - src_c.biodegradation_capacity_mg_l
    got = bs.concentrations(c, src_c, COMPAT, [0.0], [0.0]).instantaneous[0]
    assert close(got, max(expect, 0.0))
    src_n = bs.source_terms(c, cons)
    expect_n = sum(b * math.exp(-src_n.ks_inst_per_yr * t) for b in src_n.bioc_mg_l) - src_n.biodegradation_capacity_mg_l
    got_n = bs.concentrations(c, src_n, cons, [0.0], [0.0]).instantaneous[0]
    assert close(got_n, max(expect_n, 0.0))
    assert 0 < got_n < got


def test_q2_first_order_mass_retardation():
    c = case_from_workbook(CASES["three_zones_retarded"])
    assert c.retardation > 1.5
    rc, rn = bs.run(c, COMPAT), bs.run(c, bs.Settings(compatible=False))
    ratio_c = rc.grid_masses_kg["first_order"] / rc.grid_masses_kg["no_decay"]
    ratio_n = rn.grid_masses_kg["first_order"] / rn.grid_masses_kg["no_decay"]
    assert ratio_n == pytest.approx(ratio_c * c.retardation, rel=2e-3)


def test_q3_biodegradation_capacity_added_once_when_consistent():
    c = case_from_workbook(CASES["equal_zone_concentrations"])
    sc, sn = bs.source_terms(c, COMPAT), bs.source_terms(c, bs.Settings(compatible=False))
    bc = sc.biodegradation_capacity_mg_l
    assert sc.increments_mg_l[1] == 0
    assert sum(sc.bioc_mg_l) == pytest.approx(sum(sc.increments_mg_l) + 2 * bc)
    assert sum(sn.bioc_mg_l) == pytest.approx(sum(sn.increments_mg_l) + bc)


# --- Analytical checks ---------------------------------------------------------------------------
def test_no_electron_acceptors_instantaneous_equals_no_decay():
    rec = CASES["no_electron_acceptors"]
    r = bs.run(case_from_workbook(rec), COMPAT)
    assert_close(r.arrays["instantaneous"], r.arrays["no_decay"], rel=1e-12)


def test_steady_state_centerline_domenico():
    """t -> infinity, infinite source, alpha_z = 0: C = C0 exp(x/(2ax)(1-sqrt(1+4 lam ax/v))) erf(Y/(4 sqrt(ay x)))
    (limits of erfc -> 2 and Fz -> 2 in BSCR14 p. 8 [72])."""
    from scipy.special import erf
    s = bs.Settings(compatible=False)
    c = bs.Case(vs_ft_yr=100, porosity=0.3, alpha_x_ft=10, alpha_y_ft=1, alpha_z_ft=0, retardation=1.5,
                lambda_per_yr=0.5, length_ft=500, width_ft=200, time_yr=1e5, thickness_ft=10,
                zones=[(0, 0), (0, 0), (40, 5.0)])
    x = np.array([25.0, 100.0, 300.0])
    got = bs.concentrations(c, bs.source_terms(c, s), s, x, np.zeros_like(x)).first_order
    v = c.vs_ft_yr / c.retardation
    expect = 5.0 * np.exp(x / 20 * (1 - np.sqrt(1 + 4 * 0.5 * 10 / v))) * erf(40 / (4 * np.sqrt(x)))
    assert_close(got, expect, rel=1e-12)


def test_source_flux_equals_darcy_flux_times_concentration():
    s = bs.Settings(compatible=False)
    c = case_from_workbook(CASES["single_zone_steady"])
    src = bs.source_terms(c, s)
    assert src.flux_kg_yr == pytest.approx(c.vs_ft_yr * c.porosity * c.thickness_ft * 50 * 10.0 * bs.L_PER_FT3_EXACT / 1e6)
