"""vadose.vleach against EPA's VLEACH 2.2a executable.

Reference: data/reference/vleach/2.2a/<case>/ holds the input (.inp) and the output
files of EPA's executable (VleachForGUI.exe from the VLEACH 2.2a setup):
  epa_sample           EPA's own sample problem shipped with the setup (TCE, 2 polygons)
  case1..case5         generated with modules/vadose/tests/generate_vleach_reference.py
                       (Wine), covering recharge concentration, atmosphere and groundwater
                       gas boundaries, impermeable boundaries, three polygons, diffusion
                       only and advection only.
Compared quantities: phase masses and boundary fluxes at every printout (.OUT,
5 significant figures), liquid/gas/sorbed profiles (.PRF, 5 figures) and the
groundwater mass rate of every time step (GWIMP.DAT, 7 figures).
Tolerance 5e-4 (relative): the executable computes in single precision
(about 7 digits); in diffusion-dominated cases the mass decays by orders of
magnitude, which amplifies the round-off (case 4: 3e-4; other cases <= 1e-4).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from modules.vadose.vleach import VleachModule
from modules.vadose.vleach_inp import read_inp

REF = Path(__file__).resolve().parents[3] / "data" / "reference" / "vleach" / "2.2a"
CASES = sorted(p.name for p in REF.iterdir() if p.is_dir())
FT2 = 0.3048 ** 2                     # m2 per ft2
FT3_L = 0.3048 ** 3 * 1000.0          # L per ft3
TOL = 5e-4
NUM = r"([-+]?\d*\.?\d+(?:E[-+]?\d+)?)"


def _close(a, b, scale):
    """Relative agreement; values below 1e-6 of the largest value of the series are compared absolutely."""
    return abs(a - b) <= TOL * max(abs(b), 1e-6 * scale)


def _out_blocks(text):
    """Per polygon, per printout: total, gas, liquid, sorbed, then the 'since last printout' fluxes."""
    states = re.findall(r"Polygon\s+(\d+)\s*\nAt time =\s*" + NUM + r", total mass in vadose zone =\s*" + NUM
                        + r".*?\nMass in gas phase\s*=\s*" + NUM + r".*?\nMass in liquid phase\s*=\s*" + NUM
                        + r".*?\nMass sorbed\s*=\s*" + NUM, text)
    fluxes = re.findall(r"Since last printout at time =\s*" + NUM + r".*?\n.*?Change in Total Mass =\s*" + NUM
                        + r".*?\n.*?Advection in from atmosphere =\s*" + NUM
                        + r".*?\n.*?Advection in from water table =\s*" + NUM
                        + r".*?\n.*?Diffusion in from atmosphere =\s*" + NUM
                        + r".*?\n.*?Diffusion in from water table =\s*" + NUM, text)
    return states, fluxes


def _profiles(text):
    """List (polygon index, time) -> rows (Cgas, Cliq, Csol)."""
    out, poly, t, last_t = {}, -1, None, None
    for line in text.splitlines():
        m = re.match(r"\s*Time:\s*" + NUM, line)
        if m:
            t = float(m.group(1))
            if last_t is None or t <= last_t:
                poly += 1
            last_t = t
            out[(poly, t)] = []
            continue
        m = re.match(r"\s*\d+\s+" + NUM + r"\s+" + NUM + r"\s+" + NUM + r"\s*$", line)
        if m and t is not None:
            out[(poly, t)].append(tuple(float(x) for x in m.groups()))
    return out


@pytest.fixture(scope="module", params=CASES)
def case(request):
    d = REF / request.param
    m = VleachModule()
    o = m.run(m.validate(read_inp(d / "input.inp")))
    return request.param, d, o


def test_reference_cases_present():
    assert "epa_sample" in CASES and len(CASES) >= 6


def test_mass_and_fluxes(case):
    name, d, o = case
    states, fluxes = _out_blocks((d / "output.out").read_text())
    mine = [mb for p in o.polygons for mb in p.mass_balance]
    assert len(states) == len(mine), name
    scale = max(float(s[2]) for s in states)
    for s, mb in zip(states, mine):
        assert float(s[1]) == pytest.approx(mb.time_yr)
        for ref, val in zip(s[2:], (mb.total_g_m2, mb.gas_g_m2, mb.liquid_g_m2, mb.sorbed_g_m2)):
            assert _close(val * FT2, float(ref), scale), (name, mb.time_yr, ref, val * FT2)
    since_last = fluxes                      # the regex matches only the 'since last printout' blocks
    mine_f = [mb for p in o.polygons for mb in p.mass_balance[1:]]
    assert len(since_last) == len(mine_f)
    for f, mb in zip(since_last, mine_f):
        for ref, val in zip(f[2:], (mb.advection_in_atmosphere_g_m2, mb.advection_in_water_table_g_m2,
                                    mb.diffusion_in_atmosphere_g_m2, mb.diffusion_in_water_table_g_m2)):
            assert _close(val * FT2, float(ref), scale), (name, mb.time_yr, ref, val * FT2)


def test_profiles(case):
    name, d, o = case
    ref = _profiles((d / "profiles.prf").read_text())
    n = 0
    for k, p in enumerate(o.polygons):
        for pr in p.profiles:
            rows = ref.get((k, round(pr.time_yr, 3)))
            assert rows is not None, (name, k, pr.time_yr)
            scale = max(r[1] for r in rows) or 1.0
            for r, cg, cl, cs in zip(rows, pr.gas_mg_l, pr.liquid_mg_l, pr.sorbed_mg_kg):
                assert _close(cl / 1000.0 * FT3_L, r[1], scale), (name, pr.time_yr)      # g/ft3
                assert _close(cg / 1000.0 * FT3_L, r[0], scale * 10), (name, pr.time_yr)
                assert _close(cs * 1e-6, r[2], scale), (name, pr.time_yr)                # g/g
                n += 1
    assert n > 0


def test_groundwater_rate_every_step(case):
    name, d, o = case
    ref = [tuple(map(float, line.split())) for line in (d / "gwimp.dat").read_text().splitlines() if line.strip()]
    mine = [x for p in o.polygons for x in p.groundwater_series]
    assert len(ref) == len(mine), name
    scale = max(abs(r[1]) for r in ref) or 1.0
    for (t, r), (tm, v) in zip(ref, mine):
        assert t == pytest.approx(tm, rel=1e-6)
        assert _close(v, r, scale), (name, t, r, v)


def test_total_impact_table(case):
    name, d, o = case
    text = (d / "output.out").read_text()
    rows = re.findall(r"^\s*" + NUM + r"\s+" + NUM + r"\s+" + NUM + r"\s*$",
                      text.split("TOTAL GROUNDWATER IMPACT")[1], re.M)
    assert len(rows) == len(o.total_impact)
    scale = max(float(r[2]) for r in rows) or 1.0
    for r, t in zip(rows, o.total_impact):
        assert _close(t["rate_g_yr"], float(r[1]), scale) and _close(t["cumulative_g"], float(r[2]), scale), name


def test_mass_balance_closes():
    m = VleachModule()
    o = m.run(m.validate(read_inp(REF / "case1_boundaries" / "input.inp")))
    for mb in o.polygons[0].mass_balance:
        assert abs(mb.discrepancy_g_m2) < 1e-12


def test_input_checks():
    m = VleachModule()
    inp = read_inp(REF / "epa_sample" / "input.inp")
    inp["polygons"][0]["water_content"] = 0.5          # > porosity
    with pytest.raises(Exception):
        m.validate(inp)
    inp = read_inp(REF / "epa_sample" / "input.inp")
    inp["output_interval_yr"] = 150                   # not a multiple of the 100-yr step
    with pytest.raises(Exception):
        m.validate(inp)


def test_courant_warning():
    m = VleachModule()
    o = m.run(m.validate(read_inp(REF / "epa_sample" / "input.inp")))
    assert any("Courant" in w for w in o.warnings)    # EPA sample: q*dt/dz = 100*0.3048/0.3048 = 100
