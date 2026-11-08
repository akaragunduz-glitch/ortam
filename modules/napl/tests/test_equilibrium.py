"""napl.equilibrium against published values.

References:
  HSSM Vol. 1 (docs/references/hssm_vol1_users_guide.pdf), Appendix 6, RAOULT utility output:
    Table 105, p. 185 [200]: default gasoline (Baehr and Corapcioglu, 1987): density 0.7168 g/cm3,
      average MW 104.0458 g/mol; benzene Ko 311.6757 (Eq. 44) and 307.0647 (Eq. 45).
    Table 107, p. 187 [202]: benzene solubility changed to 1780 mg/L: 306.4227 and 301.8895.
  HSSM Vol. 1 Table 100, p. 174 [189]: Ko calculated from the same mixture, in parentheses:
    benzene (312), toluene (1202), o-xylene (4440).
  BIOSCREEN 1.4 manual (docs/references/bioscreen_users_manual_1996.pdf), pp. 4-5 [68-70]:
    effective solubility = mole fraction x pure-phase solubility for fresh gasoline, two
    weathered gasolines and JP-4 (printed to 1 mg/L; mole fractions printed to 2-4 digits).
  RSL Summary Table Nov 2024: Csat of 206 organics (SSG Eq. 4-9, RSL default soil); the
    single-component limit of the soil basis must reproduce it.
  Mass balance: a NAPL of known composition and saturation is placed in soil with HSSM2
    Eq. 25-26 (written out in this file); the module must recover saturation and composition.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest
import yaml

from modules.napl.equilibrium import NaplEquilibriumModule
from modules.napl import equations as eq
from modules.registry import get_module
from ortam_core.chemicals.database import lookup

ROOT = Path(__file__).resolve().parents[3]
EXAMPLES = Path(__file__).resolve().parents[1] / "examples"

# HSSM Vol. 1 Table 106, p. 186 [201]: default RAOULT.DAT (coj g/cm3, S mg/L, MW g/mol)
RAOULT_DAT = [("benzene", 0.0082, 1750.0, 78.0), ("toluene", 0.0426, 535.0, 92.0),
              ("xylenes", 0.0718, 167.0, 106.0), ("1-hexene", 0.0159, 0.0, 84.0),
              ("cyclohexane", 0.0021, 0.0, 84.0), ("n-hexane", 0.0204, 0.0, 86.0),
              ("other_aromatics", 0.0740, 0.0, 106.0), ("other_paraffins", 0.3367, 0.0, 97.2),
              ("heavy_ends", 0.1451, 0.0, 128.0)]


def run(**kw):
    m = NaplEquilibriumModule()
    return m.run(m.validate(kw))


def raoult_components(benzene_s=1750.0):
    return [{"name": n, "napl_g_cm3": c, "solubility_mg_l": benzene_s if n == "benzene" else s, "mw_g_mol": t}
            for n, c, s, t in RAOULT_DAT]


# ------------------------------------------------------------ HSSM RAOULT
def test_raoult_table_105():
    o = run(basis="napl", chemical_table="none", components=raoult_components())
    assert o.napl_density_g_cm3 == pytest.approx(0.7168, abs=5e-5)
    assert o.average_mw_g_mol == pytest.approx(104.0458, abs=5e-5)
    b = o.components[0]
    assert b.ko_composition == pytest.approx(311.6757, abs=5e-5)
    assert b.ko_average_mw == pytest.approx(307.0647, abs=5e-5)


def test_raoult_table_107():
    b = run(basis="napl", chemical_table="none", components=raoult_components(1780.0)).components[0]
    assert b.ko_composition == pytest.approx(306.4227, abs=5e-5)
    assert b.ko_average_mw == pytest.approx(301.8895, abs=5e-5)


@pytest.mark.parametrize("i,printed", [(0, 312.0), (1, 1202.0), (2, 4440.0)])
def test_hssm_table_100_calculated_ko(i, printed):
    ko = run(basis="napl", chemical_table="none", components=raoult_components()).components[i].ko_composition
    assert float(f"{ko:.3g}") == pytest.approx(float(f"{printed:.3g}"))


def test_effective_solubility_equals_co_over_ko():
    o = run(basis="napl", chemical_table="none", components=raoult_components())
    for c in o.components[:3]:
        assert c.effective_solubility_mg_l == pytest.approx(c.napl_mg_l / c.ko_composition, rel=1e-12)  # Eq. 43


def test_mass_fractions_equal_concentrations():
    rho = sum(c for _, c, _, _ in RAOULT_DAT)
    comps = [{"name": n, "napl_mass_fraction": c / rho, "solubility_mg_l": s, "mw_g_mol": t}
             for n, c, s, t in RAOULT_DAT]
    b = run(basis="napl", chemical_table="none", napl_density_g_cm3=rho, components=comps).components[0]
    assert b.ko_composition == pytest.approx(311.6757, abs=5e-5)


# ------------------------------------------------------------ BIOSCREEN effective solubility
BIOSCREEN = [  # fuel, constituent, mole fraction, S (mg/L), effective solubility (mg/L)
    ("fresh", "benzene", 0.0093, 1780, 17), ("fresh", "toluene", 0.0568, 515, 29),
    ("fresh", "ethylbenzene", 0.0, 152, 0), ("fresh", "xylenes", 0.0858, 198, 17),
    ("weathered1", "benzene", 0.0137, 1780, 24), ("weathered1", "toluene", 0.1216, 515, 63),
    ("weathered1", "ethylbenzene", 0.0, 152, 0), ("weathered1", "xylenes", 0.1247, 198, 25),
    ("weathered2", "benzene", 0.003, 1780, 5), ("weathered2", "toluene", 0.043, 515, 22),
    ("weathered2", "ethylbenzene", 0.014, 152, 2), ("weathered2", "xylenes", 0.084, 198, 15),
    ("jp4", "benzene", 0.023, 1780, 42), ("jp4", "toluene", 0.053, 515, 27),
    ("jp4", "ethylbenzene", 0.013, 152, 2), ("jp4", "xylenes", 0.080, 198, 16),
]


BIOSCREEN_KNOWN = {("weathered2", "xylenes")}


@pytest.mark.parametrize("fuel,name,x,s,printed", BIOSCREEN, ids=[f"{r[0]}-{r[1]}" for r in BIOSCREEN])
def test_bioscreen_effective_solubility(fuel, name, x, s, printed):
    # tolerance: rounding of the printed result (0.5 mg/L) + half a unit in the last printed digit of x, times S
    digits = len(f"{x}".split(".")[1]) if x else 4
    tol = 0.5 + 0.5 * 10 ** -digits * s
    if (fuel, name) in BIOSCREEN_KNOWN:
        # BIOSCREEN prints x = 0.084 but 15 mg/L (= 0.076 x 198); the 15 is consistent with the
        # printed BTEX total of 44 mg/L. Inconsistency in the source table, documented in
        # docs/methods/napl.md.
        assert abs(eq.effective_solubility(x, s) - printed) > tol
        return
    assert abs(eq.effective_solubility(x, s) - printed) <= tol


def test_bioscreen_btex_totals():
    # printed BTEX totals = sum of the printed effective solubilities (pp. 4-5 [68-70])
    for fuel, total in (("fresh", 63), ("weathered1", 112), ("weathered2", 44), ("jp4", 87)):
        assert sum(r[4] for r in BIOSCREEN if r[0] == fuel) == pytest.approx(total, abs=1)


# ------------------------------------------------------------ Csat (RSL)
def _rsl_csat_rows():
    tox = csv.DictReader((ROOT / "data/reference/rsl/2024-11/toxicity.csv").open(encoding="utf-8"))
    out = []
    for r in tox:
        if not r["csat_mg_kg"]:
            continue
        c = lookup("rsl-2024-11", cas=r["cas"], name=r["name"])
        if c.koc_l_kg and c.solubility_mg_l and c.henry_dimensionless is not None:
            out.append((r["cas"], r["name"], float(r["csat_mg_kg"])))
    return out


CSAT = _rsl_csat_rows()


def test_csat_reference_count():
    assert len(CSAT) == 206


@pytest.mark.parametrize("cas,name,ref", CSAT, ids=[r[1] for r in CSAT])
def test_single_component_limit_is_rsl_csat(cas, name, ref):
    o = run(basis="soil", components=[{"cas": cas, "name": name, "soil_mg_kg": ref}])
    c = o.components[0]
    assert float(f"{c.csat_mg_kg:.3g}") == pytest.approx(ref, rel=1e-9) or c.csat_mg_kg == pytest.approx(ref, rel=5e-3)
    below = run(basis="soil", components=[{"cas": cas, "name": name, "soil_mg_kg": c.csat_mg_kg * 0.999}])
    above = run(basis="soil", components=[{"cas": cas, "name": name, "soil_mg_kg": c.csat_mg_kg * 1.001}])
    assert below.napl_present is False and above.napl_present is True
    assert above.components[0].mole_fraction == pytest.approx(1.0)
    assert above.components[0].water_mg_l == pytest.approx(lookup("rsl-2024-11", cas=cas, name=name).solubility_mg_l)


# ------------------------------------------------------------ mass balance (HSSM2 Eq. 25-26)
BTX = {"benzene": ("71-43-2", None), "toluene": ("108-88-3", None), "xylenes": ("1330-20-7", None)}


def _forward(theta_o: float):
    """Place the RAOULT gasoline at NAPL content theta_o in RSL default soil; return soil mg/kg."""
    rb, rs, tw, foc = 1.5, 2.65, 0.15, 0.006
    n, ta = eq.porosities(rb, rs, tw)
    co = [c for _, c, _, _ in RAOULT_DAT]
    mw = [t for _, _, _, t in RAOULT_DAT]
    x = eq.mole_fractions(co, mw)
    comps = []
    for (name, c_o, s, t), xi in zip(RAOULT_DAT, x):
        d = {"name": name, "mw_g_mol": t, "solubility_mg_l": s}
        kd = h = 0.0
        if name in BTX:
            chem = lookup("rsl-2024-11", cas=BTX[name][0])
            d["cas"] = BTX[name][0]
            kd, h = chem.koc_l_kg * foc, chem.henry_dimensionless
        cw = xi * s                                                   # BIOSCREEN p. 3
        total = cw * (tw + rb * kd + ta * h) + theta_o * c_o * 1e6    # HSSM2 Eq. 25: phi*B*cw, mg/L soil
        d["soil_mg_kg"] = total / rb
        comps.append(d)
    return comps, n


@pytest.mark.parametrize("theta_o", [1e-4, 0.01, 0.05])
def test_mass_balance_round_trip(theta_o):
    comps, n = _forward(theta_o)
    o = run(basis="soil", napl_density_g_cm3=0.7168, components=comps)
    assert o.napl_present
    assert o.napl_volume_fraction == pytest.approx(theta_o, rel=1e-9)
    assert o.napl_saturation == pytest.approx(theta_o / n, rel=1e-9)
    co = [c for _, c, _, _ in RAOULT_DAT]
    for r, c_o in zip(o.components, co):
        assert r.napl_mg_l == pytest.approx(c_o * 1e6, rel=1e-9)
    assert o.components[0].ko_composition == pytest.approx(311.6757, abs=5e-4)
    # phases add up to the total
    for r in o.components:
        total = (r.water_mg_l * 0.15 + r.sorbed_mg_kg * 1.5 + r.soil_gas_mg_l * (n - 0.15)
                 + r.napl_mg_l * o.napl_volume_fraction) / 1.5
        assert total == pytest.approx(r.soil_mg_kg, rel=1e-9)


def test_no_napl_below_saturation():
    o = run(basis="soil", components=[{"cas": "71-43-2", "soil_mg_kg": 10.0},
                                      {"cas": "108-88-3", "soil_mg_kg": 20.0}])
    assert o.napl_present is False and o.saturation_index < 1
    b = o.components[0]
    assert b.in_napl_pct == 0.0 and b.mole_fraction is None
    assert b.water_mg_l * 0.15 + b.sorbed_mg_kg * 1.5 + b.soil_gas_mg_l * o.soil["theta_a"] == \
        pytest.approx(10.0 * 1.5)


def test_mixture_forms_napl_below_single_csat():
    # each below its own Csat, the sum Ct/Csat above 1: NAPL forms (Raoult's law)
    b, t = lookup("rsl-2024-11", cas="71-43-2"), lookup("rsl-2024-11", cas="108-88-3")
    o = run(basis="soil", components=[{"cas": "71-43-2", "soil_mg_kg": 0.6 * b.csat_mg_kg},
                                      {"cas": "108-88-3", "soil_mg_kg": 0.6 * t.csat_mg_kg}])
    assert o.napl_present and o.saturation_index == pytest.approx(1.2, rel=5e-3)
    assert sum(c.mole_fraction for c in o.components) == pytest.approx(1.0)


def test_insoluble_component_always_forms_napl():
    o = run(basis="soil", chemical_table="none", napl_density_g_cm3=0.8, residual_saturation=0.05,
            components=[{"name": "heavy_ends", "mw_g_mol": 128.0, "solubility_mg_l": 0.0, "soil_mg_kg": 100.0}])
    n = 1 - 1.5 / 2.65
    assert o.napl_present and o.saturation_index is None
    assert o.napl_volume_fraction == pytest.approx(100.0 * 1.5 / 0.8e6)
    assert o.napl_saturation == pytest.approx(100.0 * 1.5 / 0.8e6 / n)
    assert o.potentially_mobile is False


# ------------------------------------------------------------ indicators and input checks
def test_lnapl_indicators():
    o = run(basis="napl", chemical_table="none", tph_groundwater_mg_l=12.0, tph_soil_mg_kg=300.0, fuel="diesel",
            components=[{"name": "benzene", "cas": "71-43-2", "napl_g_cm3": 0.0082, "solubility_mg_l": 1750.0,
                         "mw_g_mol": 78.0, "groundwater_mg_l": 6.0}] + raoult_components()[1:])
    ind = {i.indicator: i for i in o.indicators}
    assert ind["benzene in groundwater"].exceeded and ind["benzene in groundwater"].threshold == 5.0
    assert not ind["TPH in groundwater"].exceeded and ind["TPH in groundwater"].threshold == 30.0
    assert ind["TPH in soil (diesel)"].exceeded and ind["TPH in soil (diesel)"].threshold == 250.0
    seff = o.components[0].effective_solubility_mg_l
    r = ind["benzene in groundwater / effective solubility"]
    assert r.value == pytest.approx(6.0 / seff) and r.exceeded == (6.0 / seff > 0.2)


def test_indicator_thresholds_as_printed():
    d = yaml.safe_load((ROOT / "data/reference/napl/lnapl_indicators.yaml").read_text())
    assert (d["groundwater_mg_l"]["benzene"]["value"], d["groundwater_mg_l"]["tph"]["value"]) == (5.0, 30.0)
    assert d["soil_mg_kg"]["benzene"]["value"] == 10.0
    assert [d["soil_mg_kg"]["tph"][k]["value"] for k in ("unweathered-gasoline", "weathered-gasoline", "diesel")] \
        == [100.0, 250.0, 250.0]
    assert d["effective_solubility_fraction"]["value"] == 0.2


def test_non_detect_cannot_be_entered_as_zero():
    m = NaplEquilibriumModule()
    with pytest.raises(Exception):
        m.validate({"basis": "soil", "components": [{"cas": "71-43-2", "soil_mg_kg": 0.0}]})


def test_input_checks():
    m = NaplEquilibriumModule()
    with pytest.raises(Exception):   # mass fractions without density
        m.validate({"basis": "napl", "components": [{"cas": "71-43-2", "napl_mass_fraction": 0.01}]})
    with pytest.raises(Exception):   # mixed composition inputs
        m.validate({"basis": "napl", "components": [{"cas": "71-43-2", "napl_mass_fraction": 0.01},
                                                    {"cas": "108-88-3", "napl_g_cm3": 0.04}]})
    with pytest.raises(Exception):   # pseudo-component without properties
        m.run(m.validate({"basis": "napl", "components": [{"name": "paraffins", "napl_g_cm3": 0.3}]}))


@pytest.mark.parametrize("example", sorted(EXAMPLES.glob("*.yaml")), ids=lambda p: p.name)
def test_examples_run(example):
    m = get_module("napl.equilibrium")
    o = m.run(m.validate(yaml.safe_load(example.read_text())))
    assert m.summary(o)
