"""Resident soil screening levels and risk.

References:
  - RSL User's Guide (2026-03), Section 4.1.1 (resident soil equations).
  - RSL Summary Table, November 2024 (TR = 1E-06, THQ = 1): reference values.
  - RSL User's Guide Section 2.6 and RAGS Part A (1989) Chapter 8: risk and HQ.
RSL tables report two significant figures, so comparisons use 2 s.f.
"""

import csv
import math
import subprocess
import sys
from pathlib import Path

import pytest

from modules.risk.resident_soil import ResidentSoilModule, age_adjusted_factors, load_exposure_defaults
from modules.risk.rsl_import import DATA

ROOT = Path(__file__).resolve().parents[3]


def sig2(x: float) -> float:
    return float(f"{x:.2g}")


def run(chemicals, **kw):
    m = ResidentSoilModule()
    return m.run(m.validate({"chemicals": chemicals, **kw}))


def test_age_adjusted_factors():
    # RSL UG 4.1.1.3 (p. 40-41): IFSres-adj = 36750 mg/kg, DFSres-adj = 103390 mg/kg
    ifs, dfs = age_adjusted_factors(load_exposure_defaults())
    assert ifs == pytest.approx(36750)
    assert dfs == pytest.approx(103390)


def test_acephate_route_values_by_hand():
    # Acephate: RfD = 3E-04 mg/kg-day, GIABS = 1, ABSd = 0.1, no RfC (RSL Nov 2024)
    # Child ingestion, RSL UG 4.1.1.1:
    #   SL = THQ*AT*BW / ((RBA/RfD)*1e-6*EF*ED*IRS) = 1*2190*15 / ((1/3e-4)*1e-6*350*6*200)
    ing = 2190 * 15 / ((1 / 3e-4) * 1e-6 * 350 * 6 * 200)
    # Child dermal: SL = THQ*AT*BW / ((1/(RfD*GIABS))*1e-6*EF*ED*SA*AF*ABSd)
    derm = 2190 * 15 / ((1 / 3e-4) * 1e-6 * 350 * 6 * 2373 * 0.2 * 0.1)
    c = run([{"cas": "30560-19-1"}]).chemicals[0]
    assert c.sl_nc_child_mg_kg.ingestion == pytest.approx(ing)
    assert c.sl_nc_child_mg_kg.dermal == pytest.approx(derm)
    assert c.sl_nc_child_mg_kg.inhalation is None
    assert c.sl_nc_child_mg_kg.total == pytest.approx(1 / (1 / ing + 1 / derm))
    assert sig2(c.sl_mg_kg) == 19  # RSL Summary Table, Nov 2024


# Ten chemicals, each exercising a different combination of routes and effects.
TEN_CHEMICALS = [
    # cas,         name,                  RSL (mg/kg), basis
    ("7440-50-8", "Copper",               3100,   "n"),  # ingestion only
    ("30560-19-1", "Acephate",            19,     "n"),  # ingestion + dermal
    ("7440-36-0", "Antimony (metallic)",  31,     "n"),  # GIABS 0.15, PEF inhalation
    ("7440-39-3", "Barium",               15000,  "n"),  # GIABS 0.07, PEF inhalation
    ("7440-41-7", "Beryllium and compounds", 160, "n"),  # nc governs over IUR
    ("7440-02-0", "Nickel Soluble Salts", 1400,   "n"),  # nc governs over IUR
    ("7440-48-4", "Cobalt",               23,     "n"),
    ("1912-24-9", "Atrazine",             2.4,    "c"),  # cancer: ingestion + dermal
    ("60-57-1",   "Dieldrin",             0.034,  "c"),  # cancer: all three routes
    ("87-86-5",   "Pentachlorophenol",    1.0,    "c"),  # cancer: ABSd 0.25, IUR
]


@pytest.mark.parametrize("cas,name,expected,basis", TEN_CHEMICALS)
def test_ten_chemicals_match_rsl(cas, name, expected, basis):
    c = run([{"cas": cas}]).chemicals[0]
    assert c.name == name
    assert sig2(c.sl_mg_kg) == expected
    assert c.sl_basis == basis


def _regression_rows():
    with (DATA / "2024-11" / "toxicity.csv").open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            has_tox = any(r[k] for k in ("sfo_per_mg_kg_day", "iur_per_ug_m3", "rfd_mg_kg_day", "rfc_mg_m3"))
            if (r["cas"] and r["sl_resident_soil_mg_kg"] and has_tox
                    and r["volatile"] != "V" and r["mutagen"] != "M"):
                yield r


def test_all_nonvolatile_nonmutagenic_chemicals_match_rsl():
    rows = list(_regression_rows())
    assert len(rows) > 500
    m = ResidentSoilModule()
    misses = []
    for r in rows:
        o = m.run(m.validate({"chemicals": [{"cas": r["cas"], "name": r["name"]}]}))
        got, ref = sig2(o.chemicals[0].sl_mg_kg), float(r["sl_resident_soil_mg_kg"])
        if got != sig2(ref):
            misses.append((r["name"], ref, got))
    assert not misses, misses[:10]


def test_forward_risk_at_screening_level_equals_targets():
    # RSL UG 2.6.1: Risk = C*TR/RSL, HQ = C*THQ/RSL
    sl = run([{"cas": "60-57-1"}]).chemicals[0]
    o = run([{"cas": "60-57-1", "concentration_mg_kg": sl.sl_ca_mg_kg.total}])
    assert o.chemicals[0].risk.total == pytest.approx(1e-6)
    assert o.total_risk == pytest.approx(1e-6)
    o = run([{"cas": "7440-50-8", "concentration_mg_kg": sl_of("7440-50-8")}])
    assert o.chemicals[0].hq_child.total == pytest.approx(1.0)
    assert o.hazard_index_child == pytest.approx(1.0)


def sl_of(cas):
    return run([{"cas": cas}]).chemicals[0].sl_nc_child_mg_kg.total


def test_one_hit_rule_above_0_01():
    # RSL UG 2.6.2; RAGS A (1989) p. 8-11: Risk = 1 - exp(-CDI x SF) when linear risk > 0.01
    sl = run([{"cas": "60-57-1"}]).chemicals[0].sl_ca_mg_kg.total
    o = run([{"cas": "60-57-1", "concentration_mg_kg": sl * 5e4}])  # linear risk = 0.05
    assert o.chemicals[0].risk.total == pytest.approx(1 - math.exp(-0.05))
    assert o.total_risk == pytest.approx(1 - math.exp(-0.05))


def test_hazard_index_and_total_risk_sum_over_chemicals():
    # RAGS A (1989) p. 8-12 (Risk_T = sum Risk_i) and p. 8-13 (HI = sum HQ)
    a = sl_of("7440-50-8")
    b = sl_of("7440-39-3")
    o = run([{"cas": "7440-50-8", "concentration_mg_kg": 0.5 * a},
             {"cas": "7440-39-3", "concentration_mg_kg": 0.25 * b}])
    assert o.hazard_index_child == pytest.approx(0.75)


def test_override_is_flagged():
    o = run([{"cas": "7440-50-8", "rfd_mg_kg_day": 0.08}])
    c = o.chemicals[0]
    assert "rfd_mg_kg_day" in c.overridden
    assert sig2(c.sl_mg_kg) == 6300  # twice the RfD -> twice the SL (ingestion only)
    assert any("overrid" in w.lower() for w in o.warnings)


def test_exposure_override_is_flagged():
    o = run([{"cas": "7440-50-8"}], exposure={"irs_child_mg_d": 100})
    assert "irs_child_mg_d" in o.exposure_overridden
    assert sig2(o.chemicals[0].sl_mg_kg) == 6300


def test_standalone_without_table():
    o = run([{"name": "Test chemical", "rfd_mg_kg_day": 0.04, "giabs": 1}], toxicity_table="none")
    assert sig2(o.chemicals[0].sl_mg_kg) == 3100  # same as copper


def test_arsenic_uses_rba_0_6():
    c = run([{"cas": "7440-38-2"}]).chemicals[0]
    assert c.rba == 0.6
    assert sig2(c.sl_mg_kg) == 0.68  # RSL Summary Table, Nov 2024


def test_volatile_without_vf_rejected():
    with pytest.raises(ValueError, match="volatil"):
        run([{"cas": "71-43-2"}])  # benzene


def test_volatile_with_user_vf_accepted():
    o = run([{"cas": "71-43-2", "vf_m3_kg": 3000.0}])
    assert o.chemicals[0].sl_ca_mg_kg.inhalation is not None
    assert "vf_m3_kg" in o.chemicals[0].overridden


def test_mutagen_rejected():
    with pytest.raises(ValueError, match="mutagen"):
        run([{"cas": "50-32-8"}])  # benzo(a)pyrene


def test_ambiguous_cas_requires_name():
    with pytest.raises(ValueError, match="name"):
        run([{"cas": "7440-43-9"}])  # Cadmium (Diet) and Cadmium (Water)
    c = run([{"cas": "7440-43-9", "name": "Cadmium (Diet)"}]).chemicals[0]
    assert sig2(c.sl_mg_kg) == 7.1


def test_ceiling_limit_flagged_not_applied():
    # RSL UG 5.13: SLs above 1E+05 mg/kg are flagged, not replaced
    o = run([{"cas": "57-55-6"}])  # propylene glycol, RSL 1.3E+06 mg/kg (key nm)
    assert sig2(o.chemicals[0].sl_mg_kg) == 1.3e6
    assert any("ceiling" in w for w in o.warnings)


def test_cli_example_runs():
    r = subprocess.run([sys.executable, "-m", "ortam_core.cli", "run", "risk.resident-soil",
                        str(ROOT / "modules/risk/examples/resident_soil.yaml")],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr
    assert "Hazard index" in r.stdout
