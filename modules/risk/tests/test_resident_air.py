"""Resident air. References: RSL UG 4.1.3 (PDF pp. 96-101); RSL Summary Table, Nov 2024."""

import csv

import pytest

from modules.risk.resident_air import ResidentAirModule
from modules.risk.rsl_import import DATA

SPECIAL = {"Refractory Ceramic Fibers (units in fibers)"}  # fibre units (RSL UG 5.23)


def sig2(x):
    return float(f"{x:.2g}")


def run(chemicals, **kw):
    m = ResidentAirModule()
    return m.run(m.validate({"chemicals": chemicals, **kw}))


def rows():
    with (DATA / "2024-11" / "toxicity.csv").open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            # Rows without inhalation toxicity values (e.g. lead: NAAQS-based, RSL UG 5.2) are not
            # risk-based and are excluded.
            has_tox = r["rfc_mg_m3"] or r["iur_per_ug_m3"]
            if r["cas"] and r["sl_resident_air_ug_m3"] and has_tox and r["name"] not in SPECIAL:
                yield r


def test_all_chemicals_match_rsl_resident_air():
    rs = list(rows())
    assert len(rs) > 350
    misses = []
    for r in rs:
        o = run([{"cas": r["cas"], "name": r["name"]}])
        got = o.chemicals[0].sl_ug_m3
        if got is None or sig2(got) != sig2(float(r["sl_resident_air_ug_m3"])):
            misses.append((r["name"], r["sl_resident_air_ug_m3"], got))
    assert not misses, misses[:10]


@pytest.mark.parametrize("cas,expected,method", [
    ("71-43-2", 0.36, "standard"),    # benzene
    ("50-32-8", 0.0017, "mutagen"),   # benzo(a)pyrene
    ("75-01-4", 0.17, "vinyl chloride"),
    ("79-01-6", 0.48, "trichloroethylene"),
])
def test_cancer_methods(cas, expected, method):
    c = run([{"cas": cas}]).chemicals[0]
    assert c.cancer_method == method
    assert sig2(c.sl_ug_m3) == expected


def test_forward_risk_and_hq():
    sl = run([{"cas": "71-43-2"}]).chemicals[0]
    o = run([{"cas": "71-43-2", "concentration_ug_m3": sl.sl_ca_ug_m3}])
    assert o.chemicals[0].risk == pytest.approx(1e-6)
    o = run([{"cas": "108-88-3", "concentration_ug_m3": run([{"cas": "108-88-3"}]).chemicals[0].sl_nc_ug_m3}])
    assert o.chemicals[0].hq == pytest.approx(1.0)
    assert o.hazard_index == pytest.approx(1.0)
