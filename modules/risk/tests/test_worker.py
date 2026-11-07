"""Worker soil and air (composite, outdoor, indoor).

References: RSL UG 4.2-4.4 (PDF pp. 101-116); RSL Summary Table, Nov 2024
(Industrial Soil and Industrial Air = composite worker).
"""

import csv

import pytest

from ortam_core.chemicals.rsl_import import DATA
from modules.risk.worker_air import WorkerAirModule
from modules.risk.worker_soil import WorkerSoilModule


def sig2(x):
    return float(f"{x:.2g}")


def run(module, chemicals, **kw):
    return module.run(module.validate({"chemicals": chemicals, **kw}))


def table():
    with (DATA / "2024-11" / "toxicity.csv").open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def has_tox(r, keys):
    return any(r[k] for k in keys)


def test_composite_soil_matches_industrial_soil_table():
    m = WorkerSoilModule()
    rows = [r for r in table() if r["cas"] and r["sl_industrial_soil_mg_kg"] and r["volatile"] != "V"
            and has_tox(r, ("sfo_per_mg_kg_day", "iur_per_ug_m3", "rfd_mg_kg_day", "rfc_mg_m3"))]
    assert len(rows) > 550
    misses = [(r["name"], r["sl_industrial_soil_mg_kg"]) for r in rows
              if sig2(run(m, [{"cas": r["cas"], "name": r["name"]}]).chemicals[0].sl_mg_kg)
              != sig2(float(r["sl_industrial_soil_mg_kg"]))]
    assert not misses, misses[:10]


def test_composite_air_matches_industrial_air_table():
    m = WorkerAirModule()
    rows = [r for r in table() if r["cas"] and r["sl_industrial_air_ug_m3"]
            and has_tox(r, ("iur_per_ug_m3", "rfc_mg_m3")) and not r["name"].startswith("Refractory")]
    assert len(rows) > 350
    misses = [(r["name"], r["sl_industrial_air_ug_m3"]) for r in rows
              if sig2(run(m, [{"cas": r["cas"], "name": r["name"]}]).chemicals[0].sl_ug_m3)
              != sig2(float(r["sl_industrial_air_ug_m3"]))]
    assert not misses, misses[:10]


def test_outdoor_worker_by_hand():
    # RSL UG 4.3.1.1 (PDF p. 107): copper, RfD 0.04, ingestion only, EF 225, ED 25, IRS 100, BW 80
    expected = 1.0 * 365 * 25 * 80 / ((1 / 0.04) * 1e-6 * 225 * 25 * 100)
    c = run(WorkerSoilModule(), [{"cas": "7440-50-8"}], receptor="outdoor").chemicals[0]
    assert c.sl_nc_mg_kg.ingestion == pytest.approx(expected)


def test_indoor_worker_has_no_dermal_route():
    # RSL UG 4.4.1 (PDF p. 112): no dermal contact; IRS 50 mg/day
    c = run(WorkerSoilModule(), [{"cas": "30560-19-1"}], receptor="indoor").chemicals[0]
    assert c.sl_nc_mg_kg.dermal is None
    expected = 1.0 * 365 * 25 * 80 / ((1 / 3e-4) * 1e-6 * 250 * 25 * 50)
    assert c.sl_nc_mg_kg.ingestion == pytest.approx(expected)


def test_outdoor_air_by_hand():
    # RSL UG 4.3.2.2 (PDF p. 111): benzene IUR 7.8E-06
    expected = 1e-6 * 365 * 70 / (7.8e-6 * 225 * 25 * 8 / 24)
    c = run(WorkerAirModule(), [{"cas": "71-43-2"}], receptor="outdoor").chemicals[0]
    assert c.sl_ca_ug_m3 == pytest.approx(expected)


def test_volatile_soil_needs_vf():
    with pytest.raises(ValueError, match="volatil"):
        run(WorkerSoilModule(), [{"cas": "71-43-2"}])


def test_forward_risk():
    m = WorkerSoilModule()
    sl = run(m, [{"cas": "60-57-1"}]).chemicals[0].sl_ca_mg_kg.total
    assert run(m, [{"cas": "60-57-1", "concentration_mg_kg": sl}]).total_risk == pytest.approx(1e-6)
