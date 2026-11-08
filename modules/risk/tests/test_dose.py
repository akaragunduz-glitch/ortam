"""Doses (CDI, DAD, EC): modules/risk/equations/dose.py and the dose output of the risk tools.

1. External reference: RAGS Part E (2004) Exhibit B-3 (209 organic chemicals in water;
   modules/risk/tests/data/rags_e_exhibit_b3.csv). The exhibit prints its inputs
   (Kp, B, tau, t*, FA) and results (DA_event, DAD) to two significant figures, so the
   recomputed values can differ by the rounding of the inputs: tolerance 10 %,
   median deviation below 2 % (two-digit rounding alone gives a median of about 1-2 %).
2. Consistency: for every RSL chemical and every risk tool, the risk and HQ obtained
   from the doses (Risk = dose*toxicity, HQ = dose/toxicity; RAGS A pp. 8-6, 8-11;
   RAGS E Eq. 4.2-4.3; RAGS F Eq. 11-12) equal the risk and HQ of the tool, which
   are verified against the RSL tables. This checks the dose equations by an
   independent route (doses are built from intake equations, the tools from SLs).
"""

from __future__ import annotations

import csv
import statistics
from pathlib import Path

import pytest

from modules.risk.common import one_hit
from modules.risk.construction_soil import ConstructionSoilModule
from modules.risk.equations import dose as dq
from modules.risk.fish import FishModule
from modules.risk.recreator_soil import RecreatorSoilModule
from modules.risk.recreator_surface_water import RecreatorSurfaceWaterModule
from modules.risk.resident_air import ResidentAirModule
from modules.risk.resident_soil import ResidentSoilModule
from modules.risk.resident_tapwater import ResidentTapwaterModule
from modules.risk.worker_air import WorkerAirModule
from modules.risk.worker_soil import WorkerSoilModule
from ortam_core.chemicals.database import lookup, subchronic
from ortam_core.chemicals.rsl_import import DATA

HERE = Path(__file__).parent
TABLE = "rsl-2024-11"
REL = 1e-9


# ------------------------------------------------------------------ 1. RAGS E Exhibit B-3
def _b3():
    with (HERE / "data" / "rags_e_exhibit_b3.csv").open(encoding="utf-8") as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def _b3_compute(r):
    """Exhibit B-3 conditions: Cw = 1 mg/L = 1e-3 mg/cm3, t_event 0.58 h, EV 1, EF 350, ED 30, SA 18000, BW 70, AT 25550."""
    da = dq.da_event_water(1e-3, float(r["kp_cm_h"]), 0.58, False, float(r["fa"]), float(r["b"]),
                           float(r["tau_h"]), float(r["t_star_h"]))
    return da, dq.dad(da, 1.0, 30.0, 350.0, 18000.0, 70.0, 25550.0)


def test_exhibit_b3_benzene_by_hand():
    r = next(x for x in _b3() if x["chemical"] == "Benzene")
    da, dad = _b3_compute(r)
    assert float(f"{da:.1e}") == 1.7e-5 and float(f"{dad:.1e}") == 1.8e-3      # RAGS E p. B-13


def test_exhibit_b3_all_chemicals():
    rows = _b3()
    assert len(rows) == 209
    dev_da, dev_dad = [], []
    for r in rows:
        da, dad = _b3_compute(r)
        dev_da.append(abs(da / float(r["da_event_mg_cm2"]) - 1))
        dev_dad.append(abs(dad / float(r["dad_mg_kg_d"]) - 1))
    assert max(dev_da) < 0.10 and max(dev_dad) < 0.10
    assert statistics.median(dev_da) < 0.02 and statistics.median(dev_dad) < 0.02


def test_eq_3_1_from_printed_da_event():
    # DAD from the printed DA_event: Eq. 3.1 alone; 2-digit rounding of both columns -> 10 %
    for r in _b3():
        dad = dq.dad(float(r["da_event_mg_cm2"]), 1.0, 30.0, 350.0, 18000.0, 70.0, 25550.0)
        assert dad == pytest.approx(float(r["dad_mg_kg_d"]), rel=0.10), r["chemical"]


def test_inorganic_eq_3_4():
    # Eq. 3.4: DA_event = Kp*Cw*t_event; chromium(VI) Kp = 2e-3 cm/h (Exhibit 3-1, p. 3-5)
    assert dq.da_event_water(1e-3, 2e-3, 0.58, True) == pytest.approx(2e-3 * 1e-3 * 0.58)


def test_soil_dermal_eq_3_12_and_rags_f_eq_6():
    assert dq.da_event_soil(100.0, 0.07, 0.1) == pytest.approx(100.0 * 1e-6 * 0.07 * 0.1)
    # EC = CA*ET*EF*ED/AT(h): 10 ug/m3, 8 h/d, 250 d/yr, 25 yr, AT = 70 yr
    assert dq.exposure_concentration(10.0, 8.0, 250.0, 25.0, 70 * 365.0) == pytest.approx(
        10.0 * 8 * 250 * 25 / (70 * 365 * 24))


# ------------------------------------------------------------------ 2. consistency with the risk tools
def _rows():
    with (DATA / "2024-11" / "toxicity.csv").open(encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if r["cas"]
                and any(r[k] for k in ("sfo_per_mg_kg_day", "iur_per_ug_m3", "rfd_mg_kg_day", "rfc_mg_m3"))]
    return rows


def _check(tool_routes, from_dose: dict[str, float], risk: bool, who: str):
    for route in ("ingestion", "dermal", "inhalation"):
        v = getattr(tool_routes, route) if tool_routes is not None else None
        d = from_dose.get(route)
        if v is None:
            assert not d, (who, route, d)
        else:
            assert (one_hit(d) if risk else d) == pytest.approx(v, rel=REL), (who, route)


def _check_scalar(v, from_dose: dict[str, float], risk: bool, who: str):
    d = from_dose.get("inhalation")
    if v is None:
        assert not d, who
    else:
        assert (one_hit(d) if risk else d) == pytest.approx(v, rel=REL), who


def _chem(r, vf_field="vf_m3_kg", conc_field="concentration_mg_kg"):
    ch = {"cas": r["cas"], "name": r["name"], conc_field: 1.0}
    if r["volatile"] == "V" and vf_field:
        ch[vf_field] = 5000.0
    return ch


@pytest.mark.parametrize("module,kw,conc,vf", [
    (ResidentSoilModule, {}, "concentration_mg_kg", "vf_m3_kg"),
    (RecreatorSoilModule, {"exposure": {"ef_d_yr": [40] * 4, "et_h_d": [2] * 4, "irs_mg_d": [200, 200, 100, 100],
                                        "sa_soil_cm2_d": [2373, 2373, 6032, 6032],
                                        "af_mg_cm2": [0.2, 0.2, 0.07, 0.07]}}, "concentration_mg_kg", "vf_m3_kg"),
    (ResidentTapwaterModule, {}, "concentration_ug_l", None),
    (RecreatorSurfaceWaterModule, {"exposure": {"ef_d_yr": [52] * 4, "ev_per_d": [1] * 4, "et_event_h": [2] * 4,
                                                "irw_l_h": [0.12] * 4}}, "concentration_ug_l", None),
], ids=["resident-soil", "recreator-soil", "resident-tapwater", "recreator-surface-water"])
def test_resident_type_tools(module, kw, conc, vf):
    m, n = module(), 0
    for r in _rows():
        o = m.run(m.validate({"chemicals": [_chem(r, vf, conc)], **kw}))
        res = o.chemicals[0]
        c = lookup(TABLE, r["cas"], r["name"])
        assert res.doses, r["name"]
        _check(res.risk, dq.risk_by_route(res.doses, c), True, r["name"])
        _check(res.hq_child, dq.hq_by_route(res.doses, c, "noncancer-child"), False, r["name"])
        _check(res.hq_adult, dq.hq_by_route(res.doses, c, "noncancer-adult"), False, r["name"])
        n += 1
    assert n > 700


@pytest.mark.parametrize("receptor", ["composite", "outdoor", "indoor"])
def test_worker_soil(receptor):
    m = WorkerSoilModule()
    for r in _rows():
        res = m.run(m.validate({"receptor": receptor, "chemicals": [_chem(r)]})).chemicals[0]
        c = lookup(TABLE, r["cas"], r["name"])
        _check(res.risk, dq.risk_by_route(res.doses, c), True, r["name"])
        _check(res.hq, dq.hq_by_route(res.doses, c, "noncancer"), False, r["name"])


def test_construction_soil_subchronic():
    m = ConstructionSoilModule()
    for r in _rows():
        res = m.run(m.validate({"pef_sc_m3_kg": 1.0e6, "chemicals": [_chem(r, "vf_sc_m3_kg")]})).chemicals[0]
        c = lookup(TABLE, r["cas"], r["name"])
        srfd, srfc = subchronic(TABLE, c)
        _check(res.risk, dq.risk_by_route(res.doses, c), True, r["name"])
        _check(res.hq, dq.hq_by_route(res.doses, c, "noncancer", srfd or c.rfd_mg_kg_day, srfc or c.rfc_mg_m3),
               False, r["name"])


@pytest.mark.parametrize("module,kw", [(ResidentAirModule, {}), (WorkerAirModule, {"receptor": "indoor"})],
                         ids=["resident-air", "worker-air"])
def test_air_tools(module, kw):
    m = module()
    for r in _rows():
        res = m.run(m.validate({"chemicals": [{"cas": r["cas"], "name": r["name"], "concentration_ug_m3": 1.0}],
                                **kw})).chemicals[0]
        c = lookup(TABLE, r["cas"], r["name"])
        _check_scalar(res.risk, dq.risk_by_route(res.doses, c), True, r["name"])
        _check_scalar(res.hq, dq.hq_by_route(res.doses, c, "noncancer"), False, r["name"])


def test_fish():
    m = FishModule()
    for r in _rows():
        res = m.run(m.validate({"irf_mg_d": 54000.0, "chemicals": [
            {"cas": r["cas"], "name": r["name"], "concentration_mg_kg": 0.01}]})).chemicals[0]
        c = lookup(TABLE, r["cas"], r["name"])
        rk, hq = dq.risk_by_route(res.doses, c).get("ingestion"), dq.hq_by_route(res.doses, c, "noncancer").get("ingestion")
        assert (res.risk is None and not rk) or one_hit(rk) == pytest.approx(res.risk, rel=REL), r["name"]
        assert (res.hq is None and not hq) or hq == pytest.approx(res.hq, rel=REL), r["name"]


def test_no_doses_without_concentration():
    o = ResidentSoilModule().run(ResidentSoilModule().validate({"chemicals": [{"cas": "7440-38-2"}]}))
    assert o.chemicals[0].doses is None


def test_benzo_a_pyrene_mutagen_rows_carry_adaf():
    o = ResidentSoilModule().run(ResidentSoilModule().validate(
        {"chemicals": [{"cas": "50-32-8", "concentration_mg_kg": 1.0}]}))
    seg = {d.segment: d.weight for d in o.chemicals[0].doses if d.basis == "cancer" and d.route == "ingestion"}
    assert seg == {"0-2": 10.0, "2-6": 3.0, "6-16": 3.0, "16-26": 1.0}
