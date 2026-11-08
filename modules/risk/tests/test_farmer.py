"""risk.farmer against the RAIS Chemical Risk Calculator (farmer, site-specific mode).

Reference: docs/references/rais_farmer_soil_output_2026-11-08.xlsx, saved by the owner from
the RAIS calculator on 2026-11-08: 8 chemicals (arsenic, benzene, benzo[a]pyrene, cadmium,
dieldrin, TCE, vinyl chloride, zinc), two scenarios:
  "Contaminated Food Products"  produce, beef and dairy concentrations given (4.8.1)
  "Contaminated Soil"           soil concentration with RAIS transfer factors (4.8.3)
The RAIS toxicity values (toxicity_table "none"), transfer factors, VF, ABSd, GIABS and RBA printed in the output are
passed to the module, so that the comparison tests the exposure equations. RAIS prints 3
significant figures; agreement is required to 3 significant figures (soil scenario 1 %, because
the transfer factors are printed with 3 figures). Finding: in the soil scenario the RAIS
'Total HI' omits the soil ingestion HQ; ORTAM reports the full sum and the test checks that
RAIS total = ORTAM total - soil ingestion HQ.
"""

from __future__ import annotations

from pathlib import Path

import openpyxl
import pytest

from modules.risk.farmer import FarmerModule

XLSX = Path(__file__).resolve().parents[3] / "docs" / "references" / "rais_farmer_soil_output_2026-11-08.xlsx"


def _sheet(prefix: str) -> list[dict]:
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    ws = next(w for w in wb.worksheets if w.title.startswith(prefix))
    rows = list(ws.iter_rows(values_only=True))
    i = next(k for k, r in enumerate(rows) if r[0] == "Chemical")
    hdr = [" ".join(str(h).split()) if h else "" for h in rows[i]]
    out = []
    for r in rows[i + 1:]:
        if not r[0] or str(r[0]).startswith("*") or str(r[0]).startswith("Output"):
            continue
        out.append({h: v for h, v in zip(hdr, r)})
    return out


def _num(v):
    return None if v in (None, "-", "") else float(v)


def sig3(x: float) -> float:
    return float(f"{x:.3g}")


def _same(mine, ref, what, tol=6e-3):
    ref = _num(ref)
    if ref is None:
        assert mine is None or mine == 0, what
        return
    assert mine is not None, what
    assert sig3(mine) == pytest.approx(ref, rel=1e-9) or mine == pytest.approx(ref, rel=tol), (what, mine, ref)


def _tox(row: dict, soil: bool) -> dict:
    t = {"cas": row["CAS Number"], "name": row["Chemical"]}
    for key, col in (("sfo_per_mg_kg_day", "Oral Slope Factor (SFO) (mg/kg-day)-1"),
                     ("iur_per_ug_m3", "Inhalation Unit Risk (IUR) (ug/m3)-1"),
                     ("rfd_mg_kg_day", "Oral Reference Dose (RfD) (mg/kg-day)")):
        t[key] = _num(row.get(col))
    t["mutagen"] = row["Mutagen?"] == "Yes"
    t["volatile"] = row["Volatile?"] == "Yes"
    if soil:
        t["rfc_mg_m3"] = _num(row.get("Inhalation Reference Concentration (RfC) (mg/m3)"))
        t["giabs"] = _num(row["RAGS Part E Gastrointestinal Absorption Factor (GIABS) (unitless)"])
        t["absd"] = _num(row["RAGS Part E Dermal Absorption Factor (ABSd) (unitless)"])
        t["rba"] = _num(row["Relative Bioavailability (RBA)"])
        t["bv_wet"] = _num(row["Soil-to-Plant Transfer Factor (BVwet) (kg-dry soil per kg-fresh plant)"])
        t["bv_dry"] = _num(row["Soil-to-Plant Transfer Factor (BVdry) (kg-dry soil per kg-dry plant)"])
        t["tf_dairy_d_kg"] = _num(row["Milk Transfer Factor (TFdairy) (day/kg)"])
        t["tf_beef_d_kg"] = _num(row["Beef Transfer Factor (TFbeef) (day/kg)"])
        t["vf_m3_kg"] = _num(row["Volatilization Factor Selected (m3/kg)"])
        t["soil_mg_kg"] = _num(row["Soil Concentration (mg/kg)"])
    else:
        t["produce_mg_kg"] = _num(row["Produce Concentration (mg/kg)"])
        t["beef_mg_kg"] = _num(row["Beef Concentration (mg/kg)"])
        t["dairy_mg_kg"] = _num(row["Dairy Concentration (mg/kg)"])
    return {k: v for k, v in t.items() if v is not None}


FOOD = _sheet("Contaminated Food Products R")
SOIL = _sheet("Contaminated Soil Risk")
LABEL = {"produce": "Produce Ingestion", "beef": "Beef Ingestion", "dairy": "Dairy Ingestion",
         "ingestion": "Ingestion", "dermal": "Dermal", "inhalation": "Inhalation"}
WHO = {"child": "Child", "adult": "Adult", "age_adjusted": "Age-adjusted"}


def _check(row, r, routes, tol=6e-3):
    for route in routes:
        L = LABEL[route]
        s = r.cdi.get(route)
        unit_ca = "(ug/m3)" if route == "inhalation" else "(mg/kg-day)"
        unit_nc = "(mg/m3)" if route == "inhalation" else "(mg/kg-day)"
        _same(s.cancer if s else None, row[f"{L} Carcinogenic CDI {unit_ca}"], (r.name, route, "cancer CDI"), tol)
        for who, W in WHO.items():
            col = f"{W} {L} Noncarcinogenic CDI {unit_nc}"
            if who == "age_adjusted" and col not in row:
                col = f"Adjusted {L} Noncarcinogenic CDI {unit_nc}"
            _same(getattr(s, who) if s else None, row[col], (r.name, route, who), tol)
        rk = next((k for k in row if k.replace("  ", " ") == f"{L} Risk"), None)
        _same(r.risk.get(route), row[rk], (r.name, route, "risk"), tol)
        for who, W in WHO.items():
            hk = next(k for k in row if k.replace("  ", " ") == f"{W} {L} HQ")
            _same(r.hq.get(who, {}).get(route), row[hk], (r.name, route, who, "HQ"), tol)
    _same(r.total_risk, row["Total Risk"], (r.name, "total risk"), tol)
    for who, W in WHO.items():
        assert r.hi.get(who, 0.0) == pytest.approx(sum(r.hq.get(who, {}).values()))    # RAGS A p. 8-13
        ref = _num(row[f"{W} Total HI"])
        if ref is None:
            continue
        if "ingestion" in routes:
            # Finding: the RAIS 'Total HI' of the soil scenario omits the soil ingestion HQ (the total
            # risk includes it). Reproduced here; ORTAM reports the full sum.
            assert sig3(r.hi[who] - r.hq[who].get("ingestion", 0.0)) == pytest.approx(ref, rel=1e-9) or \
                r.hi[who] - r.hq[who].get("ingestion", 0.0) == pytest.approx(ref, rel=tol), (r.name, who)
        else:
            _same(r.hi.get(who), ref, (r.name, who, "HI"), tol)


def test_reference_present():
    assert len(FOOD) == 8 and len(SOIL) == 8


@pytest.mark.parametrize("row", FOOD, ids=[r["Chemical"] for r in FOOD])
def test_food_products_scenario(row):
    m = FarmerModule()
    r = m.run(m.validate({"scenario": "food", "toxicity_table": "none",
                          "chemicals": [_tox(row, soil=False)]})).chemicals[0]
    _check(row, r, ("produce", "beef", "dairy"))


@pytest.mark.parametrize("row", SOIL, ids=[r["Chemical"] for r in SOIL])
def test_soil_scenario(row):
    m = FarmerModule()
    r = m.run(m.validate({"scenario": "soil", "toxicity_table": "none", "chemicals": [_tox(row, soil=True)]})).chemicals[0]
    # transfer factors are printed with 3 significant figures (RAIS uses more digits): 1 %
    _check(row, r, ("produce", "beef", "dairy", "ingestion", "dermal", "inhalation"), tol=1e-2)
