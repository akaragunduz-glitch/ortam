"""Workflow vapor.johnson-ettinger -> risk.resident-air against the J&E spreadsheet risk output.

Reference: data/reference/je/v6-2017/reference_risk.csv, cells MODEL!G153 (cancer risk)
and G154 (HQ) of EPA's J&E spreadsheet v6.0 for the 16 cases of reference_cases.csv
(residential defaults: ED 26 yr, EF 350 d/yr, ET 24 h/d, ATc 70 yr, TR 1e-6).
The spreadsheet uses its 2017 toxicity values; ORTAM uses RSL Nov 2024. Two documented
differences remain; everything else must agree to 1e-6:
  - Vinyl chloride RfC: J&E 0.1 mg/m3 (IRIS), RSL Nov 2024 0.0511 mg/m3 (key A);
    the HQ scales with 0.1 / 0.0511.
  - Trichloroethylene: the J&E spreadsheet splits the IUR into rounded parts,
    1.0e-6 (mutagenic, x MMOAF 72) + 3.1e-6 (x ED 26) (ToxSummary!E829:E830);
    the RSL uses 4.1e-6 x CAF 0.244 / MAF 0.756 (RSL UG 4.1.3.5, p. 101),
    i.e. 1.0004e-6 + 3.0996e-6. The risks differ by 1.2e-4 (relative).
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from workflows.vapor_intrusion import je_risk

ROOT = Path(__file__).resolve().parents[3]
REF = ROOT / "data/reference/je/v6-2017"
CASES = {c["case"]: c for c in csv.DictReader((REF / "reference_cases.csv").open(encoding="utf-8"))}
RISK = list(csv.DictReader((REF / "reference_risk.csv").open(encoding="utf-8")))
SOURCE = {"Groundwater": "groundwater", "Exterior Soil Gas": "exterior-soil-gas",
          "Sub-slab Soil Gas": "subslab-soil-gas"}


def _je_inputs(case: dict) -> dict:
    return {
        "source": SOURCE[case["source"]],
        "depth_to_source_m": float(case["ls_m"]),
        "temperature_c": float(case["ts_c"]),
        "foundation": case["foundation"],
        "strata": [{"soil": case[f"soil_{k}"], "thickness_m": float(case[f"h_{k}_m"])}
                   for k in "abc" if case[f"soil_{k}"] != "Not Present"],
        "chemicals": [{"name": case["chemical"], "concentration": float(case["c_medium"])}],
    }


@pytest.mark.parametrize("ref", RISK, ids=lambda r: f"case{r['case']}-{r['chemical']}")
def test_risk_matches_je_spreadsheet(ref):
    out = je_risk(_je_inputs(CASES[ref["case"]]))
    c = out["air"].chemicals[0]
    assert out["je"].chemicals[0].cia_ug_m3 == pytest.approx(float(ref["cia_ug_m3"]), rel=1e-9)
    risk_rel = 2e-4 if ref["cas"] == "79-01-6" else 1e-6
    hq_ref = float(ref["hazard_quotient"])
    if ref["cas"] == "75-01-4":
        hq_ref *= 0.1 / 0.0511
    assert c.risk == pytest.approx(float(ref["cancer_risk"]), rel=risk_rel)
    assert c.hq == pytest.approx(hq_ref, rel=1e-6)
