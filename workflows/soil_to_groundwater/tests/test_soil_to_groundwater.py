"""Risk-based SSL (DAF = 1) against the RSL Summary Table, Nov 2024.

Known deviations: the five tap water deviations documented in
modules/risk/tests/test_resident_tapwater.py propagate here, plus
3,3'-dichlorobenzidine (open question, docs/methods/vadose.md).
"""

import csv

from ortam_core.chemicals.rsl_import import DATA
from workflows.soil_to_groundwater import risk_based_ssl

KNOWN = {"Auramine", "Chlordane (alpha)", "Methylene-bis(N,N-dimethyl) Aniline, 4,4'-", "~Aroclor 1254",
         "Tetrachlorotoluene, p- alpha, alpha, alpha-", "Dichlorobenzidine, 3,3'-"}


def sig2(x):
    return float(f"{x:.2g}")


def test_risk_based_ssl_matches_rsl():
    with (DATA / "2024-11" / "toxicity.csv").open(encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if r["cas"] and r["ssl_risk_based_mg_kg"]]
    assert len(rows) > 740
    out = risk_based_ssl([{"cas": r["cas"], "name": r["name"]} for r in rows])
    got = {c.name: c.ssl_mg_kg for c in out["ssl"].chemicals}
    misses, known = [], []
    for r in rows:
        g, ref = got.get(r["name"]), float(r["ssl_risk_based_mg_kg"])
        if g is not None and sig2(g) == sig2(ref):
            continue
        (known if r["name"] in KNOWN else misses).append((r["name"], ref, g))
    assert not misses, misses[:10]
    for name, ref, g in known:
        assert abs(g - ref) / ref < 0.15, name
