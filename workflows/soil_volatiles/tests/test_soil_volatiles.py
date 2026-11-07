"""Volatile chemicals in soil: VF (air.ssl-vf) + soil SLs (risk tools) against the RSL tables.

References: RSL Summary Table, Nov 2024 (Resident Soil, Industrial Soil);
RSL UG 4.1.1, 4.2.1, 4.9.4, 5.12.
"""

import csv
import math

import pytest

from ortam_core.chemicals.rsl_import import DATA
from workflows.soil_volatiles import soil_screening


def sig2(x):
    return float(f"{x:.2g}")


def matches(got, ref):
    """Equal to two significant figures. Values that sit on a rounding boundary (e.g. 72.501 vs a
    table value of 72) are accepted when within half a unit of the last digit plus 1 % of it."""
    if sig2(got) == sig2(ref):
        return True
    unit = 10 ** (math.floor(math.log10(ref)) - 1)
    return abs(got - ref) <= 0.5 * unit * 1.01


def rows(column):
    with (DATA / "2024-11" / "toxicity.csv").open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            tox = any(r[k] for k in ("sfo_per_mg_kg_day", "iur_per_ug_m3", "rfd_mg_kg_day", "rfc_mg_m3"))
            if r["cas"] and r["volatile"] == "V" and r[column] and tox:
                yield r


@pytest.mark.parametrize("receptor,column", [("resident", "sl_resident_soil_mg_kg"),
                                             ("composite", "sl_industrial_soil_mg_kg")])
def test_volatile_soil_sls_match_rsl(receptor, column):
    rs = list(rows(column))
    assert len(rs) > 290
    chems = [{"cas": r["cas"], "name": r["name"]} for r in rs]
    out = soil_screening(chems, receptor=receptor)
    misses = [(r["name"], r[column], c.sl_mg_kg) for r, c in zip(rs, out["risk"].chemicals)
              if c.sl_mg_kg is None or not matches(c.sl_mg_kg, float(r[column]))]
    assert not misses, misses[:10]


def test_csat_flag_for_benzene_concentration():
    out = soil_screening([{"cas": "71-43-2"}])
    assert out["vf"].chemicals[0].liquid is True
    assert out["risk"].chemicals[0].sl_mg_kg < out["vf"].chemicals[0].csat_mg_kg


def test_pef_only_when_vf_not_calculable():
    out = soil_screening([{"cas": "7803-51-2", "name": "Phosphine"}])
    assert any("PEF only" in n for n in out["notes"])
    assert sig2(out["risk"].chemicals[0].sl_mg_kg) == 23
