"""Soil screening levels protective of groundwater (SSL).

References: RSL UG 4.8.3 (Method 1, PDF p. 190), 4.8.4 (Method 2, p. 193),
4.8.5 (DAF, p. 194); SSG 2002 Eq. 4-10 to 4-14; RSL Summary Table, Nov 2024
(MCL-based SSL, DAF = 1). The risk-based SSL column is verified in
workflows/soil_to_groundwater.
"""

import csv
import math
import re

import pytest

from modules.vadose.ssl_leaching import SslLeachingModule
from ortam_core.chemicals.rsl_import import DATA


def sig2(x):
    return float(f"{x:.2g}")


def run(chemicals, **kw):
    m = SslLeachingModule()
    return m.run(m.validate({"chemicals": chemicals, **kw}))


def _num(text):
    m = re.match(r"[0-9.Ee+-]+", text or "")
    return float(m.group()) if m else None


def test_mcl_based_ssl_matches_rsl_table():
    with (DATA / "2024-11" / "toxicity.csv").open(encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if r["cas"] and _num(r["mcl_ug_l"]) and _num(r["ssl_mcl_based_mg_kg"])]
    assert len(rows) > 85
    misses = []
    for r in rows:
        c = run([{"cas": r["cas"], "name": r["name"], "target": "mcl"}]).chemicals[0]
        if c.ssl_mg_kg is None or sig2(c.ssl_mg_kg) != sig2(_num(r["ssl_mcl_based_mg_kg"])):
            misses.append((r["name"], r["ssl_mcl_based_mg_kg"], c.ssl_mg_kg))
    assert not misses, misses[:10]


def test_method1_benzene_by_hand():
    # Benzene: Koc 145.8, H' 0.2269 (RSL Nov 2024); Cw = 5 ug/L (MCL)
    h = 0.22690106295994
    n = 1 - 1.5 / 2.65
    expected = 5 / 1000 * (145.8 * 0.002 + (0.3 + (n - 0.3) * h) / 1.5)
    c = run([{"cas": "71-43-2", "target_ug_l": 5.0}]).chemicals[0]
    assert c.ssl_mg_kg == pytest.approx(expected)


def test_inorganic_uses_table_kd_and_zero_henry():
    # Arsenic: Kd 29 L/kg (SSG 1996 via RSL UG 4.8, p. 185), H' = 0 for inorganics
    expected = 10 / 1000 * (29 + 0.3 / 1.5)
    c = run([{"cas": "7440-38-2", "target_ug_l": 10.0}]).chemicals[0]
    assert c.ssl_mg_kg == pytest.approx(expected)
    assert c.henry_used == 0.0


def test_dilution_factor_scales_linearly():
    a = run([{"cas": "71-43-2", "target_ug_l": 5.0}]).chemicals[0].ssl_mg_kg
    b = run([{"cas": "71-43-2", "target_ug_l": 5.0}], daf=20).chemicals[0].ssl_mg_kg
    assert b == pytest.approx(20 * a)


def test_site_specific_daf_by_hand():
    # SSG Eq. 4-11, 4-12: K 20 m/d -> 7300 m/yr, i 0.01, L 45 m, da 10 m, I 0.18 m/yr
    k, i, length, da, inf = 7300.0, 0.01, 45.0, 10.0, 0.18
    d = (0.0112 * length ** 2) ** 0.5 + da * (1 - math.exp(-length * inf / (k * i * da)))
    daf = 1 + k * i * d / (inf * length)
    o = run([{"cas": "71-43-2", "target_ug_l": 5.0}],
            site={"hydraulic_conductivity_m_yr": k, "hydraulic_gradient": i, "source_length_m": length,
                  "aquifer_thickness_m": da})
    assert o.mixing_zone_depth_m == pytest.approx(d)
    assert o.daf == pytest.approx(daf)


def test_method2_mass_limit_by_hand():
    # SSG Eq. 4-14: SSL = Cw*I*ED/(rho_b*ds); Cw 0.005 mg/L, I 0.18, ED 70, ds 3 m
    expected = 0.005 * 0.18 * 70 / (1.5 * 3.0)
    c = run([{"cas": "71-43-2", "target_ug_l": 5.0}], method="mass-limit", source_depth_m=3.0).chemicals[0]
    assert c.ssl_mg_kg == pytest.approx(expected)


def test_mcl_annotation_reported():
    with (DATA / "2024-11" / "toxicity.csv").open(encoding="utf-8") as f:
        r = next(r for r in csv.DictReader(f) if "(" in r["mcl_ug_l"])
    c = run([{"cas": r["cas"], "name": r["name"], "target": "mcl"}]).chemicals[0]
    assert c.target_ug_l == pytest.approx(_num(r["mcl_ug_l"]))
    assert any("MCL" in n for n in c.notes)


def test_target_required():
    with pytest.raises(Exception):
        run([{"cas": "71-43-2"}])
