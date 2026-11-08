"""VISL workflow. Indoor air targets equal the verified RSL air SLs; subsurface targets follow
VISL UG (2014) pp. 8-9 with the OSWER (2015) attenuation factors."""

import pytest

from workflows.vapor_intrusion import visl_screening


def test_residential_benzene():
    out = visl_screening([{"cas": "71-43-2"}])
    c = out["visl"].chemicals[0]
    assert float(f"{c.target_ia_ug_m3:.2g}") == 0.36            # RSL Resident Air, Nov 2024
    assert c.target_subslab_ug_m3 == pytest.approx(c.target_ia_ug_m3 / 0.03)
    assert c.target_groundwater_ug_l == pytest.approx(c.target_ia_ug_m3 / (c.henry_used * 0.001 * 1000))


def test_commercial_uses_worker_air():
    out = visl_screening([{"cas": "79-01-6"}], land_use="commercial")
    assert float(f"{out['visl'].chemicals[0].target_ia_ug_m3:.2g}") == 3.0   # RSL Industrial Air, Nov 2024


# Reference: EPA VISL calculator output, resident, defaults (TR 1E-06, THQ 0.1, AF_ss 0.03,
# AF_gw 0.001, 25 C), generated 08 Nov 2026: docs/references/visl_calculator_resident_2026-11-08.xlsx
VISL_REF = {
    # cas: (target indoor air ug/m3, sub-slab ug/m3, groundwater ug/L, Cvp ug/m3, Chc ug/m3)
    "71-43-2": (0.36, 12.0, 1.59, 3.98e8, 4.06e8),
    "127-18-4": (4.17, 139.0, 5.76, 1.65e8, 1.49e8),
    "79-01-6": (0.209, 6.95, 0.518, 4.88e8, 5.15e8),
    "75-01-4": (0.168, 5.59, 0.147, 1.0e10, 1.0e10),
}


def sig3(x):
    return float(f"{x:.3g}")


@pytest.mark.parametrize("cas", list(VISL_REF))
def test_matches_visl_calculator(cas):
    ia, ss, gw, cvp, chc = VISL_REF[cas]
    c = visl_screening([{"cas": cas}], air_options={"target_hq": 0.1})["visl"].chemicals[0]
    assert sig3(c.target_ia_ug_m3) == ia
    assert sig3(c.target_subslab_ug_m3) == ss
    assert sig3(c.target_groundwater_ug_l) == gw
    assert c.pure_vapor_ug_m3 == pytest.approx(cvp, rel=0.006)


# VISL calculator 'VISL Risk' sheet: site sub-slab concentration -> indoor air -> risk and HQ
VISL_RISK = {"71-43-2": (100, 8.33e-06, 0.0959), "127-18-4": (200, 5.56e-07, 0.144),
             "79-01-6": (150, 9.41e-06, 2.16), "75-01-4": (300, 5.37e-05, 0.169)}


@pytest.mark.parametrize("cas", list(VISL_RISK))
def test_risk_from_subslab_matches_visl_calculator(cas):
    from workflows.vapor_intrusion import visl_risk_from_subsurface
    csg, cr, hq = VISL_RISK[cas]
    c = visl_risk_from_subsurface([{"cas": cas, "subslab_ug_m3": csg}])["air"].chemicals[0]
    assert sig3(c.risk) == cr
    assert sig3(c.hq) == hq
