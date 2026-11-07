"""Construction worker soil. References: RSL UG 4.5 (PDF pp. 116-124). Hand calculations."""

import pytest

from modules.risk.construction_soil import ConstructionSoilModule

PEF_SC = 1.0e6  # m3/kg, arbitrary user value: PEFsc calculation is task D1.2


def run(chemicals, **kw):
    m = ConstructionSoilModule()
    return m.run(m.validate({"pef_sc_m3_kg": PEF_SC, "chemicals": chemicals, **kw}))


def test_subchronic_noncancer_by_hand():
    # Acephate: subchronic RfDo 0.004 mg/kg-day (RSL Nov 2024 subchronic table), ABSd 0.1, GIABS 1
    # RSL UG 4.5.1.1 (p. 117): AT = EW*7*ED = 350 d; EF = EW*DW = 250 d/yr; ED 1 yr; IRS 330; BW 80
    at, ef = 50 * 7 * 1, 50 * 5
    ing = 1.0 * at * 80 / ((1 / 0.004) * 1e-6 * ef * 1 * 330)
    derm = 1.0 * at * 80 / ((1 / (0.004 * 1.0)) * 1e-6 * ef * 1 * 3527 * 0.3 * 0.1)
    c = run([{"cas": "30560-19-1"}]).chemicals[0]
    assert c.toxicity_duration == "subchronic"
    assert c.sl_nc_mg_kg.ingestion == pytest.approx(ing)
    assert c.sl_nc_mg_kg.dermal == pytest.approx(derm)


def test_chronic_option():
    at, ef = 350, 250
    ing = 1.0 * at * 80 / ((1 / 3e-4) * 1e-6 * ef * 1 * 330)
    c = run([{"cas": "30560-19-1"}], toxicity_duration="chronic").chemicals[0]
    assert c.sl_nc_mg_kg.ingestion == pytest.approx(ing)


def test_fallback_to_chronic_is_flagged():
    # Inorganic arsenic: not in the RSL subchronic table -> chronic RfD with a warning
    o = run([{"cas": "7440-38-2"}])
    assert o.chemicals[0].toxicity_duration == "chronic (no subchronic value)"
    assert any("subchronic" in w for w in o.warnings)


def test_cancer_by_hand():
    # Dieldrin: SFO 16, IUR 4.6E-03, ABSd 0.1 (RSL UG 4.5.1.2, p. 119)
    at = 365 * 70
    ing = 1e-6 * at * 80 / (16 * 1e-6 * 250 * 1 * 330)
    derm = 1e-6 * at * 80 / ((16 / 1.0) * 1e-6 * 250 * 1 * 3527 * 0.3 * 0.1)
    inh = 1e-6 * at / (4.6e-3 * 1000 * 250 * 1 * 8 / 24 * (1 / PEF_SC))
    c = run([{"cas": "60-57-1"}]).chemicals[0]
    assert c.sl_ca_mg_kg.ingestion == pytest.approx(ing)
    assert c.sl_ca_mg_kg.dermal == pytest.approx(derm)
    assert c.sl_ca_mg_kg.inhalation == pytest.approx(inh)


def test_pef_sc_required():
    m = ConstructionSoilModule()
    with pytest.raises(Exception):
        m.validate({"chemicals": [{"cas": "60-57-1"}]})


def test_volatile_needs_vf_sc():
    with pytest.raises(ValueError, match="vf_sc"):
        run([{"cas": "71-43-2"}])
