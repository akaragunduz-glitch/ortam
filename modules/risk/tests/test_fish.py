"""Ingestion of fish (tissue concentration). References: RSL UG 4.7 (PDF pp. 182-183)."""

import pytest

from modules.risk.fish import FishModule

IRF = 54000.0  # mg/day, example site-specific fish ingestion rate


def run(chemicals, **kw):
    m = FishModule()
    return m.run(m.validate({"irf_mg_d": IRF, "chemicals": chemicals, **kw}))


def test_noncancer_by_hand():
    # SL = THQ*AT*BW / ((1/RfD)*1e-6*EF*ED*IRF), AT = 365*EDres; copper RfD 0.04
    expected = 1.0 * 365 * 26 * 80 / ((1 / 0.04) * 1e-6 * 350 * 26 * IRF)
    assert run([{"cas": "7440-50-8"}]).chemicals[0].sl_nc_mg_kg == pytest.approx(expected)


def test_cancer_by_hand():
    # SL = TR*AT*BW / (SFO*1e-6*EF*ED*IRF), AT = 365*LT; dieldrin SFO 16
    expected = 1e-6 * 365 * 70 * 80 / (16 * 1e-6 * 350 * 26 * IRF)
    c = run([{"cas": "60-57-1"}]).chemicals[0]
    assert c.sl_ca_mg_kg == pytest.approx(expected)
    assert c.sl_basis == "c"


def test_forward_risk():
    sl = run([{"cas": "60-57-1"}]).chemicals[0].sl_ca_mg_kg
    o = run([{"cas": "60-57-1", "concentration_mg_kg": sl}])
    assert o.total_risk == pytest.approx(1e-6)


def test_irf_required():
    with pytest.raises(Exception):
        FishModule().validate({"chemicals": [{"cas": "60-57-1"}]})
