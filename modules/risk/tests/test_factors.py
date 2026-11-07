"""Age-adjusted intake factors against the values printed in the RSL User's Guide."""

import pytest

from modules.risk.common import load_defaults
from modules.risk.equations import factors

D = load_defaults()


@pytest.mark.parametrize("name,expected,source", [
    ("ifs_mg_kg", 36750, "RSL UG 4.1.1.3, p. 40"),
    ("dfs_mg_kg", 103390, "RSL UG 4.1.1.3, p. 41"),
    ("ifsm_mg_kg", 166833.3, "RSL UG 4.1.1.4, p. 42"),
    ("dfsm_mg_kg", 428260, "RSL UG 4.1.1.4, p. 43"),
    ("ifw_l_kg", 327.95, "RSL UG 4.1.2.3, p. 71"),
    ("ifwm_l_kg", 1019.9, "RSL UG 4.1.2.4, p. 76"),
    ("dfw_cm2_event_kg", 2610650, "RSL UG 4.1.2.3, p. 74"),
    ("dfwm_cm2_event_kg", 8191633, "RSL UG 4.1.2.4, p. 78"),
    ("et_event_adj_h", 0.67077, "RSL UG 4.1.2.3, p. 72"),
    ("et_event_madj_h", 0.67077, "RSL UG 4.1.2.4, p. 78"),
])
def test_factor(name, expected, source):
    f = factors.compute(D)
    assert getattr(f, name) == pytest.approx(expected, rel=1e-5), source


@pytest.mark.parametrize("mw,log_kow,inside", [
    # RAGS Part E (2004), p. 3-7:
    #   Eq. 3.9:  -0.06831 <= 5.103e-4*MW + 0.05616*logKow <= 0.5577
    #   Eq. 3.10: -0.3010  <= -5.103e-4*MW + 0.05616*logKow <= 0.1758
    (78.115, 2.13, True),    # benzene
    (100.0, 4.0, True),      # Eq. 3.10: 0.1736 <= 0.1758
    (100.0, 4.5, False),     # Eq. 3.10: 0.2017 > 0.1758
    (100.0, -3.0, False),    # Eq. 3.9: -0.1174 < -0.06831
    (900.0, 2.0, False),     # Eq. 3.10: -0.3470 < -0.3010
])
def test_epd_boundaries_from_rags_e(mw, log_kow, inside):
    from modules.risk.equations.dermal_water import in_epd
    assert in_epd(mw, log_kow) is inside


@pytest.mark.parametrize("b,tau,expected", [
    # RSL UG 4.9.8, PDF p. 221: regression, clipped to [0, 1], rounded to 0.1
    (0.0507, 0.288, 1.0),    # benzene-like (small B, short lag): FA >= 1 -> 1
    (0.739, 23.62, None),    # fenvalerate: computed below by hand
])
def test_fraction_absorbed(b, tau, expected):
    import math
    from modules.risk.equations.dermal_water import fraction_absorbed
    if expected is None:
        lb, lt = math.log10(b), math.log10(tau)
        raw = 1.051232292 + 0.091016187 * lb - 0.286735467 * lt - 0.180504367 * lb * lt
        expected = round(min(max(raw, 0.0), 1.0), 1)
    assert fraction_absorbed(b, tau) == pytest.approx(expected)
