"""Three-point tests. Reference: hand calculation (plane equation)."""

import pytest

from modules.groundwater.three_point import ThreePointModule


def run(wells, **kw):
    m = ThreePointModule()
    data = {"wells": [dict(zip(("name", "x_m", "y_m", "head_m"), w)) for w in wells], **kw}
    return m.run(m.validate(data))


def test_eastward_flow():
    # Head drops 1 m over 100 m to the east -> i = 0.01, azimuth 90 deg
    o = run([("A", 0, 0, 10), ("B", 100, 0, 9), ("C", 0, 100, 10)])
    assert o.gradient == pytest.approx(0.01)
    assert o.flow_azimuth_deg == pytest.approx(90.0)


def test_southwest_flow_and_velocities():
    # Head rises 0.005 m/m to the north and east -> flow to the southwest (225 deg)
    o = run([("A", 0, 0, 10), ("B", 200, 0, 11), ("C", 0, 200, 11)], K_m_d=10, n_e=0.2)
    i = (0.005**2 * 2) ** 0.5
    assert o.gradient == pytest.approx(i)
    assert o.flow_azimuth_deg == pytest.approx(225.0)
    assert o.darcy_flux_m_d == pytest.approx(10 * i)
    assert o.seepage_velocity_m_d == pytest.approx(10 * i / 0.2)


def test_collinear_wells_raise():
    with pytest.raises(ValueError):
        run([("A", 0, 0, 10), ("B", 50, 0, 9), ("C", 100, 0, 8)])


def test_missing_well_rejected():
    with pytest.raises(Exception):
        ThreePointModule().validate({"wells": [{"name": "A", "x_m": 0, "y_m": 0, "head_m": 1}]})
