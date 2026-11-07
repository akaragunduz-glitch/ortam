"""Resident tap water.

References: RSL UG 4.1.2 (PDF pp. 62-95), 4.9.8 (pp. 217-225); RAGS Part E (2004)
p. 3-7, Eq. 3.9-3.10 (EPD); RSL Summary Table, November 2024.
"""

import csv
import math

import pytest

from modules.risk.resident_tapwater import ResidentTapwaterModule
from modules.risk.rsl_import import DATA

# Known deviations from the RSL Summary Table (Nov 2024), investigated and documented
# in docs/methods/risk.md. Each entry: name -> (RSL value ug/L, reason).
KNOWN_DEVIATIONS = {
    "~Aroclor 1242": (0.0078, "PCB special case, RSL UG 5.8"),
    "~Aroclor 1248": (0.0078, "PCB special case, RSL UG 5.8"),
    "~Aroclor 1254": (0.0078, "PCB special case, RSL UG 5.8"),
    "~Aroclor 1260": (0.0078, "PCB special case, RSL UG 5.8"),
    "~Mercuric Chloride (and other Mercury salts)": (5.7, "mercury salt special case, RSL UG 5.19"),
    "Auramine": (0.078, "dermal parameters; open question"),
    "Chlordane (alpha)": (3.6, "dermal parameters; open question"),
    "Methylene-bis(N,N-dimethyl) Aniline, 4,4'-": (0.7, "dermal parameters; open question"),
    "Tetrachlorotoluene, p- alpha, alpha, alpha-": (0.0017, "dermal parameters; open question"),
}


def sig2(x):
    return float(f"{x:.2g}")


def run(chemicals, **kw):
    m = ResidentTapwaterModule()
    return m.run(m.validate({"chemicals": chemicals, **kw}))


def rows():
    with (DATA / "2024-11" / "toxicity.csv").open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            # Rows without toxicity values (e.g. lead: action level, RSL UG 5.2) are excluded.
            has_tox = any(r[k] for k in ("sfo_per_mg_kg_day", "iur_per_ug_m3", "rfd_mg_kg_day", "rfc_mg_m3"))
            if r["cas"] and r["sl_tapwater_ug_l"] and has_tox:
                yield r


def test_all_chemicals_match_rsl_tapwater():
    rs = list(rows())
    assert len(rs) > 800
    misses, known = [], []
    for r in rs:
        got = run([{"cas": r["cas"], "name": r["name"]}]).chemicals[0].sl_ug_l
        if got is not None and sig2(got) == sig2(float(r["sl_tapwater_ug_l"])):
            continue
        (known if r["name"] in KNOWN_DEVIATIONS else misses).append((r["name"], r["sl_tapwater_ug_l"], got))
    assert not misses, misses[:10]
    # Known deviations stay within 30 % of the RSL value
    for name, ref, got in known:
        assert abs(got - float(ref)) / float(ref) < 0.30, name


@pytest.mark.parametrize("cas,expected,method", [
    ("7440-50-8", 800, None),             # copper: noncancer, inorganic dermal
    ("71-43-2", 0.46, "standard"),        # benzene: organic dermal + volatile inhalation
    ("7440-38-2", 0.052, "standard"),     # arsenic
    ("60-57-1", 0.0018, "standard"),      # dieldrin: long lag time (ET <= t*)
    ("18540-29-9", 0.11, "mutagen"),      # chromium(VI)
    ("50-32-8", 0.025, "mutagen"),        # benzo(a)pyrene
    ("75-01-4", 0.019, "vinyl chloride"),
    ("79-01-6", 0.49, "trichloroethylene"),
])
def test_selected_chemicals(cas, expected, method):
    c = run([{"cas": cas}]).chemicals[0]
    assert sig2(c.sl_ug_l) == expected
    if method:
        assert c.cancer_method == method


def test_benzene_dermal_by_hand():
    # RSL UG 4.1.2.3 (PDF pp. 71-74). Benzene is organic; ET_event-adj = 0.67077 h is
    # below t* = 0.691 h, so the short-event (ET <= t*) equation applies:
    #   SL = DA_event*1000 / (2*FA*Kp*sqrt(6*tau*ET/pi)),  DA_event = TR*AT / ((SFO/GIABS)/1000*DFW)
    sfo, kp, tau = 0.055, 0.0149, 0.28793386679266
    da = 1e-6 * 365 * 70 / ((sfo / 1.0) / 1000 * 2610650)
    et = (6 * 0.54 + 20 * 0.71) / 26
    fa = 1.0  # RSL UG 4.9.8 regression gives 1.02 -> clipped to 1
    expected = da * 1000 / (2 * fa * kp * math.sqrt(6 * tau * et / math.pi))
    c = run([{"cas": "71-43-2"}]).chemicals[0]
    assert c.sl_ca_ug_l.dermal == pytest.approx(expected, rel=1e-9)


def test_dermal_excluded_outside_epd():
    # Caprolactam is inside the EPD; a hypothetical organic with logKow 8, MW 100 is outside
    o = run([{"name": "Lipophilic test", "rfd_mg_kg_day": 0.01, "mw_g_mol": 100, "log_kow": 8.0,
              "kp_cm_h": 0.1, "b_dimensionless": 0.5, "tau_event_h": 1.0, "t_star_h": 2.4}],
            toxicity_table="none")
    assert o.chemicals[0].sl_nc_child_ug_l.dermal is None
    assert any("EPD" in n for n in o.chemicals[0].notes)


def test_forward_risk():
    sl = run([{"cas": "71-43-2"}]).chemicals[0].sl_ca_ug_l.total
    o = run([{"cas": "71-43-2", "concentration_ug_l": sl}])
    assert o.total_risk == pytest.approx(1e-6)
