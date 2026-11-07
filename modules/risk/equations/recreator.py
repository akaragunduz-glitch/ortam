"""Recreator parameters from age-segment values.

The recreator soil/sediment and surface water equations (RSL UG 4.6.1 and
4.6.2, PDF pp. 124-182) have the same form as the resident equations
(4.1.1 and 4.1.2) with recreator parameters. Child, adult and age-adjusted
parameters are ED-weighted means of the age-segment values:
  X_child = sum_{0-2, 2-6} ED_i*X_i / sum ED_i     (4.6.1.7, p. 149; 4.6.2.7, p. 178)
  X_adult = sum_{6-16, 16-26} ED_i*X_i / sum ED_i  (p. 150; p. 180)
  X_res   = sum_all ED_i*X_i / sum ED_i            (p. 151-152; p. 181-182)
For surface water, the daily ingestion IRW_day = EV * ET_event * IRW (L/h)
(4.6.2.1, p. 153); there is no inhalation route (p. 153).
This module maps recreator inputs onto the resident parameter structure
(modules.risk.common.Defaults), so that equations/soil.py and
equations/tapwater.py are reused unchanged.
"""

from __future__ import annotations

from modules.risk.common import Defaults

CHILD, ADULT = (0, 1), (2, 3)


def _wmean(values: list[float], ed: list[float], idx) -> float:
    return sum(ed[i] * values[i] for i in idx) / sum(ed[i] for i in idx)


def build(base: Defaults, seg: dict[str, list[float]]) -> Defaults:
    """seg: per-parameter lists of four age-segment values (0-2, 2-6, 6-16, 16-26)."""
    ed = seg["ed_yr"]
    allidx = (0, 1, 2, 3)
    n = len(ed)
    get = lambda k: seg.get(k, [0.0] * n)
    irw_day = [ev * et * irw for ev, et, irw in zip(get("ev_per_d"), get("et_event_h"), get("irw_l_h"))]

    p = dict(base.p)
    for key, vals in (("ef", seg["ef_d_yr"]), ("et", get("et_h_d"))):
        unit = "d_yr" if key == "ef" else "h_d"
        p[f"{key}_child_{unit}"] = _wmean(vals, ed, CHILD)
        p[f"{key}_adult_{unit}"] = _wmean(vals, ed, ADULT)
        p[f"{key}_res_{unit}"] = _wmean(vals, ed, allidx)
    p["ed_child_yr"] = ed[0] + ed[1]
    p["ed_adult_yr"] = ed[2] + ed[3]
    p["ed_res_yr"] = sum(ed)
    for name, vals in (("bw_{}_kg", seg["bw_kg"]), ("irs_{}_mg_d", get("irs_mg_d")),
                       ("sa_{}_cm2_d", get("sa_soil_cm2_d")), ("af_{}_mg_cm2", get("af_mg_cm2")),
                       ("sa_water_{}_cm2", seg["sa_water_cm2"]), ("ev_{}_per_d", get("ev_per_d")),
                       ("et_event_{}_h", get("et_event_h"))):
        p[name.format("child")] = _wmean(vals, ed, CHILD)
        p[name.format("adult")] = _wmean(vals, ed, ADULT)
    # IRW_day (child/adult) = EV * ET_event * IRW with the ED-weighted means (4.6.2.1, p. 153; 4.6.2.2, p. 156)
    for who, idx in (("child", CHILD), ("adult", ADULT)):
        p[f"irw_{who}_l_d"] = (_wmean(get("ev_per_d"), ed, idx) * _wmean(get("et_event_h"), ed, idx)
                               * _wmean(get("irw_l_h"), ed, idx))

    segments = []
    for i, label in enumerate(("0-2", "2-6", "6-16", "16-26")):
        segments.append({"label": label, "ed_yr": ed[i], "ef_d_yr": seg["ef_d_yr"][i], "et_h_d": get("et_h_d")[i],
                         "bw_kg": seg["bw_kg"][i], "irs_mg_d": get("irs_mg_d")[i],
                         "sa_soil_cm2_d": get("sa_soil_cm2_d")[i], "af_mg_cm2": get("af_mg_cm2")[i],
                         "irw_l_d": irw_day[i], "sa_water_cm2": seg["sa_water_cm2"][i],
                         "ev_per_d": get("ev_per_d")[i], "et_event_h": get("et_event_h")[i],
                         "adaf": seg["adaf"][i]})
    return Defaults(p=p, segments=segments, tce=base.tce, rba_by_cas=base.rba_by_cas)
