"""Doses: chronic daily intake (CDI), dermal absorbed dose (DAD), exposure concentration (EC).

Sources (printed page; PDF page in brackets):
  CDI, generic intake   I = C*CR*EF*ED/(BW*AT)          RAGS A (1989) Exhibit 6-9, p. 6-21 [130]
  CDI, soil ingestion   I = CS*IR*CF*FI*EF*ED/(BW*AT)   RAGS A Exhibit 6-14, p. 6-40 [149]
  DAD                   DAD = DA_event*EV*ED*EF*SA/(BW*AT)   RAGS E (2004) Eq. 3.1, p. 3-3 [24]
  DA_event, water       Eq. 3.2 (t_event <= t*), 3.3 (t_event > t*), p. 3-4 [25];
                        Eq. 3.4 (inorganics), p. 3-5 [26]
  DA_event, soil        DA_event = C_soil*CF*AF*ABS_d   RAGS E Eq. 3.12, p. 3-9 [30]
  EC                    EC = CA*ET*EF*ED/AT (AT in hours)   RAGS F (2009) Eq. 6, p. 14 [22]; Eq. 8, p. 18 [26]
  Risk and hazard       Risk = CDI*SF               RAGS A p. 8-6 [195]
                        HQ = E/RfD                  RAGS A p. 8-11 [200]
                        SF_ABS = SFO/ABS_GI, RfD_ABS = RfD*ABS_GI   RAGS E Eq. 4.2, 4.3, pp. 4-3, 4-4 [45-46]
                        Risk = IUR*EC; HQ = EC/(RfC*1000)   RAGS F Eq. 11, p. 22 [30]; Eq. 12, p. 24 [32]
Averaging times: noncancer AT = ED*365 d; cancer AT = LT*365 d (RAGS A Exhibit 6-14;
RAGS E Eq. 3.1). For the construction worker the noncancer AT is EW*7*ED (RSL UG 4.5.1.1).

Age segments, ADAF, vinyl chloride and TCE follow the RSL User's Guide equations
that the risk tools use (modules/risk/equations/soil.py, tapwater.py, air.py), so that
  risk_route = sum(dose * weight * toxicity)   and   HQ_route = sum(dose * weight) / toxicity
reproduce the tools' risk and HQ exactly. The weight of a dose row carries
RBA (soil ingestion), ADAF (mutagens; RSL UG 4.1.1.4), CAF/MAF (TCE; RSL UG 4.1.1.6).
Vinyl chloride has an early-life term without averaging time (RSL UG 4.1.1.5,
4.1.2.5, 4.1.3.4); it is reported as a separate row with segment 'early-life'.

Units: CDI and DAD in mg/kg-day; EC in ug/m3.
Fortran analogy: each builder fills a table of rows (like a COMMON block of
arrays); the tools only pass parameters in.
"""

from __future__ import annotations

import math
from typing import Literal

from pydantic import BaseModel, Field

from modules.risk.common import Defaults
from modules.risk.equations import dermal_water as dw
from modules.risk.equations.factors import Factors
from modules.risk.equations.worker import Worker
from ortam_core.chemicals.database import Chemical

Route = Literal["ingestion", "dermal", "inhalation"]
Basis = Literal["noncancer", "noncancer-child", "noncancer-adult", "cancer"]


class Dose(BaseModel):
    route: Route
    quantity: Literal["CDI", "DAD", "EC"]
    basis: Basis
    segment: str = Field(description="child, adult, resident, worker, an age segment (0-2 ...) or early-life")
    value: float = Field(description="CDI/DAD in mg/kg-day; EC in ug/m3")
    unit: str
    weight: float = Field(1.0, description="Multiplier applied with the toxicity value (RBA, ADAF, CAF/MAF)")


# ------------------------------------------------------------------ basic equations
def intake(c: float, cr: float, cf: float, ef: float, ed: float, bw: float, at: float) -> float:
    """CDI (mg/kg-day) = C*CR*CF*EF*ED/(BW*AT). RAGS A Exhibit 6-9 and 6-14."""
    return c * cr * cf * ef * ed / (bw * at)


def da_event_soil(c_mg_kg: float, af_mg_cm2: float, absd: float) -> float:
    """DA_event (mg/cm2-event) = C_soil*1e-6 kg/mg*AF*ABS_d. RAGS E Eq. 3.12."""
    return c_mg_kg * 1e-6 * af_mg_cm2 * absd


def da_event_water(cw_mg_cm3: float, kp: float, t_event_h: float, inorganic: bool, fa: float | None = None,
                   b: float | None = None, tau_h: float | None = None, t_star_h: float | None = None) -> float:
    """DA_event (mg/cm2-event). RAGS E Eq. 3.2 and 3.3 (organics), Eq. 3.4 (inorganics)."""
    if inorganic:
        return kp * cw_mg_cm3 * t_event_h
    if t_event_h <= t_star_h:
        return 2.0 * fa * kp * cw_mg_cm3 * math.sqrt(6.0 * tau_h * t_event_h / math.pi)
    return fa * kp * cw_mg_cm3 * (t_event_h / (1.0 + b) + 2.0 * tau_h * (1.0 + 3.0 * b + 3.0 * b * b) / (1.0 + b) ** 2)


def dad(da_event: float, ev: float, ed: float, ef: float, sa: float, bw: float, at: float) -> float:
    """DAD (mg/kg-day) = DA_event*EV*ED*EF*SA/(BW*AT). RAGS E Eq. 3.1."""
    return da_event * ev * ed * ef * sa / (bw * at)


def exposure_concentration(ca_ug_m3: float, et_h_d: float, ef_d_yr: float, ed_yr: float, at_d: float) -> float:
    """EC (ug/m3) = CA*ET*EF*ED/AT with AT in hours (= days * 24). RAGS F Eq. 6 and 8."""
    return ca_ug_m3 * et_h_d * ef_d_yr * ed_yr / (at_d * 24.0)


# ------------------------------------------------------------------ risk and hazard from doses
def toxicity(c: Chemical, route: Route, cancer: bool, rfd: float | None = None, rfc: float | None = None):
    """Toxicity value matching a dose row; None if not available.

    Cancer: SFO (CDI), SFO/GIABS (DAD; RAGS E Eq. 4.2), IUR (EC, (ug/m3)^-1).
    Noncancer: RfD (CDI), RfD*GIABS (DAD; Eq. 4.3), RfC*1000 (EC, ug/m3; RAGS F Eq. 12).
    rfd/rfc replace the chronic values (e.g. subchronic, construction worker).
    """
    rfd = rfd if rfd is not None else c.rfd_mg_kg_day
    rfc = rfc if rfc is not None else c.rfc_mg_m3
    if cancer:
        if route == "inhalation":
            return c.iur_per_ug_m3
        if not c.sfo_per_mg_kg_day:
            return None
        return c.sfo_per_mg_kg_day / c.giabs if route == "dermal" else c.sfo_per_mg_kg_day
    if route == "inhalation":
        return rfc * 1000.0 if rfc else None
    if not rfd:
        return None
    return rfd * c.giabs if route == "dermal" else rfd


def risk_by_route(doses: list[Dose], c: Chemical) -> dict[str, float]:
    """Linear risk per route = sum(dose*weight*toxicity) (RAGS A p. 8-6; RAGS E; RAGS F Eq. 11)."""
    out: dict[str, float] = {}
    for d in doses:
        if d.basis != "cancer":
            continue
        t = toxicity(c, d.route, True)
        if t:
            out[d.route] = out.get(d.route, 0.0) + d.value * d.weight * t
    return out


def hq_by_route(doses: list[Dose], c: Chemical, basis: Basis, rfd: float | None = None,
                rfc: float | None = None) -> dict[str, float]:
    """HQ per route = sum(dose*weight)/toxicity (RAGS A p. 8-11; RAGS E Eq. 4.3; RAGS F Eq. 12)."""
    out: dict[str, float] = {}
    for d in doses:
        if d.basis != basis:
            continue
        t = toxicity(c, d.route, False, rfd, rfc)
        if t:
            out[d.route] = out.get(d.route, 0.0) + d.value * d.weight / t
    return out


# ------------------------------------------------------------------ helpers
_MGKGD, _UGM3 = "mg/kg-day", "ug/m3"


def _cdi(rows, basis, seg, value, w=1.0):
    rows.append(Dose(route="ingestion", quantity="CDI", basis=basis, segment=seg, value=value, unit=_MGKGD, weight=w))


def _dad(rows, basis, seg, value, w=1.0):
    rows.append(Dose(route="dermal", quantity="DAD", basis=basis, segment=seg, value=value, unit=_MGKGD, weight=w))


def _ec(rows, basis, seg, value, w=1.0):
    rows.append(Dose(route="inhalation", quantity="EC", basis=basis, segment=seg, value=value, unit=_UGM3, weight=w))


def _ages(d: Defaults) -> dict[str, dict[str, float]]:
    """Child and adult scalar parameters of the resident (or recreator) structure."""
    p = d.p
    return {
        "child": dict(ef=p["ef_child_d_yr"], ed=p["ed_child_yr"], et=p["et_child_h_d"], bw=p["bw_child_kg"],
                      irs=p["irs_child_mg_d"], sa=p["sa_child_cm2_d"], af=p["af_child_mg_cm2"],
                      irw=p["irw_child_l_d"], ev=p["ev_child_per_d"], saw=p["sa_water_child_cm2"],
                      et_ev=p["et_event_child_h"]),
        "adult": dict(ef=p["ef_adult_d_yr"], ed=p["ed_adult_yr"], et=p["et_adult_h_d"], bw=p["bw_adult_kg"],
                      irs=p["irs_adult_mg_d"], sa=p["sa_adult_cm2_d"], af=p["af_adult_mg_cm2"],
                      irw=p["irw_adult_l_d"], ev=p["ev_adult_per_d"], saw=p["sa_water_adult_cm2"],
                      et_ev=p["et_event_adult_h"]),
    }


def _cancer_weights(method: str, d: Defaults) -> tuple[float, float]:
    """(weight of the child/adult rows, factor applied to the age-segment ADAF rows)."""
    if method == "trichloroethylene":
        return d.tce["caf_oral"], d.tce["maf_oral"]
    return 1.0, 1.0


# ------------------------------------------------------------------ resident and recreator soil
def resident_soil(conc: float, c: Chemical, rba: float, vf: float | None, d: Defaults, method: str,
                  use_pef: bool = True) -> list[Dose]:
    """Soil ingestion (CDI), dermal contact (DAD), inhalation of vapors and particulates (EC).

    Same exposure terms as modules/risk/equations/soil.py (RSL UG 4.1.1).
    Air concentration CA = C*(1/VF + 1/PEF)*1000 ug/m3 (RSL UG 4.1.1.1).
    """
    p, rows = d.p, []
    ages = _ages(d)
    inv_air = (1.0 / vf if vf else 0.0) + (1.0 / p["pef_m3_kg"] if use_pef else 0.0)
    ca = conc * inv_air * 1000.0
    absd = c.absd
    # noncancer: child and adult separately, AT = 365*ED
    for who, a in ages.items():
        ed = a["ed"] if who == "child" else p["ed_res_yr"]       # adult noncancer uses EDres (RSL UG 4.1.1.2)
        at = 365.0 * ed
        b = f"noncancer-{who}"
        _cdi(rows, b, who, intake(conc, a["irs"], 1e-6, a["ef"], ed, a["bw"], at), rba)
        if absd:
            _dad(rows, b, who, dad(da_event_soil(conc, a["af"], absd), 1.0, ed, a["ef"], a["sa"], a["bw"], at))
        if ca:
            _ec(rows, b, who, exposure_concentration(ca, a["et"], a["ef"], ed, at))
    # cancer: AT = 365*LT
    at = 365.0 * p["lt_yr"]
    w_std, w_seg = _cancer_weights(method, d)
    if method in ("standard", "vinyl chloride", "trichloroethylene"):
        for who, a in ages.items():
            _cdi(rows, "cancer", who, intake(conc, a["irs"], 1e-6, a["ef"], a["ed"], a["bw"], at), rba * w_std)
            if absd:
                _dad(rows, "cancer", who, dad(da_event_soil(conc, a["af"], absd), 1.0, a["ed"], a["ef"], a["sa"],
                                              a["bw"], at), w_std)
    if method in ("mutagen", "trichloroethylene"):
        for g in d.segments:
            w = g["adaf"] * w_seg
            _cdi(rows, "cancer", g["label"], intake(conc, g["irs_mg_d"], 1e-6, g["ef_d_yr"], g["ed_yr"], g["bw_kg"], at),
                 rba * w)
            if absd:
                _dad(rows, "cancer", g["label"], dad(da_event_soil(conc, g["af_mg_cm2"], absd), 1.0, g["ed_yr"],
                                                     g["ef_d_yr"], g["sa_soil_cm2_d"], g["bw_kg"], at), w)
    if method == "vinyl chloride":
        # early-life term without averaging time: C*1e-6*IRS_c/BW_c (RSL UG 4.1.1.5)
        a = ages["child"]
        _cdi(rows, "cancer", "early-life", conc * 1e-6 * a["irs"] / a["bw"], rba)
        if absd:
            _dad(rows, "cancer", "early-life", da_event_soil(conc, a["af"], absd) * a["sa"] / a["bw"])
    # cancer inhalation
    if method == "vinyl chloride":
        if vf:      # only the VF term (RSL UG 4.1.1.5, p. 49)
            ca_v = conc / vf * 1000.0
            _ec(rows, "cancer", "resident", exposure_concentration(ca_v, p["et_res_h_d"], p["ef_res_d_yr"],
                                                                  p["ed_res_yr"], at))
            _ec(rows, "cancer", "early-life", ca_v)
    elif ca:
        if method in ("standard", "trichloroethylene"):
            w = d.tce["caf_inh"] if method == "trichloroethylene" else 1.0
            _ec(rows, "cancer", "resident", exposure_concentration(ca, p["et_res_h_d"], p["ef_res_d_yr"],
                                                                  p["ed_res_yr"], at), w)
        if method in ("mutagen", "trichloroethylene"):
            f = d.tce["maf_inh"] if method == "trichloroethylene" else 1.0
            for g in d.segments:
                _ec(rows, "cancer", g["label"], exposure_concentration(ca, g["et_h_d"], g["ef_d_yr"], g["ed_yr"], at),
                    g["adaf"] * f)
    return rows


# ------------------------------------------------------------------ resident tap water, recreator surface water
def resident_water(conc_ug_l: float, c: Chemical, d: Defaults, f: Factors, method: str,
                   inhalation: bool = True) -> list[Dose]:
    """Water ingestion (CDI), dermal contact while bathing (DAD), inhalation of volatiles (EC).

    Same exposure terms as modules/risk/equations/tapwater.py (RSL UG 4.1.2).
    Cw (mg/cm3) = C (ug/L) * 1e-6. CA = C*K (Andelman, K = 0.5 L/m3; RSL UG 4.1.2).
    For the cancer dermal dose the event time is the ED-weighted ET_event-adj (standard,
    VC) or ET_event-madj (mutagen, TCE), as in RSL UG 4.1.2.3-4.1.2.6.
    """
    p, rows = d.p, []
    ages = _ages(d)
    applies, inorganic, fa, _ = dw.applicability(c)
    cw = conc_ug_l * 1e-6
    ca = conc_ug_l * p["k_andelman_l_m3"] if (inhalation and c.volatile) else 0.0

    def da(t_event):
        return da_event_water(cw, c.kp_cm_h, t_event, inorganic, fa, c.b_dimensionless, c.tau_event_h, c.t_star_h)

    for who, a in ages.items():
        ed = a["ed"] if who == "child" else p["ed_res_yr"]
        at = 365.0 * ed
        b = f"noncancer-{who}"
        _cdi(rows, b, who, intake(conc_ug_l, a["irw"], 1e-3, a["ef"], ed, a["bw"], at))
        if applies:
            _dad(rows, b, who, dad(da(a["et_ev"]), a["ev"], ed, a["ef"], a["saw"], a["bw"], at))
        if ca:
            _ec(rows, b, who, exposure_concentration(ca, a["et"], a["ef"], ed, at))

    at = 365.0 * p["lt_yr"]
    w_std, w_seg = _cancer_weights(method, d)
    et_std = f.et_event_adj_h if method in ("standard", "vinyl chloride") else f.et_event_madj_h
    if method in ("standard", "vinyl chloride", "trichloroethylene"):
        for who, a in ages.items():
            _cdi(rows, "cancer", who, intake(conc_ug_l, a["irw"], 1e-3, a["ef"], a["ed"], a["bw"], at), w_std)
            if applies:
                _dad(rows, "cancer", who, dad(da(et_std), a["ev"], a["ed"], a["ef"], a["saw"], a["bw"], at), w_std)
    if method in ("mutagen", "trichloroethylene"):
        for g in d.segments:
            w = g["adaf"] * w_seg
            _cdi(rows, "cancer", g["label"], intake(conc_ug_l, g["irw_l_d"], 1e-3, g["ef_d_yr"], g["ed_yr"], g["bw_kg"],
                                                    at), w)
            if applies:
                _dad(rows, "cancer", g["label"], dad(da(f.et_event_madj_h), g["ev_per_d"], g["ed_yr"], g["ef_d_yr"],
                                                     g["sa_water_cm2"], g["bw_kg"], at), w)
    if method == "vinyl chloride":
        a = ages["child"]
        _cdi(rows, "cancer", "early-life", conc_ug_l * 1e-3 * a["irw"] / a["bw"])
        if applies:
            _dad(rows, "cancer", "early-life", da(f.et_event_adj_h) * a["ev"] * a["saw"] / a["bw"])
    if ca:
        if method in ("standard", "vinyl chloride", "trichloroethylene"):
            w = d.tce["caf_inh"] if method == "trichloroethylene" else 1.0
            _ec(rows, "cancer", "resident", exposure_concentration(ca, p["et_res_h_d"], p["ef_res_d_yr"],
                                                                  p["ed_res_yr"], at), w)
        if method in ("mutagen", "trichloroethylene"):
            m = d.tce["maf_inh"] if method == "trichloroethylene" else 1.0
            for g in d.segments:
                _ec(rows, "cancer", g["label"], exposure_concentration(ca, g["et_h_d"], g["ef_d_yr"], g["ed_yr"], at),
                    g["adaf"] * m)
        if method == "vinyl chloride":
            _ec(rows, "cancer", "early-life", ca)
    return rows


# ------------------------------------------------------------------ resident air
def resident_air(conc_ug_m3: float, d: Defaults, method: str) -> list[Dose]:
    """EC for the resident (RSL UG 4.1.3; RAGS F Eq. 6 and 8)."""
    p, rows = d.p, []
    ef, ed, et = p["ef_res_d_yr"], p["ed_res_yr"], p["et_res_h_d"]
    _ec(rows, "noncancer", "resident", exposure_concentration(conc_ug_m3, et, ef, ed, 365.0 * ed))
    at = 365.0 * p["lt_yr"]
    if method in ("standard", "vinyl chloride", "trichloroethylene"):
        w = d.tce["caf_inh"] if method == "trichloroethylene" else 1.0
        _ec(rows, "cancer", "resident", exposure_concentration(conc_ug_m3, et, ef, ed, at), w)
    if method in ("mutagen", "trichloroethylene"):
        m = d.tce["maf_inh"] if method == "trichloroethylene" else 1.0
        for g in d.segments:
            _ec(rows, "cancer", g["label"], exposure_concentration(conc_ug_m3, g["et_h_d"], g["ef_d_yr"], g["ed_yr"],
                                                                  at), g["adaf"] * m)
    if method == "vinyl chloride":
        _ec(rows, "cancer", "early-life", conc_ug_m3)
    return rows


# ------------------------------------------------------------------ workers
def worker_soil(conc: float, absd: float | None, rba: float, vf: float | None, pef: float, w: Worker) -> list[Dose]:
    """Adult worker (composite, outdoor, indoor, construction): RSL UG 4.2-4.5.

    No ADAF, VC or TCE adjustments for adult workers (see equations/worker.py).
    """
    rows: list[Dose] = []
    ca = conc * ((1.0 / vf if vf else 0.0) + 1.0 / pef) * 1000.0
    for basis, at in (("noncancer", w.at_nc_d), ("cancer", 365.0 * w.lt_yr)):
        _cdi(rows, basis, "worker", intake(conc, w.irs_mg_d, 1e-6, w.ef_d_yr, w.ed_yr, w.bw_kg, at), rba)
        if absd and w.sa_cm2_d:
            _dad(rows, basis, "worker", dad(da_event_soil(conc, w.af_mg_cm2, absd), 1.0, w.ed_yr, w.ef_d_yr,
                                            w.sa_cm2_d, w.bw_kg, at))
        _ec(rows, basis, "worker", exposure_concentration(ca, w.et_h_d, w.ef_d_yr, w.ed_yr, at))
    return rows


def worker_air(conc_ug_m3: float, w: Worker) -> list[Dose]:
    """Worker EC (RSL UG 4.2.2; RAGS F Eq. 6 and 8)."""
    rows: list[Dose] = []
    for basis, at in (("noncancer", w.at_nc_d), ("cancer", 365.0 * w.lt_yr)):
        _ec(rows, basis, "worker", exposure_concentration(conc_ug_m3, w.et_h_d, w.ef_d_yr, w.ed_yr, at))
    return rows


# ------------------------------------------------------------------ fish
def fish(conc_mg_kg: float, irf_mg_d: float, ef: float, ed: float, bw: float, lt: float) -> list[Dose]:
    """Fish tissue ingestion (RSL UG 4.7): CDI = C*IRF*1e-6*EF*ED/(BW*AT)."""
    rows: list[Dose] = []
    _cdi(rows, "noncancer", "adult", intake(conc_mg_kg, irf_mg_d, 1e-6, ef, ed, bw, 365.0 * ed))
    _cdi(rows, "cancer", "adult", intake(conc_mg_kg, irf_mg_d, 1e-6, ef, ed, bw, 365.0 * lt))
    return rows
