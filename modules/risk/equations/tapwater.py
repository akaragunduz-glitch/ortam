"""Resident tap water screening levels, ug/L (RSL User's Guide 4.1.2).

  Noncancer child          4.1.2.1, PDF pp. 63-66
  Noncancer adult          4.1.2.2, PDF pp. 66-70
  Cancer (age-adjusted)    4.1.2.3, PDF pp. 70-75
  Mutagenic                4.1.2.4, PDF pp. 75-81
  Vinyl chloride           4.1.2.5, PDF pp. 81-84
  Trichloroethylene        4.1.2.6, PDF pp. 84-90
Inhalation of volatiles during household use is calculated only for volatile
chemicals (4.1.2, p. 62), with the Andelman factor K (L/m3).
Dermal: see dermal_water.py.
"""

from __future__ import annotations

from modules.risk.common import Defaults
from modules.risk.equations import dermal_water as dw
from modules.risk.equations.factors import Factors, inh_tce_d
from ortam_core.chemicals.database import Chemical

Triple = tuple[float | None, float | None, float | None]


def noncancer(c: Chemical, d: Defaults, age: str, thq: float) -> Triple:
    p = d.p
    if age == "child":
        ef, ed, et, bw = p["ef_child_d_yr"], p["ed_child_yr"], p["et_child_h_d"], p["bw_child_kg"]
        irw, ev, sa, et_ev = p["irw_child_l_d"], p["ev_child_per_d"], p["sa_water_child_cm2"], p["et_event_child_h"]
    else:
        ef, ed, et, bw = p["ef_adult_d_yr"], p["ed_res_yr"], p["et_adult_h_d"], p["bw_adult_kg"]
        irw, ev, sa, et_ev = p["irw_adult_l_d"], p["ev_adult_per_d"], p["sa_water_adult_cm2"], p["et_event_adult_h"]
    at = 365.0 * ed
    ing = derm = inh = None
    if c.rfd_mg_kg_day:
        # SL_ing = THQ*AT*BW / ((1/RfD) * (1 mg/1000 ug) * EF * ED * IRW)
        ing = thq * at * bw / ((1.0 / c.rfd_mg_kg_day) / 1000.0 * ef * ed * irw)
        applies, inorganic, fa, _ = dw.applicability(c)
        if applies:
            # DA_event = THQ*AT*BW / ((1/(RfD*GIABS)) * (1/1000) * EF * ED * EV * SA)
            da = thq * at * bw / ((1.0 / (c.rfd_mg_kg_day * c.giabs)) / 1000.0 * ef * ed * ev * sa)
            derm = dw.sl_from_da_event(da, c, et_ev, inorganic, fa)
    if c.rfc_mg_m3 and c.volatile:
        # SL_inh = THQ*AT / ((1/RfC) * (1/1000) * EF * ED * ET * (1 day/24 h) * K)
        inh = thq * at / ((1.0 / c.rfc_mg_m3) / 1000.0 * ef * ed * et / 24.0 * p["k_andelman_l_m3"])
    return ing, derm, inh


def cancer(c: Chemical, d: Defaults, f: Factors, tr: float, method: str) -> Triple:
    p = d.p
    at = 365.0 * p["lt_yr"]
    k = p["k_andelman_l_m3"]
    sfo, iur = c.sfo_per_mg_kg_day, c.iur_per_ug_m3
    applies, inorganic, fa, _ = dw.applicability(c)
    ing = derm = inh = None

    if method in ("standard", "mutagen"):
        ifw = f.ifw_l_kg if method == "standard" else f.ifwm_l_kg
        dfw = f.dfw_cm2_event_kg if method == "standard" else f.dfwm_cm2_event_kg
        et_ev = f.et_event_adj_h if method == "standard" else f.et_event_madj_h
        inh_days = f.inh_d if method == "standard" else f.inh_m_d
        if sfo:
            # SL_ing = TR*AT / (SFO * (1/1000) * IFW)          (4.1.2.3 p. 71; 4.1.2.4 p. 76)
            ing = tr * at / (sfo / 1000.0 * ifw)
            if applies:
                # DA_event = TR*AT / ((SFO/GIABS) * (1/1000) * DFW)  (p. 74; p. 78)
                da = tr * at / ((sfo / c.giabs) / 1000.0 * dfw)
                derm = dw.sl_from_da_event(da, c, et_ev, inorganic, fa)
        if iur and c.volatile:
            # SL_inh = TR*AT / (IUR * EFres*EDres*ETres/24 * K)  (p. 74; ADAF sum p. 80)
            inh = tr * at / (iur * inh_days * k)

    elif method == "vinyl chloride":
        # 4.1.2.5 (pp. 81-84): lifetime term plus early-life term without averaging time
        if sfo:
            ing = tr / (sfo / 1000.0 * f.ifw_l_kg / at + sfo / 1000.0 * p["irw_child_l_d"] / p["bw_child_kg"])
            if applies:
                da = tr / ((sfo / c.giabs) / 1000.0 * (f.dfw_cm2_event_kg / at
                           + p["ev_child_per_d"] * p["sa_water_child_cm2"] / p["bw_child_kg"]))
                derm = dw.sl_from_da_event(da, c, f.et_event_adj_h, inorganic, fa)
        if iur and c.volatile:
            inh = tr / (iur * f.inh_d * k / at + iur * k)

    elif method == "trichloroethylene":
        # 4.1.2.6 (pp. 84-90)
        caf, maf = d.tce["caf_oral"], d.tce["maf_oral"]
        if sfo:
            ing = tr * at / (sfo / 1000.0 * (caf * f.ifw_l_kg + maf * f.ifwm_l_kg))
            if applies:
                da = tr * at / ((sfo / c.giabs) / 1000.0 * (caf * f.dfw_cm2_event_kg + maf * f.dfwm_cm2_event_kg))
                derm = dw.sl_from_da_event(da, c, f.et_event_madj_h, inorganic, fa)
        if iur and c.volatile:
            inh = tr * at / (iur * k * inh_tce_d(d, f))
    else:
        raise ValueError(f"Unknown cancer method: {method}")
    return ing, derm, inh
