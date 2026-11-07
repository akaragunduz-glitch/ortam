"""Age-adjusted intake and contact factors.

Sources (RSL User's Guide, PDF pages):
  IFSres-adj, DFSres-adj         4.1.1.3, pp. 40-41
  IFSMres-adj, DFSMres-adj       4.1.1.4, pp. 42-44
  IFWres-adj, DFWres-adj         4.1.2.3, pp. 71-74
  IFWMres-adj, DFWMres-adj       4.1.2.4, pp. 76-78
  ETevent-res-adj, -madj         4.1.2.3, p. 72; 4.1.2.4, p. 78
  Mutagenic inhalation sum       4.1.1.4, p. 44; 4.1.3.3, p. 99
"""

from __future__ import annotations

from dataclasses import dataclass

from modules.risk.common import Defaults


@dataclass(frozen=True)
class Factors:
    ifs_mg_kg: float
    dfs_mg_kg: float
    ifsm_mg_kg: float
    dfsm_mg_kg: float
    ifw_l_kg: float
    ifwm_l_kg: float
    dfw_cm2_event_kg: float
    dfwm_cm2_event_kg: float
    et_event_adj_h: float
    et_event_madj_h: float
    inh_d: float      # EFres*EDres*ETres*(1 day/24 h), days
    inh_m_d: float    # sum over age segments of EF*ED*ET/24*ADAF, days


def compute(d: Defaults) -> Factors:
    p, s = d.p, d.segments
    efc, eda, efa, edc = p["ef_child_d_yr"], p["ed_adult_yr"], p["ef_adult_d_yr"], p["ed_child_yr"]
    bwc, bwa = p["bw_child_kg"], p["bw_adult_kg"]
    return Factors(
        # IFS = EFc*EDc*IRSc/BWc + EFa*EDa*IRSa/BWa
        ifs_mg_kg=efc * edc * p["irs_child_mg_d"] / bwc + efa * eda * p["irs_adult_mg_d"] / bwa,
        # DFS = EFc*EDc*SAc*AFc/BWc + EFa*EDa*SAa*AFa/BWa
        dfs_mg_kg=(efc * edc * p["sa_child_cm2_d"] * p["af_child_mg_cm2"] / bwc
                   + efa * eda * p["sa_adult_cm2_d"] * p["af_adult_mg_cm2"] / bwa),
        # IFSM = sum_i EF_i*ED_i*IRS_i*ADAF_i/BW_i
        ifsm_mg_kg=sum(g["ef_d_yr"] * g["ed_yr"] * g["irs_mg_d"] * g["adaf"] / g["bw_kg"] for g in s),
        # DFSM = sum_i EF_i*ED_i*SA_i*AF_i*ADAF_i/BW_i
        dfsm_mg_kg=sum(g["ef_d_yr"] * g["ed_yr"] * g["sa_soil_cm2_d"] * g["af_mg_cm2"] * g["adaf"] / g["bw_kg"]
                       for g in s),
        # IFW = EFc*EDc*IRWc/BWc + EFa*EDa*IRWa/BWa
        ifw_l_kg=efc * edc * p["irw_child_l_d"] / bwc + efa * eda * p["irw_adult_l_d"] / bwa,
        # IFWM = sum_i EF_i*ED_i*IRW_i*ADAF_i/BW_i
        ifwm_l_kg=sum(g["ef_d_yr"] * g["ed_yr"] * g["irw_l_d"] * g["adaf"] / g["bw_kg"] for g in s),
        # DFW = EFc*EDc*EVc*SAc/BWc + EFa*EDa*EVa*SAa/BWa
        dfw_cm2_event_kg=(efc * edc * p["ev_child_per_d"] * p["sa_water_child_cm2"] / bwc
                          + efa * eda * p["ev_adult_per_d"] * p["sa_water_adult_cm2"] / bwa),
        # DFWM = sum_i EF_i*ED_i*EV_i*SA_i*ADAF_i/BW_i
        dfwm_cm2_event_kg=sum(g["ef_d_yr"] * g["ed_yr"] * g["ev_per_d"] * g["sa_water_cm2"] * g["adaf"]
                              / g["bw_kg"] for g in s),
        # ETevent-adj = (EDc*ETevent-c + EDa*ETevent-a) / EDres
        et_event_adj_h=(edc * p["et_event_child_h"] + eda * p["et_event_adult_h"]) / p["ed_res_yr"],
        # ETevent-madj = sum_i ED_i*ETevent_i / sum_i ED_i
        et_event_madj_h=sum(g["ed_yr"] * g["et_event_h"] for g in s) / sum(g["ed_yr"] for g in s),
        inh_d=p["ef_res_d_yr"] * p["ed_res_yr"] * p["et_res_h_d"] / 24.0,
        inh_m_d=sum(g["ef_d_yr"] * g["ed_yr"] * g["et_h_d"] / 24.0 * g["adaf"] for g in s),
    )


def inh_tce_d(d: Defaults, f: Factors) -> float:
    """TCE inhalation exposure term, days: EFres*EDres*ETres/24*CAFi + sum_i(...)*MAFi*ADAF_i.

    RSL UG 4.1.1.6, p. 54; 4.1.2.6, p. 90; 4.1.3.5, p. 101.
    """
    return f.inh_d * d.tce["caf_inh"] + d.tce["maf_inh"] * f.inh_m_d
