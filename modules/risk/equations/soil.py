"""Resident soil screening levels, mg/kg (RSL User's Guide 4.1.1).

  Noncancer child          4.1.1.1, PDF pp. 36-37
  Noncancer adult          4.1.1.2, PDF pp. 38-39
  Cancer (age-adjusted)    4.1.1.3, PDF pp. 40-41
  Mutagenic                4.1.1.4, PDF pp. 42-45
  Vinyl chloride           4.1.1.5, PDF pp. 46-49
  Trichloroethylene        4.1.1.6, PDF pp. 50-56
"""

from __future__ import annotations

from modules.risk.common import Defaults
from modules.risk.equations.factors import Factors, inh_tce_d
from modules.risk.toxicity import Chemical

Triple = tuple[float | None, float | None, float | None]


def _inv_vf_pef(vf: float | None, pef: float, use_pef: bool = True) -> float:
    """(1/VF + 1/PEF), kg/m3; the VF term is omitted when no VF applies."""
    return (1.0 / vf if vf else 0.0) + (1.0 / pef if use_pef else 0.0)


def noncancer(c: Chemical, rba: float, vf: float | None, d: Defaults, age: str, thq: float) -> Triple:
    """Child: ED = EDres-c, AT = 365*EDres-c. Adult: ED = EDres, AT = 365*EDres."""
    p = d.p
    if age == "child":
        ef, ed, et = p["ef_child_d_yr"], p["ed_child_yr"], p["et_child_h_d"]
        bw, irs, sa, af = p["bw_child_kg"], p["irs_child_mg_d"], p["sa_child_cm2_d"], p["af_child_mg_cm2"]
    else:
        ef, ed, et = p["ef_adult_d_yr"], p["ed_res_yr"], p["et_adult_h_d"]
        bw, irs, sa, af = p["bw_adult_kg"], p["irs_adult_mg_d"], p["sa_adult_cm2_d"], p["af_adult_mg_cm2"]
    at = 365.0 * ed
    ing = derm = inh = None
    if c.rfd_mg_kg_day:
        # SL_ing = THQ*AT*BW / ((RBA/RfD) * 1e-6 kg/mg * EF * ED * IRS)
        ing = thq * at * bw / ((rba / c.rfd_mg_kg_day) * 1e-6 * ef * ed * irs)
        if c.absd:
            # SL_derm = THQ*AT*BW / ((1/(RfD*GIABS)) * 1e-6 * EF * ED * SA * AF * ABSd)
            derm = thq * at * bw / ((1.0 / (c.rfd_mg_kg_day * c.giabs)) * 1e-6 * ef * ed * sa * af * c.absd)
    if c.rfc_mg_m3:
        # SL_inh = THQ*AT / ((1/RfC) * EF * ED * ET * (1 day/24 h) * (1/VF + 1/PEF))
        inh = thq * at / ((1.0 / c.rfc_mg_m3) * ef * ed * et / 24.0 * _inv_vf_pef(vf, p["pef_m3_kg"]))
    return ing, derm, inh


def cancer(c: Chemical, rba: float, vf: float | None, d: Defaults, f: Factors, tr: float,
           method: str) -> Triple:
    p = d.p
    at = 365.0 * p["lt_yr"]
    sfo, iur, absd = c.sfo_per_mg_kg_day, c.iur_per_ug_m3, c.absd
    inv_air = _inv_vf_pef(vf, p["pef_m3_kg"])
    ing = derm = inh = None

    if method == "standard":
        # 4.1.1.3: SL_ing = TR*AT/(SFO*1e-6*RBA*IFS); SL_derm = TR*AT/((SFO/GIABS)*1e-6*DFS*ABSd)
        # SL_inh = TR*AT/(IUR*1000*EFres*EDres*ETres/24*(1/VF+1/PEF))
        if sfo:
            ing = tr * at / (sfo * 1e-6 * rba * f.ifs_mg_kg)
            if absd:
                derm = tr * at / ((sfo / c.giabs) * 1e-6 * f.dfs_mg_kg * absd)
        if iur:
            inh = tr * at / (iur * 1000.0 * f.inh_d * inv_air)

    elif method == "mutagen":
        # 4.1.1.4: IFSM, DFSM and the ADAF-weighted inhalation sum replace the standard factors
        if sfo:
            ing = tr * at / (sfo * 1e-6 * rba * f.ifsm_mg_kg)
            if absd:
                derm = tr * at / ((sfo / c.giabs) * 1e-6 * f.dfsm_mg_kg * absd)
        if iur:
            inh = tr * at / (iur * 1000.0 * inv_air * f.inh_m_d)

    elif method == "vinyl chloride":
        # 4.1.1.5: adult-lifetime term plus an early-life term without averaging time
        if sfo:
            ing = tr / (sfo * 1e-6 * rba * f.ifs_mg_kg / at
                        + sfo * 1e-6 * rba * p["irs_child_mg_d"] / p["bw_child_kg"])
            if absd:
                derm = tr / ((sfo / c.giabs) * 1e-6 * f.dfs_mg_kg * absd / at
                             + (sfo / c.giabs) * 1e-6 * p["sa_child_cm2_d"] * p["af_child_mg_cm2"] * absd
                             / p["bw_child_kg"])
        if iur and vf:
            # p. 49: only the VF term appears (no PEF)
            inh = tr / (iur * 1000.0 * f.inh_d / (at * vf) + iur * 1000.0 / vf)

    elif method == "trichloroethylene":
        # 4.1.1.6: SL_ing = TR*AT/(SFO*1e-6*RBA*(CAFo*IFS + MAFo*IFSM)), etc.
        caf, maf = d.tce["caf_oral"], d.tce["maf_oral"]
        if sfo:
            ing = tr * at / (sfo * 1e-6 * rba * (caf * f.ifs_mg_kg + maf * f.ifsm_mg_kg))
            if absd:
                derm = tr * at / ((sfo / c.giabs) * 1e-6 * (caf * f.dfs_mg_kg * absd + maf * f.dfsm_mg_kg * absd))
        if iur:
            inh = tr * at / (iur * 1000.0 * inv_air * inh_tce_d(d, f))
    else:
        raise ValueError(f"Unknown cancer method: {method}")
    return ing, derm, inh
