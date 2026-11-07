"""Resident air screening levels, ug/m3 (RSL User's Guide 4.1.3, PDF pp. 96-101).

  Noncancer       4.1.3.1, p. 97:  SL = THQ*AT / ((1/RfC)*(1/1000)*EFres*EDres*ETres/24), AT = 365*EDres
  Cancer          4.1.3.2, p. 98:  SL = TR*AT / (IUR*EFres*EDres*ETres/24), AT = 365*LT
  Mutagenic       4.1.3.3, p. 99:  ADAF-weighted sum over age segments
  Vinyl chloride  4.1.3.4, p. 100: SL = TR / (IUR + IUR*EFres*EDres*ETres/24/AT)
  TCE             4.1.3.5, p. 101: CAFi and MAFi fractions
"""

from __future__ import annotations

from modules.risk.common import Defaults
from modules.risk.equations.factors import Factors, inh_tce_d
from ortam_core.chemicals.database import Chemical


def noncancer(c: Chemical, d: Defaults, thq: float) -> float | None:
    p = d.p
    if not c.rfc_mg_m3:
        return None
    at = 365.0 * p["ed_res_yr"]
    return thq * at / ((1.0 / c.rfc_mg_m3) / 1000.0 * p["ef_res_d_yr"] * p["ed_res_yr"] * p["et_res_h_d"] / 24.0)


def cancer(c: Chemical, d: Defaults, f: Factors, tr: float, method: str) -> float | None:
    iur = c.iur_per_ug_m3
    if not iur:
        return None
    at = 365.0 * d.p["lt_yr"]
    if method == "standard":
        return tr * at / (iur * f.inh_d)
    if method == "mutagen":
        return tr * at / (iur * f.inh_m_d)
    if method == "vinyl chloride":
        return tr / (iur + iur * f.inh_d / at)
    if method == "trichloroethylene":
        return tr * at / (iur * inh_tce_d(d, f))
    raise ValueError(f"Unknown cancer method: {method}")
