"""Adult worker receptors: composite, outdoor, indoor and construction worker.

All four use the same equation form with different parameters
(RSL User's Guide 4.2-4.5, PDF pp. 101-124). Mutagenic, VC and TCE
adjustments do not apply to adult workers (verified against the RSL
Industrial Soil and Industrial Air tables).

Soil, mg/kg (example: composite, 4.2.1.1-4.2.1.2, pp. 101-104):
  nc ingestion  SL = THQ*ATnc*BW / ((RBA/RfD)*1e-6*EF*ED*IRS)
  nc dermal     SL = THQ*ATnc*BW / ((1/(RfD*GIABS))*1e-6*EF*ED*SA*AF*ABSd)
  nc inhal.     SL = THQ*ATnc / ((1/RfC)*EF*ED*ET/24*(1/VF + 1/PEF))
  ca ingestion  SL = TR*AT*BW / (SFO*1e-6*RBA*EF*ED*IRS)
  ca dermal     SL = TR*AT*BW / ((SFO/GIABS)*1e-6*EF*ED*SA*AF*ABSd)
  ca inhal.     SL = TR*AT / (IUR*1000*EF*ED*ET/24*(1/VF + 1/PEF))
with ATnc = 365*ED (composite, outdoor, indoor) or EW*7*ED (construction,
4.5.1.1, p. 117) and AT = 365*LT. Construction: EF = EW*DW, PEF -> PEFsc or
PEF'sc, VF -> VFsc (4.5.1-4.5.2).
Air, ug/m3 (4.2.2, pp. 105-106):
  nc  SL = THQ*ATnc / ((1/RfC)*(1/1000)*EF*ED*ET/24);  ca  SL = TR*AT / (IUR*EF*ED*ET/24)
"""

from __future__ import annotations

from dataclasses import dataclass

Triple = tuple[float | None, float | None, float | None]


@dataclass(frozen=True)
class Worker:
    ef_d_yr: float
    ed_yr: float
    et_h_d: float
    bw_kg: float
    irs_mg_d: float
    sa_cm2_d: float | None      # None: no dermal route (indoor worker)
    af_mg_cm2: float | None
    at_nc_d: float              # noncancer averaging time, days
    lt_yr: float = 70.0


def _inv_air(vf: float | None, pef: float) -> float:
    return (1.0 / vf if vf else 0.0) + 1.0 / pef


def soil_noncancer(rfd, rfc, giabs, absd, rba, vf, pef, w: Worker, thq: float) -> Triple:
    ing = derm = inh = None
    if rfd:
        ing = thq * w.at_nc_d * w.bw_kg / ((rba / rfd) * 1e-6 * w.ef_d_yr * w.ed_yr * w.irs_mg_d)
        if absd and w.sa_cm2_d:
            derm = thq * w.at_nc_d * w.bw_kg / ((1.0 / (rfd * giabs)) * 1e-6 * w.ef_d_yr * w.ed_yr
                                                 * w.sa_cm2_d * w.af_mg_cm2 * absd)
    if rfc:
        inh = thq * w.at_nc_d / ((1.0 / rfc) * w.ef_d_yr * w.ed_yr * w.et_h_d / 24.0 * _inv_air(vf, pef))
    return ing, derm, inh


def soil_cancer(sfo, iur, giabs, absd, rba, vf, pef, w: Worker, tr: float) -> Triple:
    at = 365.0 * w.lt_yr
    ing = derm = inh = None
    if sfo:
        ing = tr * at * w.bw_kg / (sfo * 1e-6 * rba * w.ef_d_yr * w.ed_yr * w.irs_mg_d)
        if absd and w.sa_cm2_d:
            derm = tr * at * w.bw_kg / ((sfo / giabs) * 1e-6 * w.ef_d_yr * w.ed_yr * w.sa_cm2_d
                                        * w.af_mg_cm2 * absd)
    if iur:
        inh = tr * at / (iur * 1000.0 * w.ef_d_yr * w.ed_yr * w.et_h_d / 24.0 * _inv_air(vf, pef))
    return ing, derm, inh


def air_noncancer(rfc, w: Worker, thq: float) -> float | None:
    return thq * w.at_nc_d / ((1.0 / rfc) / 1000.0 * w.ef_d_yr * w.ed_yr * w.et_h_d / 24.0) if rfc else None


def air_cancer(iur, w: Worker, tr: float) -> float | None:
    return tr * 365.0 * w.lt_yr / (iur * w.ef_d_yr * w.ed_yr * w.et_h_d / 24.0) if iur else None
