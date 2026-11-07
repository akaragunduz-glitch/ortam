"""Dermal absorption from water (tap water).

Sources:
  - Effective predictive domain (EPD): RAGS Part E (2004), p. 3-7, Eq. 3.9 and 3.10.
    (The RSL User's Guide 4.9.8, PDF p. 217, reprints these with a sign error in
    Eq. 3.9 and the wrong lower bound in Eq. 3.10; the RAGS Part E form is used and
    reproduces the RSL tap water table - see docs/methods/risk.md.)
  - Fraction absorbed FA: RSL UG 4.9.8, PDF p. 221 (regression in log B and log tau,
    clipped to [0, 1], rounded to the nearest tenth).
  - SL from DA_event, inorganic and organic (ET <= t* and ET > t*):
    RSL UG 4.1.2.1, PDF p. 64 (same form in 4.1.2.2-4.1.2.6).
Inorganic = no log Kow in the RSL Chemical Parameters table (Kp from RAGS Part E
Exhibit 3-1, default 0.001 cm/h).
"""

from __future__ import annotations

import math

from modules.risk.toxicity import Chemical


def in_epd(mw: float, log_kow: float) -> bool:
    """RAGS Part E (2004), p. 3-7, Eq. 3.9 and 3.10."""
    eq39 = 5.103e-4 * mw + 0.05616 * log_kow
    eq310 = -5.103e-4 * mw + 0.05616 * log_kow
    return -0.06831 <= eq39 <= 0.5577 and -0.3010 <= eq310 <= 0.1758


def fraction_absorbed(b: float, tau_event_h: float) -> float:
    """RSL UG 4.9.8, PDF p. 221."""
    lb, lt = math.log10(b), math.log10(tau_event_h)
    if b <= 0.1:
        fa = 0.9589849087 - 0.0163393790 * lb - 0.1451565908 * lt - 0.0534664095 * lb * lt
    elif b <= 1.0:
        fa = 1.051232292 + 0.091016187 * lb - 0.286735467 * lt - 0.180504367 * lb * lt
    else:
        fa = 0.992336792 + 0.479643809 * lb - 0.114381522 * lt - 1.263647642 * lb * lt
    if fa >= 1.0:
        return 1.0
    if fa < 0.0:
        return 0.0
    return round(fa, 1)


def applicability(c: Chemical) -> tuple[bool, bool, float | None, str | None]:
    """Is the dermal route calculated? Returns (applies, inorganic, FA, note)."""
    if c.kp_cm_h is None:
        return False, False, None, "dermal (water) not calculated: no Kp"
    if c.log_kow is None:
        return True, True, None, None
    if c.mw_g_mol is None or not in_epd(c.mw_g_mol, c.log_kow):
        return False, False, None, "dermal (water) not calculated: outside the EPD (RAGS E Eq. 3.9-3.10)"
    if None in (c.b_dimensionless, c.tau_event_h, c.t_star_h):
        return False, False, None, "dermal (water) not calculated: B, tau_event or t* missing"
    fa = fraction_absorbed(c.b_dimensionless, c.tau_event_h)
    if fa == 0.0:
        return False, False, fa, "dermal (water) not calculated: FA = 0"
    return True, False, fa, None


def sl_from_da_event(da_ug_cm2_event: float, c: Chemical, et_event_h: float,
                     inorganic: bool, fa: float | None) -> float:
    """SL (ug/L) from the dermally absorbed dose per event DA_event (ug/cm2-event).

    RSL UG 4.1.2.1, PDF p. 64:
      inorganic:      SL = DA*1000 / (Kp*ET)
      organic ET<=t*: SL = DA*1000 / (2*FA*Kp*sqrt(6*tau*ET/pi))
      organic ET>t*:  SL = DA*1000 / (FA*Kp*(ET/(1+B) + 2*tau*(1+3B+3B^2)/(1+B)^2))
    """
    kp = c.kp_cm_h
    if inorganic:
        return da_ug_cm2_event * 1000.0 / (kp * et_event_h)
    b, tau = c.b_dimensionless, c.tau_event_h
    if et_event_h <= c.t_star_h:
        return da_ug_cm2_event * 1000.0 / (2.0 * fa * kp * math.sqrt(6.0 * tau * et_event_h / math.pi))
    return da_ug_cm2_event * 1000.0 / (
        fa * kp * (et_event_h / (1.0 + b) + 2.0 * tau * (1.0 + 3.0 * b + 3.0 * b * b) / (1.0 + b) ** 2))
