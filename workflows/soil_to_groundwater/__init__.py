"""Soil-to-groundwater workflows: risk-based SSL (RSL) and VLEACH leachate.

Connects two standalone tools (modules never call each other):
  1. risk.resident-tapwater gives the unrounded risk-based tap water SL (ug/L);
  2. vadose.ssl-leaching converts it into a soil SSL (mg/kg).
This reproduces the RSL 'Risk-based SSL' column (DAF = 1).
"""

from __future__ import annotations

from modules.risk.resident_tapwater import ResidentTapwaterModule
from modules.vadose.ssl_leaching import SslLeachingModule


def risk_based_ssl(chemicals: list[dict], ssl_options: dict | None = None,
                   tapwater_options: dict | None = None) -> dict:
    """chemicals: [{"cas": ..., "name": ...}]. Returns {"tapwater": ..., "ssl": ..., "notes": [...]}."""
    tap = ResidentTapwaterModule()
    tap_out = tap.run(tap.validate({"chemicals": chemicals, **(tapwater_options or {})}))
    notes, ssl_chems = [], []
    for c, t in zip(chemicals, tap_out.chemicals):
        if t.sl_ug_l is None:
            notes.append(f"{t.name}: no tap water SL; no SSL.")
            continue
        ssl_chems.append({**{k: c[k] for k in ("cas", "name") if c.get(k)}, "target_ug_l": t.sl_ug_l})
    ssl = SslLeachingModule()
    ssl_out = ssl.run(ssl.validate({"chemicals": ssl_chems, **(ssl_options or {})})) if ssl_chems else None
    return {"tapwater": tap_out, "ssl": ssl_out, "notes": notes}


def vleach_to_groundwater(vleach_inputs: dict, daf: float = 1.0, target_mg_l: float | None = None) -> dict:
    """Leachate and groundwater concentration from a VLEACH simulation.

    1. vadose.vleach gives, for every time step and polygon, the mass rate across the
       water table (g/yr; advection plus diffusion).
    2. Leachate concentration C_L = rate / (area * q), g/m3 = mg/L (flux-averaged).
    3. Groundwater concentration C_gw = C_L / DAF: the target leachate concentration is
       target * DAF (RSL UG 4.8.5, PDF p. 194; SSG 2002 Eq. 4-11); DAF = 1 means no dilution.
    Returns the VLEACH output, the C_L and C_gw series per polygon, their peak, and the
    comparison with target_mg_l (e.g. an MCL or a tap water SL / 1000).
    """
    from modules.vadose.vleach import VleachModule

    if daf < 1:
        raise ValueError("DAF must be >= 1.")
    m = VleachModule()
    out = m.run(m.validate(vleach_inputs))
    polygons, notes = [], []
    for p_in, p in zip(vleach_inputs["polygons"], out.polygons):
        q = p_in["recharge_m_yr"]
        if q <= 0:
            notes.append(f"{p.name}: no recharge; leachate concentration undefined.")
            continue
        series = [(t, rate / (p.area_m2 * q)) for t, rate in p.groundwater_series]
        t_max, cl_max = max(series, key=lambda x: x[1])
        res = {"name": p.name, "leachate_mg_l": series, "groundwater_mg_l": [(t, c / daf) for t, c in series],
               "peak_time_yr": t_max, "peak_leachate_mg_l": cl_max, "peak_groundwater_mg_l": cl_max / daf}
        if target_mg_l is not None:
            res["exceeds_target"] = cl_max / daf > target_mg_l
            res["first_exceedance_yr"] = next((t for t, c in series if c / daf > target_mg_l), None)
        polygons.append(res)
    return {"vleach": out, "polygons": polygons, "daf": daf, "target_mg_l": target_mg_l, "notes": notes}
