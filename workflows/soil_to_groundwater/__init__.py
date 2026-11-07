"""Risk-based soil screening levels protective of groundwater.

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
