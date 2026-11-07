"""Soil screening with volatilization: air.ssl-vf -> risk soil tools.

The modules never call each other; this workflow connects them:
  1. air.ssl-vf computes VF (and Csat) for each chemical;
  2. the VF is passed to risk.resident-soil / risk.worker-soil as vf_m3_kg.
Where VF cannot be calculated (missing H', Dia, Diw or Kd), inhalation is
evaluated with the PEF only, as in the RSL tables; this is reported.
SLs above Csat are flagged, not replaced (RSL UG 5.12).
"""

from __future__ import annotations

from modules.air.ssl_vf import SslVfModule
from modules.risk.resident_soil import ResidentSoilModule
from modules.risk.worker_soil import WorkerSoilModule


def soil_screening(chemicals: list[dict], receptor: str = "resident", vf_options: dict | None = None,
                   risk_options: dict | None = None) -> dict:
    """Run VF + soil SLs. receptor: 'resident' or 'composite' / 'outdoor' / 'indoor' (workers).

    chemicals: [{"cas": ..., "name": ..., "concentration_mg_kg": ...}, ...]
    Returns {"vf": SslVfOutput, "risk": risk output, "notes": [...]}.
    """
    vf_tool = SslVfModule()
    ids = [{k: c[k] for k in ("cas", "name") if c.get(k)} for c in chemicals]
    vf_out = vf_tool.run(vf_tool.validate({"chemicals": ids, **(vf_options or {})}))

    notes, risk_chems = [], []
    for c, v in zip(chemicals, vf_out.chemicals):
        rc = dict(c)
        if v.volatile and v.vf_m3_kg is not None:
            rc["vf_m3_kg"] = v.vf_m3_kg
        elif v.volatile:
            rc["volatile"] = False
            notes.append(f"{v.name}: VF not calculable; inhalation evaluated with PEF only (RSL practice).")
        risk_chems.append(rc)

    if receptor == "resident":
        tool, extra = ResidentSoilModule(), {}
    else:
        tool, extra = WorkerSoilModule(), {"receptor": receptor}
    risk_out = tool.run(tool.validate({"chemicals": risk_chems, **extra, **(risk_options or {})}))

    for v, r in zip(vf_out.chemicals, risk_out.chemicals):
        if v.csat_mg_kg and r.sl_mg_kg and r.sl_mg_kg > v.csat_mg_kg:
            notes.append(f"{v.name}: SL {r.sl_mg_kg:.3g} mg/kg exceeds Csat {v.csat_mg_kg:.3g} mg/kg "
                         f"({'liquid' if v.liquid else 'solid'}); flagged, not replaced (RSL UG 5.12).")
    return {"vf": vf_out, "risk": risk_out, "notes": notes}
