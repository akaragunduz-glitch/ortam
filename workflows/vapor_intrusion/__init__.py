"""Vapor intrusion workflows: risk air tools <-> vapor.visl, vapor.johnson-ettinger.

  1. The target indoor air concentration is the risk-based air SL from
     risk.resident-air (residential) or risk.worker-air (commercial/industrial
     worker), as in the VISL calculator (VISL UG 2014, pp. 5-7; RAGS F).
  2. vapor.visl converts it into sub-slab/soil gas and groundwater targets.
Modules never call each other; this workflow connects them.
"""

from __future__ import annotations

from modules.risk.resident_air import ResidentAirModule
from modules.risk.worker_air import WorkerAirModule
from modules.vapor.visl import ViSlModule


def visl_screening(chemicals: list[dict], land_use: str = "residential", visl_options: dict | None = None,
                   air_options: dict | None = None) -> dict:
    """chemicals: [{"cas": ..., "name": ...}]; land_use: 'residential' or 'commercial'."""
    if land_use == "residential":
        air, extra = ResidentAirModule(), {}
    elif land_use == "commercial":
        air, extra = WorkerAirModule(), {"receptor": "indoor"}
    else:
        raise ValueError("land_use must be 'residential' or 'commercial'.")
    ids = [{k: c[k] for k in ("cas", "name") if c.get(k)} for c in chemicals]
    air_out = air.run(air.validate({"chemicals": ids, **extra, **(air_options or {})}))
    notes, v_chems = [], []
    for c, a in zip(ids, air_out.chemicals):
        if a.sl_ug_m3 is None:
            notes.append(f"{a.name}: no inhalation toxicity value; not evaluated.")
            continue
        v_chems.append({**c, "target_ia_ug_m3": a.sl_ug_m3})
    v = ViSlModule()
    v_out = v.run(v.validate({"chemicals": v_chems, **(visl_options or {})})) if v_chems else None
    return {"air": air_out, "visl": v_out, "notes": notes}


def visl_risk_from_subsurface(chemicals: list[dict], land_use: str = "residential", af_subslab: float = 0.03,
                              af_groundwater: float = 0.001, air_options: dict | None = None) -> dict:
    """Indoor air estimated from measured sub-slab/soil gas or groundwater, then risk.

    chemicals: [{"cas": ..., "subslab_ug_m3": ...} or {"cas": ..., "groundwater_ug_l": ...}]
    Indoor air: C_ia = C_sg * AF_ss, or C_ia = C_gw * H' * AF_gw * 1000 (inverse of VISL UG pp. 8-9);
    risk and HQ from risk.resident-air / risk.worker-air (RSL UG 2.6).
    """
    from ortam_core.chemicals.database import lookup

    rows, notes = [], []
    for c in chemicals:
        if "subslab_ug_m3" in c:
            cia = c["subslab_ug_m3"] * af_subslab
        else:
            h = lookup("rsl-2024-11", c.get("cas"), c.get("name")).henry_dimensionless
            cia = c["groundwater_ug_l"] * h * af_groundwater * 1000.0
        rows.append({**{k: c[k] for k in ("cas", "name") if c.get(k)}, "concentration_ug_m3": cia})
    if land_use == "residential":
        air, extra = ResidentAirModule(), {}
    else:
        air, extra = WorkerAirModule(), {"receptor": "indoor"}
    out = air.run(air.validate({"chemicals": rows, **extra, **(air_options or {})}))
    return {"air": out, "notes": notes}


def je_risk(je_inputs: dict, land_use: str = "residential", air_options: dict | None = None) -> dict:
    """Site-specific vapor intrusion: vapor.johnson-ettinger -> risk air tool.

    1. vapor.johnson-ettinger gives the indoor air concentration Cia (ug/m3) from a
       groundwater or soil gas concentration (J&E v6, 2017).
    2. risk.resident-air (residential) or risk.worker-air with receptor 'indoor'
       (commercial) gives the cancer risk and HQ for Cia (RSL UG 2.6), with the
       current RSL toxicity values. The J&E spreadsheet uses the same equations
       (MODEL!G153-G154) with its 2017 toxicity table.
    A chemical whose result is invalid (e.g. above solubility) gets no risk.
    """
    from modules.vapor.johnson_ettinger import JohnsonEttingerModule

    j = JohnsonEttingerModule()
    je_out = j.run(j.validate(je_inputs))
    rows, notes = [], []
    for c in je_out.chemicals:
        if c.cia_ug_m3 is None:
            notes.append(f"{c.name}: J&E result invalid; risk not computed.")
            continue
        rows.append({"cas": c.cas, "concentration_ug_m3": c.cia_ug_m3})
    if land_use == "residential":
        air, extra = ResidentAirModule(), {}
    elif land_use == "commercial":
        air, extra = WorkerAirModule(), {"receptor": "indoor"}
    else:
        raise ValueError("land_use must be 'residential' or 'commercial'.")
    air_out = air.run(air.validate({"chemicals": rows, **extra, **(air_options or {})})) if rows else None
    return {"je": je_out, "air": air_out, "notes": notes}
