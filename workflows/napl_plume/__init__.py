"""LNAPL source -> dissolved plume -> receptor workflows.

  napl_plume_risk:
    1. (optional) napl.equilibrium, basis "napl": effective solubility x_k gamma_k S_k of one
       NAPL component (Raoult's law; BSCR p. 3 of the Version 1.4 revisions, HSSM1 Eq. 43-44).
       It replaces the concentration of the central source zone: the dissolved concentration in
       equilibrium with the LNAPL, an upper bound of the source-zone concentration.
    2. napl.source-plume (BIOSCREEN 1.4): concentrations at the receptors, for one of the three
       plume models.
    3. risk.resident-tapwater: receptor concentration (mg/L -> ug/L) as the exposure point
       concentration in tap water; cancer risk and hazard (RSL tap water equations).
Modules never call each other; this workflow connects them.
"""
from __future__ import annotations

import copy

from modules.napl.equilibrium import NaplEquilibriumModule
from modules.napl.source_plume import SourcePlumeModule
from modules.risk.resident_tapwater import ResidentTapwaterModule

MODELS = {"no_decay": "no_decay_mg_l", "first_order": "first_order_mg_l", "instantaneous": "instantaneous_mg_l"}


def _match(component, key: str) -> bool:
    return key in {component.cas, component.name} or (component.name or "").lower() == key.lower()


def napl_plume_risk(plume_inputs: dict, model: str = "first_order", napl_inputs: dict | None = None,
                    component: str | None = None, tapwater_options: dict | None = None) -> dict:
    """plume_inputs: napl.source-plume input with receptors and a chemical (CAS) for the risk step.
    model: plume model used for the exposure concentration ('no_decay', 'first_order', 'instantaneous').
    napl_inputs/component: optional napl.equilibrium input (basis 'napl') and the CAS or name of the
    component whose effective solubility becomes the central source-zone concentration."""
    if model not in MODELS:
        raise ValueError(f"model must be one of {', '.join(MODELS)}.")
    plume_inputs = copy.deepcopy(plume_inputs)
    notes, eq_out = [], None
    if napl_inputs is not None:
        if component is None:
            raise ValueError("component (CAS or name) is required with napl_inputs.")
        if napl_inputs.get("basis") != "napl":
            raise ValueError("napl_inputs must use basis 'napl' (known NAPL composition).")
        eqm = NaplEquilibriumModule()
        eq_out = eqm.run(eqm.validate(napl_inputs))
        comp = next((c for c in eq_out.components if _match(c, component)), None)
        if comp is None or comp.effective_solubility_mg_l is None:
            raise ValueError(f"No effective solubility for component {component}.")
        plume_inputs["source"]["zones"][0]["concentration_mg_l"] = comp.effective_solubility_mg_l
        notes.append(f"Central source zone concentration = effective solubility of {comp.name} "
                     f"({comp.effective_solubility_mg_l:.4g} mg/L) from napl.equilibrium.")
    pm = SourcePlumeModule()
    plume = pm.run(pm.validate(plume_inputs))
    receptors = []
    chem = plume_inputs.get("chemical") or {}
    tw_out = None
    if plume.receptors:
        if not chem.get("cas") and not chem.get("name"):
            raise ValueError("plume_inputs needs a chemical (CAS) for the tap water risk step.")
        tw = ResidentTapwaterModule()
        ident = {k: chem[k] for k in ("cas", "name") if chem.get(k)}
        tw_out = []
        for r in plume.receptors:      # one tap water run per receptor (separate exposure points)
            c_mg_l = getattr(r, MODELS[model])
            o = tw.run(tw.validate({"chemicals": [{**ident, "concentration_ug_l": c_mg_l * 1000.0}],
                                    **(tapwater_options or {})}))
            c = o.chemicals[0]
            tw_out.append(o)
            receptors.append({"name": r.name, "x_m": r.x_m, "y_m": r.y_m, "concentration_mg_l": c_mg_l,
                              "sl_ug_l": c.sl_ug_l, "exceeds_sl": c_mg_l * 1000.0 > c.sl_ug_l if c.sl_ug_l else None,
                              "risk": o.total_risk, "hazard_index_child": o.hazard_index_child})
    return {"equilibrium": eq_out, "plume": plume, "tapwater": tw_out, "model": model, "receptors": receptors,
            "notes": notes}
