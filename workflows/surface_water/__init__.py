"""Surface water workflows: groundwater plume or watershed -> water body -> fish and risk.

  plume_to_water_body:
    1. napl.source-plume (BIOSCREEN 1.4) with the modelled length set to the distance of the
       gaining stream or lake shore; the mass flux across that cross-section (mg/day; numerical
       integration or BIOSCREEN's 5-point estimate, BSCR14 p. 2) of the chosen plume model is
       converted to g/yr (365 d/yr, as BIOSCREEN);
    2. surfacewater.watershed with that flux as an external load (ORTAM extension of Table B-4-7).
  water_body_fish_risk:
    surfacewater.watershed -> Cdw, Csb -> workflows.foodchain (foodchain.media fish, Tables B-4-26
    to B-4-28; risk.foodchain-ingestion). For deposition-derived soil concentrations both bases are
    evaluated: Cs for cancer risk, CstD for the hazard index (Table B-4-1, p. B-168 [171]).
Modules never call each other; this workflow connects them.
"""
from __future__ import annotations

import copy

from modules.napl.source_plume import SourcePlumeModule
from modules.surfacewater.watershed import WatershedModule
from workflows.foodchain import foodchain_risk

MODELS = {"no_decay": "No degradation", "first_order": "First-order decay", "instantaneous": "Instantaneous reaction"}
MERCURY_SPECIES_CAS = {"Hg2+": "7487-94-7", "MHg": "22967-92-6"}


def plume_to_water_body(plume_inputs: dict, watershed_inputs: dict, discharge_distance_m: float,
                        model: str = "first_order", flux: str = "integrated") -> dict:
    """Groundwater plume discharging to a water body (gaining stream or lake)."""
    if model not in MODELS:
        raise ValueError(f"model must be one of {', '.join(MODELS)}.")
    if flux not in ("integrated", "bioscreen"):
        raise ValueError("flux must be 'integrated' or 'bioscreen'.")
    p_in = copy.deepcopy(plume_inputs)
    p_in["general"]["model_length_m"] = discharge_distance_m
    pm = SourcePlumeModule()
    plume = pm.run(pm.validate(p_in))
    notes = [f"Plume model length set to the discharge distance ({discharge_distance_m:g} m)."]
    if flux == "integrated":
        if plume.integrated is None:
            raise ValueError("Integrated mass flux is not available with vertical dispersion (alpha_z > 0).")
        f_mg_day = next(i for i in plume.integrated if i.model == MODELS[model]).mass_flux_mg_day[-1]
    else:
        mb = next(m for m in plume.mass_balance if m.model == MODELS[model])
        if mb.mass_flux_mg_day is None:
            raise ValueError("BIOSCREEN mass flux is not available with vertical dispersion (alpha_z > 0).")
        f_mg_day = mb.mass_flux_mg_day[-1]
        notes.append("BIOSCREEN 5-point mass flux: about +-50 % (BSCR14 p. 2).")
    load_g_yr = f_mg_day * 365.0 / 1000.0
    w_in = copy.deepcopy(watershed_inputs)
    w_in.setdefault("external_loads", []).append({"name": f"groundwater plume ({MODELS[model]})", "load_g_yr": load_g_yr})
    wm = WatershedModule()
    water = wm.run(wm.validate(w_in))
    return {"plume": plume, "mass_flux_mg_day": f_mg_day, "load_g_yr": load_g_yr, "water_body": water, "notes": notes}


def water_body_fish_risk(watershed_inputs: dict, receptor: str = "fisher", media_options: dict | None = None,
                         ingestion_options: dict | None = None) -> dict:
    """Water body concentrations -> fish concentration -> fish ingestion risk (HHRAP)."""
    wm = WatershedModule()
    water = wm.run(wm.validate(watershed_inputs))
    chem = watershed_inputs["chemical"]
    out = {"water_body": water, "risk": {}}
    for r in water.results:
        if water.mercury:
            chems = [{"cas": cas, "surface_water_dissolved_mg_l": r.dissolved_species_mg_l[s]}
                     for s, cas in MERCURY_SPECIES_CAS.items()]
        else:
            c = {k: chem[k] for k in ("cas", "name") if chem.get(k)}
            chems = [{**c, "surface_water_dissolved_mg_l": r.dissolved_mg_l, "sediment_mg_kg": r.bed_sediment_mg_kg}]
        out["risk"][r.basis] = foodchain_risk(chems, receptor=receptor, media_options=media_options,
                                              ingestion_options=ingestion_options)
    bases = out["risk"]
    if "measured" in bases:
        out["total_risk"] = bases["measured"]["total_risk"]
        out["hazard_index"] = bases["measured"]["hazard_index"]
    else:   # Table B-4-1: Cs for carcinogens, CstD for noncarcinogens
        out["total_risk"] = bases["carcinogens"]["total_risk"]
        out["hazard_index"] = bases["noncarcinogens"]["hazard_index"]
    return out
