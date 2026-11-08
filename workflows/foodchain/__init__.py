"""Food chain workflow: soil / water / sediment -> foodchain.media -> risk.foodchain-ingestion.

Modules never call each other; this workflow connects them:
  1. foodchain.media gives produce, animal product and fish concentrations (HHRAP 2005 App. B);
  2. risk.foodchain-ingestion gives intake, CDI, risk and HQ for the chosen receptor (App. C);
  3. optionally, direct soil contact (ingestion, dermal, inhalation) from risk.resident-soil
     is added for the total (the HHRAP farmer also has direct soil exposure, Table C-1-1;
     ORTAM uses the RSL resident soil equations for that part).
"""

from __future__ import annotations

from modules.foodchain.media import FoodchainMediaModule
from modules.risk.foodchain_ingestion import FIELD, FoodchainIngestionModule


def foodchain_risk(chemicals: list[dict], receptor: str = "farmer", media_options: dict | None = None,
                   ingestion_options: dict | None = None, direct_soil: bool = False) -> dict:
    """chemicals: [{"cas": ..., "soil_mg_kg": ..., "surface_water_dissolved_mg_l": ..., ...}]"""
    m = FoodchainMediaModule()
    media = m.run(m.validate({"chemicals": chemicals, **(media_options or {})}))
    rows = []
    for c_in, c in zip(chemicals, media.chemicals):
        foods = {f: getattr(c, f) for f in FIELD.values() if getattr(c, f) is not None}
        rows.append({**{k: c_in[k] for k in ("cas", "name") if c_in.get(k)}, "foods": foods})
    r = FoodchainIngestionModule()
    ingestion = r.run(r.validate({"receptor": receptor, "chemicals": rows, **(ingestion_options or {})}))
    out = {"media": media, "ingestion": ingestion, "soil": None}
    total_risk = ingestion.total_risk or 0.0
    total_hi = ingestion.hazard_index or 0.0
    if direct_soil:
        from modules.risk.resident_soil import ResidentSoilModule
        s = ResidentSoilModule()
        soil_chems = [{**{k: c[k] for k in ("cas", "name", "vf_m3_kg") if c.get(k)},
                       "concentration_mg_kg": c["soil_mg_kg"]} for c in chemicals if c.get("soil_mg_kg") is not None]
        soil = s.run(s.validate({"chemicals": soil_chems}))
        out["soil"] = soil
        total_risk += soil.total_risk or 0.0
        child = receptor.endswith("child")
        total_hi += (soil.hazard_index_child if child else soil.hazard_index_adult) or 0.0
    out["total_risk"] = total_risk
    out["hazard_index"] = total_hi
    return out
