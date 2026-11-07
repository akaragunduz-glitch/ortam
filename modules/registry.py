"""Module registry: every module the CLI and the UI can run.

Adding a module requires one line in this dictionary.
"""

from importlib import import_module

MODULES: dict[str, str] = {
    "groundwater.three-point": "modules.groundwater.three_point:ThreePointModule",
    "risk.resident-soil":      "modules.risk.resident_soil:ResidentSoilModule",
    "risk.resident-tapwater":  "modules.risk.resident_tapwater:ResidentTapwaterModule",
    "risk.resident-air":       "modules.risk.resident_air:ResidentAirModule",
    "risk.worker-soil":        "modules.risk.worker_soil:WorkerSoilModule",
    "risk.worker-air":         "modules.risk.worker_air:WorkerAirModule",
    "risk.construction-soil":  "modules.risk.construction_soil:ConstructionSoilModule",
    "risk.recreator-soil":     "modules.risk.recreator_soil:RecreatorSoilModule",
    "risk.recreator-surface-water": "modules.risk.recreator_surface_water:RecreatorSurfaceWaterModule",
    "risk.fish":               "modules.risk.fish:FishModule",
    "air.ssl-vf":              "modules.air.ssl_vf:SslVfModule",
    "air.ssl-pef":             "modules.air.ssl_pef:SslPefModule",
    # "groundwater.modflow":  "modules.groundwater.modflow:ModflowModule",       # Phase 5
    "vadose.ssl-leaching":     "modules.vadose.ssl_leaching:SslLeachingModule",
    # "vapor.johnson-ettinger": "modules.vapor.johnson_ettinger:JEModule",       # Phase 1
    # "air.aermod":           "modules.air.aermod:AermodModule",                 # Phase 5
    # "geochem.phreeqc":      "modules.geochem.phreeqc:PhreeqcModule",           # Phase 5
}


def get_module(name: str):
    if name not in MODULES:
        raise KeyError(f"Unknown module: {name}. Available: {', '.join(MODULES)}")
    path, cls = MODULES[name].split(":")
    return getattr(import_module(path), cls)()
