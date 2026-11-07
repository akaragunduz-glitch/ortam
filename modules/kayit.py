"""Modül kaydı: komut satırı ve arayüzün çalıştırabileceği tüm modüller.

Yeni bir modül eklendiğinde yalnızca bu sözlüğe bir satır eklenir.
"""

from importlib import import_module

MODULLER: dict[str, str] = {
    "groundwater.uc-nokta": "modules.groundwater.uc_nokta:UcNoktaModulu",
    # "groundwater.modflow": "modules.groundwater.modflow:ModflowModulu",   # Faz 5
    # "vadose.ssl":          "modules.vadose.ssl:SSLSiziniModulu",          # Faz 1
    # "vapor.je":            "modules.vapor.johnson_ettinger:JEModulu",     # Faz 1
    # "air.aermod":          "modules.air.aermod:AermodModulu",             # Faz 5
    # "geochem.phreeqc":     "modules.geochem.phreeqc:PhreeqcModulu",       # Faz 5
    # "risk.rags":           "modules.risk.rags:RagsModulu",                # Faz 1
}


def modul_al(ad: str):
    if ad not in MODULLER:
        raise KeyError(f"Bilinmeyen modül: {ad}. Mevcut: {', '.join(MODULLER)}")
    yol, sinif = MODULLER[ad].split(":")
    return getattr(import_module(yol), sinif)()
