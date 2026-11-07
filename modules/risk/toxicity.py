"""Toxicity and chemical data access (RSL tables imported by rsl_import.py).

Lookup is by CAS number; when one CAS number appears in several rows
(e.g. Cadmium (Diet) and Cadmium (Water)), the chemical name is required.
Values are returned exactly as in the RSL table; empty cells become None.

Fortran analogy: a read-only COMMON block filled once from a data file.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache

from modules.risk.rsl_import import DATA

TABLES = {"rsl-2024-11": "2024-11"}


@dataclass(frozen=True)
class Chemical:
    name: str
    cas: str
    sfo_per_mg_kg_day: float | None = None   # oral slope factor, (mg/kg-day)^-1
    iur_per_ug_m3: float | None = None       # inhalation unit risk, (ug/m3)^-1
    rfd_mg_kg_day: float | None = None       # chronic oral reference dose, mg/kg-day
    rfc_mg_m3: float | None = None           # chronic inhalation reference concentration, mg/m3
    giabs: float = 1.0                       # gastrointestinal absorption fraction, -
    absd: float | None = None                # dermal absorption fraction from soil, -
    volatile: bool = False
    mutagen: bool = False
    csat_mg_kg: float | None = None          # soil saturation concentration, mg/kg
    keys: tuple[tuple[str, str], ...] = ()   # source keys of toxicity values (I = IRIS, P = PPRTV, ...)


def _f(text: str) -> float | None:
    return float(text) if text not in ("", None) else None


@lru_cache(maxsize=None)
def _rows(table: str) -> tuple[dict, ...]:
    if table not in TABLES:
        raise ValueError(f"Unknown toxicity table '{table}'. Available: {', '.join(TABLES)}")
    with (DATA / TABLES[table] / "toxicity.csv").open(encoding="utf-8") as f:
        return tuple(csv.DictReader(f))


def lookup(table: str, cas: str | None = None, name: str | None = None) -> Chemical:
    rows = _rows(table)
    if cas:
        hits = [r for r in rows if r["cas"] == cas.strip()]
        if name:
            hits = [r for r in hits if r["name"].lower() == name.strip().lower()]
    elif name:
        hits = [r for r in rows if r["name"].lower() == name.strip().lower()]
    else:
        raise ValueError("Each chemical needs a CAS number or a name.")
    if not hits:
        raise ValueError(f"Chemical not found in {table}: cas={cas!r}, name={name!r}")
    if len(hits) > 1:
        names = "; ".join(r["name"] for r in hits)
        raise ValueError(f"CAS {cas} matches several entries ({names}); give the name as well.")
    r = hits[0]
    return Chemical(
        name=r["name"], cas=r["cas"],
        sfo_per_mg_kg_day=_f(r["sfo_per_mg_kg_day"]), iur_per_ug_m3=_f(r["iur_per_ug_m3"]),
        rfd_mg_kg_day=_f(r["rfd_mg_kg_day"]), rfc_mg_m3=_f(r["rfc_mg_m3"]),
        giabs=_f(r["giabs"]) or 1.0, absd=_f(r["absd"]),
        volatile=r["volatile"] == "V", mutagen=r["mutagen"] == "M",
        csat_mg_kg=_f(r["csat_mg_kg"]),
        keys=tuple((k, r[f"{k}_key"]) for k in ("sfo", "iur", "rfd", "rfc") if r[f"{k}_key"]),
    )
