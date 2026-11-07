"""Shared pieces of the resident risk tools: defaults, inputs, outputs, risk arithmetic.

Sources:
  - Risk = C*TR/SL and HQ = C*THQ/SL: RSL UG 2.6.1 (PDF p. 28).
  - One-hit rule above 0.01: RSL UG 2.6.2 (PDF p. 28); RAGS Part A (1989) p. 8-11.
  - Total risk Risk_T = sum(Risk_i): RAGS Part A p. 8-12; HI = sum(HQ_i): p. 8-13.
  - Total SL = 1 / sum(1/SL_route): RSL UG 4.1.1.1 (PDF p. 37), and each 'Total' equation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache

import yaml
from pydantic import BaseModel, Field, create_model

from modules.risk.rsl_import import DATA
from modules.risk.toxicity import Chemical, lookup
from ortam_core.model import ModelInput

ONE_HIT_THRESHOLD = 0.01     # RSL UG 2.6.2
TCE_CAS = "79-01-6"
VC_CAS = "75-01-4"


# ---------------------------------------------------------------- defaults
@dataclass
class Defaults:
    p: dict[str, float]                    # scalar parameters
    segments: list[dict[str, float]]       # age segments 0-2, 2-6, 6-16, 16-26
    tce: dict[str, float]                  # CAF/MAF fractions
    rba_by_cas: dict[str, float] = field(default_factory=dict)


@lru_cache(maxsize=None)
def _file(version: str) -> dict:
    return yaml.safe_load((DATA / version / "exposure_resident.yaml").read_text(encoding="utf-8"))


def load_defaults(version: str = "2024-11", overrides: dict[str, float] | None = None) -> Defaults:
    f = _file(version)
    p = {k: float(v["value"]) for k, v in f["parameters"].items()}
    p.update(overrides or {})
    segs = [{k: (float(v) if k != "label" else v) for k, v in s.items()} for s in f["age_segments"]]
    tce = {k: float(v["value"]) for k, v in f["tce_fractions"].items()}
    return Defaults(p=p, segments=segs, tce=tce, rba_by_cas={k: float(v) for k, v in f["rba_by_cas"].items()})


def _exposure_model():
    f = _file("2024-11")["parameters"]
    fields = {k: (float | None, Field(None, gt=0, description=f"{v['symbol']}, {v['unit']}"))
              for k, v in f.items() if k not in ("target_risk", "target_hq")}
    return create_model("ExposureOverrides", __base__=ModelInput, **fields)


ExposureOverrides = _exposure_model()


# ------------------------------------------------------------------ inputs
class ChemicalInputBase(ModelInput):
    cas: str | None = Field(None, description="CAS registry number")
    name: str | None = Field(None, description="Chemical name; required if the CAS number is ambiguous")
    # Optional overrides of table values (every override is flagged in the output)
    sfo_per_mg_kg_day: float | None = Field(None, gt=0, description="Oral slope factor, (mg/kg-day)^-1")
    iur_per_ug_m3: float | None = Field(None, gt=0, description="Inhalation unit risk, (ug/m3)^-1")
    rfd_mg_kg_day: float | None = Field(None, gt=0, description="Chronic oral reference dose, mg/kg-day")
    rfc_mg_m3: float | None = Field(None, gt=0, description="Chronic reference concentration, mg/m3")
    giabs: float | None = Field(None, gt=0, le=1, description="Gastrointestinal absorption fraction, -")
    volatile: bool | None = Field(None, description="Volatile chemical (RSL definition)")
    mutagen: bool | None = Field(None, description="Mutagenic mode of action (ADAF applied)")


OVERRIDABLE = ("sfo_per_mg_kg_day", "iur_per_ug_m3", "rfd_mg_kg_day", "rfc_mg_m3", "giabs", "volatile",
               "mutagen", "absd", "mw_g_mol", "log_kow", "kp_cm_h", "b_dimensionless", "tau_event_h", "t_star_h")


def resolve_chemical(table: str, ci: ChemicalInputBase) -> tuple[Chemical, list[str]]:
    """Table lookup plus user overrides. Returns the chemical and the overridden field names."""
    given = {k: getattr(ci, k) for k in OVERRIDABLE if getattr(ci, k, None) is not None}
    if table == "none":
        if not ci.name:
            raise ValueError("With toxicity_table 'none', each chemical needs a name.")
        return Chemical(**{**Chemical(name=ci.name, cas=ci.cas or "").__dict__, **given}), []
    base = lookup(table, ci.cas, ci.name)
    return Chemical(**{**base.__dict__, **given}), sorted(given)


def cancer_method(c: Chemical) -> str:
    """Which cancer equation set applies (RSL UG 4.1.x.3-4.1.x.6)."""
    if c.cas == TCE_CAS:
        return "trichloroethylene"
    if c.cas == VC_CAS:
        return "vinyl chloride"
    return "mutagen" if c.mutagen else "standard"


def label(c: Chemical) -> str:
    return f"{c.name} ({c.cas})" if c.cas else c.name


# ----------------------------------------------------------------- outputs
class Routes(BaseModel):
    ingestion: float | None = None
    dermal: float | None = None
    inhalation: float | None = None
    total: float | None = None


def routes(ing: float | None, derm: float | None, inh: float | None) -> Routes:
    return Routes(ingestion=ing, dermal=derm, inhalation=inh, total=total_sl(ing, derm, inh))


def total_sl(*sls: float | None) -> float | None:
    """Total SL = 1 / sum(1/SL_i) over applicable routes."""
    inv = sum(1.0 / s for s in sls if s)
    return 1.0 / inv if inv else None


def one_hit(linear_risk: float) -> float:
    """RSL UG 2.6.2 and RAGS A (1989) p. 8-11: Risk = 1 - exp(-CDI x SF) above 0.01."""
    return 1.0 - math.exp(-linear_risk) if linear_risk > ONE_HIT_THRESHOLD else linear_risk


def linear(conc: float, target: float, sl: float | None) -> float:
    """RSL UG 2.6.1: Risk = C*TR/SL or HQ = C*THQ/SL."""
    return conc * target / sl if sl else 0.0


def ratio(conc: float, target: float, sl: Routes, risk: bool = False) -> Routes:
    f = one_hit if risk else (lambda x: x)
    return Routes(**{k: (f(linear(conc, target, getattr(sl, k))) if getattr(sl, k) else None)
                     for k in ("ingestion", "dermal", "inhalation", "total")})


def governing(ca: float | None, nc: float | None) -> tuple[float | None, str | None]:
    """Reported SL: the lower of the cancer and noncancer SL (as in the RSL Generic Tables)."""
    c = [(v, b) for v, b in ((ca, "c"), (nc, "n")) if v]
    return min(c) if c else (None, None)
