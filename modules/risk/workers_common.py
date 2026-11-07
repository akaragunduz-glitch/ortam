"""Shared input/output pieces of the worker tools (worker soil/air, construction)."""

from __future__ import annotations

from functools import lru_cache
from typing import Literal

import yaml
from pydantic import BaseModel, Field

from modules.risk.common import Routes
from modules.risk.equations.worker import Worker
from ortam_core.chemicals.rsl_import import DATA
from ortam_core.model import ModelInput


@lru_cache(maxsize=None)
def _file(version: str = "2024-11") -> dict:
    return yaml.safe_load((DATA / version / "exposure_worker.yaml").read_text(encoding="utf-8"))


def worker_defaults(receptor: str) -> dict[str, float]:
    return {k: float(v["value"]) for k, v in _file()["receptors"][receptor].items()}


class WorkerExposure(ModelInput):
    """Overrides of the receptor defaults in exposure_worker.yaml (flagged in the output)."""
    ef_d_yr: float | None = Field(None, gt=0, le=365, description="Exposure frequency, day/yr")
    ed_yr: float | None = Field(None, gt=0, description="Exposure duration, yr")
    et_h_d: float | None = Field(None, gt=0, le=24, description="Exposure time, h/day")
    bw_kg: float | None = Field(None, gt=0, description="Body weight, kg")
    irs_mg_d: float | None = Field(None, gt=0, description="Soil ingestion rate, mg/day")
    sa_cm2_d: float | None = Field(None, gt=0, description="Skin surface area, cm2/day")
    af_mg_cm2: float | None = Field(None, gt=0, description="Soil adherence factor, mg/cm2")
    lt_yr: float | None = Field(None, gt=0, description="Lifetime, yr")


def make_worker(receptor: str, overrides: dict[str, float]) -> Worker:
    d = {**worker_defaults(receptor), **overrides}
    if receptor == "construction":
        ef = d["ew_wk_yr"] * d["dw_d_wk"]           # RSL UG 4.5.1.1, p. 117: EF = EW*DW
        at_nc = d["ew_wk_yr"] * 7.0 * d["ed_yr"]     # AT_con-a = EW*7*ED
    else:
        ef = d["ef_d_yr"]
        at_nc = 365.0 * d["ed_yr"]
    return Worker(ef_d_yr=ef, ed_yr=d["ed_yr"], et_h_d=d["et_h_d"], bw_kg=d["bw_kg"], irs_mg_d=d.get("irs_mg_d", 0.0),
                  sa_cm2_d=d.get("sa_cm2_d"), af_mg_cm2=d.get("af_mg_cm2"), at_nc_d=at_nc,
                  lt_yr=d.get("lt_yr", 70.0))


class SoilResult(BaseModel):
    name: str
    cas: str
    rba: float
    toxicity_duration: str = "chronic"
    sl_nc_mg_kg: Routes
    sl_ca_mg_kg: Routes
    sl_mg_kg: float | None = Field(description="Screening level: lower of cancer and noncancer SL, mg/kg")
    sl_basis: Literal["c", "n"] | None = None
    concentration_mg_kg: float | None = None
    risk: Routes | None = None
    hq: Routes | None = None
    overridden: list[str] = []
    toxicity_keys: dict[str, str] = {}


def soil_summary(title: str, o) -> str:
    lines = [title, f"{'Chemical':32s} {'CAS':12s} {'SL mg/kg':>10s} {'basis':>5s} {'C mg/kg':>10s} "
                    f"{'Risk':>9s} {'HQ':>9s}"]
    for c in o.chemicals:
        sl = f"{c.sl_mg_kg:.2g}" if c.sl_mg_kg else "-"
        conc = f"{c.concentration_mg_kg:.3g}" if c.concentration_mg_kg is not None else "-"
        risk = f"{c.risk.total:.1e}" if c.risk and c.risk.total is not None else "-"
        hq = f"{c.hq.total:.2g}" if c.hq and c.hq.total is not None else "-"
        lines.append(f"{c.name[:32]:32s} {c.cas:12s} {sl:>10s} {c.sl_basis or '-':>5s} {conc:>10s} "
                     f"{risk:>9s} {hq:>9s}{' *' if c.overridden else ''}")
    if o.total_risk is not None:
        lines.append(f"Total cancer risk : {o.total_risk:.1e}")
        lines.append(f"Hazard index      : {o.hazard_index:.2g}")
    lines += [f"WARNING: {w}" for w in o.warnings]
    return "\n".join(lines)
