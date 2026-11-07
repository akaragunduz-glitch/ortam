"""Construction worker soil: screening levels, cancer risk and hazard (subchronic exposure).

Two activities (RSL UG 4.5.1 standard vehicle traffic; 4.5.2 other construction
activities, PDF pp. 116-124) share the same equations and differ only in the
particulate emission factor: PEFsc (traffic, 4.9.2) or PEF'sc (other, 4.9.3).
PEFsc/PEF'sc and VFsc are user inputs here; their calculation is task D1.2.

Noncancer toxicity: subchronic RfD/RfC by default (exposure <= 1 yr; RSL UG
2.3.1.2, 2.3.2.2, 4.5), falling back to chronic values with a warning when the
RSL subchronic table has none; `toxicity_duration: chronic` uses chronic values.
Equations: modules/risk/equations/worker.py.
"""

from __future__ import annotations

from typing import Literal

from pydantic import Field

from modules.risk.equations.worker import Worker  # noqa: F401  (documentation)
from modules.risk.resident_soil import ChemicalInput as _SoilChemical
from ortam_core.chemicals.database import subchronic
from modules.risk.worker_soil import run_soil
from modules.risk.workers_common import SoilResult, WorkerExposure, make_worker, soil_summary
from ortam_core.model import ModelInput, ModelOutput, Module


class ChemicalInput(_SoilChemical):
    vf_sc_m3_kg: float | None = Field(None, gt=0, description="Subchronic volatilization factor VFsc, m3/kg "
                                                              "(required for volatiles)")


class ConstructionExposure(WorkerExposure):
    ew_wk_yr: float | None = Field(None, gt=0, le=52, description="Exposure weeks, wk/yr")
    dw_d_wk: float | None = Field(None, gt=0, le=7, description="Exposure days per week, day/wk")


class ConstructionSoilInput(ModelInput):
    activity: Literal["traffic", "other"] = Field(
        "traffic", description="'traffic': standard vehicle traffic (PEFsc, 4.5.1); 'other': other "
                               "construction activities (PEF'sc, 4.5.2)")
    pef_sc_m3_kg: float = Field(gt=0, description="PEFsc or PEF'sc for the chosen activity, m3/kg")
    toxicity_duration: Literal["subchronic", "chronic"] = "subchronic"
    toxicity_table: Literal["rsl-2024-11", "none"] = "rsl-2024-11"
    target_risk: float = Field(1e-6, gt=0, lt=1, description="Target cancer risk TR, -")
    target_hq: float = Field(1.0, gt=0, description="Target hazard quotient THQ, -")
    exposure: ConstructionExposure = ConstructionExposure()
    chemicals: list[ChemicalInput] = Field(min_length=1)


class ConstructionSoilOutput(ModelOutput):
    activity: str
    toxicity_duration: str
    exposure_overridden: list[str] = []
    chemicals: list[SoilResult]
    total_risk: float | None = None
    hazard_index: float | None = None


class ConstructionSoilModule(Module):
    name = "risk.construction-soil"
    title = "Risk: construction worker soil screening levels, risk and hazard"
    input_type = ConstructionSoilInput
    output_type = ConstructionSoilOutput

    def run(self, inputs: ConstructionSoilInput) -> ConstructionSoilOutput:
        ex = inputs.exposure.model_dump(exclude_none=True)
        w = make_worker("construction", ex)
        warnings = [f"Default exposure parameters overridden: {', '.join(sorted(ex))}."] if ex else []

        def tox(c, overridden):
            if inputs.toxicity_duration == "chronic":
                return c.rfd_mg_kg_day, c.rfc_mg_m3, "chronic"
            srfd, srfc = subchronic(inputs.toxicity_table, c)
            # user-supplied RfD/RfC override everything
            rfd = c.rfd_mg_kg_day if "rfd_mg_kg_day" in overridden else (srfd or c.rfd_mg_kg_day)
            rfc = c.rfc_mg_m3 if "rfc_mg_m3" in overridden else (srfc or c.rfc_mg_m3)
            missing = (c.rfd_mg_kg_day and not srfd) or (c.rfc_mg_m3 and not srfc)
            return rfd, rfc, "chronic (no subchronic value)" if missing else "subchronic"

        results, total, hi = run_soil(inputs, w, inputs.pef_sc_m3_kg, "vf_sc_m3_kg", warnings, tox)
        return ConstructionSoilOutput(activity=inputs.activity, toxicity_duration=inputs.toxicity_duration,
                                      exposure_overridden=sorted(ex), chemicals=results, total_risk=total,
                                      hazard_index=hi, warnings=warnings)

    def summary(self, o: ConstructionSoilOutput) -> str:
        return soil_summary(f"Construction worker soil ({o.activity}; {o.toxicity_duration} noncancer toxicity)", o)
