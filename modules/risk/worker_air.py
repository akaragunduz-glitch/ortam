"""Worker air (composite, outdoor, indoor): screening levels, cancer risk and hazard.

Equations: modules/risk/equations/worker.py (RSL UG 4.2.2, 4.3.2, 4.4.2; PDF pp. 105-116).
Verification: composite worker = RSL Industrial Air table (test_worker.py).
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from modules.risk.common import ChemicalInputBase, governing, label, one_hit, resolve_chemical
from modules.risk.equations import worker as eq
from modules.risk.workers_common import WorkerExposure, make_worker
from ortam_core.model import ModelInput, ModelOutput, Module


class ChemicalInput(ChemicalInputBase):
    concentration_ug_m3: float | None = Field(None, ge=0, description="Air concentration, ug/m3")


class WorkerAirInput(ModelInput):
    receptor: Literal["composite", "outdoor", "indoor"] = Field("composite", description="Worker receptor")
    toxicity_table: Literal["rsl-2024-11", "none"] = "rsl-2024-11"
    target_risk: float = Field(1e-6, gt=0, lt=1, description="Target cancer risk TR, -")
    target_hq: float = Field(1.0, gt=0, description="Target hazard quotient THQ, -")
    exposure: WorkerExposure = WorkerExposure()
    chemicals: list[ChemicalInput] = Field(min_length=1)


class AirResult(BaseModel):
    name: str
    cas: str
    sl_nc_ug_m3: float | None = None
    sl_ca_ug_m3: float | None = None
    sl_ug_m3: float | None = None
    sl_basis: Literal["c", "n"] | None = None
    concentration_ug_m3: float | None = None
    risk: float | None = None
    hq: float | None = None
    overridden: list[str] = []


class WorkerAirOutput(ModelOutput):
    receptor: str
    chemicals: list[AirResult]
    total_risk: float | None = None
    hazard_index: float | None = None


class WorkerAirModule(Module):
    name = "risk.worker-air"
    title = "Risk: worker air (composite, outdoor, indoor) screening levels, risk and hazard"
    input_type = WorkerAirInput
    output_type = WorkerAirOutput

    def run(self, inputs: WorkerAirInput) -> WorkerAirOutput:
        ex = inputs.exposure.model_dump(exclude_none=True)
        w = make_worker(inputs.receptor, ex)
        tr, thq = inputs.target_risk, inputs.target_hq
        warnings = [f"Default exposure parameters overridden: {', '.join(sorted(ex))}."] if ex else []
        results = []
        for ci in inputs.chemicals:
            c, overridden = resolve_chemical(inputs.toxicity_table, ci)
            if overridden:
                warnings.append(f"{label(c)}: table values overridden by the user: {', '.join(overridden)}.")
            nc, ca = eq.air_noncancer(c.rfc_mg_m3, w, thq), eq.air_cancer(c.iur_per_ug_m3, w, tr)
            sl, basis = governing(ca, nc)
            r = AirResult(name=c.name, cas=c.cas, sl_nc_ug_m3=nc, sl_ca_ug_m3=ca, sl_ug_m3=sl, sl_basis=basis,
                          overridden=overridden)
            if ci.concentration_ug_m3 is not None:
                C = r.concentration_ug_m3 = ci.concentration_ug_m3
                r.risk = one_hit(C * tr / ca) if ca else None
                r.hq = C * thq / nc if nc else None
            results.append(r)
        with_c = [r for r in results if r.concentration_ug_m3 is not None]
        total = hi = None
        if with_c:
            total = one_hit(sum(r.concentration_ug_m3 * tr / r.sl_ca_ug_m3 for r in with_c if r.sl_ca_ug_m3))
            hi = sum(r.hq or 0.0 for r in with_c)
        return WorkerAirOutput(receptor=inputs.receptor, chemicals=results, total_risk=total, hazard_index=hi,
                               warnings=warnings)

    def summary(self, o: WorkerAirOutput) -> str:
        lines = [f"Worker air ({o.receptor})", f"{'Chemical':32s} {'SL ug/m3':>10s} {'basis':>5s} {'Risk':>9s} {'HQ':>9s}"]
        for c in o.chemicals:
            lines.append(f"{c.name[:32]:32s} {(f'{c.sl_ug_m3:.2g}' if c.sl_ug_m3 else '-'):>10s} "
                         f"{c.sl_basis or '-':>5s} {(f'{c.risk:.1e}' if c.risk is not None else '-'):>9s} "
                         f"{(f'{c.hq:.2g}' if c.hq is not None else '-'):>9s}")
        if o.total_risk is not None:
            lines += [f"Total cancer risk : {o.total_risk:.1e}", f"Hazard index      : {o.hazard_index:.2g}"]
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
