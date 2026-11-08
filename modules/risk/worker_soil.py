"""Worker soil (composite, outdoor, indoor): screening levels, cancer risk and hazard.

Equations: modules/risk/equations/worker.py (RSL UG 4.2.1, 4.3.1, 4.4.1; PDF pp. 101-114).
Defaults: data/reference/rsl/2024-11/exposure_worker.yaml; PEF from exposure_resident.yaml.
Verification: composite worker = RSL Industrial Soil table (test_worker.py).
"""

from __future__ import annotations

from typing import Literal

from pydantic import Field

from modules.risk.common import (governing, label, linear, load_defaults, one_hit, ratio, resolve_chemical,
                                 routes)
from modules.risk.equations import dose
from modules.risk.equations import worker as eq
from modules.risk.resident_soil import ChemicalInput
from modules.risk.workers_common import SoilResult, WorkerExposure, make_worker, soil_summary
from ortam_core.model import ModelInput, ModelOutput, Module


class WorkerSoilInput(ModelInput):
    receptor: Literal["composite", "outdoor", "indoor"] = Field("composite", description="Worker receptor")
    toxicity_table: Literal["rsl-2024-11", "none"] = "rsl-2024-11"
    target_risk: float = Field(1e-6, gt=0, lt=1, description="Target cancer risk TR, -")
    target_hq: float = Field(1.0, gt=0, description="Target hazard quotient THQ, -")
    pef_m3_kg: float | None = Field(None, gt=0, description="Particulate emission factor, m3/kg (default 1.36E+09)")
    exposure: WorkerExposure = WorkerExposure()
    chemicals: list[ChemicalInput] = Field(min_length=1)


class WorkerSoilOutput(ModelOutput):
    receptor: str
    toxicity_table: str
    exposure_overridden: list[str] = []
    chemicals: list[SoilResult]
    total_risk: float | None = None
    hazard_index: float | None = None


def run_soil(inputs, w, pef, vf_field: str, warnings: list[str], noncancer_tox=None):
    """Shared loop for worker and construction soil. noncancer_tox(c) -> (rfd, rfc, label) or None."""
    d = load_defaults()
    tr, thq = inputs.target_risk, inputs.target_hq
    results = []
    for ci in inputs.chemicals:
        c, overridden = resolve_chemical(inputs.toxicity_table, ci)
        vf = getattr(ci, vf_field)
        if c.volatile and vf is None:
            raise ValueError(f"{label(c)} is volatile; supply {vf_field} (VF calculation is task D1.2).")
        if vf is not None:
            overridden.append(vf_field)
        rba = ci.rba if ci.rba is not None else d.rba_by_cas.get(c.cas, 1.0)
        if ci.rba is not None:
            overridden.append("rba")
        if overridden:
            warnings.append(f"{label(c)}: table values overridden by the user: {', '.join(overridden)}.")
        rfd, rfc, duration = (noncancer_tox(c, overridden) if noncancer_tox
                              else (c.rfd_mg_kg_day, c.rfc_mg_m3, "chronic"))
        if duration.startswith("chronic (no"):
            warnings.append(f"{label(c)}: no subchronic value in the RSL table; chronic RfD/RfC used.")
        nc = routes(*eq.soil_noncancer(rfd, rfc, c.giabs, c.absd, rba, vf, pef, w, thq))
        ca = routes(*eq.soil_cancer(c.sfo_per_mg_kg_day, c.iur_per_ug_m3, c.giabs, c.absd, rba, vf, pef, w, tr))
        sl, basis = governing(ca.total, nc.total)
        if sl is None:
            warnings.append(f"{label(c)}: no toxicity values; no screening level calculated.")
        elif sl > 1.0e5:
            warnings.append(f"{label(c)}: SL {sl:.3g} mg/kg exceeds the 1E+05 mg/kg ceiling limit "
                            "(RSL UG 5.13); flagged, not replaced.")
        r = SoilResult(name=c.name, cas=c.cas, rba=rba, toxicity_duration=duration, sl_nc_mg_kg=nc, sl_ca_mg_kg=ca,
                       sl_mg_kg=sl, sl_basis=basis, overridden=overridden, toxicity_keys=dict(c.keys))
        if ci.concentration_mg_kg is not None:
            C = r.concentration_mg_kg = ci.concentration_mg_kg
            r.risk, r.hq = ratio(C, tr, ca, risk=True), ratio(C, thq, nc)
            r.doses = dose.worker_soil(C, c.absd, rba, vf, pef, w)
        results.append(r)
    with_c = [r for r in results if r.concentration_mg_kg is not None]
    total = hi = None
    if with_c:
        total = one_hit(sum(linear(r.concentration_mg_kg, tr, r.sl_ca_mg_kg.total) for r in with_c))
        hi = sum(r.hq.total or 0.0 for r in with_c)
    return results, total, hi


class WorkerSoilModule(Module):
    name = "risk.worker-soil"
    title = "Risk: worker soil (composite, outdoor, indoor) screening levels, risk and hazard"
    input_type = WorkerSoilInput
    output_type = WorkerSoilOutput

    def run(self, inputs: WorkerSoilInput) -> WorkerSoilOutput:
        ex = inputs.exposure.model_dump(exclude_none=True)
        w = make_worker(inputs.receptor, ex)
        pef = inputs.pef_m3_kg or load_defaults().p["pef_m3_kg"]
        warnings = [f"Default exposure parameters overridden: {', '.join(sorted(ex))}."] if ex else []
        if inputs.pef_m3_kg:
            warnings.append("PEF supplied by the user.")
        results, total, hi = run_soil(inputs, w, pef, "vf_m3_kg", warnings)
        return WorkerSoilOutput(receptor=inputs.receptor, toxicity_table=inputs.toxicity_table,
                                exposure_overridden=sorted(ex), chemicals=results, total_risk=total,
                                hazard_index=hi, warnings=warnings)

    def summary(self, o: WorkerSoilOutput) -> str:
        return soil_summary(f"Worker soil ({o.receptor}; toxicity: {o.toxicity_table})", o)
