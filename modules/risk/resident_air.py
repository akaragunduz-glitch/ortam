"""Resident air: risk-based screening levels, cancer risk and hazard quotients.

Exposure route: inhalation of ambient air. No assumption is made about how
contaminants reach the air; results are compared with air samples
(RSL User's Guide 4.1.3, PDF p. 96).

Equations: modules/risk/equations/air.py (RSL UG 4.1.3, PDF pp. 96-101).
Verification: modules/risk/tests/test_resident_air.py.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from modules.risk.common import (ChemicalInputBase, ExposureOverrides, cancer_method, governing, label,
                                 load_defaults, one_hit, resolve_chemical)
from modules.risk.equations import air, factors
from ortam_core.model import ModelInput, ModelOutput, Module


class ChemicalInput(ChemicalInputBase):
    concentration_ug_m3: float | None = Field(
        None, ge=0, description="Exposure point concentration in air, ug/m3; if given, risk and HQ are calculated")


class ResidentAirInput(ModelInput):
    toxicity_table: Literal["rsl-2024-11", "none"] = Field(
        "rsl-2024-11", description="Toxicity source; 'none' requires all values in the input")
    target_risk: float = Field(1e-6, gt=0, lt=1, description="Target cancer risk TR, -")
    target_hq: float = Field(1.0, gt=0, description="Target hazard quotient THQ, -")
    exposure: ExposureOverrides = ExposureOverrides()
    chemicals: list[ChemicalInput] = Field(min_length=1)


class ChemicalResult(BaseModel):
    name: str
    cas: str
    cancer_method: str
    sl_nc_ug_m3: float | None = None
    sl_ca_ug_m3: float | None = None
    sl_ug_m3: float | None = Field(description="Screening level: lower of cancer and noncancer SL, ug/m3")
    sl_basis: Literal["c", "n"] | None = None
    concentration_ug_m3: float | None = None
    risk: float | None = Field(None, description="Cancer risk, - (one-hit rule above 0.01)")
    hq: float | None = None
    overridden: list[str] = []
    toxicity_keys: dict[str, str] = {}


class ResidentAirOutput(ModelOutput):
    toxicity_table: str
    target_risk: float
    target_hq: float
    exposure_used: dict[str, float]
    exposure_overridden: list[str] = []
    chemicals: list[ChemicalResult]
    total_risk: float | None = None
    hazard_index: float | None = None


class ResidentAirModule(Module):
    name = "risk.resident-air"
    title = "Risk: resident air screening levels, cancer risk and hazard"
    input_type = ResidentAirInput
    output_type = ResidentAirOutput

    def run(self, inputs: ResidentAirInput) -> ResidentAirOutput:
        ex = inputs.exposure.model_dump(exclude_none=True)
        d = load_defaults(overrides=ex)
        f = factors.compute(d)
        tr, thq = inputs.target_risk, inputs.target_hq
        warnings: list[str] = []
        if ex:
            warnings.append(f"Default exposure parameters overridden: {', '.join(sorted(ex))}. "
                            "Age-segment parameters (mutagen and TCE equations) are not changed.")
        results = []
        for ci in inputs.chemicals:
            c, overridden = resolve_chemical(inputs.toxicity_table, ci)
            if overridden:
                warnings.append(f"{label(c)}: table values overridden by the user: {', '.join(overridden)}.")
            method = cancer_method(c)
            nc = air.noncancer(c, d, thq)
            ca = air.cancer(c, d, f, tr, method)
            sl, basis = governing(ca, nc)
            if sl is None:
                warnings.append(f"{label(c)}: no inhalation toxicity values; no screening level calculated.")
            res = ChemicalResult(name=c.name, cas=c.cas, cancer_method=method, sl_nc_ug_m3=nc, sl_ca_ug_m3=ca,
                                 sl_ug_m3=sl, sl_basis=basis, overridden=overridden, toxicity_keys=dict(c.keys))
            if ci.concentration_ug_m3 is not None:
                C = res.concentration_ug_m3 = ci.concentration_ug_m3
                res.risk = one_hit(C * tr / ca) if ca else None   # RSL UG 2.6.1-2.6.2
                res.hq = C * thq / nc if nc else None
            results.append(res)
        with_c = [r for r in results if r.concentration_ug_m3 is not None]
        total_risk = hi = None
        if with_c:
            total_risk = one_hit(sum(r.concentration_ug_m3 * tr / r.sl_ca_ug_m3 for r in with_c if r.sl_ca_ug_m3))
            hi = sum(r.hq or 0.0 for r in with_c)
        return ResidentAirOutput(toxicity_table=inputs.toxicity_table, target_risk=tr, target_hq=thq,
                                 exposure_used=d.p, exposure_overridden=sorted(ex), chemicals=results,
                                 total_risk=total_risk, hazard_index=hi, warnings=warnings)

    def summary(self, o: ResidentAirOutput) -> str:
        lines = [f"Resident air  (toxicity: {o.toxicity_table}; TR = {o.target_risk:.0e}, THQ = {o.target_hq:g})",
                 f"{'Chemical':32s} {'CAS':12s} {'SL ug/m3':>10s} {'basis':>5s} {'C ug/m3':>10s} "
                 f"{'Risk':>9s} {'HQ':>9s}"]
        for c in o.chemicals:
            sl = f"{c.sl_ug_m3:.2g}" if c.sl_ug_m3 else "-"
            conc = f"{c.concentration_ug_m3:.3g}" if c.concentration_ug_m3 is not None else "-"
            risk = f"{c.risk:.1e}" if c.risk is not None else "-"
            hq = f"{c.hq:.2g}" if c.hq is not None else "-"
            lines.append(f"{c.name[:32]:32s} {c.cas:12s} {sl:>10s} {c.sl_basis or '-':>5s} {conc:>10s} "
                         f"{risk:>9s} {hq:>9s}{' *' if c.overridden else ''}")
        if o.total_risk is not None:
            lines.append(f"Total cancer risk : {o.total_risk:.1e}")
            lines.append(f"Hazard index      : {o.hazard_index:.2g}")
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
