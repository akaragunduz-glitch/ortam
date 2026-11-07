"""Ingestion of fish: screening levels for fish tissue, cancer risk and hazard.

RSL UG 4.7 (PDF pp. 182-183). The SL is a concentration in fish tissue
(mg/kg), not in water or sediment; the fish ingestion rate is not age-adjusted.
  nc  SL = THQ*AT*BWres-a / ((1/RfD)*1e-6*EFres*EDres*IRF),  AT = 365*EDres   (4.7.1)
  ca  SL = TR*AT*BWres-a / (SFO*1e-6*EFres*EDres*IRF),        AT = 365*LT      (4.7.2)
IRF (mg/day) is site-specific (RSL Default Exposure Parameters: IRFres-a).
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from modules.risk.common import ChemicalInputBase, governing, label, load_defaults, one_hit, resolve_chemical
from ortam_core.model import ModelInput, ModelOutput, Module


class ChemicalInput(ChemicalInputBase):
    concentration_mg_kg: float | None = Field(None, ge=0, description="Fish tissue concentration, mg/kg")


class FishInput(ModelInput):
    irf_mg_d: float = Field(gt=0, description="Fish ingestion rate, mg/day (site-specific)")
    toxicity_table: Literal["rsl-2024-11", "none"] = "rsl-2024-11"
    target_risk: float = Field(1e-6, gt=0, lt=1, description="Target cancer risk TR, -")
    target_hq: float = Field(1.0, gt=0, description="Target hazard quotient THQ, -")
    ef_d_yr: float | None = Field(None, gt=0, le=365, description="Exposure frequency, day/yr (default 350)")
    ed_yr: float | None = Field(None, gt=0, description="Exposure duration, yr (default 26)")
    bw_kg: float | None = Field(None, gt=0, description="Body weight, kg (default 80)")
    chemicals: list[ChemicalInput] = Field(min_length=1)


class FishResult(BaseModel):
    name: str
    cas: str
    sl_nc_mg_kg: float | None = None
    sl_ca_mg_kg: float | None = None
    sl_mg_kg: float | None = None
    sl_basis: Literal["c", "n"] | None = None
    concentration_mg_kg: float | None = None
    risk: float | None = None
    hq: float | None = None
    overridden: list[str] = []


class FishOutput(ModelOutput):
    irf_mg_d: float
    chemicals: list[FishResult]
    total_risk: float | None = None
    hazard_index: float | None = None


class FishModule(Module):
    name = "risk.fish"
    title = "Risk: ingestion of fish (tissue screening levels, risk and hazard)"
    input_type = FishInput
    output_type = FishOutput

    def run(self, inputs: FishInput) -> FishOutput:
        p = load_defaults().p
        ef = inputs.ef_d_yr or p["ef_res_d_yr"]
        ed = inputs.ed_yr or p["ed_res_yr"]
        bw = inputs.bw_kg or p["bw_adult_kg"]
        irf, tr, thq = inputs.irf_mg_d, inputs.target_risk, inputs.target_hq
        warnings, results = [], []
        for ci in inputs.chemicals:
            c, overridden = resolve_chemical(inputs.toxicity_table, ci)
            if overridden:
                warnings.append(f"{label(c)}: table values overridden by the user: {', '.join(overridden)}.")
            rfd, sfo = c.rfd_mg_kg_day, c.sfo_per_mg_kg_day
            nc = thq * 365.0 * ed * bw / ((1.0 / rfd) * 1e-6 * ef * ed * irf) if rfd else None
            ca = tr * 365.0 * p["lt_yr"] * bw / (sfo * 1e-6 * ef * ed * irf) if sfo else None
            sl, basis = governing(ca, nc)
            r = FishResult(name=c.name, cas=c.cas, sl_nc_mg_kg=nc, sl_ca_mg_kg=ca, sl_mg_kg=sl, sl_basis=basis,
                           overridden=overridden)
            if ci.concentration_mg_kg is not None:
                C = r.concentration_mg_kg = ci.concentration_mg_kg
                r.risk = one_hit(C * tr / ca) if ca else None
                r.hq = C * thq / nc if nc else None
            results.append(r)
        with_c = [r for r in results if r.concentration_mg_kg is not None]
        total = hi = None
        if with_c:
            total = one_hit(sum(r.concentration_mg_kg * tr / r.sl_ca_mg_kg for r in with_c if r.sl_ca_mg_kg))
            hi = sum(r.hq or 0.0 for r in with_c)
        return FishOutput(irf_mg_d=irf, chemicals=results, total_risk=total, hazard_index=hi, warnings=warnings)

    def summary(self, o: FishOutput) -> str:
        lines = [f"Fish ingestion (IRF = {o.irf_mg_d:g} mg/day)",
                 f"{'Chemical':32s} {'SL mg/kg':>10s} {'basis':>5s} {'Risk':>9s} {'HQ':>9s}"]
        for c in o.chemicals:
            lines.append(f"{c.name[:32]:32s} {(f'{c.sl_mg_kg:.2g}' if c.sl_mg_kg else '-'):>10s} "
                         f"{c.sl_basis or '-':>5s} {(f'{c.risk:.1e}' if c.risk is not None else '-'):>9s} "
                         f"{(f'{c.hq:.2g}' if c.hq is not None else '-'):>9s}")
        if o.total_risk is not None:
            lines += [f"Total cancer risk : {o.total_risk:.1e}", f"Hazard index      : {o.hazard_index:.2g}"]
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
