"""Recreator surface water: screening levels, cancer risk and hazard.

Routes: incidental ingestion and dermal contact (no inhalation: RSL UG 4.6.2,
PDF p. 153). The equations (4.6.2, pp. 153-182) have the resident tap water
form (4.1.2) with IRW_day = EV * ET_event * IRW (L/h); recreator parameters
are mapped onto the resident structure (equations/recreator.py) and the
verified tap water equations are reused (equations/tapwater.py).
Site-specific inputs (no RSL defaults): EF, EV, ET_event per age segment.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Literal

from pydantic import Field

from modules.risk.common import cancer_method, governing, label, linear, one_hit, ratio, resolve_chemical, routes
from modules.risk.equations import dermal_water, factors, tapwater
from modules.risk.recreator_soil import Four, build_defaults, four
from modules.risk.resident_tapwater import ChemicalInput, ChemicalResult
from ortam_core.model import ModelInput, ModelOutput, Module


class RecreatorWaterExposure(ModelInput):
    ef_d_yr: Four = four("Exposure frequency, day/yr (site-specific)")
    ev_per_d: Four = four("Events per day, event/day (site-specific)")
    et_event_h: Four = four("Exposure time per event, h/event (site-specific)")
    irw_l_h: Four | None = four("Water ingestion rate, L/h (default 0.12, 0.12, 0.124, 0.0985)", False)
    sa_water_cm2: Four | None = four("Skin surface area, cm2 (default 6365, 6365, 19652, 19652)", False)
    ed_yr: Four | None = four("Exposure duration, yr (default 2, 4, 10, 10)", False)
    bw_kg: Four | None = four("Body weight, kg (default 15, 15, 80, 80)", False)


class RecreatorSurfaceWaterInput(ModelInput):
    toxicity_table: Literal["rsl-2024-11", "none"] = "rsl-2024-11"
    target_risk: float = Field(1e-6, gt=0, lt=1, description="Target cancer risk TR, -")
    target_hq: float = Field(1.0, gt=0, description="Target hazard quotient THQ, -")
    exposure: RecreatorWaterExposure
    chemicals: list[ChemicalInput] = Field(min_length=1)


class RecreatorSurfaceWaterOutput(ModelOutput):
    toxicity_table: str
    exposure_used: dict[str, float]
    chemicals: list[ChemicalResult]
    total_risk: float | None = None
    hazard_index_child: float | None = None
    hazard_index_adult: float | None = None


class RecreatorSurfaceWaterModule(Module):
    name = "risk.recreator-surface-water"
    title = "Risk: recreator surface water screening levels, risk and hazard"
    input_type = RecreatorSurfaceWaterInput
    output_type = RecreatorSurfaceWaterOutput

    def run(self, inputs: RecreatorSurfaceWaterInput) -> RecreatorSurfaceWaterOutput:
        d = build_defaults(inputs.exposure.model_dump())
        f = factors.compute(d)
        tr, thq = inputs.target_risk, inputs.target_hq
        warnings: list[str] = []
        results = []
        for ci in inputs.chemicals:
            c, overridden = resolve_chemical(inputs.toxicity_table, ci)
            if overridden:
                warnings.append(f"{label(c)}: table values overridden by the user: {', '.join(overridden)}.")
            method = cancer_method(c)
            cw = replace(c, volatile=False)   # no inhalation route for surface water (4.6.2, p. 153)
            _, _, fa, note = dermal_water.applicability(c)
            nc_c = routes(*tapwater.noncancer(cw, d, "child", thq))
            nc_a = routes(*tapwater.noncancer(cw, d, "adult", thq))
            ca = routes(*tapwater.cancer(cw, d, f, tr, method))
            sl, basis = governing(ca.total, nc_c.total)
            r = ChemicalResult(name=c.name, cas=c.cas, cancer_method=method, dermal_fa=fa, sl_nc_child_ug_l=nc_c,
                               sl_nc_adult_ug_l=nc_a, sl_ca_ug_l=ca, sl_ug_l=sl, sl_basis=basis,
                               overridden=overridden, toxicity_keys=dict(c.keys), notes=[note] if note else [])
            if ci.concentration_ug_l is not None:
                C = r.concentration_ug_l = ci.concentration_ug_l
                r.risk, r.hq_child, r.hq_adult = ratio(C, tr, ca, True), ratio(C, thq, nc_c), ratio(C, thq, nc_a)
            results.append(r)
        with_c = [r for r in results if r.concentration_ug_l is not None]
        total = hc = ha = None
        if with_c:
            total = one_hit(sum(linear(r.concentration_ug_l, tr, r.sl_ca_ug_l.total) for r in with_c))
            hc = sum(r.hq_child.total or 0.0 for r in with_c)
            ha = sum(r.hq_adult.total or 0.0 for r in with_c)
        return RecreatorSurfaceWaterOutput(toxicity_table=inputs.toxicity_table, exposure_used=d.p,
                                           chemicals=results, total_risk=total, hazard_index_child=hc,
                                           hazard_index_adult=ha, warnings=warnings)

    def summary(self, o: RecreatorSurfaceWaterOutput) -> str:
        lines = [f"Recreator surface water (toxicity: {o.toxicity_table})",
                 f"{'Chemical':32s} {'SL ug/L':>10s} {'basis':>5s} {'Risk':>9s} {'HQ child':>9s}"]
        for c in o.chemicals:
            lines.append(f"{c.name[:32]:32s} {(f'{c.sl_ug_l:.2g}' if c.sl_ug_l else '-'):>10s} "
                         f"{c.sl_basis or '-':>5s} "
                         f"{(f'{c.risk.total:.1e}' if c.risk and c.risk.total is not None else '-'):>9s} "
                         f"{(f'{c.hq_child.total:.2g}' if c.hq_child and c.hq_child.total is not None else '-'):>9s}")
        if o.total_risk is not None:
            lines += [f"Total cancer risk    : {o.total_risk:.1e}",
                      f"Hazard index (child) : {o.hazard_index_child:.2g}"]
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
