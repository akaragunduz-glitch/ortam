"""Recreator soil or sediment: screening levels, cancer risk and hazard.

The recreator equations (RSL UG 4.6.1, PDF pp. 124-153) have the resident
form (4.1.1) with site-specific recreator parameters; child/adult values are
ED-weighted means of the age segments (4.6.1.7). The recreator parameters are
mapped onto the resident structure (equations/recreator.py) and the verified
resident soil equations are reused (equations/soil.py), including the
mutagenic, vinyl chloride and TCE sets.
Site-specific inputs (no RSL defaults): EF, ET, IRS, SA, AF per age segment.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Literal

import yaml
from pydantic import BaseModel, Field

from modules.risk.common import (cancer_method, governing, label, linear, load_defaults, one_hit, ratio,
                                 resolve_chemical, routes)
from modules.risk.equations import dose, factors, recreator, soil
from modules.risk.resident_soil import ChemicalInput, ChemicalResult
from ortam_core.chemicals.rsl_import import DATA
from ortam_core.model import ModelInput, ModelOutput, Module

Four = list[float]


@lru_cache(maxsize=None)
def recreator_defaults() -> dict[str, list[float]]:
    f = yaml.safe_load((DATA / "2024-11" / "exposure_recreator.yaml").read_text(encoding="utf-8"))
    return {k: [float(x) for x in v["values"]] for k, v in f["defaults"].items()}


def four(desc: str, required: bool = True):
    return Field(... if required else None, min_length=4, max_length=4,
                 description=f"{desc}; four values for ages 0-2, 2-6, 6-16, 16-26")


class RecreatorSoilExposure(ModelInput):
    ef_d_yr: Four = four("Exposure frequency, day/yr (site-specific)")
    et_h_d: Four = four("Exposure time, h/day (site-specific)")
    irs_mg_d: Four = four("Soil or sediment ingestion rate, mg/day (site-specific)")
    sa_soil_cm2_d: Four = four("Skin surface area, cm2/day (site-specific)")
    af_mg_cm2: Four = four("Soil adherence factor, mg/cm2 (site-specific)")
    ed_yr: Four | None = four("Exposure duration, yr (default 2, 4, 10, 10)", False)
    bw_kg: Four | None = four("Body weight, kg (default 15, 15, 80, 80)", False)


class RecreatorSoilInput(ModelInput):
    medium: Literal["soil", "sediment"] = "soil"
    toxicity_table: Literal["rsl-2024-11", "none"] = "rsl-2024-11"
    target_risk: float = Field(1e-6, gt=0, lt=1, description="Target cancer risk TR, -")
    target_hq: float = Field(1.0, gt=0, description="Target hazard quotient THQ, -")
    pef_m3_kg: float | None = Field(None, gt=0, description="Particulate emission factor, m3/kg (default 1.36E+09)")
    exposure: RecreatorSoilExposure
    chemicals: list[ChemicalInput] = Field(min_length=1)


class RecreatorSoilOutput(ModelOutput):
    medium: str
    toxicity_table: str
    exposure_used: dict[str, float]
    chemicals: list[ChemicalResult]
    total_risk: float | None = None
    hazard_index_child: float | None = None
    hazard_index_adult: float | None = None


def build_defaults(ex: dict, extra: dict | None = None, pef: float | None = None):
    seg = {**recreator_defaults(), **{k: v for k, v in ex.items() if v is not None}}
    base = load_defaults(overrides={"pef_m3_kg": pef} if pef else None)
    d = recreator.build(base, seg)
    return d


class RecreatorSoilModule(Module):
    name = "risk.recreator-soil"
    title = "Risk: recreator soil or sediment screening levels, risk and hazard"
    input_type = RecreatorSoilInput
    output_type = RecreatorSoilOutput

    def run(self, inputs: RecreatorSoilInput) -> RecreatorSoilOutput:
        d = build_defaults(inputs.exposure.model_dump(), pef=inputs.pef_m3_kg)
        f = factors.compute(d)
        tr, thq = inputs.target_risk, inputs.target_hq
        warnings: list[str] = []
        results = []
        for ci in inputs.chemicals:
            c, overridden = resolve_chemical(inputs.toxicity_table, ci)
            if c.volatile and ci.vf_m3_kg is None:
                raise ValueError(f"{label(c)} is volatile; supply vf_m3_kg (VF calculation is task D1.2).")
            if ci.vf_m3_kg is not None:
                overridden.append("vf_m3_kg")
            rba = ci.rba if ci.rba is not None else d.rba_by_cas.get(c.cas, 1.0)
            if ci.rba is not None:
                overridden.append("rba")
            if overridden:
                warnings.append(f"{label(c)}: table values overridden by the user: {', '.join(overridden)}.")
            method = cancer_method(c)
            nc_c = routes(*soil.noncancer(c, rba, ci.vf_m3_kg, d, "child", thq))
            nc_a = routes(*soil.noncancer(c, rba, ci.vf_m3_kg, d, "adult", thq))
            ca = routes(*soil.cancer(c, rba, ci.vf_m3_kg, d, f, tr, method))
            sl, basis = governing(ca.total, nc_c.total)
            r = ChemicalResult(name=c.name, cas=c.cas, rba=rba, cancer_method=method, sl_nc_child_mg_kg=nc_c,
                               sl_nc_adult_mg_kg=nc_a, sl_ca_mg_kg=ca, sl_mg_kg=sl, sl_basis=basis,
                               overridden=overridden, toxicity_keys=dict(c.keys))
            if ci.concentration_mg_kg is not None:
                C = r.concentration_mg_kg = ci.concentration_mg_kg
                r.risk, r.hq_child, r.hq_adult = ratio(C, tr, ca, True), ratio(C, thq, nc_c), ratio(C, thq, nc_a)
                r.doses = dose.resident_soil(C, c, rba, ci.vf_m3_kg, d, method)
            results.append(r)
        with_c = [r for r in results if r.concentration_mg_kg is not None]
        total = hc = ha = None
        if with_c:
            total = one_hit(sum(linear(r.concentration_mg_kg, tr, r.sl_ca_mg_kg.total) for r in with_c))
            hc = sum(r.hq_child.total or 0.0 for r in with_c)
            ha = sum(r.hq_adult.total or 0.0 for r in with_c)
        return RecreatorSoilOutput(medium=inputs.medium, toxicity_table=inputs.toxicity_table, exposure_used=d.p,
                                   chemicals=results, total_risk=total, hazard_index_child=hc,
                                   hazard_index_adult=ha, warnings=warnings)

    def summary(self, o: RecreatorSoilOutput) -> str:
        lines = [f"Recreator {o.medium} (toxicity: {o.toxicity_table})",
                 f"{'Chemical':32s} {'SL mg/kg':>10s} {'basis':>5s} {'Risk':>9s} {'HQ child':>9s}"]
        for c in o.chemicals:
            lines.append(f"{c.name[:32]:32s} {(f'{c.sl_mg_kg:.2g}' if c.sl_mg_kg else '-'):>10s} "
                         f"{c.sl_basis or '-':>5s} "
                         f"{(f'{c.risk.total:.1e}' if c.risk and c.risk.total is not None else '-'):>9s} "
                         f"{(f'{c.hq_child.total:.2g}' if c.hq_child and c.hq_child.total is not None else '-'):>9s}")
        if o.total_risk is not None:
            lines += [f"Total cancer risk    : {o.total_risk:.1e}",
                      f"Hazard index (child) : {o.hazard_index_child:.2g}"]
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
