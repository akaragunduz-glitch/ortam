"""Resident tap water: risk-based screening levels, cancer risk and hazard quotients.

Exposure routes: ingestion of drinking water, dermal contact while bathing
(chemicals inside the RAGS Part E effective predictive domain; inorganics with
Kp), and inhalation of volatiles during household water use (volatile chemicals).

Equations: modules/risk/equations/tapwater.py and dermal_water.py
(RSL User's Guide 4.1.2, PDF pp. 62-95; 4.9.8, pp. 217-225; RAGS Part E p. 3-7).
Verification: modules/risk/tests/test_resident_tapwater.py.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from modules.risk.common import (ChemicalInputBase, ExposureOverrides, Routes, cancer_method, governing,
                                 label, linear, load_defaults, one_hit, ratio, resolve_chemical, routes)
from modules.risk.equations import dermal_water, dose, factors, tapwater
from modules.risk.equations.dose import Dose
from ortam_core.model import ModelInput, ModelOutput, Module


class ChemicalInput(ChemicalInputBase):
    concentration_ug_l: float | None = Field(
        None, ge=0, description="Exposure point concentration in tap water, ug/L; "
                                "if given, risk and HQ are calculated")
    mw_g_mol: float | None = Field(None, gt=0, description="Molecular weight, g/mol")
    log_kow: float | None = Field(None, description="log Kow, -")
    kp_cm_h: float | None = Field(None, gt=0, description="Dermal permeability constant, cm/h")
    b_dimensionless: float | None = Field(None, gt=0, description="Permeability ratio B, -")
    tau_event_h: float | None = Field(None, gt=0, description="Lag time per event, h/event")
    t_star_h: float | None = Field(None, gt=0, description="Time to reach steady state, h")


class ResidentTapwaterInput(ModelInput):
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
    dermal_fa: float | None = Field(None, description="Fraction absorbed (organics), -")
    sl_nc_child_ug_l: Routes
    sl_nc_adult_ug_l: Routes
    sl_ca_ug_l: Routes
    sl_ug_l: float | None = Field(description="Screening level: lower of cancer and child noncancer SL, ug/L")
    sl_basis: Literal["c", "n"] | None = None
    concentration_ug_l: float | None = None
    risk: Routes | None = None
    hq_child: Routes | None = None
    hq_adult: Routes | None = None
    overridden: list[str] = []
    toxicity_keys: dict[str, str] = {}
    notes: list[str] = []
    doses: list[Dose] | None = Field(None, description="CDI/DAD (mg/kg-day) and EC (ug/m3) by route, basis and age segment")


class ResidentTapwaterOutput(ModelOutput):
    toxicity_table: str
    target_risk: float
    target_hq: float
    exposure_used: dict[str, float]
    exposure_overridden: list[str] = []
    chemicals: list[ChemicalResult]
    total_risk: float | None = None
    hazard_index_child: float | None = None
    hazard_index_adult: float | None = None


class ResidentTapwaterModule(Module):
    name = "risk.resident-tapwater"
    title = "Risk: resident tap water screening levels, cancer risk and hazard"
    input_type = ResidentTapwaterInput
    output_type = ResidentTapwaterOutput

    def run(self, inputs: ResidentTapwaterInput) -> ResidentTapwaterOutput:
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
            applies, inorganic, fa, note = dermal_water.applicability(c)
            nc_c = routes(*tapwater.noncancer(c, d, "child", thq))
            nc_a = routes(*tapwater.noncancer(c, d, "adult", thq))
            ca = routes(*tapwater.cancer(c, d, f, tr, method))
            sl, basis = governing(ca.total, nc_c.total)
            if sl is None:
                warnings.append(f"{label(c)}: no toxicity values; no screening level calculated.")
            res = ChemicalResult(name=c.name, cas=c.cas, cancer_method=method, dermal_fa=fa,
                                 sl_nc_child_ug_l=nc_c, sl_nc_adult_ug_l=nc_a, sl_ca_ug_l=ca, sl_ug_l=sl,
                                 sl_basis=basis, overridden=overridden, toxicity_keys=dict(c.keys),
                                 notes=[note] if note else [])
            if ci.concentration_ug_l is not None:
                C = res.concentration_ug_l = ci.concentration_ug_l
                res.risk = ratio(C, tr, ca, risk=True)
                res.hq_child = ratio(C, thq, nc_c)
                res.hq_adult = ratio(C, thq, nc_a)
                res.doses = dose.resident_water(C, c, d, f, method)
            results.append(res)

        with_c = [r for r in results if r.concentration_ug_l is not None]
        total_risk = hi_c = hi_a = None
        if with_c:
            total_risk = one_hit(sum(linear(r.concentration_ug_l, tr, r.sl_ca_ug_l.total) for r in with_c))
            hi_c = sum(r.hq_child.total or 0.0 for r in with_c)
            hi_a = sum(r.hq_adult.total or 0.0 for r in with_c)
        return ResidentTapwaterOutput(toxicity_table=inputs.toxicity_table, target_risk=tr, target_hq=thq,
                                      exposure_used=d.p, exposure_overridden=sorted(ex), chemicals=results,
                                      total_risk=total_risk, hazard_index_child=hi_c, hazard_index_adult=hi_a,
                                      warnings=warnings)

    def summary(self, o: ResidentTapwaterOutput) -> str:
        lines = [f"Resident tap water  (toxicity: {o.toxicity_table}; TR = {o.target_risk:.0e}, "
                 f"THQ = {o.target_hq:g})",
                 f"{'Chemical':32s} {'CAS':12s} {'SL ug/L':>10s} {'basis':>5s} {'C ug/L':>10s} "
                 f"{'Risk':>9s} {'HQ child':>9s}"]
        for c in o.chemicals:
            sl = f"{c.sl_ug_l:.2g}" if c.sl_ug_l else "-"
            conc = f"{c.concentration_ug_l:.3g}" if c.concentration_ug_l is not None else "-"
            risk = f"{c.risk.total:.1e}" if c.risk and c.risk.total is not None else "-"
            hq = f"{c.hq_child.total:.2g}" if c.hq_child and c.hq_child.total is not None else "-"
            lines.append(f"{c.name[:32]:32s} {c.cas:12s} {sl:>10s} {c.sl_basis or '-':>5s} {conc:>10s} "
                         f"{risk:>9s} {hq:>9s}{' *' if c.overridden else ''}")
        if o.total_risk is not None:
            lines.append(f"Total cancer risk    : {o.total_risk:.1e}")
            lines.append(f"Hazard index (child) : {o.hazard_index_child:.2g}")
            lines.append(f"Hazard index (adult) : {o.hazard_index_adult:.2g}")
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
