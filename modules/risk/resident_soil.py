"""Resident soil: risk-based screening levels, cancer risk and hazard quotients.

Exposure routes: incidental soil ingestion, dermal contact with soil, and
inhalation of particulates (PEF) and, if a volatilization factor is supplied,
volatiles (VF).

Equations: modules/risk/equations/soil.py (RSL User's Guide 4.1.1, PDF pp. 36-56),
risk arithmetic: modules/risk/common.py (RSL UG 2.6; RAGS Part A ch. 8).
Default parameters: data/reference/rsl/2024-11/exposure_resident.yaml.
Cancer equations: standard, mutagenic (ADAF), vinyl chloride, trichloroethylene.
Verification: modules/risk/tests/test_resident_soil.py.

Not covered (an error rather than a wrong number):
  - volatile chemicals without a user-supplied VF: VF calculation is task D1.2.

Fortran analogy: `run` is the main program that loops over chemicals and
calls the PURE FUNCTIONs in equations/soil.py.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from modules.risk.common import (ChemicalInputBase, ExposureOverrides, Routes, cancer_method, governing,
                                 label, linear, load_defaults, one_hit, ratio, resolve_chemical, routes)
from modules.risk.equations import dose, factors, soil
from modules.risk.equations.dose import Dose
from ortam_core.model import ModelInput, ModelOutput, Module

CEILING_MG_KG = 1.0e5        # RSL UG 5.13: theoretical ceiling limit (flagged, not applied)


def age_adjusted_factors(p: dict[str, float]) -> tuple[float, float]:
    """IFSres-adj and DFSres-adj, mg/kg (RSL UG 4.1.1.3, PDF pp. 40-41)."""
    d = load_defaults(overrides=p)
    f = factors.compute(d)
    return f.ifs_mg_kg, f.dfs_mg_kg


def load_exposure_defaults(version: str = "2024-11") -> dict[str, float]:
    return load_defaults(version).p


class ChemicalInput(ChemicalInputBase):
    concentration_mg_kg: float | None = Field(
        None, ge=0, description="Exposure point concentration in soil, mg/kg (e.g. 95% UCL); "
                                "if given, risk and HQ are calculated")
    absd: float | None = Field(None, ge=0, le=1, description="Dermal absorption fraction from soil, -")
    rba: float | None = Field(None, gt=0, le=1, description="Relative bioavailability, -")
    vf_m3_kg: float | None = Field(None, gt=0, description="Volatilization factor, m3/kg (required for volatiles)")


class ResidentSoilInput(ModelInput):
    toxicity_table: Literal["rsl-2024-11", "none"] = Field(
        "rsl-2024-11", description="Toxicity source; 'none' requires all values in the input")
    target_risk: float = Field(1e-6, gt=0, lt=1, description="Target cancer risk TR, -")
    target_hq: float = Field(1.0, gt=0, description="Target hazard quotient THQ, -")
    exposure: ExposureOverrides = ExposureOverrides()
    chemicals: list[ChemicalInput] = Field(min_length=1)


class ChemicalResult(BaseModel):
    name: str
    cas: str
    rba: float
    cancer_method: str
    sl_nc_child_mg_kg: Routes
    sl_nc_adult_mg_kg: Routes
    sl_ca_mg_kg: Routes
    sl_mg_kg: float | None = Field(description="Screening level: lower of cancer and child noncancer SL, mg/kg")
    sl_basis: Literal["c", "n"] | None = None
    concentration_mg_kg: float | None = None
    risk: Routes | None = Field(None, description="Cancer risk, - (one-hit rule above 0.01)")
    hq_child: Routes | None = None
    hq_adult: Routes | None = None
    overridden: list[str] = []
    toxicity_keys: dict[str, str] = {}
    notes: list[str] = []
    doses: list[Dose] | None = Field(None, description="CDI/DAD (mg/kg-day) and EC (ug/m3) by route, basis and age segment")


class ResidentSoilOutput(ModelOutput):
    toxicity_table: str
    target_risk: float
    target_hq: float
    exposure_used: dict[str, float]
    exposure_overridden: list[str] = []
    chemicals: list[ChemicalResult]
    total_risk: float | None = Field(None, description="Total cancer risk, - (RAGS A p. 8-12; one-hit rule)")
    hazard_index_child: float | None = Field(None, description="Hazard index, child, - (RAGS A p. 8-13)")
    hazard_index_adult: float | None = None


class ResidentSoilModule(Module):
    name = "risk.resident-soil"
    title = "Risk: resident soil screening levels, cancer risk and hazard"
    input_type = ResidentSoilInput
    output_type = ResidentSoilOutput

    def run(self, inputs: ResidentSoilInput) -> ResidentSoilOutput:
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
            lab = label(c)
            if c.volatile and ci.vf_m3_kg is None:
                raise ValueError(f"{lab} is volatile; supply vf_m3_kg (VF calculation is task D1.2).")
            if ci.vf_m3_kg is not None:
                overridden.append("vf_m3_kg")
            rba = ci.rba if ci.rba is not None else d.rba_by_cas.get(c.cas, 1.0)
            if ci.rba is not None:
                overridden.append("rba")
            if overridden:
                warnings.append(f"{lab}: table values overridden by the user: {', '.join(overridden)}.")
            method = cancer_method(c)

            nc_c = routes(*soil.noncancer(c, rba, ci.vf_m3_kg, d, "child", thq))
            nc_a = routes(*soil.noncancer(c, rba, ci.vf_m3_kg, d, "adult", thq))
            ca = routes(*soil.cancer(c, rba, ci.vf_m3_kg, d, f, tr, method))
            sl, basis = governing(ca.total, nc_c.total)
            if sl is None:
                warnings.append(f"{lab}: no toxicity values; no screening level calculated.")
            elif sl > CEILING_MG_KG:
                warnings.append(f"{lab}: SL {sl:.3g} mg/kg exceeds the 1E+05 mg/kg ceiling limit "
                                "(RSL UG 5.13); flagged, not replaced.")
            if sl and c.csat_mg_kg and sl > c.csat_mg_kg:
                warnings.append(f"{lab}: SL exceeds Csat ({c.csat_mg_kg:.3g} mg/kg).")

            res = ChemicalResult(name=c.name, cas=c.cas, rba=rba, cancer_method=method, sl_nc_child_mg_kg=nc_c,
                                 sl_nc_adult_mg_kg=nc_a, sl_ca_mg_kg=ca, sl_mg_kg=sl, sl_basis=basis,
                                 overridden=overridden, toxicity_keys=dict(c.keys))
            if ci.concentration_mg_kg is not None:
                C = res.concentration_mg_kg = ci.concentration_mg_kg
                res.risk = ratio(C, tr, ca, risk=True)
                res.hq_child = ratio(C, thq, nc_c)
                res.hq_adult = ratio(C, thq, nc_a)
                res.doses = dose.resident_soil(C, c, rba, ci.vf_m3_kg, d, method)
            results.append(res)

        with_c = [r for r in results if r.concentration_mg_kg is not None]
        total_risk = hi_c = hi_a = None
        if with_c:
            total_risk = one_hit(sum(linear(r.concentration_mg_kg, tr, r.sl_ca_mg_kg.total) for r in with_c))
            hi_c = sum(r.hq_child.total or 0.0 for r in with_c)
            hi_a = sum(r.hq_adult.total or 0.0 for r in with_c)

        return ResidentSoilOutput(toxicity_table=inputs.toxicity_table, target_risk=tr, target_hq=thq,
                                  exposure_used=d.p, exposure_overridden=sorted(ex), chemicals=results,
                                  total_risk=total_risk, hazard_index_child=hi_c, hazard_index_adult=hi_a,
                                  warnings=warnings)

    def summary(self, o: ResidentSoilOutput) -> str:
        lines = [f"Resident soil  (toxicity: {o.toxicity_table}; TR = {o.target_risk:.0e}, THQ = {o.target_hq:g})",
                 f"{'Chemical':32s} {'CAS':12s} {'SL mg/kg':>10s} {'basis':>5s} {'C mg/kg':>10s} "
                 f"{'Risk':>9s} {'HQ child':>9s}"]
        for c in o.chemicals:
            sl = f"{c.sl_mg_kg:.2g}" if c.sl_mg_kg else "-"
            conc = f"{c.concentration_mg_kg:.3g}" if c.concentration_mg_kg is not None else "-"
            risk = f"{c.risk.total:.1e}" if c.risk and c.risk.total is not None else "-"
            hq = f"{c.hq_child.total:.2g}" if c.hq_child and c.hq_child.total is not None else "-"
            flag = " *" if c.overridden else ""
            lines.append(f"{c.name[:32]:32s} {c.cas:12s} {sl:>10s} {c.sl_basis or '-':>5s} {conc:>10s} "
                         f"{risk:>9s} {hq:>9s}{flag}")
        if o.total_risk is not None:
            lines.append(f"Total cancer risk    : {o.total_risk:.1e}")
            lines.append(f"Hazard index (child) : {o.hazard_index_child:.2g}")
            lines.append(f"Hazard index (adult) : {o.hazard_index_adult:.2g}")
        if any(c.overridden for c in o.chemicals) or o.exposure_overridden:
            lines.append("* user-supplied values replace table defaults (see warnings)")
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
