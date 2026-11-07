"""Soil-to-air volatilization factor (VF) and soil saturation limit (Csat).

Chronic (residential / worker) and subchronic (construction) VFs, each with
the unlimited-source or the mass-limit model.
Equations: modules/air/equations.py (RSL UG 4.9.4-4.9.7, 5.12; SSG 2002 Eq. 4-8,
4-9, 4-13, 5-14, 5-15, 5-17). Defaults: data/reference/rsl/2024-11/air_emission_factors.yaml.
Kd = Koc*foc for organics; table Kd for inorganics (RSL UG 4.9.4, p. 205).
VF is calculated only for volatile chemicals (RSL definition: HLC > 1E-05
atm-m3/mol or VP > 1 mmHg; RSL UG 4.9.4, p. 204).
The output VF is an input to the risk tools (risk.resident-soil: vf_m3_kg;
risk.construction-soil: vf_sc_m3_kg); the workflows connect them.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from modules.air import equations as eq
from modules.air.common import defaults
from ortam_core.chemicals.database import Chemical, lookup
from ortam_core.model import ModelInput, ModelOutput, Module

_PARAMS = ("henry_dimensionless", "dia_cm2_s", "diw_cm2_s", "koc_l_kg", "kd_l_kg", "solubility_mg_l", "mp_c")


class ChemicalInput(ModelInput):
    cas: str | None = None
    name: str | None = None
    henry_dimensionless: float | None = Field(None, ge=0, description="H', -")
    dia_cm2_s: float | None = Field(None, gt=0, description="Diffusivity in air, cm2/s")
    diw_cm2_s: float | None = Field(None, gt=0, description="Diffusivity in water, cm2/s")
    koc_l_kg: float | None = Field(None, ge=0, description="Koc, L/kg (organics)")
    kd_l_kg: float | None = Field(None, ge=0, description="Kd, L/kg (inorganics; overrides Koc*foc)")
    solubility_mg_l: float | None = Field(None, gt=0, description="Water solubility, mg/L")
    mp_c: float | None = Field(None, description="Melting point, C")
    volatile: bool | None = Field(None, description="Treat as volatile (default: RSL flag)")


class SoilProperties(ModelInput):
    rho_b_g_cm3: float | None = Field(None, gt=0, description="Dry bulk density, g/cm3 (default 1.5)")
    rho_s_g_cm3: float | None = Field(None, gt=0, description="Particle density, g/cm3 (default 2.65)")
    theta_w: float | None = Field(None, ge=0, lt=1, description="Water-filled porosity, - (default 0.15)")
    foc: float | None = Field(None, ge=0, lt=1, description="Organic carbon fraction, - (default 0.006)")


class Dispersion(ModelInput):
    """Either Q/C directly, or the constants A, B, C (SSG Appendix D / E) with the source area."""
    q_over_c: float | None = Field(None, gt=0, description="Q/C, g/m2-s per kg/m3")
    a: float | None = None
    b: float | None = None
    c: float | None = None
    area_acres: float | None = Field(None, ge=0.5, le=500, description="Source area, acres (0.5-500)")


class SslVfInput(ModelInput):
    scenario: Literal["chronic", "subchronic"] = Field("chronic", description="chronic: resident/worker; "
                                                       "subchronic: construction worker")
    source_model: Literal["unlimited", "mass-limit"] = "unlimited"
    source_depth_m: float | None = Field(None, gt=0, description="Average source depth ds, m (mass-limit model)")
    exposure_years: float | None = Field(None, gt=0, description="Chronic: T in yr (default 26); "
                                                                  "subchronic: ED in yr (default 1)")
    construction_weeks_per_yr: float | None = Field(None, gt=0, le=52, description="EW, wk/yr (default 50)")
    toxicity_table: Literal["rsl-2024-11", "none"] = "rsl-2024-11"
    soil: SoilProperties = SoilProperties()
    dispersion: Dispersion = Dispersion()
    chemicals: list[ChemicalInput] = Field(min_length=1)


class VfResult(BaseModel):
    name: str
    cas: str
    volatile: bool
    kd_l_kg: float | None = None
    da_cm2_s: float | None = None
    vf_m3_kg: float | None = None
    csat_mg_kg: float | None = None
    liquid: bool | None = Field(None, description="Liquid at soil temperature (MP < 20 C, RSL UG 5.12)")
    overridden: list[str] = []


class SslVfOutput(ModelOutput):
    scenario: str
    source_model: str
    q_over_c: float
    exposure_interval_s: float
    dispersion_correction: float | None = None
    chemicals: list[VfResult]


class SslVfModule(Module):
    name = "air.ssl-vf"
    title = "Air: soil-to-air volatilization factor (VF) and Csat (EPA SSL)"
    input_type = SslVfInput
    output_type = SslVfOutput

    def run(self, inputs: SslVfInput) -> SslVfOutput:
        dflt = defaults()
        soil = {**dflt["soil"], **inputs.soil.model_dump(exclude_none=True)}
        rb, rs, tw, foc = soil["rho_b_g_cm3"], soil["rho_s_g_cm3"], soil["theta_w"], soil["foc"]
        warnings: list[str] = []
        if inputs.source_model == "mass-limit" and not inputs.source_depth_m:
            raise ValueError("The mass-limit model needs source_depth_m (RSL UG 4.9.5).")

        disp = inputs.dispersion
        fd = None
        if inputs.scenario == "chronic":
            t_yr = inputs.exposure_years or dflt["vf_chronic"]["t_yr_mass_limit"]
            t_s = dflt["vf_chronic"]["t_s"] if inputs.exposure_years is None else t_yr * eq.SECONDS_PER_YEAR
            qc = disp.q_over_c or (eq.q_over_c(disp.a, disp.b, disp.c, disp.area_acres)
                                   if disp.a is not None else dflt["vf_chronic"]["q_over_c"])
        else:
            con = dflt["construction"]
            ed = inputs.exposure_years or con["ed_yr"]
            ew = inputs.construction_weeks_per_yr or con["ew_wk_yr"]
            tc_h = ed * ew * 7.0 * 24.0                      # RSL UG 4.9.6, p. 214
            t_s = tc_h * 3600.0
            fd = eq.dispersion_correction(tc_h)
            sc = dflt["vf_subchronic"]
            qc = disp.q_over_c or eq.q_over_c(disp.a if disp.a is not None else sc["a"],
                                               disp.b if disp.b is not None else sc["b"],
                                               disp.c if disp.c is not None else sc["c"],
                                               disp.area_acres or 0.5)
        if disp.q_over_c or disp.a is not None:
            warnings.append(f"Site-specific dispersion factor used: Q/C = {qc:.4g}.")

        results = []
        for ci in inputs.chemicals:
            given = {k: getattr(ci, k) for k in _PARAMS if getattr(ci, k) is not None}
            if inputs.toxicity_table == "none":
                c = Chemical(name=ci.name or "", cas=ci.cas or "", volatile=True, **given)
                overridden = []
            else:
                base = lookup(inputs.toxicity_table, ci.cas, ci.name)
                c = Chemical(**{**base.__dict__, **given})
                overridden = sorted(given)
            volatile = ci.volatile if ci.volatile is not None else c.volatile
            label = f"{c.name} ({c.cas})" if c.cas else c.name
            kd = c.kd_l_kg if c.kd_l_kg is not None else (c.koc_l_kg * foc if c.koc_l_kg is not None else None)
            r = VfResult(name=c.name, cas=c.cas, volatile=volatile, kd_l_kg=kd, overridden=overridden,
                         liquid=(c.mp_c < 20.0) if c.mp_c is not None else None)
            if not volatile:
                warnings.append(f"{label}: not volatile (RSL definition); no VF calculated.")
                results.append(r)
                continue
            missing = [k for k, v in (("H'", c.henry_dimensionless), ("Dia", c.dia_cm2_s), ("Diw", c.diw_cm2_s),
                                      ("Kd/Koc", kd)) if v is None]
            if missing:
                warnings.append(f"{label}: missing {', '.join(missing)}; no VF calculated.")
                results.append(r)
                continue
            r.da_cm2_s = eq.apparent_diffusivity(c.henry_dimensionless, c.dia_cm2_s, c.diw_cm2_s, kd, rb, rs, tw)
            if inputs.scenario == "chronic":
                r.vf_m3_kg = (eq.vf_unlimited(qc, r.da_cm2_s, t_s, rb) if inputs.source_model == "unlimited"
                              else eq.vf_mass_limit(qc, t_yr, rb, inputs.source_depth_m))
            else:
                r.vf_m3_kg = (eq.vf_subchronic_unlimited(qc, fd, r.da_cm2_s, t_s, rb)
                              if inputs.source_model == "unlimited"
                              else eq.vf_subchronic_mass_limit(qc, fd, t_s, rb, inputs.source_depth_m))
            if c.solubility_mg_l is not None:
                r.csat_mg_kg = eq.csat(c.solubility_mg_l, kd, c.henry_dimensionless, rb, rs, tw)
            results.append(r)
        return SslVfOutput(scenario=inputs.scenario, source_model=inputs.source_model, q_over_c=qc,
                           exposure_interval_s=t_s, dispersion_correction=fd, chemicals=results, warnings=warnings)

    def summary(self, o: SslVfOutput) -> str:
        lines = [f"Volatilization factor ({o.scenario}, {o.source_model}; Q/C = {o.q_over_c:.4g})",
                 f"{'Chemical':32s} {'CAS':12s} {'VF m3/kg':>11s} {'Csat mg/kg':>11s} {'liquid':>7s}"]
        for c in o.chemicals:
            lines.append(f"{c.name[:32]:32s} {c.cas:12s} {(f'{c.vf_m3_kg:.4g}' if c.vf_m3_kg else '-'):>11s} "
                         f"{(f'{c.csat_mg_kg:.4g}' if c.csat_mg_kg else '-'):>11s} "
                         f"{('yes' if c.liquid else 'no') if c.liquid is not None else '-':>7s}")
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
