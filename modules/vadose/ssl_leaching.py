"""Soil screening level protective of groundwater (SSL): migration to groundwater.

Method 1 (partitioning, infinite source) and Method 2 (mass-limit), with a
default, fixed or site-specific dilution attenuation factor (DAF).
Equations: modules/vadose/equations.py (RSL UG 4.8; SSG 2002 Eq. 4-10 to 4-14).
Defaults: data/reference/rsl/2024-11/soil_to_groundwater.yaml.

Target groundwater concentration (per chemical): a user value in ug/L
(e.g. a risk-based tap water SL from risk.resident-tapwater, connected in
workflows/soil_to_groundwater) or the MCL from the RSL table (target: mcl).
Kd: table Kd for inorganics (SSG 1996 Exhibit C-4 at pH 6.8, or Baes 1984;
RSL UG 4.8, p. 185), Koc*foc for organics. H' = 0 for inorganics except
mercury (SSG Eq. 4-10, note). No Csat substitution (RSL UG 4.8, p. 184).
"""

from __future__ import annotations

from functools import lru_cache
from typing import Literal

import yaml
from pydantic import BaseModel, Field, model_validator

from modules.vadose import equations as eq
from ortam_core.chemicals.database import Chemical, lookup
from ortam_core.chemicals.rsl_import import DATA
from ortam_core.model import ModelInput, ModelOutput, Module

MERCURY_CAS = "7439-97-6"


@lru_cache(maxsize=None)
def defaults() -> dict[str, float]:
    f = yaml.safe_load((DATA / "2024-11" / "soil_to_groundwater.yaml").read_text(encoding="utf-8"))
    return {k: float(v["value"]) for k, v in f["defaults"].items()}


class ChemicalInput(ModelInput):
    cas: str | None = None
    name: str | None = None
    target_ug_l: float | None = Field(None, gt=0, description="Target groundwater concentration, ug/L")
    target: Literal["mcl"] | None = Field(None, description="'mcl': use the MCL from the RSL table")
    kd_l_kg: float | None = Field(None, ge=0, description="Kd override, L/kg")
    koc_l_kg: float | None = Field(None, ge=0, description="Koc override, L/kg")
    henry_dimensionless: float | None = Field(None, ge=0, description="H' override, -")
    inorganic: bool | None = Field(None, description="Override the inorganic classification (no log Kow)")

    @model_validator(mode="after")
    def _one_target(self):
        if (self.target_ug_l is None) == (self.target is None):
            raise ValueError("Give exactly one of target_ug_l or target: mcl.")
        return self


class Site(ModelInput):
    """Site-specific DAF inputs (SSG Eq. 4-11, 4-12); all four are needed."""
    hydraulic_conductivity_m_yr: float = Field(gt=0, description="Aquifer hydraulic conductivity K, m/yr")
    hydraulic_gradient: float = Field(gt=0, description="Hydraulic gradient i, m/m")
    source_length_m: float = Field(gt=0, description="Source length parallel to flow L, m")
    aquifer_thickness_m: float = Field(gt=0, description="Aquifer thickness da, m")


class SslLeachingInput(ModelInput):
    method: Literal["partitioning", "mass-limit"] = "partitioning"
    source_depth_m: float | None = Field(None, gt=0, description="Source depth ds, m (mass-limit)")
    daf: float | None = Field(None, ge=1, description="Dilution attenuation factor (default 1; 20 for a "
                                                      "0.5 acre source); ignored if 'site' is given")
    site: Site | None = None
    theta_w: float | None = Field(None, gt=0, lt=1, description="Water-filled porosity, - (default 0.3)")
    rho_b_kg_l: float | None = Field(None, gt=0, description="Dry bulk density, kg/L (default 1.5)")
    rho_s_kg_l: float | None = Field(None, gt=0, description="Particle density, kg/L (default 2.65)")
    foc: float | None = Field(None, ge=0, lt=1, description="Organic carbon fraction, - (default 0.002)")
    infiltration_m_yr: float | None = Field(None, gt=0, description="Infiltration rate I, m/yr (default 0.18)")
    ed_yr: float | None = Field(None, gt=0, description="Exposure duration, yr (mass-limit; default 70)")
    toxicity_table: Literal["rsl-2024-11", "none"] = "rsl-2024-11"
    chemicals: list[ChemicalInput] = Field(min_length=1)


class SslResult(BaseModel):
    name: str
    cas: str
    target_ug_l: float | None = None
    target_basis: str
    cw_mg_l: float | None = None
    kd_l_kg: float | None = None
    henry_used: float | None = None
    ssl_mg_kg: float | None = None
    notes: list[str] = []


class SslLeachingOutput(ModelOutput):
    method: str
    daf: float
    mixing_zone_depth_m: float | None = None
    chemicals: list[SslResult]


class SslLeachingModule(Module):
    name = "vadose.ssl-leaching"
    title = "Vadose zone: soil screening level protective of groundwater (EPA SSL)"
    input_type = SslLeachingInput
    output_type = SslLeachingOutput

    def run(self, inputs: SslLeachingInput) -> SslLeachingOutput:
        d = defaults()
        tw = inputs.theta_w or d["theta_w"]
        rb = inputs.rho_b_kg_l or d["rho_b_kg_l"]
        rs = inputs.rho_s_kg_l or d["rho_s_kg_l"]
        foc = inputs.foc if inputs.foc is not None else d["foc"]
        inf = inputs.infiltration_m_yr or d["infiltration_m_yr"]
        ed = inputs.ed_yr or d["ed_yr"]
        n = 1.0 - rb / rs                     # RSL UG 4.8.3, p. 190
        ta = n - tw
        warnings: list[str] = []
        if inputs.method == "mass-limit" and not inputs.source_depth_m:
            raise ValueError("The mass-limit method needs source_depth_m (RSL UG 4.8.4).")

        mz = None
        if inputs.site:
            s = inputs.site
            mz = eq.mixing_zone_depth(s.source_length_m, inf, s.hydraulic_conductivity_m_yr, s.hydraulic_gradient,
                                      s.aquifer_thickness_m)
            daf = eq.dilution_attenuation_factor(s.hydraulic_conductivity_m_yr, s.hydraulic_gradient, mz, inf,
                                                 s.source_length_m)
            if mz > s.aquifer_thickness_m:
                warnings.append(f"Mixing zone depth {mz:.3g} m exceeds the aquifer thickness "
                                f"{s.aquifer_thickness_m:.3g} m; check site inputs.")
        else:
            daf = inputs.daf or d["daf"]

        results = []
        for ci in inputs.chemicals:
            given = {k: getattr(ci, k) for k in ("kd_l_kg", "koc_l_kg", "henry_dimensionless")
                     if getattr(ci, k) is not None}
            if inputs.toxicity_table == "none":
                c = Chemical(name=ci.name or "", cas=ci.cas or "", **given)
            else:
                c = Chemical(**{**lookup(inputs.toxicity_table, ci.cas, ci.name).__dict__, **given})
            label = f"{c.name} ({c.cas})" if c.cas else c.name
            r = SslResult(name=c.name, cas=c.cas, target_basis="user" if ci.target_ug_l else "mcl")
            if ci.target == "mcl":
                r.target_ug_l = c.mcl_ug_l
                if c.mcl_ug_l is None:
                    warnings.append(f"{label}: no MCL in the RSL table; no SSL calculated.")
                    results.append(r)
                    continue
                if c.mcl_text and "(" in c.mcl_text:
                    r.notes.append(f"MCL cell printed as '{c.mcl_text}'; numeric value {c.mcl_ug_l:g} ug/L used.")
            else:
                r.target_ug_l = ci.target_ug_l
            inorganic = ci.inorganic if ci.inorganic is not None else (c.log_kow is None)
            kd = c.kd_l_kg if c.kd_l_kg is not None else (c.koc_l_kg * foc if c.koc_l_kg is not None else None)
            h = c.henry_dimensionless or 0.0
            if inorganic and c.cas != MERCURY_CAS and "henry_dimensionless" not in given:
                h = 0.0                       # SSG Eq. 4-10 note: H' = 0 for inorganics except mercury
            r.kd_l_kg, r.henry_used = kd, h
            r.cw_mg_l = r.target_ug_l / 1000.0 * daf
            if inputs.method == "partitioning":
                if kd is None:
                    warnings.append(f"{label}: no Kd or Koc; no SSL calculated.")
                    results.append(r)
                    continue
                r.ssl_mg_kg = eq.ssl_partitioning(r.cw_mg_l, kd, tw, ta, h, rb)
            else:
                r.ssl_mg_kg = eq.ssl_mass_limit(r.cw_mg_l, inf, ed, rb, inputs.source_depth_m)
            results.append(r)
        return SslLeachingOutput(method=inputs.method, daf=daf, mixing_zone_depth_m=mz, chemicals=results,
                                 warnings=warnings)

    def summary(self, o: SslLeachingOutput) -> str:
        lines = [f"Soil screening level protective of groundwater ({o.method}; DAF = {o.daf:.4g})",
                 f"{'Chemical':32s} {'CAS':12s} {'target ug/L':>12s} {'basis':>6s} {'Kd L/kg':>9s} {'SSL mg/kg':>10s}"]
        for c in o.chemicals:
            lines.append(f"{c.name[:32]:32s} {c.cas:12s} {(f'{c.target_ug_l:.3g}' if c.target_ug_l else '-'):>12s} "
                         f"{c.target_basis:>6s} {(f'{c.kd_l_kg:.3g}' if c.kd_l_kg is not None else '-'):>9s} "
                         f"{(f'{c.ssl_mg_kg:.2g}' if c.ssl_mg_kg is not None else '-'):>10s}")
        if o.mixing_zone_depth_m is not None:
            lines.append(f"Mixing zone depth : {o.mixing_zone_depth_m:.3g} m")
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
