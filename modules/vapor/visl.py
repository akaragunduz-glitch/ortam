"""Vapor intrusion screening levels (VISL): sub-slab/soil gas and groundwater targets.

Input per chemical: a target indoor air concentration (ug/m3), e.g. the risk-based
air SL from risk.resident-air or risk.worker-air (connected in
workflows/vapor_intrusion). The module converts it with the generic attenuation
factors and checks whether the chemical can reach the target (NVT rules).

Defaults (OSWER VI Technical Guide 2015, Table 6-1, PDF p. 132): sub-slab and
near-source exterior soil gas AF = 0.03; groundwater AF = 0.001 (0.0005 for
laterally extensive fine-grained layers). The VISL User's Guide (2014) still
lists 0.1 for soil gas; the 2015 guide supersedes it.
Sub-slab and groundwater targets are given for volatile chemicals only.
Equations: modules/vapor/equations.py.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from modules.vapor import equations as eq
from ortam_core.chemicals.database import Chemical, lookup
from ortam_core.model import ModelInput, ModelOutput, Module

_PROPS = ("henry_dimensionless", "henry_atm_m3_mol", "vp_mmhg", "mw_g_mol", "solubility_mg_l")


class ChemicalInput(ModelInput):
    cas: str | None = None
    name: str | None = None
    target_ia_ug_m3: float = Field(gt=0, description="Target indoor air concentration, ug/m3")
    henry_dimensionless: float | None = Field(None, ge=0, description="H' at 25 C, -")
    henry_atm_m3_mol: float | None = Field(None, ge=0, description="HLC at 25 C, atm-m3/mol")
    vp_mmhg: float | None = Field(None, gt=0, description="Vapor pressure at 25 C, mmHg")
    mw_g_mol: float | None = Field(None, gt=0, description="Molecular weight, g/mol")
    solubility_mg_l: float | None = Field(None, gt=0, description="Water solubility, mg/L")
    boiling_point_k: float | None = Field(None, gt=0, description="Normal boiling point Tb, K")
    critical_temp_k: float | None = Field(None, gt=0, description="Critical temperature Tc, K")
    enthalpy_vap_cal_mol: float | None = Field(None, gt=0, description="Enthalpy of vaporization at Tb, cal/mol")
    volatile: bool | None = Field(None, description="Volatile (default: RSL flag)")


class ViSlInput(ModelInput):
    af_subslab: float = Field(0.03, gt=0, le=1, description="Sub-slab / soil gas attenuation factor, -")
    af_groundwater: float = Field(0.001, gt=0, le=1, description="Groundwater attenuation factor, -")
    groundwater_temp_c: float = Field(25.0, ge=0, le=40, description="Average groundwater temperature, C")
    toxicity_table: Literal["rsl-2024-11", "none"] = "rsl-2024-11"
    chemicals: list[ChemicalInput] = Field(min_length=1)


class ViSlResult(BaseModel):
    name: str
    cas: str
    volatile: bool
    target_ia_ug_m3: float
    target_subslab_ug_m3: float | None = None
    target_groundwater_ug_l: float | None = None
    henry_used: float | None = None
    henry_temp_c: float | None = None
    pure_vapor_ug_m3: float | None = None
    nvt_indoor_air: bool = False
    nvt_subslab: bool = False
    nvt_groundwater: bool = False
    notes: list[str] = []


class ViSlOutput(ModelOutput):
    af_subslab: float
    af_groundwater: float
    groundwater_temp_c: float
    chemicals: list[ViSlResult]


class ViSlModule(Module):
    name = "vapor.visl"
    title = "Vapor intrusion: screening levels for sub-slab/soil gas and groundwater (VISL)"
    input_type = ViSlInput
    output_type = ViSlOutput

    def run(self, inputs: ViSlInput) -> ViSlOutput:
        warnings: list[str] = []
        if inputs.af_subslab != 0.03 or inputs.af_groundwater != 0.001:
            warnings.append(f"Attenuation factors differ from OSWER (2015) defaults: AF_ss = {inputs.af_subslab:g}, "
                            f"AF_gw = {inputs.af_groundwater:g}.")
        results = []
        for ci in inputs.chemicals:
            given = {k: getattr(ci, k) for k in _PROPS if getattr(ci, k) is not None}
            if inputs.toxicity_table == "none":
                c = Chemical(name=ci.name or "", cas=ci.cas or "", **given)
            else:
                c = Chemical(**{**lookup(inputs.toxicity_table, ci.cas, ci.name).__dict__, **given})
            label = f"{c.name} ({c.cas})" if c.cas else c.name
            volatile = ci.volatile if ci.volatile is not None else c.volatile
            r = ViSlResult(name=c.name, cas=c.cas, volatile=volatile, target_ia_ug_m3=ci.target_ia_ug_m3)
            if not volatile:
                r.notes.append("Not volatile: indoor air target only (no sub-slab or groundwater target).")
                results.append(r)
                continue

            # Henry's constant at the groundwater temperature
            h25 = c.henry_dimensionless
            if h25 is None and c.henry_atm_m3_mol is not None:
                h25 = eq.henry_dimensionless_25c(c.henry_atm_m3_mol)
            h, t_used = h25, 25.0
            if inputs.groundwater_temp_c != 25.0:
                props = (c.henry_atm_m3_mol, ci.boiling_point_k, ci.critical_temp_k, ci.enthalpy_vap_cal_mol)
                if None in props:
                    warnings.append(f"{label}: Tb, Tc or dHv,b not available; H' at 25 C used (VISL UG p. 9).")
                else:
                    h = eq.henry_at_temperature(c.henry_atm_m3_mol, inputs.groundwater_temp_c, *props[1:])
                    t_used = inputs.groundwater_temp_c
            r.henry_used, r.henry_temp_c = h, t_used

            # NVT checks (VISL UG pp. 7-9)
            if c.vp_mmhg is not None and c.mw_g_mol is not None:
                r.pure_vapor_ug_m3 = eq.pure_vapor_concentration(c.vp_mmhg, c.mw_g_mol)
                if ci.target_ia_ug_m3 > r.pure_vapor_ug_m3:
                    r.nvt_indoor_air = True
                    r.notes.append("NVT: target indoor air exceeds the pure-component vapor concentration.")
                    results.append(r)
                    continue
            sg = eq.subslab_target(ci.target_ia_ug_m3, inputs.af_subslab)
            if r.pure_vapor_ug_m3 is not None and sg > r.pure_vapor_ug_m3:
                r.nvt_subslab = True
                r.notes.append("NVT: soil gas target exceeds the pure-component vapor concentration.")
            else:
                r.target_subslab_ug_m3 = sg
            if h:
                gw = eq.groundwater_target(ci.target_ia_ug_m3, h, inputs.af_groundwater)
                if c.solubility_mg_l is not None and gw > c.solubility_mg_l * 1000.0:
                    r.nvt_groundwater = True
                    r.notes.append("NVT: groundwater target exceeds the water solubility.")
                else:
                    r.target_groundwater_ug_l = gw
            else:
                r.notes.append("No Henry's law constant: no groundwater target.")
            results.append(r)
        return ViSlOutput(af_subslab=inputs.af_subslab, af_groundwater=inputs.af_groundwater,
                          groundwater_temp_c=inputs.groundwater_temp_c, chemicals=results, warnings=warnings)

    def summary(self, o: ViSlOutput) -> str:
        lines = [f"Vapor intrusion screening levels (AF_ss = {o.af_subslab:g}, AF_gw = {o.af_groundwater:g}, "
                 f"T_gw = {o.groundwater_temp_c:g} C)",
                 f"{'Chemical':30s} {'Indoor air ug/m3':>16s} {'Sub-slab ug/m3':>15s} {'Groundwater ug/L':>17s}"]
        f = lambda v, nvt: "NVT" if nvt else (f"{v:.2g}" if v is not None else "-")
        for c in o.chemicals:
            lines.append(f"{c.name[:30]:30s} {f(c.target_ia_ug_m3, c.nvt_indoor_air):>16s} "
                         f"{f(c.target_subslab_ug_m3, c.nvt_subslab or c.nvt_indoor_air):>15s} "
                         f"{f(c.target_groundwater_ug_l, c.nvt_groundwater or c.nvt_indoor_air):>17s}")
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
