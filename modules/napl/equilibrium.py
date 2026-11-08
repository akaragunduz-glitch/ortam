"""napl.equilibrium: NAPL/water partitioning (Raoult's law), effective solubility, NAPL saturation.

Two bases:
  basis "napl": the NAPL composition is known (coj in g/cm3, as in the HSSM RAOULT utility, or
      mass fractions with the NAPL density). Results: density, average molecular weight, mole
      fractions, effective solubilities, Ko (HSSM1 Eq. 44 and 45).
  basis "soil": total soil concentrations (mg/kg) of the NAPL components. Results: whether a
      NAPL phase is present (sum Ct/Csat > 1), its composition, the NAPL volume fraction and
      saturation, and the concentrations in water, soil and soil gas (derivation:
      docs/methods/napl.md; equations: modules/napl/equations.py).
Groundwater and soil data are screened against the LNAPL indicators of EPA 510-R-15-001
(data/reference/napl/lnapl_indicators.yaml).
Chemical properties: RSL Chemical Parameters table (lookup by CAS) unless given. A component
without a CAS number (e.g. "other paraffins") needs its own molecular weight and solubility;
solubility 0 marks an insoluble component (as in the HSSM RAOULT.DAT data set).
Ideal solution (activity coefficient 1) unless given. Air-phase partitioning is included in
the soil basis (SSG Eq. 4-9); HSSM neglects it (HSSM2 p. 16 [30]).
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, Field, model_validator

from modules.napl import equations as eq
from ortam_core.chemicals.database import lookup
from ortam_core.model import ModelInput, ModelOutput, Module

ROOT = Path(__file__).resolve().parents[2] / "data" / "reference"
BENZENE = "71-43-2"


@lru_cache
def _indicators() -> dict:
    return yaml.safe_load((ROOT / "napl" / "lnapl_indicators.yaml").read_text(encoding="utf-8"))


@lru_cache
def _soil_defaults() -> dict:
    return yaml.safe_load((ROOT / "rsl" / "2024-11" / "air_emission_factors.yaml").read_text(encoding="utf-8"))["soil"]


class Component(ModelInput):
    name: str | None = None
    cas: str | None = None
    mw_g_mol: float | None = Field(None, gt=0, description="Molecular weight, g/mol")
    solubility_mg_l: float | None = Field(None, ge=0, description="Pure-phase water solubility, mg/L (0 = insoluble)")
    activity_coefficient: float = Field(1.0, gt=0, description="Activity coefficient in the NAPL, - (1 = ideal)")
    koc_l_kg: float | None = Field(None, ge=0, description="Koc, L/kg")
    kd_l_kg: float | None = Field(None, ge=0, description="Kd, L/kg (overrides Koc*foc)")
    henry_dimensionless: float | None = Field(None, ge=0, description="H', -")
    mp_c: float | None = Field(None, description="Melting point, C")
    napl_g_cm3: float | None = Field(None, ge=0, description="Concentration in the NAPL coj, g/cm3 (basis napl)")
    napl_mass_fraction: float | None = Field(None, ge=0, le=1, description="Mass fraction in the NAPL, - (basis napl)")
    soil_mg_kg: float | None = Field(None, gt=0, description="Total soil concentration, mg/kg dry soil (basis soil). "
                                                         "Non-detects are not entered as 0 or DL/2: leave empty.")
    groundwater_mg_l: float | None = Field(None, gt=0, description="Measured groundwater concentration, mg/L. "
                                                               "Non-detects are not entered.")


class Soil(ModelInput):
    rho_b_g_cm3: float | None = Field(None, gt=0, description="Dry bulk density, g/cm3 (default 1.5)")
    rho_s_g_cm3: float | None = Field(None, gt=0, description="Particle density, g/cm3 (default 2.65)")
    theta_w: float | None = Field(None, ge=0, lt=1, description="Water-filled porosity, - (default 0.15)")
    foc: float | None = Field(None, ge=0, lt=1, description="Organic carbon fraction, - (default 0.006)")


class NaplInput(ModelInput):
    basis: Literal["napl", "soil"]
    chemical_table: Literal["rsl-2024-11", "none"] = "rsl-2024-11"
    napl_density_g_cm3: float | None = Field(None, gt=0, description="NAPL density, g/cm3 (default for coj input: sum coj)")
    soil: Soil = Soil()
    residual_saturation: float | None = Field(None, gt=0, lt=1, description="Residual NAPL saturation Sor, - (site data)")
    fuel: Literal["unweathered-gasoline", "weathered-gasoline", "diesel"] | None = None
    tph_soil_mg_kg: float | None = Field(None, gt=0, description="Measured TPH in soil, mg/kg")
    tph_groundwater_mg_l: float | None = Field(None, gt=0, description="Measured TPH in groundwater, mg/L")
    components: list[Component] = Field(min_length=1)

    @model_validator(mode="after")
    def _check(self):
        if self.basis == "napl":
            for c in self.components:
                if (c.napl_g_cm3 is None) == (c.napl_mass_fraction is None):
                    raise ValueError(f"basis napl: give exactly one of napl_g_cm3 or napl_mass_fraction "
                                     f"for {c.name or c.cas}.")
            kinds = {c.napl_g_cm3 is None for c in self.components}
            if len(kinds) > 1:
                raise ValueError("basis napl: use napl_g_cm3 for all components or napl_mass_fraction for all.")
            if self.components[0].napl_g_cm3 is None and self.napl_density_g_cm3 is None:
                raise ValueError("basis napl with mass fractions needs napl_density_g_cm3.")
        elif not any(c.soil_mg_kg for c in self.components):
            raise ValueError("basis soil: at least one component needs soil_mg_kg.")
        return self


class ComponentResult(BaseModel):
    name: str
    cas: str | None = None
    mw_g_mol: float
    solubility_mg_l: float | None
    mole_fraction: float | None = None
    mass_fraction: float | None = None
    effective_solubility_mg_l: float | None = Field(None, description="x * gamma * S (BIOSCREEN p. 3)")
    ko_composition: float | None = Field(None, description="HSSM1 Eq. 44")
    ko_average_mw: float | None = Field(None, description="HSSM1 Eq. 45")
    csat_mg_kg: float | None = Field(None, description="Single-chemical Csat, SSG Eq. 4-9")
    soil_mg_kg: float | None = None
    water_mg_l: float | None = None
    sorbed_mg_kg: float | None = None
    soil_gas_mg_l: float | None = None
    napl_mg_l: float | None = Field(None, description="Concentration in the NAPL co, mg/L of NAPL")
    in_napl_pct: float | None = Field(None, description="Share of the total mass in the NAPL, %")
    groundwater_mg_l: float | None = None
    groundwater_to_effective_solubility: float | None = None
    notes: list[str] = []


class Indicator(BaseModel):
    indicator: str
    value: float
    threshold: float
    unit: str
    exceeded: bool
    source: str


class NaplOutput(ModelOutput):
    basis: str
    napl_present: bool | None = None
    saturation_index: float | None = Field(None, description="sum Ct/Csat (basis soil); > 1: NAPL present")
    napl_density_g_cm3: float | None = None
    average_mw_g_mol: float | None = Field(None, description="Mass-weighted, as in HSSM RAOULT")
    total_effective_solubility_mg_l: float | None = None
    porosity: float | None = None
    napl_volume_fraction: float | None = Field(None, description="theta_o, L NAPL per L soil")
    napl_saturation: float | None = Field(None, description="So = theta_o/n, -")
    napl_mg_kg: float | None = Field(None, description="NAPL mass per dry soil mass, mg/kg")
    potentially_mobile: bool | None = Field(None, description="So > residual saturation (if given)")
    soil: dict = {}
    components: list[ComponentResult]
    indicators: list[Indicator] = []


def _props(c: Component, table: str) -> dict:
    """Chemical properties: user values override the RSL table (lookup by CAS)."""
    p = {"name": c.name, "cas": c.cas}
    if c.cas and table != "none":
        r = lookup(table, cas=c.cas, name=c.name)
        p.update(name=c.name or r.name, mw=r.mw_g_mol, s=r.solubility_mg_l, koc=r.koc_l_kg, kd=r.kd_l_kg,
                 h=r.henry_dimensionless, mp=r.mp_c)
    for key, val in (("mw", c.mw_g_mol), ("s", c.solubility_mg_l), ("koc", c.koc_l_kg), ("kd", c.kd_l_kg),
                     ("h", c.henry_dimensionless), ("mp", c.mp_c)):
        if val is not None:
            p[key] = val
    p["name"] = p.get("name") or c.cas or "component"
    if not p.get("mw"):
        raise ValueError(f"{p['name']}: molecular weight is required (no CAS number or not in the table).")
    if p.get("s") is None:
        raise ValueError(f"{p['name']}: water solubility is required (0 for an insoluble component).")
    return p


class NaplEquilibriumModule(Module):
    name = "napl.equilibrium"
    title = "NAPL: Raoult's law partitioning, effective solubility and NAPL saturation"
    input_type = NaplInput
    output_type = NaplOutput

    def run(self, inp: NaplInput) -> NaplOutput:
        props = [_props(c, inp.chemical_table) for c in inp.components]
        warnings: list[str] = []
        results = [ComponentResult(name=p["name"], cas=p.get("cas"), mw_g_mol=p["mw"], solubility_mg_l=p["s"])
                   for p in props]
        for p, r in zip(props, results):
            if p.get("mp") is not None and p["mp"] >= 20.0:
                r.notes.append("Solid at soil temperature (MP >= 20 C, RSL UG 5.12): the pure-phase solubility "
                               "of the solid is used in Raoult's law.")
        out = NaplOutput(basis=inp.basis, components=results, warnings=warnings)
        if inp.basis == "napl":
            self._napl_basis(inp, props, out)
        else:
            self._soil_basis(inp, props, out)
        out.indicators = self._indicators(inp, out)
        return out

    # ------------------------------------------------------------ basis napl
    def _napl_basis(self, inp: NaplInput, props: list[dict], out: NaplOutput) -> None:
        cs = inp.components
        mw = [p["mw"] for p in props]
        if cs[0].napl_g_cm3 is not None:
            co = [c.napl_g_cm3 for c in cs]
            rho = inp.napl_density_g_cm3 or eq.napl_density(co)
            if inp.napl_density_g_cm3 and abs(inp.napl_density_g_cm3 / eq.napl_density(co) - 1) > 0.01:
                out.warnings.append("napl_density_g_cm3 differs from sum coj by more than 1 %.")
        else:
            rho = inp.napl_density_g_cm3
            wsum = sum(c.napl_mass_fraction for c in cs)
            if abs(wsum - 1.0) > 0.01:
                out.warnings.append(f"Mass fractions sum to {wsum:.3f}: mole fractions and Ko are computed over "
                                    "the listed components only; add pseudo-components for the rest.")
            co = [c.napl_mass_fraction * rho for c in cs]
        out.napl_density_g_cm3 = rho
        out.average_mw_g_mol = eq.average_mw(co, mw)
        x = eq.mole_fractions(co, mw)
        tot = sum(co)
        for p, r, c, xi, comp in zip(props, out.components, co, x, cs):
            r.mole_fraction, r.mass_fraction = xi, c / tot
            r.napl_mg_l = c * eq.MG_PER_G * eq.CM3_PER_L
            g = comp.activity_coefficient
            if p["s"] > 0:
                r.effective_solubility_mg_l = eq.effective_solubility(xi, p["s"], g)
                r.ko_composition = eq.ko_composition(p["mw"], p["s"], g, co, mw)
                r.ko_average_mw = eq.ko_average_mw(rho, out.average_mw_g_mol, p["s"], p["mw"])
        out.total_effective_solubility_mg_l = sum(r.effective_solubility_mg_l or 0.0 for r in out.components)

    # ------------------------------------------------------------ basis soil
    def _soil_basis(self, inp: NaplInput, props: list[dict], out: NaplOutput) -> None:
        d = _soil_defaults()
        rb = inp.soil.rho_b_g_cm3 or d["rho_b_g_cm3"]["value"]
        rs = inp.soil.rho_s_g_cm3 or d["rho_s_g_cm3"]["value"]
        tw = inp.soil.theta_w if inp.soil.theta_w is not None else d["theta_w"]["value"]
        foc = inp.soil.foc if inp.soil.foc is not None else d["foc"]["value"]
        n, ta = eq.porosities(rb, rs, tw)
        if ta <= 0:
            raise ValueError("theta_w must be smaller than the total porosity.")
        out.porosity = n
        out.soil = {"rho_b_g_cm3": rb, "rho_s_g_cm3": rs, "theta_w": tw, "theta_a": ta, "foc": foc}
        idx = [i for i, c in enumerate(inp.components) if c.soil_mg_kg is not None]
        for i, c in enumerate(inp.components):
            if c.soil_mg_kg is None:
                out.components[i].notes.append("No soil concentration (not detected or not measured): not in the "
                                               "mass balance.")
        nk, ak, bp, kd = [], [], {}, {}
        si = 0.0
        for i in idx:
            p, c, r = props[i], inp.components[i], out.components[i]
            kd[i] = p["kd"] if p.get("kd") is not None else (p["koc"] * foc if p.get("koc") is not None else 0.0)
            if p.get("kd") is None and p.get("koc") is None:
                r.notes.append("No Koc or Kd: sorption neglected (Kd = 0).")
            h = p.get("h") or 0.0
            if p.get("h") is None:
                r.notes.append("No Henry's law constant: soil-gas partitioning neglected.")
            bp[i] = eq.bulk_partition(tw, rb, kd[i], ta, h)
            r.soil_mg_kg = c.soil_mg_kg
            nk.append(c.soil_mg_kg * rb / p["mw"])                               # mmol/L soil
            ak.append(c.activity_coefficient * p["s"] * bp[i] / p["mw"])          # mmol/L soil per unit x
            if p["s"] > 0:
                r.csat_mg_kg = eq.csat(p["s"], rb, bp[i])
                si += c.soil_mg_kg / r.csat_mg_kg * (1.0 / c.activity_coefficient)
        has_insoluble = any(props[i]["s"] == 0 for i in idx)
        out.saturation_index = None if has_insoluble else si
        no, x = eq.solve_napl_moles(nk, ak)
        out.napl_present = no > 0.0
        if not out.napl_present:
            for i in idx:
                r, p = out.components[i], props[i]
                r.water_mg_l = inp.components[i].soil_mg_kg * rb / bp[i]
                r.in_napl_pct = 0.0
                self._phases(r, p, kd[i])
            return
        mass_o = no * sum(xi * props[i]["mw"] for xi, i in zip(x, idx))          # mg NAPL per L soil
        out.napl_mg_kg = mass_o / rb
        out.average_mw_g_mol = sum(xi * props[i]["mw"] ** 2 for xi, i in zip(x, idx)) / \
            sum(xi * props[i]["mw"] for xi, i in zip(x, idx))                      # mass-weighted (RAOULT)
        for xi, i, n_i in zip(x, idx, nk):
            r, p, c = out.components[i], props[i], inp.components[i]
            r.mole_fraction = xi
            r.mass_fraction = xi * no * p["mw"] / mass_o
            r.in_napl_pct = 100.0 * xi * no / n_i
            r.water_mg_l = eq.effective_solubility(xi, p["s"], c.activity_coefficient)
            if p["s"] > 0:
                r.effective_solubility_mg_l = r.water_mg_l
            self._phases(r, p, kd[i])
        out.total_effective_solubility_mg_l = sum(r.effective_solubility_mg_l or 0.0 for r in out.components)
        rho = inp.napl_density_g_cm3
        if rho is None:
            out.warnings.append("napl_density_g_cm3 not given: NAPL volume fraction, saturation, Ko and the "
                                "concentration in the NAPL are not computed.")
            return
        out.napl_density_g_cm3 = rho
        out.napl_volume_fraction = mass_o / (rho * eq.MG_PER_G * eq.CM3_PER_L)
        out.napl_saturation = out.napl_volume_fraction / n
        if out.napl_volume_fraction > 0.1 * ta:
            out.warnings.append("NAPL volume exceeds 10 % of the air-filled porosity: the water and air contents "
                                "are not reduced by the NAPL in this calculation.")
        co = [r.mass_fraction * rho for r in (out.components[i] for i in idx)]     # g/cm3
        mw = [props[i]["mw"] for i in idx]
        for i, c_o in zip(idx, co):
            r, p = out.components[i], props[i]
            r.napl_mg_l = c_o * eq.MG_PER_G * eq.CM3_PER_L
            if p["s"] > 0:
                g = inp.components[i].activity_coefficient
                r.ko_composition = eq.ko_composition(p["mw"], p["s"], g, co, mw)
                r.ko_average_mw = eq.ko_average_mw(rho, out.average_mw_g_mol, p["s"], p["mw"])
        if inp.residual_saturation is not None:
            out.potentially_mobile = out.napl_saturation > inp.residual_saturation

    @staticmethod
    def _phases(r: ComponentResult, p: dict, kd: float) -> None:
        r.sorbed_mg_kg = kd * r.water_mg_l                    # Kd = cs/cw, HSSM2 Eq. 26
        r.soil_gas_mg_l = (p.get("h") or 0.0) * r.water_mg_l    # H' = cg/cw (SSG Eq. 4-9 term)

    # ------------------------------------------------------------ indicators
    def _indicators(self, inp: NaplInput, out: NaplOutput) -> list[Indicator]:
        d = _indicators()
        res = []

        def add(what, value, thr, unit, src):
            res.append(Indicator(indicator=what, value=value, threshold=thr, unit=unit, exceeded=value > thr,
                                 source=f"{d['source']}; {src}"))

        for c, r in zip(inp.components, out.components):
            if c.cas == BENZENE and c.groundwater_mg_l is not None:
                g = d["groundwater_mg_l"]["benzene"]
                add("benzene in groundwater", c.groundwater_mg_l, g["value"], "mg/L", g["source"])
            if c.cas == BENZENE and c.soil_mg_kg is not None:
                g = d["soil_mg_kg"]["benzene"]
                add("benzene in soil", c.soil_mg_kg, g["value"], "mg/kg", g["source"])
            if c.groundwater_mg_l is not None and r.effective_solubility_mg_l:
                r.groundwater_mg_l = c.groundwater_mg_l
                r.groundwater_to_effective_solubility = c.groundwater_mg_l / r.effective_solubility_mg_l
                f = d["effective_solubility_fraction"]
                add(f"{r.name} in groundwater / effective solubility", r.groundwater_to_effective_solubility,
                    f["value"], "-", f["source"])
            elif c.groundwater_mg_l is not None:
                r.groundwater_mg_l = c.groundwater_mg_l
        if inp.tph_groundwater_mg_l is not None:
            g = d["groundwater_mg_l"]["tph"]
            add("TPH in groundwater", inp.tph_groundwater_mg_l, g["value"], "mg/L", g["source"])
        if inp.tph_soil_mg_kg is not None:
            if inp.fuel is None:
                out.warnings.append("tph_soil_mg_kg given without fuel: TPH soil indicator not evaluated.")
            else:
                g = d["soil_mg_kg"]["tph"][inp.fuel]
                add(f"TPH in soil ({inp.fuel})", inp.tph_soil_mg_kg, g["value"], "mg/kg", g["source"])
        return res

    def summary(self, o: NaplOutput) -> str:
        lines = [f"NAPL equilibrium (basis {o.basis}); Raoult's law (HSSM Eq. 44-45)"]
        if o.basis == "soil":
            lines.append(f"NAPL present      : {'yes' if o.napl_present else 'no'}"
                         + (f"   (sum Ct/Csat = {o.saturation_index:.3g})" if o.saturation_index is not None else ""))
            if o.napl_saturation is not None:
                lines.append(f"NAPL saturation So: {o.napl_saturation:.3g}  (theta_o = {o.napl_volume_fraction:.3g})")
            if o.napl_mg_kg:
                lines.append(f"NAPL in soil      : {o.napl_mg_kg:.4g} mg/kg")
        if o.napl_density_g_cm3:
            lines.append(f"NAPL density      : {o.napl_density_g_cm3:.4g} g/cm3")
        if o.average_mw_g_mol:
            lines.append(f"Avg. mol. weight  : {o.average_mw_g_mol:.4g} g/mol")
        lines.append(f"{'Component':22s} {'x':>9s} {'Seff mg/L':>10s} {'Ko(44)':>10s} {'Ko(45)':>10s} "
                     f"{'Cw mg/L':>10s} {'Csat mg/kg':>10s}")

        def f(v, fmt):
            return format(v, fmt) if v is not None else "-"
        for c in o.components:
            lines.append(f"{c.name[:22]:22s} {f(c.mole_fraction, '9.4f'):>9s} {f(c.effective_solubility_mg_l, '10.4g'):>10s} "
                         f"{f(c.ko_composition, '10.4g'):>10s} {f(c.ko_average_mw, '10.4g'):>10s} "
                         f"{f(c.water_mg_l, '10.4g'):>10s} {f(c.csat_mg_kg, '10.4g'):>10s}")
        if o.total_effective_solubility_mg_l:
            lines.append(f"Total effective solubility: {o.total_effective_solubility_mg_l:.4g} mg/L")
        for i in o.indicators:
            lines.append(f"LNAPL indicator: {i.indicator} = {i.value:.3g} {i.unit} (threshold {i.threshold:g}) "
                         f"-> {'EXCEEDED' if i.exceeded else 'not exceeded'}")
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
