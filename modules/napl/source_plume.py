"""napl.source-plume: dissolved plume from an LNAPL source zone (BIOSCREEN 1.4).

Domenico (1987) analytical model with three source layers ("layer cake"), declining source
(first-order source mass depletion), and three plume models: no degradation, first-order decay,
instantaneous reaction with the electron acceptors O2, NO3, Fe2+, SO4 and CH4 (superposition).
Outputs: centreline and 11 x 5 array concentrations, BIOSCREEN's plume and source mass balance,
mass flux across cross-sections (BIOSCREEN 1.4), concentrations at receptors (optionally as a
time series), and a plume mass and mass flux from a fine-grid numerical integration.

Input and output are SI (m, yr, mg/L, kg, m3); the engine (modules/napl/bioscreen.py) works in
BIOSCREEN's units (ft, yr) and the conversion is exact (1 ft = 0.3048 m, 1 yr = 365 d).
bioscreen_compatible (default true) reproduces BIOSCREEN 1.4 exactly, including its constants
and three workbook inconsistencies; false gives the corrected, consistent form
(docs/methods/napl.md §4).
Sources: BIOSCREEN User's Manual 1.3 (EPA/600/R-96/087) with the Version 1.4 revisions;
BIOSCRN4.xls. instantaneous_source (default "mass-balance") replaces BIOSCREEN's source depletion
of the instantaneous model by a mass-balanced one (docs/methods/napl.md §4.8); "bioscreen" keeps it.
Verification: modules/napl/tests/test_bioscreen.py, test_source_plume.py.
"""
from __future__ import annotations

import math
from typing import Literal

import numpy as np
from pydantic import BaseModel, Field, model_validator

from modules.napl import bioscreen as bs
from ortam_core.chemicals.database import lookup
from ortam_core.model import ModelInput, ModelOutput, Module

FT = 0.3048                 # m per ft (exact)
FT3 = FT ** 3               # m3 per ft3
MODEL_NAMES = {"no_decay": "No degradation", "first_order": "First-order decay",
               "instantaneous": "Instantaneous reaction"}


class Hydrogeology(ModelInput):
    seepage_velocity_m_yr: float | None = Field(None, gt=0, description="Seepage velocity Vs, m/yr")
    hydraulic_conductivity_m_s: float | None = Field(None, gt=0, description="Hydraulic conductivity K, m/s")
    hydraulic_gradient: float | None = Field(None, gt=0, description="Hydraulic gradient i, m/m")
    effective_porosity: float = Field(gt=0, lt=1, description="Effective porosity n, -")

    @model_validator(mode="after")
    def _check(self):
        k = self.hydraulic_conductivity_m_s is not None and self.hydraulic_gradient is not None
        if (self.seepage_velocity_m_yr is None) == (not k):
            raise ValueError("Give either seepage_velocity_m_yr or hydraulic_conductivity_m_s with hydraulic_gradient.")
        return self


class Dispersion(ModelInput):
    alpha_x_m: float | None = Field(None, gt=0, description="Longitudinal dispersivity, m")
    plume_length_m: float | None = Field(None, gt=0, description="Estimated plume length Lp, m "
                                         "(alpha_x from Xu and Eckstein 1995 if alpha_x is not given)")
    alpha_y_m: float | None = Field(None, gt=0, description="Transverse dispersivity, m (default 0.1 alpha_x)")
    alpha_z_m: float = Field(0.0, ge=0, description="Vertical dispersivity, m (default 0: no vertical dispersion)")

    @model_validator(mode="after")
    def _check(self):
        if self.alpha_x_m is None and self.plume_length_m is None:
            raise ValueError("Give alpha_x_m or plume_length_m.")
        return self


class Adsorption(ModelInput):
    retardation_factor: float | None = Field(None, ge=1, description="Retardation factor R, -")
    bulk_density_kg_l: float | None = Field(None, gt=0, description="Aquifer bulk density, kg/L")
    koc_l_kg: float | None = Field(None, ge=0, description="Koc, L/kg (RSL table if a chemical is given)")
    foc: float | None = Field(None, ge=0, lt=1, description="Fraction organic carbon, -")


class Biodegradation(ModelInput):
    decay_rate_per_yr: float | None = Field(None, ge=0, description="First-order decay coefficient lambda, 1/yr")
    solute_half_life_yr: float | None = Field(None, gt=0, description="Dissolved plume half-life, yr")
    delta_oxygen_mg_l: float = Field(0.0, ge=0, description="Background minus source-zone minimum DO, mg/L")
    delta_nitrate_mg_l: float = Field(0.0, ge=0, description="Background minus source-zone minimum NO3, mg/L")
    ferrous_iron_mg_l: float = Field(0.0, ge=0, description="Average observed Fe2+ in the source area, mg/L")
    delta_sulfate_mg_l: float = Field(0.0, ge=0, description="Background minus source-zone minimum SO4, mg/L")
    methane_mg_l: float = Field(0.0, ge=0, description="Average observed CH4 in the source area, mg/L")
    utilization_factors: dict[Literal["oxygen", "nitrate", "ferrous_iron", "sulfate", "methane"], float] | None = \
        Field(None, description="Override of the BTEX utilization factors (mg/mg)")

    @model_validator(mode="after")
    def _check(self):
        if self.decay_rate_per_yr is not None and self.solute_half_life_yr is not None:
            raise ValueError("Give decay_rate_per_yr or solute_half_life_yr, not both.")
        return self


class SourceZone(ModelInput):
    width_m: float = Field(gt=0, description="Width of the zone, m (central zone: full width; outer zones: "
                                              "width of the band on each side)")
    concentration_mg_l: float = Field(gt=0, description="Average dissolved concentration in the zone, mg/L")


class Source(ModelInput):
    thickness_m: float = Field(gt=0, description="Source thickness in the saturated zone Z, m")
    zones: list[SourceZone] = Field(min_length=1, max_length=3,
                                    description="From the centre outward: central zone, then up to two side "
                                                "bands (BIOSCREEN zones 3, 2/4, 1/5)")
    soluble_mass_kg: float | None = Field(None, gt=0, description="Soluble (dissolvable) mass in NAPL and soil, kg; "
                                          "empty: infinite (constant) source")


class General(ModelInput):
    model_length_m: float = Field(gt=0, description="Modelled area length L, m")
    model_width_m: float = Field(gt=0, description="Modelled area width W, m")
    simulation_time_yr: float = Field(gt=0, description="Simulation time t, yr")


class Receptor(ModelInput):
    name: str
    x_m: float = Field(ge=0, description="Distance downgradient of the source, m")
    y_m: float = Field(0.0, description="Distance from the plume centreline, m")


class Chemical(ModelInput):
    cas: str | None = None
    name: str | None = None


class SourcePlumeInput(ModelInput):
    chemical: Chemical | None = Field(None, description="Constituent (or lumped BTEX); used for Koc lookup")
    chemical_table: Literal["rsl-2024-11", "none"] = "rsl-2024-11"
    hydrogeology: Hydrogeology
    dispersion: Dispersion
    adsorption: Adsorption = Adsorption(retardation_factor=1.0)
    biodegradation: Biodegradation = Biodegradation()
    general: General
    source: Source
    receptors: list[Receptor] = []
    receptor_times_yr: list[float] = Field([], description="Times for receptor time series, yr")
    volume_threshold_mg_l: float = Field(0.005, gt=0, description="Concentration defining the plume volume, mg/L")
    integration_points: int = Field(201, ge=21, le=2001, description="Grid points per direction for the "
                                    "numerical plume mass and mass flux")
    bioscreen_compatible: bool = True
    instantaneous_source: Literal["mass-balance", "bioscreen"] = Field(
        "mass-balance", description="Source depletion of the instantaneous model: 'mass-balance' (BC constant, "
        "source mass and concentration vanish together) or 'bioscreen' (BIOSCREEN 1.4)")


class Derived(BaseModel):
    seepage_velocity_m_yr: float
    retardation_factor: float
    contaminant_velocity_m_yr: float
    alpha_x_m: float
    alpha_y_m: float
    alpha_z_m: float
    decay_rate_per_yr: float
    biodegradation_capacity_mg_l: float
    source_width_m: float
    source_flow_m3_yr: float
    source_concentration_avg_mg_l: float
    source_flux_kg_yr: float
    source_decay_rate_per_yr: float
    source_half_life_yr: float | None = Field(None, description="No degradation and first-order models; None: infinite")
    source_concentration_avg_inst_mg_l: float
    source_flux_inst_kg_yr: float
    instantaneous_source: str
    source_decay_rate_inst_per_yr: float | None = Field(None, description="'bioscreen' source only")
    source_half_life_inst_yr: float | None = Field(None, description="'bioscreen' source only; None: infinite or n/a")
    source_depletion_time_inst_yr: float | None = Field(None, description="'mass-balance' source: time at which "
                                                        "source mass and concentration reach zero; None: never")
    layer_widths_m: list[float]
    layer_increments_mg_l: list[float]
    layer_increments_inst_mg_l: list[float]


class CenterlinePoint(BaseModel):
    x_m: float
    no_decay_mg_l: float
    first_order_mg_l: float
    instantaneous_mg_l: float


class PlumeArray(BaseModel):
    x_m: list[float]
    y_m: list[float]
    no_decay_mg_l: list[list[float]]
    first_order_mg_l: list[list[float]]
    instantaneous_mg_l: list[list[float]]


class MassBalance(BaseModel):
    model: str
    mass_left_source_kg: float | None = Field(None, description="'Plume mass if no biodegradation'")
    actual_plume_mass_kg: float | None = None
    mass_removed_kg: float | None = None
    fraction_removed: float | None = None
    electron_acceptor_change_kg: dict[str, float] | None = None
    source_mass_initial_kg: float | None = None
    source_mass_now_kg: float | None = None
    plume_volume_m3: float | None = None
    source_flow_m3_yr: float | None = None
    mass_flux_mg_day: list[float] | None = Field(None, description="At the 11 array x positions (BIOSCREEN 1.4)")
    cant_calc: bool = Field(False, description="Plume mass on the grid outside 50-150 % of the mass left in the source")


class Integrated(BaseModel):
    model: str
    plume_mass_kg: float = Field(description="n Z R integral of C over the modelled area")
    mass_flux_mg_day: list[float] = Field(description="q Z integral of C dy across the model width, at the array x")


class ReceptorResult(BaseModel):
    name: str
    x_m: float
    y_m: float
    no_decay_mg_l: float
    first_order_mg_l: float
    instantaneous_mg_l: float
    time_series: list[dict[str, float]] = []


class SourcePlumeOutput(ModelOutput):
    bioscreen_compatible: bool
    chemical: str | None = None
    derived: Derived
    centerline: list[CenterlinePoint]
    array: PlumeArray
    mass_balance: list[MassBalance]
    integrated: list[Integrated] | None = None
    receptors: list[ReceptorResult] = []
    notes: list[str] = []


def _hl(v: float) -> float | None:
    return None if math.isinf(v) else v


class SourcePlumeModule(Module):
    name = "napl.source-plume"
    title = "NAPL: dissolved plume from a source zone (BIOSCREEN 1.4, Domenico model)"
    input_type = SourcePlumeInput
    output_type = SourcePlumeOutput

    def case(self, inp: SourcePlumeInput) -> tuple[bs.Case, bs.Settings, list[str], str | None]:
        """SI input -> BIOSCREEN case in ft and yr (exact conversion)."""
        s = bs.Settings(compatible=inp.bioscreen_compatible)
        notes: list[str] = []
        h, d, a, b = inp.hydrogeology, inp.dispersion, inp.adsorption, inp.biodegradation
        n = h.effective_porosity
        if h.seepage_velocity_m_yr is not None:
            vs_ft = h.seepage_velocity_m_yr / FT
        else:   # K (m/s) -> cm/s for BSCR p. 16 [22] / XLS Input!E5
            vs_ft = bs.seepage_velocity_ft_yr(h.hydraulic_conductivity_m_s * 100.0, h.hydraulic_gradient, n)
        ax = d.alpha_x_m / FT if d.alpha_x_m is not None else bs.alpha_x_xu_eckstein_ft(d.plume_length_m / FT, s)
        ay = d.alpha_y_m / FT if d.alpha_y_m is not None else bs.alpha_y_default(ax)
        chem_name = None
        if a.retardation_factor is not None:
            r = a.retardation_factor
        else:
            koc = a.koc_l_kg
            if koc is None:
                if inp.chemical is None or inp.chemical_table == "none":
                    raise ValueError("Give retardation_factor, or koc_l_kg (or a chemical) with bulk density and foc.")
                ch = lookup(inp.chemical_table, inp.chemical.cas, inp.chemical.name)
                koc, chem_name = ch.koc_l_kg, ch.name
                if koc is None:
                    raise ValueError(f"No Koc in {inp.chemical_table} for {ch.name}; give koc_l_kg.")
                notes.append(f"Koc = {koc:g} L/kg from {inp.chemical_table} ({ch.name}).")
            if a.bulk_density_kg_l is None or a.foc is None:
                raise ValueError("Retardation from Koc needs bulk_density_kg_l and foc.")
            r = bs.retardation_factor(a.bulk_density_kg_l, koc, a.foc, n)
        if b.decay_rate_per_yr is not None:
            lam = b.decay_rate_per_yr
        elif b.solute_half_life_yr is not None:
            lam = bs.decay_rate_per_yr(b.solute_half_life_yr, s)
        else:
            lam = 0.0
            notes.append("No decay coefficient given: the first-order model equals the no-degradation model.")
        uf = dict(bs.UTILIZATION_FACTORS)
        if b.utilization_factors:
            uf.update(b.utilization_factors)
        ea = {"oxygen": b.delta_oxygen_mg_l, "nitrate": b.delta_nitrate_mg_l, "ferrous_iron": b.ferrous_iron_mg_l,
              "sulfate": b.delta_sulfate_mg_l, "methane": b.methane_mg_l}
        # Zones: input from the centre outward -> BIOSCREEN zones 3, 2, 1
        z = [(zz.width_m / FT, zz.concentration_mg_l) for zz in inp.source.zones]
        zones = [(0.0, 0.0)] * (3 - len(z)) + list(reversed(z))
        g = inp.general
        c = bs.Case(vs_ft_yr=vs_ft, porosity=n, alpha_x_ft=ax, alpha_y_ft=ay, alpha_z_ft=d.alpha_z_m / FT,
                    retardation=r, lambda_per_yr=lam, length_ft=g.model_length_m / FT, width_ft=g.model_width_m / FT,
                    time_yr=g.simulation_time_yr, thickness_ft=inp.source.thickness_m / FT, zones=zones,
                    mass_kg=inp.source.soluble_mass_kg, electron_acceptors=ea, utilization_factors=uf,
                    instantaneous_source=inp.instantaneous_source)
        if inp.chemical is not None and chem_name is None:
            chem_name = inp.chemical.name or inp.chemical.cas
        return c, s, notes, chem_name

    def run(self, inp: SourcePlumeInput) -> SourcePlumeOutput:
        c, s, notes, chem = self.case(inp)
        warnings: list[str] = []
        res = bs.run(c, s, inp.volume_threshold_mg_l)
        src = res.source
        v_ft = c.vs_ft_yr / c.retardation
        derived = Derived(
            seepage_velocity_m_yr=c.vs_ft_yr * FT, retardation_factor=c.retardation, contaminant_velocity_m_yr=v_ft * FT,
            alpha_x_m=c.alpha_x_ft * FT, alpha_y_m=c.alpha_y_ft * FT, alpha_z_m=c.alpha_z_ft * FT,
            decay_rate_per_yr=c.lambda_per_yr, biodegradation_capacity_mg_l=src.biodegradation_capacity_mg_l,
            source_width_m=src.source_width_ft * FT, source_flow_m3_yr=src.flow_l_yr / s.l_per_ft3 * FT3,
            source_concentration_avg_mg_l=src.c_avg_mg_l, source_flux_kg_yr=src.flux_kg_yr,
            source_decay_rate_per_yr=src.ks_per_yr, source_half_life_yr=_hl(src.half_life_yr),
            source_concentration_avg_inst_mg_l=src.c_avg_inst_mg_l, source_flux_inst_kg_yr=src.flux_inst_kg_yr,
            instantaneous_source=c.instantaneous_source,
            source_decay_rate_inst_per_yr=src.ks_inst_per_yr if c.instantaneous_source == "bioscreen" else None,
            source_half_life_inst_yr=_hl(src.half_life_inst_yr) if c.instantaneous_source == "bioscreen" else None,
            source_depletion_time_inst_yr=_hl(src.depletion_time_inst_yr) if c.instantaneous_source == "mass-balance"
            else None,
            layer_widths_m=[w * FT for w in src.widths_ft], layer_increments_mg_l=list(src.increments_mg_l),
            layer_increments_inst_mg_l=list(src.bioc_mg_l))
        x_m = [float(x) * FT for x in res.x_ft]
        x_m[0] = 0.0
        centerline = [CenterlinePoint(x_m=x, no_decay_mg_l=float(a), first_order_mg_l=float(b),
                                      instantaneous_mg_l=float(cc))
                      for x, a, b, cc in zip(x_m, res.centerline["no_decay"], res.centerline["first_order"],
                                            res.centerline["instantaneous"])]
        array = PlumeArray(x_m=x_m, y_m=[float(y) * FT for y in res.y_ft],
                           **{f"{k}_mg_l": res.arrays[k].tolist() for k in bs.MODELS})
        mb = []
        for k in bs.MODELS:
            m = res.models[k]
            mb.append(MassBalance(
                model=MODEL_NAMES[k], mass_left_source_kg=m.plume_mass_if_no_biodegradation_kg,
                actual_plume_mass_kg=m.actual_plume_mass_kg, mass_removed_kg=m.mass_removed_kg,
                fraction_removed=m.fraction_removed, electron_acceptor_change_kg=m.electron_acceptor_change_kg,
                source_mass_initial_kg=c.mass_kg, source_mass_now_kg=m.source_mass_now_kg,
                plume_volume_m3=None if m.plume_volume_ft3 is None else m.plume_volume_ft3 * FT3,
                source_flow_m3_yr=None if m.source_flow_ft3_yr is None else m.source_flow_ft3_yr * FT3,
                mass_flux_mg_day=m.mass_flux_mg_day, cant_calc=bool(m.cant_calc)))
        if any(m.cant_calc for m in res.models.values()):
            warnings.append("Mass balance: the no-degradation plume mass on the 11 x 5 array is outside 50-150 % of "
                            "the mass that has left the source ('Can't Calc', BSCR p. 40 [46]); enlarge the modelled "
                            "area so that the no-degradation plume lies within it.")
        integrated = None
        if c.alpha_z_ft == 0:
            integrated = self._integrate(c, s, src, inp.integration_points)
        else:
            warnings.append("Vertical dispersivity > 0: mass flux and integrated plume mass are not calculated "
                            "(the vertical concentration profile is not uniform; BSCR14 p. 2 [66]).")
        receptors = []
        for rc in inp.receptors:
            xr, yr = np.array([rc.x_m / FT]), np.array([rc.y_m / FT])
            cc = bs.concentrations(c, src, s, xr, yr)
            rr = ReceptorResult(name=rc.name, x_m=rc.x_m, y_m=rc.y_m, no_decay_mg_l=float(cc.no_decay[0]),
                                first_order_mg_l=float(cc.first_order[0]), instantaneous_mg_l=float(cc.instantaneous[0]))
            for t in inp.receptor_times_yr:
                ct = bs.concentrations(c, src, s, xr, yr, t_yr=t)
                rr.time_series.append({"time_yr": t, "no_decay_mg_l": float(ct.no_decay[0]),
                                       "first_order_mg_l": float(ct.first_order[0]),
                                       "instantaneous_mg_l": float(ct.instantaneous[0])})
            receptors.append(rr)
            if rc.x_m > inp.general.model_length_m:
                notes.append(f"Receptor {rc.name} lies beyond the modelled length.")
        if c.retardation > 6:
            warnings.append(f"R = {c.retardation:.3g} > 6: the instantaneous reaction (superposition) model should be "
                            "used with caution and checked with a numerical model (BSCR p. 21 [27], p. 49 [55]).")
        if c.mass_kg is not None and src.biodegradation_capacity_mg_l > 0:
            if c.instantaneous_source == "bioscreen":
                notes.append("Declining source with the instantaneous model (BIOSCREEN): the source of C0 + BC decays "
                             "with the shorter half-life while BC stays constant, so the source concentration can "
                             "reach zero before the source mass is depleted (BSCR p. 35 [41]; docs/methods/napl.md §4.8).")
            else:
                notes.append(f"Instantaneous model, mass-balance source: source mass and concentration reach zero "
                             f"at t = {src.depletion_time_inst_yr:.4g} yr (docs/methods/napl.md §4.8).")
        if s.compatible:
            notes.append("bioscreen_compatible: BIOSCREEN 1.4 constants and workbook logic reproduced exactly "
                         "(docs/methods/napl.md §4.6, Q1-Q3).")
        return SourcePlumeOutput(bioscreen_compatible=s.compatible, chemical=chem, derived=derived,
                                 centerline=centerline, array=array, mass_balance=mb, integrated=integrated,
                                 receptors=receptors, notes=notes, warnings=warnings)

    @staticmethod
    def _integrate(c: bs.Case, s: bs.Settings, src: bs.Source, npts: int) -> list[Integrated]:
        """Trapezoidal integration of C(x, y) on an npts x npts grid over the modelled area (0 <= x <= L,
        |y| <= W/2); mass = n Z R integral (2-D, as BIOSCREEN, BSCR p. 39 [45]); flux = Vs n Z integral C dy."""
        x = np.linspace(0.0, c.length_ft, npts)
        y = np.linspace(-c.width_ft / 2.0, c.width_ft / 2.0, npts)
        X, Y = np.meshgrid(x, y)
        cc = bs.concentrations(c, src, s, X, Y)
        xs = np.array([j * c.length_ft / 10.0 for j in range(11)])
        XS, YS = np.meshgrid(xs, y)
        cs = bs.concentrations(c, src, s, XS, YS)
        out = []
        l_per_ft3 = bs.L_PER_FT3_EXACT     # integration results use the exact conversion
        for k in bs.MODELS:
            a = getattr(cc, k)
            mass = np.trapezoid(np.trapezoid(a, y, axis=0), x) * c.thickness_ft * c.porosity * c.retardation \
                * l_per_ft3 / 1e6
            flux = np.trapezoid(getattr(cs, k), y, axis=0) * c.vs_ft_yr * c.porosity * c.thickness_ft * l_per_ft3 \
                / bs.DAYS_PER_YEAR
            out.append(Integrated(model=MODEL_NAMES[k], plume_mass_kg=float(mass), mass_flux_mg_day=flux.tolist()))
        return out

    def summary(self, out: SourcePlumeOutput) -> str:
        lines = [f"napl.source-plume ({'BIOSCREEN 1.4 compatible' if out.bioscreen_compatible else 'consistent'})"]
        if out.chemical:
            lines.append(f"Chemical: {out.chemical}")
        d = out.derived
        hl = "infinite" if d.source_half_life_yr is None else f"{d.source_half_life_yr:.3g} yr"
        if d.instantaneous_source == "bioscreen":
            hli = "infinite" if d.source_half_life_inst_yr is None else f"half-life {d.source_half_life_inst_yr:.3g} yr"
        else:
            hli = ("not depleted" if d.source_depletion_time_inst_yr is None
                   else f"depleted at {d.source_depletion_time_inst_yr:.3g} yr")
        lines.append(f"Vs = {d.seepage_velocity_m_yr:.4g} m/yr, R = {d.retardation_factor:.4g}, alpha x/y/z = "
                     f"{d.alpha_x_m:.3g}/{d.alpha_y_m:.3g}/{d.alpha_z_m:.3g} m, lambda = {d.decay_rate_per_yr:.4g} 1/yr, "
                     f"BC = {d.biodegradation_capacity_mg_l:.4g} mg/L")
        lines.append(f"Source: Q = {d.source_flow_m3_yr:.4g} m3/yr, flux = {d.source_flux_kg_yr:.4g} kg/yr, "
                     f"half-life {hl} (inst. reaction, {d.instantaneous_source} source: {hli})")
        lines.append(f"{'x (m)':>10} {'no decay':>11} {'1st order':>11} {'inst.':>11}   (mg/L, centreline)")
        for p in out.centerline:
            lines.append(f"{p.x_m:10.2f} {p.no_decay_mg_l:11.4g} {p.first_order_mg_l:11.4g} {p.instantaneous_mg_l:11.4g}")
        for m in out.mass_balance:
            if m.cant_calc:
                lines.append(f"{m.model}: mass balance can't be calculated (enlarge model area)")
            else:
                lines.append(f"{m.model}: left source {m.mass_left_source_kg:.4g} kg, in plume "
                             f"{m.actual_plume_mass_kg:.4g} kg, removed {100 * m.fraction_removed:.3g} %")
        for r in out.receptors:
            lines.append(f"Receptor {r.name} (x = {r.x_m:g} m, y = {r.y_m:g} m): {r.no_decay_mg_l:.4g} / "
                         f"{r.first_order_mg_l:.4g} / {r.instantaneous_mg_l:.4g} mg/L")
        lines += [f"Warning: {w}" for w in out.warnings] + [f"Note: {n}" for n in out.notes]
        return "\n".join(lines)
