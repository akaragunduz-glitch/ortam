"""VLEACH 2.2a: one-dimensional vadose zone leaching of a sorbing, volatile contaminant.

Source: Ravi, V. and Johnson, J.A., "VLEACH, A One-Dimensional Finite Difference
Vadose Zone Leaching Model, Version 2.2", U.S. EPA Center for Subsurface Modeling
Support (docs/references/vleach_2.2_users_guide.pdf; cited "VLEACH UG", printed pages).

Physics (VLEACH UG section 3, pp. 9-10):
  - three-phase linear equilibrium in every cell: total concentration
        C_T = theta*C_l + (n - theta)*C_g + rho_b*C_s,
        C_g = H*C_l, C_s = K_d*C_l, K_d = K_oc*f_oc          (Eq. 4-8)
  - advection of the liquid phase with the recharge (Darcy) flux q (Eq. 2)
  - diffusion of the gas phase with D = D_air*(n - theta)^(10/3)/n^2  (Eq. 9)
  - no dispersion, no degradation, homogeneous soil per polygon (pp. 6-7)
  - initial condition: the measured soil concentration is the total
    concentration of the cell, distributed over the three phases.

Numerical scheme. The guide only says "implicit finite difference" (p. 10). The
details below were identified by numerical experiments with EPA's own executable
(VleachForGUI.exe of the VLEACH 2.2a setup, run under Wine) and are verified
against its output files to their printed precision (modules/vadose/tests/test_vleach.py):
  1. liquid advection, upwind Crank-Nicolson in the liquid phase only:
        (1 + w/2) C_l,i' = (1 - w/2) C_l,i + (w/2)(C_l,i-1 + C_l,i-1'),  w = q*dt/(theta*dz)
     with C_l,0 = concentration in recharge water;
  2. gas diffusion, fully implicit in the gas phase only, storage (n - theta):
        (n - theta) dz (C_g' - C_g)/dt = D (C_g,i+1' - 2 C_g,i' + C_g,i-1')/dz,
     boundary cells exchange with the atmosphere (gas concentration C_atm) or the
     water table (gas concentration H*C_gw) over a distance dz; a negative input
     value makes the boundary impermeable;
  3. both mass changes are computed from the same old state and added
     (additive operator splitting); the cell's new total mass is then
     redistributed over the three phases (Eq. 5-7).
  Note: the D formula uses the exponent 10/3 of the guide (Eq. 9); the 2.2b
  release notes mention 7/3, but the 2.2a executable uses 10/3.
Groundwater impact (VLEACH UG p. 21): the mass flux across the water table
(advection plus diffusion) in the last time step before each printout, per year;
the cumulative mass is the sum over all time steps (as in VLEACH); the flux of
every time step is also given (VLEACH file GWIMP.DAT).

Units here: SI (m, yr, g, g/m3 = mg/L), soil concentrations in ug/kg as in VLEACH.
D_air in m2/day is converted with 365 day/yr, as in VLEACH.
Fortran analogy: one DO loop over time steps; inside, two tridiagonal sweeps
(advection: forward substitution; diffusion: Thomas algorithm).
"""

from __future__ import annotations

import math

from pydantic import BaseModel, Field, model_validator

from ortam_core.model import ModelInput, ModelOutput, Module

DAYS_PER_YEAR = 365.0          # VLEACH conversion of D_air (m2/day -> per year)
# VLEACH works in feet and converts D_air with 3.2809 ft/m (the .PRM file prints
# 0.7 m2/day as 2750.3 ft2/yr). In metres this is a factor (3.2809*0.3048)^2 = 1.0000366
# on D, kept so that results reproduce the EPA executable (effect < 0.004 % on D).
VLEACH_D_FACTOR = (3.2809 * 0.3048) ** 2


# ------------------------------------------------------------------ inputs
class Layer(ModelInput):
    top_m: float = Field(ge=0, description="Top of the interval below ground surface, m")
    bottom_m: float = Field(gt=0, description="Bottom of the interval, m")
    soil_ug_kg: float = Field(ge=0, description="Initial total soil concentration, ug/kg (dry weight). "
                                                "A modelled clean interval is 0; a non-detect is not.")


class Polygon(ModelInput):
    name: str = "Polygon 1"
    area_m2: float = Field(gt=0, description="Polygon area, m2")
    cell_thickness_m: float = Field(gt=0, description="Vertical cell size dz, m")
    depth_to_water_m: float = Field(gt=0, description="Thickness of the vadose zone, m (multiple of dz)")
    recharge_m_yr: float = Field(ge=0, description="Recharge (Darcy) flux q, m/yr")
    bulk_density_g_cm3: float = Field(gt=0, description="Dry bulk density, g/cm3")
    porosity: float = Field(gt=0, lt=1, description="Porosity n, -")
    water_content: float = Field(gt=0, lt=1, description="Volumetric water content theta, -")
    foc: float = Field(ge=0, lt=1, description="Fraction of organic carbon, -")
    recharge_conc_mg_l: float = Field(0.0, ge=0, description="Concentration in recharge water, mg/L")
    atmosphere_gas_mg_l: float | None = Field(0.0, ge=0, description="Gas concentration above the surface, mg/L; "
                                                                     "null = surface impermeable to gas")
    groundwater_mg_l: float | None = Field(0.0, ge=0, description="Concentration in groundwater, mg/L; "
                                                                  "null = water table impermeable to gas")
    initial: list[Layer] = Field(min_length=1, description="Initial soil concentration intervals, from the top")

    @model_validator(mode="after")
    def _check(self):
        if self.water_content >= self.porosity:
            raise ValueError("water_content must be smaller than porosity (VLEACH UG p. 17).")
        n = self.depth_to_water_m / self.cell_thickness_m
        if abs(n - round(n)) > 1e-6:
            raise ValueError("depth_to_water_m must be a multiple of cell_thickness_m.")
        return self


class Chemical(ModelInput):
    name: str
    koc_ml_g: float = Field(ge=0, description="Organic carbon partition coefficient Koc, mL/g")
    henry_dimensionless: float = Field(ge=0, description="Dimensionless Henry's law constant H, -")
    solubility_mg_l: float = Field(gt=0, description="Water solubility, mg/L (checked, not enforced, as in VLEACH)")
    dair_m2_d: float = Field(ge=0, description="Free-air diffusion coefficient, m2/day")


class VleachInput(ModelInput):
    title: str = ""
    time_step_yr: float = Field(gt=0, description="Time step, yr")
    simulation_time_yr: float = Field(gt=0, description="Simulation length, yr (multiple of the time step)")
    output_interval_yr: float = Field(gt=0, description="Mass balance / groundwater impact interval, yr")
    profile_interval_yr: float = Field(gt=0, description="Concentration profile interval, yr")
    chemical: Chemical
    polygons: list[Polygon] = Field(min_length=1)

    @model_validator(mode="after")
    def _times(self):
        for name in ("simulation_time_yr", "output_interval_yr"):
            k = getattr(self, name) / self.time_step_yr
            if abs(k - round(k)) > 1e-6:
                raise ValueError(f"{name} must be a multiple of time_step_yr (VLEACH requirement).")
        return self


# ------------------------------------------------------------------ outputs
class MassBalance(BaseModel):
    time_yr: float
    total_g_m2: float
    gas_g_m2: float
    liquid_g_m2: float
    sorbed_g_m2: float
    advection_in_atmosphere_g_m2: float = Field(description="Since the previous output")
    advection_in_water_table_g_m2: float
    diffusion_in_atmosphere_g_m2: float
    diffusion_in_water_table_g_m2: float
    discrepancy_g_m2: float
    groundwater_flux_g_m2_yr: float = Field(description="Flux to groundwater in the last time step (VLEACH 'impact')")


class Profile(BaseModel):
    time_yr: float
    depth_m: list[float] = Field(description="Cell centre depth, m")
    gas_mg_l: list[float]
    liquid_mg_l: list[float]
    sorbed_mg_kg: list[float]


class PolygonResult(BaseModel):
    name: str
    area_m2: float
    retardation: float = Field(description="1 + (rho_b*Kd + (n - theta)*H)/theta, -")
    deff_m2_yr: float
    mass_balance: list[MassBalance]
    profiles: list[Profile]
    cumulative_to_groundwater_g: float = Field(description="Cumulative mass leached to groundwater, g")
    groundwater_series: list[tuple[float, float]] = Field(
        description="(time yr, mass rate to groundwater g/yr) for every time step (VLEACH GWIMP.DAT)")


class VleachOutput(ModelOutput):
    title: str
    polygons: list[PolygonResult]
    total_impact: list[dict] = Field(description="time_yr, rate_g_yr (last step), cumulative_g (all polygons)")


# ------------------------------------------------------------------ numerics
def _advect(cl: list[float], w: float, c_in: float) -> list[float]:
    """Upwind Crank-Nicolson step of the liquid phase (forward substitution)."""
    new = [0.0] * len(cl)
    prev_old = prev_new = c_in
    for i, c in enumerate(cl):
        new[i] = ((1.0 - w / 2.0) * c + (w / 2.0) * (prev_old + prev_new)) / (1.0 + w / 2.0)
        prev_old, prev_new = c, new[i]
    return new


def _diffuse(cg: list[float], lam: float, top: float | None, bottom: float | None) -> list[float]:
    """Fully implicit diffusion step of the gas phase (Thomas algorithm). lam = D*dt/((n-theta)*dz^2)."""
    n = len(cg)
    a = [0.0] * n; b = [0.0] * n; c = [0.0] * n; d = list(cg)
    for i in range(n):
        b[i] = 1.0
        if i > 0:
            a[i] = -lam; b[i] += lam
        if i < n - 1:
            c[i] = -lam; b[i] += lam
    if top is not None:
        b[0] += lam; d[0] += lam * top
    if bottom is not None:
        b[-1] += lam; d[-1] += lam * bottom
    for i in range(1, n):                       # forward elimination
        m = a[i] / b[i - 1]
        b[i] -= m * c[i - 1]; d[i] -= m * d[i - 1]
    x = [0.0] * n
    x[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):              # back substitution
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]
    return x


def _initial_profile(p: Polygon, ncell: int) -> list[float]:
    conc = [None] * ncell
    dz = p.cell_thickness_m
    for lay in p.initial:
        for i in range(ncell):
            mid = (i + 0.5) * dz
            if lay.top_m <= mid < lay.bottom_m:
                conc[i] = lay.soil_ug_kg
    if any(c is None for c in conc):
        raise ValueError(f"{p.name}: the initial intervals must cover the whole vadose zone.")
    return conc


def run_polygon(p: Polygon, ch: Chemical, dt: float, nsteps: int, out_every: int, prf_interval: float) -> PolygonResult:
    dz, q = p.cell_thickness_m, p.recharge_m_yr
    th, n_ = p.water_content, p.porosity
    a = n_ - th                                       # air-filled porosity
    rho = p.bulk_density_g_cm3 * 1.0e6                # g/m3
    kd = ch.koc_ml_g * p.foc * 1.0e-6                 # m3/g  (Eq. 8)
    H = ch.henry_dimensionless
    B = th + a * H + rho * kd                         # C_T / C_l  (Eq. 4)
    D = ch.dair_m2_d * DAYS_PER_YEAR * VLEACH_D_FACTOR * a ** (10.0 / 3.0) / n_ ** 2   # m2/yr (Eq. 9)
    ncell = round(p.depth_to_water_m / dz)
    c_in = p.recharge_conc_mg_l                       # mg/L = g/m3
    top = p.atmosphere_gas_mg_l
    bot = p.groundwater_mg_l * H if p.groundwater_mg_l is not None else None
    w = q * dt / (th * dz)
    lam = D * dt / (a * dz * dz) if a > 0 else 0.0

    soil = _initial_profile(p, ncell)
    M = [s * rho * 1.0e-9 for s in soil]              # ug/kg * g/m3 * (1e-6 g/ug)(1e-3 kg/g) -> g/m3
    cl = [m / B for m in M]

    def masses():
        tot = sum(M) * dz
        lq = sum(cl) * th * dz
        return tot, lq * a * H / th, lq, sum(cl) * rho * kd * dz

    def profile(t):
        return Profile(time_yr=t, depth_m=[(i + 0.5) * dz for i in range(ncell)],
                       gas_mg_l=[c * H for c in cl], liquid_mg_l=list(cl),
                       sorbed_mg_kg=[c * kd * 1.0e6 for c in cl])

    tot0 = masses()[0]
    acc = dict(adv_top=0.0, adv_bot=0.0, dif_top=0.0, dif_bot=0.0)
    last_total = tot0
    mb = [MassBalance(time_yr=0.0, total_g_m2=tot0, gas_g_m2=masses()[1], liquid_g_m2=masses()[2],
                      sorbed_g_m2=masses()[3], advection_in_atmosphere_g_m2=0.0, advection_in_water_table_g_m2=0.0,
                      diffusion_in_atmosphere_g_m2=0.0, diffusion_in_water_table_g_m2=0.0, discrepancy_g_m2=0.0,
                      groundwater_flux_g_m2_yr=0.0)]
    profiles = [profile(0.0)]
    cum_gw = 0.0
    series: list[tuple[float, float]] = []
    for k in range(1, nsteps + 1):
        d_adv = [0.0] * ncell
        step_bot = 0.0
        if q > 0:
            new = _advect(cl, w, c_in)
            d_adv = [th * (x - y) for x, y in zip(new, cl)]
            acc["adv_top"] += q * dt * c_in
            out = q * dt * (cl[-1] + new[-1]) / 2.0
            acc["adv_bot"] -= out
            step_bot -= out
        d_dif = [0.0] * ncell
        if lam > 0:
            cg = [c * H for c in cl]
            cg2 = _diffuse(cg, lam, top, bot)
            d_dif = [a * (x - y) for x, y in zip(cg2, cg)]
            if top is not None:
                acc["dif_top"] += D * (top - cg2[0]) / dz * dt
            if bot is not None:
                f = D * (bot - cg2[-1]) / dz * dt
                acc["dif_bot"] += f
                step_bot += f
        M = [m + x + y for m, x, y in zip(M, d_adv, d_dif)]
        cl = [m / B for m in M]
        cum_gw -= step_bot * p.area_m2
        t = k * dt
        series.append((t, -step_bot * p.area_m2 / dt))
        if k % out_every == 0:
            tot, gas, liq, sor = masses()
            inflow = sum(acc.values())
            mb.append(MassBalance(time_yr=t, total_g_m2=tot, gas_g_m2=gas, liquid_g_m2=liq, sorbed_g_m2=sor,
                                  advection_in_atmosphere_g_m2=acc["adv_top"],
                                  advection_in_water_table_g_m2=acc["adv_bot"],
                                  diffusion_in_atmosphere_g_m2=acc["dif_top"],
                                  diffusion_in_water_table_g_m2=acc["dif_bot"],
                                  discrepancy_g_m2=(tot - last_total) - inflow,
                                  groundwater_flux_g_m2_yr=-step_bot / dt))
            last_total = tot
            acc = dict.fromkeys(acc, 0.0)
        r = t / prf_interval
        if abs(r - round(r)) < 1e-6 or k == nsteps:
            profiles.append(profile(t))
    return PolygonResult(name=p.name, area_m2=p.area_m2, retardation=B / th, deff_m2_yr=D, mass_balance=mb,
                         profiles=profiles, cumulative_to_groundwater_g=cum_gw, groundwater_series=series)


class VleachModule(Module):
    name = "vadose.vleach"
    title = "Vadose zone: VLEACH 2.2a one-dimensional leaching (sorption, volatilization)"
    input_type = VleachInput
    output_type = VleachOutput

    def run(self, inp: VleachInput) -> VleachOutput:
        dt = inp.time_step_yr
        nsteps = round(inp.simulation_time_yr / dt)
        out_every = round(inp.output_interval_yr / dt)
        warnings: list[str] = []
        k = inp.profile_interval_yr / dt
        if abs(k - round(k)) > 1e-6:
            warnings.append("profile_interval_yr is not a multiple of the time step: profiles are given at the "
                            "steps that fall on a multiple of it and at the end of the run (as in VLEACH).")
        prf_times = inp.profile_interval_yr
        results = []
        for p in inp.polygons:
            if p.recharge_m_yr * dt / p.cell_thickness_m > 1.0:
                warnings.append(f"{p.name}: Courant number q*dt/dz = {p.recharge_m_yr * dt / p.cell_thickness_m:.3g} "
                                "> 1; reduce the time step (numerical dispersion, possible negative values).")
            r = run_polygon(p, inp.chemical, dt, nsteps, out_every, prf_times)
            cmax = max(max(pr.liquid_mg_l) for pr in r.profiles)
            if cmax > inp.chemical.solubility_mg_l:
                warnings.append(f"{p.name}: liquid concentration {cmax:.3g} mg/L exceeds the solubility "
                                f"{inp.chemical.solubility_mg_l:g} mg/L (free product is not modelled).")
            results.append(r)
        total = []
        for j, mb in enumerate(results[0].mass_balance[1:], start=1):
            rate = sum(r.mass_balance[j].groundwater_flux_g_m2_yr * r.area_m2 for r in results)
            k_end = round(mb.time_yr / dt)
            cum = sum(rr * dt for r in results for _, rr in r.groundwater_series[:k_end])
            total.append({"time_yr": mb.time_yr, "rate_g_yr": rate, "cumulative_g": cum})
        return VleachOutput(title=inp.title, polygons=results, total_impact=total, warnings=warnings)

    def summary(self, o: VleachOutput) -> str:
        lines = [f"VLEACH 2.2a  {o.title}"]
        for r in o.polygons:
            lines.append(f"{r.name}: R = {r.retardation:.3g}, D_eff = {r.deff_m2_yr:.3g} m2/yr; "
                         f"cumulative to groundwater = {r.cumulative_to_groundwater_g:.4g} g")
            lines.append(f"  {'t (yr)':>8s} {'total g/m2':>11s} {'to GW g/m2/yr':>14s}")
            for m in r.mass_balance:
                lines.append(f"  {m.time_yr:8.4g} {m.total_g_m2:11.4g} {m.groundwater_flux_g_m2_yr:14.4g}")
        lines.append("Total groundwater impact: " + ", ".join(f"{t['time_yr']:g} yr: {t['rate_g_yr']:.4g} g/yr"
                                                            for t in o.total_impact))
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
