"""surfacewater.watershed: COPC loads from a watershed and concentrations in a water body.

HHRAP 2005 (EPA530-R-05-006) Appendix B, Tables B-4-1 to B-4-25 (equations: equations.py):
  watershed soil  Cs measured, or from air deposition (B-4-1 with ks from B-4-2 to B-4-6)
  loads           LT = LDEP + Ldif + LRI + LR + LE (B-4-7 to B-4-14) + external loads
  water body      Cwtot, fwc/fbs, kwt, kv, Kv, KL/KG (river or lake), kb (B-4-15 to B-4-22)
  concentrations  Cwctot, Cdw (dissolved), Csb (bed sediment) (B-4-23 to B-4-25)
Fish concentrations (B-4-26 to B-4-28) are in foodchain.media.

External loads (g/yr, e.g. groundwater discharge from napl.source-plume, a point discharge) are an
ORTAM extension: they are added to the sum of Table B-4-7 (docs/methods/surfacewater.md).
Mercury (total Hg emission): 0.48 Q and Fv = 0.85; soil split 98 % Hg2+ / 2 % methyl mercury; runoff
and erosion loads of both species; water body with mercuric chloride parameters; Cdw split
85 % / 15 % (Tables B-4-1, -8, -9, -11, -12, -24).
Chemical parameters: HHRAP Companion Database (data/reference/hhrap/2005/parameters.csv) unless given.
Verification: modules/surfacewater/tests (independent LibreOffice hand calculation).
"""
from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, Field, model_validator

from modules.surfacewater import equations as eq
from ortam_core.model import ModelInput, ModelOutput, Module

DATA = Path(__file__).resolve().parents[2] / "data" / "reference" / "hhrap" / "2005"
DB_KEYS = {"Kds": "kds_ml_g", "Kdsw": "kdsw_l_kg", "Kdbs": "kdbs_l_kg", "H": "henry_atm_m3_mol",
           "Dw": "dw_cm2_s", "Da": "da_cm2_s", "ksg": "ksg_per_yr", "Fv": "fv"}
MERCURY_CAS = "7439-97-6"


@lru_cache
def defaults() -> dict:
    return yaml.safe_load((DATA / "defaults.yaml").read_text(encoding="utf-8"))["surface_water"]


@lru_cache
def hhrap_params() -> dict[str, dict]:
    out: dict[str, dict] = {}
    with (DATA / "parameters.csv").open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            d = out.setdefault(r["cas"], {"name": r["name"]})
            if r["parameter"] in DB_KEYS and r["value"] != "":
                d[DB_KEYS[r["parameter"]]] = float(r["value"])
    return out


def _d(key: str) -> float:
    return defaults()[key]["value"]


class ChemicalInput(ModelInput):
    cas: str | None = None
    name: str | None = None
    organic: bool | None = Field(None, description="Organic COPC (enrichment ratio 3) or inorganic (1); "
                                 "default from the RSL table (Koc: organic, Kd: inorganic)")
    kds_ml_g: float | None = Field(None, ge=0, description="Soil-water partition coefficient Kds, mL/g")
    kdsw_l_kg: float | None = Field(None, ge=0, description="Suspended sediment-water partition coefficient Kdsw, L/kg")
    kdbs_l_kg: float | None = Field(None, ge=0, description="Bed sediment-pore water partition coefficient Kdbs, L/kg")
    henry_atm_m3_mol: float | None = Field(None, ge=0, description="Henry's law constant H, atm-m3/mol")
    dw_cm2_s: float | None = Field(None, gt=0, description="Diffusivity in water Dw, cm2/s")
    da_cm2_s: float | None = Field(None, gt=0, description="Diffusivity in air Da, cm2/s")
    ksg_per_yr: float | None = Field(None, ge=0, description="Soil degradation constant ksg, 1/yr")
    fv: float | None = Field(None, ge=0, le=1, description="Fraction in the vapor phase Fv, -")


class Usle(ModelInput):
    rainfall_rf_per_yr: float = Field(gt=0, description="USLE rainfall (erosivity) factor RF, 1/yr (50-300)")
    erodibility_k_ton_acre: float | None = Field(None, gt=0, description="USLE erodibility K, ton/acre (default 0.39)")
    length_slope_ls: float | None = Field(None, gt=0, description="USLE length-slope factor LS (default 1.5)")
    cover_management_c: float | None = Field(None, gt=0, le=1, description="USLE cover management C (default 0.1)")
    supporting_practice_pf: float | None = Field(None, gt=0, le=1, description="USLE supporting practice PF (default 1.0)")


class Watershed(ModelInput):
    area_m2: float = Field(gt=0, description="Total watershed area receiving deposition AL, m2")
    impervious_area_m2: float = Field(0.0, ge=0, description="Impervious watershed area AI, m2")
    soil_mixing_depth_cm: float | None = Field(None, gt=0, description="Zs, cm (default 2 untilled; 20 tilled)")
    bulk_density_g_cm3: float | None = Field(None, gt=0, description="BD, g/cm3 (default 1.5)")
    soil_water_content: float | None = Field(None, gt=0, lt=1, description="theta_sw, mL/cm3 (default 0.2)")
    particle_density_g_cm3: float | None = Field(None, gt=0, description="rho_soil, g/cm3 (default 2.7)")
    air_temperature_k: float | None = Field(None, gt=0, description="Ta, K (default 298)")
    runoff_cm_yr: float = Field(ge=0, description="Average annual surface runoff from pervious areas RO, cm/yr")
    precipitation_cm_yr: float = Field(ge=0, description="Average annual precipitation P, cm/yr")
    irrigation_cm_yr: float | None = Field(None, ge=0, description="Average annual irrigation I, cm/yr (default 0)")
    evapotranspiration_cm_yr: float = Field(ge=0, description="Average annual evapotranspiration Ev, cm/yr")
    unit_soil_loss_kg_m2_yr: float | None = Field(None, ge=0, description="Xe, kg/m2-yr (or give usle)")
    usle: Usle | None = None
    sediment_delivery_ratio: float | None = Field(None, gt=0, le=1, description="SD (default a AL^-b, Table B-4-14)")
    sd_b: float | None = Field(None, gt=0, description="Empirical slope coefficient b (default 0.125)")
    enrichment_ratio: float | None = Field(None, gt=0, description="ER (default 3 organic, 1 inorganic)")
    erosion_loss_constant: Literal["zero", "calculated"] = Field(
        "zero", description="kse: 'zero' (HHRAP default) or 'calculated' (Table B-4-3)")

    @model_validator(mode="after")
    def _check(self):
        if (self.unit_soil_loss_kg_m2_yr is None) == (self.usle is None):
            raise ValueError("Give unit_soil_loss_kg_m2_yr or usle (not both).")
        if self.impervious_area_m2 > self.area_m2:
            raise ValueError("impervious_area_m2 exceeds area_m2.")
        return self


class Soil(ModelInput):
    concentration_mg_kg: float | None = Field(None, ge=0, description="Measured watershed soil concentration, mg/kg "
                                              "(empty: from air deposition, Table B-4-1). Non-detects are not entered.")
    deposition_period_yr: float | None = Field(None, gt=0, description="tD, yr (default 30)")
    t1_yr: float | None = Field(None, ge=0, description="T1, yr (default 0)")
    exposure_duration_yr: float | None = Field(None, gt=0, description="T2, yr (default 30)")


class Air(ModelInput):
    emission_rate_g_s: float = Field(ge=0, description="COPC emission rate Q, g/s (mercury: total Hg)")
    watershed_dytwv_s_m2_yr: float = Field(0.0, ge=0, description="Unitized deposition, vapor, watershed, s/m2-yr")
    watershed_dytwp_s_m2_yr: float = Field(0.0, ge=0, description="Unitized deposition, particle, watershed, s/m2-yr")
    water_body_dytwv_s_m2_yr: float = Field(0.0, ge=0, description="Unitized deposition, vapor, water body, s/m2-yr")
    water_body_dytwp_s_m2_yr: float = Field(0.0, ge=0, description="Unitized deposition, particle, water body, s/m2-yr")
    water_body_cywv_ug_s_g_m3: float = Field(0.0, ge=0, description="Unitized air concentration, vapor, over the "
                                             "water body, ug-s/g-m3")


class WaterBody(ModelInput):
    type: Literal["river", "lake"]
    area_m2: float = Field(gt=0, description="Water body surface area Aw, m2")
    flow_m3_yr: float = Field(ge=0, description="Average volumetric flow rate Vfx, m3/yr")
    water_column_depth_m: float = Field(gt=0, description="dwc, m")
    benthic_depth_m: float | None = Field(None, gt=0, description="dbs, m (default 0.03)")
    tss_mg_l: float | None = Field(None, gt=0, description="Total suspended solids, mg/L (default 10)")
    bed_sediment_concentration_g_cm3: float | None = Field(None, gt=0, description="CBS, g/cm3 (default 1.0)")
    bed_sediment_porosity: float | None = Field(None, gt=0, lt=1, description="theta_bs (default 0.6)")
    temperature_k: float | None = Field(None, gt=0, description="Twk, K (default 298)")
    theta: float | None = Field(None, gt=0, description="Temperature correction factor (default 1.026)")
    current_velocity_m_s: float | None = Field(None, gt=0, description="u, m/s (river)")
    wind_speed_m_s: float | None = Field(None, gt=0, description="W, m/s (lake; default 3.9)")
    drag_coefficient: float | None = Field(None, gt=0, description="Cd (default 0.0011)")
    von_karman: float | None = Field(None, gt=0, description="k (default 0.4)")
    viscous_sublayer: float | None = Field(None, gt=0, description="lambda_z (default 4)")
    air_density_g_cm3: float | None = Field(None, gt=0, description="rho_a (default 0.0012)")
    water_density_g_cm3: float | None = Field(None, gt=0, description="rho_w (default 1.0)")
    air_viscosity_g_cm_s: float | None = Field(None, gt=0, description="mu_a (default 1.81e-4)")
    water_viscosity_g_cm_s: float | None = Field(None, gt=0, description="mu_w (default 1.69e-2)")

    @model_validator(mode="after")
    def _check(self):
        if self.type == "river" and self.current_velocity_m_s is None:
            raise ValueError("A river needs current_velocity_m_s (Table B-4-20).")
        return self


class ExternalLoad(ModelInput):
    name: str
    load_g_yr: float = Field(ge=0, description="Load to the water body, g/yr")


class WatershedInput(ModelInput):
    chemical: ChemicalInput
    chemical_table: Literal["hhrap-2005", "none"] = "hhrap-2005"
    watershed: Watershed
    soil: Soil = Soil()
    air: Air | None = None
    water_body: WaterBody
    external_loads: list[ExternalLoad] = []

    @model_validator(mode="after")
    def _check(self):
        if self.soil.concentration_mg_kg is None and self.air is None:
            raise ValueError("Give soil.concentration_mg_kg or air deposition (Table B-4-1).")
        return self


class SoilResult(BaseModel):
    species: str
    kse_per_yr: float
    ksr_per_yr: float
    ksl_per_yr: float
    ksv_per_yr: float
    ksg_per_yr: float
    ks_per_yr: float
    ds_mg_kg_yr: float | None = None
    cs_td_mg_kg: float | None = Field(None, description="CstD, noncarcinogens (Table B-4-1)")
    cs_average_mg_kg: float | None = Field(None, description="Cs averaged over T2, carcinogens (Table B-4-1)")
    cs_measured_mg_kg: float | None = None


class Loads(BaseModel):
    deposition_g_yr: float
    diffusion_g_yr: float
    impervious_runoff_g_yr: float
    pervious_runoff_g_yr: float
    erosion_g_yr: float
    external_g_yr: float
    total_g_yr: float


class Concentrations(BaseModel):
    basis: str = Field(description="'measured', 'carcinogens' (Cs) or 'noncarcinogens' (CstD)")
    soil_mg_kg: dict[str, float]
    loads: Loads
    total_water_body_mg_l: float = Field(description="Cwtot, water column and bed sediment (Table B-4-15)")
    total_water_column_mg_l: float = Field(description="Cwctot (Table B-4-23)")
    dissolved_mg_l: float = Field(description="Cdw (Table B-4-24)")
    bed_sediment_mg_kg: float = Field(description="Csb (Table B-4-25)")
    dissolved_species_mg_l: dict[str, float] | None = Field(None, description="Mercury: Hg2+ 85 %, MHg 15 %")


class Transfer(BaseModel):
    unit_soil_loss_kg_m2_yr: float
    sd_intercept_a: float | None = None
    sediment_delivery_ratio: float
    kl_m_yr: float
    kg_m_yr: float
    kv_overall_m_yr: float
    kv_per_yr: float
    kb_calculated_per_yr: float
    kb_per_yr: float
    fwc: float
    fbs: float
    kwt_per_yr: float


class WatershedOutput(ModelOutput):
    chemical: str
    mercury: bool
    parameters: dict[str, dict[str, float]]
    soil: list[SoilResult]
    transfer: Transfer
    results: list[Concentrations]
    notes: list[str] = []


def _params(ci: ChemicalInput, cas: str | None, table: str) -> dict[str, float]:
    p = dict(hhrap_params().get(cas, {})) if (table == "hhrap-2005" and cas) else {}
    p.pop("name", None)
    for k in DB_KEYS.values():
        v = getattr(ci, k)
        if v is not None:
            p[k] = v
    missing = [k for k in DB_KEYS.values() if k not in p]
    if missing:
        raise ValueError(f"Missing chemical parameters for {ci.name or cas}: {', '.join(missing)}.")
    return p


def _organic(ci: ChemicalInput, cas: str | None) -> bool:
    if ci.organic is not None:
        return ci.organic
    if cas:
        try:
            from ortam_core.chemicals.database import lookup
            c = lookup("rsl-2024-11", cas, None)
            if c.koc_l_kg is not None:
                return True
            if c.kd_l_kg is not None:
                return False
        except Exception:
            pass
    raise ValueError("Set chemical.organic (enrichment ratio 3 for organic, 1 for inorganic COPCs).")


class WatershedModule(Module):
    name = "surfacewater.watershed"
    title = "Surface water: watershed loads and water body concentrations (HHRAP 2005, Tables B-4-1 to B-4-25)"
    input_type = WatershedInput
    output_type = WatershedOutput

    def run(self, inp: WatershedInput) -> WatershedOutput:
        D, ws, wb, sl, ci = defaults(), inp.watershed, inp.water_body, inp.soil, inp.chemical
        notes, warnings = [], []
        mercury = ci.cas == MERCURY_CAS
        g = lambda v, key: _d(key) if v is None else v  # noqa: E731
        zs = ws.soil_mixing_depth_cm if ws.soil_mixing_depth_cm is not None \
            else D["soil_mixing_depth_cm"]["untilled"]
        bd = g(ws.bulk_density_g_cm3, "soil_bulk_density_g_cm3")
        tsw = g(ws.soil_water_content, "soil_water_content")
        rhos = g(ws.particle_density_g_cm3, "soil_particle_density_g_cm3")
        ta = g(ws.air_temperature_k, "air_temperature_k")
        r = _d("gas_constant_atm_m3_mol_k")
        irr = g(ws.irrigation_cm_yr, "irrigation_cm_yr")
        al, ai = ws.area_m2, ws.impervious_area_m2
        # Erosion: Xe (Table B-4-13), SD (Table B-4-14)
        if ws.usle is not None:
            u = ws.usle
            U = D["usle"]
            pick = lambda v, k: U[k]["value"] if v is None else v  # noqa: E731
            xe = eq.unit_soil_loss(u.rainfall_rf_per_yr, pick(u.erodibility_k_ton_acre, "erodibility_k_ton_acre"),
                                   pick(u.length_slope_ls, "length_slope_ls"),
                                   pick(u.cover_management_c, "cover_management_c"),
                                   pick(u.supporting_practice_pf, "supporting_practice_pf"))
        else:
            xe = ws.unit_soil_loss_kg_m2_yr
        sdd = D["sediment_delivery"]
        a_sd = None
        if ws.sediment_delivery_ratio is not None:
            sd = ws.sediment_delivery_ratio
        else:
            a_sd = eq.sd_intercept(al, sdd["a_by_area_sq_mi"], sdd["m2_per_sq_mi"])
            sd = eq.sediment_delivery_ratio(a_sd, al, ws.sd_b or sdd["b"])
        # Chemical parameters; mercury: species in soil, mercuric chloride in the water body
        if mercury:
            sp = D["mercury"]
            species = {s: _params(ChemicalInput(), c, inp.chemical_table) for s, c in sp["species_cas"].items()}
            wpar = _params(ci, sp["water_body_parameters_cas"], inp.chemical_table)
            for s in species:
                species[s]["fv"] = sp["vapor_fraction"]["value"]
            wpar["fv"] = sp["vapor_fraction"]["value"]
            organic = False
            qfac = sp["emission_fraction"]["value"]
            split = {s: sp["soil_split"][s] for s in species}
            notes.append("Mercury: 0.48 Q, Fv = 0.85; soil 98 % Hg2+ / 2 % MHg; water body with mercuric chloride "
                         "parameters; Cdw 85 % Hg2+ / 15 % MHg (Tables B-4-1, B-4-8, B-4-9, B-4-12, B-4-24).")
        else:
            wpar = _params(ci, ci.cas, inp.chemical_table)
            species = {"COPC": wpar}
            organic = _organic(ci, ci.cas)
            qfac, split = 1.0, {"COPC": 1.0}
        er = ws.enrichment_ratio or D["enrichment_ratio"]["organic" if organic else "inorganic"]
        # Watershed soil (Tables B-4-1 to B-4-6)
        td = g(sl.deposition_period_yr, "deposition_period_yr")
        t1 = g(sl.t1_yr, "t1_yr")
        t2 = g(sl.exposure_duration_yr, "exposure_duration_yr")
        q = inp.air.emission_rate_g_s * qfac if inp.air else 0.0
        soil_res, cs_by_basis = [], {}
        for s, p in species.items():
            kse = eq.k_erosion(xe, sd, er, bd, zs, p["kds_ml_g"], tsw) if ws.erosion_loss_constant == "calculated" \
                else _d("soil_erosion_loss_constant")
            ksr = eq.k_runoff(ws.runoff_cm_yr, tsw, zs, p["kds_ml_g"], bd)
            ksl = eq.k_leaching(ws.precipitation_cm_yr, irr, ws.runoff_cm_yr, ws.evapotranspiration_cm_yr, tsw, zs, bd,
                                p["kds_ml_g"])
            ksv = eq.k_volatilization(p["henry_atm_m3_mol"], zs, p["kds_ml_g"], r, ta, bd, p["da_cm2_s"], rhos, tsw)
            ks = eq.k_soil(p["ksg_per_yr"], kse, ksr, ksl, ksv)
            if ksl < 0:
                warnings.append(f"{s}: P + I - RO - Ev < 0 gives a negative leaching constant ksl (Table B-4-5).")
            res = SoilResult(species=s, kse_per_yr=kse, ksr_per_yr=ksr, ksl_per_yr=ksl, ksv_per_yr=ksv,
                             ksg_per_yr=p["ksg_per_yr"], ks_per_yr=ks)
            if sl.concentration_mg_kg is not None:
                res.cs_measured_mg_kg = sl.concentration_mg_kg * split[s]
                cs_by_basis.setdefault("measured", {})[s] = res.cs_measured_mg_kg
            else:
                if ks <= 0:
                    raise ValueError(f"{s}: soil loss constant ks <= 0; Table B-4-1 needs ks > 0.")
                a = inp.air
                ds = eq.deposition_term(q * split[s], zs, bd, p["fv"], a.watershed_dytwv_s_m2_yr, a.watershed_dytwp_s_m2_yr)
                res.ds_mg_kg_yr = ds
                res.cs_td_mg_kg = eq.soil_conc_at_td(ds, ks, td)
                res.cs_average_mg_kg = eq.soil_conc_average(ds, ks, td, t1, t2)
                cs_by_basis.setdefault("carcinogens", {})[s] = res.cs_average_mg_kg
                cs_by_basis.setdefault("noncarcinogens", {})[s] = res.cs_td_mg_kg
            soil_res.append(res)
        if mercury and sl.concentration_mg_kg is not None:
            notes.append("Mercury, measured soil: total Hg split 98 % Hg2+ / 2 % MHg as in Table B-4-1.")
        # Water body transfer (Tables B-4-16 to B-4-22)
        dbs = g(wb.benthic_depth_m, "benthic_depth_m")
        tss = g(wb.tss_mg_l, "tss_mg_l")
        cbs = g(wb.bed_sediment_concentration_g_cm3, "bed_sediment_concentration_g_cm3")
        tbs = g(wb.bed_sediment_porosity, "bed_sediment_porosity")
        twk = g(wb.temperature_k, "water_temperature_k")
        theta = g(wb.theta, "temperature_correction_theta")
        cd = g(wb.drag_coefficient, "drag_coefficient")
        w = g(wb.wind_speed_m_s, "wind_speed_m_s")
        vk = g(wb.von_karman, "von_karman")
        lz = g(wb.viscous_sublayer, "viscous_sublayer")
        rho_a = g(wb.air_density_g_cm3, "air_density_g_cm3")
        rho_w = g(wb.water_density_g_cm3, "water_density_g_cm3")
        mu_a = g(wb.air_viscosity_g_cm_s, "air_viscosity_g_cm_s")
        mu_w = g(wb.water_viscosity_g_cm_s, "water_viscosity_g_cm_s")
        dwc, aw, vfx = wb.water_column_depth_m, wb.area_m2, wb.flow_m3_yr
        dz = dwc + dbs
        if wb.type == "river":
            kl = eq.kl_river(wpar["dw_cm2_s"], wb.current_velocity_m_s, dz)
            kg = _d("gas_transfer_river_m_yr")
        else:
            kl = eq.kl_lake(cd, w, rho_a, rho_w, vk, lz, mu_w, wpar["dw_cm2_s"])
            kg = eq.kg_lake(cd, w, vk, lz, mu_a, rho_a, wpar["da_cm2_s"])
        h = wpar["henry_atm_m3_mol"]
        kv_m = eq.overall_transfer(kl, kg, h, r, twk, theta) if h > 0 else 0.0
        if h <= 0:
            notes.append("H = 0: no volatilization from the water body and no vapor diffusion load.")
        fwc = eq.fraction_water_column(wpar["kdsw_l_kg"], tss, dwc, dbs, tbs, wpar["kdbs_l_kg"], cbs)
        kv = eq.k_volatilization_water(kv_m, dz, wpar["kdsw_l_kg"], tss)
        kb_calc = eq.k_burial(xe, al, sd, vfx, tss, aw, cbs, dbs)
        kb = kb_calc
        if kb_calc < 0:
            kb = _d("burial_rate_negative")
            warnings.append(f"kb = {kb_calc:.4g} 1/yr < 0 (flow large compared with the watershed); kb = 0 used "
                            "(HHRAP Chapter 5, p. 5-81).")
        elif kb_calc > 1:
            warnings.append(f"kb = {kb_calc:.4g} 1/yr > 1: re-evaluate the unit soil loss Xe (HHRAP Chapter 5, p. 5-81).")
        kwt = eq.k_water_total(fwc, kv, kb)
        transfer = Transfer(unit_soil_loss_kg_m2_yr=xe, sd_intercept_a=a_sd, sediment_delivery_ratio=sd, kl_m_yr=kl,
                            kg_m_yr=kg, kv_overall_m_yr=kv_m, kv_per_yr=kv, kb_calculated_per_yr=kb_calc, kb_per_yr=kb,
                            fwc=fwc, fbs=1.0 - fwc, kwt_per_yr=kwt)
        # Loads (Tables B-4-7 to B-4-12)
        a = inp.air
        if a is not None:
            fvw = wpar["fv"]
            ldep = eq.load_deposition(q, fvw, a.water_body_dytwv_s_m2_yr, a.water_body_dytwp_s_m2_yr, aw)
            lri = eq.load_impervious(q, fvw, a.watershed_dytwv_s_m2_yr, a.watershed_dytwp_s_m2_yr, ai)
            ldif = eq.load_diffusion(kv_m, q, fvw, a.water_body_cywv_ug_s_g_m3, aw, h, r, twk) if h > 0 else 0.0
        else:
            ldep = lri = ldif = 0.0
        lext = sum(e.load_g_yr for e in inp.external_loads)
        results = []
        for basis, cs in cs_by_basis.items():
            lr = sum(eq.load_runoff(ws.runoff_cm_yr, al, ai, cs[s], bd, tsw, species[s]["kds_ml_g"]) for s in cs)
            le = sum(eq.load_erosion(xe, al, ai, sd, er, cs[s], species[s]["kds_ml_g"], bd, tsw) for s in cs)
            lt = ldep + ldif + lri + lr + le + lext              # Table B-4-7, p. B-199 [202] + external loads
            cwtot = eq.total_water_body(lt, vfx, fwc, kwt, aw, dwc, dbs)
            cwctot = eq.total_water_column(fwc, cwtot, dwc, dbs)
            cdw = eq.dissolved(cwctot, wpar["kdsw_l_kg"], tss)
            csb = eq.bed_sediment(fwc, cwtot, wpar["kdbs_l_kg"], tbs, cbs, dwc, dbs)
            spl = None
            if mercury:
                ws_ = D["mercury"]["water_split"]
                spl = {"Hg2+": ws_["Hg2+"] * cdw, "MHg": ws_["MHg"] * cdw}
            results.append(Concentrations(
                basis=basis, soil_mg_kg=cs, loads=Loads(deposition_g_yr=ldep, diffusion_g_yr=ldif,
                                                        impervious_runoff_g_yr=lri, pervious_runoff_g_yr=lr,
                                                        erosion_g_yr=le, external_g_yr=lext, total_g_yr=lt),
                total_water_body_mg_l=cwtot, total_water_column_mg_l=cwctot, dissolved_mg_l=cdw,
                bed_sediment_mg_kg=csb, dissolved_species_mg_l=spl))
        if lext:
            notes.append(f"External loads ({lext:.4g} g/yr) added to LT (ORTAM extension of Table B-4-7).")
        name = ci.name or (hhrap_params().get(ci.cas, {}).get("name") if ci.cas else None) or ci.cas or "COPC"
        params = {s: p for s, p in species.items()}
        if mercury:
            params["water body (mercuric chloride)"] = wpar
        return WatershedOutput(chemical=name, mercury=mercury, parameters=params, soil=soil_res, transfer=transfer,
                               results=results, notes=notes, warnings=warnings)

    def summary(self, out: WatershedOutput) -> str:
        t = out.transfer
        lines = [f"surfacewater.watershed: {out.chemical}",
                 f"Xe = {t.unit_soil_loss_kg_m2_yr:.4g} kg/m2-yr, SD = {t.sediment_delivery_ratio:.4g}, "
                 f"Kv = {t.kv_overall_m_yr:.4g} m/yr, kv = {t.kv_per_yr:.4g} 1/yr, kb = {t.kb_per_yr:.4g} 1/yr, "
                 f"fwc = {t.fwc:.4g}, kwt = {t.kwt_per_yr:.4g} 1/yr"]
        for s in out.soil:
            lines.append(f"Soil {s.species}: ks = {s.ks_per_yr:.4g} 1/yr" +
                         (f", Cs = {s.cs_average_mg_kg:.4g}, CstD = {s.cs_td_mg_kg:.4g} mg/kg" if s.ds_mg_kg_yr is not None
                          else f", Cs = {s.cs_measured_mg_kg:.4g} mg/kg (measured)"))
        for r in out.results:
            L = r.loads
            lines.append(f"[{r.basis}] LT = {L.total_g_yr:.4g} g/yr (dep {L.deposition_g_yr:.3g}, dif {L.diffusion_g_yr:.3g}, "
                         f"imp {L.impervious_runoff_g_yr:.3g}, runoff {L.pervious_runoff_g_yr:.3g}, erosion "
                         f"{L.erosion_g_yr:.3g}, external {L.external_g_yr:.3g})")
            lines.append(f"[{r.basis}] Cwtot = {r.total_water_body_mg_l:.4g} mg/L, Cwctot = {r.total_water_column_mg_l:.4g}"
                         f" mg/L, Cdw = {r.dissolved_mg_l:.4g} mg/L, Csb = {r.bed_sediment_mg_kg:.4g} mg/kg")
            if r.dissolved_species_mg_l:
                lines.append(f"[{r.basis}] Cdw Hg2+ = {r.dissolved_species_mg_l['Hg2+']:.4g}, MHg = "
                             f"{r.dissolved_species_mg_l['MHg']:.4g} mg/L")
        lines += [f"Warning: {w}" for w in out.warnings] + [f"Note: {n}" for n in out.notes]
        return "\n".join(lines)
