"""Johnson & Ettinger vapor intrusion model (EPA implementation v6.0, 2017).

Computes the steady-state attenuation coefficient alpha and the indoor air
concentration Cia = alpha * CS from a groundwater, exterior soil gas or sub-slab
soil gas source, for the five residential foundation types of EPA's spreadsheet,
1-3 soil strata (12 US SCS soil texture classes) and, for groundwater, the
capillary zone above the water table.

Data:      data/reference/je/v6-2017/ (imported from EPA's spreadsheet, modules/vapor/je_import.py)
Equations: modules/vapor/je_equations.py (J&E v6 equation numbers and pages)
Checked:   against EPA's spreadsheet, 16 cases (modules/vapor/tests/test_johnson_ettinger.py)

Errors and warnings follow the checks of the spreadsheet (MODEL!J11-J85): a
geometry that the spreadsheet marks "ERROR" is rejected; a value outside the
spreadsheet's recommended range gives a warning.
Indoor air risk is not computed here; workflows/vapor_intrusion passes Cia to
the risk air tools.
"""

from __future__ import annotations

import csv
import math
from functools import lru_cache
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field

from modules.vapor import je_equations as je
from ortam_core.model import ModelInput, ModelOutput, Module

DATA = Path(__file__).resolve().parents[2] / "data" / "reference" / "je" / "v6-2017"

Foundation = Literal["Slab-on-grade", "Basement w/ slab", "Basement w/ dirt floor",
                     "Closed crawl space w/ slab", "Closed crawl space w/ dirt floor"]
SoilName = Literal["Clay", "Clay Loam", "Loam", "Loamy Sand", "Sand", "Sandy Clay", "Sandy Clay Loam",
                   "Sandy Loam", "Silt", "Silt Loam", "Silty Clay", "Silty Clay Loam"]
DIRT = ("Basement w/ dirt floor", "Closed crawl space w/ dirt floor")


@lru_cache
def _table(name: str) -> list[dict]:
    with (DATA / name).open(encoding="utf-8") as f:
        return [{k: (float(v) if v not in ("",) and k not in ("name", "cas", "soil", "scs_symbol", "foundation")
                     else (v if v != "" else None)) for k, v in row.items()} for row in csv.DictReader(f)]


def _chemical(cas: str | None, name: str | None) -> dict:
    for row in _table("chemicals.csv"):
        if (cas and row["cas"] == cas) or (not cas and name and row["name"].lower() == name.lower()):
            return dict(row)
    raise KeyError(f"Chemical not in the J&E v6 table: {cas or name}")


class Stratum(ModelInput):
    soil: SoilName
    thickness_m: float = Field(gt=0, description="Stratum thickness, m")
    total_porosity: float | None = Field(None, gt=0, lt=1, description="n, cm3/cm3 (default: SCS class)")
    water_filled_porosity: float | None = Field(None, ge=0, lt=1, description="nw, cm3/cm3 (default: SCS class)")


class Building(ModelInput):
    lb_m: float | None = Field(None, gt=0, description="Depth below grade to base of foundation, m")
    lf_m: float | None = Field(None, ge=0, description="Foundation thickness, m")
    eta: float | None = Field(None, gt=0, le=1, description="Fraction of foundation area with cracks, -")
    abf_m2: float | None = Field(None, gt=0, description="Enclosed space floor area, m2")
    hb_m: float | None = Field(None, gt=0, description="Enclosed space mixing height, m")
    ach_per_h: float | None = Field(None, gt=0, description="Indoor air exchange rate, 1/h")
    qsoil_qb: float | None = Field(None, gt=0, lt=1, description="Qsoil/Qb, - (default: Eq. 17)")


class CapillaryZone(ModelInput):
    height_m: float | None = Field(None, gt=0, description="hcz, m (default: Eq. 7b)")
    total_porosity: float | None = Field(None, gt=0, lt=1, description="ncz, cm3/cm3")
    water_filled_porosity: float | None = Field(None, ge=0, lt=1, description="nw,cz, cm3/cm3 (default: Eq. 4b)")


class ChemicalInput(ModelInput):
    cas: str | None = None
    name: str | None = None
    concentration: float = Field(gt=0, description="Measured concentration: ug/L (groundwater) or ug/m3 (soil gas). "
                                                   "Non-detects are not entered as 0 or DL/2.")
    henry_atm_m3_mol: float | None = Field(None, gt=0, description="Hc at 25 C, atm-m3/mol")
    dair_cm2_s: float | None = Field(None, gt=0, description="Diffusivity in air, cm2/s")
    dwater_cm2_s: float | None = Field(None, gt=0, description="Diffusivity in water, cm2/s")
    tboil_k: float | None = Field(None, gt=0, description="Normal boiling point, K")
    tcrit_k: float | None = Field(None, gt=0, description="Critical temperature, K")
    dhvb_cal_mol: float | None = Field(None, gt=0, description="Enthalpy of vaporization at Tb, cal/mol")
    solubility_mg_l: float | None = Field(None, gt=0, description="Water solubility, mg/L")
    vp_mmhg: float | None = Field(None, gt=0, description="Vapor pressure, mmHg")
    mw_g_mol: float | None = Field(None, gt=0, description="Molecular weight, g/mol")


class JeInput(ModelInput):
    source: Literal["groundwater", "exterior-soil-gas", "subslab-soil-gas"]
    depth_to_source_m: float = Field(gt=0, description="Ls: depth below grade to water table or soil gas sample, m")
    temperature_c: float = Field(25.0, description="Average soil/groundwater temperature Ts, C")
    foundation: Foundation = "Slab-on-grade"
    building: Building = Building()
    strata: list[Stratum] = Field(min_length=1, max_length=3, description="Strata A, B, C from the surface down")
    capillary_zone: CapillaryZone = CapillaryZone()
    chemicals: list[ChemicalInput] = Field(min_length=1)


class JeResult(BaseModel):
    name: str
    cas: str
    concentration: float
    hs: float = Field(description="H'S at Ts, -")
    cs_ug_m3: float
    cs_fraction_of_saturation: float | None = None
    deff_strata_cm2_s: list[float]
    deff_cz_cm2_s: float | None = None
    deff_t_cm2_s: float
    a_param: float
    b_param: float
    c_param: float
    alpha: float
    cia_ug_m3: float | None
    css_ug_m3: float | None
    valid: bool = True
    notes: list[str] = []


class JeOutput(ModelOutput):
    source: str
    foundation: str
    lb_m: float
    lf_m: float
    eta: float
    abf_m2: float
    hb_m: float
    ach_per_h: float
    qsoil_qb: float
    qb_m3_h: float
    qsoil_m3_h: float | None
    capillary_zone_height_m: float | None = None
    capillary_zone_total_porosity: float | None = None
    capillary_zone_water_filled_porosity: float | None = None
    chemicals: list[JeResult]


class JohnsonEttingerModule(Module):
    name = "vapor.johnson-ettinger"
    title = "Vapor intrusion: Johnson & Ettinger model (EPA v6.0, 2017)"
    input_type = JeInput
    output_type = JeOutput

    def run(self, inp: JeInput) -> JeOutput:
        warnings: list[str] = []
        gw = inp.source == "groundwater"
        subslab = inp.source == "subslab-soil-gas"
        dirt = inp.foundation in DIRT
        ls, ts = inp.depth_to_source_m, inp.temperature_c

        # Building (BLDG_DATA residential defaults, user values override)
        bdef = next(b for b in _table("buildings.csv") if b["foundation"] == inp.foundation)
        bu = inp.building
        lb = bu.lb_m if bu.lb_m is not None else bdef["lb_m"]
        lf = bu.lf_m if bu.lf_m is not None else bdef["lf_m"]
        eta = bu.eta if bu.eta is not None else bdef["eta"]
        abf = bu.abf_m2 if bu.abf_m2 is not None else bdef["abf_m2"]
        hb = bu.hb_m if bu.hb_m is not None else bdef["hb_m"]
        ach = bu.ach_per_h if bu.ach_per_h is not None else bdef["ach_per_h"]
        if bu.qsoil_qb is not None:
            qsoil_qb = bu.qsoil_qb
        else:
            qsoil_qb = je.qsoil_qb_adjusted(ach)
            if ach != bdef["ach_per_h"]:
                warnings.append(f"Qsoil/Qb adjusted for ach = {ach:g} 1/h: {qsoil_qb:.4g} (J&E v6 Eq. 17, p. 29).")
        qb = je.ventilation_rate(abf, hb, ach)
        qsoil = None if dirt else qsoil_qb * qb
        if dirt and lf != 0:
            raise ValueError("Dirt floor: foundation thickness Lf must be 0 (MODEL!J39).")

        # Strata
        soils = {s["soil"]: s for s in _table("soils.csv")}
        n = [s.total_porosity if s.total_porosity is not None else soils[s.soil]["n_total_porosity"]
             for s in inp.strata]
        nw = [s.water_filled_porosity if s.water_filled_porosity is not None
              else soils[s.soil]["nw_water_filled_porosity"] for s in inp.strata]
        h = [s.thickness_m for s in inp.strata]
        if any(w >= p for w, p in zip(nw, n)):
            raise ValueError("Water-filled porosity must be smaller than total porosity in every stratum.")
        if not math.isclose(sum(h), ls, rel_tol=1e-9):
            raise ValueError(f"Sum of stratum thicknesses ({sum(h):g} m) must equal the depth to the source "
                             f"Ls = {ls:g} m (MODEL!J62, J68, J72).")

        # Capillary zone (groundwater only), soil of the deepest stratum (MODEL!F73-F76)
        hcz = ncz = nwcz = None
        if gw:
            deep = soils[inp.strata[-1].soil]
            cz = inp.capillary_zone
            hcz = cz.height_m if cz.height_m is not None else je.capillary_zone_height_m(deep["mean_grain_diameter_cm"])
            ncz = cz.total_porosity if cz.total_porosity is not None else deep["n_total_porosity"]
            nwcz = cz.water_filled_porosity if cz.water_filled_porosity is not None else \
                je.capillary_water_filled_porosity(deep["n_total_porosity"], deep["theta_r"], deep["vg_m"])

        # Geometry checks (MODEL!J12, J38, J56, J62, J68)
        if len(h) == 1 and gw and h[0] < lb + hcz:
            raise ValueError(f"Stratum A ({h[0]:g} m) is thinner than Lb + capillary zone height "
                             f"({lb:g} + {hcz:.3g} m) (MODEL!J56).")
        if h[0] < lb:
            raise ValueError(f"Stratum A ({h[0]:g} m) must extend below the foundation base Lb = {lb:g} m (MODEL!J56).")
        if gw and len(h) == 2 and h[1] < hcz:
            warnings.append(f"Stratum B ({h[1]:g} m) is thinner than the capillary zone ({hcz:.3g} m) (MODEL!J62).")
        if gw and len(h) == 3 and h[2] < hcz:
            raise ValueError(f"Stratum C ({h[2]:g} m) is thinner than the capillary zone ({hcz:.3g} m) (MODEL!J68).")
        if ls <= lb:
            raise ValueError("The source must lie below the base of the foundation (Ls > Lb).")
        if subslab and not (lb <= ls <= lb + 1):
            warnings.append("Sub-slab sample depth should lie between Lb and Lb + 1 m (MODEL!J12).")
        elif gw and ls < max(lb + hcz, lb + 1):
            warnings.append("Ls < max(Lb + hcz, Lb + 1 m): source very close to the foundation (MODEL!J12).")
        elif inp.source == "exterior-soil-gas" and ls < lb + 1:
            warnings.append("Exterior soil gas sample less than 1 m below the foundation (MODEL!J12).")
        t_lo, t_hi = (3, 25) if gw else (3, 30)
        if not t_lo <= ts <= t_hi:
            warnings.append(f"Ts = {ts:g} C outside the spreadsheet range {t_lo}-{t_hi} C (MODEL!J13).")

        results = []
        for ci in inp.chemicals:
            c = _chemical(ci.cas, ci.name)
            for k in ("henry_atm_m3_mol", "dair_cm2_s", "dwater_cm2_s", "tboil_k", "tcrit_k", "dhvb_cal_mol",
                      "solubility_mg_l", "vp_mmhg", "mw_g_mol"):
                if getattr(ci, k) is not None:
                    c[k] = getattr(ci, k)
            label = f"{c['name']} ({c['cas']})"
            notes: list[str] = []
            if c["henry_atm_m3_mol"] is None or c["dair_cm2_s"] is None or c["dwater_cm2_s"] is None:
                raise ValueError(f"{label}: Hc, Da and Dw are required.")
            hr = je.henry_reference(c["henry_atm_m3_mol"])
            if None in (c["tboil_k"], c["tcrit_k"], c["dhvb_cal_mol"]):
                hs = hr                    # IFERROR(..., Hr) in MODEL!F27
                if ts != 25:
                    notes.append("Tb, Tc or dHv,b missing: H' at 25 C used (MODEL!F27).")
            else:
                hs = je.henry_at_temperature(c["henry_atm_m3_mol"], ts, c["tboil_k"], c["tcrit_k"],
                                             c["dhvb_cal_mol"])
            cs = je.source_vapor_concentration(ci.concentration, hs, gw)
            deff = [je.effective_diffusivity(c["dair_cm2_s"], c["dwater_cm2_s"], ni, nwi, hs)
                    for ni, nwi in zip(n, nw)]
            deff_cz = je.effective_diffusivity(c["dair_cm2_s"], c["dwater_cm2_s"], ncz, nwcz, hs) if gw else None
            dt = je.total_effective_diffusivity(deff, h, lb, deff_cz, hcz)
            a = je.a_param(dt, abf, lb, qb, ls)
            b = je.b_param(qsoil_qb, qb, lf, deff[0], eta, abf, lb)
            if subslab:
                alpha = qsoil_qb             # J&E v6 p. 60; MODEL!G93
            else:
                alpha = je.attenuation_factor(a, b, qsoil_qb, dirt)
            valid = True
            if gw and c["solubility_mg_l"] is not None and ci.concentration > 1000.0 * c["solubility_mg_l"]:
                valid = False
                notes.append("Groundwater concentration exceeds the water solubility; result invalid (MODEL!J11).")
            frac = None
            if c["vp_mmhg"] is not None and c["mw_g_mol"] is not None:
                frac = cs / je.saturated_vapor_concentration(c["vp_mmhg"], c["mw_g_mol"])
                if frac > 1:
                    notes.append("CS exceeds the saturated vapor concentration (MODEL!F15).")
            cia = cs * alpha if valid else None                      # J&E v6 p. 60, Cbuilding = alpha CS
            css = None if (dirt or cia is None) else cia * qb / qsoil  # MODEL!G100
            results.append(JeResult(name=c["name"], cas=c["cas"], concentration=ci.concentration, hs=hs, cs_ug_m3=cs,
                                    cs_fraction_of_saturation=frac, deff_strata_cm2_s=deff, deff_cz_cm2_s=deff_cz,
                                    deff_t_cm2_s=dt, a_param=a, b_param=b, c_param=qsoil_qb, alpha=alpha,
                                    cia_ug_m3=cia, css_ug_m3=css, valid=valid, notes=notes))
        return JeOutput(source=inp.source, foundation=inp.foundation, lb_m=lb, lf_m=lf, eta=eta, abf_m2=abf, hb_m=hb,
                        ach_per_h=ach, qsoil_qb=qsoil_qb, qb_m3_h=qb, qsoil_m3_h=qsoil,
                        capillary_zone_height_m=hcz, capillary_zone_total_porosity=ncz,
                        capillary_zone_water_filled_porosity=nwcz, chemicals=results, warnings=warnings)

    def summary(self, o: JeOutput) -> str:
        unit = "ug/L" if o.source == "groundwater" else "ug/m3"
        lines = [f"Johnson & Ettinger (EPA v6.0): source {o.source}, {o.foundation}",
                 f"Qb = {o.qb_m3_h:.4g} m3/h, Qsoil/Qb = {o.qsoil_qb:.3g}"
                 + (f", hcz = {o.capillary_zone_height_m:.3g} m" if o.capillary_zone_height_m else ""),
                 f"{'Chemical':28s} {'C (' + unit + ')':>13s} {'alpha':>10s} {'Cia ug/m3':>11s} {'Css ug/m3':>11s}"]
        g = lambda v: f"{v:.3g}" if v is not None else "-"
        for c in o.chemicals:
            lines.append(f"{c.name[:28]:28s} {c.concentration:>13.4g} {c.alpha:>10.3g} {g(c.cia_ug_m3):>11s} "
                         f"{g(c.css_ug_m3):>11s}")
            lines += [f"  note: {x}" for x in c.notes]
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
