"""Particulate emission factors (PEF): wind erosion and construction activities.

  wind                 PEF   (chronic; RSL UG 4.9.1, p. 196; SSG Eq. 4-5)
  construction-traffic PEFsc (unpaved road traffic; SSG Eq. 5-5, 5-6; RSL UG 4.9.2)
  construction-other   PEF'sc (wind, excavation, dozing, grading, tilling; RSL UG 4.9.3)
Equations: modules/air/equations.py. Defaults: air_emission_factors.yaml.

PEFsc has two published forms: SSG 2002 Eq. 5-5 (constant 556; the source
cited by the RSL parameter workbook; default here) and the RSL User's Guide
4.9.2 form with silt content s and moisture M_dry (no defaults published;
`form: rsl-ug` requires both).
"""

from __future__ import annotations

from typing import Literal

from pydantic import Field, model_validator

from modules.air import equations as eq
from modules.air.common import defaults
from ortam_core.model import ModelInput, ModelOutput, Module


class WindInput(ModelInput):
    q_over_c: float | None = Field(None, gt=0, description="Q/Cwind, g/m2-s per kg/m3 (default 93.77)")
    vegetative_cover: float | None = Field(None, ge=0, lt=1, description="V, - (default 0.5)")
    um_m_s: float | None = Field(None, gt=0, description="Mean annual wind speed, m/s (default 4.69)")
    ut_m_s: float | None = Field(None, gt=0, description="Threshold wind speed at 7 m, m/s (default 11.32)")


class TrafficInput(ModelInput):
    road_length_ft: float = Field(gt=0, description="Length of contaminated road segment LR, ft")
    road_width_ft: float | None = Field(None, gt=0, description="WR, ft (default 20)")
    mean_vehicle_weight_t: float = Field(gt=0, description="W, tons")
    rain_days_per_yr: float = Field(ge=0, lt=365, description="p, days/yr with >= 0.01 in precipitation")
    vkt_km: float = Field(gt=0, description="Sum of fleet vehicle kilometers traveled, km")
    area_acres: float | None = Field(None, ge=0.5, le=500, description="Site area for Q/Csr, acres (default 0.5)")
    form: Literal["ssg-2002", "rsl-ug"] = "ssg-2002"
    silt_pct: float | None = Field(None, gt=0, description="Road surface silt content s, % (rsl-ug form)")
    moisture_dry_pct: float | None = Field(None, gt=0, description="M_dry, % (rsl-ug form)")

    @model_validator(mode="after")
    def _rsl_needs_s_m(self):
        if self.form == "rsl-ug" and (self.silt_pct is None or self.moisture_dry_pct is None):
            raise ValueError("form 'rsl-ug' requires silt_pct and moisture_dry_pct (no published defaults).")
        return self


class OtherInput(ModelInput):
    area_acres: float = Field(ge=0.5, le=500, description="Area of surface soil contamination, acres")
    excavation_area_m2: float = Field(0.0, ge=0, description="A_excav, m2")
    excavation_depth_m: float = Field(0.0, ge=0, description="d_excav, m")
    dozing_area_acres: float = Field(0.0, ge=0)
    dozer_blade_m: float = Field(1.0, gt=0)
    dozing_passes: float = Field(0.0, ge=0)
    grading_area_acres: float = Field(0.0, ge=0)
    grader_blade_m: float = Field(1.0, gt=0)
    grading_passes: float = Field(0.0, ge=0)
    tilling_area_acres: float = Field(0.0, ge=0)
    tilling_passes: float | None = Field(None, ge=0, description="NA_till (default 2)")


class SslPefInput(ModelInput):
    method: Literal["wind", "construction-traffic", "construction-other"] = "wind"
    wind: WindInput = WindInput()
    traffic: TrafficInput | None = None
    other: OtherInput | None = None


class SslPefOutput(ModelOutput):
    method: str
    pef_m3_kg: float
    q_over_c: float
    dispersion_correction: float | None = None
    details: dict[str, float] = {}


class SslPefModule(Module):
    name = "air.ssl-pef"
    title = "Air: particulate emission factor (PEF, PEFsc, PEF'sc; EPA SSL)"
    input_type = SslPefInput
    output_type = SslPefOutput

    def run(self, inputs: SslPefInput) -> SslPefOutput:
        d = defaults()
        con = d["construction"]
        tc_h = con["ed_yr"] * con["ew_wk_yr"] * 7.0 * 24.0
        fd = eq.dispersion_correction(tc_h)
        t_s = con["ed_yr"] * con["ef_d_yr"] * con["et_h_d"] * 3600.0   # Tt, RSL UG 4.9.2, p. 200

        if inputs.method == "wind":
            w, pw = inputs.wind, d["pef_wind"]
            qc = w.q_over_c or pw["q_over_c"]
            um, ut = w.um_m_s or pw["um_m_s"], w.ut_m_s or pw["ut_m_s"]
            v = w.vegetative_cover if w.vegetative_cover is not None else pw["v"]
            # F(x): workbook value 0.194 with the default wind speeds (reproduces the RSL PEF of
            # 1.36E+09 m3/kg); the RSL UG 4.9.4 polynomial (0.1931 at the defaults) when Um or Ut is site-specific.
            fx = pw["f_x"] if (w.um_m_s is None and w.ut_m_s is None) else eq.f_of_x(0.886 * ut / um)
            pef = eq.pef_wind(qc, v, um, ut, fx)
            return SslPefOutput(method="wind", pef_m3_kg=pef, q_over_c=qc, details={"f_x": fx})

        if inputs.method == "construction-traffic":
            tr, pt = inputs.traffic, d["pef_traffic"]
            if tr is None:
                raise ValueError("method 'construction-traffic' needs a 'traffic' block.")
            qc = eq.q_over_c(pt["a"], pt["b"], pt["c"], tr.area_acres or 0.5)
            ar = tr.road_length_ft * (tr.road_width_ft or pt["road_width_ft"]) * 0.092903
            if tr.form == "ssg-2002":
                pef = eq.pef_traffic_ssg(qc, fd, t_s, ar, tr.mean_vehicle_weight_t, tr.rain_days_per_yr, tr.vkt_km)
            else:
                pef = eq.pef_traffic_rsl(qc, fd, t_s, ar, tr.mean_vehicle_weight_t, tr.rain_days_per_yr,
                                         tr.vkt_km, tr.silt_pct, tr.moisture_dry_pct)
            return SslPefOutput(method=inputs.method, pef_m3_kg=pef, q_over_c=qc, dispersion_correction=fd,
                                details={"road_area_m2": ar, "t_s": t_s})

        o, po = inputs.other, d["pef_other"]
        if o is None:
            raise ValueError("method 'construction-other' needs an 'other' block.")
        qc = eq.q_over_c(po["a"], po["b"], po["c"], o.area_acres)
        asurf = o.area_acres * eq.ACRE_M2
        fx = po["f_x"]   # workbook value with the default wind speeds (see wind method)
        m = {
            "m_wind_g": eq.m_wind(po["v"], po["um_m_s"], po["ut_m_s"], fx, asurf, con["ed_yr"]),
            "m_excav_g": eq.m_excavation(po["um_m_s"], po["moisture_excav_pct"], po["rho_soil_mg_m3"],
                                         o.excavation_area_m2, o.excavation_depth_m, po["n_dump"]),
            "m_doz_g": eq.m_dozing(po["silt_doz_pct"], po["moisture_doz_pct"],
                                   eq.vkt_blade(o.dozing_area_acres, o.dozer_blade_m, o.dozing_passes),
                                   po["speed_doz_kph"]),
            "m_grade_g": eq.m_grading(po["speed_grade_kph"],
                                      eq.vkt_blade(o.grading_area_acres, o.grader_blade_m, o.grading_passes)),
            "m_till_g": eq.m_tilling(po["silt_till_pct"], o.tilling_area_acres,
                                     o.tilling_passes if o.tilling_passes is not None else po["n_till"]),
        }
        jt = sum(m.values()) / (asurf * t_s)          # <J'T>, g/m2-s (RSL UG 4.9.3, p. 202)
        pef = eq.pef_other(qc, fd, jt)
        return SslPefOutput(method=inputs.method, pef_m3_kg=pef, q_over_c=qc, dispersion_correction=fd,
                            details={**m, "j_t_g_m2_s": jt, "f_x": fx})

    def summary(self, o: SslPefOutput) -> str:
        lines = [f"Particulate emission factor ({o.method})", f"PEF        : {o.pef_m3_kg:.4g} m3/kg",
                 f"Q/C        : {o.q_over_c:.4g} g/m2-s per kg/m3"]
        if o.dispersion_correction:
            lines.append(f"FD         : {o.dispersion_correction:.5f}")
        lines += [f"{k:10s} : {v:.4g}" for k, v in o.details.items()]
        return "\n".join(lines)
