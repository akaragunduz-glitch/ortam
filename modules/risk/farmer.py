"""risk.farmer: RAIS farmer scenario (consumption of farm products; soil and biota).

Source: RAIS Chemical Risk Calculator User's Guide, section 4.8 (docs/references/
rais_chemical_risk_users_guide_2026-11.pdf; PDF pages, the guide is a saved web page):
  4.8.1  direct consumption of agricultural products (produce, dairy, beef concentrations given),
         PDF pp. 69-76
  4.8.3  direct exposure to soil and consumption of agricultural products (soil concentration),
         PDF pp. 94-108
Parameters: data/reference/rais/2026-11/farmer_defaults.yaml.

Food concentrations in the soil scenario (PDF pp. 94-95):
  produce = C_soil * (BV_wet + MLF_produce)
  dairy, beef: TF * [Qp*fp*fs*(BV_dry + MLF_pasture) + Qs*fp] * C_soil   (mg/kg FW)
CDI equations, all with C in mg/kg, 1e-6 kg/mg, AT = 365 d/yr:
  noncancer child   C*1e-6*EF*ED_c*CF*IR_c/(365*ED_c*BW_c)
  noncancer adult   C*1e-6*EF*ED_far*CF*IR_a/(365*ED_far*BW_a)
  age-adjusted      C*1e-6*CF*IF_adj/(365*ED_far),  IF_adj = EF*ED_c*IR_c/BW_c + EF*ED_a*IR_a/BW_a
  cancer            C*1e-6*CF*IF_adj/(365*LT); mutagens IF_M (ADAF segments); TCE CAF*IF + MAF*IF_M;
                    vinyl chloride IF/(365*LT) + IR_c/BW_c
  direct soil       ingestion (RBA, IRS), dermal (SA, AF, ABSd), inhalation of particulates and
                    vapors (1/VF + 1/PEF; cancer CDI in ug/m3 with 1000 ug/mg, noncancer in mg/m3)
Risk = CDI*SFO (dermal SFO/GIABS; inhalation IUR), one-hit above 0.01 per route and for the
total; HQ = CDI/RfD (dermal RfD*GIABS; inhalation RfC).
RAIS internal factor values that differ from its guide equations are kept and documented
in farmer_defaults.yaml (verified against the RAIS output, tests/test_farmer.py).
"""

from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, Field

from modules.risk.common import ChemicalInputBase, cancer_method, label, one_hit, resolve_chemical
from ortam_core.model import ModelInput, ModelOutput, Module

DEFAULTS = Path(__file__).resolve().parents[2] / "data" / "reference" / "rais" / "2026-11" / "farmer_defaults.yaml"
FOODS = ("produce", "beef", "dairy")
SOIL_ROUTES = ("ingestion", "dermal", "inhalation")


@lru_cache
def defaults() -> dict:
    return yaml.safe_load(DEFAULTS.read_text(encoding="utf-8"))


class ChemicalInput(ChemicalInputBase):
    soil_mg_kg: float | None = Field(None, ge=0, description="Soil concentration (soil scenario), mg/kg")
    produce_mg_kg: float | None = Field(None, ge=0, description="Produce concentration (food scenario), mg/kg FW")
    beef_mg_kg: float | None = Field(None, ge=0)
    dairy_mg_kg: float | None = Field(None, ge=0)
    bv_wet: float | None = Field(None, ge=0, description="Soil-to-plant factor, kg dry soil/kg fresh plant")
    bv_dry: float | None = Field(None, ge=0, description="Soil-to-plant factor, kg dry soil/kg dry plant")
    tf_beef_d_kg: float | None = Field(None, ge=0)
    tf_dairy_d_kg: float | None = Field(None, ge=0)
    vf_m3_kg: float | None = Field(None, gt=0, description="Volatilization factor (volatiles), m3/kg")
    rba: float | None = Field(None, gt=0, le=1)
    absd: float | None = Field(None, ge=0, le=1)


class FarmerInput(ModelInput):
    scenario: Literal["soil", "food"] = "soil"
    toxicity_table: Literal["rsl-2024-11", "none"] = "rsl-2024-11"
    chemicals: list[ChemicalInput] = Field(min_length=1)


class CdiSet(BaseModel):
    cancer: float | None = None
    child: float | None = None
    adult: float | None = None
    age_adjusted: float | None = None


class FarmerResult(BaseModel):
    name: str
    cas: str
    cancer_method: str
    food_mg_kg: dict[str, float] = Field(description="Produce, beef, dairy concentrations, mg/kg FW")
    cdi: dict[str, CdiSet] = Field(description="By route; mg/kg-day, inhalation cancer ug/m3, noncancer mg/m3")
    risk: dict[str, float] = {}
    total_risk: float | None = None
    hq: dict[str, dict[str, float]] = Field({}, description="child / adult / age_adjusted -> route -> HQ")
    hi: dict[str, float] = {}
    notes: list[str] = []


class FarmerOutput(ModelOutput):
    scenario: str
    chemicals: list[FarmerResult]


# ------------------------------------------------------------------ factors
def _if_adj(d, ir: dict) -> float:
    return (d["ef_d_yr"] * d["ed_child_yr"] * ir["child"] / d["bw_child_kg"]
            + d["ef_d_yr"] * d["ed_adult_yr"] * ir["adult"] / d["bw_adult_kg"])


def _if_mut(d, ir: dict, seg: dict) -> float:
    return sum(d["ef_d_yr"] * ed * ir[who] * adaf / bw
               for ed, bw, adaf, who in zip(seg["ed_yr"], seg["bw_kg"], seg["adaf"], seg.get("rates", seg.get("ages"))))


def _food_ir(d, food: str) -> dict:
    if food == "produce":
        return {k: d["ir_fruit"][k] + d["ir_vegetable"][k] for k in ("child", "adult")}
    return d[f"ir_{food}"]


def _oral_cdi(c_mg_kg: float, cf: float, ir: dict, if_adj: float, if_mut: float, d: dict, method: str,
              bw_child: float) -> CdiSet:
    """Food or soil ingestion CDIs (RAIS 4.8.1.1-4.8.1.7 and 4.8.3)."""
    ef, edf, lt = d["ef_d_yr"], d["ed_farmer_yr"], d["lt_yr"]
    k = c_mg_kg * 1e-6 * cf
    s = CdiSet(child=k * ef * d["ed_child_yr"] * ir["child"] / (365.0 * d["ed_child_yr"] * bw_child),
               adult=k * ef * edf * ir["adult"] / (365.0 * edf * d["bw_adult_kg"]),
               age_adjusted=k * if_adj / (365.0 * edf))
    if method == "mutagen":
        s.cancer = k * if_mut / (365.0 * lt)
    elif method == "trichloroethylene":
        s.cancer = k * (d["tce"]["caf_oral"] * if_adj + d["tce"]["maf_oral"] * if_mut) / (365.0 * lt)
    elif method == "vinyl chloride":
        # early-life term without averaging time (4.8.1.7, PDF p. 75; 4.8.3.7, p. 107)
        s.cancer = k * if_adj / (365.0 * lt) + k * ir["child"] / bw_child
    else:
        s.cancer = k * if_adj / (365.0 * lt)
    return s


class FarmerModule(Module):
    name = "risk.farmer"
    title = "Risk: RAIS farmer (farm products; soil and biota)"
    input_type = FarmerInput
    output_type = FarmerOutput

    def run(self, inp: FarmerInput) -> FarmerOutput:
        d = defaults()
        warnings, results = [], []
        for ci in inp.chemicals:
            c, overridden = resolve_chemical(inp.toxicity_table, ci)
            if overridden:
                warnings.append(f"{label(c)}: table values overridden by the user: {', '.join(overridden)}.")
            method = cancer_method(c)
            r = FarmerResult(name=c.name, cas=c.cas, cancer_method=method, food_mg_kg={}, cdi={})
            cf = d["contaminated_fraction"]
            # food concentrations
            if inp.scenario == "food":
                for f in FOODS:
                    v = getattr(ci, f"{f}_mg_kg")
                    if v is not None:
                        r.food_mg_kg[f] = v
            else:
                if ci.soil_mg_kg is None:
                    raise ValueError(f"{label(c)}: soil scenario needs soil_mg_kg.")
                cs = ci.soil_mg_kg
                if ci.bv_wet is None or ci.bv_dry is None or ci.tf_beef_d_kg is None or ci.tf_dairy_d_kg is None:
                    raise ValueError(f"{label(c)}: supply bv_wet, bv_dry, tf_beef_d_kg and tf_dairy_d_kg "
                                     "(RAIS transfer factors).")
                r.food_mg_kg["produce"] = cs * (ci.bv_wet + d["mlf_produce"])
                for f, tf in (("dairy", ci.tf_dairy_d_kg), ("beef", ci.tf_beef_d_kg)):
                    a = d["animals"][f]
                    r.food_mg_kg[f] = tf * (a["qp_kg_d"] * a["fp"] * a["fs"] * (ci.bv_dry + d["mlf_pasture"])
                                            + a["qs_kg_d"] * a["fp"]) * cs
            seg = d["food_segments"]
            for f, conc in r.food_mg_kg.items():
                ir = _food_ir(d, f)
                r.cdi[f] = _oral_cdi(conc, cf[f], ir, _if_adj(d, ir), _if_mut(d, ir, seg), d, method,
                                     d["bw_child_kg"])
            # direct soil contact
            if inp.scenario == "soil":
                cs = ci.soil_mg_kg
                rba = ci.rba if ci.rba is not None else 1.0
                ss = d["soil_mutagen_segments"]
                ifs = _if_adj(d, d["irs_mg_d"])
                ifsm = _if_mut(d, d["irs_mg_d"], ss)
                r.cdi["ingestion"] = _oral_cdi(cs * rba, 1.0, d["irs_mg_d"], ifs, ifsm, d, method, d["bw_child_kg"])
                absd = ci.absd if ci.absd is not None else c.absd
                if absd:
                    dfs_c = d["sa_cm2_d"]["child"] * d["af_mg_cm2"]["child"]
                    dfs_a = d["sa_cm2_d"]["adult"] * d["af_mg_cm2"]["adult"]
                    ir = {"child": dfs_c * absd, "adult": dfs_a * absd}
                    dfs = absd * (d["ef_d_yr"] * d["ed_child_yr"] * dfs_c / d["bw_child_kg"]
                                  + d["ef_d_yr"] * d["soil_dermal_adult_ed_yr"] * dfs_a / d["bw_adult_kg"])
                    dfsm = absd * _if_mut(d, {"child": dfs_c, "adult": dfs_a}, ss)
                    r.cdi["dermal"] = _oral_cdi(cs, 1.0, ir, dfs, dfsm, d, method, d["bw_child_kg"])
                inv = (1.0 / ci.vf_m3_kg if ci.vf_m3_kg else 0.0) + 1.0 / d["pef_m3_kg"]
                if c.volatile and not ci.vf_m3_kg:
                    r.notes.append("Volatile chemical without vf_m3_kg: particulates only.")
                ef, edf, et, lt = d["ef_d_yr"], d["ed_farmer_yr"], d["et_h_d"], d["lt_yr"]
                nc = cs * ef * edf * et / 24.0 * inv / (365.0 * edf)                       # mg/m3
                ist = d["inhalation_mutagen_segments"]
                mut = sum(ef * ed * et / 24.0 * adaf for ed, adaf in zip(ist["ed_yr"], ist["adaf"]))
                std = ef * edf * et / 24.0
                if method == "mutagen":
                    days = mut
                elif method == "trichloroethylene":
                    days = d["tce"]["caf_inh"] * std + d["tce"]["maf_inh"] * mut
                else:
                    days = std
                ca = cs * 1000.0 * inv * days / (365.0 * lt)                                # ug/m3
                if method == "vinyl chloride" and ci.vf_m3_kg:                           # 4.8.3.7, PDF p. 107
                    ca = cs * 1000.0 / ci.vf_m3_kg * (std / (365.0 * lt) + 1.0)
                r.cdi["inhalation"] = CdiSet(cancer=ca, child=nc, adult=nc, age_adjusted=nc)
            # risk and hazard
            lin = {}
            for route, s in r.cdi.items():
                if route == "inhalation":
                    tox_c, tox_n = c.iur_per_ug_m3, c.rfc_mg_m3
                elif route == "dermal":
                    tox_c = c.sfo_per_mg_kg_day / c.giabs if c.sfo_per_mg_kg_day else None
                    tox_n = c.rfd_mg_kg_day * c.giabs if c.rfd_mg_kg_day else None
                else:
                    tox_c, tox_n = c.sfo_per_mg_kg_day, c.rfd_mg_kg_day
                if tox_c and s.cancer is not None:
                    lin[route] = s.cancer * tox_c
                    r.risk[route] = one_hit(lin[route])
                if tox_n:
                    for who in ("child", "adult", "age_adjusted"):
                        r.hq.setdefault(who, {})[route] = getattr(s, who) / tox_n
            if lin:
                r.total_risk = one_hit(sum(lin.values()))
            r.hi = {who: sum(v.values()) for who, v in r.hq.items()}
            results.append(r)
        return FarmerOutput(scenario=inp.scenario, chemicals=results, warnings=warnings)

    def summary(self, o: FarmerOutput) -> str:
        g = lambda v, f: format(v, f) if v is not None else "-"
        lines = [f"RAIS farmer ({o.scenario})", f"{'Chemical':28s} {'Risk':>9s} {'HI child':>9s} {'HI adult':>9s}"]
        for c in o.chemicals:
            lines.append(f"{c.name[:28]:28s} {g(c.total_risk, '.2e'):>9s} {g(c.hi.get('child'), '.3g'):>9s} "
                         f"{g(c.hi.get('adult'), '.3g'):>9s}")
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
