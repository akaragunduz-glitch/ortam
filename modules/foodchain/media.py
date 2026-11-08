"""foodchain.media: concentrations in produce, animal products and fish (HHRAP 2005).

Inputs per chemical: soil concentration Cs (mg/kg, e.g. the exposure point
concentration), optionally the dissolved surface water concentration C_dw (mg/L) and
the bed sediment concentration C_sb (mg/kg) for fish, and optionally air-deposition
contributions (Pd + Pv, mg/kg DW) to exposed produce and forage/silage from an air model.
Transfer factors: HHRAP Companion Database (data/reference/hhrap/2005), or the
Appendix A-2 correlations (modules/foodchain/btf.py) from log Kow, H and Koc of the RSL
table for organic chemicals that are not in it; every value can be overridden.
Equations: modules/foodchain/equations.py (HHRAP Appendix B tables and pages).
Defaults: data/reference/hhrap/2005/defaults.yaml.

Fortran analogy: a chain of assignment statements; plant concentrations feed the
animal intake sums, which feed the tissue concentrations.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import yaml
from pydantic import BaseModel, Field

from modules.foodchain import btf
from modules.foodchain import equations as eq
from ortam_core.chemicals.database import lookup
from ortam_core.model import ModelInput, ModelOutput, Module

DATA = Path(__file__).resolve().parents[2] / "data" / "reference" / "hhrap" / "2005"
# HHRAP database parameter name -> ORTAM name
DB_KEYS = {"Brag": "br_ag", "Brf/s": "br_forage", "Brgrain": "br_grain", "Brroot": "br_rootveg",
           "Babeef": "ba_beef", "Bamilk": "ba_milk", "Bapork": "ba_pork", "Bapoult": "ba_chicken",
           "Baeggs": "ba_eggs", "BCFfish": "bcf_fish", "BAF": "baf_fish", "BSAFfish": "bsaf_fish",
           "Kow (log)": "log_kow"}
FACTORS = ["br_ag", "br_forage", "br_grain", "br_rootveg", "ba_beef", "ba_milk", "ba_pork", "ba_chicken",
           "ba_eggs", "bcf_fish", "baf_fish", "bsaf_fish"]


@lru_cache
def defaults() -> dict:
    return yaml.safe_load((DATA / "defaults.yaml").read_text(encoding="utf-8"))


@lru_cache
def hhrap_db() -> dict[str, dict]:
    out: dict[str, dict] = defaultdict(dict)
    with (DATA / "parameters.csv").open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            out[r["cas"]]["name"] = r["name"]
            if r["parameter"] in DB_KEYS and r["value"] != "":
                out[r["cas"]][DB_KEYS[r["parameter"]]] = float(r["value"])
    return dict(out)


class AirContribution(ModelInput):
    """Pd + Pv from an air model (HHRAP Tables B-2-7, B-2-8, B-3-7, B-3-8), mg/kg DW."""
    exposed_produce_mg_kg: float = Field(0.0, ge=0)
    forage_mg_kg: float = Field(0.0, ge=0)
    silage_mg_kg: float = Field(0.0, ge=0)


class ChemicalInput(ModelInput):
    cas: str | None = None
    name: str | None = None
    soil_mg_kg: float | None = Field(None, ge=0, description="Cs, mg/kg dry soil")
    surface_water_dissolved_mg_l: float | None = Field(None, ge=0, description="C_dw, mg/L")
    sediment_mg_kg: float | None = Field(None, ge=0, description="C_sb, mg/kg bed sediment")
    air: AirContribution = AirContribution()
    log_kow: float | None = None
    br_ag: float | None = Field(None, ge=0); br_forage: float | None = Field(None, ge=0)
    br_grain: float | None = Field(None, ge=0); br_rootveg: float | None = Field(None, ge=0)
    ba_beef: float | None = Field(None, ge=0); ba_milk: float | None = Field(None, ge=0)
    ba_pork: float | None = Field(None, ge=0); ba_chicken: float | None = Field(None, ge=0)
    ba_eggs: float | None = Field(None, ge=0)
    bcf_fish: float | None = Field(None, ge=0); baf_fish: float | None = Field(None, ge=0)
    bsaf_fish: float | None = Field(None, ge=0)
    metabolism_factor: float | None = Field(None, gt=0, le=1)


class MediaInput(ModelInput):
    kds_foc: float = Field(0.01, gt=0, le=1, description="foc for Kds = Koc*foc when Br_rootveg is estimated")
    fish_lipid_fraction: float | None = Field(None, gt=0, lt=1)
    sediment_oc_fraction: float | None = Field(None, gt=0, lt=1)
    chemicals: list[ChemicalInput] = Field(min_length=1)


class MediaResult(BaseModel):
    name: str
    cas: str
    factors: dict[str, float | None]
    factor_source: dict[str, str]
    produce_exposed_aboveground_mg_kg_dw: float | None = None
    produce_protected_aboveground_mg_kg_dw: float | None = None
    produce_belowground_mg_kg_dw: float | None = None
    forage_mg_kg_dw: float | None = None
    silage_mg_kg_dw: float | None = None
    grain_mg_kg_dw: float | None = None
    beef_mg_kg_fw: float | None = None
    milk_mg_kg_fw: float | None = None
    pork_mg_kg_fw: float | None = None
    eggs_mg_kg_fw: float | None = None
    chicken_mg_kg_fw: float | None = None
    fish_mg_kg_fw: float | None = None
    fish_method: str | None = None
    notes: list[str] = []


class MediaOutput(ModelOutput):
    chemicals: list[MediaResult]


def _factors(ci: ChemicalInput, foc: float) -> tuple[dict, dict, str, str]:
    db = hhrap_db()
    rsl = None
    cas = ci.cas
    if cas is None or cas not in db:
        try:
            rsl = lookup("rsl-2024-11", ci.cas, ci.name)
            cas = rsl.cas
        except Exception:
            rsl = None
    row = db.get(cas or "", {})
    name = row.get("name") or (rsl.name if rsl else (ci.name or cas or ""))
    f, src = {}, {}
    for k in FACTORS + ["log_kow"]:
        if getattr(ci, k, None) is not None:
            f[k], src[k] = getattr(ci, k), "user"
        elif k in row:
            f[k], src[k] = row[k], "HHRAP 2005 database"
    if f.get("log_kow") is None and rsl is not None and rsl.log_kow is not None:
        f["log_kow"], src["log_kow"] = rsl.log_kow, "RSL 2024-11"
    if not row and rsl is not None and rsl.log_kow is not None:     # organic chemical not in HHRAP
        lk = rsl.log_kow
        est = {"br_ag": btf.br_aboveground(lk), "br_forage": btf.br_aboveground(lk),
               "br_grain": btf.br_aboveground(lk)}
        if rsl.koc_l_kg:
            est["br_rootveg"] = btf.br_rootveg(btf.rcf_dw(lk), rsl.koc_l_kg * foc)
        for p in ("beef", "milk", "pork", "chicken", "eggs"):
            est[f"ba_{p}"] = btf.ba(lk, p)
        for k, v in est.items():
            if k not in f:
                f[k], src[k] = v, "HHRAP App. A-2 correlation (RSL log Kow, Koc)"
    return f, src, name, cas or ""


class FoodchainMediaModule(Module):
    name = "foodchain.media"
    title = "Food chain: produce, animal product and fish concentrations (HHRAP 2005)"
    input_type = MediaInput
    output_type = MediaOutput

    def run(self, inp: MediaInput) -> MediaOutput:
        d = defaults()
        an = d["animals"]
        F, Bs = an["fraction_contaminated_feed"]["value"], an["soil_bioavailability"]["value"]
        lipid = inp.fish_lipid_fraction or d["fish"]["lipid_fraction"]["value"]
        ocsed = inp.sediment_oc_fraction or d["fish"]["oc_sediment"]["value"]
        warnings, results = [], []
        for ci in inp.chemicals:
            f, src, name, cas = _factors(ci, inp.kds_foc)
            r = MediaResult(name=name, cas=cas, factors={k: f.get(k) for k in FACTORS + ["log_kow"]},
                            factor_source=src)
            if any(s.startswith("HHRAP App.") for s in src.values()):
                r.notes.append("Not in the HHRAP database: transfer factors from the Appendix A-2 correlations.")
            cs = ci.soil_mg_kg
            if cs is not None:
                need = ["br_ag", "br_forage", "br_grain", "br_rootveg"]
                missing = [k for k in need if f.get(k) is None]
                if missing:
                    raise ValueError(f"{name}: transfer factors missing ({', '.join(missing)}); supply them.")
                pr_ag = eq.produce_root_uptake(cs, f["br_ag"])
                r.produce_protected_aboveground_mg_kg_dw = pr_ag
                r.produce_exposed_aboveground_mg_kg_dw = pr_ag + ci.air.exposed_produce_mg_kg
                r.produce_belowground_mg_kg_dw = eq.produce_belowground(cs, f["br_rootveg"],
                                                                        eq.vg_rootveg(f.get("log_kow"), d))
                feeds = {"forage": eq.produce_root_uptake(cs, f["br_forage"]) + ci.air.forage_mg_kg,
                         "silage": eq.produce_root_uptake(cs, f["br_forage"]) + ci.air.silage_mg_kg,
                         "grain": eq.produce_root_uptake(cs, f["br_grain"])}
                r.forage_mg_kg_dw, r.silage_mg_kg_dw, r.grain_mg_kg_dw = feeds["forage"], feeds["silage"], feeds["grain"]
                mf_def = an["metabolism_factor"]
                mf = ci.metabolism_factor or mf_def["by_cas"].get(cas, mf_def["default"])
                for prod, attr in (("beef", "beef_mg_kg_fw"), ("milk", "milk_mg_kg_fw"), ("pork", "pork_mg_kg_fw"),
                                   ("eggs", "eggs_mg_kg_fw"), ("chicken", "chicken_mg_kg_fw")):
                    ba = f.get(f"ba_{prod}")
                    if ba is None:
                        r.notes.append(f"No Ba for {prod}.")
                        continue
                    a = an[prod]
                    intake = eq.animal_intake(feeds, a["qp_kg_dw_d"], F, a["qs_kg_d"], cs, Bs)
                    m = mf if prod in mf_def["applies_to"] else 1.0
                    setattr(r, attr, eq.animal_tissue(intake, ba, m))
            # fish (HHRAP App. A-2-2.13.4: BSAF for PCDD/F and PCBs, BAF for log Kow >= 4, BCF otherwise)
            if f.get("bsaf_fish") is not None and ci.sediment_mg_kg is not None:
                r.fish_mg_kg_fw = eq.fish_bsaf(ci.sediment_mg_kg, lipid, f["bsaf_fish"], ocsed)
                r.fish_method = "BSAF (Table B-4-28)"
            elif ci.surface_water_dissolved_mg_l is not None:
                if f.get("baf_fish") is not None:
                    r.fish_mg_kg_fw = eq.fish_baf(ci.surface_water_dissolved_mg_l, f["baf_fish"])
                    r.fish_method = "BAF (Table B-4-27)"
                elif f.get("bcf_fish") is not None:
                    r.fish_mg_kg_fw = eq.fish_bcf(ci.surface_water_dissolved_mg_l, f["bcf_fish"])
                    r.fish_method = "BCF (Table B-4-26)"
                else:
                    r.notes.append("No BCF/BAF: fish concentration not calculated.")
            elif f.get("bsaf_fish") is not None:
                r.notes.append("BSAF chemical: supply sediment_mg_kg for the fish concentration.")
            if cs is None and ci.surface_water_dissolved_mg_l is None and ci.sediment_mg_kg is None:
                warnings.append(f"{name}: no medium concentration given.")
            results.append(r)
        return MediaOutput(chemicals=results, warnings=warnings)

    def summary(self, o: MediaOutput) -> str:
        g = lambda v: f"{v:.3g}" if v is not None else "-"
        lines = ["Food chain media (HHRAP 2005); produce mg/kg DW, animal products and fish mg/kg FW",
                 f"{'Chemical':26s} {'AG prod':>9s} {'BG prod':>9s} {'beef':>9s} {'milk':>9s} {'pork':>9s} "
                 f"{'eggs':>9s} {'chicken':>9s} {'fish':>9s}"]
        for c in o.chemicals:
            lines.append(f"{c.name[:26]:26s} {g(c.produce_exposed_aboveground_mg_kg_dw):>9s} "
                         f"{g(c.produce_belowground_mg_kg_dw):>9s} {g(c.beef_mg_kg_fw):>9s} {g(c.milk_mg_kg_fw):>9s} "
                         f"{g(c.pork_mg_kg_fw):>9s} {g(c.eggs_mg_kg_fw):>9s} {g(c.chicken_mg_kg_fw):>9s} "
                         f"{g(c.fish_mg_kg_fw):>9s}")
            lines += [f"  note: {n}" for n in c.notes]
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
