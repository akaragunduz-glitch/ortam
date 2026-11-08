"""risk.foodchain-ingestion: intake, cancer risk and hazard from eating home-grown food (HHRAP 2005).

Source: HHRAP 2005 Appendix C (docs/references/hhrap_2005_appendix_c.pdf), printed page [PDF page].
  Table C-1-2, p. C-6 [7]   I_ag = [(Pd+Pv+Pr)*CR_ag + Pr*CR_pp + Pr_bg*CR_bg] * F_ag      (mg/kg-day)
  Table C-1-3, p. C-10 [11] I_j  = A_j * CR_j * F_j      (beef, milk, pork, poultry, eggs)
  Table C-1-4, p. C-14 [15] I_fish = C_fish * CR_fish * F_fish
  Table C-1-7, p. C-23 [24] Cancer risk = I * ED * EF * CSF / (AT * 365), AT = 70 yr
  Table C-1-8, p. C-26 [27] HQ = I * ED * EF / (RfD * AT * 365), AT = ED
Consumption rates are body-weight normalized (kg/kg-day): produce DW, animal tissue and
fish FW (HHRAP Tables C-1-2 to C-1-4; data/reference/hhrap/2005/defaults.yaml).
Food concentrations come from foodchain.media (connected in workflows/foodchain) or
from measurements. Toxicity values: RSL table (as in all ORTAM risk tools). HHRAP 2005
predates the ADAF guidance; mutagens are flagged, no ADAF is applied here.
Doses are reported as CDI rows (modules/risk/equations/dose.py), one per food.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Literal

import yaml

from pydantic import BaseModel, Field

from modules.risk.common import ChemicalInputBase, label, one_hit, resolve_chemical
from modules.risk.equations.dose import Dose
from ortam_core.model import ModelInput, ModelOutput, Module

DEFAULTS = Path(__file__).resolve().parents[2] / "data" / "reference" / "hhrap" / "2005" / "defaults.yaml"


@lru_cache
def _defaults() -> dict:
    return yaml.safe_load(DEFAULTS.read_text(encoding="utf-8"))


Receptor = Literal["farmer", "farmer_child", "resident", "resident_child", "fisher", "fisher_child"]
FOODS = ("ag", "pp", "bg", "beef", "milk", "pork", "chicken", "eggs", "fish")


class FoodConcentrations(ModelInput):
    produce_exposed_aboveground_mg_kg_dw: float | None = Field(None, ge=0, description="Pd + Pv + Pr")
    produce_protected_aboveground_mg_kg_dw: float | None = Field(None, ge=0, description="Pr")
    produce_belowground_mg_kg_dw: float | None = Field(None, ge=0, description="Pr_bg")
    beef_mg_kg_fw: float | None = Field(None, ge=0)
    milk_mg_kg_fw: float | None = Field(None, ge=0)
    pork_mg_kg_fw: float | None = Field(None, ge=0)
    chicken_mg_kg_fw: float | None = Field(None, ge=0)
    eggs_mg_kg_fw: float | None = Field(None, ge=0)
    fish_mg_kg_fw: float | None = Field(None, ge=0)


FIELD = {"ag": "produce_exposed_aboveground_mg_kg_dw", "pp": "produce_protected_aboveground_mg_kg_dw",
         "bg": "produce_belowground_mg_kg_dw", "beef": "beef_mg_kg_fw", "milk": "milk_mg_kg_fw",
         "pork": "pork_mg_kg_fw", "chicken": "chicken_mg_kg_fw", "eggs": "eggs_mg_kg_fw", "fish": "fish_mg_kg_fw"}


class ChemicalInput(ChemicalInputBase):
    foods: FoodConcentrations


class Exposure(ModelInput):
    consumption_kg_kg_d: dict[str, float] = Field({}, description="Overrides of CR by food: ag, pp, bg, beef, ...")
    fraction_contaminated: float | None = Field(None, gt=0, le=1)
    ed_yr: float | None = Field(None, gt=0)
    ef_d_yr: float | None = Field(None, gt=0, le=365)
    at_cancer_yr: float | None = Field(None, gt=0)


class IngestionInput(ModelInput):
    receptor: Receptor = "farmer"
    toxicity_table: Literal["rsl-2024-11", "none"] = "rsl-2024-11"
    exposure: Exposure = Exposure()
    chemicals: list[ChemicalInput] = Field(min_length=1)


class IngestionResult(BaseModel):
    name: str
    cas: str
    intake_mg_kg_d: dict[str, float] = Field(description="Daily intake I by food (Tables C-1-2 to C-1-4)")
    risk_by_food: dict[str, float] = {}
    hq_by_food: dict[str, float] = {}
    risk: float | None = None
    hq: float | None = None
    doses: list[Dose] = []
    notes: list[str] = []


class IngestionOutput(ModelOutput):
    receptor: str
    parameters: dict
    chemicals: list[IngestionResult]
    total_risk: float | None = None
    hazard_index: float | None = None


class FoodchainIngestionModule(Module):
    name = "risk.foodchain-ingestion"
    title = "Risk: ingestion of home-grown produce, animal products and fish (HHRAP 2005)"
    input_type = IngestionInput
    output_type = IngestionOutput

    def run(self, inp: IngestionInput) -> IngestionOutput:
        d = _defaults()
        rec = inp.receptor
        cr = {**d["consumption"][rec], **inp.exposure.consumption_kg_kg_d}
        F = inp.exposure.fraction_contaminated or d["consumption"]["fraction_contaminated"]["value"]
        ed = inp.exposure.ed_yr or d["exposure"]["ed_yr"][rec]
        ef = inp.exposure.ef_d_yr or d["exposure"]["ef_d_yr"]
        at_ca = inp.exposure.at_cancer_yr or d["exposure"]["at_cancer_yr"]
        warnings = []
        if inp.exposure.model_dump(exclude_defaults=True):
            warnings.append("HHRAP default exposure parameters overridden.")
        results, total, hi = [], 0.0, 0.0
        any_ca = any_nc = False
        for ci in inp.chemicals:
            c, overridden = resolve_chemical(inp.toxicity_table, ci)
            r = IngestionResult(name=c.name, cas=c.cas, intake_mg_kg_d={})
            if overridden:
                warnings.append(f"{label(c)}: table values overridden by the user: {', '.join(overridden)}.")
            if c.mutagen:
                r.notes.append("Mutagenic mode of action: HHRAP 2005 applies no ADAF; risk may be underestimated "
                               "for early-life exposure.")
            for food in FOODS:
                conc = getattr(ci.foods, FIELD[food])
                if conc is None:
                    continue
                if food not in cr:
                    r.notes.append(f"{food}: no consumption rate for receptor {rec}; not included.")
                    continue
                r.intake_mg_kg_d[food] = conc * cr[food] * F          # C-1-2, C-1-3, C-1-4
            for food, I in r.intake_mg_kg_d.items():
                cdi_ca = I * ed * ef / (at_ca * 365.0)                  # C-1-7 without CSF
                cdi_nc = I * ed * ef / (ed * 365.0)                     # C-1-8, AT = ED
                r.doses.append(Dose(route="ingestion", quantity="CDI", basis="cancer", segment=food,
                                    value=cdi_ca, unit="mg/kg-day"))
                r.doses.append(Dose(route="ingestion", quantity="CDI", basis="noncancer", segment=food,
                                    value=cdi_nc, unit="mg/kg-day"))
                if c.sfo_per_mg_kg_day:
                    r.risk_by_food[food] = cdi_ca * c.sfo_per_mg_kg_day
                if c.rfd_mg_kg_day:
                    r.hq_by_food[food] = cdi_nc / c.rfd_mg_kg_day
            if r.risk_by_food:
                r.risk = one_hit(sum(r.risk_by_food.values()))
                total += sum(r.risk_by_food.values()); any_ca = True
            if r.hq_by_food:
                r.hq = sum(r.hq_by_food.values())
                hi += r.hq; any_nc = True
            results.append(r)
        return IngestionOutput(receptor=rec, parameters={"consumption_kg_kg_d": cr, "fraction": F, "ed_yr": ed,
                                                         "ef_d_yr": ef, "at_cancer_yr": at_ca},
                               chemicals=results, total_risk=one_hit(total) if any_ca else None,
                               hazard_index=hi if any_nc else None, warnings=warnings)

    def summary(self, o: IngestionOutput) -> str:
        lines = [f"Food chain ingestion, receptor {o.receptor} (HHRAP 2005 App. C)",
                 f"{'Chemical':30s} {'Intake mg/kg-d':>15s} {'Risk':>9s} {'HQ':>9s}"]
        for c in o.chemicals:
            lines.append(f"{c.name[:30]:30s} {sum(c.intake_mg_kg_d.values()):15.3g} "
                         f"{(f'{c.risk:.1e}' if c.risk is not None else '-'):>9s} "
                         f"{(f'{c.hq:.2g}' if c.hq is not None else '-'):>9s}")
        if o.total_risk is not None:
            lines.append(f"Total cancer risk : {o.total_risk:.1e}")
        if o.hazard_index is not None:
            lines.append(f"Hazard index      : {o.hazard_index:.2g}")
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
