"""Resident soil: risk-based screening levels, cancer risk and hazard quotients.

Exposure routes: incidental soil ingestion, dermal contact with soil, and
inhalation of particulates (PEF) and, if a volatilization factor is supplied,
volatiles (VF).

Equations (no equation is written from memory; every one cites its source):
  - Screening levels: EPA RSL User's Guide (web page of 2026-03-13, saved as
    docs/references/rsl_users_guide_2026-03.pdf), Section 4.1.1.1-4.1.1.3,
    PDF pp. 36-41. The RSL guide does not number its equations; comments cite
    section and page.
  - Risk and HQ from a concentration: RSL User's Guide Section 2.6.1 (PDF p. 28);
    one-hit rule Section 2.6.2 (PDF p. 28) and RAGS Part A (1989) p. 8-11.
  - Total risk and hazard index: RAGS Part A (1989) p. 8-12 and p. 8-13.
Default parameters: data/reference/rsl/2024-11/exposure_resident.yaml.
Verification: modules/risk/tests/test_resident_soil.py (RSL Summary Table,
November 2024; all non-volatile, non-mutagenic chemicals).

Not covered (raise an error rather than return a wrong number):
  - mutagens (age-dependent adjustment factors): planned task D1.1b;
  - volatiles without a user-supplied VF: VF calculation is task D1.2.

Fortran analogy: every SL function below is a pure FUNCTION of the
parameters; `run` is the main program that loops over chemicals.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Literal

import yaml
from pydantic import BaseModel, Field

from modules.risk.rsl_import import DATA
from modules.risk.toxicity import Chemical, lookup
from ortam_core.model import ModelInput, ModelOutput, Module

CEILING_MG_KG = 1.0e5        # RSL UG 5.13: theoretical ceiling limit (flagged, not applied)
ONE_HIT_THRESHOLD = 0.01     # RSL UG 2.6.2: linear risk equation valid below 0.01


# ---------------------------------------------------------------- parameters
@lru_cache(maxsize=None)
def _defaults_file(version: str = "2024-11") -> dict:
    return yaml.safe_load((DATA / version / "exposure_resident.yaml").read_text(encoding="utf-8"))


def load_exposure_defaults(version: str = "2024-11") -> dict[str, float]:
    return {k: float(v["value"]) for k, v in _defaults_file(version)["parameters"].items()}


def age_adjusted_factors(p: dict[str, float]) -> tuple[float, float]:
    """Age-adjusted soil ingestion and dermal factors, mg/kg.

    RSL UG 4.1.1.3, PDF p. 40 (IFSres-adj) and p. 41 (DFSres-adj):
      IFS = EFc*EDc*IRSc/BWc + EFa*EDa*IRSa/BWa
      DFS = EFc*EDc*SAc*AFc/BWc + EFa*EDa*SAa*AFa/BWa
    """
    ifs = (p["ef_child_d_yr"] * p["ed_child_yr"] * p["irs_child_mg_d"] / p["bw_child_kg"]
           + p["ef_adult_d_yr"] * p["ed_adult_yr"] * p["irs_adult_mg_d"] / p["bw_adult_kg"])
    dfs = (p["ef_child_d_yr"] * p["ed_child_yr"] * p["sa_child_cm2_d"] * p["af_child_mg_cm2"] / p["bw_child_kg"]
           + p["ef_adult_d_yr"] * p["ed_adult_yr"] * p["sa_adult_cm2_d"] * p["af_adult_mg_cm2"] / p["bw_adult_kg"])
    return ifs, dfs


def _inv_vf_pef(vf_m3_kg: float | None, pef_m3_kg: float) -> float:
    """(1/VF + 1/PEF), kg/m3; the VF term is omitted when no VF applies."""
    return (1.0 / vf_m3_kg if vf_m3_kg else 0.0) + 1.0 / pef_m3_kg


def _total(*sls: float | None) -> float | None:
    """Total SL = 1 / sum(1/SL_i) over applicable routes (RSL UG 4.1.1.1, p. 37, 'Total')."""
    inv = sum(1.0 / s for s in sls if s)
    return 1.0 / inv if inv else None


# ----------------------------------------------- screening levels (mg/kg)
def sl_noncancer(c: Chemical, rba: float, vf: float | None, p: dict[str, float], age: str,
                 thq: float) -> tuple[float | None, float | None, float | None]:
    """Noncancer SLs by route for age = 'child' or 'adult', mg/kg.

    Child: RSL UG 4.1.1.1, PDF pp. 36-37 (ED = EDres-c, AT = 365*EDres-c).
    Adult: RSL UG 4.1.1.2, PDF pp. 38-39 (ED = EDres, AT = 365*EDres).
    """
    if age == "child":
        ef, ed, et = p["ef_child_d_yr"], p["ed_child_yr"], p["et_child_h_d"]
        bw, irs, sa, af = p["bw_child_kg"], p["irs_child_mg_d"], p["sa_child_cm2_d"], p["af_child_mg_cm2"]
    else:
        ef, ed, et = p["ef_adult_d_yr"], p["ed_res_yr"], p["et_adult_h_d"]
        bw, irs, sa, af = p["bw_adult_kg"], p["irs_adult_mg_d"], p["sa_adult_cm2_d"], p["af_adult_mg_cm2"]
    at = 365.0 * ed  # days

    ing = derm = inh = None
    if c.rfd_mg_kg_day:
        # SL_ing = THQ*AT*BW / ((RBA/RfD) * 1e-6 kg/mg * EF * ED * IRS)
        ing = thq * at * bw / ((rba / c.rfd_mg_kg_day) * 1e-6 * ef * ed * irs)
        if c.absd:
            # SL_derm = THQ*AT*BW / ((1/(RfD*GIABS)) * 1e-6 kg/mg * EF * ED * SA * AF * ABSd)
            derm = thq * at * bw / ((1.0 / (c.rfd_mg_kg_day * c.giabs)) * 1e-6 * ef * ed * sa * af * c.absd)
    if c.rfc_mg_m3:
        # SL_inh = THQ*AT / ((1/RfC) * EF * ED * ET * (1 day/24 h) * (1/VF + 1/PEF))
        inh = thq * at / ((1.0 / c.rfc_mg_m3) * ef * ed * et / 24.0 * _inv_vf_pef(vf, p["pef_m3_kg"]))
    return ing, derm, inh


def sl_cancer(c: Chemical, rba: float, vf: float | None, p: dict[str, float],
              tr: float) -> tuple[float | None, float | None, float | None]:
    """Cancer SLs by route (age-adjusted), mg/kg. RSL UG 4.1.1.3, PDF pp. 40-41."""
    ifs, dfs = age_adjusted_factors(p)
    at = 365.0 * p["lt_yr"]  # days
    ing = derm = inh = None
    if c.sfo_per_mg_kg_day:
        # SL_ing = TR*AT / (SFO * 1e-6 kg/mg * RBA * IFSres-adj)
        ing = tr * at / (c.sfo_per_mg_kg_day * 1e-6 * rba * ifs)
        if c.absd:
            # SL_derm = TR*AT / ((SFO/GIABS) * 1e-6 kg/mg * DFSres-adj * ABSd)
            derm = tr * at / ((c.sfo_per_mg_kg_day / c.giabs) * 1e-6 * dfs * c.absd)
    if c.iur_per_ug_m3:
        # SL_inh = TR*AT / (IUR * 1000 ug/mg * EFres * EDres * ETres * (1 day/24 h) * (1/VF + 1/PEF))
        inh = tr * at / (c.iur_per_ug_m3 * 1000.0 * p["ef_res_d_yr"] * p["ed_res_yr"]
                         * p["et_res_h_d"] / 24.0 * _inv_vf_pef(vf, p["pef_m3_kg"]))
    return ing, derm, inh


def one_hit(linear_risk: float) -> float:
    """RSL UG 2.6.2 and RAGS A (1989) p. 8-11: Risk = 1 - exp(-CDI x SF) above 0.01."""
    return 1.0 - math.exp(-linear_risk) if linear_risk > ONE_HIT_THRESHOLD else linear_risk


# --------------------------------------------------------------- input/output
class ChemicalInput(ModelInput):
    cas: str | None = Field(None, description="CAS registry number")
    name: str | None = Field(None, description="Chemical name; required if the CAS number is ambiguous")
    concentration_mg_kg: float | None = Field(
        None, ge=0, description="Exposure point concentration in soil, mg/kg (e.g. 95% UCL); "
                                "if given, risk and HQ are calculated")
    # Optional overrides of table values (every override is flagged in the output)
    sfo_per_mg_kg_day: float | None = Field(None, gt=0, description="Oral slope factor, (mg/kg-day)^-1")
    iur_per_ug_m3: float | None = Field(None, gt=0, description="Inhalation unit risk, (ug/m3)^-1")
    rfd_mg_kg_day: float | None = Field(None, gt=0, description="Chronic oral reference dose, mg/kg-day")
    rfc_mg_m3: float | None = Field(None, gt=0, description="Chronic reference concentration, mg/m3")
    giabs: float | None = Field(None, gt=0, le=1, description="Gastrointestinal absorption fraction, -")
    absd: float | None = Field(None, ge=0, le=1, description="Dermal absorption fraction from soil, -")
    rba: float | None = Field(None, gt=0, le=1, description="Relative bioavailability, -")
    vf_m3_kg: float | None = Field(None, gt=0, description="Volatilization factor, m3/kg (required for volatiles)")


class ExposureOverrides(ModelInput):
    """Any default in exposure_resident.yaml may be overridden (flagged in the output)."""
    lt_yr: float | None = Field(None, gt=0, description="Lifetime, yr")
    ef_child_d_yr: float | None = Field(None, gt=0, le=365, description="Exposure frequency, child, day/yr")
    ef_adult_d_yr: float | None = Field(None, gt=0, le=365, description="Exposure frequency, adult, day/yr")
    ef_res_d_yr: float | None = Field(None, gt=0, le=365, description="Exposure frequency, resident, day/yr")
    ed_child_yr: float | None = Field(None, gt=0, description="Exposure duration, child, yr")
    ed_adult_yr: float | None = Field(None, gt=0, description="Exposure duration, adult, yr")
    ed_res_yr: float | None = Field(None, gt=0, description="Exposure duration, resident, yr")
    et_child_h_d: float | None = Field(None, gt=0, le=24, description="Exposure time, child, h/day")
    et_adult_h_d: float | None = Field(None, gt=0, le=24, description="Exposure time, adult, h/day")
    et_res_h_d: float | None = Field(None, gt=0, le=24, description="Exposure time, resident, h/day")
    bw_child_kg: float | None = Field(None, gt=0, description="Body weight, child, kg")
    bw_adult_kg: float | None = Field(None, gt=0, description="Body weight, adult, kg")
    irs_child_mg_d: float | None = Field(None, gt=0, description="Soil ingestion rate, child, mg/day")
    irs_adult_mg_d: float | None = Field(None, gt=0, description="Soil ingestion rate, adult, mg/day")
    sa_child_cm2_d: float | None = Field(None, gt=0, description="Skin surface area, child, cm2/day")
    sa_adult_cm2_d: float | None = Field(None, gt=0, description="Skin surface area, adult, cm2/day")
    af_child_mg_cm2: float | None = Field(None, gt=0, description="Soil adherence factor, child, mg/cm2")
    af_adult_mg_cm2: float | None = Field(None, gt=0, description="Soil adherence factor, adult, mg/cm2")
    pef_m3_kg: float | None = Field(None, gt=0, description="Particulate emission factor, m3/kg")


class ResidentSoilInput(ModelInput):
    toxicity_table: Literal["rsl-2024-11", "none"] = Field(
        "rsl-2024-11", description="Toxicity source; 'none' requires all values in the input")
    target_risk: float = Field(1e-6, gt=0, lt=1, description="Target cancer risk TR, -")
    target_hq: float = Field(1.0, gt=0, description="Target hazard quotient THQ, -")
    exposure: ExposureOverrides = ExposureOverrides()
    chemicals: list[ChemicalInput] = Field(min_length=1)


class Routes(BaseModel):
    ingestion: float | None = None
    dermal: float | None = None
    inhalation: float | None = None
    total: float | None = None


class ChemicalResult(BaseModel):
    name: str
    cas: str
    rba: float
    sl_nc_child_mg_kg: Routes
    sl_nc_adult_mg_kg: Routes
    sl_ca_mg_kg: Routes
    sl_mg_kg: float | None = Field(description="Screening level: lower of cancer and child noncancer SL, mg/kg")
    sl_basis: Literal["c", "n"] | None = None
    concentration_mg_kg: float | None = None
    risk: Routes | None = Field(None, description="Cancer risk, - (one-hit rule above 0.01)")
    hq_child: Routes | None = None
    hq_adult: Routes | None = None
    overridden: list[str] = []
    toxicity_keys: dict[str, str] = {}


class ResidentSoilOutput(ModelOutput):
    toxicity_table: str
    target_risk: float
    target_hq: float
    exposure_used: dict[str, float]
    exposure_overridden: list[str] = []
    chemicals: list[ChemicalResult]
    total_risk: float | None = Field(None, description="Total cancer risk, - (RAGS A p. 8-12; one-hit rule)")
    hazard_index_child: float | None = Field(None, description="Hazard index, child, - (RAGS A p. 8-13)")
    hazard_index_adult: float | None = None


# --------------------------------------------------------------------- module
_OVERRIDABLE = ("sfo_per_mg_kg_day", "iur_per_ug_m3", "rfd_mg_kg_day", "rfc_mg_m3", "giabs", "absd")


class ResidentSoilModule(Module):
    name = "risk.resident-soil"
    title = "Risk: resident soil screening levels, cancer risk and hazard"
    input_type = ResidentSoilInput
    output_type = ResidentSoilOutput

    def _chemical(self, inp: ResidentSoilInput, ci: ChemicalInput) -> tuple[Chemical, list[str]]:
        given = {k: getattr(ci, k) for k in _OVERRIDABLE if getattr(ci, k) is not None}
        if inp.toxicity_table == "none":
            if not ci.name:
                raise ValueError("With toxicity_table 'none', each chemical needs a name.")
            base = Chemical(name=ci.name, cas=ci.cas or "")
            return Chemical(**{**base.__dict__, **given}), []
        base = lookup(inp.toxicity_table, ci.cas, ci.name)
        return Chemical(**{**base.__dict__, **given}), sorted(given)

    def run(self, inputs: ResidentSoilInput) -> ResidentSoilOutput:
        defaults_file = _defaults_file()
        p = load_exposure_defaults()
        ex = inputs.exposure.model_dump(exclude_none=True)
        p.update(ex)
        tr, thq = inputs.target_risk, inputs.target_hq
        warnings: list[str] = []
        if ex:
            warnings.append(f"Default exposure parameters overridden: {', '.join(sorted(ex))}.")

        results = []
        for ci in inputs.chemicals:
            c, overridden = self._chemical(inputs, ci)
            label = f"{c.name} ({c.cas})" if c.cas else c.name
            if c.mutagen:
                raise ValueError(f"{label} is a mutagen; mutagenic equations (ADAF) are not yet "
                                 "implemented (task D1.1b).")
            if c.volatile and ci.vf_m3_kg is None:
                raise ValueError(f"{label} is volatile; supply vf_m3_kg (VF calculation is task D1.2).")
            if ci.vf_m3_kg is not None:
                overridden.append("vf_m3_kg")
            rba = ci.rba if ci.rba is not None else float(defaults_file["rba_by_cas"].get(c.cas, 1.0))
            if ci.rba is not None:
                overridden.append("rba")
            if overridden:
                warnings.append(f"{label}: table values overridden by the user: {', '.join(overridden)}.")

            nc_c = sl_noncancer(c, rba, ci.vf_m3_kg, p, "child", thq)
            nc_a = sl_noncancer(c, rba, ci.vf_m3_kg, p, "adult", thq)
            ca = sl_cancer(c, rba, ci.vf_m3_kg, p, tr)
            r_nc_c, r_nc_a, r_ca = (Routes(ingestion=s[0], dermal=s[1], inhalation=s[2], total=_total(*s))
                                    for s in (nc_c, nc_a, ca))

            candidates = [(v, b) for v, b in ((r_ca.total, "c"), (r_nc_c.total, "n")) if v]
            sl, basis = min(candidates) if candidates else (None, None)
            if sl is None:
                warnings.append(f"{label}: no toxicity values; no screening level calculated.")
            elif sl > CEILING_MG_KG:
                warnings.append(f"{label}: SL {sl:.3g} mg/kg exceeds the 1E+05 mg/kg ceiling limit "
                                "(RSL UG 5.13); flagged, not replaced.")
            if sl and c.csat_mg_kg and sl > c.csat_mg_kg:
                warnings.append(f"{label}: SL exceeds Csat ({c.csat_mg_kg:.3g} mg/kg).")

            res = ChemicalResult(name=c.name, cas=c.cas, rba=rba, sl_nc_child_mg_kg=r_nc_c,
                                 sl_nc_adult_mg_kg=r_nc_a, sl_ca_mg_kg=r_ca, sl_mg_kg=sl, sl_basis=basis,
                                 overridden=overridden, toxicity_keys=dict(c.keys))
            if ci.concentration_mg_kg is not None:
                C = ci.concentration_mg_kg
                res.concentration_mg_kg = C
                # RSL UG 2.6.1 (PDF p. 28): Risk = C*TR/SL, HQ = C*THQ/SL, per route and in total
                res.risk = _ratio(C, tr, r_ca, risk=True)
                res.hq_child = _ratio(C, thq, r_nc_c)
                res.hq_adult = _ratio(C, thq, r_nc_a)
            results.append(res)

        with_c = [r for r in results if r.concentration_mg_kg is not None]
        total_risk = hi_c = hi_a = None
        if with_c:
            # RAGS A (1989) p. 8-12: Risk_T = sum(Risk_i), from linear values, then one-hit (RSL UG 2.6.2)
            total_risk = one_hit(sum(_linear(r.concentration_mg_kg, tr, r.sl_ca_mg_kg.total) for r in with_c))
            # RAGS A (1989) p. 8-13: HI = sum(HQ_i)
            hi_c = sum(r.hq_child.total or 0.0 for r in with_c)
            hi_a = sum(r.hq_adult.total or 0.0 for r in with_c)

        return ResidentSoilOutput(toxicity_table=inputs.toxicity_table, target_risk=tr, target_hq=thq,
                                  exposure_used=p, exposure_overridden=sorted(ex), chemicals=results,
                                  total_risk=total_risk, hazard_index_child=hi_c, hazard_index_adult=hi_a,
                                  warnings=warnings)

    def summary(self, o: ResidentSoilOutput) -> str:
        lines = [f"Resident soil  (toxicity: {o.toxicity_table}; TR = {o.target_risk:.0e}, THQ = {o.target_hq:g})",
                 f"{'Chemical':32s} {'CAS':12s} {'SL mg/kg':>10s} {'basis':>5s} {'C mg/kg':>10s} "
                 f"{'Risk':>9s} {'HQ child':>9s}"]
        for c in o.chemicals:
            sl = f"{c.sl_mg_kg:.2g}" if c.sl_mg_kg else "-"
            conc = f"{c.concentration_mg_kg:.3g}" if c.concentration_mg_kg is not None else "-"
            risk = f"{c.risk.total:.1e}" if c.risk and c.risk.total is not None else "-"
            hq = f"{c.hq_child.total:.2g}" if c.hq_child and c.hq_child.total is not None else "-"
            flag = " *" if c.overridden else ""
            lines.append(f"{c.name[:32]:32s} {c.cas:12s} {sl:>10s} {c.sl_basis or '-':>5s} {conc:>10s} "
                         f"{risk:>9s} {hq:>9s}{flag}")
        if o.total_risk is not None:
            lines.append(f"Total cancer risk    : {o.total_risk:.1e}")
            lines.append(f"Hazard index (child) : {o.hazard_index_child:.2g}")
            lines.append(f"Hazard index (adult) : {o.hazard_index_adult:.2g}")
        if any(c.overridden for c in o.chemicals) or o.exposure_overridden:
            lines.append("* user-supplied values replace table defaults (see warnings)")
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)


def _linear(C: float, target: float, sl: float | None) -> float:
    return C * target / sl if sl else 0.0


def _ratio(C: float, target: float, sl: Routes, risk: bool = False) -> Routes:
    f = (lambda x: one_hit(x)) if risk else (lambda x: x)
    vals = {k: (f(_linear(C, target, getattr(sl, k))) if getattr(sl, k) else None)
            for k in ("ingestion", "dermal", "inhalation", "total")}
    return Routes(**vals)
