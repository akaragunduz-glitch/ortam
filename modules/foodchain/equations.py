"""Food chain concentrations from soil, water and sediment (HHRAP 2005 Appendix B).

Source: HHRAP 2005 Appendix B (docs/references/hhrap_2005_appendix_b.pdf), printed page [PDF page].
ORTAM starts from measured (or modelled) media concentrations; the HHRAP air-deposition
terms Pd (Table B-2-7) and Pv (Table B-2-8) are optional inputs here (from an air model).

  Table B-2-9,  p. B-80 [83]   Pr_ag = Cs * Br_ag                        (aboveground produce, root uptake)
  Table B-2-10, p. B-83 [86]   Pr_bg = Cs * Br_rootveg * VG_rootveg     (belowground produce)
  Table B-3-9,  p. B-134 [137] Pr_i  = Cs * Br_i                        (forage, silage: Br_forage; grain: Br_grain)
  Table B-3-10, p. B-137 [140] A_beef = (sum_i F_i*Qp_i*P_i + Qs*Cs*Bs) * Ba_beef * MF
  Table B-3-11, p. B-144 [147] A_milk = (...) * Ba_milk * MF
  Table B-3-12, p. B-151 [154] A_pork = (...) * Ba_pork * MF
  Table B-3-13, p. B-158 [161] A_egg = (sum_i F_i*Qp_i*P_i + Qs*Cs*Bs) * Ba_egg
  Table B-3-14, p. B-163 [166] A_chicken = (...) * Ba_chicken
  Table B-4-26, p. B-267 [270] C_fish = C_dw * BCF_fish
  Table B-4-27, p. B-270 [273] C_fish = C_dw * BAF_fish
  Table B-4-28, p. B-273 [276] C_fish = C_sb * f_lipid * BSAF / OC_sed
Plant concentrations are dry weight (mg/kg DW); animal tissue and fish fresh weight.
"""

from __future__ import annotations


def produce_root_uptake(cs_mg_kg: float, br: float) -> float:
    """Pr = Cs * Br, mg/kg DW (Tables B-2-9 and B-3-9)."""
    return cs_mg_kg * br


def produce_belowground(cs_mg_kg: float, br_rootveg: float, vg_rootveg: float) -> float:
    """Pr_bg = Cs * Br_rootveg * VG_rootveg, mg/kg DW (Table B-2-10)."""
    return cs_mg_kg * br_rootveg * vg_rootveg


def vg_rootveg(log_kow: float | None, defaults: dict) -> float:
    """0.01 for log Kow > 4, otherwise 1.0 (Table B-2-10, p. B-85 [88])."""
    v = defaults["plants"]["vg_rootveg"]
    return v["value_log_kow_gt_4"] if log_kow is not None and log_kow > 4.0 else v["value_otherwise"]


def animal_intake(feeds: dict[str, float], qp: dict[str, float], f: float, qs: float, cs: float,
                  bs: float) -> float:
    """Daily COPC intake of the animal, mg/day: sum_i(F_i*Qp_i*P_i) + Qs*Cs*Bs (Tables B-3-10 to B-3-14)."""
    return sum(f * q * feeds[k] for k, q in qp.items()) + qs * cs * bs


def animal_tissue(intake_mg_d: float, ba_day_kg: float, mf: float = 1.0) -> float:
    """A = intake * Ba * MF, mg/kg FW. MF applies to beef, milk and pork only (Table B-3-10, p. B-140)."""
    return intake_mg_d * ba_day_kg * mf


def fish_bcf(cdw_mg_l: float, bcf_l_kg: float) -> float:
    """C_fish = C_dw * BCF, mg/kg FW (Table B-4-26)."""
    return cdw_mg_l * bcf_l_kg


def fish_baf(cdw_mg_l: float, baf_l_kg: float) -> float:
    """C_fish = C_dw * BAF, mg/kg FW (Table B-4-27)."""
    return cdw_mg_l * baf_l_kg


def fish_bsaf(csb_mg_kg: float, f_lipid: float, bsaf: float, oc_sed: float) -> float:
    """C_fish = C_sb * f_lipid * BSAF / OC_sed, mg/kg FW (Table B-4-28)."""
    return csb_mg_kg * f_lipid * bsaf / oc_sed
