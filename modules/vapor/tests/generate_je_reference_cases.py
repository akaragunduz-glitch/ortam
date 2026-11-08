"""Regenerate data/reference/je/v6-2017/reference_cases.csv with EPA's J&E spreadsheet.

Not part of the test run. Requires LibreOffice (soffice) and a recalculation
script that opens a workbook headless, recalculates all formulas and saves it
(e.g. the recalc.py of the xlsx tools; path given with --recalc).
For each case the script writes the inputs into worksheet MODEL of
docs/references/je_spreadsheet_v6_2017.xlsm (F10 source, F11 Cmedium, F12 Ls,
F13 Ts, F17 chemical, F37 foundation, F55-F68 strata; the thickness of a
"Not Present" stratum is left blank, as the spreadsheet requires), recalculates,
and reads the results (F14 Cs ... G100 Css; K93/L93 error flags).

Run:  python modules/vapor/tests/generate_je_reference_cases.py --recalc /path/to/recalc.py
"""
import argparse, csv, subprocess, sys, tempfile
from pathlib import Path
import openpyxl

ROOT = Path(__file__).resolve().parents[3]
XLSM = ROOT / "docs/references/je_spreadsheet_v6_2017.xlsm"
OUTCSV = ROOT / "data/reference/je/v6-2017/reference_cases.csv"
NP = "Not Present"
CASES = [
 # id, source, chem, Cmed, Ls, Ts, foundation, (SCS_A,hSA), (SCS_B,hSB), (SCS_C,hSC)
 (1, "Groundwater", "Benzene", 100, 3, 25, "Slab-on-grade", ("Sand", 3), (NP, 0), (NP, 0)),
 (2, "Groundwater", "Trichloroethylene", 50, 5, 15, "Basement w/ slab", ("Loam", 5), (NP, 0), (NP, 0)),
 (3, "Groundwater", "Tetrachloroethylene", 20, 6, 12, "Slab-on-grade", ("Clay", 6), (NP, 0), (NP, 0)),
 (4, "Groundwater", "Vinyl Chloride", 10, 8, 20, "Closed crawl space w/ slab", ("Silty Clay", 8), (NP, 0), (NP, 0)),
 (5, "Groundwater", "Benzene", 100, 4, 25, "Basement w/ dirt floor", ("Sand", 4), (NP, 0), (NP, 0)),
 (6, "Groundwater", "Trichloroethylene", 50, 3, 25, "Closed crawl space w/ dirt floor", ("Loamy Sand", 3), (NP, 0), (NP, 0)),
 (7, "Groundwater", "Trichloroethylene", 50, 5, 18, "Slab-on-grade", ("Sandy Loam", 2), ("Clay", 3), (NP, 0)),
 (8, "Groundwater", "Tetrachloroethylene", 20, 6, 16, "Slab-on-grade", ("Loam", 1), ("Sand", 2), ("Silt", 3)),
 (9, "Exterior Soil Gas", "Trichloroethylene", 1000, 3, 25, "Slab-on-grade", ("Sand", 3), (NP, 0), (NP, 0)),
 (10, "Exterior Soil Gas", "Benzene", 500, 4, 25, "Basement w/ slab", ("Loam", 4), (NP, 0), (NP, 0)),
 (11, "Exterior Soil Gas", "Tetrachloroethylene", 2000, 4, 10, "Slab-on-grade", ("Silt Loam", 1.5), ("Sand", 2.5), (NP, 0)),
 (12, "Exterior Soil Gas", "Vinyl Chloride", 300, 5, 25, "Basement w/ dirt floor", ("Clay", 5), (NP, 0), (NP, 0)),
 (13, "Sub-slab Soil Gas", "Trichloroethylene", 1000, 0.5, 25, "Slab-on-grade", ("Sand", 0.5), (NP, 0), (NP, 0)),
 (14, "Sub-slab Soil Gas", "Benzene", 500, 2.5, 25, "Basement w/ slab", ("Loam", 2.5), (NP, 0), (NP, 0)),
 (15, "Groundwater", "Naphthalene", 30, 4, 14, "Slab-on-grade", ("Silt Loam", 4), (NP, 0), (NP, 0)),
 (16, "Groundwater", "Acetaldehyde", 500, 2, 25, "Slab-on-grade", ("Sand", 2), (NP, 0), (NP, 0)),
]
OUT = ["F14", "F26", "F27", "F38", "F39", "F40", "F41", "F42", "F43", "F44", "F45", "F46", "F74", "F75", "F76",
       "G104", "G105", "G106", "G107", "G108", "G111", "G112", "G113", "G93", "G96", "G100", "K93", "L93"]
NAMES = ["Cs_ug_m3", "Hr", "Hs", "Lb_m", "Lf_m", "eta", "Abf_m2", "Hb_m", "ach_per_h", "Qsoil_Qb", "Qb_m3_h",
         "Qsoil_m3_h", "hcz_m", "ncz", "nwcz", "DeffA_cm2_s", "DeffB_cm2_s", "DeffC_cm2_s", "DeffCZ_cm2_s",
         "DeffT_cm2_s", "A_param", "B_param", "C_param", "alpha", "Cia_ug_m3", "Css_ug_m3", "flag", "comment"]
ap = argparse.ArgumentParser()
ap.add_argument("--recalc", required=True)
RECALC = ap.parse_args().recalc
tmp = Path(tempfile.mkdtemp())
rows = []
for case in CASES:
    cid, src, chem, cmed, ls, ts, found, a, b, c = case
    wb = openpyxl.load_workbook(XLSM)
    m = wb["MODEL"]
    m["F10"], m["F11"], m["F12"], m["F13"], m["F17"], m["F37"] = src, cmed, ls, ts, chem, found
    m["F55"], m["F56"], m["F61"], m["F67"] = a[0], a[1], b[0], c[0]
    m["F62"] = b[1] if b[0] != NP else None
    m["F68"] = c[1] if c[0] != NP else None
    fn = str(tmp / f"case{cid}.xlsx")
    wb.save(fn)
    r = subprocess.run([sys.executable, RECALC, fn, "180"], capture_output=True, text=True)
    v = openpyxl.load_workbook(fn, data_only=True)["MODEL"]
    rec = {"case": cid, "source": src, "chemical": chem, "c_medium": cmed, "ls_m": ls, "ts_c": ts, "foundation": found,
           "soil_a": a[0], "h_a_m": a[1], "soil_b": b[0], "h_b_m": b[1], "soil_c": c[0], "h_c_m": c[1]}
    for cell, name in zip(OUT, NAMES):
        rec[name] = v[cell].value
    rows.append(rec)
    print(cid, chem, src, found, "alpha=", rec["alpha"], "Cia=", rec["Cia_ug_m3"], "flag=", rec["flag"], rec["comment"])
with open(OUTCSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
