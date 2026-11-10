"""Regenerate data/reference/bioscreen/v1.4/reference_cases.json with EPA/AFCEE's BIOSCREEN 1.4 workbook.

Not part of the test run. Requires LibreOffice and its Python-UNO bridge (run with the system
python3 that can `import uno`, not with the project environment).

BIOSCRN4.xls (in docs/references/bioscreen_1.4.zip) computes everything in worksheet formulas.
The error function is not Excel's: the formulas call the user functions gerf/gerfC of the
workbook's VBA module "Macro Module". LibreOffice imports that module but does not resolve the
functions from worksheet formulas, so this script copies the two functions verbatim from the
workbook's own module into the document's Standard Basic library (with Option VBASupport 1),
and then recalculates. For the distributed data set (Keesler AFB) this reproduces the values
cached by Excel in the .xls file (checked in the tests).

For each case the script loads a fresh copy of the workbook, writes the inputs into worksheet
"Input" (cells listed in INPUT_CELLS; cells not given keep their value or formula), selects
each of the three models in 'Plume Output'!P4 ("t" no degradation, "f" first-order decay,
"i" instantaneous reaction), recalculates and reads the outputs listed in OUTPUTS.

Run:  python3 modules/napl/tests/generate_bioscreen_reference_cases.py
"""
from __future__ import annotations

import hashlib, json, os, re, subprocess, tempfile, time, zipfile
from pathlib import Path

import uno  # noqa: F401  (LibreOffice Python-UNO)
from com.sun.star.beans import PropertyValue

ROOT = Path(__file__).resolve().parents[3]
ZIP = ROOT / "docs/references/bioscreen_1.4.zip"
OUT = ROOT / "data/reference/bioscreen/v1.4/reference_cases.json"

# Input cells of worksheet "Input" (BIOSCREEN 1.4 units: ft, yr, cm/s, mg/L, kg).
INPUT_CELLS = {
    "vs_ft_yr": "E5", "k_cm_s": "E7", "gradient": "E8", "porosity": "E9",
    "alpha_x_ft": "E12", "alpha_y_ft": "E13", "alpha_z_ft": "E14", "plume_length_ft": "E16",
    "retardation": "E19", "rho_kg_l": "E21", "koc_l_kg": "E22", "foc": "E23",
    "lambda_per_yr": "E26", "half_life_yr": "E28",
    "do_mg_l": "E30", "no3_mg_l": "E31", "fe_mg_l": "E32", "so4_mg_l": "E33", "ch4_mg_l": "E34",
    "length_ft": "J5", "width_ft": "J6", "time_yr": "J7", "thickness_ft": "K10",
    "y1_ft": "H13", "c1_mg_l": "I13", "y2_ft": "H14", "c2_mg_l": "I14", "y3_ft": "H15", "c3_mg_l": "I15",
    "mass_kg": "I21",
}
# Model-independent outputs.
OUTPUTS = {
    "Input": ["E5", "E12", "E13", "E19", "E26", "E73", "F52", "F59", "F60", "F61", "F62", "F63",
              "F65", "F66", "F67", "F68", "L45", "L46", "L47", "Q45", "Q46", "Q47", "T45", "T46", "T47"],
    "Centerline Output": ["B35:L35", "B51:L51", "B60:L60", "B83:L83"],
    "Plume Output": ["C57:M57", "B88:B92", "C88:M92", "C103:M107", "C111:M115", "C119:M123",
                     "P88", "P91", "P101", "P104", "P107", "P111", "P116", "P121", "P125"],
}
# Outputs that depend on the model selected in 'Plume Output'!P4.
MODEL_OUTPUTS = ["B5:L9", "B11:L11", "M17", "M19", "M21", "N22", "I25:M25", "M28", "M30", "M31"]

KEESLER_CAL = {"alpha_x_ft": 32.5, "alpha_y_ft": 3.25}
CASES = [
    {"id": "keesler_as_distributed", "note": "Workbook data set as distributed (Keesler AFB, alpha_x from Lp)", "inputs": {}},
    {"id": "keesler_example1", "note": "BIOSCREEN 1.4 Revisions, Example 1 (Fig. 2-4), W = 200 ft",
     "inputs": dict(KEESLER_CAL)},
    {"id": "keesler_example1_w50", "note": "BIOSCREEN 1.4 Revisions, Example 1 (Fig. 5-6), W = 50 ft",
     "inputs": {**KEESLER_CAL, "width_ft": 50}},
    {"id": "hill_example2", "note": "BIOSCREEN 1.3 manual, Example 2 (Fig. 6-8), Hill AFB Site 870; K = 8.1E-03 cm/s as entered in Fig. 6 (Vs 1609.1 ft/yr; the data table gives 8.05E-03)",
     "inputs": {"k_cm_s": 8.1e-3, "gradient": 0.048, "porosity": 0.25, "alpha_x_ft": 28.5, "alpha_y_ft": 2.85,
                "alpha_z_ft": 0, "plume_length_ft": 1450, "rho_kg_l": 1.7, "koc_l_kg": 38, "foc": 8e-4,
                "half_life_yr": 0.1, "do_mg_l": 5.78, "no3_mg_l": 17, "fe_mg_l": 11.3, "so4_mg_l": 100,
                "ch4_mg_l": 0.414, "length_ft": 1450, "width_ft": 320, "time_yr": 5, "thickness_ft": 10,
                "y1_ft": 50, "c1_mg_l": 0.07, "y2_ft": 25, "c2_mg_l": 2.8, "y3_ft": 100, "c3_mg_l": 9,
                "mass_kg": "Infinite"}},
    {"id": "single_zone_steady", "note": "One zone, infinite source, direct Vs/R/lambda, 100 yr",
     "inputs": {"vs_ft_yr": 100, "porosity": 0.3, "plume_length_ft": 500, "alpha_z_ft": 0, "retardation": 1.0,
                "lambda_per_yr": 1.0, "do_mg_l": 4, "no3_mg_l": 2, "fe_mg_l": 5, "so4_mg_l": 20, "ch4_mg_l": 3,
                "length_ft": 1000, "width_ft": 300, "time_yr": 100, "thickness_ft": 10,
                "y1_ft": 0, "c1_mg_l": 0, "y2_ft": 0, "c2_mg_l": 0, "y3_ft": 50, "c3_mg_l": 10, "mass_kg": "Infinite"}},
    {"id": "three_zones_retarded", "note": "Three zones, finite mass, R from Koc*foc", 
     "inputs": {"koc_l_kg": 200, "foc": 0.002, "alpha_x_ft": 20, "alpha_y_ft": 2, "length_ft": 600, "width_ft": 200,
                "time_yr": 20, "mass_kg": 500}},
    {"id": "vertical_dispersion", "note": "Example 1 with alpha_z = 0.5 ft", "inputs": {**KEESLER_CAL, "alpha_z_ft": 0.5}},
    {"id": "no_electron_acceptors", "note": "BC = 0: instantaneous model equals no degradation",
     "inputs": {**KEESLER_CAL, "do_mg_l": 0, "no3_mg_l": 0, "fe_mg_l": 0, "so4_mg_l": 0, "ch4_mg_l": 0,
                "mass_kg": 100}},
    {"id": "short_time_front", "note": "Plume front inside the model area (t = 1 yr)",
     "inputs": {"vs_ft_yr": 200, "alpha_x_ft": 10, "alpha_y_ft": 1, "length_ft": 400, "width_ft": 100,
                "time_yr": 1, "mass_kg": 50}},
    {"id": "two_zones", "note": "Zones 2 and 3 only (zone 1 empty)",
     "inputs": {"y1_ft": 0, "c1_mg_l": 0, "y2_ft": 20, "c2_mg_l": 1, "y3_ft": 30, "c3_mg_l": 5, "mass_kg": 300,
                "time_yr": 10, "alpha_x_ft": 15, "alpha_y_ft": 1.5}},
    {"id": "equal_zone_concentrations", "note": "Co1 = Co2: zero layer increment, BC added twice (workbook logic)",
     "inputs": {"y1_ft": 10, "c1_mg_l": 2, "y2_ft": 15, "c2_mg_l": 2, "y3_ft": 20, "c3_mg_l": 8, "mass_kg": 400}},
    {"id": "large_plume_length", "note": "Dispersivities from Lp = 1000 ft, lambda from half-life, finite mass",
     "inputs": {"plume_length_ft": 1000, "half_life_yr": 1.0, "length_ft": 1000, "width_ft": 300, "time_yr": 15,
                "mass_kg": 1e4}},
    {"id": "decreasing_inner_zone", "note": "Middle zone lower than outer zone (negative increment)",
     "inputs": {"y1_ft": 20, "c1_mg_l": 5, "y2_ft": 20, "c2_mg_l": 3, "y3_ft": 20, "c3_mg_l": 8, "mass_kg": 800,
                **KEESLER_CAL}},
    {"id": "high_retardation", "note": "R > 6 (Koc 1000, foc 0.01)",
     "inputs": {"koc_l_kg": 1000, "foc": 0.01, "time_yr": 50, "length_ft": 200, "width_ft": 100, "mass_kg": 2000,
                **KEESLER_CAL}},
]


def pv(name, value):
    p = PropertyValue(); p.Name = name; p.Value = value
    return p


def cells(sheet, rng):
    r = sheet.getCellRangeByName(rng).getRangeAddress()
    return [[sheet.getCellByPosition(c, row) for c in range(r.StartColumn, r.EndColumn + 1)]
            for row in range(r.StartRow, r.EndRow + 1)]


def value(cell):
    if cell.getError():
        return f"#ERR{cell.getError()}"
    kind = cell.getType().value
    if kind == "FORMULA":
        kind = {1: "VALUE", 2: "TEXT"}.get(cell.FormulaResultType2, "VALUE")
    if kind == "EMPTY":
        return None
    return cell.getString() if kind == "TEXT" else cell.getValue()


def read(sheet, rng):
    v = [[value(c) for c in row] for row in cells(sheet, rng)]
    if len(v) == 1 and len(v[0]) == 1:
        return v[0][0]
    return v[0] if len(v) == 1 else v


def install_erf(doc):
    """Copy gerf and gerfC verbatim from the workbook's VBA module into the Standard library."""
    libs = doc.BasicLibraries
    src = None
    for lname in libs.getElementNames():
        lib = libs.getByName(lname)
        for mname in lib.getElementNames():
            code = lib.getByName(mname)
            if "Function gerf(" in code:
                src = code
    if src is None:
        raise RuntimeError("gerf not found in the workbook's VBA modules")
    funcs = re.findall(r"^Function gerfC?\(.*?^End Function", src, flags=re.S | re.M)
    funcs = ["\n".join(l for l in f.splitlines() if not l.strip().startswith("Attribute")) for f in funcs]
    code = "Option VBASupport 1\n\n" + "\n\n".join(funcs) + "\n"
    if not libs.hasByName("Standard"):
        libs.createLibrary("Standard")
    std = libs.getByName("Standard")
    if std.hasByName("BioscreenErf"):
        std.removeByName("BioscreenErf")
    std.insertByName("BioscreenErf", code)
    return code


def excel_cached_values(xls: Path) -> dict | None:
    """Values cached by Excel in the distributed .xls (last saved state; needs xlrd)."""
    try:
        import xlrd
    except ImportError:
        return None
    b = xlrd.open_workbook(str(xls))
    cl, po = b.sheet_by_name("Centerline Output"), b.sheet_by_name("Plume Output")
    row = lambda sh, r: [sh.cell_value(r - 1, c) for c in range(1, 12)]
    return {"Centerline Output!B51:L51": row(cl, 51), "Centerline Output!B60:L60": row(cl, 60),
            "Centerline Output!B83:L83": row(cl, 83), "Plume Output!P4": po.cell_value(3, 15),
            **{f"Plume Output!{a}": po.cell_value(r - 1, 12) for a, r in
               (("M17", 17), ("M19", 19), ("M21", 21), ("M28", 28), ("M30", 30), ("M31", 31))}}


def main():
    tmp = Path(tempfile.mkdtemp())
    with zipfile.ZipFile(ZIP) as z:
        z.extract("BIOSCRN4.xls", tmp)
    xls = tmp / "BIOSCRN4.xls"
    sha_zip = hashlib.sha256(ZIP.read_bytes()).hexdigest()
    sha_xls = hashlib.sha256(xls.read_bytes()).hexdigest()
    port = 2002
    proc = subprocess.Popen(["soffice", "--headless", "--invisible", "--norestore", "--nologo",
                             f"--accept=socket,host=localhost,port={port};urp;"],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    local = uno.getComponentContext()
    resolver = local.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local)
    for _ in range(90):
        try:
            ctx = resolver.resolve(f"uno:socket,host=localhost,port={port};urp;StarOffice.ComponentContext")
            break
        except Exception:
            time.sleep(1)
    desktop = ctx.ServiceManager.createInstanceWithContext("com.sun.star.frame.Desktop", ctx)
    results, erf_code = [], None
    try:
        for case in CASES:
            doc = desktop.loadComponentFromURL(uno.systemPathToFileUrl(str(xls)), "_blank", 0,
                                               (pv("Hidden", True), pv("MacroExecutionMode", 4)))
            erf_code = install_erf(doc)
            sh = doc.Sheets
            inp = sh.getByName("Input")
            for key, val in case["inputs"].items():
                c = inp.getCellRangeByName(INPUT_CELLS[key])
                if isinstance(val, str):
                    c.setString(val)
                else:
                    c.setValue(float(val))
            doc.calculateAll()
            rec = {"id": case["id"], "note": case["note"], "inputs": case["inputs"],
                   "input_cells": {k: read(inp, a) for k, a in INPUT_CELLS.items()}, "outputs": {}, "models": {}}
            for sname, rngs in OUTPUTS.items():
                s = sh.getByName(sname)
                for rng in rngs:
                    rec["outputs"][f"{sname}!{rng}"] = read(s, rng)
            po = sh.getByName("Plume Output")
            for model in ("t", "f", "i"):
                po.getCellRangeByName("P4").setString(model)
                doc.calculateAll()
                rec["models"][model] = {rng: read(po, rng) for rng in MODEL_OUTPUTS}
            results.append(rec)
            print(case["id"], "centerline no-decay x=L:", rec["outputs"]["Centerline Output!B51:L51"][-1])
            doc.close(True)
    finally:
        try:
            desktop.terminate()
        except Exception:
            pass
        proc.wait(timeout=60)
    cached = excel_cached_values(xls)
    meta = {"source": "docs/references/bioscreen_1.4.zip (BIOSCRN4.xls)", "sha256_zip": sha_zip,
            "sha256_xls": sha_xls, "engine": "LibreOffice (headless, UNO) with gerf/gerfC copied from the "
            "workbook's 'Macro Module'", "erf_basic_module": erf_code, "units": "ft, yr, cm/s, mg/L, kg, L",
            "excel_cached_values": cached,
            "cases": results}
    OUT.write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print("written", OUT)


if __name__ == "__main__":
    main()
