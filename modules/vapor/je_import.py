"""Import the data tables of EPA's J&E spreadsheet (v6.0, September 2017) into ORTAM.

Source: docs/references/je_spreadsheet_v6_2017.xlsm (see README.md there).
Worksheets read:
  CHEM_DATA  -> chemicals.csv  (rows 6-288, 283 chemicals; physical-chemical properties only)
  SOIL_DATA  -> soils.csv      (rows 6-17, the 12 US SCS soil texture classes)
  BLDG_DATA  -> buildings.csv  (rows 10-14, residential default buildings)

Output folder: data/reference/je/v6-2017/ with metadata.yaml (source file and SHA-256).
Values are copied, not converted; units are carried in the column names.
Derived quantities that the spreadsheet computes with cell formulas (capillary
zone porosity and height, building ventilation rate, saturated vapor concentration)
are NOT imported: ORTAM computes them in modules/vapor/je_equations.py, where each
equation cites its source.
Toxicity values of CHEM_DATA are not imported; ORTAM uses the RSL tables for risk.

Fortran analogy: a one-off pre-processing program that converts the vendor's
input deck into fixed-format data files read by the solver.

Run:  uv run python -m modules.vapor.je_import
"""

from __future__ import annotations

import csv
import hashlib
from datetime import date
from pathlib import Path

import openpyxl
import yaml

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "docs" / "references" / "je_spreadsheet_v6_2017.xlsm"
OUT = ROOT / "data" / "reference" / "je" / "v6-2017"

# CHEM_DATA: column letter -> output name (header rows 3-5 of the worksheet)
CHEM_COLUMNS = {
    "A": "name",
    "B": "cas",
    "C": "mw_g_mol",
    "E": "vp_mmhg",
    "I": "solubility_mg_l",
    "K": "henry_atm_m3_mol",
    "M": "henry_dimensionless_25c",
    "N": "dair_cm2_s",
    "P": "dwater_cm2_s",
    "R": "tboil_k",
    "T": "tcrit_k",
    "V": "dhvb_cal_mol",
    "X": "koc_cm3_g",
}
# SOIL_DATA: column letter -> output name (header rows 2-4)
SOIL_COLUMNS = {
    "A": "soil",
    "B": "scs_symbol",
    "C": "n_total_porosity",
    "D": "nw_water_filled_porosity",
    "E": "rho_bulk_density_g_cm3",
    "F": "ks_cm_h",
    "I": "vg_alpha1_per_cm",
    "J": "vg_n",
    "K": "vg_m",
    "L": "theta_r",
    "M": "mean_grain_diameter_cm",
}
# BLDG_DATA residential rows: column letter -> output name (header rows 2-4)
BLDG_COLUMNS = {
    "B": "foundation",
    "C": "lb_m",
    "D": "lf_m",
    "E": "eta",
    "F": "abf_m2",
    "G": "hb_m",
    "H": "ach_per_h",
    "J": "qsoil_qb",
}


def _value(v):
    """Numbers stay numbers; 'Not Available', 'NA' and blanks become empty fields."""
    if isinstance(v, (int, float)):
        return v
    if isinstance(v, str):
        s = v.strip()
        if s.startswith("="):            # e.g. '=2.44' in BLDG_DATA!G14
            try:
                return float(s[1:])
            except ValueError:
                raise ValueError(f"Unexpected formula: {s}")
        if s in ("", "NA", "Not Available"):
            return ""
        try:
            return float(s)
        except ValueError:
            return s
    return ""


def _rows(ws, first: int, last: int, columns: dict[str, str]) -> list[dict]:
    rows = []
    for r in range(first, last + 1):
        if ws[f"A{r}"].value in (None, ""):
            continue
        rows.append({name: _value(ws[f"{col}{r}"].value) for col, name in columns.items()})
    return rows


def _write(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    wb = openpyxl.load_workbook(SOURCE, keep_vba=False)
    chems = _rows(wb["CHEM_DATA"], 6, 288, CHEM_COLUMNS)
    for c in chems:
        c["name"], c["cas"] = str(c["name"]).strip(), str(c["cas"]).strip()   # text columns
        for k, v in c.items():            # 'No Tcrit', 'No DHv,b', ... -> missing value
            if k not in ("name", "cas") and isinstance(v, str):
                c[k] = ""
    soils = _rows(wb["SOIL_DATA"], 6, 17, SOIL_COLUMNS)
    bldg = _rows(wb["BLDG_DATA"], 10, 14, BLDG_COLUMNS)
    OUT.mkdir(parents=True, exist_ok=True)
    _write(OUT / "chemicals.csv", chems)
    _write(OUT / "soils.csv", soils)
    _write(OUT / "buildings.csv", bldg)
    sha = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    meta = {
        "je_version": "6.0 (September 2017)",
        "imported_on": date.today().isoformat(),
        "source": SOURCE.name,
        "source_sha256": sha,
        "files": {
            "chemicals.csv": {"rows": len(chems), "worksheet": "CHEM_DATA", "source_rows": "6-288"},
            "soils.csv": {"rows": len(soils), "worksheet": "SOIL_DATA", "source_rows": "6-17"},
            "buildings.csv": {"rows": len(bldg), "worksheet": "BLDG_DATA", "source_rows": "10-14 (Residential)"},
            "reference_cases.csv": {
                "rows": 16,
                "note": "Spreadsheet results for 16 cases, recalculated with LibreOffice; "
                        "see modules/vapor/tests/generate_je_reference_cases.py",
            },
        },
    }
    (OUT / "metadata.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True), encoding="utf-8")
    print(f"{len(chems)} chemicals, {len(soils)} soils, {len(bldg)} buildings -> {OUT}")


if __name__ == "__main__":
    main()
