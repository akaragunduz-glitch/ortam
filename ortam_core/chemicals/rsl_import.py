"""Import the EPA Regional Screening Level (RSL) tables into ORTAM data files.

Source files (docs/references/, see README.md there):
  rsl_summary_table_<version>.xlsx        -> toxicity.csv
  rsl_chemical_parameters_<version>.xlsx  -> chemical_parameters.csv
  rsl_subchronic_toxicity_<version>.xlsx  -> subchronic_toxicity.csv

Output folder: data/reference/rsl/<version>/ with a metadata.yaml that records
the source files and their SHA-256 checksums, so that every result can be
traced back to the exact table version.

The import only copies values; it performs no unit conversion. Units are
those of the RSL tables and are carried in the column names.

Fortran analogy: this is a one-off pre-processing program that converts a
vendor input deck into the fixed-format data file the solver reads.

Run:  uv run python -m ortam_core.chemicals.rsl_import --version 2024-11
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from datetime import date
from pathlib import Path

import openpyxl
import yaml

ROOT = Path(__file__).resolve().parents[2]
REFERENCES = ROOT / "docs" / "references"
DATA = ROOT / "data" / "reference" / "rsl"

# Summary table: column index (0-based) -> output column name.
# Header row 3 of worksheet "Summary 1124" (RSL Summary Table, Nov 2024).
SUMMARY_COLUMNS = {
    13: "name",
    14: "cas",
    0: "sfo_per_mg_kg_day",
    1: "sfo_key",
    2: "iur_per_ug_m3",
    3: "iur_key",
    4: "rfd_mg_kg_day",
    5: "rfd_key",
    6: "rfc_mg_m3",
    7: "rfc_key",
    8: "volatile",
    9: "mutagen",
    10: "giabs",
    11: "absd",
    12: "csat_mg_kg",
    15: "sl_resident_soil_mg_kg",
    16: "sl_resident_soil_key",
    17: "sl_industrial_soil_mg_kg",
    18: "sl_industrial_soil_key",
    19: "sl_resident_air_ug_m3",
    20: "sl_resident_air_key",
    21: "sl_industrial_air_ug_m3",
    22: "sl_industrial_air_key",
    23: "sl_tapwater_ug_l",
    24: "sl_tapwater_key",
    25: "mcl_ug_l",
    26: "ssl_risk_based_mg_kg",
    27: "ssl_risk_based_key",
    28: "ssl_mcl_based_mg_kg",
}

# Chemical-specific parameters: header text (normalized) -> output column name.
PARAMETER_COLUMNS = {
    "Analyte": "name",
    "CAS No.": "cas",
    "MW": "mw_g_mol",
    "MW Ref": "mw_ref",
    "H` (unitless)": "henry_dimensionless",
    "HLC (atm-m3/mole)": "henry_atm_m3_mol",
    "H` and HLC Ref": "henry_ref",
    "VP (mmHg)": "vp_mmhg",
    "VP Ref": "vp_ref",
    "MP (C)": "mp_c",
    "MP Ref": "mp_ref",
    "Density (g/cm3)": "density_g_cm3",
    "Density Ref": "density_ref",
    "Dia (cm2/s)": "dia_cm2_s",
    "Diw (cm2/s)": "diw_cm2_s",
    "Dia and Diw Ref": "dia_diw_ref",
    "Kd (L/kg)": "kd_l_kg",
    "Kd Ref": "kd_ref",
    "Koc (L/kg)": "koc_l_kg",
    "Koc Ref": "koc_ref",
    "log Kow (unitless)": "log_kow",
    "log Kow Ref": "log_kow_ref",
    "S (mg/L)": "solubility_mg_l",
    "S Ref": "solubility_ref",
    "B (unitless)": "b_dimensionless",
    "τevent (hr/event)": "tau_event_hr",
    "t* (hr)": "t_star_hr",
    "Kp (cm/hr)": "kp_cm_hr",
    "Kp Ref": "kp_ref",
}


def _norm(text) -> str:
    return " ".join(str(text).replace("\n", " ").split()) if text is not None else ""


def _cell(value) -> str:
    """Excel cell -> CSV text. Empty cells and blank strings become ''."""
    if value is None:
        return ""
    if isinstance(value, float):
        return repr(value)
    return str(value).strip()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(path: Path, first_data_row: int):
    ws = openpyxl.load_workbook(path, read_only=True, data_only=True).active
    return ws.title, list(ws.iter_rows(min_row=first_data_row, values_only=True))


def _write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def import_summary(src: Path, out: Path) -> tuple[str, int]:
    sheet, rows = _rows(src, first_data_row=4)
    cols = list(SUMMARY_COLUMNS.items())
    data = [[_cell(r[i]) for i, _ in cols] for r in rows if r[13]]
    _write_csv(out, [name for _, name in cols], data)
    return sheet, len(data)


def import_by_header(src: Path, out: Path, header_row: int, mapping: dict[str, str]) -> tuple[str, int]:
    """Import a worksheet whose columns are identified by the header row."""
    ws = openpyxl.load_workbook(src, read_only=True, data_only=True).active
    rows = list(ws.iter_rows(values_only=True))
    header = [_norm(h) for h in rows[header_row - 1]]
    index = {}
    for i, h in enumerate(header):
        if not h:
            continue
        index[i] = mapping.get(h, _slug(h))
    data = [[_cell(r[i]) if i < len(r) else "" for i in index] for r in rows[header_row:] if r and r[0]]
    _write_csv(out, list(index.values()), data)
    return ws.title, len(data)


def _slug(header: str) -> str:
    keep = "".join(c.lower() if c.isalnum() else "_" for c in header)
    return "_".join(p for p in keep.split("_") if p)


def import_rsl(version: str, references: Path = REFERENCES, out_dir: Path | None = None) -> Path:
    out_dir = out_dir or DATA / version
    out_dir.mkdir(parents=True, exist_ok=True)
    sources = {
        "toxicity.csv": references / f"rsl_summary_table_{version}.xlsx",
        "chemical_parameters.csv": references / f"rsl_chemical_parameters_{version}.xlsx",
        "subchronic_toxicity.csv": references / f"rsl_subchronic_toxicity_{version}.xlsx",
    }
    meta = {"rsl_version": version, "imported_on": date.today().isoformat(), "files": {}}

    sheet, n = import_summary(sources["toxicity.csv"], out_dir / "toxicity.csv")
    meta["files"]["toxicity.csv"] = {"rows": n, "worksheet": sheet}
    sheet, n = import_by_header(sources["chemical_parameters.csv"], out_dir / "chemical_parameters.csv",
                                header_row=3, mapping=PARAMETER_COLUMNS)
    meta["files"]["chemical_parameters.csv"] = {"rows": n, "worksheet": sheet}
    sheet, n = import_by_header(sources["subchronic_toxicity.csv"], out_dir / "subchronic_toxicity.csv",
                                header_row=_find_header_row(sources["subchronic_toxicity.csv"]), mapping={
                                    "Analyte": "name", "CAS No.": "cas"})
    meta["files"]["subchronic_toxicity.csv"] = {"rows": n, "worksheet": sheet}

    for out_name, src in sources.items():
        meta["files"][out_name].update({"source": src.name, "source_sha256": _sha256(src)})
    (out_dir / "metadata.yaml").write_text(yaml.safe_dump(meta, sort_keys=False), encoding="utf-8")
    return out_dir


def _find_header_row(path: Path) -> int:
    ws = openpyxl.load_workbook(path, read_only=True, data_only=True).active
    for i, r in enumerate(ws.iter_rows(max_row=10, values_only=True), start=1):
        if r and any(_norm(c) == "CAS No." for c in r):
            return i
    raise ValueError(f"No 'CAS No.' header found in {path.name}")


def main() -> None:
    p = argparse.ArgumentParser(description="Import EPA RSL tables into data/reference/rsl/<version>/")
    p.add_argument("--version", default="2024-11", help="RSL table version, e.g. 2024-11")
    a = p.parse_args()
    print(f"Written: {import_rsl(a.version)}")


if __name__ == "__main__":
    main()
