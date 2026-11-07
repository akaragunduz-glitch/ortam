"""RSL import tests. Reference: the EPA RSL workbooks in docs/references/."""

import csv
import re
from pathlib import Path

import openpyxl
import pytest
import yaml

from ortam_core.chemicals.rsl_import import DATA, REFERENCES, import_rsl

VERSION = "2024-11"
DATA_DIR = DATA / VERSION


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


@pytest.fixture(scope="module")
def reimported(tmp_path_factory) -> Path:
    return import_rsl(VERSION, out_dir=tmp_path_factory.mktemp("rsl"))


@pytest.mark.parametrize("name", ["toxicity.csv", "chemical_parameters.csv", "subchronic_toxicity.csv"])
def test_committed_data_is_reproducible(reimported, name):
    # The CSV files in data/ must be exactly what the importer produces.
    assert (reimported / name).read_text(encoding="utf-8") == (DATA_DIR / name).read_text(encoding="utf-8")


def test_metadata_records_source_checksums():
    meta = yaml.safe_load((DATA_DIR / "metadata.yaml").read_text(encoding="utf-8"))
    assert meta["rsl_version"] == VERSION
    assert meta["files"]["toxicity.csv"]["rows"] == 890
    for f in meta["files"].values():
        assert re.fullmatch(r"[0-9a-f]{64}", f["source_sha256"])


# Ten chemicals checked value by value against the workbook (read independently).
SPOT_CHECK = ["30560-19-1", "7440-36-0", "1912-24-9", "7440-39-3", "7440-41-7",
              "7440-48-4", "7440-50-8", "60-57-1", "7440-02-0", "87-86-5"]


@pytest.mark.parametrize("cas", SPOT_CHECK)
def test_toxicity_values_match_workbook(cas):
    ws = openpyxl.load_workbook(REFERENCES / f"rsl_summary_table_{VERSION}.xlsx",
                                read_only=True, data_only=True).active
    xl = next(r for r in ws.iter_rows(min_row=4, values_only=True) if r[14] == cas)
    row = next(r for r in read_csv(DATA_DIR / "toxicity.csv") if r["cas"] == cas)

    def same(csv_value: str, xl_value):
        if isinstance(xl_value, (int, float)):
            return float(csv_value) == pytest.approx(xl_value, rel=1e-12)
        return csv_value == ("" if xl_value is None else str(xl_value).strip())

    for col, idx in [("sfo_per_mg_kg_day", 0), ("iur_per_ug_m3", 2), ("rfd_mg_kg_day", 4),
                     ("rfc_mg_m3", 6), ("giabs", 10), ("absd", 11), ("sl_resident_soil_mg_kg", 15)]:
        assert same(row[col], xl[idx]), (cas, col)


def _workbook():
    return openpyxl.load_workbook(REFERENCES / f"rsl_default_exposure_parameters_{VERSION}.xlsx", data_only=True)


def _cell_value(wb, source: str) -> float | str:
    """Read the 'Default' column (C) of '<worksheet>, row <n>' and return its number."""
    sheet, row = re.fullmatch(r"(.+), row (\d+)", source).groups()
    raw = wb[sheet].cell(row=int(row), column=3).value
    if isinstance(raw, (int, float)):
        return float(raw)
    m = re.search(r"(?:Resident\s*-?\s*)?([0-9.]+)", str(raw))  # e.g. "Resident - 350", "Resident 2373"
    return float(m.group(1)) if m else str(raw)


def test_exposure_defaults_match_workbook():
    params = yaml.safe_load((DATA_DIR / "exposure_resident.yaml").read_text(encoding="utf-8"))
    wb = _workbook()
    for key, p in params["parameters"].items():
        src = p["source"]
        if p["symbol"] == "THQ":
            continue  # summary-table basis, documented in the YAML
        sheet, row = re.fullmatch(r"(.+), row (\d+)", src).groups()
        raw = wb[sheet].cell(row=int(row), column=3).value
        if p["symbol"] == "TR":
            assert raw == "1 x 10-6"
        elif p["symbol"] == "PEF":
            assert str(raw).startswith("1.36 x 109")
        else:
            assert p["value"] == pytest.approx(_cell_value(wb, src)), key
        assert str(wb[sheet].cell(row=int(row), column=1).value).strip() == p["symbol"], key
    assert "Arsenic = 0.6" in str(wb["Miscellaneous Variables"].cell(row=5, column=3).value)
    assert params["rba_by_cas"]["7440-38-2"] == 0.6


def test_age_segments_match_workbook():
    params = yaml.safe_load((DATA_DIR / "exposure_resident.yaml").read_text(encoding="utf-8"))
    wb = _workbook()
    for field, sources in params["age_segment_sources"].items():
        if isinstance(sources, str):
            continue  # taken from the RSL User's Guide equation images (checked via factors)
        for seg, src in zip(params["age_segments"], sources):
            assert seg[field] == pytest.approx(_cell_value(wb, src)), (field, seg["label"])


def test_worker_defaults_match_workbook():
    params = yaml.safe_load((DATA_DIR / "exposure_worker.yaml").read_text(encoding="utf-8"))
    wb = _workbook()
    for receptor, fields in params["receptors"].items():
        for key, p in fields.items():
            if p["source"].startswith("RSL UG"):
                continue  # EWcon, DWcon: equation images; EF = EW*DW checked below
            assert p["value"] == pytest.approx(_cell_value(wb, p["source"])), (receptor, key)
    con = params["receptors"]["construction"]
    assert con["ew_wk_yr"]["value"] * con["dw_d_wk"]["value"] == con["ef_d_yr"]["value"]


def test_recreator_defaults_match_workbook():
    params = yaml.safe_load((DATA_DIR / "exposure_recreator.yaml").read_text(encoding="utf-8"))
    wb = _workbook()
    for key, p in params["defaults"].items():
        if isinstance(p["sources"], str):
            continue
        for v, src in zip(p["values"], p["sources"]):
            assert v == pytest.approx(_cell_value(wb, src)), (key, src)


def test_air_emission_defaults_match_workbook():
    params = yaml.safe_load((DATA_DIR / "air_emission_factors.yaml").read_text(encoding="utf-8"))
    wb = _workbook()
    for group, fields in params.items():
        if not isinstance(fields, dict):
            continue
        for key, p in fields.items():
            if p["source"].startswith("RSL UG"):
                continue
            raw = _cell_value(wb, p["source"])
            if key == "t_s":
                assert p["value"] == pytest.approx(8.2e8)   # cell text "8.2 x 10^8 (used for unlimited ...)"
                continue
            assert p["value"] == pytest.approx(raw), (group, key)


def test_soil_to_groundwater_defaults_match_workbook():
    params = yaml.safe_load((DATA_DIR / "soil_to_groundwater.yaml").read_text(encoding="utf-8"))
    wb = _workbook()
    for key, p in params["defaults"].items():
        assert p["value"] == pytest.approx(_cell_value(wb, p["source"])), key
