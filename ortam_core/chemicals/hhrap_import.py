"""Import the HHRAP 2005 Companion Database (MS Access) into ORTAM data files.

Source: docs/references/hhrap_2005_chemical_database.mdb (HHRAP Companion Database,
U.S. EPA 2005, EPA530-R-05-006; see docs/references/README.md).
Tables read: ParamVals (values), Definitions (parameter names and units),
ChemNames (names), ParamExplanations (how a value was derived), References.

Output folder: data/reference/hhrap/2005/
  parameters.csv    long format: cas, name, parameter, value, unit, explanation, reference
  explanations.csv  explanation code -> text (e.g. 60 = Travis and Arms correlation for Br_ag)
  references.csv    reference id -> citation
  metadata.yaml     source file and SHA-256
Values are copied without conversion. The database's toxicity values (2005)
are imported for completeness but ORTAM uses the RSL tables for risk.

The .mdb is read with the pure-Python package access-parser (no MS Access needed).
Fortran analogy: a one-off conversion of a vendor binary file into a plain data deck.

Run:  uv run python -m ortam_core.chemicals.hhrap_import
"""

from __future__ import annotations

import contextlib
import csv
import hashlib
import io
from datetime import date
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "docs" / "references" / "hhrap_2005_chemical_database.mdb"
OUT = ROOT / "data" / "reference" / "hhrap" / "2005"


def _table(db, name: str) -> list[dict]:
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        t = db.parse_table(name)            # access-parser prints harmless overflow notices
    keys = list(t)
    return [{k: t[k][i] for k in keys} for i in range(len(t[keys[0]]))]


def _clean(v):
    if v is None:
        return ""
    if isinstance(v, str):
        return v.strip().replace("\r\n", " ").replace("\n", " ")
    return v


def main() -> None:
    from access_parser import AccessParser

    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        db = AccessParser(str(SOURCE))
    defs = {r["DefRef"]: r for r in _table(db, "Definitions")}
    names: dict[str, str] = {}
    for r in sorted(_table(db, "ChemNames"), key=lambda r: r["SynID"]):
        names.setdefault(r["CAS#"], r["Name"])
    rows = []
    for r in _table(db, "ParamVals"):
        d = defs[r["Param"]]
        rows.append({"cas": r["CAS"], "name": names.get(r["CAS"], ""), "parameter": d["Parameter"],
                     "value": _clean(r["ParamVal"]), "unit": _clean(d["Units"]),
                     "explanation": _clean(r["ParamExp"]), "reference": _clean(r["ParamRef"])})
    rows.sort(key=lambda r: (r["cas"], r["parameter"]))
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "parameters.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with (OUT / "explanations.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["code", "explanation"])
        for r in sorted(_table(db, "ParamExplanations"), key=lambda r: r["NotesRef"]):
            w.writerow([r["NotesRef"], _clean(r["Explanation"])])
    with (OUT / "references.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id", "short", "citation"])
        for r in sorted(_table(db, "References"), key=lambda r: r["RefID"]):
            w.writerow([r["RefID"], _clean(r["ShortCitation"]), _clean(r["LongCitation"])])
    meta = {
        "source": SOURCE.name,
        "document": "HHRAP Companion Database, U.S. EPA 2005 (EPA530-R-05-006)",
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "imported_on": date.today().isoformat(),
        "chemicals": len({r["cas"] for r in rows}),
        "rows": len(rows),
    }
    (OUT / "metadata.yaml").write_text(yaml.safe_dump(meta, sort_keys=False), encoding="utf-8")
    print(f"{meta['chemicals']} chemicals, {len(rows)} values -> {OUT}")


if __name__ == "__main__":
    main()
