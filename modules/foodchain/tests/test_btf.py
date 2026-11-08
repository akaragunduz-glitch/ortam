"""foodchain.btf against the HHRAP 2005 Companion Database.

The database records, for every value, how it was derived (explanation code; see
data/reference/hhrap/2005/explanations.csv). For the codes that name an Appendix A-2
correlation, the value is recomputed from the database's own log Kow (and H, Koc) and
compared. The database stores values rounded to 3 significant figures: tolerance 1 %.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import pytest

from modules.foodchain import btf

DATA = Path(__file__).resolve().parents[3] / "data" / "reference" / "hhrap" / "2005"


def _db():
    out = defaultdict(dict)
    with (DATA / "parameters.csv").open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            out[r["cas"]][r["parameter"]] = (float(r["value"]) if r["value"] else None, r["explanation"])
    return out


DB = _db()

# parameter, explanation code (correlation used), function of the database row, minimum count
CHECKS = [
    ("Brag", "60", lambda c: btf.br_aboveground(c["Kow (log)"][0]), 180),
    ("Brf/s", "58", lambda c: btf.br_aboveground(c["Kow (log)"][0]), 180),
    ("Brgrain", "56", lambda c: btf.br_aboveground(c["Kow (log)"][0]), 180),
    ("RCF", "42", lambda c: btf.rcf_dw(c["Kow (log)"][0]), 180),
    # Kds values derived from Koc (explanation 48) are stored rounded; Br_rootveg was computed
    # with the unrounded Kds = Koc * 0.01 (foc = 0.01, explanation 2)
    ("Brroot", "55", lambda c: btf.br_rootveg(c["RCF"][0], c["Koc"][0] * 0.01 if c["Kds"][1] == "48"
                                              else c["Kds"][0]), 180),
    ("Bvag", "53", lambda c: btf.bv(c["Kow (log)"][0], c["H"][0]), 160),
    ("Bvforage", "53", lambda c: btf.bv(c["Kow (log)"][0], c["H"][0]), 160),
    ("Babeef", "77", lambda c: btf.ba(c["Kow (log)"][0], "beef"), 170),
    ("Bamilk", "79", lambda c: btf.ba(c["Kow (log)"][0], "milk"), 170),
    ("Bapork", "75", lambda c: btf.ba(c["Kow (log)"][0], "pork"), 170),
    ("Bapoult", "76", lambda c: btf.ba(c["Kow (log)"][0], "chicken"), 170),
    ("Baeggs", "96", lambda c: btf.ba(c["Kow (log)"][0], "eggs"), 170),
]
# Values that differ by more than the rounding (documented in docs/methods/foodchain.md)
KNOWN = {("Brf/s", "118-74-1"),        # hexachlorobenzene: 0.033 printed, 0.0335 computed (2-digit rounding)
         ("RCF", "39001-02-0"),        # OCDF: database RCF corresponds to log Kow 8.2; its log Kow is 8.0
         ("Brroot", "75-34-3")}        # 1,1-dichloroethane: Kds stored as 0.06; Br_rootveg 189 implies Kds 0.0635


@pytest.mark.parametrize("param,code,fn,nmin", CHECKS, ids=[c[0] for c in CHECKS])
def test_correlation_matches_database(param, code, fn, nmin):
    n, misses = 0, []
    for cas, c in DB.items():
        v, e = c.get(param, (None, ""))
        if e != code or not v or c.get("Kow (log)", (None,))[0] is None:
            continue
        n += 1
        got = fn(c)
        # Br_rootveg: quotient of two rounded table values (RCF, Kds of 1-2 digits) -> 5 %
        tol = 0.05 if param == "Brroot" else 0.01
        if abs(got / v - 1) > tol and (param, cas) not in KNOWN:
            misses.append((cas, v, got))
    assert n >= nmin, n
    assert not misses, misses[:5]


def test_benzene_by_hand():
    # log Kow 2.1 (database): Br = 10^(1.588 - 0.578*2.1) = 2.367; database 2.37
    assert btf.br_aboveground(2.1) == pytest.approx(2.367, rel=1e-3)
    # RCF = 10^(0.77*2.1 - 1.52)/(1 - 0.87) = 9.62
    assert btf.rcf_dw(2.1) == pytest.approx(9.62, rel=1e-3)
    # Ba_beef = 10^(-0.099*2.1^2 + 1.07*2.1 - 3.56) * 0.19 = 3.382e-3
    assert btf.ba(2.1, "beef") == pytest.approx(3.382e-3, rel=1e-3)
