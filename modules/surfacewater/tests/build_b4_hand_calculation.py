"""Independent hand calculation of HHRAP 2005 Appendix B, Tables B-4-1 to B-4-25, as a LibreOffice workbook.

Not part of the test run. Run with the system python3 that can `import uno`:
    python3 modules/surfacewater/tests/build_b4_hand_calculation.py
Writes data/reference/hhrap/2005/b4_hand_calculation.ods (formulas, one column per case) and
data/reference/hhrap/2005/b4_hand_calculation.json (values recalculated by LibreOffice).

The worksheet formulas are transcribed directly from the equations printed in the tables
(hhrap_2005_appendix_b.pdf; table and printed page [PDF page] in column C), independently of the
Python module modules/surfacewater/. Chemical parameters are copied from the HHRAP Companion
Database (data/reference/hhrap/2005/parameters.csv) into the input rows. Units are those of the
tables. Row formulas use {name} for the cell of row `name` in the same column.
"""
from __future__ import annotations

import csv, json, os, subprocess, time
from pathlib import Path

import uno  # noqa: F401
from com.sun.star.beans import PropertyValue

ROOT = Path(__file__).resolve().parents[3]
OUT_ODS = ROOT / "data/reference/hhrap/2005/b4_hand_calculation.ods"
OUT_JSON = ROOT / "data/reference/hhrap/2005/b4_hand_calculation.json"


def db(cas):
    keys = {"Kds": "kds", "Kdsw": "kdsw", "Kdbs": "kdbs", "H": "h", "Dw": "dw", "Da": "da", "ksg": "ksg", "Fv": "fv"}
    out = {}
    with (ROOT / "data/reference/hhrap/2005/parameters.csv").open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["cas"] == cas and r["parameter"] in keys:
                out[keys[r["parameter"]]] = float(r["value"])
    return out


BENZ, CD, TCDD, NAPH = db("71-43-2"), db("7440-43-9"), db("1746-01-6"), db("91-20-3")
HGCL2, MEHG = db("7487-94-7"), db("22967-92-6")

COMMON = dict(zs=2, bd=1.5, tsw=0.2, rhos=2.7, ta=298, R=8.205e-5, irr=0, kse_mode=0, k_usle=0.39, ls=1.5,
              c_usle=0.1, pf=1.0, b_sd=0.125, dbs=0.03, tss=10, cbs=1.0, tbs=0.6, twk=298, theta=1.026,
              cd=0.0011, wind=3.9, vk=0.4, lz=4, rho_a=0.0012, rho_w=1.0, mu_a=1.81e-4, mu_w=1.69e-2,
              q=0, dytwv_ws=0, dytwp_ws=0, dytwv_wb=0, dytwp_wb=0, cywv=0, td=30, t1=0, t2=30, l_ext=0,
              basis=0, hg_factor=1)
RIVER = dict(river=1, aw=5.0e4, vfx=3.0e7, dwc=1.5, u=0.5)
LAKE = dict(river=0, aw=2.0e5, vfx=1.0e6, dwc=4.0, u=0)
WATERSHED = dict(al=5.0e6, ai=2.0e5, ro=15, p=100, ev=60, rf=150)
AIR = dict(q=1.0e-4, dytwv_ws=0.05, dytwp_ws=0.4, dytwv_wb=0.06, dytwp_wb=0.5, cywv=0.8)
CASES = [
    ("river_measured_organic", {**COMMON, **RIVER, **WATERSHED, **BENZ, "er": 3, "cs_in": 5.0}),
    ("lake_measured_inorganic", {**COMMON, **LAKE, **WATERSHED, **CD, "er": 1, "cs_in": 20.0, "al": 8.0e7, "ai": 0}),
    ("river_deposition_t2_le_td", {**COMMON, **RIVER, **WATERSHED, **TCDD, **AIR, "er": 3, "cs_in": None}),
    ("lake_deposition_t2_gt_td_ext", {**COMMON, **LAKE, **WATERSHED, **NAPH, **AIR, "er": 3, "cs_in": None,
                                      "t2": 40, "l_ext": 50.0, "zs": 20}),
    ("river_negative_kb_kse_irr", {**COMMON, **RIVER, **WATERSHED, **BENZ, "er": 3, "cs_in": 2.0, "vfx": 5.0e9,
                                   "kse_mode": 1, "irr": 20, "al": 2.0e5, "ai": 0}),
    ("lake_deposition_cstd_basis", {**COMMON, **LAKE, **WATERSHED, **NAPH, **AIR, "er": 3, "cs_in": None, "basis": 1}),
    # Mercury: two soil species columns and one water body column (Tables B-4-1, -8, -9, -10, -11, -12, -24)
    ("hg_soil_hg2", {**COMMON, **LAKE, **WATERSHED, **HGCL2, **AIR, "er": 1, "cs_in": None, "fv": 0.85,
                     "hg_factor": 0.48 * 0.98}),
    ("hg_soil_mhg", {**COMMON, **LAKE, **WATERSHED, **MEHG, **AIR, "er": 1, "cs_in": None, "fv": 0.85,
                     "hg_factor": 0.48 * 0.02}),
    ("hg_water_body", {**COMMON, **LAKE, **WATERSHED, **HGCL2, **AIR, "er": 1, "cs_in": None, "fv": 0.85,
                       "hg_factor": 0.48}),
]
INPUTS = [  # key, unit, source
    ("river", "1 river / 0 lake", "Tables B-4-20, B-4-21"), ("cs_in", "mg/kg", "measured Cs (blank: Table B-4-1)"),
    ("basis", "0 Cs / 1 CstD", "Table B-4-1 description p. B-168 [171]"),
    ("kds", "mL/g", "App. A-2 database"), ("kdsw", "L/kg", "database"), ("kdbs", "L/kg", "database"),
    ("h", "atm-m3/mol", "database"), ("dw", "cm2/s", "database"), ("da", "cm2/s", "database"),
    ("ksg", "1/yr", "database"), ("fv", "-", "database (Hg: 0.85)"), ("er", "-", "B-4-3 p. B-181 [184]"),
    ("zs", "cm", "B-4-1 p. B-172 [175]"), ("bd", "g/cm3", "B-4-1 p. B-172 [175]"), ("tsw", "-", "B-4-3 p. B-182 [185]"),
    ("rhos", "g/cm3", "B-4-6 p. B-196 [199]"), ("ta", "K", "B-4-6 p. B-196 [199]"), ("R", "atm-m3/mol-K", "B-4-6"),
    ("p", "cm/yr", "B-4-5"), ("irr", "cm/yr", "B-4-5"), ("ro", "cm/yr", "B-4-4"), ("ev", "cm/yr", "B-4-5"),
    ("kse_mode", "0 default / 1 B-4-3", "B-4-2 p. B-177 [180]"),
    ("rf", "1/yr", "B-4-13"), ("k_usle", "ton/acre", "B-4-13"), ("ls", "-", "B-4-13"), ("c_usle", "-", "B-4-13"),
    ("pf", "-", "B-4-13"), ("al", "m2", "B-4-10"), ("ai", "m2", "B-4-9"), ("b_sd", "-", "B-4-14"),
    ("q", "g/s", "B-4-1"), ("hg_factor", "-", "B-4-1/-8/-9/-12: 0.48 (x 0.98 / 0.02 soil split)"),
    ("dytwv_ws", "s/m2-yr", "watershed"), ("dytwp_ws", "s/m2-yr", "watershed"),
    ("dytwv_wb", "s/m2-yr", "water body"), ("dytwp_wb", "s/m2-yr", "water body"), ("cywv", "ug-s/g-m3", "water body"),
    ("td", "yr", "B-4-1"), ("t1", "yr", "B-4-1"), ("t2", "yr", "B-4-1"),
    ("aw", "m2", "B-4-8"), ("vfx", "m3/yr", "B-4-15"), ("dwc", "m", "B-4-15"), ("dbs", "m", "B-4-15"),
    ("tss", "mg/L", "B-4-16"), ("cbs", "g/cm3", "B-4-16"), ("tbs", "-", "B-4-16"), ("twk", "K", "B-4-19"),
    ("theta", "-", "B-4-19"), ("u", "m/s", "B-4-20"), ("cd", "-", "B-4-20"), ("wind", "m/s", "B-4-20"),
    ("vk", "-", "B-4-20"), ("lz", "-", "B-4-20"), ("rho_a", "g/cm3", "B-4-20"), ("rho_w", "g/cm3", "B-4-20"),
    ("mu_a", "g/cm-s", "B-4-21"), ("mu_w", "g/cm-s", "B-4-20"), ("l_ext", "g/yr", "ORTAM extension: external load"),
]
F = [  # key, unit, source, formula (LibreOffice syntax with ; separators)
    ("qe", "g/s", "Q (x 0.48 x split for Hg)", "={q}*{hg_factor}"),
    ("xe", "kg/m2-yr", "B-4-13 p. B-219 [222]", "={rf}*{k_usle}*{ls}*{c_usle}*{pf}*907.18/4047"),
    ("a_sd", "-", "B-4-14 p. B-223 [226]",
     "=IF({al}/2590000<=0.1;2.1;IF({al}/2590000<=1;1.9;IF({al}/2590000<=10;1.4;IF({al}/2590000<=100;1.2;0.6))))"),
    ("sd", "-", "B-4-14 p. B-223 [226]", "={a_sd}*{al}^(-{b_sd})"),
    ("kse", "1/yr", "B-4-3 p. B-180 [183]",
     "=IF({kse_mode}=1;0.1*{xe}*{sd}*{er}/({bd}*{zs})*({kds}*{bd}/({tsw}+{kds}*{bd}));0)"),
    ("ksr", "1/yr", "B-4-4 p. B-185 [188]", "={ro}/({tsw}*{zs})*(1/(1+({kds}*{bd}/{tsw})))"),
    ("ksl", "1/yr", "B-4-5 p. B-189 [192]", "=({p}+{irr}-{ro}-{ev})/({tsw}*{zs}*(1+({bd}*{kds}/{tsw})))"),
    ("ksv", "1/yr", "B-4-6 p. B-194 [197]",
     "=(31536000*{h}/({zs}*{kds}*{R}*{ta}*{bd}))*({da}/{zs})*(1-({bd}/{rhos})-{tsw})"),
    ("ks", "1/yr", "B-4-2 p. B-176 [179]", "={ksg}+{kse}+{ksr}+{ksl}+{ksv}"),
    ("ds", "mg/kg-yr", "B-4-1 p. B-169 [172]", "=100*{qe}/({zs}*{bd})*({fv}*{dytwv_ws}+{dytwp_ws}*(1-{fv}))"),
    ("cstd", "mg/kg", "B-4-1 p. B-169 [172]", "={ds}*(1-EXP(-{ks}*{td}))/{ks}"),
    ("cs_dep", "mg/kg", "B-4-1 p. B-168 [171]",
     "=IF({t2}<={td};{ds}/({ks}*({td}-{t1}))*(({td}+EXP(-{ks}*{td})/{ks})-({t1}+EXP(-{ks}*{t1})/{ks}));"
     "((({ds}*{td}-{cstd})/{ks})+({cstd}/{ks}*(1-EXP(-{ks}*({t2}-{td})))))/({t2}-{t1}))"),
    ("cs", "mg/kg", "measured, or Cs / CstD (basis)", "=IF(ISBLANK({cs_in});IF({basis}=1;{cstd};{cs_dep});{cs_in})"),
    ("ldep", "g/yr", "B-4-8 p. B-201 [204]", "={qe}*({fv}*{dytwv_wb}+(1-{fv})*{dytwp_wb})*{aw}"),
    ("lri", "g/yr", "B-4-9 p. B-204 [207]", "={qe}*({fv}*{dytwv_ws}+(1-{fv})*{dytwp_ws})*{ai}"),
    ("lr", "g/yr", "B-4-10 p. B-207 [210]", "={ro}*({al}-{ai})*({cs}*{bd}/({tsw}+{kds}*{bd}))*0.01"),
    ("le", "g/yr", "B-4-11 p. B-211 [214]",
     "={xe}*({al}-{ai})*{sd}*{er}*({cs}*{kds}*{bd}/({tsw}+{kds}*{bd}))*0.001"),
    ("dz", "m", "B-4-16 p. B-232 [235]", "={dwc}+{dbs}"),
    ("kl", "m/yr", "B-4-20 p. B-244 [247]",
     "=IF({river}=1;SQRT(0.0001*{dw}*{u}/{dz})*31536000;"
     "({cd}^0.5*{wind})*({rho_a}/{rho_w})^0.5*({vk}^0.33/{lz})*({mu_w}/({rho_w}*{dw}))^(-0.67)*31536000)"),
    ("kg", "m/yr", "B-4-21 p. B-249 [252]",
     "=IF({river}=1;36500;({cd}^0.5*{wind})*({vk}^0.33/{lz})*({mu_a}/({rho_a}*{da}))^(-0.67)*31536000)"),
    ("kv_overall", "m/yr", "B-4-19 p. B-241 [244]", "=((1/{kl})+1/({kg}*{h}/({R}*{twk})))^(-1)*{theta}^({twk}-293)"),
    ("ldif", "g/yr", "B-4-12 p. B-215 [218]", "={kv_overall}*{qe}*{fv}*{cywv}*{aw}*0.000001/({h}/({R}*{twk}))"),
    ("lr_sum", "g/yr", "B-4-10 (Hg: Hg2+ + MHg)", "={lr}"),
    ("le_sum", "g/yr", "B-4-11 (Hg: Hg2+ + MHg)", "={le}"),
    ("lt", "g/yr", "B-4-7 p. B-199 [202] + external load", "={ldep}+{ldif}+{lri}+{lr_sum}+{le_sum}+{l_ext}"),
    ("fwc", "-", "B-4-16 p. B-230 [233]",
     "=(1+{kdsw}*{tss}*0.000001)*{dwc}/{dz}/((1+{kdsw}*{tss}*0.000001)*{dwc}/{dz}+({tbs}+{kdbs}*{cbs})*{dbs}/{dz})"),
    ("fbs", "-", "B-4-16 p. B-230 [233]", "=1-{fwc}"),
    ("kv", "1/yr", "B-4-18 p. B-237 [240]", "={kv_overall}/({dz}*(1+{kdsw}*{tss}*0.000001))"),
    ("kb_raw", "1/yr", "B-4-22 p. B-253 [256]",
     "=(({xe}*{al}*{sd}*1000-{vfx}*{tss})/({aw}*{tss}))*({tss}*0.000001/({cbs}*{dbs}))"),
    ("kb", "1/yr", "Chapter 5 p. 5-81 [82]: negative -> 0", "=MAX(0;{kb_raw})"),
    ("kwt", "1/yr", "B-4-17 p. B-235 [238]", "={fwc}*{kv}+{fbs}*{kb}"),
    ("cwtot", "mg/L", "B-4-15 p. B-226 [229]", "={lt}/({vfx}*{fwc}+{kwt}*{aw}*({dwc}+{dbs}))"),
    ("cwctot", "mg/L", "B-4-23 p. B-257 [260]", "={fwc}*{cwtot}*({dwc}+{dbs})/{dwc}"),
    ("cdw", "mg/L", "B-4-24 p. B-260 [263]", "={cwctot}/(1+{kdsw}*{tss}*0.000001)"),
    ("csb", "mg/kg", "B-4-25 p. B-263 [266]", "={fbs}*{cwtot}*({kdbs}/({tbs}+{kdbs}*{cbs}))*({dwc}+{dbs})/{dbs}"),
]
HG_OVERRIDES = {"lr_sum": "={lr@hg_soil_hg2}+{lr@hg_soil_mhg}", "le_sum": "={le@hg_soil_hg2}+{le@hg_soil_mhg}"}
HG_EXTRA = [("cdw_hg2", "mg/L", "B-4-24 p. B-260 [263]: 0.85 Cdw", "=0.85*{cdw}"),
            ("cdw_mhg", "mg/L", "B-4-24 p. B-260 [263]: 0.15 Cdw", "=0.15*{cdw}")]


def pv(n, v):
    p = PropertyValue(); p.Name = n; p.Value = v
    return p


def col(i):
    return chr(ord("D") + i)


def main():
    proc = subprocess.Popen(["soffice", "--headless", "--invisible", "--norestore",
                             "--accept=socket,host=localhost,port=2005;urp;"],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    local = uno.getComponentContext()
    res = local.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local)
    for _ in range(90):
        try:
            ctx = res.resolve("uno:socket,host=localhost,port=2005;urp;StarOffice.ComponentContext"); break
        except Exception:
            time.sleep(1)
    desk = ctx.ServiceManager.createInstanceWithContext("com.sun.star.frame.Desktop", ctx)
    doc = desk.loadComponentFromURL("private:factory/scalc", "_blank", 0, (pv("Hidden", True),))
    sh = doc.Sheets.getByIndex(0); sh.Name = "B-4 hand calculation"
    rows = INPUTS + [(k, u, s) for k, u, s, _ in F + HG_EXTRA]
    rowof = {k: i + 2 for i, (k, _, _) in enumerate(rows)}          # 1-based row numbers
    casecol = {name: col(j) for j, (name, _) in enumerate(CASES)}
    hdr = ["symbol", "unit", "source"] + [n for n, _ in CASES]
    for j, h in enumerate(hdr):
        sh.getCellByPosition(j, 0).setString(h)
    for i, (k, u, s) in enumerate(rows):
        for j, v in enumerate((k, u, s)):
            sh.getCellByPosition(j, i + 1).setString(v)

    def expand(formula, c):
        import re
        def rep(m):
            key = m.group(1)
            if "@" in key:
                key, case = key.split("@")
                return f"{casecol[case]}{rowof[key]}"
            return f"{c}{rowof[key]}"
        return re.sub(r"\{([A-Za-z0-9_@]+)\}", rep, formula)

    for j, (name, vals) in enumerate(CASES):
        c = col(j)
        for k, _, _ in INPUTS:
            v = vals.get(k)
            cell = sh.getCellRangeByName(f"{c}{rowof[k]}")
            if v is not None:
                cell.setValue(float(v))
        for k, _, _, f in F:
            ff = HG_OVERRIDES.get(k, f) if name == "hg_water_body" else f
            sh.getCellRangeByName(f"{c}{rowof[k]}").setFormula(expand(ff, c))
        if name == "hg_water_body":
            for k, _, _, f in HG_EXTRA:
                sh.getCellRangeByName(f"{c}{rowof[k]}").setFormula(expand(f, c))
    doc.calculateAll()
    out = {"source": "hhrap_2005_appendix_b.pdf, Tables B-4-1 to B-4-25 (formulas in b4_hand_calculation.ods)",
           "units": {k: u for k, u, _ in rows}, "cases": {}}
    errors = []
    for j, (name, vals) in enumerate(CASES):
        c = col(j); rec = {}
        for k, _, _ in rows:
            cell = sh.getCellRangeByName(f"{c}{rowof[k]}")
            if cell.getError():
                errors.append((name, k, cell.getError()))
            elif cell.getType().value != "EMPTY":
                rec[k] = cell.getValue()
        out["cases"][name] = {"inputs": {k: vals.get(k) for k, _, _ in INPUTS}, "values": rec}
    doc.storeToURL(uno.systemPathToFileUrl(str(OUT_ODS)), (pv("FilterName", "calc8"), pv("Overwrite", True)))
    doc.close(True)
    try:
        desk.terminate()
    except Exception:
        pass
    proc.wait(timeout=60)
    OUT_JSON.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("errors:", errors)
    for name in out["cases"]:
        v = out["cases"][name]["values"]
        print(name, {k: f"{v.get(k, float('nan')):.4g}" for k in ("ks", "cs", "lt", "kb_raw", "fwc", "cwtot", "cdw", "csb")})


if __name__ == "__main__":
    main()
