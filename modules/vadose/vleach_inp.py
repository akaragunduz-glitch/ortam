"""Read an original VLEACH 2.2a input file (.INP, feet units) into ORTAM vadose.vleach input.

Card layout: VLEACH UG section 6, pp. 15-18, Figure 4. Fixed format:
  card 1  TITLE (A80)            card 2  NPOLY (I3)
  card 3  DELT, STIME, PTIME, PRTIME (4G10.0), years
  card 4  KOC (mL/g), KH (-), CMAX (mg/L), DAIR (m2/day) (4G10.0)
  per polygon:
  card 1  TITLE (A80)
  card 2  AREA (ft2), DELZ (ft), Q (ft/yr), RHOB (g/mL), POR, THETA, FOC (7G10.0)
  card 3  CINF, CATM, CGW (mg/L; negative = impermeable to gas) (3G10.0)
  card 4  NCELL (I5), PLT (A1), PLTIME (G10.0)
  card 5  J1, J2, XCON (ug/kg) (2I5, G10.0), repeated until J2 = NCELL.
Lengths are converted with 1 ft = 0.3048 m.

Run:  uv run python -m modules.vadose.vleach_inp FILE.INP > input.yaml
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml

FT = 0.3048


def _f(s: str) -> float:
    s = s.strip()
    return float(s) if s else 0.0


def _fields(line: str, width: int, n: int) -> list[float]:
    line = line.rstrip("\r\n").ljust(width * n)
    return [_f(line[i * width:(i + 1) * width]) for i in range(n)]


def read_inp(path: str | Path, chemical_name: str = "chemical") -> dict:
    lines = Path(path).read_text(encoding="latin-1").splitlines()
    it = iter(lines)
    title = next(it).strip()
    npoly = int(next(it)[:3])
    delt, stime, ptime, prtime = _fields(next(it), 10, 4)
    koc, kh, cmax, dair = _fields(next(it), 10, 4)
    polygons = []
    for _ in range(npoly):
        name = next(it).strip()
        area, delz, q, rhob, por, theta, foc = _fields(next(it), 10, 7)
        cinf, catm, cgw = _fields(next(it), 10, 3)
        card4 = next(it)
        ncell = int(card4[:5])
        layers, j2 = [], 0
        while j2 < ncell:
            card5 = next(it)
            j1, j2 = int(card5[:5]), int(card5[5:10])
            xcon = _f(card5[10:20])
            layers.append({"top_m": (j1 - 1) * delz * FT, "bottom_m": j2 * delz * FT, "soil_ug_kg": xcon})
        polygons.append({
            "name": name, "area_m2": area * FT * FT, "cell_thickness_m": delz * FT,
            "depth_to_water_m": ncell * delz * FT, "recharge_m_yr": q * FT, "bulk_density_g_cm3": rhob,
            "porosity": por, "water_content": theta, "foc": foc, "recharge_conc_mg_l": cinf,
            "atmosphere_gas_mg_l": catm if catm >= 0 else None,
            "groundwater_mg_l": cgw if cgw >= 0 else None, "initial": layers,
        })
    return {"title": title, "time_step_yr": delt, "simulation_time_yr": stime, "output_interval_yr": ptime,
            "profile_interval_yr": prtime,
            "chemical": {"name": chemical_name, "koc_ml_g": koc, "henry_dimensionless": kh, "solubility_mg_l": cmax,
                         "dair_m2_d": dair},
            "polygons": polygons}


if __name__ == "__main__":
    print(yaml.safe_dump(read_inp(sys.argv[1]), sort_keys=False))
