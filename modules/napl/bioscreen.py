"""BIOSCREEN 1.4 engine: Domenico (1987) model with source decay, first-order decay and
instantaneous reaction (electron acceptor superposition), as implemented in BIOSCRN4.xls.

Units inside this file are BIOSCREEN's: ft, yr, mg/L, kg, L (the module napl.source-plume
converts SI input at the boundary, as a Fortran driver converts units before CALLing a solver).

Sources (printed page [PDF page]; full references in docs/methods/napl.md):
  BSCR    Newell, McLeod, Gonzales (1996), BIOSCREEN User's Manual 1.3, EPA/600/R-96/087
          (bioscreen_users_manual_1996.pdf, PDF pages 1-63)
  BSCR14  BIOSCREEN Version 1.4 Revisions, July 1997 (same PDF, pages 64-85)
  XLS     BIOSCRN4.xls (bioscreen_1.4.zip): worksheet formulas, cited by sheet!cell; the error
          function gerf of its VBA "Macro Module".

Two settings (Settings.compatible):
  True   reproduce BIOSCREEN 1.4 exactly (verified against the workbook to 1e-9), including
         its constants (28.3 L/ft3, 0.693, 3.28 ft/m, the gerf approximation) and three workbook
         inconsistencies (documented in docs/methods/napl.md §4.6):
           Q1 instantaneous model: the middle source layer decays with the source rate of the
              no-degradation model (ks) instead of ks_inst (XLS 'Centerline Output'!B80,
              'Plume Output'!C113: the zone-2 denominator uses ksourceinst);
           Q2 first-order plume mass omits the retardation factor ('Plume Output'!P104), while
              the other plume masses include it (P91, P116, P121);
           Q3 the biodegradation capacity is added to a source layer whenever the next outer
              layer increment is zero ('Input'!T46-T47), so it can be added twice.
  False  consistent form: exact erf, 0.3048 m/ft, 28.316846592 L/ft3, ln 2, and Q1-Q3 corrected.

Source depletion of the instantaneous model (Case.instantaneous_source):
  "bioscreen"     BIOSCREEN 1.4: the effective source C0 + BC decays as e^(-ks_inst t),
                  ks_inst = Q (Cavg + BC)/M0, and BC is subtracted (BSCR p. 35 [41]); the observed
                  source concentration reaches zero while most of the mass is still in the source.
  "mass-balance"  (default) BC is constant (property of the upgradient groundwater, BSCR p. 39 [45]);
                  the proportionality between source concentration and source mass (BSCR App. A.3,
                  p. 51 [57]) applies to the observed concentration only. With BSCR App. A.3 Eq. 4,
                  dM/dt = -Q (Cavg M/M0 + BC), so M = M0 f(t), f = (1 + b) e^(-ks t) - b, b = BC/Cavg,
                  ks = Q Cavg/M0, every zone concentration is Co,k f(t), and source mass and source
                  concentration reach zero together at t_end = ln((1 + b)/b)/ks (derivation:
                  docs/methods/napl.md §4.8). Equal to "bioscreen" for BC = 0 or an infinite source.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.special import erf as _erf_exact

X_SOURCE_FT = 1e-4          # XLS 'Centerline Output'!B35, 'Plume Output'!C57: x of the source column
ALPHA_Z_ZERO_FT = 1e-99     # XLS Input!E73: alpha_z = 0 is replaced by 1e-99 ft
L_PER_FT3_BIOSCREEN = 28.3  # XLS Input!F59, 'Plume Output'!P91, B11
L_PER_FT3_EXACT = 28.316846592
FT2_PER_ACRE = 43560.0      # XLS 'Plume Output'!M30-M31 (acre-ft)
DAYS_PER_YEAR = 365.0       # XLS Input!E5 (86400*365 s/yr), 'Plume Output'!B11 (/365)

# Utilization factors, mg electron acceptor (or by-product) per mg BTEX: BSCR p. 25-27 [31-33]
# (O2 3.14, NO3 4.9, Fe2+ 21.8, SO4 4.7, CH4 0.78); XLS Input!E47-E51.
UTILIZATION_FACTORS = {"oxygen": 3.14, "nitrate": 4.9, "ferrous_iron": 21.8, "sulfate": 4.7, "methane": 0.78}
# Sign of the change in mass: electron acceptors are consumed (-), by-products produced (+)
# (XLS 'Plume Output'!I25-M25).
EA_SIGN = {"oxygen": -1, "nitrate": -1, "ferrous_iron": 1, "sulfate": -1, "methane": 1}


def erf_bioscreen(x):
    """Error function of BIOSCREEN 1.4 (XLS 'Macro Module', Function gerf; the Chebyshev
    approximation of erfc with fractional error < 1.2e-7). Vectorized."""
    x = np.asarray(x, dtype=float)
    g = np.abs(x)
    t = 1.0 / (1.0 + 0.5 * g)
    a = t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))
    a = t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + a)))
    a = t * np.exp(-g * g - 1.26551223 + t * (1.00002368 + t * (0.37409196 + a)))
    return np.where(x >= 0.0, 1.0 - a, -1.0 + a)


@dataclass(frozen=True)
class Settings:
    compatible: bool = True

    @property
    def l_per_ft3(self) -> float:
        return L_PER_FT3_BIOSCREEN if self.compatible else L_PER_FT3_EXACT

    @property
    def ln2(self) -> float:
        # BSCR p. 23 [29]: k = 0.693 / t1/2
        return 0.693 if self.compatible else math.log(2.0)

    @property
    def ft_per_m(self) -> float:
        # BSCR p. 19 [25], Xu and Eckstein (1995) written for Lp in ft with 3.28 ft/m
        return 3.28 if self.compatible else 1.0 / 0.3048

    def erf(self, x):
        return erf_bioscreen(x) if self.compatible else _erf_exact(np.asarray(x, dtype=float))


# ---------------------------------------------------------------------------------------------
# Input relations (BSCR data entry section)
# ---------------------------------------------------------------------------------------------
def seepage_velocity_ft_yr(k_cm_s: float, gradient: float, porosity: float) -> float:
    """Vs = K i / n. BSCR p. 16 [22]; XLS Input!E5 (= K*86400*365*i/2.54/12/n, K in cm/s)."""
    return k_cm_s * 86400.0 * DAYS_PER_YEAR * gradient / 2.54 / 12.0 / porosity


def alpha_x_xu_eckstein_ft(plume_length_ft: float, s: Settings) -> float:
    """alpha_x = 3.28 * 0.83 * (log10(Lp/3.28))^2.414 (Lp in ft), Xu and Eckstein (1995).
    BSCR p. 19 [25] and Appendix A.4 p. 53 [59]; XLS Input!E12."""
    return s.ft_per_m * 0.83 * math.log10(plume_length_ft / s.ft_per_m) ** 2.414


def alpha_y_default(alpha_x: float) -> float:
    """alpha_y = 0.10 alpha_x. BSCR p. 19 [25]; XLS Input!E13 (= alpha_x/10)."""
    return alpha_x / 10.0


def retardation_factor(rho_kg_l: float, koc_l_kg: float, foc: float, porosity: float) -> float:
    """R = 1 + Kd rho_b / n, Kd = Koc foc. BSCR p. 21 [27]; XLS Input!E19."""
    return 1.0 + rho_kg_l / porosity * koc_l_kg * foc


def decay_rate_per_yr(half_life_yr: float, s: Settings) -> float:
    """lambda = 0.693 / t1/2. BSCR p. 23 [29]; XLS Input!E26."""
    return s.ln2 / half_life_yr


# ---------------------------------------------------------------------------------------------
# Case description and source terms
# ---------------------------------------------------------------------------------------------
@dataclass
class Case:
    """One BIOSCREEN run in BIOSCREEN units. zones: [(Yo1, Co1), (Yo2, Co2), (Yo3, Co3)] with
    BIOSCREEN's numbering (zone 3 central, zones 2 and 1 on both sides; Yo1/Yo2 are the widths of
    each side band, BSCR p. 30 [36]); a missing zone has width 0 and concentration 0."""
    vs_ft_yr: float
    porosity: float
    alpha_x_ft: float
    alpha_y_ft: float
    alpha_z_ft: float
    retardation: float
    lambda_per_yr: float
    length_ft: float
    width_ft: float
    time_yr: float
    thickness_ft: float
    zones: list[tuple[float, float]]
    mass_kg: float | None = None                      # None: infinite source
    electron_acceptors: dict[str, float] = field(default_factory=dict)   # mg/L, keys of UTILIZATION_FACTORS
    utilization_factors: dict[str, float] = field(default_factory=lambda: dict(UTILIZATION_FACTORS))
    instantaneous_source: str = "mass-balance"          # or "bioscreen"


@dataclass
class Source:
    widths_ft: tuple[float, float, float]        # Y1, Y2, Y3: full widths of the layers (XLS Input!L45-L47)
    increments_mg_l: tuple[float, float, float]  # C1, C2, C3: layer increments (XLS Input!Q45-Q47)
    bioc_mg_l: tuple[float, float, float]        # layer increments of the instantaneous model (XLS Input!T45-T47)
    biodegradation_capacity_mg_l: float          # BC (XLS Input!F52)
    ea_share: dict[str, float]                   # (C_ea/UF)/BC (XLS Input!H47-H51)
    source_width_ft: float                       # max(Y1, Y2, Y3)
    flow_l_yr: float                             # Q (XLS Input!F59)
    c_avg_mg_l: float                            # width-weighted source concentration (Input!F60)
    flux_kg_yr: float                            # Q * Cavg (Input!F61)
    ks_per_yr: float                             # source decay, no degradation / first order (Input!F62)
    c_avg_inst_mg_l: float                       # Cavg + BC (Input!F65)
    flux_inst_kg_yr: float                       # (Input!F66)
    ks_inst_per_yr: float                        # source decay, instantaneous reaction (Input!F67)
    half_life_yr: float                          # 0.693/ks (Input!F63), inf for an infinite source
    half_life_inst_yr: float                     # 0.693/ks_inst (Input!F68)
    b_ratio: float = 0.0                         # BC/Cavg (mass-balance source)
    depletion_time_inst_yr: float = math.inf     # t_end of the mass-balance source

    def remaining_fraction(self, t):
        """f(t) = M(t)/M0 of the mass-balance source of the instantaneous model (docs/methods/napl.md §4.8)."""
        t = np.asarray(t, dtype=float)
        if self.ks_per_yr == 0:
            return np.ones_like(t)
        return np.maximum(0.0, (1.0 + self.b_ratio) * np.exp(-self.ks_per_yr * t) - self.b_ratio)


def source_terms(c: Case, s: Settings) -> Source:
    (yo1, co1), (yo2, co2), (yo3, co3) = c.zones
    # Layer-cake source: three superimposed Domenico sources (BSCR p. 47 [53], p. 30 [36]).
    y1 = 0.0 if yo1 == 0 else 2 * yo1 + 2 * yo2 + yo3       # XLS Input!L45
    y2 = 0.0 if yo2 == 0 else 2 * yo2 + yo3                 # Input!L46
    y3 = yo3                                                # Input!L47
    c1 = 0.0 if co1 == 0 else co1                           # Input!Q45
    c2 = 0.0 if co2 == 0 else co2 - co1                     # Input!Q46
    c3 = 0.0 if co3 == 0 else co3 - co2                     # Input!Q47
    # Biodegradation capacity BC = sum(C_ea,n / UF_n), BSCR14 p. 8 [72]; XLS Input!F47-F52.
    terms = {k: c.electron_acceptors.get(k, 0.0) / c.utilization_factors[k] for k in UTILIZATION_FACTORS}
    bc = sum(terms.values())
    share = {k: (v / bc if bc > 0 else 0.0) for k, v in terms.items()}
    # Effective source concentration C0 + BC for the instantaneous model (BSCR p. 47 [53]):
    # BC enters the outermost layer only, since the layers are increments.
    if s.compatible:   # XLS Input!T45-T47 (Q3: BC again if the next outer increment is zero)
        b1 = 0.0 if c1 == 0 else c1 + bc
        b2 = 0.0 if c2 == 0 else (c2 + bc if c1 == 0 else c2)
        b3 = 0.0 if c3 == 0 else (c3 + bc if c2 == 0 else c3)
    else:
        inc = [c1, c2, c3]
        first = next((k for k, v in enumerate(inc) if v != 0), None)
        b1, b2, b3 = (v + bc if k == first else v for k, v in enumerate(inc))
    w = max(y1, y2, y3)
    if w <= 0:
        raise ValueError("Source width is zero.")
    # Q = Vs n Z W (Darcy velocity x source area), BSCR p. 41 [47]; XLS Input!F59 (L/yr).
    q = c.vs_ft_yr * c.porosity * c.thickness_ft * w * s.l_per_ft3
    cavg = (yo1 * co1 * 2 + yo2 * co2 * 2 + yo3 * co3) / w   # Input!F60
    flux = q * cavg / 1e6                                       # Input!F61, kg/yr
    cavg_i = cavg + bc
    flux_i = q * cavg_i / 1e6
    # Source decay ks = Q C0 / M0, BSCR Appendix A.3 p. 50-51 [56-57]; XLS Input!F62, F67.
    if c.mass_kg is None:
        ks = ksi = 0.0
    else:
        ks, ksi = flux / c.mass_kg, flux_i / c.mass_kg
    hl = math.inf if ks == 0 else s.ln2 / ks
    hli = math.inf if ks == 0 else s.ln2 / ksi     # Input!F68 tests F62 (ks), not F67
    if c.instantaneous_source not in ("bioscreen", "mass-balance"):
        raise ValueError("instantaneous_source must be 'bioscreen' or 'mass-balance'.")
    b = bc / cavg if cavg > 0 else 0.0
    t_end = math.inf if (ks == 0 or b == 0) else math.log((1.0 + b) / b) / ks
    return Source((y1, y2, y3), (c1, c2, c3), (b1, b2, b3), bc, share, w, q, cavg, flux, ks, cavg_i, flux_i,
                  ksi, hl, hli, b, t_end)


# ---------------------------------------------------------------------------------------------
# Concentrations
# ---------------------------------------------------------------------------------------------
@dataclass
class Concentrations:
    no_decay: np.ndarray
    first_order: np.ndarray
    inst_before: np.ndarray     # no-decay plume of the effective source (C0 + BC), before subtracting BC
    instantaneous: np.ndarray


def concentrations(c: Case, src: Source, s: Settings, x_ft, y_ft, t_yr=None) -> Concentrations:
    """C(x, y, z=0, t) of the three models; x, y broadcast against each other.

    Domenico model with source decay (BSCR14 p. 8 [72]; BSCR Appendix A.1 p. 46 [52]):
      C = sum_k C_k exp(-ks (t - x/v)) / 8 * Fx * Fy_k * Fz
      Fx (no decay) = erfc((x - v t) / (2 (alpha_x v t)^0.5))
      Fx (first order) = exp(x/(2 alpha_x) (1 - (1 + 4 lambda alpha_x / v)^0.5))
                         * erfc((x - v t (1 + 4 lambda alpha_x / v)^0.5) / (2 (alpha_x v t)^0.5))
      Fy_k = erf((y + Y_k/2) / (2 (alpha_y x)^0.5)) - erf((y - Y_k/2) / (2 (alpha_y x)^0.5))
      Fz = erf(Z / (2 (alpha_z x)^0.5)) - erf(-Z / (2 (alpha_z x)^0.5))  (downward dispersion
           only, BSCR14 p. 7 [71])
      v = Vs / R; the source time is t - x/v, not below 0 (BSCR p. 35 [41]; XLS MAX(0, t - x/Vc)).
    Instantaneous reaction (superposition, BSCR Appendix A.2 p. 48-49 [54-55]): no-decay model
    with source C0 + BC and source decay ks_inst, minus BC, not below 0 (XLS 'Plume Output'!C119).
    Source column (x <= 1e-4 ft): C = sum_k C_k exp(-ks t) over the layers that contain y
    (XLS 'Plume Output'!C63: IF(Fy<>0, C_k exp(-ks (t - x/v)), 0)).
    """
    t = c.time_yr if t_yr is None else t_yr
    x = np.asarray(x_ft, dtype=float)
    y = np.asarray(y_ft, dtype=float)
    x, y = np.broadcast_arrays(x, y)
    at_source = x <= X_SOURCE_FT
    xe = np.where(at_source, X_SOURCE_FT, x)
    v = c.vs_ft_yr / c.retardation                      # XLS name Vc = Vs/Ret
    tau = np.maximum(0.0, t - xe / v)
    az = c.alpha_z_ft if c.alpha_z_ft != 0 else ALPHA_Z_ZERO_FT
    erf = s.erf
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        if t > 0:
            den = 2.0 * math.sqrt(c.alpha_x_ft * v * t)                 # 'Centerline Output'!E34
            fx0 = 1.0 - erf((xe - v * t) / den)                          # gerfC, row 39
            root = math.sqrt(1.0 + 4.0 * c.lambda_per_yr * c.alpha_x_ft / v)   # erfcconst3
            fxl = np.exp(xe / 2.0 / c.alpha_x_ft * (1.0 - root)) * (1.0 - erf((xe - v * t * root) / den))
        else:
            fx0 = fxl = np.zeros_like(xe)
        fz = 2.0 * erf(c.thickness_ft / 2.0 / np.sqrt(az * xe))         # rows 40-41
        sq = 2.0 * np.sqrt(c.alpha_y_ft * xe)
        fys, geo = [], []
        for yk in src.widths_ft:
            fy = erf((y + yk / 2.0) / sq) - erf((y - yk / 2.0) / sq)     # rows 42-44
            fys.append(fy)
            geo.append(np.where(at_source, np.where(fy != 0, 1.0, 0.0), fx0 * fz * fy / 8.0))
        dec = np.exp(-src.ks_per_yr * tau)
        deci = np.exp(-src.ks_inst_per_yr * tau)
        no_decay = sum(ck * dec * g for ck, g in zip(src.increments_mg_l, geo))
        if s.compatible:
            # XLS 'Centerline Output'!C60, 'Plume Output'!D103: C_nodecay * Fx(lambda) / Fx(0), 0 if Fx(0) = 0
            ratio = np.where(fx0 == 0, 0.0, fxl / np.where(fx0 == 0, 1.0, fx0))
            first = np.where(at_source, no_decay, no_decay * ratio)
        else:
            # direct form: sum_k C_k exp(-ks tau) Fx(lambda) Fz Fy_k / 8
            direct = sum(ck * dec * fxl * fz * fy / 8.0 for ck, fy in zip(src.increments_mg_l, fys))
            first = np.where(at_source, no_decay, direct)
        before = 0.0
        if c.instantaneous_source == "mass-balance":
            # effective source = observed zone concentrations x f(tau) + BC (BC once, outermost layer)
            f = src.remaining_fraction(tau)
            nonzero = [k for k, ck in enumerate(src.increments_mg_l) if ck != 0]
            for k in nonzero:
                before = before + src.increments_mg_l[k] * f * geo[k]
            if nonzero:
                before = before + src.biodegradation_capacity_mg_l * geo[nonzero[0]]
        else:
            for k, (ck, bk, g) in enumerate(zip(src.increments_mg_l, src.bioc_mg_l, geo)):
                if ck == 0:
                    continue
                # Q1: XLS zone-2 term uses ks (not ks_inst) in compatible mode
                d = dec if (s.compatible and k == 1) else deci
                before = before + bk * d * g
        before = np.broadcast_to(np.asarray(before, dtype=float), x.shape).copy()
        inst = np.maximum(0.0, before - src.biodegradation_capacity_mg_l)
    nd = np.broadcast_to(np.asarray(no_decay, dtype=float), x.shape).copy()
    return Concentrations(nd, np.asarray(first, dtype=float), before, inst)


# ---------------------------------------------------------------------------------------------
# BIOSCREEN output screens
# ---------------------------------------------------------------------------------------------
MODELS = ("no_decay", "first_order", "instantaneous")


def grid(c: Case):
    """11 x 5 array of the output screens: x = 0 (1e-4 ft), L/10, ..., L; y = 2W/4, W/4, 0, -W/4, -2W/4
    (XLS names delx = L/10, dely = W/4; BSCR p. 38 [44])."""
    dx, dy = c.length_ft / 10.0, c.width_ft / 4.0
    x = np.array([X_SOURCE_FT] + [j * dx for j in range(1, 11)])
    y = np.array([2 * dy, dy, 0.0, -dy, -2 * dy])
    return x, y, dx, dy


@dataclass
class ModelMass:
    plume_mass_if_no_biodegradation_kg: float | None
    actual_plume_mass_kg: float | None
    mass_removed_kg: float | None
    fraction_removed: float | None
    electron_acceptor_change_kg: dict[str, float] | None
    source_mass_now_kg: float | None          # None: infinite source
    plume_volume_ft3: float | None
    source_flow_ft3_yr: float | None
    mass_flux_mg_day: list[float] | None       # at the 11 x positions; None if alpha_z != 0
    cant_calc: bool


@dataclass
class Result:
    source: Source
    x_ft: np.ndarray
    y_ft: np.ndarray
    centerline: dict[str, np.ndarray]
    arrays: dict[str, np.ndarray]              # 5 x 11, rows y = +2dy ... -2dy
    inst_before_array: np.ndarray
    grid_masses_kg: dict[str, float]           # XLS P91 (no decay), P104 (first order), P116, P121
    mass_left_source_kg: float                 # P88 (no decay and first order)
    mass_left_source_inst_kg: float            # P111
    fraction_biodegraded: dict[str, float]     # P107, P125
    models: dict[str, ModelMass]


def run(c: Case, s: Settings, volume_threshold_mg_l: float = 0.005) -> Result:
    src = source_terms(c, s)
    x, y, dx, dy = grid(c)
    cl = concentrations(c, src, s, x, np.zeros_like(x))
    X, Y = np.meshgrid(x, np.abs(y))       # rows 4-5 repeat +dy, +2dy (XLS rows 91-92 = rows 89, 88)
    ar = concentrations(c, src, s, X, Y)
    arrays = {"no_decay": ar.no_decay, "first_order": ar.first_order, "instantaneous": ar.instantaneous}
    centerline = {"no_decay": cl.no_decay, "first_order": cl.first_order, "instantaneous": cl.instantaneous}
    # Plume masses on the 11 x 5 grid, first column half width (BSCR p. 39 [45]; XLS P91, P104, P116, P121)
    cell = dx * dy * c.thickness_ft * c.porosity * s.l_per_ft3 / 1e6

    def gmass(a, ret):
        return (a.sum() - a[:, 0].sum() / 2.0) * cell * ret

    r = c.retardation
    m_nd = gmass(ar.no_decay, r)
    m_fo = gmass(ar.first_order, r if not s.compatible else 1.0)     # Q2
    m_ib = gmass(ar.inst_before, r)
    m_in = gmass(ar.instantaneous, r)
    t = c.time_yr
    # Mass that has left the source (BSCR p. 39 [45]; XLS P88, P111)
    if c.mass_kg is None:
        left, left_i = src.flux_kg_yr * t, src.flux_inst_kg_yr * t
    else:
        left = c.mass_kg - c.mass_kg * math.exp(-src.ks_per_yr * t)
        if c.instantaneous_source == "mass-balance":
            left_i = c.mass_kg * (1.0 - float(src.remaining_fraction(t)))
        else:
            left_i = c.mass_kg - c.mass_kg * math.exp(-src.ks_inst_per_yr * t)
    frac_fo = (m_nd - m_fo) / m_nd if m_nd else math.nan     # P107
    frac_in = (m_ib - m_in) / m_ib if m_ib else math.nan     # P125
    cant = m_nd < 0.5 * left or m_nd > 1.5 * left             # M17 test (also used for the inst. model)
    models = {}
    for name in MODELS:
        a = arrays[name]
        if name == "no_decay":
            m17, m19, fr = left, left, 0.0
        elif name == "first_order":
            m17, fr = left, frac_fo
            m19 = left * (1 - frac_fo)
        else:
            m17, fr = left_i, frac_in
            m19 = left_i * (1 - frac_in)
        removed = m17 - m19
        ea = None
        if name == "instantaneous" and not cant and src.biodegradation_capacity_mg_l > 0:
            ea = {k: EA_SIGN[k] * removed * src.ea_share[k] * c.utilization_factors[k] for k in UTILIZATION_FACTORS}
        ks = src.ks_inst_per_yr if name == "instantaneous" else src.ks_per_yr
        now = None if c.mass_kg is None else c.mass_kg * math.exp(-ks * t)     # M28
        if now is not None and name == "instantaneous" and c.instantaneous_source == "mass-balance":
            now = c.mass_kg * float(src.remaining_fraction(t))
        # M30: cells above 0.005 mg/L x cell volume; the cell length is D57 - C57 = dx - 1e-4 ft
        dxv = dx - X_SOURCE_FT if s.compatible else dx
        vol = None if cant else float((a > volume_threshold_mg_l).sum()) * dxv * dy * c.thickness_ft * c.porosity
        flow = None if cant else src.flow_l_yr / s.l_per_ft3                   # M31 (ft3/yr)
        # Mass flux across the cross-sections, BSCR14 p. 2 [66]; XLS 'Plume Output'!B11-L11:
        # Vs n Z dy (C(0) + 2 C(dy) + 2 C(2dy)) / 365; x = 0: source flux Q*Cavg.
        if c.alpha_z_ft != 0:
            flux = None
        else:
            f = [src.flux_kg_yr * 1e6 / DAYS_PER_YEAR]
            for j in range(1, 11):
                f.append(c.vs_ft_yr * c.porosity * c.thickness_ft * s.l_per_ft3
                         * (dy * (a[2, j] + 2 * a[1, j] + 2 * a[0, j])) / DAYS_PER_YEAR)
            flux = f
        models[name] = ModelMass(None if cant else m17, None if cant else m19, None if cant else removed,
                                 None if cant else fr, ea, now, vol, flow, flux, cant)
    return Result(src, x, y, centerline, arrays, ar.inst_before, {"no_decay": m_nd, "first_order": m_fo,
                  "inst_before": m_ib, "instantaneous": m_in}, left, left_i,
                  {"first_order": frac_fo, "instantaneous": frac_in}, models)
