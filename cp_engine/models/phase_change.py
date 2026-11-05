"""Phase change: melting fractions and solid–liquid equilibrium (solubility).

When a component melts, dissolves or crystallises while T changes, the
enthalpy changes not only through Cp but also through the *amount* that
changes phase.  The measured ("apparent") heat capacity is then

    cp_app = dH/dT = Σ n_i Cp_i  +  Cp^E  +  Σ_i Δh_i (dn_i/dT)
                     (sensible)    (excess)  (latent / apparent part)

Two ways of describing how much has changed phase are provided here:

1. Smooth melting of a pure substance (PCM slurry, ice in pure water).
   A pure substance melts at a single T_m, which would give an infinite cp
   spike.  Real samples melt over a finite range (impurities, size
   distribution, instrument heating rate), so the melt fraction is modelled
   with a smooth step:
        f(T) = ½ [1 + erf((T − T_m) / (σ √2))],   σ = width / 4
   (≈ 95 % of the melting happens within ±width/2 of T_m).  df/dT is a
   Gaussian whose area is exactly 1. The latent part is f'(T)[H_L(T) − H_S(T)] with
   H_L − H_S = ΔH_fus + ∫_{Tm}^{T} ΔCp dT (Kirchhoff), so its area is n[ΔH_fus + ½σ² dΔCp/dT] for a
   linear ΔCp: ΔH_fus exactly only if ΔCp is constant (the correction is < 0.001 % for the PCM data).
   The step is cut at ±8σ (f = 0 or 1 exactly; neglected tail mass 1.2e-15), so the property
   correlations of a phase that is absent are not evaluated far outside their ranges.

2. Solid–liquid equilibrium (SLE) in a solution (ice slurry in a brine or
   an alcohol solution, salt crystals in their saturated solution).
   The solubility of solid k in the liquid follows from equal chemical
   potentials of pure solid and dissolved k (ideal solubility equation):
        ln(x_k γ_k) = −(ΔH_fus/R)(1/T − 1/T_m)
   With γ_k = 1 this is the ideal solubility / freezing-point-depression
   equation.  For electrolytes it is useless (NaCl would be "insoluble"),
   so a tabulated solubility (g solute/100 g solvent) can be used instead.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.interpolate import PchipInterpolator

from ..database import R


# ---------------------------------------------------------------------------
# 1. Smooth melting of a pure substance
# ---------------------------------------------------------------------------
TAIL_SIGMAS = 8.0      # Gaussian cut-off: 1 − erf(8/√2) = 1.2e-15


def melt_fraction(T: float, Tm: float, width: float) -> float:
    """Liquid fraction f(T) = ½[1 + erf((T − Tm)/(σ√2))], σ = width/4."""
    sigma = width / 4.0
    z = (T - Tm) / sigma
    if z <= -TAIL_SIGMAS:
        return 0.0
    if z >= TAIL_SIGMAS:
        return 1.0
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def melt_fraction_derivative(T: float, Tm: float, width: float) -> float:
    """df/dT = Gaussian with standard deviation σ = width/4 (area 1)."""
    sigma = width / 4.0
    if abs(T - Tm) >= TAIL_SIGMAS * sigma:
        return 0.0
    return math.exp(-0.5 * ((T - Tm) / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi))


# ---------------------------------------------------------------------------
# 2. Solid–liquid equilibrium
# ---------------------------------------------------------------------------
def sle_rhs(T: float, Tm: float, dHfus: float, dCp: float = 0.0) -> float:
    """Right-hand side of the solid–liquid equilibrium equation, ln(x·γ) at temperature T:

        ln(x γ) = −(ΔH_fus/R)(1/T − 1/T_m) + (ΔCp/R)(T_m/T − 1 − ln(T_m/T))

    ΔCp = Cp,liquid − Cp,solid at T_m. The second term accounts for the temperature dependence of
    the enthalpy of fusion (Kirchhoff's law); with ΔCp = 0 it is the textbook ideal equation.
    """
    return -dHfus / R * (1.0 / T - 1.0 / Tm) + dCp / R * (Tm / T - 1.0 - math.log(Tm / T))


def ideal_solubility(T: float, Tm: float, dHfus: float, dCp: float = 0.0) -> float:
    """Ideal (γ = 1) mole-fraction solubility x_sat = exp[sle_rhs(T)], capped at 1.

    Above T_m the 'solid' cannot exist, so x_sat = 1 (fully liquid / miscible).
    """
    if T >= Tm:
        return 1.0
    return math.exp(sle_rhs(T, Tm, dHfus, dCp))


def ideal_freezing_point(x: float, Tm: float, dHfus: float, dCp: float = 0.0) -> float:
    """Inverse of the ideal solubility equation: liquidus temperature for mole fraction x.

    With ΔCp = 0 this is explicit, 1/T = 1/Tm − (R/ΔH_fus) ln x; otherwise it is solved numerically.
    """
    if x >= 1.0:
        return Tm
    if dCp == 0.0:
        return 1.0 / (1.0 / Tm - R / dHfus * math.log(x))
    return _root_below(lambda T: sle_rhs(T, Tm, dHfus, dCp) - math.log(x), Tm)


def _root_below(f, Tm: float) -> float:
    """First root of f(T) = 0 below T_m, found by stepping down from T_m and refining with brentq.

    (With the ΔCp term the equation is not monotonic far below T_m, so a wide fixed bracket can fail;
    the physically relevant liquidus is the first crossing below T_m.)
    """
    from scipy.optimize import brentq
    hi = Tm * (1 - 1e-12)
    f_hi = f(hi)
    T = hi
    while T > 0.3 * Tm:
        lo = T - 0.005 * Tm
        if f_hi * f(lo) <= 0:
            return brentq(f, lo, T)
        T = lo
    raise ValueError("no liquidus temperature found above 0.3 T_m")


def liquidus_temperature(x: float, Tm: float, dHfus: float, ln_gamma=None, dCp: float = 0.0) -> float:
    """Freezing (liquidus) temperature of a solution with solvent mole fraction x.

    Solves ln(x) + ln γ(T) = sle_rhs(T) for T; `ln_gamma(T)` is the activity coefficient of
    the crystallising component at fixed liquid composition (None → ideal solution).
    """
    if x >= 1.0:
        return Tm
    f = lambda T: math.log(x) + (ln_gamma(T) if ln_gamma else 0.0) - sle_rhs(T, Tm, dHfus, dCp)
    return _root_below(f, Tm)


class TableSolubility:
    """Tabulated solubility s(T) in g solute per 100 g solvent (monotone PCHIP interpolation)."""

    def __init__(self, entry: dict):
        self.T = np.asarray(entry["T"], dtype=float)
        self.s = np.asarray(entry["values"], dtype=float)
        self.entry = entry
        self._f = PchipInterpolator(self.T, self.s, extrapolate=True)
        dh = entry.get("dH_sol")
        self.dH_sol = float(dh["value"]) * (1000.0 if "kJ" in dh.get("unit", "kJ/mol") else 1.0) if dh else None

    @property
    def Tmin(self) -> float:
        return float(self.T.min())

    @property
    def Tmax(self) -> float:
        return float(self.T.max())

    def mass_ratio(self, T: float) -> float:
        """kg solute that can dissolve per kg solvent at T."""
        return float(self._f(T)) / 100.0
