"""Relation between Cp and Cv.

From the fundamental property relations (Maxwell relations) one obtains,
exactly, for any homogeneous substance or mixture at fixed composition:

    Cp − Cv = T v α² / κ_T

    v    specific volume (m^3/kg) → gives Cp, Cv per kg
    α    = (1/v)(∂v/∂T)_P   volumetric thermal expansion coefficient  [1/K]
    κ_T  = −(1/v)(∂v/∂P)_T  isothermal compressibility                [1/Pa]

For an ideal gas α = 1/T and κ_T = 1/P, which gives Cp − Cv = R (per mole).

Mixture α and κ_T (assumption)
------------------------------
For a mixture we assume ideal volume of mixing (V^E = 0).  Then the
mixture volume is the sum of component volumes, and its α and κ_T are the
*volume-fraction* averages:

    α_mix   = Σ φ_i α_i ,     κ_T,mix = Σ φ_i κ_T,i ,     v_mix = Σ w_i / ρ_i
"""
from __future__ import annotations

from collections.abc import Mapping


def cp_minus_cv(T: float, v: float, alpha: float, kappa_T: float) -> float:
    """Cp − Cv = T v α² / κ_T  (same basis as v: per kg if v in m^3/kg)."""
    return T * v * alpha ** 2 / kappa_T


def cp_from_cv(cv: float, T: float, v: float, alpha: float, kappa_T: float) -> float:
    """Cp = Cv + T v α² / κ_T."""
    return cv + cp_minus_cv(T, v, alpha, kappa_T)


def cv_from_cp(cp: float, T: float, v: float, alpha: float, kappa_T: float) -> float:
    """Cv = Cp − T v α² / κ_T."""
    return cp - cp_minus_cv(T, v, alpha, kappa_T)


def mixture_volumetric_properties(w: Mapping[str, float], rho: Mapping[str, float],
                                  alpha: Mapping[str, float], kappa: Mapping[str, float]
                                  ) -> tuple[float, float, float]:
    """Return (rho_mix, alpha_mix, kappa_mix) assuming ideal volume of mixing.

    w: mass fractions; rho, alpha, kappa: pure-component values at T, P.
    """
    v = sum(w[k] / rho[k] for k in w)              # m^3 per kg of mixture
    phi = {k: (w[k] / rho[k]) / v for k in w}      # volume fractions
    a = sum(phi[k] * alpha[k] for k in w)
    kt = sum(phi[k] * kappa[k] for k in w)
    return 1.0 / v, a, kt
