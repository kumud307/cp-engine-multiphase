"""Ideal-mixture heat capacity.

For an ideal mixture the enthalpy of mixing is zero (H^E = 0), so the
mixture enthalpy is the mole-weighted sum of pure-component enthalpies,
H = Σ x_i H_i.  Differentiating with respect to T at constant P and x:

    molar basis:  Cp_mix = Σ x_i Cp_i          [J/(mol K)]
    mass  basis:  cp_mix = Σ w_i cp_i          [J/(kg K)]

Both forms are the same physics: the second follows from the first by
dividing by the mixture molar mass M = Σ x_i M_i, because w_i = x_i M_i / M.
"""
from __future__ import annotations

from collections.abc import Mapping


def cp_ideal_molar(x: Mapping[str, float], cp_molar: Mapping[str, float]) -> float:
    """Cp = Σ x_i Cp_i  (x: mole fractions, Cp_i in J/(mol K))."""
    return sum(x[k] * cp_molar[k] for k in x)


def cp_ideal_mass(w: Mapping[str, float], cp_mass: Mapping[str, float]) -> float:
    """cp = Σ w_i cp_i  (w: mass fractions, cp_i in J/(kg K))."""
    return sum(w[k] * cp_mass[k] for k in w)


def mole_to_mass_fractions(x: Mapping[str, float], MW: Mapping[str, float]) -> dict[str, float]:
    """w_i = x_i M_i / Σ_j x_j M_j."""
    M = sum(x[k] * MW[k] for k in x)
    return {k: x[k] * MW[k] / M for k in x}


def mass_to_mole_fractions(w: Mapping[str, float], MW: Mapping[str, float]) -> dict[str, float]:
    """x_i = (w_i / M_i) / Σ_j (w_j / M_j)."""
    s = sum(w[k] / MW[k] for k in w)
    return {k: (w[k] / MW[k]) / s for k in w}
