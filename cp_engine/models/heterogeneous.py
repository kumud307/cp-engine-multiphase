"""Heterogeneous mixtures: suspensions, nanofluids, emulsions.

In a heterogeneous mixture each phase keeps its own properties, and
enthalpy is extensive, so the total enthalpy is the sum over phases:

    H = Σ_phases m_p h_p      ⇒      cp_mix = Σ_phases w_p cp_p       (mass basis)

where w_p is the *mass* fraction of phase p.  (Interfacial enthalpy is
neglected – an assumption that can matter for very small nanoparticles.)

Nanofluid models (particles s in base fluid f, particle volume fraction φ)
---------------------------------------------------------------------------
Model I  (volume-weighted cp; Pak & Cho 1998):
    cp = φ cp_s + (1 − φ) cp_f
    Thermodynamically INCONSISTENT: cp is energy per unit *mass*, so a volume
    fraction cannot weight it. The model ignores the density difference.

Model II (thermal equilibrium; Xuan & Roetzel 2000):
    (ρ cp)_mix = φ (ρ cp)_s + (1 − φ)(ρ cp)_f,   ρ_mix = φ ρ_s + (1 − φ) ρ_f
    ⇒ cp = [φ ρ_s cp_s + (1 − φ) ρ_f cp_f] / [φ ρ_s + (1 − φ) ρ_f]
    This is exactly the mass-weighted rule because w_s = φ ρ_s / ρ_mix.
"""
from __future__ import annotations

from collections.abc import Mapping

from ..database import Database


def cp_mass_weighted(w_phase: Mapping[str, float], cp_phase: Mapping[str, float]) -> float:
    """cp = Σ w_p cp_p over phases (J/(kg K))."""
    return sum(w_phase[p] * cp_phase[p] for p in w_phase)


def volume_to_mass_fraction(phi: float, rho_s: float, rho_f: float) -> float:
    """Particle mass fraction w_s = φ ρ_s / [φ ρ_s + (1 − φ) ρ_f]."""
    return phi * rho_s / (phi * rho_s + (1.0 - phi) * rho_f)


def mass_to_volume_fraction(w: float, rho_s: float, rho_f: float) -> float:
    """Particle volume fraction φ = (w/ρ_s) / [w/ρ_s + (1 − w)/ρ_f]."""
    return (w / rho_s) / (w / rho_s + (1.0 - w) / rho_f)


def cp_model_I(phi: float, cp_s: float, cp_f: float) -> float:
    """Model I (volume-weighted, inconsistent): cp = φ cp_s + (1 − φ) cp_f."""
    return phi * cp_s + (1.0 - phi) * cp_f


def cp_model_II(phi: float, rho_s: float, cp_s: float, rho_f: float, cp_f: float) -> float:
    """Model II (thermal equilibrium): cp = [φρ_s cp_s + (1−φ)ρ_f cp_f]/[φρ_s + (1−φ)ρ_f]."""
    return (phi * rho_s * cp_s + (1.0 - phi) * rho_f * cp_f) / (phi * rho_s + (1.0 - phi) * rho_f)


def nanofluid_cp(db: Database, particle: str, fluid: str, phi: float, T: float,
                 model: str = "II") -> float:
    """Specific heat (J/(kg K)) of a single-component base fluid + one particle type."""
    s, f = db.get(particle), db.get(fluid)
    cp_s, cp_f = s.cp_mass("solid", T), f.cp_mass("liquid", T)
    if model.upper() == "I":
        return cp_model_I(phi, cp_s, cp_f)
    rho_s, rho_f = s.density("solid", T), f.density("liquid", T)
    return cp_model_II(phi, rho_s, cp_s, rho_f, cp_f)
