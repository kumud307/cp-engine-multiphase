"""Validation of normalised input.

Checks performed here (before any thermodynamics is attempted):

* temperatures > 0 K, pressure > 0, amounts ≥ 0 (fractions are checked in
  the normaliser while they are being read)
* T1, T2, dT mutually consistent when all three are given
* declared phases physically possible at the given T and P
  (e.g. water declared solid at 350 K, 1 atm → PhaseError)
* process description consistent (adiabatic with Q ≠ 0 → error)
* positive Cv, γ ≥ 1, positive κ_T and density

Correlation temperature ranges are checked when the correlations are
evaluated (database.check_range), and consistency between over-specified
energy inputs (e.g. Q vs ΔH at constant P) is checked by the solver,
because it requires applying the thermodynamic relations.

All problems raise subclasses of CpEngineError with a specific message.
"""
from __future__ import annotations

import math

from .database import Database
from .errors import InputError, PhaseError
from .inputs import NormalisedInput

PHASE_TOL_K = 0.5   # tolerance for phase-consistency checks


def validate(inp: NormalisedInput, db: Database) -> None:
    k, log = inp.known, inp.log

    # ---- every number must be finite (NaN passes "x <= 0" tests silently) -----
    for sym, v in k.items():
        if isinstance(v, (int, float)) and not isinstance(v, bool) and not math.isfinite(v):
            raise InputError(f"{sym} = {v} is not a finite number.")
    for comp, n in (inp.amounts or {}).items():
        if not math.isfinite(n):
            raise InputError(f"The amount of {comp} is not a finite number ({n}).")

    # ---- state ----------------------------------------------------------------
    for sym in ("T", "T1", "T2"):
        if sym in k and k[sym] <= 0:
            raise InputError(f"{sym} = {k[sym]:.4g} K is not above absolute zero.")
    if "P" in k and k["P"] <= 0:
        raise InputError(f"P = {k['P']:.4g} Pa must be positive.")
    if all(s in k for s in ("T1", "T2", "dT")) and abs(k["T2"] - k["T1"] - k["dT"]) > 1e-6:
        raise InputError(f"Inconsistent temperatures: T2 − T1 = {k['T2'] - k['T1']:.4g} K but dT = {k['dT']:.4g} K.")
    if "dT" in k and k["dT"] == 0:
        raise InputError("dT = 0: a heat capacity cannot be obtained from a zero temperature change.")
    if "T1" in k and "T2" in k and k["T1"] == k["T2"]:
        raise InputError("T1 = T2: a heat capacity cannot be obtained from a zero temperature change.")

    # ---- amounts --------------------------------------------------------------
    for sym in ("m", "n"):
        if sym in k and k[sym] <= 0:
            raise InputError(f"Total {'mass' if sym == 'm' else 'amount'} must be positive.")

    # ---- process --------------------------------------------------------------
    if "adiabatic" in inp.flags and "Q" in k and abs(k["Q"]) > 0:
        raise InputError(f"The process is adiabatic but Q = {k['Q']:.6g} J was given; adiabatic means Q = 0.")
    if "W_pv" in k and "constant_V" in inp.flags and abs(k["W_pv"]) > 0:
        raise InputError("Expansion work W_pv must be zero at constant volume.")

    # ---- properties -----------------------------------------------------------
    for sym in ("Cv_mass", "Cv_molar", "Cv_total", "kappa_T", "rho"):
        if sym in k and k[sym] <= 0:
            raise InputError(f"{sym} must be positive.")
    if "gamma" in k and k["gamma"] < 1.0:
        raise InputError(f"gamma = Cp/Cv = {k['gamma']:.4g} < 1 is impossible (Cp ≥ Cv for stable matter).")

    # ---- phase consistency of declared phases ---------------------------------
    if inp.amounts:
        temps = [k[s] for s in ("T", "T1", "T2") if s in k]
        if temps:
            check_phases(inp, db, min(temps), max(temps))


def check_phases(inp: NormalisedInput, db: Database, T_lo: float, T_hi: float) -> None:
    """Compare declared phases with pure-component transition temperatures."""
    P = inp.known.get("P", 101325.0)
    log = inp.log
    for key, phase in inp.phase_spec.items():
        c = db.get(key)
        Tm, Tb = c.Tm, c.Tb_at(P)
        if phase == "solid" and Tm is not None and T_hi > Tm + PHASE_TOL_K:
            raise PhaseError(f"{c.name} was declared solid, but at {T_hi:.2f} K it is above its melting point "
                             f"({Tm:.2f} K). Declare it liquid, or use phase 'auto' to model melting.")
        if phase == "liquid":
            if Tm is not None and T_lo < Tm - PHASE_TOL_K:
                log.warn(f"{c.name} was declared liquid at {T_lo:.2f} K, below its pure melting point "
                         f"({Tm:.2f} K): treated as sub-cooled liquid (valid only if dissolved / freezing-point "
                         "depressed). Use phase 'auto' to let the engine model freezing.")
            if Tb is not None and T_hi > Tb + PHASE_TOL_K:
                log.warn(f"{c.name} was declared liquid at {T_hi:.2f} K, above its boiling point at this "
                         f"pressure ({Tb:.2f} K): treated as superheated / dissolved liquid.")
        if phase == "gas" and Tb is not None and T_lo < Tb - PHASE_TOL_K:
            log.warn(f"{c.name} was declared gas at {T_lo:.2f} K, below its boiling point ({Tb:.2f} K): "
                     "valid only as a vapour below its saturation partial pressure (ideal gas assumed).")
