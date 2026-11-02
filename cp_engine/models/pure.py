"""Pure-component heat capacities and enthalpies.

Heat capacities come straight from the database correlations
(Component.cp).  This module adds *enthalpies*, which are needed whenever a
phase change occurs (latent heat) or when a mean Cp over an interval is
computed from ΔH.

Reference state (H = 0) for each component
------------------------------------------
* components with liquid data and a melting point: liquid at T_m
* otherwise: the stable phase that has data (e.g. a solid oxide) at 298.15 K

Then, by the path rule for a state function (Kirchhoff's law),

    H_L(T) = ∫_{Tm}^{T} Cp_L dT
    H_S(T) = −ΔH_fus + ∫_{Tm}^{T} Cp_S dT
    H_G(T) = H_L(T_b) + ΔH_vap + ∫_{Tb}^{T} Cp_G dT

so the enthalpy of melting at any T is H_L(T) − H_S(T) = ΔH_fus + ∫(Cp_L − Cp_S)dT.
"""
from __future__ import annotations

from ..database import Component, check_range
from ..errors import DataError
from ..report import Log, null_log

T_REF_DEFAULT = 298.15


def reference_temperature(comp: Component) -> tuple[str, float]:
    """Phase and temperature at which H = 0 for this component."""
    if comp.has_phase("liquid") and comp.Tm is not None:
        return "liquid", comp.Tm
    for ph in ("liquid", "solid", "gas"):
        if comp.has_phase(ph):
            return ph, T_REF_DEFAULT
    raise DataError(f"{comp.name} has no Cp correlation at all.")


def enthalpy(comp: Component, phase: str, T: float, log: Log | None = None) -> float:
    """Molar enthalpy H_i^phase(T) in J/mol, relative to the component reference state."""
    log = log or null_log()
    if not comp.has_phase(phase):
        raise DataError(f"No {phase}-phase Cp correlation for {comp.name}; cannot compute its enthalpy.")
    corr = comp.cp_corr[phase]
    check_range(T, corr.Tmin, corr.Tmax, f"Cp of {comp.name} ({phase})", log)
    ref_phase, T0 = reference_temperature(comp)

    if phase == ref_phase:
        return corr.integral(T0, T)

    if phase == "solid":
        if ref_phase == "liquid":
            if comp.dHfus is None:
                raise DataError(f"ΔH_fus of {comp.name} is needed for solid↔liquid enthalpies.")
            return -comp.dHfus + corr.integral(T0, T)          # T0 = Tm
        return corr.integral(T0, T)

    if phase == "liquid":           # reference is solid (no Tm-based liquid reference)
        if comp.Tm is None or comp.dHfus is None:
            raise DataError(f"Tm and ΔH_fus of {comp.name} are needed for its liquid enthalpy.")
        H_s_Tm = comp.cp_corr[ref_phase].integral(T0, comp.Tm)
        return H_s_Tm + comp.dHfus + corr.integral(comp.Tm, T)

    if phase == "gas":
        if comp.Tb is None or comp.dHvap is None:
            raise DataError(f"Tb and ΔH_vap of {comp.name} are needed for its vapour enthalpy.")
        H_L_Tb = enthalpy(comp, "liquid", comp.Tb, log) if comp.has_phase("liquid") else 0.0
        return H_L_Tb + comp.dHvap + corr.integral(comp.Tb, T)

    raise DataError(f"Unknown phase '{phase}'.")
