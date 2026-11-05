"""System model: classification of the mixture and Route A (property route).

Route A predicts Cp from state (T, P) and composition only.  The steps are

1.  Decide how each component behaves in the temperature window of interest
    (fixed solid / liquid / gas, dissolved solute, melting, crystallising
    from solution, boiling).  User-declared phases are respected; 'auto'
    components are inferred from T_m, T_b(P), miscibility and solubility.
2.  Group the liquids into liquid phases (miscible liquids → one solution;
    immiscible liquids → separate phases = emulsion).  Every solid
    component is a separate particulate phase (suspension / slurry).
3.  At any temperature T, 'split' the components between phases:
    smooth melting f(T), solid–liquid equilibrium, tabulated solubility.
4.  Total enthalpy of the system (state function, extensive):
        H(T) = Σ_phases Σ_i n_i^p H_i^p(T)  +  Σ_liquid phases n^p H^E,p(T, x^p)
5.  Apparent heat capacity = dH/dT, split into three parts by the chain rule:
        dH/dT = Σ n_i^p Cp_i^p             (ideal / sensible part)
              + Σ n^p Cp^E,p              (excess part, real solutions, NRTL)
              + Σ (∂H/∂n_i^p)(dn_i^p/dT)   (latent / apparent part, phase change)
    The latent part is evaluated as [H(T; n(T+h)) − H(T; n(T−h))]/(2h), i.e.
    the change of H caused *only* by the change of phase amounts.
6.  Mean Cp over [T1, T2] = [H(T2) − H(T1)] / (m ΔT), which includes all
    latent heat released/absorbed in the interval.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy import integrate

from .database import P_ATM, Component, Database
from .errors import DataError, InputError, PhaseError
from .models import cv_cp, phase_change
from .models.heterogeneous import cp_model_I
from .models.nrtl import NRTL
from .models.pure import enthalpy as pure_enthalpy
from .report import Log

DEFAULT_MELT_WIDTH = 2.0   # K; width of the smoothed melting step (≈95 % of melting)

BEHAVIOUR_TEXT = {
    "fixed_solid": "solid (no phase change)",
    "fixed_liquid": "liquid",
    "fixed_gas": "gas (ideal)",
    "dissolved": "dissolved solute",
    "melt_smooth": "melting / freezing (pure, smoothed over a finite range)",
    "sle_ideal": "crystallises from / dissolves into the solution (solid–liquid equilibrium)",
    "sle_table": "partly dissolved (tabulated solubility, saturated slurry)",
    "boil_smooth": "boiling (pure, smoothed over a finite range)",
}
PHASE_CHANGE_BEHAVIOURS = {"melt_smooth", "sle_ideal", "sle_table", "boil_smooth"}


def pure_phase_at(comp: Component, T: float, P: float = P_ATM) -> str:
    """Stable phase of the *pure* component at (T, P) from T_m and T_b(P)."""
    Tb = comp.Tb_at(P)
    if comp.Tm is not None and T < comp.Tm and comp.has_phase("solid"):
        return "solid"
    if Tb is not None and T > Tb and comp.has_phase("gas"):
        return "gas"
    if comp.has_phase("liquid"):
        return "liquid"
    for ph in ("solid", "gas"):
        if comp.has_phase(ph):
            return ph
    raise DataError(f"{comp.name} has no Cp data for any phase.")


@dataclass
class Plan:
    """How one component behaves in the temperature window."""
    key: str
    comp: Component
    n: float                      # mol (in the calculation basis)
    behaviour: str
    group: str | None = None      # liquid phase it belongs to
    solvent: str | None = None    # for sle_table
    reason: str = ""


@dataclass
class CpBreakdown:
    """Heat capacity of the whole basis amount, split into contributions (J/K)."""
    ideal: float
    excess: float
    latent: float
    mass: float                   # kg in basis
    moles: float                  # mol in basis
    kind: str                     # "point" or "mean"
    T: float | None = None
    T1: float | None = None
    T2: float | None = None
    phase_amounts: dict[str, dict[str, float]] = field(default_factory=dict)
    extras: dict[str, Any] = field(default_factory=dict)

    @property
    def total(self) -> float:
        return self.ideal + self.excess + self.latent

    def per_mass(self) -> dict[str, float]:
        return {k: getattr(self, k) / self.mass for k in ("ideal", "excess", "latent")} | {"total": self.total / self.mass}

    def per_mole(self) -> dict[str, float]:
        return {k: getattr(self, k) / self.moles for k in ("ideal", "excess", "latent")} | {"total": self.total / self.moles}


class SystemModel:
    """Mixture of components with given amounts; computes Route-A heat capacities."""

    def __init__(self, db: Database, amounts: dict[str, float], phase_spec: dict[str, str],
                 P: float = P_ATM, options: dict[str, Any] | None = None, log: Log | None = None):
        if not amounts:
            raise InputError("No components given.")
        self.db, self.P = db, P
        self.amounts = {k: float(v) for k, v in amounts.items() if v > 0}
        if not self.amounts:
            raise InputError("All component amounts are zero.")
        self.comps = {k: db.get(k) for k in self.amounts}
        self.phase_spec = {k: phase_spec.get(k, "auto") for k in self.amounts}
        self.options = dict(options or {})
        self.width = float(self.options.get("melt_width", DEFAULT_MELT_WIDTH))
        if self.width <= 0:
            raise InputError("melt_width must be positive (K).")
        self.solution_model = self.options.get("solution_model", "auto")
        # "auto": use NRTL activity coefficients in the solid–liquid equilibrium only if the parameter
        # set was fitted to freezing-point data (then ONE liquid model gives both H^E and the SLE);
        # otherwise ideal (γ = 1). "ideal" / "nrtl" force the choice.
        self.sle_activity = self.options.get("sle_activity", "auto")
        # which NRTL parameter set to use for pairs that have several: "HE" (excess-enthalpy
        # fitted – default, appropriate for caloric properties) or "VLE"
        self.nrtl_parameters = self.options.get("nrtl_parameters", "HE")
        self.log = log or Log()
        if abs(P - P_ATM) > 1.0:
            self.log.assume(f"P = {P / 1e5:.4g} bar: pure-component Cp and H corrected from 1 atm with the exact "
                            "relations (∂Cp/∂P)_T = −T(∂²v/∂T²)_P and (∂H/∂P)_T = v − T(∂v/∂T)_P, v(T,P) from the "
                            "Tait equation (validated against NIST to 1.2 % up to 100 MPa).")
            if abs(P - P_ATM) > 1.0e8:
                self.log.warn(f"P = {P / 1e6:.0f} MPa is beyond the 100 MPa range over which the pressure "
                              "correction was validated against NIST isotherms.")
        self.mass = sum(n * self.comps[k].MW for k, n in self.amounts.items())
        self.moles = sum(self.amounts.values())
        self.plans: dict[str, Plan] = {}
        self.groups: dict[str, list[str]] = {}
        self._nrtl_cache: dict[tuple[str, ...], NRTL | None] = {}
        self._H_cache: dict[tuple, float] = {}
        self._stability_checked: set = set()
        self.window: tuple[float, float] | None = None

    # ==================================================================
    # 1. Behaviour of each component in the window [T_lo, T_hi]
    # ==================================================================
    def configure(self, T_lo: float, T_hi: float) -> None:
        if T_lo > T_hi:
            T_lo, T_hi = T_hi, T_lo
        self.window = (T_lo, T_hi)
        margin = 0.75 * self.width                    # 3σ of the smooth step
        lo, hi = T_lo - margin, T_hi + margin
        keys = list(self.amounts)

        def liquid_capable(k: str) -> bool:
            c, s = self.comps[k], self.phase_spec[k]
            if s == "liquid":
                return True
            if s == "auto" and c.has_phase("liquid"):
                Tb = c.Tb_at(self.P)
                return not (Tb is not None and lo > Tb and c.has_phase("gas"))
            return False

        def soluble_in(k: str) -> str | None:
            """Liquid solvent present in which solid k has tabulated solubility."""
            for j in keys:
                if j != k and liquid_capable(j) and self.db.solubility_entry(k, j):
                    return j
            return None

        plans: dict[str, Plan] = {}
        for k in keys:
            c, spec, n = self.comps[k], self.phase_spec[k], self.amounts[k]
            if spec in ("solid", "liquid", "gas"):
                if not c.has_phase(spec):
                    raise DataError(f"{c.name} was declared {spec}, but the database has no {spec}-phase "
                                    f"Cp correlation for it (available: {', '.join(c.cp_corr)}).")
                plans[k] = Plan(k, c, n, f"fixed_{spec}", reason="declared by user")
                continue
            if spec == "dissolved":
                plans[k] = Plan(k, c, n, "dissolved", reason="declared by user")
                continue
            if spec != "auto":
                raise InputError(f"Unknown phase '{spec}' for {c.name}. Use solid, liquid, gas, dissolved or auto.")

            # ---- automatic inference --------------------------------------
            if not c.has_phase("liquid"):
                solvent = soluble_in(k)
                if solvent and c.has_phase("solid"):
                    plans[k] = Plan(k, c, n, "sle_table", solvent=solvent,
                                    reason=f"soluble in {solvent} (tabulated solubility)")
                elif c.has_phase("solid"):
                    if c.Tm is not None and T_hi > c.Tm:
                        raise PhaseError(f"{c.name} melts at {c.Tm:.1f} K but has no liquid-phase data; "
                                         f"cannot handle T up to {T_hi:.1f} K.")
                    plans[k] = Plan(k, c, n, "fixed_solid", reason="no liquid data; insoluble solid")
                else:
                    plans[k] = Plan(k, c, n, "fixed_gas", reason="only gas-phase data")
                continue

            Tm, Tb = c.Tm, c.Tb_at(self.P)
            partners = [j for j in keys if j != k and (
                (liquid_capable(j) and self._miscible(k, j))
                or self.phase_spec[j] == "dissolved"
                or (self.phase_spec[j] == "auto" and self.db.solubility_entry(j, k) is not None))]
            if Tb is not None and lo > Tb and c.has_phase("gas"):
                plans[k] = Plan(k, c, n, "fixed_gas", reason=f"T above boiling point {Tb:.1f} K")
            elif Tm is None or lo > Tm:
                if Tb is not None and lo <= Tb <= hi:
                    if not partners and c.has_phase("gas"):
                        plans[k] = Plan(k, c, n, "boil_smooth", reason=f"boiling point {Tb:.1f} K inside the range")
                    else:
                        self.log.warn(f"The normal boiling point of {c.name} ({Tb:.1f} K) lies in the temperature range, "
                                      "but vapour–liquid equilibrium of mixtures is not modelled; it is treated as liquid.")
                        plans[k] = Plan(k, c, n, "fixed_liquid", reason="liquid (boiling not modelled in a mixture)")
                else:
                    plans[k] = Plan(k, c, n, "fixed_liquid", reason=f"above melting point {Tm:.1f} K" if Tm else "liquid")
            elif partners:
                if not c.has_phase("solid") or c.dHfus is None:
                    self.log.warn(f"{c.name} might crystallise below {Tm:.1f} K but has no solid data; treated as liquid.")
                    plans[k] = Plan(k, c, n, "fixed_liquid", reason="no solid data")
                else:
                    plans[k] = Plan(k, c, n, "sle_ideal",
                                    reason=f"may freeze out of solution with {', '.join(partners)} (T_m = {Tm:.1f} K)")
            elif hi < Tm:
                if c.has_phase("solid"):
                    plans[k] = Plan(k, c, n, "fixed_solid", reason=f"below melting point {Tm:.1f} K")
                else:
                    self.log.warn(f"{c.name} is below its melting point ({Tm:.1f} K) but has no solid data; "
                                  "treated as sub-cooled liquid.")
                    plans[k] = Plan(k, c, n, "fixed_liquid", reason="sub-cooled liquid (no solid data)")
            else:
                if c.has_phase("solid") and c.dHfus is not None:
                    plans[k] = Plan(k, c, n, "melt_smooth", reason=f"melting point {Tm:.1f} K inside the range")
                else:
                    self.log.warn(f"{c.name} melts in the range but lacks solid data; treated as liquid.")
                    plans[k] = Plan(k, c, n, "fixed_liquid", reason="no solid data")
        self.plans = plans
        self._group_liquids()
        self._nrtl_cache.clear()
        self._H_cache.clear()

    # ------------------------------------------------------------------
    def _pair_state(self, a: str, b: str) -> tuple[float | None, float]:
        """(temperature, mole fraction of a in the a–b pair) used to decide miscibility."""
        T = 0.5 * (self.window[0] + self.window[1]) if self.window else None
        na, nb = self.amounts.get(a, 0.0), self.amounts.get(b, 0.0)
        return T, (na / (na + nb) if na + nb > 0 else 0.5)

    def _miscible(self, a: str, b: str) -> bool:
        T, x = self._pair_state(a, b)
        ok, why = self.db.miscibility(a, b, T, x)
        if why.startswith("database binodal") and self.window and self.window[1] - self.window[0] > 0:
            lo = self.db.miscibility(a, b, self.window[0], x)[0]
            hi = self.db.miscibility(a, b, self.window[1], x)[0]
            if lo != hi:
                n1, n2 = sorted((self.comps[a].name, self.comps[b].name))
                self.log.warn(f"{n1} + {n2} crosses its liquid–liquid coexistence "
                              f"curve inside [{self.window[0]:.1f}, {self.window[1]:.1f}] K; the phase split is "
                              f"evaluated at the mid-point ({why}). Use a narrower interval on either side.")
        if why.startswith("database binodal"):
            n1, n2 = sorted((self.comps[a].name, self.comps[b].name))
            self.log.note(f"Miscibility of {n1} + {n2}: {why}.")
        return ok

    def _group_liquids(self) -> None:
        """Assign liquid-capable components to liquid phases."""
        solvents = [k for k, p in self.plans.items()
                    if p.behaviour in ("fixed_liquid", "sle_ideal", "melt_smooth", "boil_smooth")]
        solutes = [k for k, p in self.plans.items() if p.behaviour in ("dissolved", "sle_table")]
        groups: dict[str, list[str]] = {}
        user = self.options.get("liquid_phases")
        if user:
            for gname, members in user.items():
                keys = [self.db.get(m).key for m in members]
                groups[gname] = [k for k in keys if k in solvents]
            missing = [k for k in solvents if not any(k in g for g in groups.values())]
            if missing:
                raise InputError(f"liquid_phases does not assign: {', '.join(missing)}.")
        else:
            raw: list[list[str]] = []
            for k in sorted(solvents, key=lambda k: -self.amounts[k]):
                placed = None
                fits = [g for g in raw if all(self._miscible(k, j) for j in g)]
                if fits:
                    placed = fits[0]
                    placed.append(k)
                    if len(fits) > 1:
                        self.log.warn(f"{self.comps[k].name} is miscible with more than one liquid phase; "
                                      f"it was placed entirely in the phase rich in {self.comps[fits[0][0]].name} "
                                      "(partitioning between phases is not modelled).")
                else:
                    raw.append([k])
            for g in raw:
                for i, a in enumerate(g):
                    for b in g[i + 1:]:
                        ok, why = self.db.miscibility(a, b, *self._pair_state(a, b))
                        if why.startswith("heuristic"):
                            self.log.assume(f"Miscibility of {a}–{b}: {why}.")
            for g in raw:
                name = "liquid" if len(raw) == 1 else f"liquid ({self.comps[g[0]].name}-rich)"
                groups[name] = g
        for k in solutes:
            p = self.plans[k]
            target = None
            for gname, members in groups.items():
                if (p.solvent and p.solvent in members) or any(self.db.solubility_entry(k, m) for m in members):
                    target = gname
                    break
            if target is None and groups:
                target = max(groups, key=lambda g: sum(self.amounts[m] for m in groups[g]))
            if target is None:
                raise InputError(f"{self.comps[k].name} is a dissolved solute but there is no liquid solvent.")
            groups[target].append(k)
        self.groups = groups
        for gname, members in groups.items():
            for k in members:
                self.plans[k].group = gname

    # ==================================================================
    # 2. Phase split at temperature T
    # ==================================================================
    def _nrtl(self, members: tuple[str, ...]) -> NRTL | None:
        if self.solution_model == "ideal" or len(members) < 2:
            return None
        if members not in self._nrtl_cache:
            m = NRTL(list(members), self.db, log=self.log, prefer=self.nrtl_parameters)
            self._nrtl_cache[members] = m if m.has_any_parameters else None
            if self.solution_model == "nrtl" and not m.has_any_parameters:
                self.log.warn(f"solution_model='nrtl' requested but no NRTL parameters exist for "
                              f"{', '.join(members)}; ideal solution used.")
        return self._nrtl_cache[members]

    def _sle_uses_gamma(self, members, k: str) -> bool:
        """Should the solid–liquid equilibrium of component k use NRTL activity coefficients?"""
        if self.sle_activity == "ideal":
            return False
        model = self._nrtl(tuple(members))
        if model is None:
            return False
        if self.sle_activity == "nrtl":
            return True
        pairs = [p for (a, b), p in model.pair_sources.items() if k in (a, b)]
        use = bool(pairs) and all("SLE" in str(p.get("fitted_to", "")) for p in pairs)
        if use:
            self.log.assume(f"Solid–liquid equilibrium of {self.comps[k].name} uses NRTL activity coefficients "
                            "from the parameter set fitted jointly to H^E and freezing points (one consistent "
                            "liquid model for the heat of mixing and the freezing equilibrium).")
        return use

    def _dCp_fus(self, k: str) -> float:
        """ΔCp of fusion at T_m (Cp,liquid − Cp,solid), used in the solid–liquid equilibrium
        equation so that it is consistent with the Kirchhoff-law enthalpies used for H(T)."""
        if self.options.get("sle_dcp", True) is False:
            return 0.0
        c = self.comps[k]
        if c.Tm is None or not (c.has_phase("liquid") and c.has_phase("solid")):
            return 0.0
        return c.cp_corr["liquid"].raw(c.Tm) - c.cp_corr["solid"].raw(c.Tm)

    def _check_stability(self, model: NRTL, T: float) -> None:
        """Warn if, at T, the NRTL parameters of a pair would predict a (false) liquid–liquid split.

        Every pair in a liquid phase was judged miscible, so a non-convex ΔG_mix/RT means the binary
        parameters are being used outside their reliable range (Cp^E is then less trustworthy).
        """
        key = (tuple(model.keys), round(T / 5.0) * 5.0)
        if key in self._stability_checked:
            return
        self._stability_checked.add(key)
        for i in range(len(model.keys)):
            for j in range(i + 1, len(model.keys)):
                if (model.keys[i], model.keys[j]) not in model.pair_sources:
                    continue
                if not model.sub_binary(i, j).is_stable_binary(key[1], 400):
                    a, b = sorted((self.comps[model.keys[i]].name, self.comps[model.keys[j]].name))
                    self.log.warn(f"The NRTL parameters of {a} + {b} predict a spurious liquid–liquid split near "
                                  f"{key[1]:.0f} K (outside their reliable range); Cp^E for this pair is uncertain.")

    def _ions(self, k: str) -> float:
        return float(self.comps[k].raw.get("ions_in_solution", {}).get("value", 1))

    def _solve_sle(self, gname: str, T: float) -> dict[str, float]:
        """Amounts (mol) of each group member in the liquid at T (SLE iteration)."""
        members = self.groups[gname]
        nL = {k: self.amounts[k] for k in members}
        sle = [k for k in members if self.plans[k].behaviour in ("sle_ideal", "sle_table")]
        if not sle:
            return nL
        for _ in range(200):
            change = 0.0
            for k in sle:
                c, ntot = self.comps[k], self.amounts[k]
                p = self.plans[k]
                if p.behaviour == "sle_ideal":
                    x_sat = phase_change.ideal_solubility(T, c.Tm, c.dHfus, self._dCp_fus(k))
                    if x_sat < 1.0 and self._sle_uses_gamma(members, k):
                        model = self._nrtl(tuple(members))
                        if model is not None:
                            tot = sum(nL.values())
                            x = [nL[m] / tot for m in members]
                            x_sat = min(1.0, x_sat / float(model.gamma(x, T)[members.index(k)]))
                    if x_sat >= 1.0:
                        new = ntot
                    else:
                        others = sum(nL[j] * self._ions(j) for j in members if j != k)
                        new = min(ntot, x_sat / (1.0 - x_sat) * others)
                else:
                    sol = phase_change.TableSolubility(self.db.solubility_entry(k, p.solvent))
                    if not (sol.Tmin <= T <= sol.Tmax):
                        self.log.warn(f"Solubility of {c.name} in {p.solvent}: T = {T:.1f} K outside the table "
                                      f"[{sol.Tmin:.0f}, {sol.Tmax:.0f}] K (extrapolated).")
                    m_solvent = nL[p.solvent] * self.comps[p.solvent].MW
                    new = min(ntot, sol.mass_ratio(T) * m_solvent / c.MW)
                change = max(change, abs(new - nL[k]) / max(ntot, 1e-30))
                nL[k] = new
            if change < 1e-13:
                break
        else:
            self.log.warn(f"Solid–liquid equilibrium in {gname} did not converge at {T:.2f} K "
                          "(possibly below a eutectic point).")
        crystallising = [k for k in sle if nL[k] < self.amounts[k] * (1 - 1e-9)]
        if len(crystallising) > 1:
            self.log.warn(f"At {T:.1f} K more than one component crystallises from {gname} "
                          f"({', '.join(crystallising)}): below a eutectic, results are approximate.")
        return nL

    def split(self, T: float) -> list[tuple[str, str, str, float]]:
        """Phase state at T: list of (phase name, component, property state, mol)."""
        entries: list[tuple[str, str, str, float]] = []
        liquid_amounts: dict[str, float] = {}
        for gname in self.groups:
            liquid_amounts.update(self._solve_sle(gname, T))
        for k, p in self.plans.items():
            n = p.n
            if p.behaviour.startswith("fixed_"):
                ph = p.behaviour[6:]
                name = p.group if ph == "liquid" else ("gas" if ph == "gas" else f"solid ({p.comp.name})")
                entries.append((name, k, ph, n))
            elif p.behaviour == "dissolved":
                state = "liquid" if p.comp.has_phase("liquid") else "dissolved"
                entries.append((p.group, k, state, n))
            elif p.behaviour == "melt_smooth":
                f = phase_change.melt_fraction(T, p.comp.Tm, self.width)
                entries += [(f"solid ({p.comp.name})", k, "solid", (1 - f) * n), (p.group, k, "liquid", f * n)]
            elif p.behaviour == "boil_smooth":
                f = phase_change.melt_fraction(T, p.comp.Tb_at(self.P), self.width)
                entries += [(p.group, k, "liquid", (1 - f) * n), ("gas", k, "gas", f * n)]
            elif p.behaviour == "sle_ideal":
                nL = liquid_amounts[k]
                entries += [(p.group, k, "liquid", nL), (f"solid ({p.comp.name})", k, "solid", n - nL)]
            elif p.behaviour == "sle_table":
                nL = liquid_amounts[k]
                entries += [(p.group, k, "dissolved", nL), (f"solid ({p.comp.name})", k, "solid", n - nL)]
        return entries

    # ==================================================================
    # 3. Enthalpy and heat capacity
    # ==================================================================
    def _dissolution_enthalpy(self, k: str) -> float:
        p = self.plans[k]
        if p.solvent:
            entry = self.db.solubility_entry(k, p.solvent)
            sol = phase_change.TableSolubility(entry)
            if sol.dH_sol is not None:
                return sol.dH_sol
        self.log.warn(f"No enthalpy of solution for {self.comps[k].name}; dissolution heat taken as zero.")
        return 0.0

    def _pure_H(self, k: str, state: str, T: float) -> float:
        key = (k, state, round(T, 10))
        if key not in self._H_cache:
            c = self.comps[k]
            if state == "dissolved":
                h = pure_enthalpy(c, "solid", T, self.log) + self._dissolution_enthalpy(k)
            else:
                h = pure_enthalpy(c, state, T, self.log)
            h += self._p_corr(c, state, T, "H")
            self._H_cache[key] = h
        return self._H_cache[key]

    def _p_corr(self, c: Component, state: str, T: float, what: str) -> float:
        """Pressure correction of a pure-component property (correlations are for 1 atm)."""
        if abs(self.P - P_ATM) < 1.0:
            return 0.0
        phase = "solid" if state == "dissolved" else state
        try:
            if what == "H":
                return c.h_pressure_correction(phase, T, self.P)
            return c.cp_pressure_correction(phase, T, self.P)
        except DataError as exc:
            self.log.warn(f"No pressure correction for {c.name} ({phase}): {exc} Cp taken at 1 atm.")
            return 0.0

    def _pure_cp(self, k: str, state: str, T: float) -> float:
        c = self.comps[k]
        if state == "dissolved":
            if not any(self.db.apparent_cp_entry(k, j) for j in self.amounts if j != k):
                self.log.assume(f"Dissolved {c.name} contributes its solid-state Cp (Neumann–Kopp approximation); "
                                "for electrolytes the true partial molar Cp can differ strongly.")
            return c.cp("solid", T, self.log) + self._p_corr(c, "solid", T, "Cp")
        return c.cp(state, T, self.log) + self._p_corr(c, state, T, "Cp")

    def _liquid_groups(self, entries):
        groups: dict[str, dict[str, float]] = {}
        for name, k, state, n in entries:
            if state in ("liquid", "dissolved") and n > 0:
                groups.setdefault(name, {})
                groups[name][k] = groups[name].get(k, 0.0) + n
        return groups

    def enthalpy(self, T: float, entries=None) -> float:
        """Total enthalpy (J) of the basis amount at T for a given phase state."""
        entries = self.split(T) if entries is None else entries
        H = sum(n * self._pure_H(k, state, T) for _, k, state, n in entries if n > 0)
        for gname, comp_n in self._liquid_groups(entries).items():
            model = self._nrtl(tuple(self.groups[gname]))
            if model is not None:
                tot = sum(comp_n.values())
                x = [comp_n.get(m, 0.0) / tot for m in self.groups[gname]]
                H += tot * model.hE(x, T)
        H += self._electrolyte_terms(T, entries)[1]
        return H

    def _electrolyte_terms(self, T: float, entries) -> tuple[float, float]:
        """(Cp, H) corrections for dissolved electrolytes with an apparent-molar-Cp model.

        Cp(solution) = n_w Cp_w + n_s φCp(m). The ideal part already counts the salt with its
        solid Cp (Neumann–Kopp), so the correction is n_s (φCp − Cp_s). The enthalpy gets the
        matching term n_s [φCp (T − T_ref) − (H_s(T) − H_s(T_ref))], whose T-derivative is that
        correction, so point and interval-mean Cp stay consistent.
        """
        dcp = dh = 0.0
        for name, k, state, n in entries:
            if state != "dissolved" or n <= 0:
                continue
            for solvent in self.groups.get(name, []):
                e = self.db.apparent_cp_entry(k, solvent)
                if e is None:
                    continue
                n_solv = sum(m for nm, kk, st, m in entries if nm == name and kk == solvent and st == "liquid")
                if n_solv <= 0:
                    continue
                mol = n / (n_solv * self.comps[solvent].MW)          # molality, mol/kg
                c = e["coeffs"]
                phi = c["c0"] + c["d0"] * mol ** 0.5 + c["e"] * mol
                Tref = e.get("T_ref", 298.15)
                comp = self.comps[k]
                self.log.assume(f"Dissolved {comp.name}: apparent molar heat capacity φCp(m) fitted to measured "
                                f"Cp of {comp.name}(aq) at {Tref} K ({e.get('source', '')[:60]}...).")
                if abs(T - Tref) > 10.0:
                    self.log.warn(f"φCp of {comp.name}(aq) is known only at {Tref} K; used at {T:.1f} K without a "
                                  "temperature correction.")
                if mol > e.get("m_max", 1e9) * 1.05:
                    self.log.warn(f"{comp.name} molality {mol:.2f} mol/kg exceeds the fitted range "
                                  f"(≤ {e['m_max']:.1f} mol/kg); φCp extrapolated.")
                dcp += n * (phi - comp.cp("solid", T))
                dh += n * (phi * (T - Tref) - (pure_enthalpy(comp, "solid", T) - pure_enthalpy(comp, "solid", Tref)))
                break
        return dcp, dh

    def cp_frozen(self, T: float, entries=None) -> tuple[float, float]:
        """(ideal, excess) heat capacity in J/K at fixed phase amounts."""
        entries = self.split(T) if entries is None else entries
        ideal = sum(n * self._pure_cp(k, state, T) for _, k, state, n in entries if n > 0)
        excess = 0.0
        for gname, comp_n in self._liquid_groups(entries).items():
            model = self._nrtl(tuple(self.groups[gname]))
            if model is not None:
                model.check_temperature(T)
                self._check_stability(model, T)
                tot = sum(comp_n.values())
                x = [comp_n.get(m, 0.0) / tot for m in self.groups[gname]]
                excess += tot * model.cpE(x, T)
        excess += self._electrolyte_terms(T, entries)[0]
        return ideal, excess

    @property
    def has_phase_change(self) -> bool:
        return any(p.behaviour in PHASE_CHANGE_BEHAVIOURS for p in self.plans.values())

    def _ensure(self, T_lo: float, T_hi: float) -> None:
        if self.window is None or not (self.window[0] <= T_lo and T_hi <= self.window[1]):
            self.configure(T_lo, T_hi)

    def cp_point(self, T: float) -> CpBreakdown:
        """Apparent heat capacity at a single temperature T (J/K for the basis)."""
        self._ensure(T, T)
        entries = self.split(T)
        ideal, excess = self.cp_frozen(T, entries)
        latent = 0.0
        if self.has_phase_change:
            h = min(1e-3, self.width / 200.0)
            latent = (self.enthalpy(T, self.split(T + h)) - self.enthalpy(T, self.split(T - h))) / (2 * h)
        bd = CpBreakdown(ideal, excess, latent, self.mass, self.moles, "point", T=T,
                         phase_amounts=self._amounts_table(entries))
        bd.extras.update(self._extras(T, entries))
        # A point value inside a smoothed melting/boiling step depends on the assumed width (peak ∝ 1/width).
        # Width-independent quantities: the latent heat of the transition, and the mean Cp over any fixed
        # interval that contains the whole transition. Report both.
        for p in self.plans.values():
            Tt = p.comp.Tm if p.behaviour == "melt_smooth" else (p.comp.Tb_at(self.P) if p.behaviour == "boil_smooth" else None)
            if Tt is not None and abs(T - Tt) < self.width:
                dh = p.comp.dHfus if p.behaviour == "melt_smooth" else p.comp.dHvap
                latent_per_kg = p.n * dh / self.mass
                half = max(5.0, self.width)                   # fixed window (contains the whole step if width ≤ 5 K)
                T1, T2 = Tt - half, Tt + half
                mean = (self.enthalpy(T2) - self.enthalpy(T1)) / (T2 - T1) / self.mass
                bd.extras.update(transition_latent_J_per_kg=latent_per_kg, transition_mean_cp_mass=mean,
                                 transition_range_K=(T1, T2))
                self.log.warn(f"T = {T:.2f} K lies inside the smoothed phase-change range of {p.comp.name} "
                              f"(width {self.width:g} K), so this point apparent Cp depends on the assumed width. "
                              f"Width-independent: latent heat {latent_per_kg / 1000:,.1f} kJ per kg of mixture; "
                              f"mean Cp over [{T1:.2f}, {T2:.2f}] K = {mean:,.0f} J/(kg K). For energy balances "
                              "give T1 and T2 spanning the whole transition.")
        return bd

    def cp_mean(self, T1: float, T2: float) -> CpBreakdown:
        """Mean apparent heat capacity over [T1, T2] = ΔH/ΔT (J/K for the basis)."""
        if abs(T2 - T1) < 1e-9:
            raise InputError("T1 and T2 are equal; a mean Cp needs a finite temperature interval.")
        self._ensure(min(T1, T2), max(T1, T2))
        dT = T2 - T1
        dH = self.enthalpy(T2) - self.enthalpy(T1)
        a, b = min(T1, T2), max(T1, T2)
        breaks = []
        for p in self.plans.values():
            for Tt in (p.comp.Tm, p.comp.Tb_at(self.P)):
                if Tt is not None and a < Tt < b:
                    breaks.append(Tt)
        # absolute tolerance relative to the size of the integral (J); the excess
        # integrand can be ~0, where a purely relative tolerance is unreachable
        scale = abs(dH) + self.mass * 1000.0 * abs(dT)
        opts = dict(limit=400, epsrel=1e-8, epsabs=1e-10 * scale, points=sorted(breaks) or None)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", integrate.IntegrationWarning)
            ideal_int, _ = integrate.quad(lambda T: self.cp_frozen(T)[0], a, b, **opts)
            excess_int, _ = integrate.quad(lambda T: self.cp_frozen(T)[1], a, b, **opts)
        if caught:
            self.log.note("Numerical integration of the ideal/excess parts reported limited accuracy "
                          "(kinks at phase boundaries); the total mean Cp = ΔH/ΔT is unaffected.")
        sign = 1.0 if dT > 0 else -1.0
        ideal_m, excess_m = sign * ideal_int / dT, sign * excess_int / dT
        latent_m = dH / dT - ideal_m - excess_m
        if not self.has_phase_change:
            latent_m = 0.0 if abs(latent_m) < 1e-6 * abs(dH / dT) else latent_m
        bd = CpBreakdown(ideal_m, excess_m, latent_m, self.mass, self.moles, "mean", T1=T1, T2=T2,
                         phase_amounts={f"at T1 = {T1:.2f} K": self._flat_amounts(self.split(T1)),
                                        f"at T2 = {T2:.2f} K": self._flat_amounts(self.split(T2))})
        bd.extras["dH_J"] = dH
        return bd

    # ------------------------------------------------------------------
    def _flat_amounts(self, entries) -> dict[str, float]:
        out: dict[str, float] = {}
        for name, k, state, n in entries:
            key = f"{name}: {self.comps[k].name}"
            out[key] = out.get(key, 0.0) + n * self.comps[k].MW / self.mass
        return out

    def _amounts_table(self, entries) -> dict[str, dict[str, float]]:
        """Mass fraction (of whole system) of each component in each phase."""
        out: dict[str, dict[str, float]] = {}
        for name, k, state, n in entries:
            out.setdefault(name, {})
            out[name][self.comps[k].name] = out[name].get(self.comps[k].name, 0.0) + n * self.comps[k].MW / self.mass
        return out

    def _extras(self, T: float, entries) -> dict[str, Any]:
        """Extra diagnostics: Model I comparison for suspensions, Cv estimate."""
        ex: dict[str, Any] = {}
        try:
            rho, alpha, kappa = self.volumetric_properties(T, entries)
            ex["rho_mix"], ex["alpha_mix"], ex["kappa_T_mix"] = rho, alpha, kappa
        except Exception:                         # missing data is not fatal for extras
            pass
        solids = {name for name, _, st, n in entries if st == "solid" and n > 0}
        liquids = {name for name, _, st, n in entries if st in ("liquid", "dissolved") and n > 0}
        if solids and len(liquids) == 1 and not self.has_phase_change:
            try:
                vol, cpm = {}, {}
                for name, k, state, n in entries:
                    if n <= 0:
                        continue
                    c = self.comps[k]
                    ph = "solid" if state in ("solid", "dissolved") else state
                    vol[name] = vol.get(name, 0.0) + n * c.MW / c.density(ph, T, self.P)
                    cpm.setdefault(name, [0.0, 0.0])
                    cpm[name][0] += n * self._pure_cp(k, state, T)
                    cpm[name][1] += n * c.MW
                V = sum(vol.values())
                cp_I = sum(vol[p] / V * cpm[p][0] / cpm[p][1] for p in vol)
                ex["model_I_cp_mass"] = cp_I
                ex["phase_volume_fractions"] = {p: vol[p] / V for p in vol}
            except Exception:
                pass
        return ex

    def volumetric_properties(self, T: float, entries=None) -> tuple[float, float, float]:
        """(ρ, α, κ_T) of the mixture assuming ideal volume of mixing."""
        entries = self.split(T) if entries is None else entries
        w, rho, alpha, kappa = {}, {}, {}, {}
        for i, (name, k, state, n) in enumerate(entries):
            if n <= 0:
                continue
            c = self.comps[k]
            ph = "solid" if state == "dissolved" else state
            tag = f"{i}"
            w[tag] = n * c.MW / self.mass
            rho[tag] = c.density(ph, T, self.P, self.log)
            alpha[tag] = c.alpha(ph, T, self.P, self.log)
            kappa[tag] = c.kappa_T(ph, T, self.P)
        s = sum(w.values())
        w = {k: v / s for k, v in w.items()}
        return cv_cp.mixture_volumetric_properties(w, rho, alpha, kappa)

    # ------------------------------------------------------------------
    def classification(self) -> str:
        """Human-readable system type (ideal / real solution / heterogeneous ...)."""
        if not self.plans:
            raise RuntimeError("configure() first")
        n_liq = len(self.groups)
        has_solid = any(p.behaviour in ("fixed_solid", "sle_ideal", "sle_table", "melt_smooth")
                        for p in self.plans.values())
        gases = [k for k, p in self.plans.items() if p.behaviour in ("fixed_gas", "boil_smooth")]
        parts = []
        if len(self.plans) == 1:
            parts.append("Pure component")
        else:
            liquid_desc = []
            for gname, members in self.groups.items():
                if len(members) == 1:
                    kind = "pure liquid"
                elif self._nrtl(tuple(members)) is not None:
                    kind = "real liquid solution (NRTL Cp^E)"
                else:
                    kind = "ideal liquid mixture"
                liquid_desc.append(f"{gname}: {kind}")
            if n_liq == 1 and not has_solid and not gases:
                parts.append("Homogeneous liquid solution")
            elif n_liq == 0 and not has_solid:
                parts.append("Ideal-gas mixture")
            else:
                features = []
                if has_solid:
                    features.append("solid particles (suspension / slurry)" if n_liq else "solid phases")
                if n_liq >= 2:
                    features.append("immiscible liquid phases (emulsion)")
                if gases and (n_liq or has_solid):
                    features.append("gas phase")
                parts.append("Heterogeneous mixture: " + " + ".join(features))
            parts += liquid_desc
        if self.has_phase_change:
            parts.append("with phase change (apparent Cp)")
        return "; ".join(parts)

    def behaviour_table(self) -> list[dict[str, str]]:
        return [{"component": p.comp.name, "behaviour": BEHAVIOUR_TEXT[p.behaviour],
                 "phase": p.group or "-", "reason": p.reason} for p in self.plans.values()]
