"""Thermodynamic relations as rule objects.

Every relation the engine knows is a `Rule`:

    inputs (symbols that must be known)  ──►  outputs (symbols it produces)

together with its name, the governing equation and the route it belongs to.
The solver (solver.py) applies rules by *forward chaining* until no rule can
fire.  Keeping each relation small and separate makes the derivation fully
traceable: the solver can print exactly which equations led to Cp.

Routes
------
A  property route   : Cp from T, P, composition and pure-component data
B  calorimetric route: Cp from heat Q and work W (first law)
C  derived routes    : Cp from ΔH, ΔU, ΔS, Cv or γ
helper              : unit-free bookkeeping (ΔT from T1 and T2, mass from moles, ...)

Sign convention (IUPAC): Q > 0 heat added to the system, W > 0 work done on
the system, ΔU = Q + W.  W is shaft/stirring (non-expansion) work; expansion
work against the surroundings is W_pv = −∫P dV.  At constant P the first law
becomes ΔH = Q + W (because ΔH = ΔU + PΔV and W_pv = −PΔV).
"""
from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

from .models import cv_cp
from .system import SystemModel

# ---------------------------------------------------------------------------
# Symbols: name -> (description shown to the user, SI unit)
# ---------------------------------------------------------------------------
SYMBOLS: dict[str, tuple[str, str]] = {
    "T": ("temperature T", "K"),
    "T1": ("initial temperature T1", "K"),
    "T2": ("final temperature T2", "K"),
    "dT": ("temperature change ΔT", "K"),
    "P": ("pressure P", "Pa"),
    "composition": ("composition (components with fractions or amounts)", "-"),
    "m": ("total mass of the system", "kg"),
    "n": ("total amount of the system (moles)", "mol"),
    "M": ("mean molar mass of the mixture", "kg/mol"),
    "Q": ("heat Q added to the system", "J"),
    "W": ("shaft / stirring work W done on the system", "J"),
    "W_pv": ("expansion work W_pv done on the system", "J"),
    "dH": ("enthalpy change ΔH", "J"),
    "dU": ("internal-energy change ΔU", "J"),
    "dS": ("entropy change ΔS", "J/K"),
    "dV": ("volume change ΔV", "m^3"),
    "Cv_mass": ("isochoric heat capacity Cv (per kg)", "J/(kg K)"),
    "Cv_molar": ("isochoric heat capacity Cv (per mol)", "J/(mol K)"),
    "Cv_total": ("isochoric heat capacity Cv of the whole system", "J/K"),
    "gamma": ("heat-capacity ratio γ = Cp/Cv", "-"),
    "alpha": ("thermal expansion coefficient α", "1/K"),
    "kappa_T": ("isothermal compressibility κ_T", "1/Pa"),
    "rho": ("density ρ", "kg/m^3"),
    "constant_P": ("process = constant pressure", "-"),
    "constant_V": ("process = constant volume", "-"),
    "adiabatic": ("process = adiabatic", "-"),
    "system": ("classified mixture model (phases, solution models)", "-"),
}

# Symbols a user can supply directly (used for 'what is missing' suggestions)
BASE_SYMBOLS = {"T", "T1", "T2", "dT", "P", "composition", "m", "n", "Q", "W", "W_pv", "dH", "dU",
                "dS", "dV", "Cv_mass", "Cv_molar", "Cv_total", "gamma", "alpha", "kappa_T", "rho",
                "constant_P", "constant_V", "adiabatic"}

CP_TAGS = {
    "A:point": "Route A – property route (pure-component data + mixing rules), at T",
    "A:mean": "Route A – property route, mean over [T1, T2] = ΔH_model/ΔT",
    "dH": "mean Cp from ΔH/ΔT",
    "Cv": "Cp from Cv",
    "dS": "Route C – from ΔS over [T1, T2]",
    "gamma": "Route C – from γ = Cp/Cv",
}


def cp_symbol(tag: str, basis: str = "mass") -> str:
    return f"Cp_{basis}[{tag}]"


@dataclass
class Rule:
    """One thermodynamic relation."""
    name: str
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    func: Callable[[dict[str, Any], Any], Any]
    equation: str
    route: str
    description: str = ""
    check: bool = False          # evaluate even if outputs are known, to test consistency
    check_level: str = "warn"    # "warn" or "note" for inconsistencies

    def apply(self, values: dict[str, Any], ctx: Any) -> tuple:
        out = self.func(values, ctx)
        return out if isinstance(out, tuple) else (out,)


# ===========================================================================
# Helper rules
# ===========================================================================
def _delta_T(v, ctx):
    """ΔT = T2 − T1"""
    return v["T2"] - v["T1"]


def _mean_T(v, ctx):
    """T = (T1 + T2)/2  (mid-point temperature of the interval)"""
    return 0.5 * (v["T1"] + v["T2"])


def _interval(v, ctx):
    """T1 = T − ΔT/2,  T2 = T + ΔT/2  (interval centred on T)"""
    return v["T"] - 0.5 * v["dT"], v["T"] + 0.5 * v["dT"]


def _molar_mass(v, ctx):
    """M = Σ n_i M_i / Σ n_i"""
    a = ctx.inp.amounts
    return sum(a[k] * ctx.db.get(k).MW for k in a) / sum(a.values())


def _m_from_n(v, ctx):
    """m = n M"""
    return v["n"] * v["M"]


def _n_from_m(v, ctx):
    """n = m / M"""
    return v["m"] / v["M"]


def _cv_mass_from_molar(v, ctx):
    """cv = Cv,m / M"""
    return v["Cv_molar"] / v["M"]


def _cv_mass_from_total(v, ctx):
    """cv = Cv,total / m"""
    return v["Cv_total"] / v["m"]


def _classify(v, ctx):
    """Build the mixture model: infer phases, group liquids, choose ideal / NRTL."""
    k = ctx.known_values()
    temps = [k[s] for s in ("T", "T1", "T2") if s in k]
    s = SystemModel(ctx.db, ctx.inp.amounts, ctx.inp.phase_spec, k.get("P", 101325.0),
                    ctx.inp.options, ctx.log)
    s.configure(min(temps), max(temps))
    return s


# ===========================================================================
# Route B – calorimetric (first law)
# ===========================================================================
def _adiabatic_Q(v, ctx):
    """Adiabatic process: Q = 0"""
    return 0.0


def _dH_constP(v, ctx):
    """Constant P: ΔH = Q + W  (W = shaft/stirring work; expansion work is absorbed in H)"""
    return v["Q"] + v["W"]


def _dU_constV(v, ctx):
    """Constant V: ΔU = Q + W  (no expansion work)"""
    return v["Q"] + v["W"]


def _dU_first_law(v, ctx):
    """First law: ΔU = Q + W + W_pv"""
    return v["Q"] + v["W"] + v["W_pv"]


def _cp_from_dH_mass(v, ctx):
    """cp = ΔH / (m ΔT)   (mean over the interval; initial and final P equal)"""
    ctx.log.assume("ΔH-based Cp: initial and final states at the same pressure; result is the mean Cp over ΔT.")
    return v["dH"] / (v["m"] * v["dT"])


def _cp_from_dH_molar(v, ctx):
    """Cp,m = ΔH / (n ΔT)"""
    return v["dH"] / (v["n"] * v["dT"])


def _cv_from_dU(v, ctx):
    """Constant V: cv = ΔU / (m ΔT)"""
    return v["dU"] / (v["m"] * v["dT"])


def _cv_from_dU_molar(v, ctx):
    """Constant V: Cv,m = ΔU / (n ΔT)"""
    return v["dU"] / (v["n"] * v["dT"])


# ===========================================================================
# Route C – derived properties
# ===========================================================================
def _dH_from_dU(v, ctx):
    """Constant P: ΔH = ΔU + P ΔV"""
    return v["dU"] + v["P"] * v["dV"]


def _cp_from_dS(v, ctx):
    """Constant P, constant cp, no phase change: ΔS = m cp ln(T2/T1)  ⇒  cp = ΔS / (m ln(T2/T1))"""
    ctx.log.assume("ΔS route: constant pressure, Cp constant over [T1, T2], no phase change in the interval.")
    if ctx.is_known("system") and ctx.value("system").has_phase_change:
        ctx.log.warn("ΔS route assumes no phase change, but the mixture model predicts one in [T1, T2]; "
                     "the ΔS-based Cp is not meaningful there.")
    return v["dS"] / (v["m"] * math.log(v["T2"] / v["T1"]))


def _cp_from_dS_molar(v, ctx):
    """Cp,m = ΔS / (n ln(T2/T1))"""
    return v["dS"] / (v["n"] * math.log(v["T2"] / v["T1"]))


def _cp_from_cv(v, ctx):
    """cp = cv + T v α² / κ_T,  v = 1/ρ   (exact thermodynamic identity)"""
    return cv_cp.cp_from_cv(v["Cv_mass"], v["T"], 1.0 / v["rho"], v["alpha"], v["kappa_T"])


def _cp_from_gamma(v, ctx):
    """cp = γ cv"""
    return v["gamma"] * v["Cv_mass"]


def _volumetric(v, ctx):
    """Mixture ρ, α, κ_T with ideal volume of mixing: α = Σφ_i α_i, κ_T = Σφ_i κ_i, 1/ρ = Σ w_i/ρ_i"""
    ctx.log.assume("Mixture α, κ_T and ρ from volume-fraction averaging (ideal volume of mixing, V^E = 0).")
    rho, alpha, kappa = v["system"].volumetric_properties(v["T"])
    return alpha, kappa, rho


# ===========================================================================
# Route A – property route
# ===========================================================================
def _cp_point(v, ctx):
    """Cp(T) = Σ n_i Cp_i + Σ n_p Cp^E_p + Σ Δh (dn/dT)   (ideal + excess + latent)"""
    bd = v["system"].cp_point(v["T"])
    ctx.breakdowns["A:point"] = bd
    return bd.total / bd.mass, bd.total / bd.moles


def _cp_mean(v, ctx):
    """Cp_mean = [H(T2) − H(T1)] / (m (T2 − T1))   (includes latent heat in the interval)"""
    bd = v["system"].cp_mean(v["T1"], v["T2"])
    ctx.breakdowns["A:mean"] = bd
    return bd.total / bd.mass, bd.total / bd.moles


# ===========================================================================
def build_rules() -> list[Rule]:
    """All relations known to the engine (order = firing priority)."""
    m, mo = cp_symbol, lambda t: cp_symbol(t, "molar")
    return [
        # helpers
        Rule("delta_T", ("T1", "T2"), ("dT",), _delta_T, "ΔT = T2 − T1", "helper"),
        Rule("mean_temperature", ("T1", "T2"), ("T",), _mean_T, "T = (T1 + T2)/2", "helper"),
        Rule("interval_from_T_and_dT", ("T", "dT"), ("T1", "T2"), _interval,
             "T1 = T − ΔT/2 ;  T2 = T + ΔT/2", "helper"),
        Rule("mean_molar_mass", ("composition",), ("M",), _molar_mass, "M = Σ x_i M_i", "helper"),
        Rule("mass_from_moles", ("n", "M"), ("m",), _m_from_n, "m = n M", "helper"),
        Rule("moles_from_mass", ("m", "M"), ("n",), _n_from_m, "n = m / M", "helper"),
        Rule("cv_per_mass_from_molar", ("Cv_molar", "M"), ("Cv_mass",), _cv_mass_from_molar, "cv = Cv,m / M", "helper"),
        Rule("cv_per_mass_from_total", ("Cv_total", "m"), ("Cv_mass",), _cv_mass_from_total, "cv = Cv / m", "helper"),
        Rule("classify_system", ("composition", "T"), ("system",), _classify,
             "phases from Tm, Tb(P), miscibility, solubility; ideal vs NRTL per liquid phase", "A"),
        # Route B
        Rule("adiabatic_no_heat", ("adiabatic",), ("Q",), _adiabatic_Q, "Q = 0", "B", check=True),
        Rule("first_law_constant_P", ("Q", "W", "constant_P"), ("dH",), _dH_constP, "ΔH = Q + W", "B", check=True),
        Rule("first_law_constant_V", ("Q", "W", "constant_V"), ("dU",), _dU_constV, "ΔU = Q + W", "B", check=True),
        Rule("first_law_general", ("Q", "W", "W_pv"), ("dU",), _dU_first_law, "ΔU = Q + W + W_pv", "B", check=True),
        Rule("enthalpy_from_internal_energy", ("dU", "P", "dV", "constant_P"), ("dH",), _dH_from_dU,
             "ΔH = ΔU + P ΔV", "C", check=True),
        Rule("cv_from_dU", ("dU", "m", "dT", "constant_V"), ("Cv_mass",), _cv_from_dU, "cv = ΔU / (m ΔT)", "B", check=True),
        Rule("cv_from_dU_molar", ("dU", "n", "dT", "constant_V"), ("Cv_molar",), _cv_from_dU_molar,
             "Cv,m = ΔU / (n ΔT)", "B", check=True),
        Rule("cp_from_dH", ("dH", "m", "dT"), (m("dH"),), _cp_from_dH_mass, "cp = ΔH / (m ΔT)", "B/C"),
        Rule("cp_from_dH_molar", ("dH", "n", "dT"), (mo("dH"),), _cp_from_dH_molar, "Cp,m = ΔH / (n ΔT)", "B/C"),
        # Route C
        Rule("mixture_volumetric_properties", ("system", "T"), ("alpha", "kappa_T", "rho"), _volumetric,
             "α = Σφ_i α_i ; κ_T = Σφ_i κ_T,i ; 1/ρ = Σ w_i/ρ_i", "C", check=True, check_level="note"),
        Rule("cp_from_cv", ("Cv_mass", "T", "alpha", "kappa_T", "rho"), (m("Cv"),), _cp_from_cv,
             "cp = cv + T v α² / κ_T", "C"),
        Rule("cp_from_gamma", ("Cv_mass", "gamma"), (m("gamma"),), _cp_from_gamma, "cp = γ cv", "C"),
        Rule("cp_from_dS", ("dS", "m", "T1", "T2"), (m("dS"),), _cp_from_dS, "cp = ΔS / (m ln(T2/T1))", "C"),
        Rule("cp_from_dS_molar", ("dS", "n", "T1", "T2"), (mo("dS"),), _cp_from_dS_molar,
             "Cp,m = ΔS / (n ln(T2/T1))", "C"),
        # Route A
        Rule("cp_property_point", ("system", "T"), (m("A:point"), mo("A:point")), _cp_point,
             "Cp = Σ n_i Cp_i + Cp^E + Σ Δh dn/dT", "A"),
        Rule("cp_property_mean", ("system", "T1", "T2"), (m("A:mean"), mo("A:mean")), _cp_mean,
             "Cp_mean = [H(T2) − H(T1)] / (m ΔT)", "A"),
    ]
