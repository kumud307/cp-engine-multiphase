"""Input layer and normaliser.

The user describes the system with a plain Python dict (or the SystemInput
dataclass).  Every key is optional – the solver later decides which
thermodynamic routes the given information allows.

Accepted keys (synonyms in brackets)
------------------------------------
State          T [temperature], T1 [T_initial], T2 [T_final], dT [delta_T], P [pressure]
Components     components: list of names, or dict name -> phase
               phases: dict name -> "solid" | "liquid" | "gas" | "dissolved" | "auto"
Composition    exactly ONE of
                 mole_fractions [x], mass_fractions [w], volume_fractions,
                 moles, masses               (dict name -> value / quantity)
               optionally phi (solid volume fractions of the WHOLE mixture);
               with phi, the fraction dict describes only the fluid.
               One missing fraction (value None or not listed) is filled by Σ = 1.
Amount         total_mass [m, mass], total_moles [n]
Process        process: "constant_P" [isobaric], "constant_V" [isochoric],
               "adiabatic", "adiabatic_constant_P", or a list of these
Energies       Q [heat], W [work, W_shaft], W_pv [W_expansion], dH, dU, dS, dV
               (totals, or per kg / per mol – e.g. "2.1 kJ/kg")
Properties     Cv, gamma [Cp/Cv], alpha, kappa_T, density [rho]
Options        melt_width, solution_model ("auto"|"ideal"|"nrtl"),
               sle_activity ("auto"|"ideal"|"nrtl"), liquid_phases, tolerance,
               nrtl_parameters ("HE" = H^E-fitted set, default | "VLE"),
               vle_parameters ("VLE", default | "HE": NRTL set for boiling of mixtures),
               vapour_model ("eos", default: φ from gas_eos | "ideal" | "lewis_randall": virial φ),
               gas_eos ("PR", default | "SRK" | "RK" | "vdW" | "virial" | "ideal": residual Cp of gases)

Sign convention (IUPAC): Q > 0 is heat ADDED TO the system; W > 0 is work
DONE ON the system.  First law: ΔU = Q + W.

Quantities are given as numbers (SI assumed), (value, "unit") tuples or
"value unit" strings – see cp_engine.units.

Canonical form produced here
----------------------------
* every quantity in SI
* composition as moles of each component (per declared phase);
  if only fractions are known the amounts refer to a calculation basis
  (1 mol, 1 kg or 1 m^3) and `amount_known` is False
"""
from __future__ import annotations

import difflib
from dataclasses import asdict, dataclass, field
from typing import Any

from .database import P_ATM, Database
from .errors import InputError, UnitError
from .report import Log
from .system import pure_phase_at
from .units import to_si

SYNONYMS = {
    "temperature": "T", "t": "T", "T": "T",
    "T1": "T1", "t1": "T1", "T_initial": "T1", "T_start": "T1",
    "T2": "T2", "t2": "T2", "T_final": "T2", "T_end": "T2",
    "dT": "dT", "delta_T": "dT", "deltaT": "dT", "ΔT": "dT",
    "P": "P", "pressure": "P", "p": "P",
    "components": "components", "phases": "phases",
    "mole_fractions": "mole_fractions", "x": "mole_fractions",
    "mass_fractions": "mass_fractions", "w": "mass_fractions",
    "volume_fractions": "volume_fractions",
    "moles": "moles", "masses": "masses", "phi": "phi",
    "total_mass": "total_mass", "m": "total_mass", "mass": "total_mass",
    "total_moles": "total_moles", "n": "total_moles",
    "process": "process",
    "Q": "Q", "heat": "Q",
    "W": "W", "work": "W", "W_shaft": "W", "W_stirring": "W",
    "W_pv": "W_pv", "W_expansion": "W_pv",
    "dH": "dH", "delta_H": "dH", "ΔH": "dH",
    "dU": "dU", "delta_U": "dU", "ΔU": "dU",
    "dS": "dS", "delta_S": "dS", "ΔS": "dS",
    "dV": "dV", "delta_V": "dV",
    "Cv": "Cv", "cv": "Cv",
    "gamma": "gamma", "Cp/Cv": "gamma",
    "alpha": "alpha", "kappa_T": "kappa_T",
    "density": "density", "rho": "density",
    "melt_width": "melt_width", "solution_model": "solution_model",
    "sle_activity": "sle_activity", "liquid_phases": "liquid_phases",
    "nrtl_parameters": "nrtl_parameters", "sle_dcp": "sle_dcp", "gas_eos": "gas_eos", "eos": "gas_eos",
    "vapour_model": "vapour_model", "vle_parameters": "vle_parameters",
    "tolerance": "tolerance",
}
COMPOSITION_KEYS = ("mole_fractions", "mass_fractions", "volume_fractions", "moles", "masses")
PROCESS_SYNONYMS = {
    "constant_p": {"constant_P"}, "isobaric": {"constant_P"}, "constant pressure": {"constant_P"},
    "constant_v": {"constant_V"}, "isochoric": {"constant_V"}, "constant volume": {"constant_V"},
    "adiabatic": {"adiabatic"},
    "adiabatic_constant_p": {"adiabatic", "constant_P"},
}
FRACTION_TOL = 1e-3


@dataclass
class SystemInput:
    """Dataclass alternative to the dict input (same field names)."""
    T: Any = None
    T1: Any = None
    T2: Any = None
    dT: Any = None
    P: Any = None
    components: Any = None
    phases: dict | None = None
    mole_fractions: dict | None = None
    mass_fractions: dict | None = None
    volume_fractions: dict | None = None
    moles: dict | None = None
    masses: dict | None = None
    phi: dict | None = None
    total_mass: Any = None
    total_moles: Any = None
    process: Any = None
    Q: Any = None
    W: Any = None
    W_pv: Any = None
    dH: Any = None
    dU: Any = None
    dS: Any = None
    dV: Any = None
    Cv: Any = None
    gamma: Any = None
    alpha: Any = None
    kappa_T: Any = None
    density: Any = None
    melt_width: Any = None
    solution_model: str | None = None
    sle_activity: str | None = None
    nrtl_parameters: str | None = None
    vle_parameters: str | None = None
    vapour_model: str | None = None
    gas_eos: str | None = None
    liquid_phases: dict | None = None
    tolerance: float | None = None

    def to_dict(self) -> dict:
        return {k: v for k, v in asdict(self).items() if v is not None}


@dataclass
class NormalisedInput:
    known: dict[str, float] = field(default_factory=dict)      # SI scalars
    flags: set[str] = field(default_factory=set)               # process flags
    amounts: dict[str, float] | None = None                    # mol per component
    phase_spec: dict[str, str] = field(default_factory=dict)
    amount_known: bool = False
    basis: str | None = None
    options: dict[str, Any] = field(default_factory=dict)
    given: dict[str, str] = field(default_factory=dict)        # symbol -> how the user gave it
    log: Log = field(default_factory=Log)


# ---------------------------------------------------------------------------
def _canonical_keys(raw: dict) -> dict:
    out = {}
    for k, v in raw.items():
        if v is None:
            continue
        key = SYNONYMS.get(k) or SYNONYMS.get(str(k).strip())
        if key is None:
            close = difflib.get_close_matches(str(k), list(SYNONYMS), n=3)
            hint = f" Did you mean {', '.join(repr(c) for c in close)}?" if close else ""
            raise InputError(f"Unknown input '{k}'.{hint} Accepted inputs: "
                             f"{', '.join(sorted(set(SYNONYMS.values())))}.")
        if key in out:
            raise InputError(f"'{k}' given twice (as '{key}' and a synonym).")
        out[key] = v
    return out


def _fill_fractions(name: str, fr: dict[str, Any], expected: list[str] | None, log: Log) -> dict[str, float]:
    """Validate a fraction dict; fill exactly one missing value by Σ = 1."""
    fr = dict(fr)
    if expected:
        for k in expected:
            fr.setdefault(k, None)
    missing = [k for k, v in fr.items() if v is None]
    vals = {}
    for k, v in fr.items():
        if v is None:
            continue
        try:
            vals[k] = float(v)
        except (TypeError, ValueError):
            raise InputError(f"{name}: value for '{k}' is not a number ({v!r}).") from None
    if vals and all(v >= 1.0 for v in vals.values()) and abs(sum(vals.values()) - 100.0) < 0.1 and not missing:
        log.assume(f"{name} appear to be given in percent (sum = 100); divided by 100.")
        vals = {k: v / 100.0 for k, v in vals.items()}
    for k, v in vals.items():
        if v < 0 or v > 1:
            raise InputError(f"{name}: '{k}' = {v} is outside [0, 1].")
    if len(missing) > 1:
        raise InputError(f"{name}: {len(missing)} values are missing ({', '.join(missing)}); "
                         "at most one can be inferred from Σ = 1.")
    s = sum(vals.values())
    if len(missing) == 1:
        rest = 1.0 - s
        if rest < -FRACTION_TOL:
            raise InputError(f"{name}: the given fractions already sum to {s:.4f} > 1, "
                             f"so '{missing[0]}' cannot be filled.")
        vals[missing[0]] = max(rest, 0.0)
        log.note(f"{name}: '{missing[0]}' was not given and was set to {max(rest, 0.0):.6g} so that Σ = 1.")
        return vals
    if abs(s - 1.0) > FRACTION_TOL:
        raise InputError(f"{name} sum to {s:.6g}, not 1 (tolerance {FRACTION_TOL}). "
                         "Correct them, or leave exactly one out to have it filled automatically.")
    if abs(s - 1.0) > 1e-12:
        log.note(f"{name} summed to {s:.6g}; renormalised to 1.")
        vals = {k: v / s for k, v in vals.items()}
    return vals


def _process_flags(p: Any) -> set[str]:
    items = [p] if isinstance(p, str) else list(p)
    flags: set[str] = set()
    for it in items:
        key = str(it).strip().lower().replace("-", "_")
        if key not in PROCESS_SYNONYMS:
            raise InputError(f"Unknown process '{it}'. Use constant_P, constant_V, adiabatic "
                             "or adiabatic_constant_P (or a list).")
        flags |= PROCESS_SYNONYMS[key]
    return flags


# ---------------------------------------------------------------------------
def normalise(raw: dict | SystemInput, db: Database) -> NormalisedInput:
    """Convert raw user input into the canonical SI form (see module docstring)."""
    if isinstance(raw, SystemInput):
        raw = raw.to_dict()
    if not isinstance(raw, dict):
        raise InputError("Input must be a dict (or SystemInput).")
    d = _canonical_keys(raw)
    out = NormalisedInput()
    log = out.log

    # ---- state ------------------------------------------------------------
    for sym in ("T", "T1", "T2"):
        if sym in d:
            v, _, assumed = to_si(d[sym], "temperature", absolute_temperature=True, name=sym)
            out.known[sym] = v
            out.given[sym] = str(d[sym])
            if assumed:
                log.assume(f"{sym} given without unit; taken as kelvin.")
    if "dT" in d:
        v, _, assumed = to_si(d["dT"], "temperature", absolute_temperature=False, name="dT")
        out.known["dT"], out.given["dT"] = v, str(d["dT"])
        if assumed:
            log.assume("dT given without unit; taken as kelvin.")
    if "P" in d:
        v, _, assumed = to_si(d["P"], "pressure", name="P")
        out.known["P"], out.given["P"] = v, str(d["P"])
        if assumed:
            log.assume("P given without unit; taken as Pa.")
    else:
        out.known["P"] = P_ATM
        log.assume("Pressure not given; P = 1 atm (101325 Pa) assumed.")

    # temperature used for densities in volume-fraction conversions
    T_ref = out.known.get("T")
    if T_ref is None and "T1" in out.known and "T2" in out.known:
        T_ref = 0.5 * (out.known["T1"] + out.known["T2"])
    if T_ref is None:
        T_ref = out.known.get("T1", out.known.get("T2"))

    # ---- process ----------------------------------------------------------
    if "process" in d:
        out.flags = _process_flags(d["process"])
        if {"constant_P", "constant_V"} <= out.flags:
            raise InputError("A process cannot be both constant-P and constant-V "
                             "(for condensed phases choose the one that was actually controlled).")

    # ---- options ----------------------------------------------------------
    if "melt_width" in d:
        out.options["melt_width"], _, _ = to_si(d["melt_width"], "temperature", name="melt_width")
    for k in ("solution_model", "sle_activity", "liquid_phases", "tolerance", "nrtl_parameters", "sle_dcp",
              "gas_eos", "vapour_model", "vle_parameters"):
        if k in d:
            out.options[k] = d[k]
    if out.options.get("solution_model", "auto") not in ("auto", "ideal", "nrtl"):
        raise InputError("solution_model must be 'auto', 'ideal' or 'nrtl'.")
    if out.options.get("sle_activity", "auto") not in ("auto", "ideal", "nrtl"):
        raise InputError("sle_activity must be 'auto', 'ideal' or 'nrtl'.")
    if out.options.get("nrtl_parameters", "HE") not in ("HE", "VLE"):
        raise InputError("nrtl_parameters must be 'HE' (excess-enthalpy fitted) or 'VLE'.")
    if out.options.get("vle_parameters", "VLE") not in ("HE", "VLE"):
        raise InputError("vle_parameters must be 'VLE' (fitted to phase equilibrium) or 'HE'.")
    if out.options.get("vapour_model", "eos") not in ("eos", "ideal", "lewis_randall"):
        raise InputError("vapour_model must be 'eos', 'ideal' or 'lewis_randall'.")
    if out.options.get("gas_eos", "PR") not in ("PR", "SRK", "RK", "vdW", "virial", "ideal"):
        raise InputError("gas_eos must be one of PR, SRK, RK, vdW, virial, ideal.")

    # ---- components & composition -----------------------------------------
    _composition(d, db, out, T_ref)

    # ---- energies and other measured quantities ---------------------------
    _energies(d, db, out)
    return out


def _species_key(db: Database, name: str, phases: dict[str, str]) -> tuple[str, str]:
    comp, implied = db.resolve(name)
    phase = phases.get(name) or phases.get(comp.key) or implied or "auto"
    phase = str(phase).lower()
    if phase not in ("solid", "liquid", "gas", "dissolved", "auto"):
        raise InputError(f"Unknown phase '{phase}' for {name}. Use solid, liquid, gas, dissolved or auto.")
    return comp.key, phase


def _composition(d: dict, db: Database, out: NormalisedInput, T_ref: float | None) -> None:
    log = out.log
    phases_in = dict(d.get("phases") or {})
    comps_in = d.get("components")
    if isinstance(comps_in, dict):
        for k, v in comps_in.items():
            if v is not None:
                phases_in.setdefault(k, v)
        comp_names = list(comps_in)
    elif comps_in is not None:
        comp_names = [comps_in] if isinstance(comps_in, str) else list(comps_in)
    else:
        comp_names = []

    given = [k for k in COMPOSITION_KEYS if k in d]
    phi = d.get("phi")
    if len(given) > 1:
        raise InputError(f"Composition over-specified: {', '.join(given)} given. Use exactly one "
                         "(plus optionally 'phi' for solid particles).")
    if not given and phi is None:
        if len(comp_names) == 1:
            given_key, comp_dict = "mole_fractions", {comp_names[0]: 1.0}
        elif comp_names:
            raise InputError(f"Components {', '.join(comp_names)} were listed but no composition was "
                             "given (mole_fractions, mass_fractions, volume_fractions, moles, masses or phi).")
        else:
            return                                   # no composition at all (calorimetric routes only)
    else:
        given_key = given[0] if given else None
        comp_dict = dict(d[given_key]) if given_key else {}

    # names expected in the fraction dict (listed components not in phi)
    phi_names = list(phi) if phi else []
    expected = [n for n in comp_names if n not in phi_names] or None

    def T_for_density() -> float:
        if T_ref is None:
            log.assume("No temperature given; densities for volume-fraction conversion evaluated at 298.15 K.")
            return 298.15
        return T_ref

    def rho(name: str) -> float:
        key, ph = _species_key(db, name, phases_in)
        comp = db.get(key)
        T = T_for_density()
        if ph in ("auto", "dissolved"):
            ph = "solid" if ph == "dissolved" else pure_phase_at(comp, T, out.known["P"])
        return comp.density(ph, T, out.known["P"], log)

    species: dict[tuple[str, str], float] = {}     # (key, phase) -> mol

    def add(name: str, mol: float) -> None:
        if mol < 0:
            raise InputError(f"Amount of '{name}' is negative.")
        sk = _species_key(db, name, phases_in)
        species[sk] = species.get(sk, 0.0) + mol

    if given_key in ("moles", "masses"):
        if phi:
            raise InputError("'phi' can only be combined with a fraction-type composition of the fluid.")
        for name, q in comp_dict.items():
            if q is None:
                raise InputError(f"{given_key}: amount of '{name}' is missing.")
            dim = "amount" if given_key == "moles" else "mass"
            v, _, assumed = to_si(q, dim, name=f"{given_key}[{name}]")
            if assumed:
                log.assume(f"{given_key}[{name}] given without unit; taken as {'mol' if dim == 'amount' else 'kg'}.")
            add(name, v if dim == "amount" else v / db.get(name).MW)
        out.amount_known = True
    else:
        fluid_names = [n for n in (list(comp_dict) or (expected or [])) if n not in phi_names]
        if not comp_dict and phi and len(fluid_names) == 0:
            fluid_names = [n for n in comp_names if n not in phi_names]
        if comp_dict or (expected and len(expected) > 1):
            fr = _fill_fractions(given_key or "fractions", comp_dict, expected, log)
        elif fluid_names:
            fr = {fluid_names[0]: 1.0}
            given_key = given_key or "mole_fractions"
        else:
            fr = {}
        # convert fractions to moles on a basis
        if given_key == "mole_fractions":
            mol = {n: v for n, v in fr.items()}
            mass_f = {n: v * db.get(n).MW for n, v in mol.items()}
        elif given_key == "mass_fractions":
            mass_f = dict(fr)
        elif given_key == "volume_fractions":
            mass_f = {n: v * rho(n) for n, v in fr.items()}
        else:
            mass_f = {}
        if phi:
            phi_vals = {}
            for n, v in phi.items():
                v = float(v)
                if not 0.0 <= v <= 1.0:
                    raise InputError(f"phi['{n}'] = {v} is outside [0, 1].")
                phi_vals[n] = v
            phi_sum = sum(phi_vals.values())
            if phi_sum > 1.0 + FRACTION_TOL:
                raise InputError(f"Solid volume fractions phi sum to {phi_sum:.4f} > 1.")
            if phi_sum < 1.0 - 1e-12 and not mass_f:
                raise InputError("phi < 1 was given but no fluid component; list the fluid in 'components'.")
            # fluid: mass per m^3 of fluid from its composition (ideal volume of mixing)
            if mass_f:
                tot = sum(mass_f.values())
                w_f = {n: v / tot for n, v in mass_f.items()}
                v_f = sum(w / rho(n) for n, w in w_f.items())          # m^3 per kg fluid
                fluid_mass = (1.0 - phi_sum) / v_f                     # kg fluid per m^3 mixture
                mass_f = {n: w * fluid_mass for n, w in w_f.items()}
            for n, v in phi_vals.items():
                mass_f[n] = mass_f.get(n, 0.0) + v * rho(n)
            out.basis = "1 m^3 of mixture"
            log.note("phi interpreted as volume fraction of solid particles in the whole mixture "
                     "(ideal volume of mixing).")
        else:
            out.basis = {"mole_fractions": "1 mol of mixture", "mass_fractions": "1 kg of mixture",
                         "volume_fractions": "1 m^3 of (unmixed) components"}.get(given_key)
        for n, m in mass_f.items():
            add(n, m / db.get(n).MW)
        if given_key == "volume_fractions":
            log.assume("Volume fractions refer to pure-component volumes before mixing.")

    # merge the same component declared in two phases (e.g. water + ice)
    amounts: dict[str, float] = {}
    phase_spec: dict[str, str] = {}
    for (key, ph), mol in species.items():
        if key in amounts and phase_spec[key] != ph:
            log.note(f"{db.get(key).name} was given in more than one phase; treated as one component whose "
                     "phase split is computed from phase equilibrium at each temperature.")
            phase_spec[key] = "auto"
        else:
            phase_spec.setdefault(key, ph)
        amounts[key] = amounts.get(key, 0.0) + mol
    if sum(amounts.values()) <= 0:
        raise InputError("Total amount of all components is zero.")
    out.amounts, out.phase_spec = amounts, phase_spec

    # ---- total amount -------------------------------------------------------
    M = sum(amounts[k] * db.get(k).MW for k in amounts) / sum(amounts.values())
    for key, dim in (("total_mass", "mass"), ("total_moles", "amount")):
        if key not in d:
            continue
        v, _, assumed = to_si(d[key], dim, name=key)
        if v <= 0:
            raise InputError(f"{key} must be positive.")
        if assumed:
            log.assume(f"{key} given without unit; taken as {'kg' if dim == 'mass' else 'mol'}.")
        current = sum(amounts.values()) * (M if dim == "mass" else 1.0)
        if out.amount_known:
            if abs(current - v) / v > 5e-3:
                raise InputError(f"{key} = {v:.6g} SI conflicts with the sum of the component amounts "
                                 f"({current:.6g} SI).")
        else:
            scale = v / current
            out.amounts = {k: n * scale for k, n in amounts.items()}
            amounts = out.amounts
            out.amount_known = True
            out.basis = None
    if out.amount_known:
        out.known["m"] = sum(amounts[k] * db.get(k).MW for k in amounts)
        out.known["n"] = sum(amounts.values())


def _energies(d: dict, db: Database, out: NormalisedInput) -> None:
    """Convert Q, W, ΔH, ΔU, ΔS, Cv, ... to SI totals (or intensive symbols)."""
    log = out.log
    # Total amount without a composition
    if out.amounts is None:
        for key, dim, sym in (("total_mass", "mass", "m"), ("total_moles", "amount", "n")):
            if key in d:
                v, _, assumed = to_si(d[key], dim, name=key)
                if v <= 0:
                    raise InputError(f"{key} must be positive.")
                out.known[sym] = v
                if assumed:
                    log.assume(f"{key} given without unit; taken as SI.")

    M = None
    if out.amounts:
        M = sum(out.amounts[k] * db.get(k).MW for k in out.amounts) / sum(out.amounts.values())

    def make_total(sym: str, value: float, dim: str) -> float:
        """Turn a specific quantity (per kg / per mol) into a total for the system."""
        if dim in ("energy", "entropy"):
            return value
        per_mass = dim.endswith("per_mass")
        if per_mass and "m" in out.known:
            return value * out.known["m"]
        if not per_mass and "n" in out.known:
            return value * out.known["n"]
        # amount unknown: choose a calculation basis (Cp is intensive)
        if per_mass:
            if "n" in out.known and M:
                out.known["m"] = out.known["n"] * M
                return value * out.known["m"]
            _set_basis(out, db, mass=1.0)
            return value * out.known["m"]
        if "m" in out.known and M:
            out.known["n"] = out.known["m"] / M
            return value * out.known["n"]
        _set_basis(out, db, moles=1.0)
        return value * out.known["n"]

    for sym in ("Q", "W", "W_pv", "dH", "dU"):
        if sym in d:
            v, dim, assumed = to_si(d[sym], ["energy", "energy_per_mass", "energy_per_amount"], name=sym)
            if assumed:
                log.assume(f"{sym} given without unit; taken as J (total).")
            out.known[sym] = make_total(sym, v, dim)
            out.given[sym] = str(d[sym])
    if "dS" in d:
        v, dim, assumed = to_si(d["dS"], ["entropy", "entropy_per_mass", "entropy_per_amount"], name="dS")
        if assumed:
            log.assume("dS given without unit; taken as J/K (total).")
        out.known["dS"] = make_total("dS", v, dim)
        out.given["dS"] = str(d["dS"])
    if "dV" in d:
        v, _, assumed = to_si(d["dV"], "volume", name="dV")
        out.known["dV"] = v
    if "Cv" in d:
        v, dim, assumed = to_si(d["Cv"], ["entropy_per_mass", "entropy_per_amount", "entropy"], name="Cv")
        if assumed:
            log.assume("Cv given without unit; taken as J/(kg K).")
        if v <= 0:
            raise InputError("Cv must be positive.")
        sym = {"entropy_per_mass": "Cv_mass", "entropy_per_amount": "Cv_molar", "entropy": "Cv_total"}[dim]
        out.known[sym] = v
        out.given[sym] = str(d["Cv"])
    for key, sym, dims in (("gamma", "gamma", "dimensionless"), ("alpha", "alpha", "per_temperature"),
                           ("kappa_T", "kappa_T", "per_pressure"), ("density", "rho", "density")):
        if key in d:
            v, _, assumed = to_si(d[key], dims, name=key)
            if assumed and key != "gamma":
                log.assume(f"{key} given without unit; taken as SI.")
            out.known[sym] = v
            out.given[sym] = str(d[key])
    # default shaft work
    if "Q" in out.known and "W" not in out.known and "adiabatic" not in out.flags:
        out.known["W"] = 0.0
        log.assume("No work W given: shaft/stirring work taken as zero (W = 0).")
    for flag in out.flags:
        out.known[flag] = 1.0


def _set_basis(out: NormalisedInput, db: Database, mass: float | None = None, moles: float | None = None) -> None:
    """Fix the system amount to a calculation basis when only intensive data were given."""
    if out.amounts:
        M = sum(out.amounts[k] * db.get(k).MW for k in out.amounts) / sum(out.amounts.values())
        n_tot = moles if moles is not None else mass / M
        scale = n_tot / sum(out.amounts.values())
        out.amounts = {k: v * scale for k, v in out.amounts.items()}
        out.known["n"], out.known["m"] = n_tot, n_tot * M
    elif mass is not None:
        out.known["m"] = mass
    else:
        out.known["n"] = moles
    out.amount_known = True
    out.basis = f"{mass} kg" if mass is not None else f"{moles} mol"
    out.log.assume(f"Energy given per unit amount and no total amount: calculation basis of {out.basis} "
                   "(Cp is intensive, so the result does not depend on this choice).")


def describe_inputs() -> str:
    """Short help text listing every accepted input key."""
    return __doc__.split("Accepted keys")[1].split("Sign convention")[0]


__all__ = ["SystemInput", "NormalisedInput", "normalise", "describe_inputs", "UnitError"]
