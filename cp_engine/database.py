"""Component and binary-interaction database.

The data live in two JSON files (data/components.json, data/binary_nrtl.json).
This module only *reads* them and converts every value to SI, so adding a new
component never requires changing code: add an entry to the JSON file.

Internal SI conventions
-----------------------
    MW            kg/mol
    Cp            J/(mol K)    (molar, per phase)
    dHfus, dHvap  J/mol
    density       kg/m^3
    alpha         1/K          (volumetric thermal expansion coefficient)
    kappa_T       1/Pa         (isothermal compressibility)
"""
from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

from scipy import integrate

from .errors import DataError, RangeError
from .report import EXTRAPOLATION_LIMIT_K, Log, null_log
from .units import convert, parse_unit

R = 8.314462618  # J/(mol K), CODATA 2018 molar gas constant
P_ATM = 101325.0  # Pa

PHASES = ("solid", "liquid", "gas")


def default_data_dir() -> Path:
    """Locate the data/ folder: $CP_ENGINE_DATA, else <project>/data."""
    env = os.environ.get("CP_ENGINE_DATA")
    if env:
        return Path(env)
    return Path(__file__).resolve().parent.parent / "data"


# ---------------------------------------------------------------------------
# Correlation evaluation
# ---------------------------------------------------------------------------
def check_range(T: float, Tmin: float | None, Tmax: float | None, what: str,
                log: Log, strict: bool = True) -> None:
    """Warn if T is outside [Tmin, Tmax]; raise RangeError if far outside.

    'Far' means more than EXTRAPOLATION_LIMIT_K kelvin beyond the range.
    """
    if Tmin is None or Tmax is None:
        return
    if Tmin - 0.05 <= T <= Tmax + 0.05:      # 0.05 K tolerance (e.g. 273.15 vs 273.16 K)
        return
    dist = (Tmin - T) if T < Tmin else (T - Tmax)
    msg = (f"{what}: T = {T:.2f} K is outside the correlation range "
           f"[{Tmin:g}, {Tmax:g}] K (extrapolating by {dist:.1f} K).")
    if strict and dist > EXTRAPOLATION_LIMIT_K:
        raise RangeError(msg + f" This is more than {EXTRAPOLATION_LIMIT_K:g} K outside; "
                         "refusing to extrapolate. Supply data valid at this temperature.")
    log.warn(msg)


def _molar_factor(unit: str, MW: float) -> float:
    """Factor converting a heat-capacity unit to J/(mol K) (uses MW for mass units)."""
    pu = parse_unit(unit)
    if pu.dims == {"E": 1, "N": -1, "T": -1}:
        return convert(1.0, unit, "J/(mol K)")
    if pu.dims == {"E": 1, "M": -1, "T": -1}:
        return convert(1.0, unit, "J/(kg K)") * MW
    raise DataError(f"Heat-capacity unit '{unit}' is neither per mole nor per mass.")


@dataclass
class CpCorrelation:
    """Heat capacity of one phase of one component as a function of T."""
    type: str
    coeffs: list[float]
    factor: float              # converts the correlation's unit to J/(mol K)
    Tmin: float | None
    Tmax: float | None
    source: str
    verify: bool

    def raw(self, T: float) -> float:
        """Evaluate the correlation (J/(mol K)) without any range check."""
        c = self.coeffs
        if self.type in ("dippr100", "poly"):
            val = sum(ck * T ** k for k, ck in enumerate(c))
        elif self.type == "constant":
            val = c[0]
        elif self.type == "shomate":
            t = T / 1000.0
            val = c[0] + c[1] * t + c[2] * t ** 2 + c[3] * t ** 3 + c[4] / t ** 2
        elif self.type == "dippr107":
            x, y = c[2] / T, c[4] / T
            val = c[0] + c[1] * (x / math.sinh(x)) ** 2 + c[3] * (y / math.cosh(y)) ** 2
        else:
            raise DataError(f"Unknown Cp correlation type '{self.type}'.")
        return val * self.factor

    def __call__(self, T: float, log: Log | None = None, what: str = "Cp") -> float:
        check_range(T, self.Tmin, self.Tmax, what, log or null_log())
        return self.raw(T)

    def integral(self, T1: float, T2: float) -> float:
        """∫ Cp dT from T1 to T2 in J/mol (used for enthalpies)."""
        if T1 == T2:
            return 0.0
        if self.type == "constant":
            return self.raw(T1) * (T2 - T1)
        val, _ = integrate.quad(self.raw, T1, T2, epsabs=1e-9, epsrel=1e-11, limit=100)
        return val


@dataclass
class DensityModel:
    """Density of one phase of one component (kg/m^3)."""
    type: str
    coeffs: list[float]
    MW: float
    Tmin: float | None
    Tmax: float | None
    source: str
    verify: bool
    alpha_fixed: float | None = None    # stored alpha (for constant density)

    def __call__(self, T: float, P: float = P_ATM, log: Log | None = None, what: str = "density") -> float:
        log = log or null_log()
        c = self.coeffs
        if self.type == "constant":
            return c[0]
        check_range(T, self.Tmin, self.Tmax, what, log)
        if self.type == "poly":
            return sum(ck * T ** k for k, ck in enumerate(c))
        if self.type == "dippr105":
            molar = c[0] / c[1] ** (1.0 + (1.0 - T / c[2]) ** c[3])   # mol/dm^3
            return molar * 1000.0 * self.MW
        if self.type == "kell_water":
            t = T - 273.15
            num = (999.83952 + 16.945176 * t - 7.9870401e-3 * t ** 2 - 46.170461e-6 * t ** 3
                   + 105.56302e-9 * t ** 4 - 280.54253e-12 * t ** 5)
            return num / (1.0 + 16.879850e-3 * t)
        if self.type == "ideal_gas":
            return P * self.MW / (R * T)
        raise DataError(f"Unknown density model '{self.type}'.")

    def alpha(self, T: float, P: float = P_ATM, log: Log | None = None) -> float | None:
        """Volumetric thermal expansion coefficient α = −(1/ρ)(∂ρ/∂T)_P  [1/K].

        For T-dependent density models α is obtained by a central difference
        of the density model itself (so ρ and α are always consistent).
        """
        if self.type == "constant":
            return self.alpha_fixed
        if self.type == "ideal_gas":
            return 1.0 / T
        h = 0.01
        rp, rm = self(T + h, P), self(T - h, P)     # range already checked by the density call itself
        return -(rp - rm) / (2 * h) / self(T, P)


@dataclass
class Component:
    """All pure-component data of one substance, in SI units."""
    key: str
    name: str
    formula: str
    MW: float                               # kg/mol
    miscibility_class: str
    Tm: float | None = None                 # K
    dHfus: float | None = None              # J/mol
    Tb: float | None = None                 # K (normal boiling point)
    dHvap: float | None = None              # J/mol
    cp_corr: dict[str, CpCorrelation] = field(default_factory=dict)
    density_models: dict[str, DensityModel] = field(default_factory=dict)
    kappa: dict[str, Any] = field(default_factory=dict)     # value or {'type': 'ln_poly', ...}
    raw: dict[str, Any] = field(default_factory=dict, repr=False)

    # ------------------------------------------------------------------
    def has_phase(self, phase: str) -> bool:
        return phase in self.cp_corr

    @property
    def critical(self):
        """Critical constants (T_c, P_c, V_c, ω) for the equations of state, or None."""
        from .models.eos import Critical
        e = self.raw.get("critical")
        if not e:
            return None
        Vc = e.get("Vc")
        return Critical(self.name, float(e["Tc"]), float(e["Pc"]), float(e["omega"]),
                        None if Vc is None else float(Vc) * 1e-6)

    def cp(self, phase: str, T: float, log: Log | None = None) -> float:
        """Molar isobaric heat capacity of the pure phase, J/(mol K)."""
        if phase not in self.cp_corr:
            have = ", ".join(self.cp_corr) or "none"
            raise DataError(f"No {phase}-phase Cp correlation for {self.name} "
                            f"(available phases: {have}). Add one to data/components.json.")
        return self.cp_corr[phase](T, log, what=f"Cp of {self.name} ({phase})")

    def cp_mass(self, phase: str, T: float, log: Log | None = None) -> float:
        """Specific (mass) isobaric heat capacity, J/(kg K)."""
        return self.cp(phase, T, log) / self.MW

    def density(self, phase: str, T: float, P: float = P_ATM, log: Log | None = None) -> float:
        if phase not in self.density_models:
            raise DataError(f"No {phase}-phase density for {self.name}.")
        return self.density_models[phase](T, P, log, what=f"density of {self.name} ({phase})")

    def alpha(self, phase: str, T: float, P: float = P_ATM, log: Log | None = None) -> float:
        dm = self.density_models.get(phase)
        a = dm.alpha(T, P, log) if dm else None
        if a is None:
            raise DataError(f"No thermal expansion coefficient for {self.name} ({phase}).")
        return a

    def kappa_T(self, phase: str, T: float, P: float = P_ATM) -> float:
        """Isothermal compressibility at 1 atm, 1/Pa (constant, or ln κ = Σ c_k T^k)."""
        if phase == "gas":
            return 1.0 / P                      # ideal gas
        if phase not in self.kappa:
            raise DataError(f"No isothermal compressibility for {self.name} ({phase}).")
        k = self.kappa[phase]
        if isinstance(k, dict):                 # {"type": "ln_poly", "coeffs": [...]}
            return math.exp(sum(c * T ** i for i, c in enumerate(k["coeffs"])))
        return k

    # ------------------------------------------------------------------
    # Pressure dependence (all Cp correlations are for ~1 atm)
    # ------------------------------------------------------------------
    TAIT_C = 0.0894     # universal Tait constant (dimensionless)

    def molar_volume(self, phase: str, T: float, P: float = P_ATM) -> float:
        """Molar volume v(T, P) in m^3/mol from the Tait equation

            v(T,P) = v0(T) [1 − C ln(1 + (P − P0)/(B(T) + P0))],   B(T) = C/κ_T(T) − P0

        so that −(1/v)(∂v/∂P)_T = κ_T at P0. v0(T) comes from the density model at 1 atm; for a
        constant-density phase v0(T) = v_ref exp[α (T − 298.15 K)].
        """
        dm = self.density_models[phase]
        if dm.type == "constant":
            a = dm.alpha_fixed or 0.0
            v0 = self.MW / dm.coeffs[0] * math.exp(a * (T - 298.15))
        else:
            v0 = self.MW / dm(T, P_ATM)
        if abs(P - P_ATM) < 1.0:
            return v0
        kT = self.kappa_T(phase, T)
        B = self.TAIT_C / kT - P_ATM
        return v0 * (1.0 - self.TAIT_C * math.log(1.0 + (P - P_ATM) / (B + P_ATM)))

    def _p_integral(self, f, P: float, n: int = 8) -> float:
        """∫_{P0}^{P} f(p) dp by n-point Gauss–Legendre quadrature."""
        import numpy as np
        x, w = np.polynomial.legendre.leggauss(n)
        mid, half = 0.5 * (P + P_ATM), 0.5 * (P - P_ATM)
        return half * sum(wi * f(mid + half * xi) for xi, wi in zip(x, w))

    def cp_pressure_correction(self, phase: str, T: float, P: float) -> float:
        """Cp(T,P) − Cp(T,1 atm) in J/(mol K), from the exact relation

            (∂Cp/∂P)_T = −T (∂²v/∂T²)_P      [ = −T v (α² + (∂α/∂T)_P) ]

        with v(T, P) from the Tait equation (molar_volume). Zero for an ideal gas.
        """
        if phase == "gas" or abs(P - P_ATM) < 1.0 or phase not in self.density_models:
            return 0.0
        h = 0.5
        d2v = lambda p: (self.molar_volume(phase, T + h, p) - 2 * self.molar_volume(phase, T, p)
                         + self.molar_volume(phase, T - h, p)) / h ** 2
        return -T * self._p_integral(d2v, P)

    def h_pressure_correction(self, phase: str, T: float, P: float) -> float:
        """H(T,P) − H(T,1 atm) in J/mol, from (∂H/∂P)_T = v − T (∂v/∂T)_P.

        Differentiating this with respect to T gives exactly −T∫(∂²v/∂T²)dP, i.e.
        cp_pressure_correction(), so point and interval-mean Cp stay consistent at any P.
        """
        if phase == "gas" or abs(P - P_ATM) < 1.0 or phase not in self.density_models:
            return 0.0
        h = 0.5
        f = lambda p: self.molar_volume(phase, T, p) - T * (self.molar_volume(phase, T + h, p)
                                                            - self.molar_volume(phase, T - h, p)) / (2 * h)
        return self._p_integral(f, P)

    def psat(self, T: float, log: Log | None = None) -> float:
        """Vapour pressure in Pa from the Antoine equation log10(P/bar) = A − B/(T + C) (NIST sets).

        The set whose range contains T is used; otherwise the nearest set, with a warning.
        """
        sets = self.raw.get("antoine")
        if not sets:
            raise DataError(f"No vapour-pressure (Antoine) data for {self.name}.")
        inside = [a for a in sets if a["Tmin"] <= T <= a["Tmax"]]
        a = inside[0] if inside else min(sets, key=lambda a: min(abs(T - a["Tmin"]), abs(T - a["Tmax"])))
        if not inside and log is not None:
            log.warn(f"Vapour pressure of {self.name}: T = {T:.1f} K outside the Antoine ranges; nearest set extrapolated.")
        return 10.0 ** (a["A"] - a["B"] / (T + a["C"])) * 1e5

    def Tb_at(self, P: float) -> float | None:
        """Boiling point at pressure P: from the Antoine equation if available, otherwise from the
        Clausius–Clapeyron equation

        1/Tb(P) = 1/Tb(1 atm) − (R/ΔHvap) ln(P / 1 atm)
        (assumes constant ΔHvap and ideal vapour).
        """
        if self.Tb is None:
            return None
        if abs(P - P_ATM) < 1.0:
            return self.Tb
        if self.raw.get("antoine"):                 # solve P_sat(T) = P with the Antoine equation
            from scipy.optimize import brentq
            try:
                return brentq(lambda T: self.psat(T) - P, 150.0, 900.0)
            except (ValueError, OverflowError, ZeroDivisionError):
                pass
        if self.dHvap is None:
            return self.Tb
        return 1.0 / (1.0 / self.Tb - R / self.dHvap * math.log(P / P_ATM))


# ---------------------------------------------------------------------------
# Database
# ---------------------------------------------------------------------------
def _val(entry: dict | None, si_unit: str) -> float | None:
    if entry is None or entry.get("value") is None:
        return None
    return convert(float(entry["value"]), entry.get("unit", si_unit), si_unit)


class Database:
    """Loads and serves component and binary data (SI units)."""

    def __init__(self, data_dir: str | Path | None = None):
        self.data_dir = Path(data_dir) if data_dir else default_data_dir()
        comp_file = self.data_dir / "components.json"
        bin_file = self.data_dir / "binary_nrtl.json"
        if not comp_file.exists():
            raise DataError(f"Component database not found at {comp_file}. "
                            "Set CP_ENGINE_DATA or pass data_dir=.")
        with open(comp_file, encoding="utf-8") as f:
            self._comp_json = json.load(f)["components"]
        self._bin_json: dict[str, Any] = {"nrtl": [], "miscibility": [], "solubility": []}
        if bin_file.exists():
            with open(bin_file, encoding="utf-8") as f:
                self._bin_json.update(json.load(f))
        self.components: dict[str, Component] = {k: self._build(k, v) for k, v in self._comp_json.items()}
        # alias table: lower-case name -> (key, phase or None)
        self._alias: dict[str, tuple[str, str | None]] = {}
        for k, v in self._comp_json.items():
            for a in [k, v.get("name", k), v.get("formula", "")] + v.get("aliases", []):
                if a:
                    self._alias.setdefault(a.lower(), (k, None))
            for a, ph in v.get("phase_aliases", {}).items():
                self._alias[a.lower()] = (k, ph)

    # ------------------------------------------------------------------
    @staticmethod
    def _build(key: str, d: dict[str, Any]) -> Component:
        try:
            MW = _val(d["MW"], "kg/mol")
        except KeyError:
            raise DataError(f"Component '{key}' has no molecular weight (MW).") from None
        comp = Component(
            key=key, name=d.get("name", key), formula=d.get("formula", ""), MW=MW,
            miscibility_class=d.get("miscibility_class", "unknown"),
            Tm=_val(d.get("Tm"), "K"), dHfus=_val(d.get("dHfus"), "J/mol"),
            Tb=_val(d.get("Tb"), "K"), dHvap=_val(d.get("dHvap"), "J/mol"), raw=d,
        )
        for phase, c in d.get("cp", {}).items():
            coeffs = c.get("coeffs", [c.get("value")])
            comp.cp_corr[phase] = CpCorrelation(
                type=c["type"], coeffs=[float(x) for x in coeffs],
                factor=_molar_factor(c.get("unit", "J/(mol K)"), MW),
                Tmin=c.get("Tmin"), Tmax=c.get("Tmax"),
                source=c.get("source", ""), verify=bool(c.get("verify", False)))
        alphas = d.get("alpha", {})
        for phase, r in d.get("density", {}).items():
            if r["type"] == "constant":
                coeffs = [convert(float(r["value"]), r.get("unit", "kg/m3"), "kg/m3")]
            else:
                coeffs = [float(x) for x in r.get("coeffs", [])]
            comp.density_models[phase] = DensityModel(
                type=r["type"], coeffs=coeffs, MW=MW, Tmin=r.get("Tmin"), Tmax=r.get("Tmax"),
                source=r.get("source", ""), verify=bool(r.get("verify", False)),
                alpha_fixed=_val(alphas.get(phase), "1/K"))
        for phase, r in d.get("kappa_T", {}).items():
            comp.kappa[phase] = dict(r) if r.get("type") == "ln_poly" else _val(r, "1/Pa")
        return comp

    # ------------------------------------------------------------------
    def resolve(self, name: str) -> tuple[Component, str | None]:
        """Find a component by key, name, formula or alias (case-insensitive).

        Returns (component, implied_phase); implied_phase is e.g. 'solid' for 'ice'.
        """
        hit = self._alias.get(str(name).strip().lower())
        if hit is None:
            close = [k for k in self.components if str(name).lower()[:3] in k.lower()]
            hint = f" Did you mean: {', '.join(close)}?" if close else ""
            raise DataError(f"Unknown component '{name}'.{hint} Available: "
                            f"{', '.join(self.components)}. Add new components to data/components.json.")
        key, phase = hit
        return self.components[key], phase

    def get(self, name: str) -> Component:
        return self.resolve(name)[0]

    def names(self) -> list[str]:
        return list(self.components)

    # ------------------------------------------------------------------
    def nrtl_pair(self, i: str, j: str, prefer: str = "HE") -> dict[str, Any] | None:
        """NRTL parameters for the pair (i, j), oriented so that '_ij' refers to i→j.

        A pair can have several parameter sets (e.g. regressed to VLE data or to excess-enthalpy
        data). `prefer` selects by the "fitted_to" field: H^E-fitted sets describe the temperature
        dependence of G^E – and hence H^E and Cp^E – far better, so they are the default for Cp.
        If the preferred set does not exist, any available set is returned.
        """
        found = []
        for p in self._bin_json.get("nrtl", []):
            if p["i"] == i and p["j"] == j:
                found.append(dict(p))
            elif p["i"] == j and p["j"] == i:
                q = dict(p)
                q["i"], q["j"] = i, j
                for s in ("a", "b", "e", "f"):
                    q[f"{s}_ij"], q[f"{s}_ji"] = p.get(f"{s}_ji", 0.0), p.get(f"{s}_ij", 0.0)
                found.append(q)
        if not found:
            return None
        # "HE" also matches "HE+SLE" (fitted jointly to excess enthalpies and freezing points)
        return next((p for p in found if str(p.get("fitted_to", "")).startswith(prefer)), found[0])

    def miscibility(self, i: str, j: str, T: float | None = None, x_i: float | None = None) -> tuple[bool, str]:
        """Are liquids i and j miscible? Returns (miscible, reason).

        If the database holds a measured liquid–liquid coexistence curve (binodal) for the pair, the
        answer depends on temperature and composition: two liquid phases form below the curve, one
        above it (upper critical solution temperature). Without T, the pair is reported as it behaves
        at room temperature (the 'miscible' flag).
        """
        for p in self._bin_json.get("miscibility", []):
            if {p["i"], p["j"]} == {i, j}:
                if "binodal" in p and T is not None:
                    import numpy as np
                    xs, Ts = zip(*p["binodal"])
                    xa = x_i if (x_i is not None and p.get("binodal_x_of", p["i"]) == i) else (
                        None if x_i is None else 1.0 - x_i)
                    if xa is None:
                        T_split = float(p.get("ucst_K", max(Ts)))
                    elif xs[0] <= xa <= xs[-1]:
                        T_split = float(np.interp(xa, xs, Ts))
                    else:                        # outside the measured curve: dilute, single phase
                        return True, f"database binodal: x = {xa:.3f} outside the two-phase region ({p.get('source', '')[:60]}...)"
                    if T > T_split:
                        return True, f"database binodal: T = {T:.1f} K above the coexistence curve ({T_split:.1f} K)"
                    return False, f"database binodal: T = {T:.1f} K below the coexistence curve ({T_split:.1f} K)"
                return bool(p["miscible"]), f"database: {p.get('source', '')}"
        ci, cj = self.components[i].miscibility_class, self.components[j].miscibility_class
        classes = {ci, cj}
        if "nonpolar" in classes and classes & {"aqueous", "polar"}:
            return False, f"heuristic: {ci} + {cj} assumed immiscible"
        return True, f"heuristic: {ci} + {cj} assumed miscible"

    def solubility_entry(self, solute: str, solvent: str) -> dict[str, Any] | None:
        for p in self._bin_json.get("solubility", []):
            if p["solute"] == solute and p["solvent"] == solvent:
                return p
        return None

    def apparent_cp_entry(self, solute: str, solvent: str) -> dict[str, Any] | None:
        """Apparent-molar-heat-capacity model φCp(m) of an electrolyte solute in a solvent, if any."""
        for p in self._bin_json.get("apparent_molar_cp", []):
            if p["solute"] == solute and p["solvent"] == solvent:
                return p
        return None

    def verify_items(self) -> list[str]:
        """List every database value flagged "verify": true (for DATA_TO_VERIFY.md)."""
        out = []
        for k, d in self._comp_json.items():
            for prop in ("MW", "Tm", "dHfus", "Tb", "dHvap"):
                e = d.get(prop)
                if isinstance(e, dict) and e.get("verify"):
                    out.append(f"{k}.{prop} = {e.get('value')} {e.get('unit', '')}  ({e.get('source', '')})")
            crit = d.get("critical")
            if isinstance(crit, dict) and crit.get("verify"):
                out.append(f"{k}.critical = Tc {crit.get('Tc')} K, Pc {crit.get('Pc')} Pa  ({crit.get('source', '')})")
            for group in ("cp", "density", "alpha", "kappa_T"):
                for ph, e in d.get(group, {}).items():
                    if isinstance(e, dict) and e.get("verify"):
                        v = e.get("coeffs", e.get("value"))
                        out.append(f"{k}.{group}[{ph}] = {v} {e.get('unit', '')}  ({e.get('source', '')})")
        for sect in ("nrtl", "miscibility", "solubility"):
            for p in self._bin_json.get(sect, []):
                if p.get("verify"):
                    who = f"{p.get('i', p.get('solute'))}-{p.get('j', p.get('solvent'))}"
                    out.append(f"binary.{sect}[{who}]  ({p.get('source', '')})")
                if sect == "solubility" and p.get("dH_sol", {}).get("verify"):
                    out.append(f"binary.solubility[{p['solute']}-{p['solvent']}].dH_sol = "
                               f"{p['dH_sol']['value']} {p['dH_sol']['unit']}")
        return out


    def eos(self, keys: list[str], kind: str = "PR"):
        """Equation of state (PR, SRK, RK, vdW or virial) for a pure gas or a gas mixture of `keys`."""
        from .models.eos import make_eos
        crit = [self.get(k).critical for k in keys]
        missing = [self.get(k).name for k, c in zip(keys, crit) if c is None]
        if missing:
            raise DataError(f"No critical constants for {', '.join(missing)}; real-gas correction unavailable.")
        return make_eos(crit, kind)


@lru_cache(maxsize=4)
def load_database(data_dir: str | None = None) -> Database:
    """Cached loader so the JSON is parsed only once per session."""
    return Database(data_dir)
