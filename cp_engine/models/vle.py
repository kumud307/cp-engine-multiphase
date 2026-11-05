"""Vapour–liquid equilibrium of liquid mixtures (γ–φ approach at low pressure).

Equal fugacities of every component in the two phases,  f̂_i^V = f̂_i^L, with
    vapour:  f̂_i^V = y_i φ̂_i P
    liquid:  f̂_i^L = x_i γ_i f_i  ≈  x_i γ_i φ_i^sat P_i^sat   (Lewis–Randall reference, Poynting factor ≈ 1)
gives  y_i P Φ_i = x_i γ_i P_i^sat,  Φ_i = φ̂_i / φ_i^sat, i.e. K_i = y_i/x_i = γ_i P_i^sat / (Φ_i P).
Vapour models:
    "eos"           φ̂_i and φ_i^sat from the SAME equation of state (Peng–Robinson by default) that gives the
                    departure enthalpy of the vapour, so phase equilibrium and enthalpy use one vapour model
                    (engine default);
    "ideal"         φ̂_i = φ_i^sat = 1: modified Raoult's law;
    "lewis_randall" φ̂_i = φ_i from the pure-component virial B_ii.  γ_i comes from the same NRTL model as the excess enthalpy; P_i^sat from the
Antoine equations of the database.  Components without vapour-pressure or gas-phase data (salts,
glycerol, heavy alkanes) are non-volatile: K_i = 0.

* bubble point: Σ_i K_i(T, P, x) x_i = 1 at the liquid composition x = z;
* dew point:    Σ_i y_i / K_i(T, P, x) = 1 with x_i ∝ y_i / K_i;
* isothermal flash (Rachford–Rice):  Σ_i z_i (K_i − 1) / [1 + V (K_i − 1)] = 0  for the vapour fraction V,
  x_i = z_i / [1 + V (K_i − 1)],  y_i = K_i x_i, iterated until γ(x) is converged.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

from ..errors import DataError
from .eos import R, VirialEOS


class VLE:
    """Vapour–liquid equilibrium of one liquid phase (components `keys`) at pressure P."""

    def __init__(self, keys: list[str], db, nrtl=None, log=None, vapour: str = "ideal", eos_kind: str = "PR"):
        self.keys, self.db, self.nrtl, self.log = list(keys), db, nrtl, log
        self.comps = [db.get(k) for k in self.keys]
        self.volatile = np.array([bool(c.raw.get("antoine")) and c.has_phase("gas") for c in self.comps])
        if not self.volatile.any():
            raise DataError("No volatile component (no vapour-pressure data).")
        if vapour not in ("ideal", "lewis_randall", "eos"):
            raise ValueError("vapour must be 'ideal', 'lewis_randall' or 'eos'")
        self._virial = None
        self._eos_mix = self._eos_pure = None
        self._vol_idx = [i for i, v in enumerate(self.volatile) if v]
        if vapour == "eos":
            vol_keys = [self.keys[i] for i in self._vol_idx]
            try:
                self._eos_mix = db.eos(vol_keys, eos_kind)
                self._eos_pure = {i: db.eos([self.keys[i]], eos_kind) for i in self._vol_idx}
            except (DataError, ValueError):
                vapour = "ideal"                  # no critical constants: fall back, and say so
                if log is not None:
                    log.assume("Vapour treated as an ideal gas in the phase equilibrium (no critical constants).")
        if vapour == "lewis_randall":
            self._virial = [VirialEOS([c.critical]) if (v and c.critical) else None
                            for c, v in zip(self.comps, self.volatile)]
        self.vapour = vapour
        self._ideal_only = False                  # True while bracketing with the cheap ideal-vapour K
        self._sat_cache: dict[float, np.ndarray] = {}

    # ------------------------------------------------------------------ K-values
    def psat(self, T: float) -> np.ndarray:
        return np.array([c.psat(T) if v else 0.0 for c, v in zip(self.comps, self.volatile)])

    def _check_range(self, T: float) -> None:
        """Warn (once per component) if T is outside the Antoine ranges of a volatile component."""
        if self.log is None:
            return
        for c, v in zip(self.comps, self.volatile):
            if v and not any(a["Tmin"] - 1.0 <= T <= a["Tmax"] + 1.0 for a in c.raw["antoine"]):
                self.log.warn(f"Vapour pressure of {c.name} extrapolated to {T:.1f} K (outside its Antoine ranges).")

    def gamma(self, x, T: float) -> np.ndarray:
        if self.nrtl is None:
            return np.ones(len(self.keys))
        return np.asarray(self.nrtl.gamma(list(x), T), float)

    @staticmethod
    def _ln_phi_vapour(eos, T: float, P: float, y) -> np.ndarray | None:
        """ln φ̂ of the vapour root, or None if the EOS has no vapour-like root at (T, P, y)."""
        try:
            if eos.Z(T, P, y, "vapour") < 0.3:     # only a liquid-like root exists here
                return None
            return np.asarray(eos.ln_phi(T, P, y, "vapour"), float)
        except (ValueError, ZeroDivisionError, OverflowError):
            return None

    def K(self, T: float, P: float, x, y=None) -> np.ndarray:
        x = np.asarray(x, float)
        ps = self.psat(T)
        K = self.gamma(x, T) * ps / P
        if self._virial is not None:          # Φ_i = exp[B_ii (P − P_i^sat)/RT]  (Lewis–Randall, no Poynting)
            for i, v in enumerate(self._virial):
                if v is not None:
                    K[i] /= math.exp(v.B_mix(T, [1.0])[0] * (P - ps[i]) / (R * T))
        if self._eos_mix is not None and not self._ideal_only:   # Φ_i = φ̂_i(T,P,y)/φ_i^sat, same EOS as the enthalpy
            base = K.copy()
            sat = self._sat_cache.get(T)
            if sat is None:
                sat = np.zeros(len(K))
                for i in self._vol_idx:
                    l = self._ln_phi_vapour(self._eos_pure[i], T, ps[i], [1.0])
                    sat[i] = 0.0 if l is None else l[0]
                if len(self._sat_cache) > 2000:
                    self._sat_cache.clear()
                self._sat_cache[T] = sat
            yy = np.asarray(y, float) if y is not None else np.where(self.volatile, base * x, 0.0)
            for _ in range(4):                # φ̂ depends weakly on y: a few fixed-point passes suffice
                yv = yy[self._vol_idx]
                if yv.sum() <= 0:
                    break
                lv = self._ln_phi_vapour(self._eos_mix, T, P, list(yv / yv.sum()))
                if lv is None:
                    K = base
                    break
                corr = np.zeros(len(K))
                corr[self._vol_idx] = sat[self._vol_idx] - lv
                K = base * np.exp(corr)
                if y is not None:
                    break
                yy = np.where(self.volatile, K * x, 0.0)
        return np.where(self.volatile, K, 0.0)

    # ------------------------------------------------------------------ bubble / dew points
    def _bracket(self, f, P: float) -> tuple[float, float]:
        """Bracket a root of f(T). With an EOS vapour, the root is first located with the ideal-vapour K
        (cheap) and the bracket is then confirmed / widened with the full model."""
        if self._eos_mix is not None and not self._ideal_only:
            self._ideal_only = True
            try:
                a, b = self._bracket(f, P)
            except ValueError:
                a, b = None, None
            finally:
                self._ideal_only = False
            if a is not None:
                for d in (1.0, 3.0, 10.0):
                    lo, hi = a - d, b + d
                    try:
                        if f(lo) * f(hi) <= 0:
                            return lo, hi
                    except (ValueError, DataError, ZeroDivisionError, OverflowError):
                        break
        tb = [c.Tb_at(P) for c, v in zip(self.comps, self.volatile) if v and c.Tb_at(P)]
        lo, hi = (min(tb) - 40.0, max(tb) + 40.0) if tb else (200.0, 700.0)
        Ts = np.concatenate([np.arange(lo, hi, 2.0), np.arange(200.0, 700.0, 5.0)])
        vals = []
        for T in Ts:
            try:
                vals.append(f(T))
            except (ValueError, DataError, ZeroDivisionError, OverflowError):
                vals.append(np.nan)
        vals = np.array(vals)
        for i in range(len(Ts) - 1):
            if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] <= 0:
                return Ts[i], Ts[i + 1]
        raise ValueError("no bubble/dew point between 200 and 700 K")

    def bubble_T(self, P: float, x) -> float:
        x = np.asarray(x, float)
        f = lambda T: float(np.sum(self.K(T, P, x) * x)) - 1.0
        T = brentq(f, *self._bracket(f, P), xtol=1e-10)
        self._check_range(T)
        return T

    def bubble_P(self, T: float, x) -> float:
        """Bubble pressure: Σ K_i(T, P, x) x_i = 1 (exact in one step for an ideal vapour, K ∝ 1/P)."""
        x = np.asarray(x, float)
        P = float(np.sum(self.K(T, 1.0, x) * x))
        if self._eos_mix is None and self._virial is None:
            return P
        for _ in range(100):
            S = float(np.sum(self.K(T, P, x) * x))
            P *= S
            if abs(S - 1.0) < 1e-12:
                break
        return P

    def _dew_x(self, T: float, P: float, y) -> np.ndarray:
        x = np.asarray(y, float).copy()
        for _ in range(200):
            K = self.K(T, P, x, y)
            new = np.divide(y, K, out=np.zeros_like(x), where=K > 0)
            s = new.sum()
            new = new / s if s > 0 else x
            if np.max(np.abs(new - x)) < 1e-12:
                return new
            x = 0.5 * (x + new)
        return x

    def dew_T(self, P: float, y) -> float:
        y = np.asarray(y, float)
        if np.any(y[~self.volatile] > 0):
            return math.inf                                 # a non-volatile component never fully evaporates
        def f(T):
            x = self._dew_x(T, P, y)
            K = self.K(T, P, x, y)
            return 1.0 - float(np.sum(np.divide(y, K, out=np.zeros_like(y), where=K > 0)))
        T = brentq(f, *self._bracket(f, P), xtol=1e-10)
        self._check_range(T)
        return T

    # ------------------------------------------------------------------ isothermal flash
    def flash(self, T: float, P: float, z) -> tuple[float, np.ndarray, np.ndarray]:
        """(vapour fraction V, liquid x, vapour y) at T, P for overall composition z."""
        z = np.asarray(z, float) / np.sum(z)
        x = z.copy()
        y = None
        V = 0.5
        for _ in range(500):
            K = self.K(T, P, x, y)
            g = lambda v: float(np.sum(z * (K - 1.0) / (1.0 + v * (K - 1.0))))
            if g(0.0) <= 0.0:
                V = 0.0                                     # sub-cooled liquid
            elif g(1.0) >= 0.0:
                V = 1.0                                     # superheated vapour
            else:
                V = brentq(g, 0.0, 1.0, xtol=1e-14)
            x_new = z / (1.0 + V * (K - 1.0))
            x_new /= x_new.sum()
            y_new = K * x_new
            y_new = y_new / y_new.sum() if y_new.sum() > 0 else y_new
            done = np.max(np.abs(x_new - x)) < 1e-13 and (y is None or np.max(np.abs(y_new - y)) < 1e-13)
            x, y = x_new, y_new
            if done:
                break
        K = self.K(T, P, x, y)
        y = K * x
        y = y / y.sum() if y.sum() > 0 else y
        return V, x, y
