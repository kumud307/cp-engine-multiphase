"""Vapour–liquid equilibrium of liquid mixtures (γ–φ approach at low pressure).

Equal fugacities of every component in the two phases,  f̂_i^V = f̂_i^L, with
    vapour:  f̂_i^V = y_i φ̂_i P            (φ̂_i = 1: ideal gas, or Lewis–Randall with the virial B_ii)
    liquid:  f̂_i^L = x_i γ_i f_i  ≈  x_i γ_i φ_i^sat P_i^sat   (Lewis–Randall reference, Poynting factor ≈ 1)
gives the modified Raoult's law  y_i P Φ_i = x_i γ_i P_i^sat,  Φ_i = φ̂_i / φ_i^sat, i.e. K_i = y_i/x_i =
γ_i P_i^sat / (Φ_i P).  γ_i comes from the same NRTL model as the excess enthalpy; P_i^sat from the
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

    def __init__(self, keys: list[str], db, nrtl=None, log=None, vapour: str = "ideal"):
        self.keys, self.db, self.nrtl, self.log = list(keys), db, nrtl, log
        self.comps = [db.get(k) for k in self.keys]
        self.volatile = np.array([bool(c.raw.get("antoine")) and c.has_phase("gas") for c in self.comps])
        if not self.volatile.any():
            raise DataError("No volatile component (no vapour-pressure data).")
        if vapour not in ("ideal", "lewis_randall"):
            raise ValueError("vapour must be 'ideal' or 'lewis_randall'")
        self.vapour = vapour
        self._virial = None
        if vapour == "lewis_randall":
            self._virial = [VirialEOS([c.critical]) if (v and c.critical) else None
                            for c, v in zip(self.comps, self.volatile)]

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

    def K(self, T: float, P: float, x) -> np.ndarray:
        ps = self.psat(T)
        K = self.gamma(x, T) * ps / P
        if self._virial is not None:          # Φ_i = exp[B_ii (P − P_i^sat)/RT]  (Lewis–Randall, no Poynting)
            for i, v in enumerate(self._virial):
                if v is not None:
                    K[i] /= math.exp(v.B_mix(T, [1.0])[0] * (P - ps[i]) / (R * T))
        return np.where(self.volatile, K, 0.0)

    # ------------------------------------------------------------------ bubble / dew points
    def _bracket(self, f, P: float) -> tuple[float, float]:
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
        x = np.asarray(x, float)
        return float(np.sum(self.K(T, 1.0, x) * x))         # K ∝ 1/P for an ideal vapour

    def _dew_x(self, T: float, P: float, y) -> np.ndarray:
        x = np.asarray(y, float).copy()
        for _ in range(200):
            K = self.K(T, P, x)
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
            K = self.K(T, P, x)
            return 1.0 - float(np.sum(np.divide(y, K, out=np.zeros_like(y), where=K > 0)))
        T = brentq(f, *self._bracket(f, P), xtol=1e-10)
        self._check_range(T)
        return T

    # ------------------------------------------------------------------ isothermal flash
    def flash(self, T: float, P: float, z) -> tuple[float, np.ndarray, np.ndarray]:
        """(vapour fraction V, liquid x, vapour y) at T, P for overall composition z."""
        z = np.asarray(z, float) / np.sum(z)
        x = z.copy()
        V = 0.5
        for _ in range(500):
            K = self.K(T, P, x)
            g = lambda v: float(np.sum(z * (K - 1.0) / (1.0 + v * (K - 1.0))))
            if g(0.0) <= 0.0:
                V = 0.0                                     # sub-cooled liquid
            elif g(1.0) >= 0.0:
                V = 1.0                                     # superheated vapour
            else:
                V = brentq(g, 0.0, 1.0, xtol=1e-14)
            x_new = z / (1.0 + V * (K - 1.0))
            x_new /= x_new.sum()
            if np.max(np.abs(x_new - x)) < 1e-13:
                x = x_new
                break
            x = x_new
        K = self.K(T, P, x)
        y = K * x
        y = y / y.sum() if y.sum() > 0 else y
        return V, x, y
