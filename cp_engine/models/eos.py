"""Equations of state for gases and vapours: residual (departure) properties and fugacity coefficients.

Two families from the CL 202 lectures:

* **Virial** (pressure form, truncated after B):  Z = 1 + B P / (R T).
  Pure B from the Pitzer/Abbott correlation  B P_c/(R T_c) = B⁰ + ω B¹,
  B⁰ = 0.083 − 0.422/T_r^1.6,  B¹ = 0.139 − 0.172/T_r^4.2.
  Mixtures: the exact quadratic rule B_mix = ΣΣ y_i y_j B_ij, with the cross coefficient B_ij from
  the Prausnitz combining rules  T_c,ij = √(T_c,i T_c,j)(1 − k_ij),  V_c,ij = [(V_c,i^⅓ + V_c,j^⅓)/2]³,
  Z_c,ij = (Z_c,i + Z_c,j)/2,  ω_ij = (ω_i + ω_j)/2,  P_c,ij = Z_c,ij R T_c,ij / V_c,ij.
  Residual properties:  H^R = P (B − T dB/dT),  C_p^R = −T P d²B/dT²,  ln φ̂_i = (P/RT)(2 Σ_j y_j B_ij − B_mix).

* **Generic cubic**  P = RT/(V − b) − a(T)/[(V + εb)(V + σb)]  with the van der Waals, Redlich–Kwong,
  Soave–Redlich–Kwong and Peng–Robinson constants (ε, σ, Ω_a, Ω_b, α(T_r, ω)).
  Mixtures: van der Waals one-fluid rules  a = ΣΣ x_i x_j √(a_i a_j)(1 − k_ij),  b = Σ x_i b_i.
  With β = bP/RT, q = a/(bRT) and I = ln[(Z + σβ)/(Z + εβ)]/(σ − ε)  (I = β/Z if ε = σ):
      ln φ   = Z − 1 − ln(Z − β) − q I                       (pure)
      ln φ̂_i = (b_i/b)(Z − 1) − ln(Z − β) − q̄_i I,  q̄_i = q (2 Σ_j x_j a_ij / a − b_i / b)
      H^R/RT = Z − 1 + (d ln a / d ln T − 1) q I
  and C_p^R = (∂H^R/∂T)_P by central differences.  The vapour is the largest root of the cubic in Z,
  the liquid the smallest.

Phase equilibrium from one equation for both phases (φ–φ): bubble_P() iterates K_i = φ̂_i^L/φ̂_i^V from
Wilson K-values (the lecture's bubble-pressure algorithm).

C_p(real gas) = C_p^ig(T) + C_p^R(T, P, y): the ideal-gas part comes from the database correlations,
the residual part from the EOS.  For an ideal gas every residual property is zero.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

R = 8.314462618

CUBICS = {
    #        ε              σ              Ω_a      Ω_b
    "vdW": (0.0, 0.0, 27.0 / 64.0, 1.0 / 8.0),
    "RK": (0.0, 1.0, 0.42748, 0.08664),
    "SRK": (0.0, 1.0, 0.42748, 0.08664),
    "PR": (1.0 - math.sqrt(2.0), 1.0 + math.sqrt(2.0), 0.45724, 0.07780),
}


@dataclass(frozen=True)
class Critical:
    """Critical constants of one component (SI: K, Pa, m³/mol)."""
    name: str
    Tc: float
    Pc: float
    omega: float
    Vc: float | None = None

    @property
    def Zc(self) -> float | None:
        return None if self.Vc is None else self.Pc * self.Vc / (R * self.Tc)


def _alpha(kind: str, Tr: float, omega: float) -> tuple[float, float]:
    """α(T_r) and d ln α / d ln T for the four cubics."""
    if kind == "vdW":
        return 1.0, 0.0
    if kind == "RK":
        return Tr ** -0.5, -0.5
    if kind == "SRK":
        k = 0.480 + 1.574 * omega - 0.176 * omega ** 2
    else:  # PR
        k = 0.37464 + 1.54226 * omega - 0.26992 * omega ** 2
    s = 1.0 + k * (1.0 - math.sqrt(Tr))
    return s * s, -k * math.sqrt(Tr) / s


class CubicEOS:
    """vdW / RK / SRK / PR equation of state for a pure fluid or a mixture."""

    def __init__(self, comps: list[Critical], kind: str = "PR", kij: dict[tuple[int, int], float] | None = None):
        if kind not in CUBICS:
            raise ValueError(f"Unknown cubic EOS '{kind}'; choose from {', '.join(CUBICS)}.")
        self.comps, self.kind = list(comps), kind
        self.eps, self.sig, self.Oa, self.Ob = CUBICS[kind]
        self.kij = kij or {}
        self.b_i = np.array([self.Ob * R * c.Tc / c.Pc for c in self.comps])

    # ---------------------------------------------------------------- parameters
    def _a_i(self, T: float) -> tuple[np.ndarray, np.ndarray]:
        """a_i(T) and d ln a_i / d ln T for each component."""
        a, dl = [], []
        for c in self.comps:
            al, d = _alpha(self.kind, T / c.Tc, c.omega)
            a.append(self.Oa * (R * c.Tc) ** 2 / c.Pc * al)
            dl.append(d)
        return np.array(a), np.array(dl)

    def _mix(self, T: float, x) -> tuple[float, float, np.ndarray, float]:
        """a, b of the mixture, Σ_j x_j a_ij for each i, and d ln a / d ln T."""
        x = np.asarray(x, float)
        a_i, dl_i = self._a_i(T)
        n = len(a_i)
        k = np.zeros((n, n))
        for (i, j), v in self.kij.items():
            k[i, j] = k[j, i] = v
        a_ij = np.sqrt(np.outer(a_i, a_i)) * (1.0 - k)
        a = float(x @ a_ij @ x)
        # d a_ij / d ln T = a_ij (dl_i + dl_j)/2
        da = float(x @ (a_ij * 0.5 * (dl_i[:, None] + dl_i[None, :])) @ x)
        return a, float(x @ self.b_i), a_ij @ x, da / a

    # ---------------------------------------------------------------- roots
    def Z_roots(self, T: float, P: float, x) -> list[float]:
        """Real roots of the cubic in Z with Z > β (physically meaningful)."""
        a, b, _, _ = self._mix(T, x)
        beta, q = b * P / (R * T), a / (b * R * T)
        e, s = self.eps, self.sig
        # Z³ + c2 Z² + c1 Z + c0 = 0 (generic form of Smith, Van Ness & Abbott)
        c2 = (e + s - 1.0) * beta - 1.0
        c1 = e * s * beta ** 2 - (e + s) * beta * (beta + 1.0) + q * beta
        c0 = -(e * s * beta ** 2 * (beta + 1.0) + q * beta ** 2)
        r = np.roots([1.0, c2, c1, c0])
        real = sorted(float(z.real) for z in r if abs(z.imag) < 1e-9 and z.real > beta)
        return real

    def Z(self, T: float, P: float, x, phase: str = "vapour") -> float:
        roots = self.Z_roots(T, P, x)
        if not roots:
            raise ValueError("no physical root of the cubic equation of state")
        return roots[-1] if phase == "vapour" else roots[0]

    def _I(self, Z: float, beta: float) -> float:
        if self.eps == self.sig:
            return beta / (Z + self.eps * beta)
        return math.log((Z + self.sig * beta) / (Z + self.eps * beta)) / (self.sig - self.eps)

    # ---------------------------------------------------------------- properties
    def ln_phi(self, T: float, P: float, x, phase: str = "vapour") -> np.ndarray:
        """ln φ̂_i of every component (pure fluid: ln φ)."""
        x = np.asarray(x, float)
        a, b, sum_xa, _ = self._mix(T, x)
        Z = self.Z(T, P, x, phase)
        beta, q = b * P / (R * T), a / (b * R * T)
        I = self._I(Z, beta)
        qbar = q * (2.0 * sum_xa / a - self.b_i / b)
        return self.b_i / b * (Z - 1.0) - math.log(Z - beta) - qbar * I

    def H_residual(self, T: float, P: float, x, phase: str = "vapour") -> float:
        """Residual molar enthalpy H − H^ig (J/mol)."""
        a, b, _, dlna = self._mix(T, x)
        Z = self.Z(T, P, x, phase)
        beta, q = b * P / (R * T), a / (b * R * T)
        return R * T * (Z - 1.0 + (dlna - 1.0) * q * self._I(Z, beta))

    def Cp_residual(self, T: float, P: float, x, phase: str = "vapour") -> float:
        """Residual molar heat capacity C_p − C_p^ig (J/(mol K)) = (∂H^R/∂T)_P."""
        h = 1e-3 * T
        return (self.H_residual(T + h, P, x, phase) - self.H_residual(T - h, P, x, phase)) / (2 * h)

    def molar_volume(self, T: float, P: float, x, phase: str = "vapour") -> float:
        return self.Z(T, P, x, phase) * R * T / P

    def bubble_P(self, T: float, x, P0: float | None = None) -> tuple[float, np.ndarray]:
        """Bubble pressure and vapour composition by the φ–φ method (one EOS for both phases):
        guess P and y from Wilson K-values, K_i = φ̂_i^L/φ̂_i^V, y_i = K_i x_i, scale P by Σ K_i x_i
        until Σ y_i = 1. Needs no vapour pressure, activity coefficient or standard state."""
        x = np.asarray(x, float)
        Kw = np.array([c.Pc * math.exp(5.373 * (1.0 + c.omega) * (1.0 - c.Tc / T)) for c in self.comps])
        P = P0 or float(np.sum(Kw * x))
        K = Kw / P
        for _ in range(500):
            y = K * x / np.sum(K * x)
            K = np.exp(self.ln_phi(T, P, x, "liquid") - self.ln_phi(T, P, y, "vapour"))
            S = float(np.sum(K * x))
            P *= S
            if abs(S - 1.0) < 1e-11:
                break
        if np.max(np.abs(np.log(K))) < 1e-3:              # all K_i = 1: both roots are the same phase
            raise ValueError("trivial solution x = y (near the mixture critical point)")
        return P, K * x / np.sum(K * x)

    def psat(self, T: float) -> float:
        """Vapour pressure of a pure fluid from the equal-fugacity condition φ_L = φ_V (Pa)."""
        if len(self.comps) != 1:
            raise ValueError("psat() is defined for a pure fluid")
        c = self.comps[0]
        if T >= c.Tc:
            raise ValueError("T ≥ T_c: no vapour pressure")
        # initial guess from the definition of ω (Lee–Kesler-type: log10 Pr = 7/3 (1 + ω)(1 − 1/Tr))
        P = c.Pc * 10 ** (7.0 / 3.0 * (1.0 + c.omega) * (1.0 - c.Tc / T))
        for _ in range(200):
            roots = self.Z_roots(T, P, [1.0])
            if len(roots) < 2:                     # only one phase at this P: move P towards the 3-root region
                P *= 0.9 if roots and roots[0] < 0.3 else 1.1
                continue
            lv = self.ln_phi(T, P, [1.0], "vapour")[0]
            ll = self.ln_phi(T, P, [1.0], "liquid")[0]
            P_new = P * math.exp(ll - lv)
            if abs(P_new / P - 1.0) < 1e-10:
                return P_new
            P = P_new
        raise RuntimeError(f"psat did not converge for {c.name} at {T} K")


class VirialEOS:
    """Truncated virial equation Z = 1 + BP/RT with Pitzer/Abbott B and Prausnitz combining rules."""

    def __init__(self, comps: list[Critical], kij: dict[tuple[int, int], float] | None = None):
        self.comps = list(comps)
        n = len(self.comps)
        if n > 1 and any(c.Vc is None for c in self.comps):
            missing = ", ".join(c.name for c in self.comps if c.Vc is None)
            raise ValueError(f"The virial cross coefficient needs V_c (missing for {missing}).")
        self.cross: dict[tuple[int, int], tuple[float, float, float]] = {}
        k = kij or {}
        for i in range(n):
            for j in range(i, n):
                ci, cj = self.comps[i], self.comps[j]
                if i == j:
                    self.cross[i, j] = (ci.Tc, ci.Pc, ci.omega)
                    continue
                Tc = math.sqrt(ci.Tc * cj.Tc) * (1.0 - k.get((i, j), k.get((j, i), 0.0)))
                Vc = ((ci.Vc ** (1 / 3) + cj.Vc ** (1 / 3)) / 2.0) ** 3
                Zc = (ci.Zc + cj.Zc) / 2.0
                self.cross[i, j] = (Tc, Zc * R * Tc / Vc, (ci.omega + cj.omega) / 2.0)

    @staticmethod
    def _B(T: float, Tc: float, Pc: float, w: float) -> tuple[float, float, float]:
        """B, dB/dT, d²B/dT² (m³/mol) from the Abbott correlation."""
        Tr = T / Tc
        B0, B1 = 0.083 - 0.422 / Tr ** 1.6, 0.139 - 0.172 / Tr ** 4.2
        dB0, dB1 = 0.6752 / Tr ** 2.6, 0.7224 / Tr ** 5.2
        d2B0, d2B1 = -1.75552 / Tr ** 3.6, -3.75648 / Tr ** 6.2
        f = R * Tc / Pc
        return f * (B0 + w * B1), f * (dB0 + w * dB1) / Tc, f * (d2B0 + w * d2B1) / Tc ** 2

    def B_matrix(self, T: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        n = len(self.comps)
        B, dB, d2B = (np.zeros((n, n)) for _ in range(3))
        for (i, j), (Tc, Pc, w) in self.cross.items():
            B[i, j], dB[i, j], d2B[i, j] = self._B(T, Tc, Pc, w)
            B[j, i], dB[j, i], d2B[j, i] = B[i, j], dB[i, j], d2B[i, j]
        return B, dB, d2B

    def B_mix(self, T: float, y) -> tuple[float, float, float]:
        y = np.asarray(y, float)
        B, dB, d2B = self.B_matrix(T)
        return float(y @ B @ y), float(y @ dB @ y), float(y @ d2B @ y)

    def Z(self, T: float, P: float, y, phase: str = "vapour") -> float:
        return 1.0 + self.B_mix(T, y)[0] * P / (R * T)

    def ln_phi(self, T: float, P: float, y, phase: str = "vapour") -> np.ndarray:
        y = np.asarray(y, float)
        B = self.B_matrix(T)[0]
        return P / (R * T) * (2.0 * B @ y - float(y @ B @ y))

    def H_residual(self, T: float, P: float, y, phase: str = "vapour") -> float:
        B, dB, _ = self.B_mix(T, y)
        return P * (B - T * dB)

    def Cp_residual(self, T: float, P: float, y, phase: str = "vapour") -> float:
        return -T * P * self.B_mix(T, y)[2]

    def molar_volume(self, T: float, P: float, y, phase: str = "vapour") -> float:
        return R * T / P + self.B_mix(T, y)[0]


def make_eos(comps: list[Critical], kind: str = "PR", kij=None):
    """'PR', 'SRK', 'RK', 'vdW' (cubic) or 'virial'."""
    if kind == "virial":
        return VirialEOS(comps, kij)
    return CubicEOS(comps, kind, kij)
