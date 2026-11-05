"""Partial molar properties, mixing functions and the Gibbs–Duhem relation (CL 202, Ch. 4).

For an extensive property B(T, P, N_1 … N_n) the partial molar property is
    B̄_i = (∂B/∂N_i)_{T,P,N_j≠i},      and Euler's theorem gives  B = Σ_i N_i B̄_i   (B/N = Σ_i x_i B̄_i).
Given the molar property B(x) the same quantity is
    B̄_k = B + (∂B/∂x_k) − Σ_{j≠i} x_j (∂B/∂x_j)        (x_i eliminated with Σ x = 1).
Gibbs–Duhem at constant T and P:  Σ_i x_i (∂B̄_i/∂x_j) = 0.

Applied to the excess Gibbs energy of the NRTL model:
    RT ln γ_i = (∂(N G^E)/∂N_i)_{T,P,N_j}              (γ_i is the partial molar G^E / RT)
    H̄_i^E = −RT² (∂ ln γ_i/∂T)_{P,x}                    (Gibbs–Helmholtz, partial molar form)
    C̄_p,i = C_p,i (pure) + C̄_p,i^E,   C̄_p,i^E = (∂H̄_i^E/∂T)_{P,x}
Mixing functions with pure liquids at the same T and P as reference:
    ΔG_mix = RT Σ x_i ln x_i + G^E,  ΔH_mix = H^E,  ΔS_mix = (ΔH_mix − ΔG_mix)/T,  ΔV_mix = V^E (= 0 here).

These are computed numerically from the molar model, so they test the model's thermodynamic
consistency rather than repeating its formulas (see tests/test_partial_molar.py).
"""
from __future__ import annotations

from typing import Callable, Sequence

import numpy as np

R = 8.314462618


def partial_molar(B_molar: Callable[[np.ndarray], float], x: Sequence[float], h: float = 1e-6) -> np.ndarray:
    """B̄_i for all components from a molar property B(x) by differentiating N·B(N/ΣN) w.r.t. N_i."""
    x = np.asarray(x, float)
    out = np.empty(len(x))
    for i in range(len(x)):
        up, dn = x.copy(), x.copy()
        up[i] += h
        dn[i] -= h
        Bup = up.sum() * B_molar(up / up.sum())
        Bdn = dn.sum() * B_molar(dn / dn.sum())
        out[i] = (Bup - Bdn) / (2 * h)
    return out


def partial_molar_from_x(B_molar: Callable[[np.ndarray], float], x: Sequence[float], k: int,
                         i_ref: int = -1, h: float = 1e-6) -> float:
    """B̄_k from B(x) with the mole-fraction formula of the lecture (x_{i_ref} eliminated)."""
    x = np.asarray(x, float)
    n = len(x)
    i_ref %= n

    def dB(j: int) -> float:            # (∂B/∂x_j) at constant x_l (l ≠ j, i_ref), x_ref = 1 − Σ others
        if j == i_ref:
            return 0.0
        up, dn = x.copy(), x.copy()
        up[j] += h
        up[i_ref] -= h
        dn[j] -= h
        dn[i_ref] += h
        return (B_molar(up) - B_molar(dn)) / (2 * h)

    d = np.array([dB(j) for j in range(n)])
    return float(B_molar(x) + d[k] - sum(x[j] * d[j] for j in range(n) if j != i_ref))


def gibbs_duhem_residual(Bbar: Callable[[np.ndarray], np.ndarray], x: Sequence[float], j: int,
                         i_ref: int = -1, h: float = 1e-6) -> float:
    """Σ_i x_i (∂B̄_i/∂x_j)_{T,P} — zero for a thermodynamically consistent model."""
    x = np.asarray(x, float)
    i_ref %= len(x)
    up, dn = x.copy(), x.copy()
    up[j] += h
    up[i_ref] -= h
    dn[j] -= h
    dn[i_ref] += h
    return float(x @ ((Bbar(up) - Bbar(dn)) / (2 * h)))


# ---------------------------------------------------------------------------
# NRTL-specific quantities
# ---------------------------------------------------------------------------
def ln_gamma_from_gE(model, x: Sequence[float], T: float) -> np.ndarray:
    """ln γ_i = ∂(N g^E/RT)/∂N_i, evaluated numerically (cross-check of the analytic NRTL ln γ)."""
    return partial_molar(lambda z: model.gE_RT(list(z), T), x)


def partial_molar_hE(model, x: Sequence[float], T: float) -> np.ndarray:
    """H̄_i^E = −RT² (∂ ln γ_i/∂T)_{P,x}  (J/mol)."""
    h = 1e-3 * T
    d = (np.asarray(model.ln_gamma(list(x), T + h)) - np.asarray(model.ln_gamma(list(x), T - h))) / (2 * h)
    return -R * T * T * d


def partial_molar_cpE(model, x: Sequence[float], T: float) -> np.ndarray:
    """C̄_p,i^E = (∂H̄_i^E/∂T)_{P,x}  (J/(mol K))."""
    h = 1e-2 * T
    return (partial_molar_hE(model, x, T + h) - partial_molar_hE(model, x, T - h)) / (2 * h)


def mixing_functions(model, x: Sequence[float], T: float) -> dict[str, float]:
    """ΔG_mix, ΔH_mix, ΔS_mix (J/mol, J/(mol K)) with pure liquids at T, P as reference states."""
    x = np.asarray(x, float)
    ideal = R * T * float(sum(xi * np.log(xi) for xi in x if xi > 0))
    gE = R * T * model.gE_RT(list(x), T) if model is not None else 0.0
    hE = model.hE(list(x), T) if model is not None else 0.0
    dG = ideal + gE
    return {"dG_mix": dG, "dH_mix": hE, "dS_mix": (hE - dG) / T, "dG_mix_ideal": ideal,
            "dS_mix_ideal": -ideal / T, "gE": gE, "sE": (hE - gE) / T}
