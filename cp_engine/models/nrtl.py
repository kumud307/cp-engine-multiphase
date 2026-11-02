"""Multicomponent NRTL model: G^E, activity coefficients, H^E and Cp^E.

NRTL (Renon & Prausnitz, AIChE J. 14 (1968) 135)
------------------------------------------------
    G^E/(RT) = Σ_i x_i [ Σ_j τ_ji G_ji x_j ] / [ Σ_k G_ki x_k ]

    τ_ij  = a_ij + b_ij/T + e_ij ln T + f_ij T          (τ_ii = 0)
    α_ij  = c_ij + d_ij (T − 273.15 K)                  (α_ij = α_ji)
    G_ij  = exp(−α_ij τ_ij)

Excess enthalpy and heat capacity (Gibbs–Helmholtz)
---------------------------------------------------
    H^E   = −R T² ∂(G^E/RT)/∂T              at constant P, x
    Cp^E  = (∂H^E/∂T)_P,x = −R [ 2T g' + T² g'' ],   g ≡ G^E/(RT)

The real-solution heat capacity is then Cp = Σ x_i Cp_i + Cp^E.

Numerical derivatives
---------------------
g' and g'' are obtained by central differences with step h:

    g'  ≈ [g(T+h) − g(T−h)] / (2h)              error O(h²)
    g'' ≈ [g(T+h) − 2g(T) + g(T−h)] / h²         error O(h²)

Truncation error grows like h², round-off error like ε/h² (ε ≈ 2e-16 is
machine precision), so the total error of g'' is smallest near
h ≈ ε^(1/4)·T ≈ 1e-4·T  (≈ 0.03 K at 300 K).  That is the default step.
`derivative_convergence()` demonstrates this by halving h; the test suite
cross-checks against exact symbolic derivatives (sympy).

If a pair has no parameters in the database, τ_ij = τ_ji = 0 for that pair,
i.e. the pair mixes ideally, and a warning is recorded.
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from ..database import R, Database, check_range
from ..report import Log, null_log

DEFAULT_REL_STEP = 1e-4     # h = 1e-4 * T  (see module docstring)


class NRTL:
    """NRTL model for an ordered list of component keys."""

    def __init__(self, keys: Sequence[str], db: Database | None = None,
                 params: dict | None = None, log: Log | None = None, prefer: str = "HE"):
        """Build parameter matrices.

        Parameters
        ----------
        keys   : component keys (order defines the index i)
        db     : database to read pair parameters from
        params : optional explicit matrices {'a','b','e','f','c','d'} (N×N arrays),
                 used by tests and by the parameter-fitting routine
        prefer : which parameter set to use when a pair has several ("HE" = fitted to
                 excess enthalpies, best for Cp; "VLE" = fitted to phase equilibria)
        """
        self.keys = list(keys)
        n = len(self.keys)
        self.log = log or null_log()
        self.a, self.b, self.e, self.f = (np.zeros((n, n)) for _ in range(4))
        self.c, self.d = np.full((n, n), 0.3), np.zeros((n, n))
        self.missing_pairs: list[tuple[str, str]] = []
        self.pair_sources: dict[tuple[str, str], dict] = {}
        if params is not None:
            for name in ("a", "b", "e", "f", "c", "d"):
                if name in params:
                    setattr(self, name, np.array(params[name], dtype=float))
            return
        if db is None:
            return
        for i in range(n):
            for j in range(i + 1, n):
                ki, kj = self.keys[i], self.keys[j]
                p = db.nrtl_pair(ki, kj, prefer)
                if p is None:
                    self.missing_pairs.append((ki, kj))
                    continue
                self.pair_sources[(ki, kj)] = p
                for s in ("a", "b", "e", "f"):
                    getattr(self, s)[i, j] = p.get(f"{s}_ij", 0.0)
                    getattr(self, s)[j, i] = p.get(f"{s}_ji", 0.0)
                self.c[i, j] = self.c[j, i] = p.get("c", 0.3)
                self.d[i, j] = self.d[j, i] = p.get("d", 0.0)
        for ki, kj in self.missing_pairs:
            self.log.warn(f"No NRTL parameters for {ki}–{kj}: this pair is treated as ideal (Cp^E contribution = 0).")

    # ------------------------------------------------------------------
    @property
    def has_any_parameters(self) -> bool:
        return bool(np.any(self.a) or np.any(self.b) or np.any(self.e) or np.any(self.f))

    def check_temperature(self, T: float) -> None:
        """Warn if T is outside the fitted range of any parameter pair."""
        for (ki, kj), p in self.pair_sources.items():
            check_range(T, p.get("Tmin"), p.get("Tmax"), f"NRTL parameters {ki}–{kj}", self.log, strict=False)

    def tau(self, T: float) -> np.ndarray:
        t = self.a + self.b / T + self.e * np.log(T) + self.f * T
        np.fill_diagonal(t, 0.0)
        return t

    def alpha(self, T: float) -> np.ndarray:
        return self.c + self.d * (T - 273.15)

    def G(self, T: float) -> np.ndarray:
        return np.exp(-self.alpha(T) * self.tau(T))

    # ------------------------------------------------------------------
    def gE_RT(self, x: Sequence[float], T: float) -> float:
        """Dimensionless excess Gibbs energy g = G^E/(RT) per mole of solution."""
        x = np.asarray(x, dtype=float)
        tau, G = self.tau(T), self.G(T)
        num = (tau * G).T @ x          # num_i = Σ_j τ_ji G_ji x_j
        den = G.T @ x                  # den_i = Σ_k G_ki x_k
        return float(np.sum(x * num / den))

    def ln_gamma(self, x: Sequence[float], T: float) -> np.ndarray:
        """Activity coefficients ln γ_i (multicomponent NRTL)."""
        x = np.asarray(x, dtype=float)
        tau, G = self.tau(T), self.G(T)
        den = G.T @ x                              # Σ_k G_ki x_k   (index i)
        num = (tau * G).T @ x                      # Σ_j τ_ji G_ji x_j
        term1 = num / den
        # term2_i = Σ_j x_j G_ij / den_j * (τ_ij − num_j / den_j)
        term2 = (G * (tau - (num / den)[None, :])) @ (x / den)
        return term1 + term2

    def gamma(self, x: Sequence[float], T: float) -> np.ndarray:
        return np.exp(self.ln_gamma(x, T))

    # ------------------------------------------------------------------
    def _derivs(self, x, T, h):
        h = h if h is not None else DEFAULT_REL_STEP * T
        gp, g0, gm = self.gE_RT(x, T + h), self.gE_RT(x, T), self.gE_RT(x, T - h)
        d1 = (gp - gm) / (2 * h)
        d2 = (gp - 2 * g0 + gm) / h ** 2
        return g0, d1, d2

    def gE(self, x: Sequence[float], T: float) -> float:
        """Excess Gibbs energy G^E, J/mol."""
        return R * T * self.gE_RT(x, T)

    def hE(self, x: Sequence[float], T: float, h: float | None = None) -> float:
        """Excess enthalpy H^E = −R T² ∂(G^E/RT)/∂T, J/mol."""
        _, d1, _ = self._derivs(x, T, h)
        return -R * T ** 2 * d1

    def cpE(self, x: Sequence[float], T: float, h: float | None = None) -> float:
        """Excess heat capacity Cp^E = −R (2 T g' + T² g''), J/(mol K)."""
        _, d1, d2 = self._derivs(x, T, h)
        return -R * (2 * T * d1 + T ** 2 * d2)

    def derivative_convergence(self, x: Sequence[float], T: float,
                               steps: Sequence[float] = (1.0, 0.3, 0.1, 0.03, 0.01, 3e-3, 1e-3, 1e-4, 1e-5)):
        """Cp^E for a sequence of step sizes h (K) – shows the O(h²) plateau
        and the round-off blow-up at very small h."""
        return [(h, self.hE(x, T, h), self.cpE(x, T, h)) for h in steps]


    def is_stable_binary(self, T: float, n: int = 2000) -> bool:
        """Phase-stability test for a binary: the Gibbs energy of mixing
        ΔG_mix/RT = x ln x + (1−x) ln(1−x) + G^E/RT must be convex in x
        (d²/dx² > 0 everywhere). A non-convex curve means the model predicts a
        liquid–liquid split – unphysical for fully miscible pairs such as water–ethanol."""
        xs = np.linspace(0.5 / n, 1 - 0.5 / n, n)
        g = [x * np.log(x) + (1 - x) * np.log(1 - x) + self.gE_RT([x, 1 - x], T) for x in xs]
        return bool(np.min(np.diff(g, 2)) > 0)


# ---------------------------------------------------------------------------
# Parameter fitting to excess-property data
# ---------------------------------------------------------------------------
def fit_nrtl_to_excess(x1: Sequence[float], T: Sequence[float], values: Sequence[float],
                       kind: str = "HE", alpha: float = 0.3,
                       fit: Sequence[str] = ("a12", "a21", "b12", "b21"),
                       initial: dict | None = None) -> dict:
    """Fit binary NRTL parameters to H^E (J/mol) or Cp^E (J/(mol K)) data.

    VLE-fitted parameters rarely reproduce H^E and Cp^E because those depend
    on the *temperature derivatives* of τ.  This routine fits the chosen
    parameters (any of a12, a21, b12, b21, e12, e21, f12, f21) by least
    squares and returns a dict ready to paste into data/binary_nrtl.json.

    Use data at SEVERAL temperatures: H^E at a single T leaves the b-parameters
    nearly degenerate (different (b12, b21) pairs fit almost equally well).

    x1     : mole fraction of component 1
    T      : temperatures (K)
    values : measured H^E or Cp^E
    kind   : "HE" or "CpE"
    """
    from scipy.optimize import least_squares

    x1, T, values = (np.asarray(v, dtype=float) for v in (x1, T, values))
    init = {"a12": 0.0, "a21": 0.0, "b12": 0.0, "b21": 0.0,
            "e12": 0.0, "e21": 0.0, "f12": 0.0, "f21": 0.0}
    init.update(initial or {})

    def model(p):
        q = dict(init)
        q.update(dict(zip(fit, p)))
        mats = {s: np.array([[0.0, q[f"{s}12"]], [q[f"{s}21"], 0.0]]) for s in ("a", "b", "e", "f")}
        mats["c"] = np.full((2, 2), alpha)
        m = NRTL(["1", "2"], params=mats)
        f = m.hE if kind == "HE" else m.cpE
        return np.array([f([xi, 1 - xi], Ti) for xi, Ti in zip(x1, T)]), q

    res = least_squares(lambda p: model(p)[0] - values, x0=[init[k] for k in fit])
    pred, q = model(res.x)
    rmse = float(np.sqrt(np.mean((pred - values) ** 2)))
    return {"a_ij": q["a12"], "a_ji": q["a21"], "b_ij": q["b12"], "b_ji": q["b21"],
            "e_ij": q["e12"], "e_ji": q["e21"], "f_ij": q["f12"], "f_ji": q["f21"],
            "c": alpha, "d": 0.0, "fitted_to": kind, "rmse": rmse, "success": bool(res.success)}


def fit_nrtl_joint(he: tuple, sle: tuple | None = None, alpha: float = 0.3,
                   fit: Sequence[str] = ("a12", "a21", "b12", "b21", "e12", "e21"),
                   initial: dict | None = None, sigma_he: float = 30.0, sigma_sle: float = 0.01,
                   cpe: tuple | None = None, sigma_cpe: float = 0.5) -> dict:
    """Fit binary NRTL parameters to excess enthalpies and (optionally) solid–liquid equilibrium data.

    he  = (x1, T, H^E in J/mol)
    cpe = (x1, T, Cp^E in J/(mol K)), optional excess heat capacities
    sle = (x1, T_freeze, Tm, dHfus, dCp, k): freezing points where pure component k (0 or 1)
          crystallises; residual = [ln x_k + ln γ_k(x, T) − sle_rhs(T)] / sigma_sle
          (sigma_sle = 0.01 in ln x corresponds to about 1 K for ice).
    Fitting both data types makes ONE liquid model describe the heat of mixing (→ Cp^E) and the
    freezing equilibrium (→ latent heat of ice slurries) consistently.
    """
    from scipy.optimize import least_squares
    from .phase_change import sle_rhs

    x1, T, H = (np.asarray(v, dtype=float) for v in he)
    init = {k: 0.0 for k in ("a12", "a21", "b12", "b21", "e12", "e21", "f12", "f21")}
    init.update(initial or {})

    def model(p):
        q = dict(init)
        q.update(dict(zip(fit, p)))
        mats = {s: np.array([[0.0, q[f"{s}12"]], [q[f"{s}21"], 0.0]]) for s in ("a", "b", "e", "f")}
        mats["c"] = np.full((2, 2), alpha)
        return NRTL(["1", "2"], params=mats), q

    def resid(p):
        m, _ = model(p)
        r = [(m.hE([x, 1 - x], t) - v) / sigma_he for x, t, v in zip(x1, T, H)]
        if cpe is not None:
            r += [(m.cpE([x, 1 - x], t) - v) / sigma_cpe for x, t, v in zip(*cpe)]
        if sle is not None:
            xs, Ts, Tm, dH, dCp, k = sle
            for x, t in zip(xs, Ts):
                xv = [x, 1 - x]
                r.append((np.log(xv[k]) + m.ln_gamma(xv, t)[k] - sle_rhs(t, Tm, dH, dCp)) / sigma_sle)
        return np.array(r)

    res = least_squares(resid, x0=[init[k] for k in fit], max_nfev=600)
    m, q = model(res.x)
    rmse = float(np.sqrt(np.mean([(m.hE([x, 1 - x], t) - v) ** 2 for x, t, v in zip(x1, T, H)])))
    return {"a_ij": q["a12"], "a_ji": q["a21"], "b_ij": q["b12"], "b_ji": q["b21"], "e_ij": q["e12"],
            "e_ji": q["e21"], "f_ij": q["f12"], "f_ji": q["f21"], "c": alpha, "d": 0.0,
            "fitted_to": "HE" + ("+CpE" if cpe is not None else "") + ("+SLE" if sle is not None else ""),
            "rmse": rmse, "success": bool(res.success),
            "model": m}


def fit_nrtl_best(he: tuple, sle: tuple | None, initial: dict, alphas=(0.2, 0.3, 0.4, 0.47),
                  forms=(("a12", "a21", "b12", "b21"), ("a12", "a21", "b12", "b21", "e12", "e21")),
                  T_stable=(220, 390), cpe: tuple | None = None) -> dict:
    """Try several α values and τ(T) forms; keep the lowest-H^E-error set that is phase-stable
    (no liquid–liquid split) over T_stable. Returns the fit dict (see fit_nrtl_joint)."""
    best = None
    for alpha in alphas:
        for form in forms:
            f = fit_nrtl_joint(he, sle, alpha, form, initial, cpe=cpe)
            if all(f["model"].is_stable_binary(T, 600) for T in range(T_stable[0], T_stable[1] + 1, 10)):
                if best is None or f["rmse"] < best["rmse"]:
                    best = f
    if best is None:
        raise ValueError("no phase-stable NRTL parameter set found")
    return best
