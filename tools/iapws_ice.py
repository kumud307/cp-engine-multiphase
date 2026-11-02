"""IAPWS R10-06(2009) equation of state for ice Ih (Feistel & Wagner 2006) – used to generate and
verify the ice data in data/components.json.

    g(T,p) = g0(p) − s0·Tt·τ + Tt·Re Σ_{k=1,2} r_k [ (t_k−τ)ln(t_k−τ) + (t_k+τ)ln(t_k+τ) − 2t_k ln t_k − τ²/t_k ]
    τ = T/Tt,  π = p/pt,  g0(p) = Σ_{k=0..4} g0k (π−π0)^k,  r2(p) = Σ_{k=0..2} r2k (π−π0)^k,  r1 constant

Properties: cp = −T g_TT, ρ = 1/g_p, α = g_Tp/g_p, κ_T = −g_pp/g_p.
Coefficients: Table 2 of the release; check values: Table 6 (https://iapws.org/technical-guidance/release/Ice-2009).

Run:  python tools/iapws_ice.py   → verifies Table 6 and prints the fit stored in components.json
"""
import numpy as np

Tt, pt, p0 = 273.16, 611.657, 101325.0
g0 = [-0.632020233335886e6, 0.655022213658955, -0.189369929326131e-7, 0.339746123271053e-14,
      -0.556464869058991e-21]
s0 = -0.332733756492168e4                       # IAPWS-95 reference
t1 = 0.368017112855051e-1 + 0.510878114959572e-1j
r1 = 0.447050716285388e2 + 0.656876847463481e2j
t2 = 0.337315741065416 + 0.335449415919309j
r2k = [-0.725974574329220e2 - 0.781008427112870e2j, -0.557107698030123e-4 + 0.464578634580806e-4j,
       0.234801409215913e-10 - 0.285651142904972e-10j]
pi0 = p0 / pt


def props(T: float, p: float) -> dict:
    tau, pi = T / Tt, p / pt
    dpi = pi - pi0
    r2 = sum(r * dpi ** k for k, r in enumerate(r2k))
    r2_p = sum(k * r * dpi ** (k - 1) for k, r in enumerate(r2k) if k) / pt
    r2_pp = 2 * r2k[2] / pt ** 2
    g0_p = sum(k * g * dpi ** (k - 1) for k, g in enumerate(g0) if k) / pt
    g0_pp = sum(k * (k - 1) * g * dpi ** (k - 2) for k, g in enumerate(g0) if k > 1) / pt ** 2
    L = lambda t: (t - tau) * np.log(t - tau) + (t + tau) * np.log(t + tau) - 2 * t * np.log(t) - tau ** 2 / t
    dL_dtau = lambda t: -np.log(t - tau) + np.log(t + tau) - 2 * tau / t
    g_TT = (1 / Tt) * np.real(sum(r * (1 / (t - tau) + 1 / (t + tau) - 2 / t) for r, t in ((r1, t1), (r2, t2))))
    g_p = g0_p + Tt * np.real(r2_p * L(t2))
    g_Tp = np.real(r2_p * dL_dtau(t2))
    g_pp = g0_pp + Tt * np.real(r2_pp * L(t2))
    return {"cp": -T * g_TT, "rho": 1 / g_p, "alpha": g_Tp / g_p, "kappa_T": -g_pp / g_p}


if __name__ == "__main__":
    # Table 6 check values
    checks = [(273.16, 611.657, 0.209678431622e4, 0.916709492200e3, 0.159863102566e-3, 0.117793449348e-9),
              (273.152519, 101325.0, 0.209671391024e4, 0.916721463419e3, 0.159841589458e-3, 0.117785291765e-9),
              (100.0, 100e6, 0.866333195517e3, 0.941678203297e3, 0.258495528207e-4, None)]
    for T, p, cp, rho, a, k in checks:
        r = props(T, p)
        print(f"T={T} p={p:g}: cp {r['cp']:.9g} (ref {cp:.9g})  rho {r['rho']:.9g} (ref {rho:.9g})  "
              f"alpha {r['alpha']:.9g} (ref {a:.9g})" + (f"  kappa {r['kappa_T']:.9g} (ref {k:.9g})" if k else ""))
    # cubic fit of cp(T) at 1 atm, 200–273.15 K, for the database ("poly" type, J/(kg K))
    T = np.linspace(200, 273.15, 200)
    cp = np.array([props(t, p0)["cp"] for t in T])
    c = np.polyfit(T, cp, 3)[::-1]
    err = np.max(np.abs(np.polyval(c[::-1], T) - cp) / cp) * 100
    print("cp fit coeffs (ascending, J/(kg K)):", [float(f"{v:.8g}") for v in c], f"max dev {err:.4f} %")
    r = props(273.15, p0)
    print(f"at 273.15 K, 1 atm: rho {r['rho']:.2f}, alpha {r['alpha']:.4e}, kappa_T {r['kappa_T']:.4e}")
    print("Fukusako (1990) linear form 185 + 6.89 T at 273.15 K:", 185 + 6.89 * 273.15, " IAPWS:", r["cp"])
