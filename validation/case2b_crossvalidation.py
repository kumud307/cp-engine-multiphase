"""Validation case 2b – HELD-OUT tests of the fitted NRTL parameters (cross-validation).

A model is only validated by data it was NOT fitted to. The stored parameter sets were fitted to all
available H^E data (and, for water + ethanol, to freezing points), so here each fit is REPEATED with
part of the data left out, and the left-out part is predicted:

  A  water + ethanol, leave-one-isotherm-out: refit (H^E of the other isotherms + freezing points),
     predict the H^E of the held-out isotherm.
  B  water + ethanol Cp^E: refit WITHOUT the Nagamachi & Francesconi paper, then predict the mean Cp^E
     over 298–323 K obtained from that paper's two isotherms (a genuine prediction).
  C  freezing points, 3-fold: refit with 2/3 of the freezing points, predict the remaining 1/3.
  D  water + methanol, leave-one-temperature-band-out (H^E only).
  E  ethanol + toluene (fitted to VLE + Cp^E): (i) a fit to the VLE data ALONE predicts Cp^E, (ii) leave one
     pressure of the VLE data out, (iii) fit with every second Cp^E composition and predict the others.

The same α and τ(T) form as the stored set are used for every refit. Results are cached in
validation/results/crossvalidation.csv (each refit takes ~10 s); run with --recompute to redo them.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from common import RESULTS, load, row
from cp_engine import CpEngine
from cp_engine.models import phase_change as pc
from cp_engine.models.nrtl import NRTL, fit_nrtl_joint

CACHE = RESULTS / "crossvalidation.csv"
NAGAMACHI = "10.1016/j.jct.2005.06.018"


def he_data(d, alc):
    g = d[(((d.comp1 == alc) & (d.comp2 == "water")) | ((d.comp1 == "water") & (d.comp2 == alc))) & (d.property == "HE")].copy()
    g["x"] = np.where(g.comp1 == alc, g.x1, 1 - g.x1)
    return g[(g.x > 0.01) & (g.x < 0.99) & (g.T_K <= 383.2)]


def freezing(db):
    fp = load("ethanol_water_freezing_points.csv")
    fp = fp[(fp.branch == "ice") & (fp.w_ethanol_percent > 0) & (fp.w_ethanol_percent <= 56.1)].copy()
    nw = (100 - fp.w_ethanol_percent) / db.get("water").MW
    ne = fp.w_ethanol_percent / db.get("ethanol").MW
    fp["x"] = ne / (nw + ne)
    return fp.reset_index(drop=True)


def gamma_inf(alc):
    """Same gamma_inf data as the stored water + ethanol fit (tools/fit_nrtl_he.py)."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
    import fit_nrtl_he as F
    if alc != "ethanol" or not F.INCLUDE_GAMMA_INF:
        return None
    g = load("thermoml_gamma_inf.csv")
    g = g[(g.solute == alc) & (g.solvent == "water")]
    return (np.full(len(g), 1e-6), g.T_K.values, np.log(g.gamma_inf.values), np.zeros(len(g)))


def refit(db, alc, he, fp=None):
    p = db.nrtl_pair(alc, "water")                      # stored set: same alpha and tau(T) form
    form = ("a12", "a21", "b12", "b21") + (("e12", "e21") if p.get("e_ij") or p.get("e_ji") else ())
    init = {"a12": p["a_ij"], "a21": p["a_ji"], "b12": p["b_ij"], "b21": p["b_ji"],
            "e12": p.get("e_ij", 0.0), "e21": p.get("e_ji", 0.0)}
    sle = None
    if fp is not None:
        w = db.get("water")
        dCp = w.cp_corr["liquid"].raw(w.Tm) - w.cp_corr["solid"].raw(w.Tm)
        sle = (fp.x.values, fp.freezing_point_K.values, w.Tm, w.dHfus, dCp, 1)
    gi = gamma_inf(alc)
    return fit_nrtl_joint((he.x, he.T_K, he.value), sle, p["c"], form, init, lng=gi, sigma_lng=0.1)["model"]


def rmsd(m, g):
    return float(np.sqrt(np.mean([(m.hE([x, 1 - x], T) - v) ** 2 for x, T, v in zip(g.x, g.T_K, g.value)])))


def cpe_from_isotherms(g):
    """Model-free mean Cp^E (298–323 K) from two measured isotherms (Redlich–Kister fit of each)."""
    A = {}
    for T in (298.15, 323.15):
        s = g[(g.T_K - T).abs() < 0.1]
        X = np.column_stack([s.x * (1 - s.x) * (2 * s.x - 1) ** k for k in range(4)])
        A[T] = np.linalg.lstsq(X, s.value, rcond=None)[0]
    xs = np.linspace(0.05, 0.95, 19)
    h = lambda T: xs * (1 - xs) * sum(a * (2 * xs - 1) ** k for k, a in enumerate(A[T]))
    return xs, (h(323.15) - h(298.15)) / 25.0


def compute() -> pd.DataFrame:
    db = CpEngine().db
    d = load("thermoml_binary_data.csv")
    out = []
    # A: ethanol leave-one-isotherm-out
    he, fp = he_data(d, "ethanol"), freezing(db)
    vle = NRTL(["ethanol", "water"], db, prefer="VLE")
    for T in sorted(he.T_K.round(0).unique()):
        test, train = he[he.T_K.round(0) == T], he[he.T_K.round(0) != T]
        m = refit(db, "ethanol", train, fp)
        out.append({"test": "A", "system": "ethanol + water", "held_out": f"H^E isotherm {T:.0f} K", "N": len(test),
                    "model": rmsd(m, test), "vle": rmsd(vle, test), "metric": "RMSD J/mol"})
        print(f"  A  ethanol, held-out {T:.0f} K: RMSD {out[-1]['model']:.1f} (VLE params {out[-1]['vle']:.1f}) J/mol")
    # B: ethanol Cp^E without the Nagamachi paper
    m = refit(db, "ethanol", he[he.doi != NAGAMACHI], fp)
    xs, data = cpe_from_isotherms(he[he.doi == NAGAMACHI])
    pred = np.array([(m.hE([x, 1 - x], 323.15) - m.hE([x, 1 - x], 298.15)) / 25.0 for x in xs])
    vpred = np.array([(vle.hE([x, 1 - x], 323.15) - vle.hE([x, 1 - x], 298.15)) / 25.0 for x in xs])
    out.append({"test": "B", "system": "ethanol + water", "held_out": "mean Cp^E 298-323 K (Nagamachi paper excluded from fit)",
                "N": len(xs), "model": float(np.sqrt(np.mean((pred - data) ** 2))),
                "vle": float(np.sqrt(np.mean((vpred - data) ** 2))), "metric": "RMSD J/(mol K)"})
    print(f"  B  ethanol Cp^E held-out: RMSD {out[-1]['model']:.2f} (VLE {out[-1]['vle']:.2f}) J/(mol K)")
    # C: freezing points 3-fold
    w = db.get("water")
    dCp = w.cp_corr["liquid"].raw(w.Tm) - w.cp_corr["solid"].raw(w.Tm)
    devs, failed = [], []
    for k in range(3):
        test, train = fp[fp.index % 3 == k], fp[fp.index % 3 != k]
        m = refit(db, "ethanol", he, train)
        for x, Tf, wt in zip(test.x, test.freezing_point_K, test.w_ethanol_percent):
            try:
                T = pc.liquidus_temperature(1 - x, w.Tm, w.dHfus, lambda t: m.ln_gamma([x, 1 - x], t)[1], dCp)
                devs.append(abs(T - Tf))
            except ValueError:                  # the held-out model predicts no ice at this composition
                failed.append(wt)
    out.append({"test": "C", "system": "ice + water + ethanol", "held_out": "freezing points, 3-fold"
                + (f" ({len(failed)} point(s) at {', '.join(f'{v:g}' for v in failed)} wt% without a liquidus)" if failed else ""),
                "N": len(devs), "model": float(np.mean(devs)), "vle": np.nan, "metric": "|ΔT| K",
                "failed": len(failed)})
    print(f"  C  freezing points held-out: mean |dT| {out[-1]['model']:.2f} K over {len(devs)} points; "
          f"{len(failed)} without a liquidus ({failed})")
    # D: methanol leave-one-temperature-band-out
    he = he_data(d, "methanol")
    vle = NRTL(["methanol", "water"], db, prefer="VLE")
    bands = [(290, 310), (320, 335), (335, 350), (350, 365), (365, 385)]
    for lo, hi in bands:
        test, train = he[(he.T_K >= lo) & (he.T_K < hi)], he[(he.T_K < lo) | (he.T_K >= hi)]
        if len(test) == 0:
            continue
        m = refit(db, "methanol", train)
        out.append({"test": "D", "system": "methanol + water", "held_out": f"H^E {lo}-{hi} K", "N": len(test),
                    "model": rmsd(m, test), "vle": rmsd(vle, test), "metric": "RMSD J/mol"})
        print(f"  D  methanol, held-out {lo}-{hi} K: RMSD {out[-1]['model']:.1f} (VLE params {out[-1]['vle']:.1f}) J/mol")
    # E: ethanol + toluene
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
    import fit_nrtl_he as F
    lng, cpe, v, cp = F.load_ethanol_toluene(db)
    stored = db.nrtl_pair("ethanol", "toluene")
    form = ("a12", "a21", "b12", "b21")
    init = {"a12": stored["a_ij"], "a21": stored["a_ji"], "b12": stored["b_ij"], "b21": stored["b_ji"]}
    cpe_rmsd = lambda m, x, T, c: float(np.sqrt(np.mean([(m.cpE([a, 1 - a], t) - y) ** 2 for a, t, y in zip(x, T, c)])))
    lng_rms = lambda m, t: float(np.sqrt(np.mean([((m.ln_gamma([x, 1 - x], T)[0] - a) ** 2 + (m.ln_gamma([x, 1 - x], T)[1] - b) ** 2) / 2
                                                  for x, T, a, b in zip(*t)])))
    only_vle = fit_nrtl_joint(None, None, stored["c"], form, init, lng=lng)["model"]
    out.append({"test": "E", "system": "ethanol + toluene", "held_out": "Cp^E predicted by a fit to VLE only",
                "N": len(cpe[0]), "model": cpe_rmsd(only_vle, *cpe), "vle": np.nan, "metric": "RMSD J/(mol K)"})
    print(f"  E  ethanol + toluene, Cp^E from a VLE-only fit: RMSD {out[-1]['model']:.2f} J/(mol K)")
    for P in sorted(v.P_kPa.unique()):
        sel = v.P_kPa.values == P
        tr = tuple(a[~sel] for a in lng)
        te = tuple(a[sel] for a in lng)
        m = fit_nrtl_joint(None, None, stored["c"], form, init, cpe=cpe, lng=tr)["model"]
        out.append({"test": "F", "system": "ethanol + toluene", "held_out": f"ln gamma at {P:g} kPa",
                    "N": int(sel.sum()), "model": lng_rms(m, te), "vle": np.nan, "metric": "RMS ln gamma"})
        print(f"  F  ethanol + toluene, held-out VLE at {P:g} kPa: RMS ln gamma {out[-1]['model']:.3f}")
    idx = np.arange(len(cpe[0]))
    for k in (0, 1):
        te = tuple(a[idx % 2 == k] for a in cpe)
        tr = tuple(a[idx % 2 != k] for a in cpe)
        m = fit_nrtl_joint(None, None, stored["c"], form, init, cpe=tr, lng=lng)["model"]
        out.append({"test": "G", "system": "ethanol + toluene", "held_out": f"Cp^E, compositions {k}::2",
                    "N": len(te[0]), "model": cpe_rmsd(m, *te), "vle": np.nan, "metric": "RMSD J/(mol K)"})
        print(f"  G  ethanol + toluene, held-out Cp^E half {k}: RMSD {out[-1]['model']:.2f} J/(mol K)")
    df = pd.DataFrame(out)
    RESULTS.mkdir(exist_ok=True)
    df.to_csv(CACHE, index=False)
    return df


def run(save: bool = True, recompute: bool = False) -> list[dict]:
    df = compute() if (recompute or not CACHE.exists()) else pd.read_csv(CACHE)
    rows = []
    for (test, system, metric), g in df.groupby(["test", "system", "metric"], sort=True):
        n = int(g.N.sum())
        if metric.startswith("RMSD"):
            mean = float(np.sqrt(np.sum(g.model ** 2 * g.N) / n))
        else:
            mean = float(np.sum(g.model * g.N) / n)
        label = {"A": "H^E, leave-one-isotherm-out (held-out)", "B": "mean Cp^E 298-323 K, source paper held out",
                 "C": "freezing points, 3-fold cross-validation (held-out)",
                 "D": "H^E, leave-one-temperature-band-out (held-out)",
                 "E": "Cp^E predicted from VLE data alone (no Cp data in the fit)",
                 "F": "VLE activity coefficients, leave-one-pressure-out (held-out)",
                 "G": "Cp^E, every second composition held out"}[test]
        fails = int(g["failed"].fillna(0).sum()) if "failed" in g else 0
        note = f" ({fails} held-out point(s) with no predicted liquidus, excluded)" if fails else ""
        tag = "" if test == "E" else ", NRTL refitted"
        rows.append(row("2b", system, f"{label}{tag}{note}", n, mean, float(g.model.max()),
                        "cross-validation", metric=metric))
        if g.vle.notna().any():
            vmean = float(np.sqrt(np.sum(g.vle ** 2 * g.N) / n)) if metric.startswith("RMSD") else float(np.mean(g.vle))
            rows.append(row("2b", system, f"{label}, VLE parameters (for comparison)", n, vmean, float(g.vle.max()),
                            "cross-validation", metric=metric))
    return rows


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    print(pd.DataFrame(run(recompute="--recompute" in sys.argv)).to_string(index=False))
