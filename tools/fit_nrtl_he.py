"""Refit NRTL parameters to measured excess enthalpies (H^E) and store them in data/binary_nrtl.json.

Why: Cp^E = (∂H^E/∂T)_P,x, so a model used for heat capacities must reproduce H^E and its temperature
dependence. VLE-regressed parameters only fix G^E near the boiling range and give H^E errors of
140-260 J/mol for water + alcohols; refitting to H^E reduces this to ~25 J/mol.

Data: validation/data/thermoml_binary_data.csv (NIST ThermoML archive; see tools/extract_thermoml.py).
Selection: 0.01 < x < 0.99 (excludes infinite-dilution titration points), T ≤ 383.2 K.

Candidates: τ_ij = a_ij + b_ij/T (+ e_ij ln T), with α ∈ {0.2, 0.3, 0.4, 0.47}. A candidate is accepted only
if the binary stays thermodynamically stable (no false liquid-liquid split) from 240 K to 390 K – the
range used by the ice-slurry calculations. The accepted candidate with the lowest RMSD is stored.

Run:  python tools/fit_nrtl_he.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from cp_engine.database import Database  # noqa: E402
from cp_engine.models.nrtl import NRTL, fit_nrtl_to_excess  # noqa: E402

PAIRS = ["ethanol", "methanol"]          # each with water
T_MAX = 383.2
FORMS = [("a12", "a21", "b12", "b21"), ("a12", "a21", "b12", "b21", "e12", "e21")]
ALPHAS = (0.2, 0.3, 0.4, 0.47)


def load_he(db_dir: Path, comp: str) -> pd.DataFrame:
    d = pd.read_csv(db_dir / "validation" / "data" / "thermoml_binary_data.csv", comment="#")
    g = d[((d.comp1 == comp) & (d.comp2 == "water")) | ((d.comp1 == "water") & (d.comp2 == comp))]
    g = g[(g.property == "HE") & (g.T_K <= T_MAX)].copy()
    g["x"] = np.where(g.comp1 == comp, g.x1, 1 - g.x1)
    return g[(g.x > 0.01) & (g.x < 0.99)]


def to_model(f: dict) -> NRTL:
    mats = {s: [[0, f[f"{s}_ij"]], [f[f"{s}_ji"], 0]] for s in "abef"}
    mats["c"] = np.full((2, 2), f["c"])
    return NRTL(["1", "2"], params=mats)


def main() -> None:
    db = Database(ROOT / "data")
    path = ROOT / "data" / "binary_nrtl.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    data["nrtl"] = [p for p in data["nrtl"] if p.get("fitted_to") != "HE"]     # replace old H^E fits
    for comp in PAIRS:
        he = load_he(ROOT, comp)
        vle = db.nrtl_pair(comp, "water", prefer="VLE")
        init = {"a12": vle["a_ij"], "a21": vle["a_ji"], "b12": vle["b_ij"], "b21": vle["b_ji"]}
        vle_rmsd = float(np.sqrt(np.mean([(NRTL([comp, "water"], db, prefer="VLE").hE([x, 1 - x], T) - v) ** 2
                                          for x, T, v in zip(he.x, he.T_K, he.value)])))
        best = None
        for alpha in ALPHAS:
            for form in FORMS:
                f = fit_nrtl_to_excess(he.x, he.T_K, he.value, "HE", alpha=alpha, fit=form, initial=init)
                m = to_model(f)
                if all(m.is_stable_binary(T, 1000) for T in range(240, 391, 10)):
                    if best is None or f["rmse"] < best["rmse"]:
                        best = f
        dois = sorted(he.doi.unique())
        entry = {"i": comp, "j": "water",
                 **{k: round(float(best[k]), 6) for k in ("a_ij", "a_ji", "b_ij", "b_ji", "e_ij", "e_ji", "f_ij", "f_ji")},
                 "c": best["c"], "d": 0.0, "Tmin": float(he.T_K.min()), "Tmax": float(he.T_K.max()),
                 "fitted_to": "HE",
                 "rmse_HE_J_per_mol": round(best["rmse"], 1), "rmse_HE_VLE_parameters": round(vle_rmsd, 1),
                 "n_points": int(len(he)),
                 "source": "Fitted (tools/fit_nrtl_he.py) to measured H^E from the NIST ThermoML archive: " + "; ".join(dois)
                           + ". Stable (no liquid-liquid split) from 240 to 390 K.",
                 "verify": False}
        data["nrtl"].append(entry)
        print(f"{comp}+water: {len(he)} points, T {he.T_K.min():.1f}-{he.T_K.max():.1f} K; "
              f"RMSD H^E: VLE params {vle_rmsd:.1f} J/mol -> refit {best['rmse']:.1f} J/mol (alpha={best['c']})")
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
