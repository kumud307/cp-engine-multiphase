"""Refit NRTL parameters to measured excess enthalpies (and freezing points) and store them in
data/binary_nrtl.json.

Why: Cp^E = (∂H^E/∂T)_P,x, so a model used for heat capacities must reproduce H^E and its temperature
dependence. VLE-regressed parameters only fix G^E near the boiling range and give H^E errors of
140-260 J/mol for water + alcohols.

water + ethanol : fitted JOINTLY to H^E (298-383 K, NIST ThermoML) and to the freezing points of
                  aqueous ethanol (Lange 1967, ≤ 56 wt%), so the SAME liquid model gives the heat of
                  mixing (→ Cp^E) and the ice–solution equilibrium (→ latent heat of ice slurries).
water + methanol: fitted to H^E only (no freezing data in the project).

A candidate is accepted only if the binary stays phase-stable (no false liquid–liquid split) from
220 K to 390 K; the stable candidate with the lowest H^E error is stored.

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
from cp_engine.models.nrtl import NRTL, fit_nrtl_best  # noqa: E402

T_MAX = 383.2


def load_he(comp: str) -> pd.DataFrame:
    d = pd.read_csv(ROOT / "validation" / "data" / "thermoml_binary_data.csv", comment="#")
    g = d[((d.comp1 == comp) & (d.comp2 == "water")) | ((d.comp1 == "water") & (d.comp2 == comp))]
    g = g[(g.property == "HE") & (g.T_K <= T_MAX)].copy()
    g["x"] = np.where(g.comp1 == comp, g.x1, 1 - g.x1)
    return g[(g.x > 0.01) & (g.x < 0.99)]


def load_freezing(db: Database, max_wt: float = 56.1) -> pd.DataFrame:
    fp = pd.read_csv(ROOT / "validation" / "data" / "ethanol_water_freezing_points.csv", comment="#")
    fp = fp[(fp.branch == "ice") & (fp.w_ethanol_percent > 0) & (fp.w_ethanol_percent <= max_wt)].copy()
    nw = (100 - fp.w_ethanol_percent) / db.get("water").MW
    ne = fp.w_ethanol_percent / db.get("ethanol").MW
    fp["x"] = ne / (nw + ne)                       # mole fraction of ethanol (component 1)
    return fp


def sle_tuple(db: Database, fp: pd.DataFrame) -> tuple:
    w = db.get("water")
    dCp = w.cp_corr["liquid"].raw(w.Tm) - w.cp_corr["solid"].raw(w.Tm)
    return (fp.x.values, fp.freezing_point_K.values, w.Tm, w.dHfus, dCp, 1)   # ice = component index 1


def main() -> None:
    db = Database(ROOT / "data")
    path = ROOT / "data" / "binary_nrtl.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    data["nrtl"] = [p for p in data["nrtl"] if not str(p.get("fitted_to", "")).startswith("HE")]
    for comp in ("ethanol", "methanol"):
        he = load_he(comp)
        vle = db.nrtl_pair(comp, "water", prefer="VLE")
        init = {"a12": vle["a_ij"], "a21": vle["a_ji"], "b12": vle["b_ij"], "b21": vle["b_ji"]}
        sle = sle_tuple(db, load_freezing(db)) if comp == "ethanol" else None
        best = fit_nrtl_best((he.x, he.T_K, he.value), sle, init)
        vle_m = NRTL([comp, "water"], db, prefer="VLE")
        vle_rmsd = float(np.sqrt(np.mean([(vle_m.hE([x, 1 - x], T) - v) ** 2
                                          for x, T, v in zip(he.x, he.T_K, he.value)])))
        refs = "; ".join(sorted(he.doi.unique()))
        extra = " and to the freezing points of aqueous ethanol 2.5-56 wt% (Lange 1967)" if sle else ""
        entry = {"i": comp, "j": "water",
                 **{k: round(float(best[k]), 6) for k in ("a_ij", "a_ji", "b_ij", "b_ji", "e_ij", "e_ji", "f_ij", "f_ji")},
                 "c": best["c"], "d": 0.0, "Tmin": 220.0 if sle else float(he.T_K.min()), "Tmax": float(he.T_K.max()),
                 "fitted_to": best["fitted_to"],
                 "rmse_HE_J_per_mol": round(best["rmse"], 1), "rmse_HE_VLE_parameters": round(vle_rmsd, 1),
                 "n_points": int(len(he)) + (len(sle[0]) if sle else 0),
                 "source": "Fitted (tools/fit_nrtl_he.py) to measured H^E from the NIST ThermoML archive: " + refs
                           + extra + ". Phase-stable (no liquid-liquid split) from 220 to 390 K.",
                 "verify": False}
        data["nrtl"].append(entry)
        print(f"{comp}+water [{best['fitted_to']}]: H^E RMSD VLE params {vle_rmsd:.1f} -> {best['rmse']:.1f} J/mol "
              f"(alpha={best['c']})")
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
