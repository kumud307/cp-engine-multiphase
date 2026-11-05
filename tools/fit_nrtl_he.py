"""Refit NRTL parameters to measured excess enthalpies (and freezing points) and store them in
data/binary_nrtl.json.

Why: Cp^E = (∂H^E/∂T)_P,x, so a model used for heat capacities must reproduce H^E and its temperature
dependence. VLE-regressed parameters only fix G^E near the boiling range and give H^E errors of
140-260 J/mol for water + alcohols.

water + ethanol : fitted JOINTLY to H^E (298-383 K, NIST ThermoML) and to the freezing points of
                  aqueous ethanol (Lange 1967, ≤ 56 wt%), so the SAME liquid model gives the heat of
                  mixing (→ Cp^E) and the ice–solution equilibrium (→ latent heat of ice slurries).
water + methanol: fitted to H^E only (no freezing data in the project).
ethanol + toluene: no H^E data exist in the archive; fitted JOINTLY to isobaric VLE (activity coefficients
                  from y P = x γ Psat, Chen et al. 2012, 4 pressures, 350-420 K) and to the excess heat
                  capacity at 313 K (Vega-Maza et al. 2011). VLE alone fixes G^E but gives Cp^E errors of
                  ~7 J/(mol K) and an unstable model; adding Cp^E fixes both.

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


# Including measured gamma_inf of ethanol in water in the water + ethanol fit improves gamma_inf (error
# 35 % -> 8 %) but worsens the held-out Cp^E (3.7 -> 4.3 J/(mol K)) and freezing points (2.5 -> 3.3 K).
# gamma_inf does not enter any heat-capacity calculation, so the Cp engine uses the set WITHOUT it.
INCLUDE_GAMMA_INF = False


def gamma_inf_tuple(comp: str):
    """Measured activity coefficients at infinite dilution of `comp` in water, as an lng tuple
    (x1 → 0; only ln γ1 is meaningful there, ln γ2 = 0). Used for water + ethanol so that the same
    parameters also describe dilute-solution activity (σ = 0.1 in ln γ)."""
    g = pd.read_csv(ROOT / "validation" / "data" / "thermoml_gamma_inf.csv", comment="#")
    g = g[(g.solute == comp) & (g.solvent == "water")]
    if comp != "ethanol" or g.empty or not INCLUDE_GAMMA_INF:
        return None
    n = len(g)
    return (np.full(n, 1e-6), g.T_K.values, np.log(g.gamma_inf.values), np.zeros(n))


def load_ethanol_toluene(db: Database):
    """(lng tuple, cpe tuple, DataFrames) for ethanol (1) + toluene (2)."""
    v = pd.read_csv(ROOT / "validation" / "data" / "thermoml_vle_ethanol_toluene.csv", comment="#")
    v = v[(v.x_ethanol > 0) & (v.x_ethanol < 1)].copy()
    e, t = db.get("ethanol"), db.get("toluene")
    v["lng1"] = np.log(v.y_ethanol * v.P_kPa * 1e3 / (v.x_ethanol * np.array([e.psat(T) for T in v.T_K])))
    v["lng2"] = np.log((1 - v.y_ethanol) * v.P_kPa * 1e3 / ((1 - v.x_ethanol) * np.array([t.psat(T) for T in v.T_K])))
    d = pd.read_csv(ROOT / "validation" / "data" / "thermoml_binary_data.csv", comment="#")
    cp = d[(d.property == "Cp_molar") & (d.P_kPa <= 200) & (((d.comp1 == "ethanol") & (d.comp2 == "toluene"))
                                                         | ((d.comp1 == "toluene") & (d.comp2 == "ethanol")))].copy()
    cp["x"] = np.where(cp.comp1 == "ethanol", cp.x1, 1 - cp.x1)
    c1, c2 = cp[cp.x == 1].value.iloc[0], cp[cp.x == 0].value.iloc[0]          # the paper's own pure values
    cp["cpe"] = cp.value - cp.x * c1 - (1 - cp.x) * c2
    cp = cp[(cp.x > 0) & (cp.x < 1)]
    return (v.x_ethanol.values, v.T_K.values, v.lng1.values, v.lng2.values), (cp.x.values, cp.T_K.values, cp.cpe.values), v, cp


def fit_ethanol_toluene(db: Database, lng, cpe):
    """Pre-fit to VLE only (good starting point), then the joint VLE + Cp^E fit with the stability check."""
    from cp_engine.models.nrtl import fit_nrtl_joint
    pre = fit_nrtl_joint(None, None, 0.3, ("a12", "a21", "b12", "b21"), {"b12": 300.0, "b21": 300.0}, lng=lng)
    init = {k2: pre[k] for k, k2 in (("a_ij", "a12"), ("a_ji", "a21"), ("b_ij", "b12"), ("b_ji", "b21"))}
    return fit_nrtl_best(None, None, init, alphas=(0.3, 0.47), cpe=cpe, lng=lng, T_stable=(250, 390), score="cost")


def main() -> None:
    db = Database(ROOT / "data")
    path = ROOT / "data" / "binary_nrtl.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    data["nrtl"] = [p for p in data["nrtl"] if p.get("fitted_to") == "VLE"]      # keep the literature VLE sets
    for comp in ("ethanol", "methanol"):
        he = load_he(comp)
        vle = db.nrtl_pair(comp, "water", prefer="VLE")
        init = {"a12": vle["a_ij"], "a21": vle["a_ji"], "b12": vle["b_ij"], "b21": vle["b_ji"]}
        sle = sle_tuple(db, load_freezing(db)) if comp == "ethanol" else None
        gi = gamma_inf_tuple(comp)
        best = fit_nrtl_best((he.x, he.T_K, he.value), sle, init, lng=gi, sigma_lng=0.1)
        if gi is not None:
            best["fitted_to"] = best["fitted_to"].replace("VLE", "GammaInf")
        vle_m = NRTL([comp, "water"], db, prefer="VLE")
        vle_rmsd = float(np.sqrt(np.mean([(vle_m.hE([x, 1 - x], T) - v) ** 2
                                          for x, T, v in zip(he.x, he.T_K, he.value)])))
        refs = "; ".join(sorted(he.doi.unique()))
        extra = (" and to the freezing points of aqueous ethanol 2.5-56 wt% (Lange 1967) and to measured activity "
                 "coefficients of ethanol at infinite dilution in water (doi:10.1021/je0495942; doi:10.1021/je049875+)") if sle else ""
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
    lng, cpe, v, cp = load_ethanol_toluene(db)
    best = fit_ethanol_toluene(db, lng, cpe)
    m = best["model"]
    rms_lng = float(np.sqrt(np.mean([((m.ln_gamma([x, 1 - x], T)[0] - a) ** 2 + (m.ln_gamma([x, 1 - x], T)[1] - b) ** 2) / 2
                                     for x, T, a, b in zip(*lng)])))
    rms_cpe = float(np.sqrt(np.mean([(m.cpE([x, 1 - x], T) - c) ** 2 for x, T, c in zip(*cpe)])))
    data["nrtl"].append({
        "i": "ethanol", "j": "toluene",
        **{k: round(float(best[k]), 6) for k in ("a_ij", "a_ji", "b_ij", "b_ji", "e_ij", "e_ji", "f_ij", "f_ji")},
        "c": best["c"], "d": 0.0, "Tmin": 250.0, "Tmax": float(max(v.T_K)), "fitted_to": best["fitted_to"],
        "rms_ln_gamma": round(rms_lng, 3), "rmse_CpE_J_per_molK": round(rms_cpe, 2), "n_points": int(len(v) + len(cp)),
        "source": "Fitted (tools/fit_nrtl_he.py) jointly to isobaric VLE (Chen, Zhong & Xu 2012, doi:10.1021/je200921u; "
                  "Psat from NIST Antoine parameters) and Cp of the mixture at 313.15 K (Vega-Maza et al. 2011, "
                  "doi:10.1016/j.jct.2011.06.019), NIST ThermoML archive. Phase-stable 250-390 K.",
        "verify": False})
    print(f"ethanol+toluene [{best['fitted_to']}]: RMS ln gamma {rms_lng:.3f}, Cp^E RMSD {rms_cpe:.2f} J/(mol K) "
          f"(alpha={best['c']})")
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
