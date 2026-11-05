"""Validation case 1 – pure components vs NIST reference data (sanity check).

Data: NIST Chemistry WebBook SRD 69 isobaric tables at 0.101325 MPa
(downloaded 2026-11-02, validation/data/nist_<name>_1atm.csv) and single
measured values quoted on the NIST condensed-phase pages
(validation/data/pure_component_points.csv; also O'Hanley et al. 2012 Table 5 for glycerol).

Checks Cp(T) of the liquid correlations, and the Cv → Cp route
(cp = cv + T v α²/κ_T) using NIST cv for water.
"""
from __future__ import annotations

import numpy as np

from common import FIGS, aard, load, row
from cp_engine import CpEngine
from cp_engine.plotting import plot_validation

FLUIDS = ["water", "methanol", "n-hexane", "toluene", "n-dodecane"]


def run(save: bool = True) -> list[dict]:
    engine = CpEngine()
    db = engine.db
    rows = []
    for name in FLUIDS:
        d = load(f"nist_{name}_1atm.csv")
        c = db.get(name)
        corr = c.cp_corr["liquid"]
        d = d[(d.T_K >= corr.Tmin - 0.05) & (d.T_K <= corr.Tmax)]
        model = [engine.solve({"T": t, "components": {name: "liquid"}}).cp_mass for t in d.T_K]
        a, mx, n = aard(model, d.cp_J_kgK)
        rows.append(row("1", name, "cp (liquid, 1 atm)", n, a, mx, "NIST SRD 69 (Lemmon et al.)"))
        if save:
            T = np.linspace(d.T_K.min(), d.T_K.max(), 200)
            plot_validation(T, [c.cp_mass("liquid", t) for t in T], d.T_K, d.cp_J_kgK,
                            "c$_p$ / J kg$^{-1}$ K$^{-1}$", f"{c.name} (liquid, 0.101 MPa)",
                            label_exp="NIST SRD 69", name=f"val1_{name}", outdir=FIGS)

    pts = load("pure_component_points.csv")
    for (comp, phase), g in pts.groupby(["component", "phase"]):
        c = db.get(comp)
        model = [c.cp(phase, t) for t in g.T_K]
        a, mx, n = aard(model, g.Cp_J_molK)
        ref = "NIST WebBook" + (" + O'Hanley et al. 2012" if g.reference.str.contains("O'Hanley").any() else "")
        rows.append(row("1", f"{comp} ({phase})", "Cp (point values)", n, a, mx, ref))

    # pressure dependence: NIST isotherms at 298.15 K, 0.1-100 MPa
    if save:
        from cp_engine.plotting import plot_data_vs_models
        groups = {}
        for name in ("water", "methanol", "n-hexane", "toluene"):
            iso = load(f"nist_{name}_298K_isotherm.csv")
            P = np.linspace(0.1, 100.0, 41)
            cp0 = db.get(name).cp_mass("liquid", 298.15)
            tait = [engine.solve({"T": 298.15, "P": (p, "MPa"), "components": {name: "liquid"}}).cp_mass / cp0 for p in P]
            groups[name] = {"data": (iso.P_MPa, iso.cp_J_kgK / cp0),
                            "models": {"Tait pressure correction": (P, tait), "no correction": (P, [1.0] * len(P))}}
        plot_data_vs_models(groups, "P / MPa", r"c$_p$(P) / c$_{p,\mathrm{model}}$(0.1 MPa)", "Pressure dependence at 298.15 K (NIST)",
                            name="fig_pressure", outdir=FIGS, w=4.8, legend_outside=True)
    for name in FLUIDS:
        iso = load(f"nist_{name}_298K_isotherm.csv")
        pred = [engine.solve({"T": t, "P": (P, "MPa"), "components": {name: "liquid"}}).cp_mass
                for t, P in zip(iso.T_K, iso.P_MPa)]
        a, mx, n = aard(pred, iso.cp_J_kgK)
        rows.append(row("1", name, "cp vs pressure 0.1-100 MPa (Tait pressure correction)", n, a, mx, "NIST SRD 69"))
        p0 = engine.solve({"T": 298.15, "components": {name: "liquid"}}).cp_mass
        a0, mx0, _ = aard([p0] * len(iso), iso.cp_J_kgK)
        rows.append(row("1", name, "cp vs pressure 0.1-100 MPa, without pressure correction (for comparison)",
                        n, a0, mx0, "NIST SRD 69"))

    # ideal-gas Cp vs NIST recommended values (200-1000 K)
    gas = load("nist_ideal_gas_cp.csv")
    for comp, g in gas.groupby("component"):
        corr = db.get(comp).cp_corr["gas"]
        a, mx, n = aard([corr.raw(t) for t in g.T_K], g.Cp_J_molK)
        rows.append(row("1", f"{comp} (ideal gas)", "Cp 200-1000 K", n, a, mx, "NIST WebBook (TRC / Chao / Scott / Draeger)"))

    # Cv -> Cp route for water with NIST cv
    d = load("nist_water_1atm.csv")
    cp_from_cv = [engine.solve({"T": t, "Cv": (cv, "J/(kg K)"), "components": ["water"]})
                  for t, cv in zip(d.T_K, d.cv_J_kgK)]
    vals = [[r for r in res.routes if r.tag == "Cv"][0].cp_mass for res in cp_from_cv]
    a, mx, n = aard(vals, d.cp_J_kgK)
    rows.append(row("1", "water", "cp from NIST cv via cp = cv + Tvα²/κ_T", n, a, mx, "NIST SRD 69"))
    return rows


if __name__ == "__main__":
    import pandas as pd
    print(pd.DataFrame(run()).to_string(index=False))
