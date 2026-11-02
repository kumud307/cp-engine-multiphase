"""Validation case 4 – phase change: freezing (liquidus) points and apparent-Cp consistency.

Experimental: freezing points of water + ethanol, 2.5-71.9 wt% ethanol (Lange's Handbook, 1967),
validation/data/ethanol_water_freezing_points.csv. The freezing point is exactly where the engine's
ice slurry loses its last ice, so it tests the solid–liquid-equilibrium model that drives the latent
part of the apparent heat capacity of ice slurries.

Models compared (ice = pure solid water; liquid mole fraction x_w):
  ideal                 ln x_w = −(ΔH_fus/R)(1/T − 1/T_m)
  ideal + ΔCp           adds (ΔCp/R)(T_m/T − 1 − ln(T_m/T)), consistent with the Kirchhoff enthalpies
  NRTL γ_w, H^E+SLE set ln(x_w γ_w) = ... with the parameters fitted JOINTLY to H^E and these freezing points
                        (engine default: one liquid model for H^E and SLE). Because these points were
                        used in the fit, its held-out accuracy is reported by case 2b (3-fold CV).
  NRTL γ_w, VLE set     same with the VLE parameters
Internal checks (no data needed): the engine's liquidus equals the model equation; ∫cp_app dT
through a melting peak equals ΔH (energy conservation) for an ice slurry and a PCM slurry.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from common import FIGS, aard, load, row
from cp_engine import CpEngine
from cp_engine.models import phase_change as pc
from cp_engine.models.nrtl import NRTL
from cp_engine.plotting import plot_apparent_peak, plot_cp_vs_T, plot_data_vs_models
from cp_engine.system import SystemModel

CSV = "ethanol_water_freezing_points.csv"
REF = "Lange's Handbook (1967) via Wikipedia 'Ethanol (data page)'"
SYSTEMS = {
    "ice_ethanol_10wt": {"mass_fractions": {"water": 0.9, "ethanol": 0.1}, "T": (250.0, 272.0),
                         "title": "Ice slurry, 10 wt% ethanol(aq)"},
    "pcm_octadecane_20wt": {"mass_fractions": {"water": 0.8, "n-octadecane": 0.2}, "T": (292.0, 312.0),
                            "title": "PCM slurry, 20 wt% n-octadecane in water"},
}


def engine_liquidus(db, w_eth: float) -> float:
    """Temperature at which the engine's ice slurry has just lost its last ice."""
    n = {"water": (1 - w_eth) / db.get("water").MW, "ethanol": w_eth / db.get("ethanol").MW}
    s = SystemModel(db, n, {})
    s.configure(200.0, 273.15)

    def ice(T):
        return sum(m for name, k, st, m in s.split(T) if st == "solid") - 1e-9 * n["water"]
    return brentq(ice, 200.0, 273.1499)


def run(save: bool = True) -> list[dict]:
    engine = CpEngine()
    db = engine.db
    w = db.get("water")
    dCp = w.cp_corr["liquid"].raw(w.Tm) - w.cp_corr["solid"].raw(w.Tm)
    data = load(CSV)
    data = data[(data.branch == "ice") & (data.w_ethanol_percent > 0)].copy()
    nw = (100 - data.w_ethanol_percent) / w.MW
    ne = data.w_ethanol_percent / db.get("ethanol").MW
    data["x_w"] = nw / (nw + ne)
    DEFAULT = "NRTL γ, H^E+SLE set (engine default; fit quality - held-out: case 2b)"
    models = {"ideal": lambda x: pc.ideal_freezing_point(x, w.Tm, w.dHfus),
              "ideal + ΔCp (option sle_activity='ideal')": lambda x: pc.ideal_freezing_point(x, w.Tm, w.dHfus, dCp)}
    for label, prefer in ((DEFAULT, "HE"), ("NRTL γ, VLE set", "VLE")):
        m = NRTL(["water", "ethanol"], db, prefer=prefer)
        models[label] = (lambda m: lambda x: pc.liquidus_temperature(
            x, w.Tm, w.dHfus, lambda T: m.ln_gamma([x, 1 - x], T)[0], dCp))(m)

    def safe(f, x):
        try:
            return f(x)
        except ValueError:
            return np.nan

    rows = []
    for label, f in models.items():
        data[label] = [safe(f, x) for x in data.x_w]
        for lim, tag in ((30.0, "≤ 30 wt% (ice-slurry range)"), (100.0, "all, 2.5-72 wt%")):
            s = data[(data.w_ethanol_percent <= lim) & data[label].notna()]
            dev = (s[label] - s.freezing_point_K).abs()
            rows.append(row("4", "ice + water + ethanol", f"freezing point, {label}, {tag}",
                            len(s), float(dev.mean()), float(dev.max()), REF, metric="|ΔT| K"))
    # engine consistency: the full engine's liquidus equals the ideal+ΔCp equation
    test = data[data.w_ethanol_percent <= 30.0]
    eng = [engine_liquidus(db, wp / 100) for wp in test.w_ethanol_percent]
    err = np.abs(np.array(eng) - test[DEFAULT].values).max()
    rows.append(row("4", "ice + water + ethanol", "engine liquidus vs SLE equation", len(test),
                    float(err), float(err), "internal check", metric="|ΔT| K"))

    # energy conservation through phase change
    for key, spec in SYSTEMS.items():
        T1, T2 = spec["T"]
        inp = {"mass_fractions": spec["mass_fractions"], "T1": T1, "T2": T2}
        T = np.linspace(T1, T2, 1601)
        curve = engine.cp_curve(inp, T)
        res = engine.solve(dict(inp, total_mass=1.0))
        mean = [r for r in res.routes if r.tag == "A:mean"][0]
        integral = np.trapezoid(curve.total, T)
        dev = abs(integral - mean.cp_mass * (T2 - T1)) / (mean.cp_mass * (T2 - T1)) * 100
        rows.append(row("4", key, "∫cp_app dT vs ΔH (energy conservation)", len(T), dev, dev, "internal check",
                        metric="rel. dev %"))
        if save:
            plot_cp_vs_T(curve, title=spec["title"], name=f"fig_cpT_{key}", outdir=FIGS)
    if save:
        curves = {}
        for wdt in (1.0, 2.0, 4.0):
            inp = {"mass_fractions": SYSTEMS["pcm_octadecane_20wt"]["mass_fractions"], "melt_width": wdt}
            curves[f"melt width {wdt:g} K"] = engine.cp_curve(inp, np.linspace(292, 310, 721))
        plot_apparent_peak(curves, "20 wt% n-octadecane in water", name="fig_pcm_peak", outdir=FIGS)
        wp = np.linspace(0.5, 72, 120)
        xw = ((100 - wp) / w.MW) / ((100 - wp) / w.MW + wp / db.get("ethanol").MW)
        shown = {DEFAULT: r"NRTL $\gamma$, H$^E$+SLE fit (engine default)",
                 "ideal + ΔCp (option sle_activity='ideal')": r"ideal + $\Delta C_p$", "ideal": "ideal"}
        curves_m = {disp: (wp, [safe(models[lab], x) - 273.15 for x in xw]) for lab, disp in shown.items()}
        plot_data_vs_models({"Lange (1967)": {"data": (data.w_ethanol_percent, data.freezing_point_K - 273.15),
                                              "models": curves_m}},
                            "ethanol / wt%", "freezing point / °C", "Ice + aqueous ethanol: liquidus",
                            name="fig_freezing_points", outdir=FIGS)
        data.to_csv(FIGS.parent / "validation" / "results" / "freezing_point_comparison.csv", index=False)
    return rows


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    print(pd.DataFrame(run()).to_string(index=False))
