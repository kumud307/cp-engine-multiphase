"""Validation case 6 – vapour–liquid equilibrium of mixtures, real-gas Cp and the cubic equations of state.

6a  Bubble temperatures of water + ethanol and water + methanol: 26 papers in the NIST ThermoML archive
    (validation/data/thermoml_vle_water_alcohols.csv, 13–101 kPa). Two papers whose data disagree with all
    the others by more than 2 K on average are excluded (listed in EXCLUDED_PAPERS).  Compared: the NRTL
    set used for phase equilibrium (VLE-fitted, engine default) and the caloric (H^E-fitted) set.
6b  Cp of superheated vapours at 0.5–10 MPa (NIST SRD 69 isobars): ideal gas vs ideal-gas Cp + residual
    Cp^R from the Peng–Robinson, Soave–Redlich–Kwong and virial (Pitzer/Abbott B) equations.
6c  Vapour pressure from each cubic equation (equal-fugacity condition) vs NIST, 0.55 ≤ T_r ≤ 0.95:
    the effect of the α(T) function (vdW, RK: none / T_r^-½; SRK, PR: fitted to vapour pressures).
6d  Energy conservation through the boiling range of a mixture: ∫ cp_app dT = ΔH.
"""
from __future__ import annotations

import numpy as np

from common import FIGS, aard, load, row
from cp_engine import CpEngine
from cp_engine.models.nrtl import NRTL
from cp_engine.models.vle import VLE

EXCLUDED_PAPERS = {"10.1016/j.fluid.2009.11.014": "water + ethanol: mean deviation −2.3 K, up to 12 K, from the other 17 papers",
                   "10.1016/j.fluid.2009.10.002": "water + methanol: mean deviation −5.6 K, up to 14 K, from the other 8 papers"}
GAS_FLUIDS = ("water", "methanol", "n-hexane", "toluene")


def bubble_rows(engine, save):
    rows = []
    db = engine.db
    d = load("thermoml_vle_water_alcohols.csv")
    d = d[(d.x_alcohol > 0) & (d.x_alcohol < 1) & ~d.doi.isin(EXCLUDED_PAPERS)]
    for alc in ("ethanol", "methanol"):
        g = d[d.alcohol == alc]
        for prefer, label in (("VLE", "NRTL VLE set (engine default for phase equilibrium)"),
                              ("HE", "NRTL H^E-fitted set (caloric set, for comparison)")):
            v = VLE([alc, "water"], db, NRTL([alc, "water"], db, prefer=prefer))
            dT = np.array([v.bubble_T(P * 1e3, [x, 1 - x]) - T for x, P, T in zip(g.x_alcohol, g.P_kPa, g.T_K)])
            rows.append(row("6", f"{alc} + water", f"bubble temperature 13-101 kPa, {label}", len(dT),
                            float(np.mean(np.abs(dT))), float(np.max(np.abs(dT))),
                            f"{g.doi.nunique()} papers, NIST ThermoML", metric="|ΔT| K"))
    if save:
        from cp_engine.plotting import plot_data_vs_models
        g = d[(d.alcohol == "ethanol") & ((d.P_kPa - 101.3).abs() < 0.5)]
        v = VLE(["ethanol", "water"], db, NRTL(["ethanol", "water"], db, prefer="VLE"))
        vh = VLE(["ethanol", "water"], db, NRTL(["ethanol", "water"], db, prefer="HE"))
        xs = np.linspace(0.002, 0.998, 120)
        groups = {"101.3 kPa": {"data": (g.x_alcohol, g.T_K - 273.15), "models": {
            "bubble point, VLE set": (xs, [v.bubble_T(101325.0, [x, 1 - x]) - 273.15 for x in xs]),
            "dew point, VLE set": (xs, [v.dew_T(101325.0, [x, 1 - x]) - 273.15 for x in xs]),
            "bubble point, H$^E$ set": (xs, [vh.bubble_T(101325.0, [x, 1 - x]) - 273.15 for x in xs])}}}
        plot_data_vs_models(groups, "x, y (ethanol)", "T / °C", "Water + ethanol, 101.3 kPa (data: bubble points)",
                            name="fig_vle_water_ethanol", outdir=FIGS, w=4.8, legend_outside=True)
    return rows


def real_gas_rows(engine, save):
    rows = []
    for name in GAS_FLUIDS:
        from pathlib import Path
        files = sorted(Path(__file__).resolve().parent.joinpath("data").glob(f"nist_{name}_vapour_*MPa.csv"))
        for f in files:
            d = load(f.name)
            P = float(d.P_MPa.iloc[0])
            for eos in ("ideal", "PR", "SRK", "virial"):
                pred = [engine.solve({"T": t, "P": (P, "MPa"), "components": {name: "gas"}, "gas_eos": eos}).cp_mass
                        for t in d.T_K]
                a, mx, n = aard(pred, d.cp_J_kgK)
                label = "ideal gas (no residual)" if eos == "ideal" else f"ideal gas + C_p^R from {eos}"
                rows.append(row("6", f"{name} vapour, {P:g} MPa", f"cp {d.T_K.min():.0f}-{d.T_K.max():.0f} K, {label}",
                                n, a, mx, "NIST SRD 69"))
    if save:
        from cp_engine.plotting import plot_data_vs_models
        groups = {}
        for P in (1.0, 5.0, 10.0):
            d = load(f"nist_water_vapour_{P:g}MPa.csv")
            T = np.linspace(d.T_K.min(), d.T_K.max(), 60)
            groups[f"{P:g} MPa"] = {"data": (d.T_K, d.cp_J_kgK / 1000), "models": {
                eos_label: (T, [engine.solve({"T": t, "P": (P, "MPa"), "components": {"water": "gas"},
                                              "gas_eos": eos}).cp_mass / 1000 for t in T])
                for eos, eos_label in (("ideal", "ideal gas"), ("PR", "ideal gas + Peng–Robinson C$_p^R$"))}}
        plot_data_vs_models(groups, "T / K", "c$_p$ / kJ kg$^{-1}$ K$^{-1}$", "Superheated steam (NIST SRD 69)",
                            name="fig_real_gas_steam", outdir=FIGS, w=4.8, legend_outside=True)
    return rows


def psat_rows(engine):
    rows = []
    for name in GAS_FLUIDS + ("n-dodecane",):
        d = load(f"nist_sat_{name}.csv").dropna(subset=["P_MPa"])
        c = engine.db.get(name).critical
        d = d[(d.T_K >= 0.55 * c.Tc) & (d.T_K <= 0.95 * c.Tc)].iloc[::20]
        for kind in ("vdW", "RK", "SRK", "PR"):
            eos = engine.db.eos([name], kind)
            pred = [eos.psat(t) / 1e6 for t in d.T_K]
            a, mx, n = aard(pred, d.P_MPa)
            rows.append(row("6", name, f"vapour pressure 0.55-0.95 T_c from the {kind} equation (φ_L = φ_V)", n, a, mx,
                            "NIST SRD 69"))
    return rows


def energy_rows(engine):
    rows = []
    for x in (0.1, 0.5):
        inp = {"mole_fractions": {"ethanol": x, "water": 1 - x}}
        s = engine.system(inp)
        T1, T2 = 340.0, 385.0
        s.configure(T1, T2)
        from scipy import integrate
        info = next(iter(s.vle_info.values()))
        integral, _ = integrate.quad(lambda t: s.cp_point(float(t)).total, T1, T2, limit=500, epsrel=1e-10,
                                     points=[info["T_bubble"], info["T_dew"]])
        dH = s.enthalpy(T2) - s.enthalpy(T1)
        rows.append(row("6", f"ethanol + water x = {x}", "∫cp_app dT vs ΔH through the boiling range 340-385 K "
                        "(adaptive quadrature, break points at the bubble and dew points)",
                        1, abs(integral / dH - 1) * 100, abs(integral / dH - 1) * 100, "internal consistency",
                        metric="rel. dev %"))
    return rows


def run(save: bool = True) -> list[dict]:
    engine = CpEngine()
    return bubble_rows(engine, save) + real_gas_rows(engine, save) + psat_rows(engine) + energy_rows(engine)


if __name__ == "__main__":
    import pandas as pd
    pd.set_option("display.width", 200)
    print(pd.DataFrame(run()).to_string(index=False))
