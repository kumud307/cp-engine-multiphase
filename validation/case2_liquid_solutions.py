"""Validation case 2 – real liquid solutions: excess enthalpy and mixture heat capacity.

All experimental data come from the NIST TRC ThermoML archive (validation/data/thermoml_binary_data.csv,
extracted by tools/extract_thermoml.py); each row carries the DOI of the original paper.

2a  water + ethanol   H^E, 298-383 K  (Nagamachi & Francesconi 2006; Brandt et al. 2014; Modaresifar & Kowalski 2017)
2b  water + methanol  H^E, 298-378 K  (Behbehani & Ghamamy 2006; Dettmann et al. 2006)
    -> NRTL with VLE parameters vs NRTL refitted to H^E (tools/fit_nrtl_he.py)
2c  water + methanol  Cp, 280 and 320 K, 0.1 MPa (Tanaka, Fujita & Uematsu 2007) – INDEPENDENT test:
    the H^E-fitted model was never fitted to Cp, so this checks the Gibbs–Helmholtz route to Cp^E.
2d  ethanol + toluene Cp, 313.15 K, 0.1 MPa (Vega-Maza et al. 2011) and isobaric VLE (Chen et al. 2012):
    NRTL fitted jointly to both (held-out tests in case 2b); the ideal mixture is shown for comparison.
2j  activity coefficients at infinite dilution in water (Straver & de Loos 2005; Atik et al. 2004) –
    independent of every fit.
2e  water + glycerol  H^E, 298.15 K (Zehioua et al. 2009) – shows the size of H^E that the ideal
    fallback ignores for this pair.
2f  water + ethanol   Cp – still PENDING: no open dataset was found (validation/data/water_ethanol_cp_exp.csv).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import FIGS, aard, load, pending, row
from cp_engine import CpEngine
from cp_engine.models.nrtl import NRTL
from cp_engine.plotting import plot_data_vs_models, plot_parity

DATA = "thermoml_binary_data.csv"


def binary(d: pd.DataFrame, a: str, b: str, prop: str) -> pd.DataFrame:
    """Rows for pair (a, b) and property, with x = mole fraction of a."""
    g = d[(((d.comp1 == a) & (d.comp2 == b)) | ((d.comp1 == b) & (d.comp2 == a))) & (d.property == prop)].copy()
    g["x"] = np.where(g.comp1 == a, g.x1, 1 - g.x1)
    return g


def refs(g: pd.DataFrame) -> str:
    return "; ".join(f"doi:{x}" for x in sorted(g.doi.unique()))


def run(save: bool = True) -> list[dict]:
    engine = CpEngine()
    db = engine.db
    d = load(DATA)
    rows = []

    # ---- 2a/2b: excess enthalpies, VLE vs H^E-fitted NRTL ------------------------------------
    for alc in ("ethanol", "methanol"):
        he = binary(d, alc, "water", "HE")
        he = he[(he.x > 0.01) & (he.x < 0.99) & (he.T_K <= 383.2)]
        groups = {}
        for label, prefer in (("NRTL (VLE parameters)", "VLE"), ("NRTL refitted (fit quality: these data were used in the fit; held-out tests in case 2b)", "HE")):
            m = NRTL([alc, "water"], db, prefer=prefer)
            pred = np.array([m.hE([x, 1 - x], T) for x, T in zip(he.x, he.T_K)])
            rmsd = float(np.sqrt(np.mean((pred - he.value) ** 2)))
            rows.append(row("2", f"{alc} + water", f"H^E, {label.replace('$^E$', '^E')}", len(he),
                            rmsd, float(np.max(np.abs(pred - he.value))), refs(he), metric="RMSD J/mol"))
        if save:
            Ts = sorted(he.T_K.unique())
            show = [T for T in Ts if any(abs(T - t) < 0.1 for t in (298.15, 323.15, 363.15))] or Ts[:3]
            xs = np.linspace(0.001, 0.999, 200)
            for T in show:
                g = he[(he.T_K - T).abs() < 0.1]
                models = {}
                for label, prefer in (("NRTL refitted to H$^E$", "HE"), ("NRTL, VLE parameters", "VLE")):
                    m = NRTL([alc, "water"], db, prefer=prefer)
                    models[label] = (xs, [m.hE([x, 1 - x], T) for x in xs])
                groups[f"{T:.2f} K"] = {"data": (g.x, g.value), "models": models}
            plot_data_vs_models(groups, f"x$_{{\\mathrm{{{alc}}}}}$", "H$^E$ / J mol$^{-1}$",
                                f"Water + {alc}: excess enthalpy", name=f"fig_HE_water_{alc}", outdir=FIGS,
                                w=4.8, legend_outside=True)

    # ---- 2c: water + methanol Cp (independent) ---------------------------------------------------
    cp = binary(d, "methanol", "water", "Cp_molar")
    cp = cp[cp.P_kPa <= 200]
    preds = {}
    for label, opt in (("ideal mixture", {"solution_model": "ideal"}),
                       ("NRTL, VLE parameters", {"nrtl_parameters": "VLE"}),
                       ("NRTL refitted to H^E", {"nrtl_parameters": "HE"})):
        p = np.array([engine.solve(dict(T=t, mole_fractions={"methanol": x, "water": 1 - x}, **opt)).cp_molar
                      for x, t in zip(cp.x, cp.T_K)])
        preds[label] = p
        a, mx, n = aard(p, cp.value)
        rows.append(row("2", "methanol + water", f"Cp, {label}", n, a, mx, refs(cp)))
    if save:
        plot_parity({k: (cp.value.values, v) for k, v in preds.items()}, "C$_p$", "J mol$^{-1}$ K$^{-1}$",
                    "Water + methanol, 0.1 MPa", name="fig_cp_parity_methanol_water", outdir=FIGS)
        xs = np.linspace(0.0, 1.0, 41)
        groups = {}
        for T in sorted(cp.T_K.unique()):
            g = cp[cp.T_K == T]
            mod = {}
            for label, opt in (("NRTL refitted to H$^E$", {"nrtl_parameters": "HE"}),
                               ("ideal mixture", {"solution_model": "ideal"})):
                mod[label] = (xs, [engine.solve(dict(T=T, mole_fractions={"methanol": x, "water": 1 - x}, **opt)).cp_molar
                                   for x in xs])
            groups[f"{T:.0f} K"] = {"data": (g.x, g.value), "models": mod}
        plot_data_vs_models(groups, "x$_{\\mathrm{methanol}}$", "C$_p$ / J mol$^{-1}$ K$^{-1}$",
                            "Water + methanol heat capacity", name="fig_cp_methanol_water", outdir=FIGS,
                            w=4.8, legend_outside=True)

    # ---- 2d: ethanol + toluene Cp (NRTL fitted to VLE + Cp^E) and VLE ---------------------------
    et = binary(d, "ethanol", "toluene", "Cp_molar")
    et = et[et.P_kPa <= 200]
    p = np.array([engine.solve({"T": t, "mole_fractions": {"ethanol": x, "toluene": 1 - x}}).cp_molar
                  for x, t in zip(et.x, et.T_K)])
    a, mx, n = aard(p, et.value)
    rows.append(row("2", "ethanol + toluene", "Cp, NRTL fitted to VLE + Cp^E (fit quality; held-out in case 2b)",
                    n, a, mx, refs(et)))
    p0 = np.array([engine.solve({"T": t, "mole_fractions": {"ethanol": x, "toluene": 1 - x}, "solution_model": "ideal"}).cp_molar
                   for x, t in zip(et.x, et.T_K)])
    a, mx, n = aard(p0, et.value)
    rows.append(row("2", "ethanol + toluene", "Cp, ideal mixture (for comparison)", n, a, mx, refs(et)))
    # VLE of ethanol + toluene: bubble temperatures predicted with the stored NRTL + NIST Antoine Psat
    from scipy.optimize import brentq
    vle = load("thermoml_vle_ethanol_toluene.csv")
    vle = vle[(vle.x_ethanol > 0) & (vle.x_ethanol < 1)]
    m_et = NRTL(["ethanol", "toluene"], engine.db)
    e_, t_ = engine.db.get("ethanol"), engine.db.get("toluene")

    def bubble_T(x, P):
        f = lambda T: (x * np.exp(m_et.ln_gamma([x, 1 - x], T)[0]) * e_.psat(T)
                       + (1 - x) * np.exp(m_et.ln_gamma([x, 1 - x], T)[1]) * t_.psat(T)) - P
        return brentq(f, 300.0, 450.0)
    dT = [abs(bubble_T(x, P * 1e3) - T) for x, P, T in zip(vle.x_ethanol, vle.P_kPa, vle.T_K)]
    rows.append(row("2", "ethanol + toluene", "VLE bubble temperature (fit quality)", len(dT), float(np.mean(dT)),
                    float(np.max(dT)), "doi:10.1021/je200921u", metric="|ΔT| K"))

    # activity coefficients at infinite dilution in water (independent: never used in a fit)
    gi = load("thermoml_gamma_inf.csv")
    for solute, prefer, label in (("ethanol", "HE", "Cp-optimised set (H^E+SLE)"), ("ethanol", "VLE", "VLE set"),
                                  ("acetone", "VLE", "VLE set (literature)")):
        g = gi[gi.solute == solute]
        m = NRTL([solute, "water"], engine.db, prefer=prefer)
        pred = [np.exp(m.ln_gamma([1e-9, 1 - 1e-9], t)[0]) for t in g.T_K]
        a, mx, n = aard(pred, g.gamma_inf)
        rows.append(row("2", f"{solute} in water", f"activity coefficient at infinite dilution, {label}", n, a, mx,
                        "; ".join(f"doi:{d}" for d in sorted(g.doi.unique()))))

    # ---- 2e: glycerol + water H^E vs the ideal fallback (H^E = 0) -------------------------------
    gw = binary(d, "glycerol", "water", "HE")
    rows.append(row("2", "glycerol + water", "H^E, ideal fallback (H^E = 0)", len(gw),
                    float(np.sqrt(np.mean(gw.value ** 2))), float(gw.value.abs().max()), refs(gw), metric="RMSD J/mol"))

    # ---- 2g: water + ethanol Cp^E derived directly from measured H^E(T) ----------------------------
    # Cp^E = (∂H^E/∂T)_P,x. Nagamachi & Francesconi (2006) measured H^E at 298.15 and 323.15 K with
    # the same apparatus; a Redlich–Kister fit of each isotherm gives H^E(x) and the finite difference
    # [H^E(323.15) − H^E(298.15)]/25 K is a model-free estimate of the mean Cp^E over that interval.
    he = binary(d, "ethanol", "water", "HE")
    he = he[he.doi == "10.1016/j.jct.2005.06.018"]
    iso = {T: he[(he.T_K - T).abs() < 0.1] for T in (298.15, 323.15)}
    if all(len(g) >= 6 for g in iso.values()):
        def rk(g):            # H^E = x(1−x) Σ A_k (2x−1)^k, 4 terms
            X = np.column_stack([g.x * (1 - g.x) * (2 * g.x - 1) ** k for k in range(4)])
            return np.linalg.lstsq(X, g.value, rcond=None)[0]
        A = {T: rk(g) for T, g in iso.items()}
        xs = np.linspace(0.05, 0.95, 19)
        h = lambda T, x: x * (1 - x) * sum(a * (2 * x - 1) ** k for k, a in enumerate(A[T]))
        cpe_data = (h(323.15, xs) - h(298.15, xs)) / 25.0
        model_cpe = {}
        for label, prefer in (("NRTL, VLE parameters", "VLE"), ("NRTL refitted to H^E", "HE")):
            m = NRTL(["ethanol", "water"], db, prefer=prefer)
            mc = np.array([(m.hE([x, 1 - x], 323.15) - m.hE([x, 1 - x], 298.15)) / 25.0 for x in xs])
            model_cpe[label] = mc
            rows.append(row("2", "ethanol + water", f"mean Cp^E 298-323 K from measured H^E(T), {label}"
                            + (" (consistency check: data used in the fit; held-out version in case 2b)" if "refitted" in label else ""), len(xs),
                            float(np.sqrt(np.mean((mc - cpe_data) ** 2))), float(np.max(np.abs(mc - cpe_data))),
                            "doi:10.1016/j.jct.2005.06.018", metric="RMSD J/(mol K)"))
        pd.DataFrame({"x_ethanol": xs, "CpE_from_HE_data": cpe_data, **model_cpe}).to_csv(
            FIGS.parent / "validation" / "results" / "ethanol_water_CpE_from_HE.csv", index=False)
        if save:
            plot_data_vs_models({"310.65 K (mean of 298–323 K)": {
                "data": (xs, cpe_data),
                "models": {"NRTL refitted to H$^E$": (xs, model_cpe["NRTL refitted to H^E"]),
                           "NRTL, VLE parameters": (xs, model_cpe["NRTL, VLE parameters"])}}},
                "x$_{\\mathrm{ethanol}}$", "C$_p^E$ / J mol$^{-1}$ K$^{-1}$",
                "Water + ethanol: C$_p^E$ from measured H$^E$(T)", name="fig_cpE_from_HE_water_ethanol",
                outdir=FIGS, w=4.8, legend_outside=True)

    # ---- 2h: NaCl(aq) – apparent molar Cp model (tools/fit_apparent_cp.py) ---------------------------
    na = binary(d, "NaCl", "water", "Cp_molar")
    na = na[na.doi != "10.1016/j.jct.2018.10.027"]           # excluded as inconsistent (see the fitting tool)
    for label, sel in (("fit points (doi:10.1021/je0502400)", na.doi == "10.1021/je0502400"),
                       ("INDEPENDENT point (doi:10.1016/j.jct.2016.09.031)", na.doi == "10.1016/j.jct.2016.09.031")):
        g = na[sel]
        p = [engine.solve({"T": t, "mole_fractions": {"NaCl": x, "water": 1 - x}}).cp_molar for x, t in zip(g.x, g.T_K)]
        a, mx, n = aard(p, g.value)
        rows.append(row("2", "NaCl + water (brine)", f"Cp, apparent-molar-Cp model, {label}", n, a, mx, refs(g)))
    nk = [engine.db.get("water").cp("liquid", t) * (1 - x) + engine.db.get("NaCl").cp("solid", t) * x
          for x, t in zip(na.x, na.T_K)]
    a, mx, n = aard(nk, na.value)
    rows.append(row("2", "NaCl + water (brine)", "Cp, previous Neumann–Kopp model (for comparison)", n, a, mx, refs(na)))

    # ---- 2i: glycerol + water – apparent molar Cp data vs the ideal fallback ------------------------
    ph = binary(d, "glycerol", "water", "phiCp")
    if len(ph):
        cpw = np.array([engine.db.get("water").cp("liquid", t) for t in ph.T_K])
        cp_data = ph.x * ph.value + (1 - ph.x) * cpw              # mixture Cp from φCp of glycerol
        p = [engine.solve({"T": t, "mole_fractions": {"glycerol": x, "water": 1 - x}}).cp_molar
             for x, t in zip(ph.x, ph.T_K)]
        a, mx, n = aard(p, cp_data)
        rows.append(row("2", "glycerol + water", "Cp (x ≤ 0.08, 278-363 K), ideal fallback", n, a, mx, refs(ph)))

    # ---- 2f: water + ethanol Cp – pending ----------------------------------------------------------
    exp = load("water_ethanol_cp_exp.csv")
    if len(exp) and exp.Cp_J_molK.notna().any():
        e = exp.dropna(subset=["Cp_J_molK"])
        p = [engine.solve({"T": t, "mole_fractions": {"ethanol": x, "water": 1 - x}}).cp_molar
             for x, t in zip(e.x_ethanol, e.T_K)]
        a, mx, n = aard(p, e.Cp_J_molK)
        rows.append(row("2", "ethanol + water", "Cp, NRTL refitted to H^E", n, a, mx, "; ".join(e.reference.unique())))
    else:
        rows.append(pending("2", "ethanol + water", "Cp (mixture)", "water_ethanol_cp_exp.csv"))
    return rows


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    print(pd.DataFrame(run()).to_string(index=False))
