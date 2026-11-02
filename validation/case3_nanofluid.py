"""Validation case 3 – nanofluids: Model I vs Model II vs DSC data.

Data: O'Hanley, Buongiorno, McKrell & Hu, Adv. Mech. Eng. 2012, 181079 (open access): alumina, copper
oxide and silica in water, 35/45/55 °C, 5-50 wt%. Values digitised from the vector figures
(tools/digitize_ohanley2012.py → validation/data/nanofluid_ohanley2012.csv).

Model I  (volume-weighted)   uses the paper's nominal volume fraction (Table 2)
Model II (thermal equilibrium) = Route A of the engine, which only needs the MASS fraction (Table 2):
         cp = w cp_s + (1 − w) cp_f  – no particle density is required.
Note: the silica in the paper is amorphous colloidal silica (Ludox); the database holds α-quartz Cp,
which differs from amorphous silica by < 1 % near room temperature.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import FIGS, aard, load, row
from cp_engine import CpEngine
from cp_engine.models.heterogeneous import cp_model_I
from cp_engine.plotting import plot_nanofluid_models, plot_parity

CSV = "nanofluid_ohanley2012.csv"
REF = "O'Hanley et al. 2012, doi:10.1155/2012/181079 (digitised)"


def run(save: bool = True) -> list[dict]:
    engine = CpEngine()
    db = engine.db
    exp = load(CSV)
    m1, m2 = [], []
    for r in exp.itertuples():
        s, f = db.get(r.particle), db.get("water")
        m1.append(cp_model_I(r.phi_table2, s.cp_mass("solid", r.T_K), f.cp_mass("liquid", r.T_K)))
        m2.append(engine.solve({"T": r.T_K, "mass_fractions": {r.particle: r.w_table2, "water": None}}).cp_mass)
    exp = exp.assign(model_I=m1, model_II=m2)
    rows = []
    for particle, g in exp.groupby("particle"):
        for col, label in (("model_I", "Model I (volume-weighted)"), ("model_II", "Model II = Route A (mass-weighted)")):
            a, mx, n = aard(g[col], g.cp_J_kgK)
            rows.append(row("3", f"{particle} + water", f"cp, {label}", n, a, mx, REF))
    for col, label in (("model_I", "Model I"), ("model_II", "Model II = Route A")):
        a, mx, n = aard(exp[col], exp.cp_J_kgK)
        rows.append(row("3", "all nanofluids", f"cp, {label}", n, a, mx, REF))
    if save:
        al = exp[(exp.particle == "Al2O3") & (exp.T_K.round(2) == 308.15)]
        plot_nanofluid_models(db, "Al2O3", "water", 308.15, 0.10,
                              exp=pd.DataFrame({"phi": al.phi_table2, "cp_J_kgK": al.cp_J_kgK}),
                              name="fig_nanofluid_models", outdir=FIGS)
        cu = exp[(exp.particle == "CuO") & (exp.T_K.round(2) == 308.15)]
        plot_nanofluid_models(db, "CuO", "water", 308.15, 0.15,
                              exp=pd.DataFrame({"phi": cu.phi_table2, "cp_J_kgK": cu.cp_J_kgK}),
                              name="fig_nanofluid_models_cuo", outdir=FIGS)
        plot_parity({"Model I (volume-weighted)": (exp.cp_J_kgK / 1000, np.array(m1) / 1000),
                     "Model II = Route A": (exp.cp_J_kgK / 1000, np.array(m2) / 1000)},
                    "c$_p$", "kJ kg$^{-1}$ K$^{-1}$", "Nanofluids (Al$_2$O$_3$, CuO, SiO$_2$ in water)",
                    name="fig_nanofluid_parity", outdir=FIGS)
    return rows


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    print(pd.DataFrame(run()).to_string(index=False))
