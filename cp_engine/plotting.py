"""Publication-quality figures.

Style rules used throughout (so all figures look like one set):
* one fixed categorical colour order (colour-blind-validated palette),
  never cycled; line style + marker as secondary encoding so figures also
  read in grayscale print
* thin lines, recessive grid, black/grey text, legend for ≥ 2 series
* one y-axis per panel (contributions of very different size get their own
  panel instead of a second y-scale)
* every figure saved as PNG and PDF at 300 dpi
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Sequence

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .database import Database
from .models.heterogeneous import nanofluid_cp
from .models.nrtl import NRTL

# Categorical slots, fixed order (blue, orange, aqua, yellow, magenta, green, violet, red)
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
STYLES = ["-", "--", "-.", ":"]
MARKERS = ["o", "s", "^", "D"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"


def formula_tex(formula: str) -> str:
    """'Al2O3' -> 'Al$_2$O$_3$' for axis labels and titles."""
    return re.sub(r"(\d+)", lambda m: "$_{" + m.group(1) + "}$", formula)


def apply_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 300, "font.size": 9, "font.family": "DejaVu Sans",
        "axes.labelsize": 9, "axes.titlesize": 9.5, "legend.fontsize": 8, "xtick.labelsize": 8,
        "ytick.labelsize": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
        "ytick.color": INK2, "text.color": INK, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.6, "lines.markersize": 5, "legend.frameon": False,
        "axes.prop_cycle": matplotlib.cycler(color=COLORS),
    })


def save(fig, name: str, outdir: str | Path | None = None) -> list[Path]:
    """Save `fig` as <name>.png and <name>.pdf (300 dpi)."""
    outdir = Path(outdir) if outdir else FIG_DIR
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("png", "pdf"):
        p = outdir / f"{name}.{ext}"
        fig.savefig(p, dpi=300, bbox_inches="tight")
        paths.append(p)
    return paths


def _new(nrows=1, ncols=1, w=3.4, h=2.6, **kw):
    apply_style()
    return plt.subplots(nrows, ncols, figsize=(w * ncols, h * nrows), constrained_layout=True, **kw)


# ---------------------------------------------------------------------------
def plot_cp_vs_T(curve: pd.DataFrame, title: str = "", name: str | None = None, outdir=None):
    """Cp(T) with ideal, excess and latent contributions (two panels, shared T axis).

    `curve` is the DataFrame returned by CpEngine.cp_curve().
    """
    has_small = (curve["excess"].abs().max() > 0) or (curve["latent"].abs().max() > 0)
    fig, axes = _new(2 if has_small else 1, 1, w=3.6, h=2.2, sharex=True)
    axes = np.atleast_1d(axes)
    T = curve["T_K"]
    ax = axes[0]
    ax.plot(T, curve["total"], color=COLORS[0], ls=STYLES[0], label="total (apparent)")
    ax.plot(T, curve["ideal"], color=COLORS[1], ls=STYLES[1], label="ideal Σ w$_i$c$_{p,i}$")
    ax.set_ylabel("c$_p$ / J kg$^{-1}$ K$^{-1}$")
    ax.legend(loc="best")
    if title:
        ax.set_title(title)
    if has_small:
        ax2 = axes[1]
        ax2.axhline(0, color=INK2, lw=0.6)
        ax2.plot(T, curve["excess"], color=COLORS[2], ls=STYLES[2], label="excess c$_p^E$")
        ax2.plot(T, curve["latent"], color=COLORS[3], ls=STYLES[3], label="latent / apparent")
        ax2.set_ylabel("contribution /\nJ kg$^{-1}$ K$^{-1}$")
        ax2.legend(loc="best")
    axes[-1].set_xlabel("T / K")
    if name:
        save(fig, name, outdir)
    return fig


def plot_cp_vs_composition(db: Database, comp1: str, comp2: str, T: float, engine, name=None, outdir=None):
    """Binary liquid mixture: real (NRTL) and ideal cp vs mole fraction of comp1."""
    xs = np.linspace(0.0, 1.0, 41)
    real, ideal = [], []
    for x in xs:
        inp = {"T": T, "mole_fractions": {comp1: x, comp2: 1 - x}}
        real.append(engine.solve(inp).cp_molar if 0 < x < 1 else None)
        ideal.append(engine.solve(dict(inp, solution_model="ideal")).cp_molar if 0 < x < 1 else None)
    c1, c2 = db.get(comp1), db.get(comp2)
    ends = {0.0: c2.cp("liquid", T), 1.0: c1.cp("liquid", T)}
    real = [ends.get(x, r) if r is None else r for x, r in zip(xs, real)]
    ideal = [ends.get(x, r) if r is None else r for x, r in zip(xs, ideal)]
    fig, ax = _new()
    ax.plot(xs, real, color=COLORS[0], ls=STYLES[0], label="real solution (NRTL)")
    ax.plot(xs, ideal, color=COLORS[1], ls=STYLES[1], label="ideal mixture")
    ax.set_xlabel(f"x$_{{\\mathrm{{{c1.name}}}}}$")
    ax.set_ylabel("C$_p$ / J mol$^{-1}$ K$^{-1}$")
    ax.set_title(f"{c1.name} + {c2.name}, T = {T:.2f} K")
    ax.legend()
    if name:
        save(fig, name, outdir)
    return fig


def plot_cpE_vs_x(db: Database, comp1: str, comp2: str, temperatures: Sequence[float],
                  exp: pd.DataFrame | None = None, name=None, outdir=None):
    """Excess heat capacity Cp^E(x) of a binary at up to four temperatures (+ data if given)."""
    m = NRTL([comp1, comp2], db)
    xs = np.linspace(0.0, 1.0, 101)
    fig, ax = _new()
    for i, T in enumerate(list(temperatures)[:4]):
        ax.plot(xs, [m.cpE([x, 1 - x], T) for x in xs], color=COLORS[i], ls=STYLES[i],
                label=f"NRTL, {T:.2f} K")
    if exp is not None and len(exp):
        for i, (T, g) in enumerate(exp.groupby("T_K")):
            ax.plot(g["x1"], g["CpE_J_molK"], ls="none", marker=MARKERS[i % 4], mfc="white",
                    color=COLORS[i % 8], label=f"exp., {T:.2f} K")
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xlabel(f"x$_{{\\mathrm{{{db.get(comp1).name}}}}}$")
    ax.set_ylabel("C$_p^E$ / J mol$^{-1}$ K$^{-1}$")
    ax.legend()
    if name:
        save(fig, name, outdir)
    return fig


def plot_nanofluid_models(db: Database, particle: str, fluid: str, T: float, phi_max: float = 0.1,
                          exp: pd.DataFrame | None = None, name=None, outdir=None):
    """Model I (volume-weighted) vs Model II (thermal equilibrium) vs experimental data."""
    phis = np.linspace(0.0, phi_max, 51)
    fig, ax = _new()
    ax.plot(phis * 100, [nanofluid_cp(db, particle, fluid, p, T, "I") for p in phis], color=COLORS[1],
            ls=STYLES[1], label="Model I: φc$_{p,s}$ + (1−φ)c$_{p,f}$")
    ax.plot(phis * 100, [nanofluid_cp(db, particle, fluid, p, T, "II") for p in phis], color=COLORS[0],
            ls=STYLES[0], label="Model II: thermal equilibrium")
    if exp is not None and len(exp):
        ax.plot(exp["phi"] * 100, exp["cp_J_kgK"], ls="none", marker="o", mfc="white", color=INK,
                label="experiment")
    ax.set_xlabel("particle volume fraction φ / %")
    ax.set_ylabel("c$_p$ / J kg$^{-1}$ K$^{-1}$")
    ax.set_title(f"{formula_tex(db.get(particle).formula)}–{db.get(fluid).name.lower()} nanofluid, {T:.2f} K")
    ax.legend()
    if name:
        save(fig, name, outdir)
    return fig


def plot_apparent_peak(curves: dict[str, pd.DataFrame], title: str = "", name=None, outdir=None):
    """Apparent cp around a phase change, for up to four cases (e.g. melt widths)."""
    fig, ax = _new(w=3.6, h=2.7)
    for i, (label, c) in enumerate(list(curves.items())[:4]):
        ax.plot(c["T_K"], c["total"] / 1000.0, color=COLORS[i], ls=STYLES[i], label=label)
    ax.set_xlabel("T / K")
    ax.set_ylabel("apparent c$_p$ / kJ kg$^{-1}$ K$^{-1}$")
    if title:
        ax.set_title(title)
    if len(curves) > 1:
        ax.legend()
    if name:
        save(fig, name, outdir)
    return fig


def plot_validation(T, model, exp_T, exp_val, ylabel: str, title: str, label_exp="reference data",
                    name=None, outdir=None):
    """Model curve vs reference points, with a deviation panel underneath."""
    fig, (ax, ax2) = _new(2, 1, w=3.6, h=2.0, sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    ax.plot(T, model, color=COLORS[0], ls=STYLES[0], label="this engine")
    ax.plot(exp_T, exp_val, ls="none", marker="o", mfc="white", color=INK, label=label_exp)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    m_at = np.interp(exp_T, T, model)
    ax2.axhline(0, color=INK2, lw=0.6)
    ax2.plot(exp_T, (m_at - np.asarray(exp_val)) / np.asarray(exp_val) * 100, ls="none",
             marker="s", color=COLORS[1], ms=4)
    ax2.set_ylabel("deviation / %")
    ax2.set_xlabel("T / K")
    if name:
        save(fig, name, outdir)
    return fig


def plot_parity(series: dict[str, tuple], quantity: str, unit: str, title: str = "", name=None, outdir=None):
    """Predicted vs measured for up to four models; series = {label: (measured, predicted)}."""
    fig, ax = _new(w=3.3, h=3.1)
    lo = min(min(np.min(m), np.min(p)) for m, p in series.values())
    hi = max(max(np.max(m), np.max(p)) for m, p in series.values())
    pad = 0.03 * (hi - lo)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=INK2, lw=0.8, label="y = x")
    for i, (label, (meas, pred)) in enumerate(list(series.items())[:4]):
        ax.plot(meas, pred, ls="none", marker=MARKERS[i], mfc="white" if i else COLORS[i], color=COLORS[i],
                label=label)
    ax.set_xlabel(f"measured {quantity} / {unit}")
    ax.set_ylabel(f"predicted {quantity} / {unit}")
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_aspect("equal")
    if title:
        ax.set_title(title)
    ax.legend(loc="upper left")
    if name:
        save(fig, name, outdir)
    return fig


def plot_data_vs_models(groups: dict[str, dict], xlabel: str, ylabel: str, title: str = "",
                        style_legend: dict[str, str] | None = None, name=None, outdir=None, w=3.6, h=2.8,
                        legend_outside: bool = False):
    """One colour per group (≤ 4, e.g. temperatures); within a group: markers = data,
    line styles = models. groups = {label: {"data": (x, y), "models": {model: (x, y)}}}."""
    fig, ax = _new(w=w, h=h)
    model_names = []
    for i, (glabel, g) in enumerate(list(groups.items())[:4]):
        if "data" in g and g["data"] is not None:
            ax.plot(*g["data"], ls="none", marker=MARKERS[i], mfc="white", color=COLORS[i], ms=4.5)
        for j, (mname, (x, y)) in enumerate(g.get("models", {}).items()):
            ax.plot(x, y, color=COLORS[i], ls=STYLES[j % 4], lw=1.4)
            if mname not in model_names:
                model_names.append(mname)
    # legend: colours = groups (with data marker), line styles = models
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=COLORS[i], marker=MARKERS[i], mfc="white", ls="none", label=gl)
               for i, gl in enumerate(list(groups)[:4])]
    handles += [Line2D([], [], color=INK2, ls=STYLES[j % 4], label=m) for j, m in enumerate(model_names)]
    if legend_outside:
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7)
    else:
        ax.legend(handles=handles, loc="best", fontsize=7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if name:
        save(fig, name, outdir)
    return fig
