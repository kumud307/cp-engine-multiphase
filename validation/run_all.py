"""Run every validation case and write the summary.

Outputs (validation/results/):
  summary.csv            every comparison (machine-readable)
  summary.md             same, as a Markdown table
  summary_table.tex      condensed LaTeX table for the manuscript (key results only)
  summary_table_full.tex full LaTeX table (supporting information)
"""
from __future__ import annotations

import pandas as pd

import case1_pure_components
import case2_liquid_solutions
import case2b_crossvalidation
import case3_nanofluid
import case4_phase_change
import case5_route_consistency
from common import RESULTS

TEX = {"%": r"\%", "_": " ", "∫": r"$\int$", "ΔH": r"$\Delta H$", "ΔCp": r"$\Delta C_p$", "|ΔT|": r"$|\Delta T|$",
       "α²": r"$\alpha^2$", "κ T": r"$\kappa_T$", "φ": r"$\phi$", "^E": r"$^E$", "γ": r"$\gamma$", "≤": r"$\le$",
       "–": "--", "≥": r"$\ge$", "Σ": r"$\Sigma$"}


def tex(s: str) -> str:
    s = str(s)
    for k in ("%", "_"):
        s = s.replace(k, TEX[k])
    for k, v in TEX.items():
        if k not in ("%", "_"):
            s = s.replace(k, v)
    return s


def fmt(v) -> str:
    if pd.isna(v):
        return "–"
    if v == 0:
        return "0"
    return f"{v:.1e}" if abs(v) < 1e-3 else (f"{v:.3g}" if abs(v) < 100 else f"{v:.0f}")


def condensed(df: pd.DataFrame) -> pd.DataFrame:
    """Key results for the manuscript (aggregates case 1 and case 5)."""
    ok = df[df.status == "ok"]
    rows = []
    c1 = ok[(ok.case == "1") & ok.quantity.str.contains("1 atm")]
    rows.append({"case": "1", "what": "5 pure liquids vs NIST SRD 69 (0.1 MPa)", "metric": "AARD %",
                 "N": c1.N.sum(), "mean": (c1["mean"] * c1.N).sum() / c1.N.sum(), "max": c1["max"].max()})
    c1p = ok[(ok.case == "1") & ok.quantity.str.contains("point")]
    rows.append({"case": "1", "what": "4 compounds, literature point values", "metric": "AARD %",
                 "N": c1p.N.sum(), "mean": (c1p["mean"] * c1p.N).sum() / c1p.N.sum(), "max": c1p["max"].max()})
    for tag, label in (("Tait pressure correction", "5 liquids, cp at 0.1-100 MPa (Tait pressure correction)"),
                       ("without pressure correction", "5 liquids, cp at 0.1-100 MPa, no pressure correction")):
        c = ok[(ok.case == "1") & ok.quantity.str.contains(tag, regex=False)]
        rows.append({"case": "1", "what": label, "metric": "AARD %", "N": c.N.sum(),
                     "mean": (c["mean"] * c.N).sum() / c.N.sum(), "max": c["max"].max()})
    for _, r in ok[(ok.case == "1") & ok.quantity.str.contains("cv")].iterrows():
        rows.append({"case": "1", "what": "water: cp from NIST cv, cp = cv + Tvα²/κ T", "metric": r.metric,
                     "N": r.N, "mean": r["mean"], "max": r["max"]})
    for _, r in ok[ok.case.isin(["2", "2b", "3"])].iterrows():
        if r.case == "3" and r.system != "all nanofluids":
            continue
        rows.append({"case": r.case, "what": f"{r.system}: {r.quantity}", "metric": r.metric, "N": r.N,
                     "mean": r["mean"], "max": r["max"]})
    for _, r in ok[(ok.case == "4")].iterrows():
        if "all, 2.5" in r.quantity:
            continue
        rows.append({"case": r.case, "what": f"{r.system}: {r.quantity}", "metric": r.metric, "N": r.N,
                     "mean": r["mean"], "max": r["max"]})
    c5 = ok[ok.case == "5"]
    rows.append({"case": "5", "what": f"Route B vs Route A, {len(c5)} systems", "metric": "rel. dev %",
                 "N": len(c5), "mean": c5["mean"].mean(), "max": c5["max"].max()})
    return pd.DataFrame(rows)


def write_tex(df: pd.DataFrame, path, what_col: str) -> None:
    lines = [r"\begin{tabular}{c p{8.6cm} c r r r}", r"\toprule",
             r"Case & Comparison & Metric & $N$ & mean & max \\", r"\midrule"]
    for _, r in df.iterrows():
        lines.append(f"{r['case']} & {tex(r[what_col])} & {tex(r['metric'])} & {r['N']} & {fmt(r['mean'])} & "
                     f"{fmt(r['max'])} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(save: bool = True) -> pd.DataFrame:
    rows = []
    for mod in (case1_pure_components, case2_liquid_solutions, case2b_crossvalidation, case3_nanofluid,
                case4_phase_change, case5_route_consistency):
        print(f"Running {mod.__name__} ...")
        rows += mod.run(save=save)
    df = pd.DataFrame(rows)
    RESULTS.mkdir(exist_ok=True)
    df.to_csv(RESULTS / "summary.csv", index=False)
    show = df.copy()
    for c in ("mean", "max"):
        show[c] = show[c].map(fmt)
    md = ["| " + " | ".join(show.columns) + " |", "|" + "---|" * len(show.columns)]
    md += ["| " + " | ".join(str(v) for v in r) + " |" for r in show.itertuples(index=False)]
    (RESULTS / "summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    full = df[df.status == "ok"].assign(what=lambda d: d.system + ": " + d.quantity)
    write_tex(full, RESULTS / "summary_table_full.tex", "what")
    write_tex(condensed(df), RESULTS / "summary_table.tex", "what")
    return df


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 70)
    print(main()[["case", "system", "quantity", "metric", "N", "mean", "max", "status"]].to_string(index=False))
