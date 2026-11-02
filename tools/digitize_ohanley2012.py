"""Extract the measured nanofluid specific heats from the vector figures of O'Hanley et al. (2012).

Reference: H. O'Hanley, J. Buongiorno, T. McKrell, L.-W. Hu, "Measurement and Model Validation of
Nanofluid Specific Heat Capacity with Differential Scanning Calorimetry", Adv. Mech. Eng. 2012,
181079, doi:10.1155/2012/181079 (open access, CC BY 2.0; PDF from https://hdl.handle.net/1721.1/72068).

The paper reports its data only in Figures 1-9 (alumina, copper oxide and silica in water at 35, 45
and 55 °C). The figures are vector graphics, so every data marker's exact position can be read from
the PDF drawing commands; the axes are calibrated from the printed tick labels. This is digitisation,
not a re-measurement: the precision is limited only by the plotting (better than ±0.002 J/(g K)).

Usage:  python tools/digitize_ohanley2012.py <path to ohanley2012.pdf>
Output: validation/data/nanofluid_ohanley2012.csv
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pymupdf

ROOT = Path(__file__).resolve().parent.parent

# Table 2 of the paper: nominal (wt%, vol%) of the four concentrations of each nanofluid
TABLE2 = {
    "Al2O3": [(20.0, 6.4), (15.0, 4.6), (10.0, 2.9), (5.0, 1.4)],
    "CuO": [(50.0, 13.7), (37.5, 8.7), (25.0, 5.0), (12.5, 2.2)],
    "SiO2": [(34.0, 19.0), (25.5, 13.5), (15.0, 8.5), (8.5, 4.1)],
}
NAMES = {"Alumina": "Al2O3", "Copper-oxide": "CuO", "Silica": "SiO2"}


def frames(drawings):
    """Plot frames = paths containing four long axis lines forming a rectangle."""
    out = []
    for d in drawings:
        long = [it for it in d["items"] if it[0] == "l" and
                (abs(it[1].x - it[2].x) > 100 or abs(it[1].y - it[2].y) > 100)]
        if len(long) >= 4:
            xs = sorted({round(v, 1) for it in long for v in (it[1].x, it[2].x)})
            ys = sorted({round(v, 1) for it in long for v in (it[1].y, it[2].y)})
            out.append((xs[0], ys[0], xs[-1], ys[-1]))
    return out


def calibrate(words, fr):
    """Linear maps x_pt -> volume fraction and y_pt -> cp from the tick labels around frame `fr`."""
    x0, y0, x1, y1 = fr
    xs, ys = [], []
    for w in words:
        try:
            val = float(w[4])
        except ValueError:
            continue
        cx, cy = (w[0] + w[2]) / 2, (w[1] + w[3]) / 2
        if x0 - 30 < w[2] <= x0 + 1 and y0 - 6 <= cy <= y1 + 6:        # y-axis labels (left)
            ys.append((cy, val))
        elif y1 <= w[1] <= y1 + 14 and x0 - 8 <= cx <= x1 + 8:          # x-axis labels (below)
            xs.append((cx, val))
    fx = np.polyfit(*zip(*xs), 1)
    fy = np.polyfit(*zip(*ys), 1)
    return (lambda x: np.polyval(fx, x)), (lambda y: np.polyval(fy, y)), len(xs), len(ys)


def markers(drawings, fr):
    """Data markers ('+' of the asterisk symbols) and error bars inside frame `fr`."""
    x0, y0, x1, y1 = fr
    lines = []
    for d in drawings:
        if d.get("dashes") not in (None, "[] 0"):
            continue
        for it in d["items"]:
            if it[0] == "l":
                a, b = it[1], it[2]
                if x0 + 1 < (a.x + b.x) / 2 < x1 and y0 <= (a.y + b.y) / 2 <= y1:
                    lines.append((a.x, a.y, b.x, b.y))
    vert = [l for l in lines if abs(l[0] - l[2]) < 1e-3]
    horiz = [l for l in lines if abs(l[1] - l[3]) < 1e-3]
    out = []
    for v in vert:
        L = abs(v[3] - v[1])
        if not 2.5 < L < 2.9:                         # '+' arm of the marker (2.69 pt)
            continue
        cx, cy = v[0], (v[1] + v[3]) / 2
        if not any(abs((h[0] + h[2]) / 2 - cx) < 0.05 and abs(h[1] - cy) < 0.05 and 2.5 < abs(h[2] - h[0]) < 2.9
                   for h in horiz):
            continue
        # error bar: longest vertical line through cx that brackets the marker
        bars = [b for b in vert if abs(b[0] - cx) < 0.05 and min(b[1], b[3]) <= cy <= max(b[1], b[3])
                and abs(b[3] - b[1]) > 2.9]
        half = max((abs(b[3] - b[1]) / 2 for b in bars), default=0.0)
        out.append((cx, cy, half))
    return out


def main(pdf: str) -> None:
    doc = pymupdf.open(pdf)
    rows = []
    for page in doc:
        words = page.get_text("words")
        # caption blocks: (top y, text, centre x); a caption belongs to the frame above it in the SAME column
        caps = [(b[1], b[4], (b[0] + b[2]) / 2) for b in page.get_text("blocks") if b[4].strip().startswith("Figure")]
        drawings = page.get_drawings()
        for fr in frames(drawings):
            same_col = [c for c in caps if c[0] > fr[3] and fr[0] - 20 <= c[2] <= fr[2] + 20]
            cap = min(same_col, key=lambda c: c[0] - fr[3], default=None)
            if cap is None:
                continue
            m = re.search(r"Figure\s*(\d+):\s*([A-Za-z-]+)-water at (\d+)", cap[1].replace("\n", " "))
            if not m:
                continue
            fig, mat, T_C = int(m.group(1)), NAMES[m.group(2)], float(m.group(3))
            fx, fy, nx, ny = calibrate(words, fr)
            for cx, cy, half in markers(drawings, fr):
                phi, cp = float(fx(cx)), float(fy(cy))
                sd = abs(float(fy(cy - half)) - cp) if half else None
                wt, vol = min(TABLE2[mat], key=lambda t: abs(t[1] / 100 - phi))
                rows.append((mat, fig, T_C + 273.15, phi, vol / 100, wt / 100, cp * 1000,
                             sd * 1000 if sd else None))
    rows.sort(key=lambda r: (r[0], r[2], r[3]))
    out = ROOT / "validation" / "data" / "nanofluid_ohanley2012.csv"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("# SOURCE: O'Hanley, Buongiorno, McKrell, Hu, Adv. Mech. Eng. 2012, 181079, doi:10.1155/2012/181079\n")
        fh.write("#   (open access, CC BY 2.0). DSC data of water-based nanofluids, Figures 1-9.\n")
        fh.write("# Values DIGITISED from the vector figures by tools/digitize_ohanley2012.py (marker centres\n")
        fh.write("#   read from the PDF drawing commands; axes calibrated from tick labels).\n")
        fh.write("# phi_digitised = abscissa of the marker; phi_table2 and w_table2 = nominal concentration from\n")
        fh.write("#   Table 2 of the paper; cp_sd = half-length of the plotted error bar (one standard deviation).\n")
        fh.write("particle,figure,T_K,phi_digitised,phi_table2,w_table2,cp_J_kgK,cp_sd_J_kgK,reference\n")
        for r in rows:
            sd = f"{r[7]:.1f}" if r[7] is not None else ""
            fh.write(f"{r[0]},{r[1]},{r[2]:.2f},{r[3]:.5f},{r[4]:.3f},{r[5]:.3f},{r[6]:.1f},{sd},"
                     f"O'Hanley et al. 2012 Fig. {r[1]}\n")
    print(f"{len(rows)} points written to {out}")


if __name__ == "__main__":
    main(sys.argv[1])
