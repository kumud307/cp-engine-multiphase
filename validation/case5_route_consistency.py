"""Validation case 5 – calorimetric-route self-consistency.

For each system: Route A gives ΔH over [T1, T2]; that ΔH is fed back as a
"measured" heat Q at constant P (Route B) and also as ΔS / Cv-type inputs.
Route B must return exactly the Route-A mean Cp.  This proves that the
first-law bookkeeping, unit handling and route reconstruction are consistent.
"""
from __future__ import annotations

import pandas as pd

from common import aard, row
from cp_engine import CpEngine

SYSTEMS = {
    "water (pure)": ({"components": ["water"]}, 290.0, 310.0),
    "water + ethanol 40 wt%": ({"mass_fractions": {"water": 0.6, "ethanol": 0.4}}, 290.0, 310.0),
    "methanol + water + ethanol + glycerol": (
        {"mass_fractions": {"water": 0.4, "methanol": 0.2, "ethanol": 0.2, "glycerol": 0.2}}, 295.0, 315.0),
    "CuO–water nanofluid φ = 4 %": ({"components": ["water", "CuO"], "phi": {"CuO": 0.04}}, 300.0, 320.0),
    "water–toluene–ethanol emulsion + SiO2": (
        {"mass_fractions": {"water": 0.45, "toluene": 0.3, "ethanol": 0.15, "SiO2": 0.1}}, 300.0, 320.0),
    "ice slurry 10 wt% ethanol": ({"mass_fractions": {"water": 0.9, "ethanol": 0.1}}, 255.0, 270.0),
    "PCM slurry 20 wt% n-octadecane": ({"mass_fractions": {"water": 0.8, "n-octadecane": 0.2}}, 295.0, 308.0),
}


def run(save: bool = True) -> list[dict]:
    engine = CpEngine()
    rows = []
    for name, (comp, T1, T2) in SYSTEMS.items():
        a = engine.solve(dict(comp, T1=T1, T2=T2, total_mass="1 kg"))
        rA = [r for r in a.routes if r.tag == "A:mean"][0]
        Q = rA.breakdown.extras["dH_J"] / rA.breakdown.mass          # J for 1 kg
        b = engine.solve({"T1": T1, "T2": T2, "total_mass": "1 kg", "process": "constant_P", "Q": (Q, "J")})
        err, _, n = aard([b.cp_mass], [rA.cp_mass])
        rows.append(row("5", name, "Route B (Q) vs Route A (ΔH_model)", n, err, err, "self-consistency",
                        metric="rel. dev %"))
    return rows


if __name__ == "__main__":
    print(pd.DataFrame(run()).to_string(index=False))
