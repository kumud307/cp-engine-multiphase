"""Download NIST SRD 69 fluid tables used for vapour pressures, critical points and real-gas Cp.

Run from the project root (needs internet and curl):  python tools/fetch_nist_fluid.py

Writes, for water, methanol, n-hexane, toluene and n-dodecane (the liquids of the database that
have a reference equation of state in the NIST Chemistry WebBook):
  validation/data/nist_sat_<name>.csv          saturation curve (T, P, liquid and vapour Cp)
  validation/data/nist_<name>_vapour_<P>MPa.csv superheated-vapour isobars (real-gas Cp)
Source: Lemmon, McLinden & Friend, "Thermophysical Properties of Fluid Systems",
NIST Chemistry WebBook, NIST SRD 69, doi:10.18434/T4D303.
"""
from __future__ import annotations

import io
import subprocess
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent.parent / "validation" / "data"
URL = ("https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID={id}&Type={type}&Digits=6"
       "&{range}&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs"
       "&VisUnit=uPa*s&STUnit=N%2Fm")

# NIST ID, saturation range (K), superheated isobars: (P / MPa, T_low, T_high)
FLUIDS = {
    "water": ("C7732185", (275, 645), [(1.0, 460, 900), (5.0, 540, 900), (10.0, 590, 900)]),
    "methanol": ("C67561", (250, 510), [(1.0, 420, 620), (2.0, 450, 620)]),
    "n-hexane": ("C110543", (250, 505), [(1.0, 450, 650), (2.0, 490, 650)]),
    "toluene": ("C108883", (260, 590), [(1.0, 520, 700)]),
    "n-dodecane": ("C112403", (300, 655), [(0.5, 600, 750)]),
}


def get(fid: str, type_: str, rng: str) -> pd.DataFrame:
    text = subprocess.run(["curl", "-s", URL.format(id=fid, type=type_, range=rng)],
                          capture_output=True, text=True, check=True, timeout=120).stdout
    d = pd.read_csv(io.StringIO(text), sep="\t")
    for c in d.columns:
        if c != "Phase":
            d[c] = pd.to_numeric(d[c], errors="coerce")     # "undefined" near the critical point -> NaN
    return d


def main() -> None:
    for name, (fid, (t0, t1), isobars) in FLUIDS.items():
        d = get(fid, "SatT", f"THigh={t1}&TLow={t0}&TInc=5")
        sat = pd.DataFrame({"T_K": d["Temperature (K)"], "P_MPa": d["Pressure (MPa)"],
                            "rho_liq_kg_m3": d["Density (l, kg/m3)"], "cp_liq_J_kgK": d["Cp (l, J/g*K)"] * 1e3,
                            "cp_vap_J_kgK": d["Cp (v, J/g*K)"] * 1e3})
        with open(OUT / f"nist_sat_{name}.csv", "w", encoding="utf-8", newline="") as f:
            f.write(f"# SOURCE: NIST Chemistry WebBook SRD 69 (doi:10.18434/T4D303), saturation table, {fid}\n")
            sat.to_csv(f, index=False)
        for P, a, b in isobars:
            d = get(fid, "IsoBar", f"P={P}&THigh={b}&TLow={a}&TInc=10")
            d = d[d["Phase"] == "vapor"]
            v = pd.DataFrame({"T_K": d["Temperature (K)"], "P_MPa": d["Pressure (MPa)"],
                              "rho_kg_m3": d["Density (kg/m3)"], "cp_J_kgK": d["Cp (J/g*K)"] * 1e3})
            with open(OUT / f"nist_{name}_vapour_{P:g}MPa.csv", "w", encoding="utf-8", newline="") as f:
                f.write(f"# SOURCE: NIST Chemistry WebBook SRD 69 (doi:10.18434/T4D303), isobar {P} MPa, vapour phase\n")
                v.to_csv(f, index=False)
        print(f"{name}: {len(sat)} saturation points, {len(isobars)} vapour isobars")


if __name__ == "__main__":
    main()
