"""Vapour pressures and critical constants for the volatile liquids of the database.

Run from the project root:  python tools/fit_vapour_pressure.py   (after tools/fetch_nist_fluid.py)

* Antoine equations  log10(P/bar) = A − B/(T + C)  fitted to the NIST SRD 69 saturation tables of
  water, methanol, n-hexane and toluene (two temperature ranges each, as NIST does for measured data).
* Critical temperature, pressure and volume: the critical point of the NIST reference equation of
  state (last row of the saturation table).  Ethanol and acetone have no NIST reference EOS:
  ethanol uses the measurement of doi:10.1016/j.fluid.2018.04.021 (514.0 K, 6.141 MPa; confirmed by
  doi:10.1016/j.jct.2013.03.002, 514.2 K and 6.135 MPa); acetone uses the NIST WebBook average
  (508 ± 2 K, 48 ± 4 bar), flagged "verify".
* Acentric factor from its definition (Pitzer): ω = −1 − log10[P_sat(T = 0.7 T_c) / P_c].

The results are written into data/components.json ("antoine" and "critical" entries).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "validation" / "data"
NIST = "NIST SRD 69 reference equation of state (saturation table, doi:10.18434/T4D303)"

# fitting ranges (K) for the two Antoine sets of each fluid
RANGES = {"water": [(273.2, 400.0), (390.0, 500.0), (490.0, 620.0)],
          "methanol": [(250.0, 360.0), (350.0, 490.0)],
          "n-hexane": [(250.0, 360.0), (350.0, 490.0)],
          "toluene": [(380.0, 425.0)]}              # fills the gap between the two measured NIST sets


def antoine(T, A, B, C):
    return A - B / (T + C)


def sat(name: str) -> pd.DataFrame:
    return pd.read_csv(DATA / f"nist_sat_{name}.csv", comment="#").dropna(subset=["P_MPa"])


def main() -> None:
    path = ROOT / "data" / "components.json"
    db = json.loads(path.read_text(encoding="utf-8"))
    comps = db["components"]
    MW = {k: comps[k]["MW"]["value"] for k in comps}

    for name, ranges in RANGES.items():
        d = sat(name)
        sets = [a for a in comps[name].get("antoine", []) if "NIST SRD 69" not in a.get("source", "")]
        for lo, hi in ranges:
            s = d[(d.T_K >= lo) & (d.T_K <= hi)]
            p, _ = curve_fit(antoine, s.T_K, np.log10(s.P_MPa * 10), p0=(4.5, 1500.0, -40.0), maxfev=20000)
            err = 100 * np.max(np.abs(10 ** antoine(s.T_K, *p) / (s.P_MPa * 10) - 1))
            sets.append({"A": round(p[0], 6), "B": round(p[1], 4), "C": round(p[2], 4), "Tmin": lo, "Tmax": hi,
                         "unit": "bar", "source": f"Fitted to the {NIST}; max deviation {err:.2f} %", "verify": False})
            print(f"{name:9s} {lo:.0f}-{hi:.0f} K: A={p[0]:.5f} B={p[1]:.3f} C={p[2]:.3f}  max dev {err:.3f} %")
        comps[name]["antoine"] = sorted(sets, key=lambda a: a["Tmin"])

    # acetone: NIST WebBook Antoine set (Ambrose, Sprake & Townsend 1974), 259.16–507.60 K
    comps["acetone"]["antoine"] = [{"A": 4.42448, "B": 1312.253, "C": -32.445, "Tmin": 259.16, "Tmax": 507.60,
                                    "unit": "bar", "source": "NIST WebBook (Ambrose, Sprake & Townsend 1974)",
                                    "verify": False}]

    def psat_db(name: str, T: float) -> float:            # bar
        sets = comps[name]["antoine"]
        inside = [a for a in sets if a["Tmin"] <= T <= a["Tmax"]] or sets
        a = inside[0]
        return 10 ** (a["A"] - a["B"] / (T + a["C"]))

    for name in ("water", "methanol", "n-hexane", "toluene", "n-dodecane"):
        d = pd.read_csv(DATA / f"nist_sat_{name}.csv", comment="#")
        crit = d.iloc[-1]
        Tc, Pc, rhoc = float(crit.T_K), float(crit.P_MPa) * 1e6, float(crit.rho_liq_kg_m3)
        s = d.dropna(subset=["P_MPa"])
        p07 = float(np.exp(np.interp(0.7 * Tc, s.T_K, np.log(s.P_MPa)))) * 1e6
        omega = -1.0 - np.log10(p07 / Pc)
        Vc = MW[name] / 1000.0 / rhoc
        comps[name]["critical"] = {"Tc": round(Tc, 3), "Pc": round(Pc, 0), "Vc": round(Vc * 1e6, 3), "omega": round(omega, 4),
                                   "units": "K, Pa, cm3/mol", "source": f"Critical point of the {NIST}; ω from P_sat at 0.7 T_c",
                                   "verify": False}
        print(f"{name:10s} Tc={Tc:.2f} K Pc={Pc / 1e6:.4f} MPa Vc={Vc * 1e6:.1f} cm3/mol omega={omega:.4f}")

    for name, Tc, Pc, src, verify in (
            ("ethanol", 514.0, 6.141e6, "Measured, doi:10.1016/j.fluid.2018.04.021 (514.0 ± 0.4 K, 6141 ± 51 kPa); "
             "doi:10.1016/j.jct.2013.03.002 gives 514.2 K, 6135 kPa", False),
            ("acetone", 508.0, 4.8e6, "NIST WebBook average of measured values (508 ± 2 K, 48 ± 4 bar)", True)):
        omega = -1.0 - np.log10(psat_db(name, 0.7 * Tc) * 1e5 / Pc)
        comps[name]["critical"] = {"Tc": Tc, "Pc": Pc, "Vc": None, "omega": round(float(omega), 4),
                                   "units": "K, Pa, cm3/mol", "source": src + "; ω from the NIST Antoine equation at 0.7 T_c",
                                   "verify": verify}
        print(f"{name:10s} Tc={Tc} K Pc={Pc / 1e6} MPa omega={omega:.4f}")

    path.write_text(json.dumps(db, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("written", path)


if __name__ == "__main__":
    main()
