"""Fit the temperature dependence of the isothermal compressibility, ln κ_T = c0 + c1 T + c2 T²,
to the NIST SRD 69 tables at 0.101325 MPa (κ_T = (Cp/Cv)/(ρ w²)) and store it in data/components.json.

κ_T(T) is needed by the Tait equation used for the pressure dependence of Cp (see
cp_engine/database.py, Component.cp_pressure_correction).

Run from the project root:  python tools/fit_kappa_nist.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
FLUIDS = ["water", "methanol", "n-hexane", "toluene", "n-dodecane"]

path = ROOT / "data" / "components.json"
data = json.loads(path.read_text(encoding="utf-8"))
for f in FLUIDS:
    t = pd.read_csv(ROOT / "validation" / "data" / f"nist_{f}_1atm.csv")
    k = t.cp_J_kgK / t.cv_J_kgK / (t.rho_kg_m3 * t.sound_speed_m_s ** 2)
    c = np.polyfit(t.T_K, np.log(k), 2)[::-1]
    err = float(np.max(np.abs(np.exp(np.polyval(c[::-1], t.T_K)) - k) / k) * 100)
    data["components"][f]["kappa_T"]["liquid"] = {
        "type": "ln_poly", "coeffs": [float(f"{v:.8g}") for v in c], "unit": "1/Pa",
        "Tmin": float(t.T_K.min()), "Tmax": float(t.T_K.max()),
        "source": f"Fit (tools/fit_kappa_nist.py) of ln(kappa_T) to NIST SRD 69 at 0.101325 MPa, "
                  f"kappa_T = (Cp/Cv)/(rho w^2); max dev {err:.2f} %",
        "verify": False}
    print(f"{f}: kappa_T(298.15 K) = {np.exp(np.polyval(c[::-1], 298.15)):.4e} 1/Pa, fit max dev {err:.2f} %")
path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
