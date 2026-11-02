"""Fit the apparent molar heat capacity of NaCl in water and store it in data/binary_nrtl.json.

For an electrolyte solution the heat capacity is written as
    Cp(solution) = n_w Cp_w°(T) + n_s φCp(m, T)
where φCp is the apparent molar heat capacity of the dissolved salt. Debye–Hückel theory gives a
√m dependence at low molality, so the empirical form used here is
    φCp = c0 + d0 √m + e·m        (298.15 K; m in mol/kg)

Data (NIST ThermoML archive), three papers report Cp of NaCl(aq):
  * doi:10.1021/je0502400 – a dedicated study of salt-solution heat capacities, m = 1, 2, 3 mol/kg
    at 298.15 K: USED FOR THE FIT.
  * doi:10.1016/j.jct.2016.09.031 – m = 2.05 mol/kg at 298.15 K: kept as the INDEPENDENT TEST.
  * doi:10.1016/j.jct.2018.10.027 – NaCl(aq) reported incidentally in a study of heterocyclic compounds
    in amino-acid buffers; its values imply φCp ≈ +17 J/(mol K) at 1 mol/kg against −40 in the
    dedicated study (a 1.1 % difference in Cp), i.e. they are inconsistent: EXCLUDED, with this reason.
The temperature dependence of φCp is therefore not available; the engine warns away from 298 K.

Run:  python tools/fit_apparent_cp.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from cp_engine.database import Database  # noqa: E402

FIT_DOI = "10.1021/je0502400"
EXCLUDED_DOI = "10.1016/j.jct.2018.10.027"


def nacl_data() -> pd.DataFrame:
    db = Database(ROOT / "data")
    d = pd.read_csv(ROOT / "validation" / "data" / "thermoml_binary_data.csv", comment="#")
    n = d[(((d.comp1 == "NaCl") & (d.comp2 == "water")) | ((d.comp1 == "water") & (d.comp2 == "NaCl")))
          & (d.property == "Cp_molar")].copy()
    n["x"] = np.where(n.comp1 == "NaCl", n.x1, 1 - n.x1)
    n["m"] = n.x / (1 - n.x) / db.get("water").MW
    cpw = np.array([db.get("water").cp("liquid", T) for T in n.T_K])
    n["phiCp"] = (n.value - (1 - n.x) * cpw) / n.x             # J/(mol K) of NaCl
    return n.drop_duplicates(subset=["doi", "T_K", "m"])


def design(T, m):
    m = np.asarray(m, dtype=float)
    return np.column_stack([np.ones_like(m), np.sqrt(m), m])


def main() -> None:
    n = nacl_data()
    fit = n[n.doi == FIT_DOI]
    X = design(fit.T_K, fit.m)
    coef, *_ = np.linalg.lstsq(X, fit.phiCp, rcond=None)
    rmse = float(np.sqrt(np.mean((X @ coef - fit.phiCp) ** 2)))
    path = ROOT / "data" / "binary_nrtl.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    data["apparent_molar_cp"] = [{
        "solute": "NaCl", "solvent": "water",
        "form": "phiCp = c0 + d0 sqrt(m) + e m ; m in mol/kg, J/(mol K), at 298.15 K (T-dependence not available)",
        "coeffs": {"c0": coef[0], "d0": coef[1], "e": coef[2]},
        "T_ref": 298.15, "m_max": float(fit.m.max()), "excluded": EXCLUDED_DOI,
        "rmse_phiCp_J_per_molK": round(rmse, 2),
        "source": f"Fitted (tools/fit_apparent_cp.py) to Cp of NaCl(aq), doi:{FIT_DOI} (NIST ThermoML archive), "
                  f"{len(fit)} points; water Cp from this database",
        "verify": False}]
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"fitted on {len(fit)} points: phiCp RMSD {rmse:.2f} J/(mol K); phiCp(298 K, m->0) = {coef[0]:.1f} J/(mol K)")
    test = n[(n.doi != FIT_DOI) & (n.doi != EXCLUDED_DOI)]
    pred = design(test.T_K, test.m) @ coef
    for (doi, T, m, p, o) in zip(test.doi, test.T_K, test.m, pred, test.phiCp):
        print(f"  test {doi}: T={T} m={m:.2f}: phiCp model {p:.1f} vs data {o:.1f}")


if __name__ == "__main__":
    main()
