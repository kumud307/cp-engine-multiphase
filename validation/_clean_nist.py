"""Convert raw NIST WebBook fluid tables (0.101325 MPa isobars) into clean CSVs
and print the derived density fits and isothermal compressibilities used in
data/components.json.  Source: NIST Chemistry WebBook, SRD 69, Thermophysical
Properties of Fluid Systems (Lemmon, McLinden & Friend), downloaded 2026-11-02.
"""
import glob, os
import numpy as np, pandas as pd
here = os.path.dirname(os.path.abspath(__file__))
for raw in sorted(glob.glob(os.path.join(here, "data", "nist_raw_*.tsv"))):
    if os.path.basename(raw).startswith("nist_raw_iso_"):
        continue
    name = os.path.basename(raw)[9:-4]
    df = pd.read_csv(raw, sep="\t")
    df = df[df["Phase"] == "liquid"]
    out = pd.DataFrame({
        "T_K": df["Temperature (K)"],
        "rho_kg_m3": df["Density (kg/m3)"],
        "cp_J_kgK": df["Cp (J/g*K)"] * 1000.0,
        "cv_J_kgK": df["Cv (J/g*K)"] * 1000.0,
        "sound_speed_m_s": df["Sound Spd. (m/s)"],
    })
    out.to_csv(os.path.join(here, "data", f"nist_{name}_1atm.csv"), index=False)
    T, rho = out.T_K.values, out.rho_kg_m3.values
    c = np.polyfit(T, rho, 2)                        # rho = c0 T^2 + c1 T + c2
    fit_err = np.max(np.abs(np.polyval(c, T) - rho) / rho) * 100
    i = np.argmin(abs(T - 300.0))
    g = out.cp_J_kgK.values[i] / out.cv_J_kgK.values[i]
    kT = g / (rho[i] * out.sound_speed_m_s.values[i] ** 2)
    print(f"{name:11s} T={T.min():.0f}-{T.max():.0f} rho2={c[0]:.6e} rho1={c[1]:.6e} rho0={c[2]:.6e} "
          f"maxerr={fit_err:.3f}%  kappa_T(300K)={kT:.4e} 1/Pa")


# ---- isotherms at 298.15 K, 0.1-100 MPa (pressure dependence of Cp) ----------------------------
for raw in sorted(glob.glob(os.path.join(here, "data", "nist_raw_iso_*.tsv"))):
    name = os.path.basename(raw)[13:-4]
    df = pd.read_csv(raw, sep="	")
    df = df[df["Phase"] == "liquid"]
    pd.DataFrame({"T_K": df["Temperature (K)"], "P_MPa": df["Pressure (MPa)"],
                  "rho_kg_m3": df["Density (kg/m3)"], "cp_J_kgK": df["Cp (J/g*K)"] * 1000.0}).to_csv(
        os.path.join(here, "data", f"nist_{name}_298K_isotherm.csv"), index=False)
    print(f"{name}: isotherm {len(df)} points")
