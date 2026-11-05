"""Replace the remaining approximate n-octadecane (PCM) properties with fits to measured data.

Data (NIST ThermoML archive, doi:10.18434/mds2-2422):
  * solid Cp: Fonseca, Gushterov & Dohrn, J. Chem. Thermodyn. 2014 (doi:10.1016/j.jct.2013.11.038),
    283.85-290.46 K, plus Messerly et al. (1967) at 298.15 K (NIST WebBook) – Cp rises toward T_m (pre-melting)
  * liquid specific volume: Espeau & Ceolin, Thermochim. Acta 2006 (doi:10.1016/j.tca.2006.04.001), 299-568 K
  * liquid compressibility check: Kariznovi et al. 2012 (doi:10.1016/j.fluid.2011.10.018) and
    Wu et al. 2011 (doi:10.1016/j.fluid.2011.08.020), densities at 323 K and 2-14 MPa
  * T_m and ΔH_fus confirmations: Wei et al. 2013/2014/2016, Fonseca et al. 2014, Parsa et al. 2015

Run:  python tools/fit_octadecane.py
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent

SOLID_CP = [(283.85, 448.6), (285.5, 455.2), (287.15, 460.7), (288.81, 467.3), (290.46, 472.4),
            (298.15, 485.64)]   # J/(mol K); last point Messerly et al. 1967 (NIST WebBook)
LIQ_V = [(299.25, 0.0012896), (301.35, 0.0012921), (303.35, 0.0012945), (305.45, 0.0012949), (307.55, 0.0012984),
         (309.55, 0.0013007), (311.65, 0.0013014), (313.65, 0.0013042), (315.75, 0.0013076), (317.75, 0.0013113),
         (320.65, 0.0013186), (323.55, 0.0013231), (325.45, 0.001324), (328.35, 0.0013296), (330.35, 0.0013306),
         (333.25, 0.0013349), (335.15, 0.0013375), (338.05, 0.0013414), (341.05, 0.001345), (342.95, 0.0013459),
         (345.85, 0.0013504), (347.85, 0.0013543), (350.75, 0.0013566), (352.65, 0.0013597), (355.55, 0.0013627),
         (357.55, 0.0013644), (360.45, 0.0013682), (362.35, 0.0013715), (365.25, 0.0013745), (367.25, 0.0013768),
         (370.15, 0.0013822), (372.05, 0.0013846), (374.95, 0.0013883), (376.95, 0.0013906), (379.85, 0.0013949),
         (381.75, 0.0013977), (384.75, 0.001402), (386.65, 0.0014039), (389.55, 0.0014093), (391.55, 0.0014115),
         (394.45, 0.0014158), (396.35, 0.0014192), (399.25, 0.0014227)]     # m^3/kg, subset 299-400 K
KAPPA_CHECK = [(2.08e6, 763.0), (3.57e6, 764.0), (5.06e6, 765.0), (6.50e6, 766.0), (8.09e6, 766.0), (9.49e6, 767.0)]


def main() -> None:
    T, cp = np.array(SOLID_CP).T
    c_s = np.polyfit(T, cp, 1)[::-1]
    err_s = np.max(np.abs(np.polyval(c_s[::-1], T) - cp) / cp) * 100
    T, v = np.array(LIQ_V).T
    rho = 1.0 / v
    c_l = np.polyfit(T, rho, 2)[::-1]
    err_l = np.max(np.abs(np.polyval(c_l[::-1], T) - rho) / rho) * 100
    P, r = np.array(KAPPA_CHECK).T
    kappa = np.polyfit(P, r, 1)[0] / r.mean()

    path = ROOT / "data" / "components.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    oc = data["components"]["n-octadecane"]
    oc["cp"]["solid"] = {"type": "poly", "coeffs": [float(f"{x:.8g}") for x in c_s], "unit": "J/(mol K)",
                         "Tmin": 280.0, "Tmax": 301.3,
                         "source": "Linear fit (tools/fit_octadecane.py) to Fonseca et al. 2014 (doi:10.1016/j.jct.2013.11.038, "
                                   f"283.85-290.46 K) and Messerly et al. 1967 at 298.15 K; max dev {err_s:.1f} %; "
                                   "rises toward T_m (pre-melting)", "verify": False}
    oc["density"]["liquid"] = {"type": "poly", "coeffs": [float(f"{x:.8g}") for x in c_l], "Tmin": 299.0, "Tmax": 400.0,
                               "source": "Quadratic fit (tools/fit_octadecane.py) to the liquid specific volume of Espeau & "
                                         f"Ceolin 2006 (doi:10.1016/j.tca.2006.04.001), 299-400 K; max dev {err_l:.2f} %",
                               "verify": False}
    oc.get("alpha", {}).pop("liquid", None)            # now follows from the density model
    oc["kappa_T"]["liquid"].update(
        value=9.0e-10, verify=False,
        source=f"Consistent with densities at 323 K, 2-9.5 MPa (Kariznovi et al. 2012, doi:10.1016/j.fluid.2011.10.018): "
               f"{kappa:.2e} 1/Pa (densities given to 1 kg/m3, so about ±25 %)")
    oc["Tm"]["source"] += "; confirmed by 300.95-301.6 K in Wei et al. 2013/2014/2016, Fonseca et al. 2014, Parsa et al. 2015 (ThermoML)"
    oc["dHfus"]["source"] += "; confirmed by 61.34-61.65 kJ/mol (Fonseca et al. 2014; Wei et al. 2013, 2014, 2016)"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"solid Cp fit max dev {err_s:.2f} %, coeffs {c_s}")
    print(f"liquid density fit max dev {err_l:.3f} %, rho(301.3 K) = {np.polyval(c_l[::-1], 301.3):.1f} kg/m3")
    print(f"kappa_T from density vs P at 323 K: {kappa:.2e} 1/Pa")


if __name__ == "__main__":
    main()
