"""One-off, documented update of data/*.json with values verified on 2026-11-02.

Every change below records its evidence in the "source" string. Kept in the repository so the
team can see exactly what was changed and why (see docs/DATA_SOURCES.md).
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
cf, bf = ROOT / "data" / "components.json", ROOT / "data" / "binary_nrtl.json"
comp = json.loads(cf.read_text(encoding="utf-8"))
binary = json.loads(bf.read_text(encoding="utf-8"))
C = comp["components"]
IAPWS = "IAPWS R10-06(2009) ice Ih equation of state (Feistel & Wagner 2006), implemented in tools/iapws_ice.py and checked against Table 6 of the release to 9 digits"

# ---- ice (water, solid phase) -------------------------------------------------------------
C["water"]["cp"]["solid"] = {
    "type": "poly", "coeffs": [393.92143, 4.3037611, 0.0099635884, -1.0608406e-05], "unit": "J/(kg K)",
    "Tmin": 200.0, "Tmax": 273.16,
    "source": f"Cubic fit (max dev 0.0004 %) to {IAPWS}; 2096.7 J/(kg K) at 273.15 K",
    "verify": False}
C["water"]["density"]["solid"] = {"type": "constant", "value": 916.72, "unit": "kg/m3",
                                  "source": f"{IAPWS}, at 273.15 K and 0.101325 MPa", "verify": False}
C["water"]["alpha"]["solid"] = {"value": 1.5984e-4, "unit": "1/K", "source": f"{IAPWS}, at 273.15 K", "verify": False}
C["water"]["kappa_T"]["solid"] = {"value": 1.1778e-10, "unit": "1/Pa", "source": f"{IAPWS}, at 273.15 K", "verify": False}

# ---- solid densities (CRC values as quoted in the Wikipedia infoboxes) ----------------------
C["Al2O3"]["density"]["solid"].update(value=3987.0, verify=False,
    source="alpha-Al2O3 (corundum) 3.987 g/cm3, CRC Handbook via Wikipedia infobox (checked 2026-11-02). Nanoparticles (often gamma-alumina) may be less dense")
C["CuO"]["density"]["solid"].update(value=6315.0, verify=False,
    source="6.315 g/cm3, CRC Handbook via Wikipedia infobox (checked 2026-11-02)")
C["SiO2"]["density"]["solid"].update(verify=False,
    source="alpha-quartz 2.648 g/cm3 (amorphous silica 2.196 g/cm3), CRC Handbook 92nd ed. via Wikipedia (checked 2026-11-02)")
C["NaCl"]["density"]["solid"].update(value=2170.0, verify=False,
    source="2.17 g/cm3, CRC Handbook via Wikipedia infobox (checked 2026-11-02)")

# ---- binary data ----------------------------------------------------------------------------
for p in binary["nrtl"]:
    if p["i"] == "acetone" and p["j"] == "water":
        # Table convention: row i, column j -> a_ij. Water row/acetone column: a = 0.0544, b = 419.97;
        # acetone row/water column: a = 6.3981, b = -1808.99. The previous entry had i and j swapped.
        p.update(a_ij=6.3981, a_ji=0.0544, b_ij=-1808.991, b_ji=419.9716,
                 source="Aspen Plus databank (DDB-regressed VLE), as tabulated in the NRTL appendix reproduced in BIT Mesra CAPE notes; orientation corrected 2026-11-02 (gamma_inf acetone in water = 5.8, water in acetone = 3.6 at 298 K)",
                 verify=True)
    elif p["i"] in ("ethanol", "methanol") and p["j"] == "water":
        p["verify"] = False
        p["source"] += "; values and orientation confirmed against the published table (2026-11-02)"

for s in binary["solubility"]:
    if s["solute"] == "NaCl":
        s.update(T=[273.15 + 10 * k for k in range(11)],
                 values=[35.65, 35.72, 35.89, 36.09, 36.37, 36.69, 37.04, 37.46, 37.93, 38.47, 38.99],
                 unit="g solute / 100 g solvent",
                 source="CRC Handbook aqueous solubility (mass %, 0-100 C) converted to g/100 g water, as tabulated in Wikipedia 'Solubility table' (checked: 26.28 wt% -> 35.65 g/100 g)",
                 verify=False)
        s["dH_sol"] = {"value": 3.87, "unit": "kJ/mol", "source": "Enthalpy of solution at infinite dilution, 298 K (CRC Handbook via Wikipedia 'Enthalpy change of solution')", "verify": False}

cf.write_text(json.dumps(comp, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
bf.write_text(json.dumps(binary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print("updated", cf.name, bf.name)
