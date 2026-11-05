"""Extract phase-equilibrium data (VLE, activity coefficients at infinite dilution, LLE) for selected
pairs from the NIST ThermoML archive (doi:10.18434/mds2-2422).

Usage:  python tools/extract_thermoml_phase_eq.py <extracted_archive_folder>

Writes to validation/data/:
  thermoml_vle_ethanol_toluene.csv  isobaric T-x-y of ethanol + toluene (Chen et al. 2012) and the pure
                                     boiling points from the same paper (x = 0 or 1)
  thermoml_gamma_inf.csv             activity coefficients at infinite dilution in water
  thermoml_lle.csv                   liquid-liquid coexistence data (methanol + n-hexane, + n-dodecane)
"""
from __future__ import annotations

import csv
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "validation" / "data"
NAMES = {"ethanol": "ethanol", "toluene": "toluene", "methylbenzene": "toluene", "water": "water",
         "acetone": "acetone", "2-propanone": "acetone", "propan-2-one": "acetone", "methanol": "methanol",
         "hexane": "n-hexane", "n-hexane": "n-hexane", "dodecane": "n-dodecane", "n-dodecane": "n-dodecane"}


def load(folder: Path, doi: str):
    r = ET.parse(folder / f"{doi}.xml").getroot()
    for el in r.iter():
        if "}" in el.tag:
            el.tag = el.tag.split("}", 1)[1]
    names = {c.findtext("RegNum/nOrgNum"): NAMES.get((c.findtext("sCommonName") or "").lower())
             for c in r.findall("Compound")}
    return r, names


def records(d, names):
    """Yield dicts {variable/property label: value} for every data point of a PureOrMixtureData block."""
    var = {v.findtext("nVarNumber"): (v.find("VariableID/VariableType")[0].text, names.get(v.findtext("VariableID/RegNum/nOrgNum")))
           for v in d.findall("Variable")}
    prop = {p.findtext("nPropNumber"): (p.findtext(".//ePropName"), names.get(p.findtext(".//RegNum/nOrgNum")))
            for p in d.findall("Property")}
    const = {(c.find("ConstraintID/ConstraintType")[0].text, names.get(c.findtext("ConstraintID/RegNum/nOrgNum"))):
             float(c.findtext("nConstraintValue")) for c in d.findall("Constraint")}
    for nv in d.findall("NumValues"):
        rec = dict(const)
        for v in nv.findall("VariableValue"):
            rec[var[v.findtext("nVarNumber")]] = float(v.findtext("nVarValue"))
        for p in nv.findall("PropertyValue"):
            rec[("PROP",) + prop[p.findtext("nPropNumber")]] = float(p.findtext("nPropValue"))
        yield rec


def comps_of(d, names):
    return [names.get(c.findtext("RegNum/nOrgNum")) for c in d.findall("Component")]


def main(folder: str) -> None:
    folder = Path(folder)
    # ---- VLE ethanol + toluene (Chen, Zhong & Xu, J. Chem. Eng. Data 2012, 57, 155) ---------------
    doi = "10.1021/je200921u"
    r, names = load(folder, doi)
    Tb, y = {}, {}
    for d in r.findall("PureOrMixtureData"):
        cs = comps_of(d, names)
        for rec in records(d, names):
            if None not in cs and sorted(cs) == ["ethanol", "toluene"]:
                x, P = rec[("Mole fraction", "ethanol")], rec[("Pressure, kPa", None)]
                if ("PROP", "Boiling temperature at pressure P, K", None) in rec:
                    Tb[(x, P)] = rec[("PROP", "Boiling temperature at pressure P, K", None)]
                if ("PROP", "Mole fraction", "ethanol") in rec:
                    y[(x, P)] = rec[("PROP", "Mole fraction", "ethanol")]
    with open(OUT / "thermoml_vle_ethanol_toluene.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write(f"# SOURCE: Chen, Zhong & Xu, J. Chem. Eng. Data 2012, 57, 155 (doi:{doi}), via NIST ThermoML archive.\n"
                 "# Isobaric VLE: liquid mole fraction x and vapour mole fraction y of ethanol, boiling temperature T.\n")
        w = csv.writer(fh)
        w.writerow(["P_kPa", "x_ethanol", "y_ethanol", "T_K", "doi"])
        for (x, P) in sorted(Tb, key=lambda k: (k[1], k[0])):
            w.writerow([P, x, y.get((x, P), ""), Tb[(x, P)], doi])
    print(f"VLE ethanol + toluene: {len(Tb)} points")

    # ---- activity coefficients at infinite dilution in water -----------------------------------------
    rows = []
    for doi in ("10.1021/je0495942", "10.1021/je049875+"):
        r, names = load(folder, doi)
        for d in r.findall("PureOrMixtureData"):
            cs = comps_of(d, names)
            if len(cs) != 2 or "water" not in cs or None in cs:
                continue
            solute = [c for c in cs if c != "water"][0]
            for rec in records(d, names):
                key = ("PROP", "Activity coefficient", solute)
                if key in rec:
                    rows.append([solute, "water", rec[("Temperature, K", None)], rec[key], doi])
    with open(OUT / "thermoml_gamma_inf.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write("# SOURCE: NIST ThermoML archive: Straver & de Loos, J. Chem. Eng. Data 2005, 50, 1171 (doi:10.1021/je0495942);\n"
                 "#   Atik, Gruber, Krummen & Gmehling, J. Chem. Eng. Data 2004, 49, 1429 (doi:10.1021/je049875+).\n"
                 "# Activity coefficient of the solute at infinite dilution in water.\n")
        w = csv.writer(fh)
        w.writerow(["solute", "solvent", "T_K", "gamma_inf", "doi"])
        w.writerows(rows)
    print(f"gamma_inf: {len(rows)} points")

    # ---- LLE ------------------------------------------------------------------------------------------
    rows = []
    r, names = load(folder, "10.1021/je0505321")            # Trejo et al. 2006: coexistence curve
    for d in r.findall("PureOrMixtureData"):
        if comps_of(d, names) == ["methanol", "n-hexane"]:
            for rec in records(d, names):
                rows.append(["methanol", "n-hexane", rec[("PROP", "Liquid-liquid equilibrium temperature, K", None)],
                             rec[("Mole fraction", "methanol")], "coexistence curve", "10.1021/je0505321"])
    r, names = load(folder, "10.1016/j.jct.2003.11.009")    # Casas et al. 2004: phase compositions
    for d in r.findall("PureOrMixtureData"):
        if comps_of(d, names) == ["methanol", "n-dodecane"]:
            for rec in records(d, names):
                rows.append(["methanol", "n-dodecane", rec[("Temperature, K", None)],
                             rec[("PROP", "Mole fraction", "methanol")], "phase composition", "10.1016/j.jct.2003.11.009"])
    with open(OUT / "thermoml_lle.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write("# SOURCE: NIST ThermoML archive: Trejo, Yanez & Eustaquio-Rincon, J. Chem. Eng. Data 2006, 51, 1070\n"
                 "#   (doi:10.1021/je0505321); Casas, Orge, Diaz & Tojo, J. Chem. Thermodyn. 2004, 36, 237\n"
                 "#   (doi:10.1016/j.jct.2003.11.009). Liquid-liquid equilibrium: T and methanol mole fraction.\n")
        w = csv.writer(fh)
        w.writerow(["comp1", "comp2", "T_K", "x1", "kind", "doi"])
        w.writerows(rows)
    print(f"LLE: {len(rows)} points")


if __name__ == "__main__":
    main(sys.argv[1])
