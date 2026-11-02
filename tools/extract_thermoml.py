"""Extract mixture heat-capacity and excess-enthalpy data from the NIST ThermoML archive.

The NIST TRC ThermoML Archive (https://doi.org/10.18434/mds2-2422, file ThermoML.v2020-09-30.tgz)
contains the experimental data of papers published in J. Chem. Eng. Data, J. Chem. Thermodyn.,
Fluid Phase Equilib., Thermochim. Acta and Int. J. Thermophys. (2003-2019) in XML form.

Usage (after downloading and extracting the archive into a folder):
    python tools/extract_thermoml.py <extracted_archive_folder>

Writes validation/data/thermoml_binary_data.csv with one row per data point:
    doi, comp1, comp2, property, x1, T_K, P_kPa, value, unit
for every binary mixture whose two components are both in data/components.json.
"""
from __future__ import annotations

import csv
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from cp_engine.database import Database  # noqa: E402

PROPS = {
    "Excess molar enthalpy (molar enthalpy of mixing), kJ/mol": ("HE", 1000.0, "J/mol"),
    "Excess molar heat capacity, J/K/mol": ("CpE", 1.0, "J/(mol K)"),
    "Molar heat capacity at constant pressure, J/K/mol": ("Cp_molar", 1.0, "J/(mol K)"),
    "Specific heat capacity at constant pressure, J/K/kg": ("cp_mass", 1.0, "J/(kg K)"),
}

# ThermoML common names -> database keys
db = Database(ROOT / "data")
NAME_MAP = {}
for key, comp in db.components.items():
    for n in [key, comp.name, comp.formula] + comp.raw.get("aliases", []):
        NAME_MAP[n.lower()] = key
NAME_MAP.update({"water": "water", "ethanol": "ethanol", "methanol": "methanol", "glycerol": "glycerol",
                 "1,2,3-propanetriol": "glycerol", "propane-1,2,3-triol": "glycerol", "acetone": "acetone",
                 "2-propanone": "acetone", "propan-2-one": "acetone", "hexane": "n-hexane",
                 "toluene": "toluene", "methylbenzene": "toluene", "dodecane": "n-dodecane",
                 "octadecane": "n-octadecane"})


def strip_ns(root):
    for el in root.iter():
        if "}" in el.tag:
            el.tag = el.tag.split("}", 1)[1]
    return root


def text(el, path):
    f = el.find(path)
    return f.text.strip() if f is not None and f.text else None


def parse(path: Path):
    root = strip_ns(ET.parse(path).getroot())
    doi = text(root, "Citation/sDOI") or path.stem
    names = {}
    for c in root.findall("Compound"):
        num = text(c, "RegNum/nOrgNum")
        nm = (text(c, "sCommonName") or "").lower()
        names[num] = NAME_MAP.get(nm)
    for d in root.findall("PureOrMixtureData"):
        comps = [text(c, "RegNum/nOrgNum") for c in d.findall("Component")]
        if len(comps) != 2 or any(names.get(c) is None for c in comps):
            continue
        keys = [names[c] for c in comps]
        props = {}
        for p in d.findall("Property"):
            pname = text(p, ".//ePropName")
            if pname in PROPS:
                props[text(p, "nPropNumber")] = PROPS[pname]
        if not props:
            continue
        # variables: temperature, pressure, mole fraction of a component
        var = {}
        for v in d.findall("Variable"):
            num = text(v, "nVarNumber")
            vt = v.find("VariableID/VariableType")
            kind = vt[0].tag if vt is not None and len(vt) else ""
            label = vt[0].text if vt is not None and len(vt) else ""
            comp = text(v, "VariableID/RegNum/nOrgNum")
            var[num] = (kind, label, comp)
        # constraints (fixed T, P or composition)
        const = {}
        for c in d.findall("Constraint"):
            ct = c.find("ConstraintID/ConstraintType")
            kind = ct[0].tag if ct is not None and len(ct) else ""
            label = ct[0].text if ct is not None and len(ct) else ""
            comp = text(c, "ConstraintID/RegNum/nOrgNum")
            const[(kind, label, comp)] = float(text(c, "nConstraintValue"))
        for nv in d.findall("NumValues"):
            vals = dict(const)
            for vv in nv.findall("VariableValue"):
                vals[var[text(vv, "nVarNumber")]] = float(text(vv, "nVarValue"))
            T = P = x1 = None
            for (kind, label, comp), val in vals.items():
                if kind == "eTemperature":
                    T = val
                elif kind == "ePressure":
                    P = val
                elif kind == "eComponentComposition" and label and label.startswith("Mole fraction"):
                    x1 = val if comp == comps[0] else 1.0 - val
            for pv in nv.findall("PropertyValue"):
                pn = text(pv, "nPropNumber")
                if pn in props and x1 is not None and T is not None:
                    tag, factor, unit = props[pn]
                    yield [doi, keys[0], keys[1], tag, x1, T, P, float(text(pv, "nPropValue")) * factor, unit]


def main(folder: str) -> None:
    rows = []
    files = sorted(Path(folder).rglob("*.xml"))
    for f in files:
        try:
            rows += list(parse(f))
        except ET.ParseError:
            continue
    out = ROOT / "validation" / "data" / "thermoml_binary_data.csv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        fh.write("# SOURCE: NIST TRC ThermoML Archive v2020-09-30 (doi:10.18434/mds2-2422), extracted by\n")
        fh.write("# tools/extract_thermoml.py. x1 = mole fraction of comp1. Cite the original paper (doi column).\n")
        w = csv.writer(fh)
        w.writerow(["doi", "comp1", "comp2", "property", "x1", "T_K", "P_kPa", "value", "unit"])
        w.writerows(rows)
    print(f"{len(files)} files scanned, {len(rows)} points written to {out}")


if __name__ == "__main__":
    main(sys.argv[1])
