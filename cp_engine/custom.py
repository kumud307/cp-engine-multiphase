"""Components and pure-component Cp values supplied by the user at run time.

The problem statement asks for mixture properties "using the pure component properties", so a user
must be able to bring their own: either a whole new component ('custom_components') or their own Cp
for a database component ('pure_cp'). Both create a modified COPY of the database for that one
calculation; data/components.json is never changed.

    engine.solve({"T": "25 C", "mass_fractions": {"water": 0.5, "ethylene glycol": 0.5},
                  "custom_components": {"ethylene glycol": {"MW": "62.07 g/mol", "cp": "2.43 kJ/(kg K)",
                                                            "density": "1110 kg/m3", "miscibility_class": "polar"}}})
    engine.solve({..., "pure_cp": {"water": "4.18 kJ/(kg K)"}})
"""
from __future__ import annotations

import copy

from .database import CpCorrelation, Database
from .errors import DataError, InputError
from .report import Log
from .units import to_si

SRC = {"source": "user input", "verify": False}


def custom_to_schema(name: str, spec: dict) -> dict:
    """Turn a short user specification into the components.json schema.

    Keys (quantities with units): MW (required), cp (required: a quantity, {phase: quantity}, or a full
    correlation dict such as {"type": "dippr100", ...}), phase (default 'liquid'), density, alpha,
    kappa_T, Tm, dHfus, Tb, dHvap, miscibility_class ('aqueous', 'polar', 'nonpolar', ...), aliases.
    """
    if not isinstance(spec, dict) or "MW" not in spec or "cp" not in spec:
        raise InputError(f"custom component '{name}' needs at least 'MW' and 'cp'.")
    MW, _, _ = to_si(spec["MW"], "molar_mass", name=f"{name}.MW")
    phase = str(spec.get("phase", "liquid"))

    def cp_entry(q):
        if isinstance(q, dict) and "type" in q:
            return dict(SRC, **q)
        v, dim, _ = to_si(q, ["entropy_per_mass", "entropy_per_amount"], name=f"{name}.cp")
        return dict(SRC, type="constant", value=v, unit="J/(kg K)" if dim == "entropy_per_mass" else "J/(mol K)")

    cp_spec = spec["cp"]
    if isinstance(cp_spec, dict) and "type" not in cp_spec and "value" not in cp_spec:
        cp = {ph: cp_entry(q) for ph, q in cp_spec.items()}
    else:
        cp = {phase: cp_entry(cp_spec)}
    out = {"name": spec.get("name", name), "formula": spec.get("formula", ""),
           "aliases": list(spec.get("aliases", [])), "miscibility_class": spec.get("miscibility_class", "unknown"),
           "MW": dict(SRC, value=MW, unit="kg/mol"), "cp": cp}
    for key in ("Tm", "Tb"):
        if key in spec:
            out[key] = dict(SRC, value=to_si(spec[key], "temperature", absolute_temperature=True, name=key)[0], unit="K")
    for key in ("dHfus", "dHvap"):
        if key in spec:
            v, dim, _ = to_si(spec[key], ["energy_per_amount", "energy_per_mass"], name=key)
            out[key] = dict(SRC, value=v * (MW if dim == "energy_per_mass" else 1.0), unit="J/mol")
    for key, dims, unit in (("density", "density", "kg/m3"), ("alpha", "per_temperature", "1/K"),
                            ("kappa_T", "per_pressure", "1/Pa")):
        if key in spec:
            q = spec[key]
            per_phase = q if isinstance(q, dict) and "value" not in q else {phase: q}
            out[key] = {}
            for ph, qq in per_phase.items():
                v = to_si(qq, dims, name=f"{name}.{key}")[0]
                entry = dict(SRC, value=v, unit=unit)
                if key == "density":
                    entry["type"] = "constant"
                out[key][ph] = entry
    return out


def extended_database(db: Database, custom: dict | None = None, pure_cp: dict | None = None,
                      log: Log | None = None) -> Database:
    """Copy of `db` with user components added and/or pure-component Cp overridden (db is unchanged)."""
    new = copy.copy(db)
    new.components = dict(db.components)
    new._alias = dict(db._alias)
    new._comp_json = dict(db._comp_json)
    for name, spec in (custom or {}).items():
        key = str(name)
        if key.lower() in new._alias:
            raise DataError(f"custom component '{name}' clashes with a database component; choose another "
                            "name, or use 'pure_cp' to override its Cp.")
        schema = custom_to_schema(key, spec)
        new._comp_json[key] = schema
        new.components[key] = Database._build(key, schema)
        for a in [key, schema["name"]] + schema["aliases"]:
            new._alias.setdefault(str(a).lower(), (key, None))
        if log:
            log.assume(f"Custom component '{key}' supplied by the user (miscibility class "
                       f"'{schema['miscibility_class']}'); its data are not checked against references.")
    for name, q in (pure_cp or {}).items():
        comp = copy.copy(new.get(name))
        comp.cp_corr = dict(comp.cp_corr)
        default_phase = "liquid" if comp.has_phase("liquid") else next(iter(comp.cp_corr))
        per_phase = q if isinstance(q, dict) and "value" not in q else {default_phase: q}
        for ph, qq in per_phase.items():
            v, dim, _ = to_si(qq, ["entropy_per_mass", "entropy_per_amount"], name=f"pure_cp[{name}]")
            molar = v * comp.MW if dim == "entropy_per_mass" else v
            comp.cp_corr[ph] = CpCorrelation("constant", [molar], 1.0, None, None, "user input (pure_cp)", False)
            if log:
                log.assume(f"Cp of {comp.name} ({ph}) set by the user to {molar:.4g} J/(mol K), constant in T.")
        new.components[comp.key] = comp
    return new
