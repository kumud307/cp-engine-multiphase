"""Unit parsing and conversion to SI.

Every number that enters the engine passes through this module, so all
internal calculations use SI units only:

    temperature  K          pressure  Pa         energy  J
    mass         kg         amount    mol        volume  m^3

A quantity can be given as
    * a plain number            -> assumed to already be SI (a note is recorded)
    * a tuple   (value, "unit") -> e.g. (25, "C"), (4.18, "kJ/(kg K)")
    * a string  "value unit"    -> e.g. "25 °C", "1 atm", "2.5 kJ"
    * a dict    {"value": v, "unit": "u"}

Compound units are parsed from their pieces, e.g. "kJ/(kg K)", "J/kg/K",
"J kg^-1 K^-1", "cal/(g*°C)" are all understood.

Temperature needs care: an *absolute* temperature in °C or °F needs an
offset (25 °C = 298.15 K) but a temperature *difference* or a unit such as
J/(kg °C) does not (1 °C of difference = 1 K of difference).
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

from .errors import UnitError

# ---------------------------------------------------------------------------
# Base unit table: symbol -> (factor to SI, {base dimension: exponent})
# Base dimensions: E energy, M mass, N amount, T temperature, P pressure, L length
# ---------------------------------------------------------------------------
_E, _M, _N, _T, _P, _L = "E", "M", "N", "T", "P", "L"

_UNITS: dict[str, tuple[float, dict[str, int]]] = {
    # energy
    "J": (1.0, {_E: 1}), "kJ": (1e3, {_E: 1}), "MJ": (1e6, {_E: 1}),
    "cal": (4.184, {_E: 1}), "kcal": (4184.0, {_E: 1}),          # thermochemical calorie
    "Wh": (3600.0, {_E: 1}), "kWh": (3.6e6, {_E: 1}),
    "BTU": (1055.05585, {_E: 1}), "Btu": (1055.05585, {_E: 1}),  # International Table BTU
    # mass
    "kg": (1.0, {_M: 1}), "g": (1e-3, {_M: 1}), "mg": (1e-6, {_M: 1}),
    "t": (1e3, {_M: 1}), "tonne": (1e3, {_M: 1}), "lb": (0.45359237, {_M: 1}),
    # amount of substance
    "mol": (1.0, {_N: 1}), "kmol": (1e3, {_N: 1}), "mmol": (1e-3, {_N: 1}),
    # temperature (as a *difference*; absolute offsets are handled separately)
    "K": (1.0, {_T: 1}), "C": (1.0, {_T: 1}), "degC": (1.0, {_T: 1}), "°C": (1.0, {_T: 1}),
    "F": (5.0 / 9.0, {_T: 1}), "degF": (5.0 / 9.0, {_T: 1}), "°F": (5.0 / 9.0, {_T: 1}),
    # pressure
    "Pa": (1.0, {_P: 1}), "kPa": (1e3, {_P: 1}), "MPa": (1e6, {_P: 1}), "GPa": (1e9, {_P: 1}),
    "bar": (1e5, {_P: 1}), "mbar": (1e2, {_P: 1}), "atm": (101325.0, {_P: 1}),
    "psi": (6894.757293168, {_P: 1}), "mmHg": (133.322387415, {_P: 1}), "Torr": (101325.0 / 760.0, {_P: 1}),
    # length / volume
    "m": (1.0, {_L: 1}), "cm": (1e-2, {_L: 1}), "mm": (1e-3, {_L: 1}),
    "L": (1e-3, {_L: 3}), "l": (1e-3, {_L: 3}), "mL": (1e-6, {_L: 3}), "ml": (1e-6, {_L: 3}),
    # dimensionless
    "1": (1.0, {}), "-": (1.0, {}), "%": (0.01, {}),
}

# Named dimensions used by the input layer (exponent dictionaries).
DIMENSIONS: dict[str, dict[str, int]] = {
    "temperature": {_T: 1},
    "pressure": {_P: 1},
    "energy": {_E: 1},
    "energy_per_mass": {_E: 1, _M: -1},
    "energy_per_amount": {_E: 1, _N: -1},
    "entropy": {_E: 1, _T: -1},                       # J/K  (also total heat capacity)
    "entropy_per_mass": {_E: 1, _M: -1, _T: -1},      # J/(kg K)
    "entropy_per_amount": {_E: 1, _N: -1, _T: -1},    # J/(mol K)
    "mass": {_M: 1},
    "amount": {_N: 1},
    "volume": {_L: 3},
    "density": {_M: 1, _L: -3},
    "per_temperature": {_T: -1},
    "per_pressure": {_P: -1},
    "molar_mass": {_M: 1, _N: -1},
    "dimensionless": {},
}

_SI_UNIT_NAME = {
    "temperature": "K", "pressure": "Pa", "energy": "J", "energy_per_mass": "J/kg",
    "energy_per_amount": "J/mol", "entropy": "J/K", "entropy_per_mass": "J/(kg K)",
    "entropy_per_amount": "J/(mol K)", "mass": "kg", "amount": "mol", "volume": "m^3",
    "density": "kg/m^3", "per_temperature": "1/K", "per_pressure": "1/Pa",
    "molar_mass": "kg/mol", "dimensionless": "-",
}


@dataclass
class ParsedUnit:
    """Result of parsing a unit string."""
    factor: float                 # multiply the value by this to get SI
    dims: dict[str, int]          # exponents of the base dimensions
    temp_scale: str | None        # "C", "F" or "K" if the unit is a bare temperature
    text: str                     # original text, for messages


def _clean(unit: str) -> str:
    u = unit.strip()
    for a, b in (("·", " "), ("*", " "), ("⋅", " "), ("−", "-"), ("º", "°"),
                 ("deg C", "degC"), ("deg F", "degF"), ("° C", "°C"), ("° F", "°F"),
                 ("³", "3"), ("²", "2"), ("⁻¹", "^-1")):
        u = u.replace(a, b)
    u = re.sub(r"\bkelvin\b", "K", u, flags=re.IGNORECASE)
    u = re.sub(r"\bcelsius\b", "C", u, flags=re.IGNORECASE)
    return u


def _parse_token(tok: str) -> tuple[float, dict[str, int]]:
    """Parse one unit token such as 'kJ', 'm3', 'K^-1' or 'cm^3'."""
    m = re.fullmatch(r"([A-Za-z°%\-1]+?)\^?(-?\d+)?", tok)
    if not m:
        raise UnitError(f"Cannot understand unit piece '{tok}'.")
    sym, exp = m.group(1), int(m.group(2)) if m.group(2) else 1
    if sym not in _UNITS:
        # 'm3' style may have been split oddly; give a helpful message.
        raise UnitError(f"Unknown unit '{sym}'. Known units: {', '.join(sorted(_UNITS))}.")
    f, d = _UNITS[sym]
    return f ** exp, {k: v * exp for k, v in d.items()}


def _add_dims(a: dict[str, int], b: dict[str, int], sign: int = 1) -> dict[str, int]:
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + sign * v
        if out[k] == 0:
            del out[k]
    return out


def parse_unit(unit: str) -> ParsedUnit:
    """Parse a (possibly compound) unit string into an SI factor and dimensions.

    Grammar (informal): numerator ['/' denominator]*, where each part is a
    space-separated product of tokens; parentheses group a denominator,
    e.g. 'kJ/(kg K)'.
    """
    text = unit
    u = _clean(unit)
    if u == "":
        return ParsedUnit(1.0, {}, None, text)
    # split "a/(b c)/d" into parts at top-level slashes
    parts, depth, cur = [], 0, ""
    for ch in u:
        if ch == "(":
            depth += 1
            continue
        if ch == ")":
            depth -= 1
            continue
        if ch == "/" and depth == 0:
            parts.append(cur)
            cur = ""
            continue
        cur += ch
    parts.append(cur)
    factor, dims = 1.0, {}
    for i, part in enumerate(parts):
        sign = 1 if i == 0 else -1
        for tok in part.split():
            f, d = _parse_token(tok)
            factor *= f ** sign
            dims = _add_dims(dims, d, sign)
    scale = None
    stripped = u.replace(" ", "")
    if dims == {_T: 1} and len(parts) == 1:
        scale = {"K": "K", "C": "C", "degC": "C", "°C": "C",
                 "F": "F", "degF": "F", "°F": "F"}.get(stripped)
    return ParsedUnit(factor, dims, scale, text)


def dimension_name(dims: dict[str, int]) -> str | None:
    """Return the name of a dimension dictionary (or None if unknown)."""
    for name, d in DIMENSIONS.items():
        if d == dims:
            return name
    return None


def split_quantity(q: Any) -> tuple[float, str | None]:
    """Split the accepted input forms into (value, unit-or-None)."""
    if isinstance(q, bool):
        raise UnitError(f"Expected a number, got {q!r}.")
    if isinstance(q, (int, float)):
        return float(q), None
    if isinstance(q, dict) and "value" in q:
        return float(q["value"]), q.get("unit") or q.get("units")
    if isinstance(q, (tuple, list)) and len(q) == 2:
        try:
            return float(q[0]), str(q[1])
        except (TypeError, ValueError) as exc:
            raise UnitError(f"Cannot read quantity {q!r}: {exc}") from None
    if isinstance(q, str):
        m = re.fullmatch(r"\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*(.*?)\s*", q)
        if not m:
            raise UnitError(f"Cannot read quantity '{q}'. Use e.g. '25 C' or (25, 'C').")
        return float(m.group(1)), (m.group(2) or None)
    raise UnitError(f"Cannot read quantity {q!r}. Use a number, (value, 'unit') or 'value unit'.")


def to_si(q: Any, allowed: list[str] | str, absolute_temperature: bool = False,
          name: str = "quantity") -> tuple[float, str, bool]:
    """Convert a quantity to SI.

    Parameters
    ----------
    q : the user's quantity (see module docstring for accepted forms)
    allowed : name(s) of acceptable dimensions (see DIMENSIONS)
    absolute_temperature : apply °C/°F offsets (True for T, T1, T2; False for ΔT)
    name : input name used in error messages

    Returns
    -------
    (value_SI, dimension_name, unit_was_assumed)
    """
    allowed = [allowed] if isinstance(allowed, str) else list(allowed)
    value, unit = split_quantity(q)
    if unit is None:
        # No unit: assume SI in the first allowed dimension.
        return value, allowed[0], True
    pu = parse_unit(unit)
    dname = dimension_name(pu.dims)
    if dname not in allowed:
        expect = " or ".join(f"{a} [{_SI_UNIT_NAME[a]}]" for a in allowed)
        got = f"{dname} [{_SI_UNIT_NAME[dname]}]" if dname else "an unrecognised dimension"
        raise UnitError(f"'{name}' was given in '{unit}', which is {got}; expected {expect}.")
    if dname == "temperature" and absolute_temperature:
        scale = pu.temp_scale or "K"
        if scale == "C":
            return value + 273.15, dname, False
        if scale == "F":
            return (value - 32.0) * 5.0 / 9.0 + 273.15, dname, False
        return value, dname, False
    return value * pu.factor, dname, False


def convert(value: float, from_unit: str, to_unit: str) -> float:
    """Convert a value between two units of the same dimension (no offsets).

    Example: convert(4.18, "kJ/(kg K)", "J/(g K)") -> 4.18
    """
    a, b = parse_unit(from_unit), parse_unit(to_unit)
    if a.dims != b.dims:
        raise UnitError(f"Cannot convert '{from_unit}' to '{to_unit}': different dimensions.")
    return value * a.factor / b.factor


def si_unit(dimension: str) -> str:
    """SI unit symbol used for reporting a dimension."""
    return _SI_UNIT_NAME[dimension]
