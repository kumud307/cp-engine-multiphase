"""Unit parsing and conversion."""
import pytest

from cp_engine.errors import UnitError
from cp_engine.units import convert, to_si


@pytest.mark.parametrize("q, expected", [
    ("25 C", 298.15), ("25 °C", 298.15), ((25, "degC"), 298.15), ("77 F", 298.15),
    ((298.15, "K"), 298.15), ({"value": 0, "unit": "C"}, 273.15),
])
def test_absolute_temperatures(q, expected):
    v, _, _ = to_si(q, "temperature", absolute_temperature=True)
    assert v == pytest.approx(expected, abs=1e-9)


def test_temperature_difference_has_no_offset():
    assert to_si("10 C", "temperature", absolute_temperature=False)[0] == pytest.approx(10.0)
    assert to_si("18 F", "temperature", absolute_temperature=False)[0] == pytest.approx(10.0)


@pytest.mark.parametrize("q, dims, expected", [
    ("1 atm", "pressure", 101325.0), ("2 bar", "pressure", 2e5), ("760 Torr", "pressure", 101325.0), ("1 mmHg", "pressure", 133.322387415),
    ("1 kJ", "energy", 1000.0), ("1 kcal", "energy", 4184.0), ("1 cal", "energy", 4.184),
    ("500 g", "mass", 0.5), ("2 kmol", "amount", 2000.0), ("1 L", "volume", 1e-3),
    ("4.18 kJ/(kg K)", "entropy_per_mass", 4180.0), ("4.18 J/g/K", "entropy_per_mass", 4180.0),
    ("1 cal/(g*°C)", "entropy_per_mass", 4184.0), ("75.3 J mol^-1 K^-1", "entropy_per_amount", 75.3),
    ("2.5e-4 1/K", "per_temperature", 2.5e-4), ("4.5e-10 1/Pa", "per_pressure", 4.5e-10),
    ("0.997 g/cm3", "density", 997.0),
])
def test_conversions(q, dims, expected):
    assert to_si(q, dims)[0] == pytest.approx(expected, rel=1e-12)


def test_unitless_number_is_assumed_si():
    v, dim, assumed = to_si(300.0, "temperature", True)
    assert (v, dim, assumed) == (300.0, "temperature", True)


def test_wrong_dimension_is_rejected():
    with pytest.raises(UnitError, match="expected"):
        to_si("5 kg", "energy", name="Q")


def test_unknown_unit_is_rejected():
    with pytest.raises(UnitError, match="Unknown unit"):
        to_si("5 furlongs", "energy")


def test_convert_helper():
    assert convert(4.184, "kJ/(kg K)", "J/(g K)") == pytest.approx(4.184)
    with pytest.raises(UnitError):
        convert(1, "J", "kg")
