"""Bad, missing and inconsistent input: clear messages, never tracebacks."""
import pytest

from cp_engine import InputError, PhaseError


def failed(res, text):
    assert not res.success
    blob = " ".join(res.errors + res.missing)
    assert text.lower() in blob.lower(), blob
    return res


def test_fractions_not_summing_to_one(engine):
    failed(engine.solve({"T": 300, "mass_fractions": {"water": 0.6, "ethanol": 0.6}}), "sum to")


def test_negative_fraction(engine):
    failed(engine.solve({"T": 300, "mole_fractions": {"water": 1.2, "ethanol": -0.2}}), "outside [0, 1]")


def test_one_missing_fraction_is_filled(engine):
    res = engine.solve({"T": 300, "components": ["water", "ethanol", "methanol"],
                        "mass_fractions": {"water": 0.5, "ethanol": 0.3}})
    assert res.success and any("methanol" in n and "0.2" in n for n in res.notes)


def test_two_missing_fractions_rejected(engine):
    failed(engine.solve({"T": 300, "components": ["water", "ethanol", "methanol"],
                         "mass_fractions": {"water": 0.5}}), "at most one")


def test_percent_fractions(engine):
    res = engine.solve({"T": 300, "mass_fractions": {"water": 60, "ethanol": 40}})
    assert res.success and any("percent" in a for a in res.assumptions)


def test_negative_amount(engine):
    failed(engine.solve({"T": 300, "masses": {"water": "-1 kg"}}), "negative")


def test_negative_temperature(engine):
    failed(engine.solve({"T": "-300 C", "components": ["water"]}), "absolute zero")


def test_unknown_input_key_suggests(engine):
    res = failed(engine.solve({"temprature": 300}), "Unknown input")
    assert "temperature" in res.errors[0]


def test_unknown_component(engine):
    failed(engine.solve({"T": 300, "mass_fractions": {"water": 0.5, "kerosene": 0.5}}), "Unknown component")


def test_phase_inconsistency_error(engine):
    failed(engine.solve({"T": "350 K", "components": {"water": "solid"}}), "melting point")


def test_liquid_below_melting_point_warns(engine):
    res = engine.solve({"T": "268 K", "components": {"water": "liquid"}})
    assert res.success and any("sub-cooled" in w for w in res.warnings)


def test_far_outside_correlation_range_is_reported(engine):
    res = engine.solve({"T": "600 K", "components": {"glycerol": "liquid"}})
    assert not res.success and any("refusing to extrapolate" in e for e in res.errors)


def test_slightly_outside_range_warns(engine):
    res = engine.solve({"T": "335 K", "components": {"n-dodecane": "liquid"}})
    assert res.success and any("outside the correlation range" in w for w in res.warnings)


def test_interval_crossing_phase_change(engine, db):
    res = engine.solve({"T1": "-5 C", "T2": "5 C", "components": ["water"], "total_mass": 1})
    r = [r for r in res.routes if r.tag == "A:mean"][0]
    latent = r.breakdown.per_mass()["latent"] * 10.0         # J/kg absorbed as latent heat
    assert latent == pytest.approx(db.get("water").dHfus / db.get("water").MW, rel=0.01)
    assert "phase change" in res.classification


def test_underspecified_reports_missing(engine):
    res = failed(engine.solve({"T": 300}), "composition")
    assert any("Route A" in m for m in res.missing)


def test_underspecified_calorimetry(engine):
    res = failed(engine.solve({"Q": "10 kJ", "process": "constant_P", "dT": 5}), "total mass")


def test_adiabatic_with_heat_is_rejected(engine):
    failed(engine.solve({"process": "adiabatic", "Q": 100, "dT": 1, "m": 1}), "adiabatic")


def test_constant_P_and_V_rejected(engine):
    failed(engine.solve({"process": ["constant_P", "constant_V"], "T": 300}), "both")


def test_inconsistent_temperatures(engine):
    failed(engine.solve({"T1": 300, "T2": 310, "dT": 5, "m": 1, "Q": 100}), "Inconsistent temperatures")


def test_composition_overspecified(engine):
    failed(engine.solve({"T": 300, "mass_fractions": {"water": 1}, "mole_fractions": {"water": 1}}),
           "over-specified")


def test_total_mass_conflicts_with_masses(engine):
    failed(engine.solve({"T": 300, "masses": {"water": "1 kg"}, "total_mass": "2 kg"}), "conflicts")


def test_wrong_unit_dimension(engine):
    failed(engine.solve({"T": 300, "components": ["water"], "Q": "5 kg"}), "expected energy")


def test_gamma_below_one(engine):
    failed(engine.solve({"Cv": 1000, "gamma": 0.9}), "impossible")


def test_raise_errors_flag(engine):
    with pytest.raises(PhaseError):
        engine.solve({"T": 350, "components": {"water": "solid"}}, raise_errors=True)
    with pytest.raises(InputError):
        engine.solve({"T": 300, "mass_fractions": {"water": 2}}, raise_errors=True)


def test_ice_alias_and_merge(engine):
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
    from iapws_ice import props          # IAPWS R10-06 ice equation of state
    res = engine.solve({"T": "-10 C", "components": ["ice"]})
    assert res.success and res.cp_mass == pytest.approx(props(263.15, 101325.0)["cp"], rel=1e-4)
    res2 = engine.solve({"T": "0 C", "masses": {"ice": "0.2 kg", "water": "0.8 kg"}})
    assert res2.success and any("more than one phase" in n for n in res2.notes)
