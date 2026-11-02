"""Solver routes (A, B, C), limits, N-component systems and route self-consistency."""
import math

import pytest

from cp_engine.database import R


def cp(res, basis="mass"):
    assert res.success, str(res)
    return res.cp_mass if basis == "mass" else res.cp_molar


# ------------------------------------------------------------------ Route A
def test_pure_water(engine, db):
    res = engine.solve({"T": "25 C", "components": ["water"]})
    assert cp(res) == pytest.approx(db.get("water").cp_mass("liquid", 298.15))
    assert cp(res, "molar") == pytest.approx(db.get("water").cp("liquid", 298.15))


def test_mass_vs_mole_basis_inputs_agree(engine, db):
    w = {"water": 0.6, "ethanol": 0.4}
    MW = {k: db.get(k).MW for k in w}
    s = sum(w[k] / MW[k] for k in w)
    x = {k: w[k] / MW[k] / s for k in w}
    r1 = engine.solve({"T": 298.15, "mass_fractions": w})
    r2 = engine.solve({"T": 298.15, "mole_fractions": x})
    assert cp(r1) == pytest.approx(cp(r2), rel=1e-10)
    M = sum(x[k] * MW[k] for k in x)
    assert cp(r1, "molar") == pytest.approx(cp(r1) * M, rel=1e-10)


def test_ideal_solution_model_gives_sum_of_wi_cpi(engine, db):
    w = {"water": 0.5, "methanol": 0.3, "ethanol": 0.2}
    res = engine.solve({"T": 300, "mass_fractions": w, "solution_model": "ideal"})
    expected = sum(w[k] * db.get(k).cp_mass("liquid", 300) for k in w)
    assert cp(res) == pytest.approx(expected, rel=1e-10)
    assert res.breakdown.excess == 0.0


def test_real_solution_adds_excess(engine):
    res = engine.solve({"T": 300, "mole_fractions": {"water": 0.7, "ethanol": 0.3}})
    assert "real liquid solution" in res.classification
    assert res.breakdown.excess != 0.0


@pytest.mark.parametrize("x_eth", [1e-9, 1 - 1e-9])
def test_composition_limits_give_pure_values(engine, db, x_eth):
    res = engine.solve({"T": 300, "mole_fractions": {"water": 1 - x_eth, "ethanol": x_eth}})
    pure = "ethanol" if x_eth > 0.5 else "water"
    assert cp(res, "molar") == pytest.approx(db.get(pure).cp("liquid", 300), rel=1e-6)


@pytest.mark.parametrize("phi", [1e-9, 0.999999999])
def test_phi_limits(engine, db, phi):
    res = engine.solve({"T": 300, "components": ["water", "Al2O3"], "phi": {"Al2O3": phi}})
    pure = db.get("Al2O3").cp_mass("solid", 300) if phi > 0.5 else db.get("water").cp_mass("liquid", 300)
    assert cp(res) == pytest.approx(pure, rel=1e-6)


def test_nanofluid_route_A_equals_model_II(engine, db):
    from cp_engine.models.heterogeneous import nanofluid_cp
    res = engine.solve({"T": 300, "components": ["water", "Al2O3"], "phi": {"Al2O3": 0.03}})
    assert cp(res) == pytest.approx(nanofluid_cp(db, "Al2O3", "water", 0.03, 300, "II"), rel=1e-6)
    assert res.breakdown.extras["model_I_cp_mass"] == pytest.approx(
        nanofluid_cp(db, "Al2O3", "water", 0.03, 300, "I"), rel=1e-6)


@pytest.mark.parametrize("comps", [
    ["water", "ethanol"],
    ["water", "ethanol", "methanol"],
    ["water", "ethanol", "glycerol", "toluene", "Al2O3"],
])
def test_n_component_systems(engine, db, comps):
    w = {c: 1.0 / len(comps) for c in comps}
    res = engine.solve({"T": 300, "mass_fractions": w, "solution_model": "ideal"})
    phase = {c: ("solid" if c == "Al2O3" else "liquid") for c in comps}
    expected = sum(w[c] * db.get(c).cp_mass(phase[c], 300) for c in comps)
    assert cp(res) == pytest.approx(expected, rel=1e-9)


def test_five_component_emulsion_suspension_classification(engine):
    res = engine.solve({"T": 300, "mass_fractions": {"water": 0.4, "ethanol": 0.1, "glycerol": 0.1,
                                                     "toluene": 0.3, "Al2O3": 0.1}})
    assert "emulsion" in res.classification and "solid particles" in res.classification


# ------------------------------------------------------------------ Route B
def test_route_B_constant_P(engine):
    res = engine.solve({"T1": "20 C", "T2": "30 C", "total_mass": "2 kg", "process": "constant_P", "Q": "83.7 kJ"})
    assert cp(res) == pytest.approx(4185.0)
    assert "Route B" in res.recommended


def test_route_B_adiabatic_stirring(engine):
    res = engine.solve({"dT": "2.5 K", "m": "0.5 kg", "process": "adiabatic_constant_P", "W": "5.2 kJ"})
    assert cp(res) == pytest.approx(5200 / (0.5 * 2.5))


def test_route_B_constant_V_needs_cv_to_cp(engine):
    res = engine.solve({"T1": 298.15, "T2": 308.15, "components": ["water"], "total_mass": 1.0,
                        "process": "constant_V", "Q": "41.3 kJ"})
    labels = [r.label for r in res.routes]
    assert any("constant V" in l for l in labels)
    cv_route = [r for r in res.routes if "constant V" in r.label][0]
    assert cv_route.cp_mass > 4130.0           # Cp > Cv


def test_specific_heat_input_without_mass_uses_basis(engine):
    res = engine.solve({"T1": 300, "T2": 310, "process": "constant_P", "Q": "41.8 kJ/kg"})
    assert cp(res) == pytest.approx(4180.0)
    assert any("basis" in a for a in res.assumptions)


# ------------------------------------------------------------------ Route C
def test_route_C_from_dH_molar_only(engine):
    res = engine.solve({"dT": 10, "n": "2 mol", "dH": "1.5 kJ"})
    assert cp(res, "molar") == pytest.approx(75.0)
    assert res.cp_mass is None          # no composition → no molar mass


def test_route_C_from_dS(engine, db):
    w = db.get("water")
    T1, T2, m = 290.0, 320.0, 1.0
    from scipy.integrate import quad
    dS = m * quad(lambda T: w.cp_mass("liquid", T) / T, T1, T2)[0]
    res = engine.solve({"T1": T1, "T2": T2, "m": m, "dS": (dS, "J/K"), "components": ["water"]})
    r = [r for r in res.routes if r.tag == "dS"][0]
    assert r.cp_mass == pytest.approx(4180, rel=0.002)


def test_route_C_from_cv_matches_nist(engine):
    # NIST SRD 69, water, 300 K, 1 atm: cv = 4130.2 J/(kg K), cp = 4180.6 J/(kg K)
    res = engine.solve({"T": 300, "Cv": "4.1302 kJ/(kg K)", "components": ["water"]})
    r = [r for r in res.routes if r.tag == "Cv"][0]
    assert r.cp_mass == pytest.approx(4180.6, rel=0.002)


def test_route_C_from_gamma(engine):
    res = engine.solve({"Cv": "20.8 J/(mol K)", "gamma": 1.4, "components": ["n-hexane"]})
    assert cp(res, "molar") == pytest.approx(29.12)


def test_user_supplied_alpha_kappa(engine):
    res = engine.solve({"T": 300, "Cv": "1000 J/(kg K)", "alpha": "1e-3 1/K", "kappa_T": "1e-9 1/Pa",
                        "density": "800 kg/m3"})
    assert cp(res) == pytest.approx(1000 + 300 * (1 / 800) * 1e-6 / 1e-9)


# ------------------------------------------------------------------ self-consistency (validation case 5)
@pytest.mark.parametrize("system", [
    {"mass_fractions": {"water": 0.6, "ethanol": 0.4}},
    {"components": ["water", "CuO"], "phi": {"CuO": 0.05}},
    {"mass_fractions": {"water": 0.9, "ethanol": 0.1}},      # ice slurry when cooled
])
def test_route_A_to_B_roundtrip(engine, system):
    T1, T2 = 280.0, 290.0
    if "ethanol" in system.get("mass_fractions", {}) and system["mass_fractions"]["ethanol"] == 0.1:
        T1, T2 = 260.0, 268.0
    base = dict(system, T1=T1, T2=T2, total_mass=1.0)
    a = engine.solve(base)
    rA = [r for r in a.routes if r.tag == "A:mean"][0]
    Q = rA.breakdown.extras["dH_J"] * 1.0 / rA.breakdown.mass
    b = engine.solve({"T1": T1, "T2": T2, "total_mass": 1.0, "process": "constant_P", "Q": Q})
    assert cp(b) == pytest.approx(rA.cp_mass, rel=1e-7)   # floating-point round-off only


def test_overspecified_consistent_is_quiet_and_inconsistent_warns(engine):
    ok = engine.solve({"dT": 10, "m": 1, "process": "constant_P", "Q": 41800, "dH": 41800})
    assert not any("INCONSISTENT" in w for w in ok.warnings)
    bad = engine.solve({"dT": 10, "m": 1, "process": "constant_P", "Q": 41800, "dH": 50000})
    assert any("INCONSISTENT" in w for w in bad.warnings)


def test_routes_compared(engine):
    res = engine.solve({"T1": 293.15, "T2": 303.15, "components": ["water"], "m": 1,
                        "process": "constant_P", "Q": 41850.0})
    assert res.discrepancies and abs(res.discrepancies[0]["relative_difference_%"]) < 0.1


def test_derivation_lists_rules(engine):
    res = engine.solve({"T1": 300, "T2": 310, "m": 1, "process": "constant_P", "Q": 41800})
    text = res.derivation()
    assert "first_law_constant_P" in text and "cp_from_dH" in text and "delta_T" in text
