"""Equations of state, vapour–liquid equilibrium, real gases, partial molar properties, uncertainty.

The EOS tests reproduce the worked CO₂ and CO₂ + CH₄ examples of the CL 202 lectures
("Equations of State for Pure Fluids", slides 30–33; "Equations of State for Mixtures", slides 23–24).
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from cp_engine import CpEngine
from cp_engine.database import load_database
from cp_engine.models import partial_molar as pm
from cp_engine.models.eos import CubicEOS, Critical, VirialEOS
from cp_engine.models.nrtl import NRTL
from cp_engine.models.vle import VLE

DATA = __import__("pathlib").Path(__file__).resolve().parent.parent / "validation" / "data"
CO2 = Critical("CO2", 304.13, 73.77e5, 0.224, 94.07e-6)      # lecture slide 30
CH4 = Critical("CH4", 190.564, 45.992e5, 0.0114, 98.6e-6)


@pytest.fixture(scope="module")
def engine():
    return CpEngine()


@pytest.fixture(scope="module")
def db():
    return load_database()


# ---------------------------------------------------------------- lecture examples
def test_abbott_B_of_co2_matches_lecture():
    assert VirialEOS([CO2]).B_mix(350.0, [1.0])[0] * 1e6 == pytest.approx(-83.73, abs=0.01)


@pytest.mark.parametrize("kind,psat,vliq", [("vdW", 52.55, 81.13), ("RK", 45.62, 61.03),
                                            ("SRK", 41.99, 58.42), ("PR", 41.60, 51.68)])
def test_cubic_vapour_pressure_and_liquid_volume_of_co2_match_lecture(kind, psat, vliq):
    e = CubicEOS([CO2], kind)
    P = e.psat(280.0)
    assert P / 1e5 == pytest.approx(psat, abs=0.02)
    assert e.molar_volume(280.0, P, [1.0], "liquid") * 1e6 == pytest.approx(vliq, abs=0.02)


@pytest.mark.parametrize("P,virial,pr0,pr10", [(10, 0.9825, 0.9782, 0.9798), (20, 0.9650, 0.9569, 0.9600),
                                               (40, 0.9300, 0.9157, 0.9223), (60, 0.8951, 0.8771, 0.8874)])
def test_co2_methane_compressibility_matches_lecture(P, virial, pr0, pr10):
    y = [0.5, 0.5]
    assert VirialEOS([CO2, CH4]).Z(350.0, P * 1e5, y) == pytest.approx(virial, abs=1e-4)
    assert CubicEOS([CO2, CH4], "PR").Z(350.0, P * 1e5, y) == pytest.approx(pr0, abs=1e-4)
    assert CubicEOS([CO2, CH4], "PR", {(0, 1): 0.10}).Z(350.0, P * 1e5, y) == pytest.approx(pr10, abs=1e-4)


def test_cross_virial_coefficient_and_fugacity_coefficient_match_lecture():
    v = VirialEOS([CO2, CH4])
    assert v.B_matrix(350.0)[0][0, 1] * 1e6 == pytest.approx(-47.0, abs=0.05)
    assert math.exp(v.ln_phi(350.0, 20e5, [0.5, 0.5])[0]) == pytest.approx(0.947, abs=0.0005)


def test_phi_phi_bubble_pressure_reproduces_lecture_co2_methane_vle():
    """Slide 26: at 273.15 K and 6.375 MPa the measured liquid has 12.5 % CH4; PR with k12 = 0 predicts 18.6 %,
    with k12 = 0.10 the data are reproduced. Pure CO2 (x_CH4 → 0) gives its vapour pressure, 3.48 MPa."""
    from scipy.optimize import brentq
    x_at = {}
    for k in (0.0, 0.10):
        e = CubicEOS([CO2, CH4], "PR", {(0, 1): k})
        x_at[k] = brentq(lambda xm: e.bubble_P(273.15, [1 - xm, xm])[0] - 6.375e6, 0.05, 0.2)
    assert x_at[0.0] == pytest.approx(0.186, abs=0.002)
    assert x_at[0.10] == pytest.approx(0.125, abs=0.006)
    P, y = CubicEOS([CO2, CH4], "PR", {(0, 1): 0.10}).bubble_P(273.15, [1 - 1e-5, 1e-5])
    assert P / 1e6 == pytest.approx(3.48, abs=0.01)


def test_lewis_randall_limit_when_delta12_is_zero():
    """δ12 = 2B12 − B11 − B22 = 0 (identical components) gives φ̂_i = φ_i (Lewis–Randall rule)."""
    v = VirialEOS([CO2, CO2])
    ln_mix = v.ln_phi(350.0, 20e5, [0.3, 0.7])
    ln_pure = VirialEOS([CO2]).ln_phi(350.0, 20e5, [1.0])
    assert np.allclose(ln_mix, ln_pure[0], atol=1e-12)


@pytest.mark.parametrize("kind", ["PR", "SRK", "virial"])
def test_residual_properties_vanish_in_the_ideal_gas_limit(kind, db):
    e = db.eos(["water"], kind)
    assert abs(e.H_residual(600.0, 1.0, [1.0])) < 1e-3
    assert abs(e.Cp_residual(600.0, 1.0, [1.0])) < 1e-5


def test_pr_reproduces_the_definition_of_omega(db):
    """Soave/PR α(T) was fitted to vapour pressures: log10(P_sat/P_c) at T_r = 0.7 ≈ −1 − ω."""
    for k in ("water", "methanol", "n-hexane", "toluene"):
        c = db.get(k).critical
        p = db.eos([k], "PR").psat(0.7 * c.Tc)
        assert -1.0 - math.log10(p / c.Pc) == pytest.approx(c.omega, abs=0.03)


# ---------------------------------------------------------------- data
def test_critical_constants_and_omega_from_data(db):
    w = db.get("water").critical
    assert w.Tc == pytest.approx(647.096, abs=0.01) and w.Pc == pytest.approx(22.064e6, rel=1e-4)
    for k, om in (("water", 0.344), ("methanol", 0.565), ("n-hexane", 0.300), ("toluene", 0.266), ("ethanol", 0.643)):
        assert db.get(k).critical.omega == pytest.approx(om, abs=0.003)


@pytest.mark.parametrize("k", ["water", "methanol", "n-hexane", "toluene"])
def test_antoine_fits_reproduce_nist_saturation_curve(k, db):
    d = pd.read_csv(DATA / f"nist_sat_{k}.csv", comment="#").dropna(subset=["P_MPa"])
    c = db.get(k)
    lo = 273.2 if k == "water" else 300.0
    d = d[(d.T_K >= lo) & (d.T_K <= 450.0)]
    dev = [abs(c.psat(T) / (P * 1e6) - 1) for T, P in zip(d.T_K, d.P_MPa)]
    assert max(dev) < 0.01


# ---------------------------------------------------------------- VLE
def test_bubble_points_of_water_ethanol_against_measurements(db):
    d = pd.read_csv(DATA / "thermoml_vle_water_alcohols.csv", comment="#")
    d = d[(d.alcohol == "ethanol") & (d.doi == "10.1021/je2008704") & (d.x_alcohol > 0) & (d.x_alcohol < 1)]
    v = VLE(["ethanol", "water"], db, NRTL(["ethanol", "water"], db, prefer="VLE"))
    dT = [abs(v.bubble_T(P * 1e3, [x, 1 - x]) - T) for x, P, T in zip(d.x_alcohol, d.P_kPa, d.T_K)]
    assert np.mean(dT) < 0.6 and max(dT) < 1.5


def test_ethanol_water_azeotrope(db):
    v = VLE(["ethanol", "water"], db, NRTL(["ethanol", "water"], db, prefer="VLE"))
    x = [0.9, 0.1]
    Tb, Td = v.bubble_T(101325.0, x), v.dew_T(101325.0, x)
    assert abs(Tb - Td) < 0.01 and Tb == pytest.approx(351.3, abs=0.3)   # minimum-boiling azeotrope


def test_flash_material_balance_and_limits(db):
    v = VLE(["ethanol", "water"], db, NRTL(["ethanol", "water"], db, prefer="VLE"))
    z = np.array([0.3, 0.7])
    V, x, y = v.flash(355.0, 101325.0, z)
    assert 0 < V < 1
    assert np.allclose((1 - V) * x + V * y, z, atol=1e-10)
    K = v.K(355.0, 101325.0, x)
    assert np.allclose(y, K * x, rtol=1e-8)                              # equal fugacities
    assert v.flash(340.0, 101325.0, z)[0] == 0.0
    assert v.flash(380.0, 101325.0, z)[0] == 1.0


def test_non_volatile_component_has_no_dew_point(db):
    v = VLE(["water", "glycerol"], db)
    assert math.isinf(v.dew_T(101325.0, [0.8, 0.2]))
    assert v.bubble_T(101325.0, [0.8, 0.2]) > 373.15                    # boiling-point elevation


# ---------------------------------------------------------------- engine
def test_mixture_boils_as_a_mixture(engine):
    # 5 mol % ethanol at 85 °C: below the bubble point (≈ 90 °C) → liquid, although T > T_b(ethanol)
    r = engine.solve({"T": "85 C", "mole_fractions": {"ethanol": 0.05, "water": 0.95}})
    assert r.classification.startswith("Homogeneous liquid")
    # 50 mol % at 90 °C: above the dew point → vapour, not "liquid water + gaseous ethanol"
    r = engine.solve({"T": "90 C", "mole_fractions": {"ethanol": 0.5, "water": 0.5}})
    assert r.classification.startswith("Gas mixture")
    # inside the two-phase region: vapour–liquid equilibrium with latent heat
    r = engine.solve({"T": "84 C", "mole_fractions": {"ethanol": 0.3, "water": 0.7}})
    assert "vapour–liquid equilibrium" in r.classification and r.breakdown.latent > 0
    assert 0 < r.breakdown.extras["vapour_mole_fraction_of_system"] < 1


def test_energy_conservation_through_the_boiling_range(engine):
    from scipy import integrate
    s = engine.system({"mole_fractions": {"ethanol": 0.3, "water": 0.7}})
    s.configure(345.0, 380.0)
    info = next(iter(s.vle_info.values()))
    integral, _ = integrate.quad(lambda t: s.cp_point(float(t)).total, 345.0, 380.0, limit=400,
                                 points=[info["T_bubble"], info["T_dew"]])
    assert integral == pytest.approx(s.enthalpy(380.0) - s.enthalpy(345.0), rel=1e-3)


def test_superheated_steam_at_10_MPa(engine):
    d = pd.read_csv(DATA / "nist_water_vapour_10MPa.csv", comment="#")
    ref = float(d[d.T_K == 770.0].cp_J_kgK.iloc[0])
    pr = engine.solve({"T": 770.0, "P": (10, "MPa"), "components": {"water": "gas"}}).cp_mass
    ideal = engine.solve({"T": 770.0, "P": (10, "MPa"), "components": {"water": "gas"}, "gas_eos": "ideal"}).cp_mass
    assert abs(pr / ref - 1) < 0.05 < abs(ideal / ref - 1)


def test_unused_heat_input_is_reported(engine):
    r = engine.solve({"T1": "20 C", "T2": "30 C", "Q": (20.0, "kJ"), "total_mass": "1 kg", "components": ["water"]})
    assert any("NOT used" in w and "constant pressure" in w for w in r.warnings)


def test_uncertainty_estimate_is_validation_based(engine):
    r = engine.solve({"T": "25 C", "mole_fractions": {"methanol": 0.3, "water": 0.7}})
    u = r.uncertainty
    assert 1.0 < u["relative_percent"] < 5.0
    assert any("methanol+water" in c["source"] for c in u["contributions"])
    assert all(c["validation"] for c in u["contributions"])
    assert engine.solve({"T": "25 C", "components": ["water"]}).uncertainty["relative_percent"] < 0.5


def test_uncertainty_table_matches_validation_results():
    """data/model_uncertainty.json must be regenerated whenever the validation results change."""
    import json
    root = DATA.parent.parent
    tab = json.loads((root / "data" / "model_uncertainty.json").read_text(encoding="utf-8"))
    d = pd.read_csv(root / "validation" / "results" / "summary.csv")
    for entry in (tab["mixtures"]["methanol+water"], tab["suspension"], tab["ideal_fallback"]):
        case, rest = entry["row"].split(": ", 1)
        system, quantity = rest.split(": ", 1)
        r = d[(d.case.astype(str) == case.replace("case ", "")) & (d.system == system) & (d.quantity == quantity)]
        assert len(r) == 1 and float(r["mean"].iloc[0]) == pytest.approx(entry["value"], abs=1e-3)


# ---------------------------------------------------------------- partial molar properties
@pytest.fixture(scope="module")
def nrtl(db):
    return NRTL(["ethanol", "water"], db)


def test_ln_gamma_is_the_partial_molar_excess_gibbs_energy(nrtl):
    x = [0.3, 0.7]
    assert np.allclose(pm.ln_gamma_from_gE(nrtl, x, 320.0), nrtl.ln_gamma(x, 320.0), atol=1e-7)


def test_euler_relation_and_lecture_formula(nrtl):
    x, T = [0.3, 0.7], 320.0
    Hbar = pm.partial_molar_hE(nrtl, x, T)
    assert np.dot(x, Hbar) == pytest.approx(nrtl.hE(x, T), rel=1e-4)                # H^E = Σ x_i H̄_i^E
    for k in (0, 1):                                                                 # B̄_k = B + ∂B/∂x_k − Σ x_j ∂B/∂x_j
        assert pm.partial_molar_from_x(lambda z: nrtl.hE(list(z), T), x, k) == pytest.approx(Hbar[k], rel=1e-4)
    assert np.dot(x, pm.partial_molar_cpE(nrtl, x, T)) == pytest.approx(nrtl.cpE(x, T), rel=1e-3)


@pytest.mark.parametrize("x", [[0.1, 0.9], [0.5, 0.5], [0.8, 0.2]])
def test_gibbs_duhem(nrtl, x):
    res = pm.gibbs_duhem_residual(lambda z: np.asarray(nrtl.ln_gamma(list(z), 320.0)), x, 0)
    assert abs(res) < 1e-8


def test_mixing_functions(nrtl):
    f = pm.mixing_functions(nrtl, [0.3, 0.7], 320.0)
    assert f["dH_mix"] == pytest.approx(nrtl.hE([0.3, 0.7], 320.0))
    assert f["dS_mix_ideal"] == pytest.approx(-8.314462618 * (0.3 * math.log(0.3) + 0.7 * math.log(0.7)))
    assert f["dG_mix"] == pytest.approx(f["dH_mix"] - 320.0 * f["dS_mix"])


def test_engine_reports_partial_molar_cp(engine):
    r = engine.solve({"T": "25 C", "mole_fractions": {"ethanol": 0.3, "water": 0.7}})
    pc = r.breakdown.extras["partial_molar_cp_J_per_molK"]["liquid"]
    assert 0.3 * pc["Ethanol"] + 0.7 * pc["Water"] == pytest.approx(r.cp_molar, rel=1e-3)
