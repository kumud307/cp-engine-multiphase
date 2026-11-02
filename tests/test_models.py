"""Thermodynamic models: ideal, NRTL, heterogeneous, phase change, Cv-Cp."""
import math

import numpy as np
import pytest
import sympy as sp

from cp_engine.database import R
from cp_engine.models import cv_cp, heterogeneous as het, phase_change as pc
from cp_engine.models.ideal import cp_ideal_mass, cp_ideal_molar, mass_to_mole_fractions, mole_to_mass_fractions
from cp_engine.models.nrtl import NRTL, fit_nrtl_to_excess
from cp_engine.system import SystemModel


# ---------------------------------------------------------------- ideal mixing
def test_ideal_mass_and_molar_forms_agree():
    x = {"a": 0.3, "b": 0.7}
    MW = {"a": 0.018, "b": 0.046}
    Cp = {"a": 75.0, "b": 112.0}
    w = mole_to_mass_fractions(x, MW)
    cp_mass = cp_ideal_mass(w, {k: Cp[k] / MW[k] for k in Cp})
    M = sum(x[k] * MW[k] for k in x)
    assert cp_mass * M == pytest.approx(cp_ideal_molar(x, Cp))
    assert mass_to_mole_fractions(w, MW)["a"] == pytest.approx(0.3)


# ---------------------------------------------------------------- NRTL
@pytest.fixture
def binary():
    p = {"a": [[0, -0.8], [3.46, 0]], "b": [[0, 246.2], [-586.1, 0]], "c": np.full((2, 2), 0.3),
         "e": [[0, 0.5], [-0.3, 0]]}
    return NRTL(["1", "2"], params=p)


def test_zero_parameters_reduce_to_ideal():
    m = NRTL(["1", "2", "3"], params={"a": np.zeros((3, 3))})
    x = [0.2, 0.3, 0.5]
    assert m.gE_RT(x, 300) == 0.0 and m.cpE(x, 300) == pytest.approx(0.0, abs=1e-9)
    assert np.allclose(m.ln_gamma(x, 300), 0.0)


def test_excess_properties_vanish_for_pure_components(binary):
    for x in ([1.0, 0.0], [0.0, 1.0]):
        assert binary.gE_RT(x, 300) == pytest.approx(0.0, abs=1e-14)
        assert binary.hE(x, 300) == pytest.approx(0.0, abs=1e-6)
        assert binary.cpE(x, 300) == pytest.approx(0.0, abs=1e-5)


def test_gibbs_duhem_consistency(binary):
    x = [0.37, 0.63]
    assert binary.gE_RT(x, 310) == pytest.approx(sum(xi * li for xi, li in zip(x, binary.ln_gamma(x, 310))))


def test_numerical_derivatives_match_sympy(binary):
    """H^E and Cp^E from central differences vs exact symbolic derivatives."""
    T, x1 = sp.symbols("T x1", positive=True)
    x2 = 1 - x1
    a, b, e = binary.a, binary.b, binary.e
    t12 = a[0, 1] + b[0, 1] / T + e[0, 1] * sp.log(T)
    t21 = a[1, 0] + b[1, 0] / T + e[1, 0] * sp.log(T)
    G12, G21 = sp.exp(-0.3 * t12), sp.exp(-0.3 * t21)
    g = x1 * x2 * (t21 * G21 / (x1 + x2 * G21) + t12 * G12 / (x2 + x1 * G12))
    HE = -R * T ** 2 * sp.diff(g, T)
    CpE = sp.diff(HE, T)
    for xv, Tv in ((0.2, 290.0), (0.5, 310.0), (0.8, 340.0)):
        subs = {x1: xv, T: Tv}
        assert binary.gE_RT([xv, 1 - xv], Tv) == pytest.approx(float(g.subs(subs)), rel=1e-12)
        assert binary.hE([xv, 1 - xv], Tv) == pytest.approx(float(HE.subs(subs)), rel=1e-6)
        assert binary.cpE([xv, 1 - xv], Tv) == pytest.approx(float(CpE.subs(subs)), rel=1e-4, abs=1e-4)


def test_step_size_convergence(binary):
    table = binary.derivative_convergence([0.4, 0.6], 300.0)
    cp = {h: c for h, _, c in table}
    assert cp[0.1] == pytest.approx(cp[0.03], rel=1e-4)       # O(h²) plateau
    assert abs(cp[1e-5] - cp[0.03]) > abs(cp[0.1] - cp[0.03])  # round-off blows up at tiny h


def test_fit_nrtl_recovers_parameters():
    true = NRTL(["1", "2"], params={"a": [[0, 0.2], [0.5, 0]], "b": [[0, 150.0], [-300.0, 0]],
                                    "c": np.full((2, 2), 0.3)})
    # H^E at several temperatures is needed: at a single T the b-parameters are nearly degenerate
    xs = np.tile(np.linspace(0.05, 0.95, 10), 3)
    Ts = np.repeat([288.15, 298.15, 318.15], 10)
    HE = [true.hE([x, 1 - x], T) for x, T in zip(xs, Ts)]
    res = fit_nrtl_to_excess(xs, Ts, HE, kind="HE", fit=("b12", "b21"), initial={"a12": 0.2, "a21": 0.5})
    assert res["b_ij"] == pytest.approx(150.0, rel=1e-3) and res["b_ji"] == pytest.approx(-300.0, rel=1e-3)


# ---------------------------------------------------------------- heterogeneous
def test_model_II_equals_mass_weighting():
    phi, rs, cs, rf, cf = 0.04, 3970.0, 773.0, 997.0, 4180.0
    ws = het.volume_to_mass_fraction(phi, rs, rf)
    assert het.cp_model_II(phi, rs, cs, rf, cf) == pytest.approx(ws * cs + (1 - ws) * cf)
    assert het.cp_model_I(phi, cs, cf) > het.cp_model_II(phi, rs, cs, rf, cf)     # model I over-predicts


@pytest.mark.parametrize("phi", [0.0, 1.0])
def test_models_reduce_to_pure_limits(phi):
    for f in (het.cp_model_I(phi, 773.0, 4180.0), het.cp_model_II(phi, 3970.0, 773.0, 997.0, 4180.0)):
        assert f == pytest.approx(773.0 if phi == 1.0 else 4180.0)


def test_volume_mass_fraction_roundtrip():
    w = het.volume_to_mass_fraction(0.1, 6310.0, 997.0)
    assert het.mass_to_volume_fraction(w, 6310.0, 997.0) == pytest.approx(0.1)


# ---------------------------------------------------------------- phase change
def test_melt_fraction_properties():
    assert pc.melt_fraction(300.0, 300.0, 2.0) == pytest.approx(0.5)
    assert pc.melt_fraction(290.0, 300.0, 2.0) < 1e-12 and pc.melt_fraction(310.0, 300.0, 2.0) > 1 - 1e-12
    T = np.linspace(290, 310, 4001)
    area = np.trapezoid([pc.melt_fraction_derivative(t, 300.0, 2.0) for t in T], T)
    assert area == pytest.approx(1.0, rel=1e-6)


def test_ideal_solubility_and_freezing_point_are_inverse():
    x = pc.ideal_solubility(263.15, 273.15, 6009.0)
    assert pc.ideal_freezing_point(x, 273.15, 6009.0) == pytest.approx(263.15)
    assert pc.ideal_solubility(280.0, 273.15, 6009.0) == 1.0
    # with the ΔCp (Kirchhoff) term the inverse is solved numerically
    x2 = pc.ideal_solubility(250.0, 273.15, 6009.0, 37.6)
    assert pc.ideal_freezing_point(x2, 273.15, 6009.0, 37.6) == pytest.approx(250.0, abs=1e-8)
    assert pc.liquidus_temperature(x2, 273.15, 6009.0, None, 37.6) == pytest.approx(250.0, abs=1e-8)
    # an activity coefficient > 1 raises the freezing point (less depression)
    assert pc.liquidus_temperature(0.9, 273.15, 6009.0, lambda T: 0.05) > pc.liquidus_temperature(0.9, 273.15, 6009.0)


# ---------------------------------------------------------------- Cv / Cp
def test_ideal_gas_cp_minus_cv_equals_R():
    T, P = 400.0, 1e5
    v_molar = R * T / P
    assert cv_cp.cp_minus_cv(T, v_molar, 1 / T, 1 / P) == pytest.approx(R)


def test_cp_cv_roundtrip():
    cv = cv_cp.cv_from_cp(4180.0, 298.15, 1 / 997.0, 2.57e-4, 4.52e-10)
    assert cv_cp.cp_from_cv(cv, 298.15, 1 / 997.0, 2.57e-4, 4.52e-10) == pytest.approx(4180.0)
    # NIST: cp - cv of water at 298 K ≈ 0.034 kJ/(kg K)
    assert 4180.0 - cv == pytest.approx(43.6, rel=0.05)


# ---------------------------------------------------------------- system model
def test_latent_heat_conservation_through_melting_peak(db):
    """∫ cp_app dT over a melting peak = sensible heat + ΔH_fus exactly (energy conservation)."""
    s = SystemModel(db, {"n-octadecane": 1.0}, {}, options={"melt_width": 2.0})
    T1, T2 = 296.0, 306.0
    mean = s.cp_mean(T1, T2)
    T = np.linspace(T1, T2, 2001)
    pts = [s.cp_point(t).total for t in T]
    assert np.trapezoid(pts, T) == pytest.approx(mean.total * (T2 - T1), rel=2e-4)
    assert mean.latent * (T2 - T1) == pytest.approx(db.get("n-octadecane").dHfus, rel=1e-3)


def test_point_cp_equals_dH_dT(db):
    s = SystemModel(db, {"water": 0.9 / 0.018015, "ethanol": 0.1 / 0.046069}, {})
    T, h = 263.15, 0.01
    s.configure(T - 1, T + 1)
    num = (s.enthalpy(T + h) - s.enthalpy(T - h)) / (2 * h)
    assert s.cp_point(T).total == pytest.approx(num, rel=1e-5)


def test_ice_slurry_obeys_freezing_point_equation(db):
    """Default: ln(x_w γ_w) = SLE equation with the same NRTL used for H^E (HE+SLE parameter set)."""
    w = db.get("water")
    dCp = w.cp_corr["liquid"].raw(w.Tm) - w.cp_corr["solid"].raw(w.Tm)
    n = {"water": 0.9 / 0.018015, "ethanol": 0.1 / 0.046069}
    s = SystemModel(db, n, {})
    s.configure(263.15, 263.15)
    liq = {k: m for name, k, st, m in s.split(263.15) if st == "liquid"}
    x_w = liq["water"] / sum(liq.values())
    g_w = NRTL(["water", "ethanol"], db).gamma([x_w, 1 - x_w], 263.15)[0]
    assert x_w * g_w == pytest.approx(pc.ideal_solubility(263.15, w.Tm, w.dHfus, dCp), rel=1e-8)
    # forcing the ideal solution gives the textbook equation
    s2 = SystemModel(db, n, {}, options={"sle_activity": "ideal"})
    s2.configure(263.15, 263.15)
    liq = {k: m for name, k, st, m in s2.split(263.15) if st == "liquid"}
    assert liq["water"] / sum(liq.values()) == pytest.approx(pc.ideal_solubility(263.15, w.Tm, w.dHfus, dCp), rel=1e-10)


def test_no_ice_above_liquidus(db):
    s = SystemModel(db, {"water": 0.9 / 0.018015, "ethanol": 0.1 / 0.046069}, {})
    s.configure(275.0, 275.0)
    assert all(st != "solid" for _, _, st, n in s.split(275.0) if n > 0)
    assert s.cp_point(275.0).latent == pytest.approx(0.0, abs=1e-6)


def test_saturated_salt_slurry_dissolves_more_when_heated(db):
    s = SystemModel(db, {"water": 0.7 / 0.018015, "NaCl": 0.3 / 0.058443}, {})
    s.configure(290, 330)
    solid = lambda T: sum(n for name, k, st, n in s.split(T) if st == "solid")
    assert solid(330) < solid(290)
    assert s.cp_point(310).latent > 0       # endothermic dissolution (ΔH_sol > 0)
