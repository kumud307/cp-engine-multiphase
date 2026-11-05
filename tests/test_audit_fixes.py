"""Regression tests for the corrections made after the independent audit (MODEL_VERIFICATION_REPORT.md)."""
from __future__ import annotations

import math

import pytest

from cp_engine import CpEngine
from cp_engine.models.phase_change import TAIL_SIGMAS, melt_fraction, melt_fraction_derivative


@pytest.fixture(scope="module")
def engine():
    return CpEngine()


# ---- G2: non-finite inputs -----------------------------------------------------------------------
@pytest.mark.parametrize("inp", [
    {"T": float("nan"), "components": ["water"]},
    {"T": 300.0, "P": float("inf"), "components": ["water"]},
    {"T1": float("nan"), "T2": 300.0, "components": ["water"]},
    {"T": 300.0, "mass_fractions": {"water": float("nan"), "ethanol": 0.5}},
])
def test_non_finite_inputs_are_rejected(engine, inp):
    r = engine.solve(inp)
    assert not r.success and "finite" in r.errors[0]


# ---- G1: saturated vapour in the VLE region keeps its departure term, no false warning -------------
def test_saturated_vapour_is_not_flagged(engine):
    s = engine.system({"mole_fractions": {"ethanol": 0.1, "water": 0.9}})
    s.configure(355.0, 372.0)
    for T in (360.0, 365.0, 370.0):
        assert s._gas_residual(T, s.split(T))[0] > 0
    assert not any("below its dew point" in w for w in s.log.warnings)


def test_supersaturated_declared_gas_is_still_flagged(engine):
    r = engine.solve({"T": 300.0, "P": "1 bar", "components": {"water": "gas"}})
    assert any("below its dew point" in w for w in r.warnings)


# ---- G3: one vapour model for equilibrium and enthalpy --------------------------------------------
def test_vle_uses_the_gas_equation_of_state_by_default(engine):
    s = engine.system({"mole_fractions": {"ethanol": 0.3, "water": 0.7}})
    s.configure(350.0, 360.0)
    v = next(iter(s._vle_cache.values()))
    assert v.vapour == "eos" and v._eos_mix is not None
    s2 = engine.system({"mole_fractions": {"ethanol": 0.3, "water": 0.7}, "gas_eos": "ideal"})
    s2.configure(350.0, 360.0)
    assert next(iter(s2._vle_cache.values())).vapour == "ideal"


def test_eos_flash_is_an_equilibrium(engine):
    import numpy as np
    s = engine.system({"mole_fractions": {"ethanol": 0.3, "water": 0.7}})
    s.configure(350.0, 360.0)
    v = next(iter(s._vle_cache.values()))
    V, x, y = v.flash(355.0, 101325.0, [0.3, 0.7])
    assert np.allclose((1 - V) * x + V * y, [0.3, 0.7], atol=1e-10)
    K = v.K(355.0, 101325.0, x, y)
    assert np.allclose(y, K * x / np.sum(K * x), rtol=1e-8)


# ---- G5: interval-mean breakdown adds up exactly to ΔH/ΔT --------------------------------------------
def test_mean_breakdown_sums_to_dH_over_dT(engine):
    s = engine.system({"mole_fractions": {"ethanol": 0.3, "water": 0.7}})
    bd = s.cp_mean(290.0, 340.0)
    assert bd.total == pytest.approx((s.enthalpy(340.0) - s.enthalpy(290.0)) / 50.0, rel=1e-13)
    assert bd.latent == 0.0 and abs(bd.extras["quadrature_mismatch_J_per_K"]) < 1e-6 * bd.total


# ---- Gaussian tails cut at ±8σ ------------------------------------------------------------------------
def test_gaussian_melting_is_cut_at_eight_sigma():
    Tm, width = 301.3, 2.0
    sig = width / 4
    assert melt_fraction(Tm - TAIL_SIGMAS * sig - 1e-9, Tm, width) == 0.0
    assert melt_fraction(Tm + TAIL_SIGMAS * sig + 1e-9, Tm, width) == 1.0
    assert melt_fraction_derivative(Tm + 5.0, Tm, width) == 0.0
    inside = melt_fraction(Tm - TAIL_SIGMAS * sig + 1e-9, Tm, width)
    assert inside < 1e-14                                      # the cut removes < 1.2e-15 of the transition


def test_pcm_energy_still_conserved_after_tail_cut(engine):
    from scipy import integrate
    s = engine.system({"mass_fractions": {"water": 0.8, "n-octadecane": 0.2}})
    s.configure(290.0, 312.0)
    val = integrate.quad(lambda t: s.cp_point(float(t)).total, 290.0, 312.0, points=[301.3], limit=400,
                         epsrel=1e-11)[0]
    assert val == pytest.approx(s.enthalpy(312.0) - s.enthalpy(290.0), rel=1e-8)
    assert math.isfinite(val)
