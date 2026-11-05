"""Pytest wrapper around the independent audit checks (audit/independent_checks.py).

Run:  pytest audit/test_independent.py -q

The defects found by the audit (G1 false condensation guard, G2 non-finite inputs) were fixed after the
audit; their tests are now ordinary tests. Tolerance-only differences that the report explains as numerical
(D1 at 300.3 K, FD at the PCM peak with ΔT = 0.01 K) are checked with tolerances justified in
MODEL_VERIFICATION_REPORT.md §E.
"""
from __future__ import annotations

import math

import pytest

import independent_checks as ic

ic.RESULTS.clear()
for fn in (ic.test_A, ic.test_B, ic.test_C, ic.test_D, ic.test_E, ic.test_F, ic.test_FD, ic.test_departure):
    fn()
BY_NAME = {r["test"]: r for r in ic.RESULTS}

NUMERICAL = {  # name: justified relative tolerance (%)
    "D1 PCM slurry apparent Cp at 300.3 K": 1e-3,                 # O(h²) of the latent finite difference, h = 1e-3 K
    "FD PCM slurry at Tm: ΔT = 0.01 K": 1e-2,                     # truncation of the REFERENCE: ΔT²/(6σ²) = 6.7e-3 %
}


@pytest.mark.parametrize("name", sorted(BY_NAME))
def test_independent_check(name):
    r = BY_NAME[name]
    if name in NUMERICAL:
        assert r["rel_error_%"] <= NUMERICAL[name], r
        return
    assert r["pass"], r


def test_nan_temperature_is_rejected():          # was defect G2
    r = ic.E.solve({"T": float("nan"), "components": ["water"]})
    assert not r.success


def test_no_false_condensation_warning_in_vle_region():          # was defect G1
    r = ic.E.solve({"T": 362.0, "mole_fractions": {"ethanol": 0.1, "water": 0.9}})
    assert not any("below its dew point" in w for w in r.warnings)


def test_reported_cp_is_finite_for_valid_inputs():
    for inp in ({"T": 300.0, "components": ["water"]},
                {"T": 357.0, "mole_fractions": {"ethanol": 0.3, "water": 0.7}},
                {"T": 301.3, "mass_fractions": {"water": 0.8, "n-octadecane": 0.2}}):
        assert math.isfinite(ic.E.solve(inp).cp_mass)
