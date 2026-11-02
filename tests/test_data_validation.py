"""Regression tests against the real experimental data in validation/data.

If a future change to the code or the database makes the engine disagree with the literature
data, these tests fail. Thresholds are set slightly above the values achieved on 2026-11-02.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cp_engine.models import phase_change as pc
from cp_engine.models.heterogeneous import cp_model_I
from cp_engine.models.nrtl import NRTL

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "validation" / "data"
sys.path.insert(0, str(ROOT / "tools"))


def binary(a, b, prop):
    d = pd.read_csv(DATA / "thermoml_binary_data.csv", comment="#")
    g = d[(((d.comp1 == a) & (d.comp2 == b)) | ((d.comp1 == b) & (d.comp2 == a))) & (d.property == prop)].copy()
    g["x"] = np.where(g.comp1 == a, g.x1, 1 - g.x1)
    return g


def test_parameter_set_selection(db):
    assert db.nrtl_pair("ethanol", "water")["fitted_to"] == "HE+SLE"        # default: H^E + freezing points
    assert db.nrtl_pair("methanol", "water")["fitted_to"] == "HE"
    assert db.nrtl_pair("ethanol", "water", prefer="VLE")["fitted_to"] == "VLE"
    assert db.nrtl_pair("water", "ethanol")["a_ij"] == db.nrtl_pair("ethanol", "water")["a_ji"]   # orientation


@pytest.mark.parametrize("alc", ["ethanol", "methanol"])
def test_he_fitted_sets_are_stable_and_accurate(db, alc):
    m = NRTL([alc, "water"], db, prefer="HE")
    assert all(m.is_stable_binary(T, 500) for T in range(240, 391, 25))     # no false liquid-liquid split
    he = binary(alc, "water", "HE")
    he = he[(he.x > 0.01) & (he.x < 0.99) & (he.T_K <= 383.2)]
    rmsd = np.sqrt(np.mean([(m.hE([x, 1 - x], T) - v) ** 2 for x, T, v in zip(he.x, he.T_K, he.value)]))
    assert rmsd < 36.0                                                     # J/mol (VLE sets: 140-260)


def test_methanol_water_cp_independent_data(engine):
    cp = binary("methanol", "water", "Cp_molar")
    cp = cp[cp.P_kPa <= 200]
    pred = [engine.solve({"T": t, "mole_fractions": {"methanol": x, "water": 1 - x}}).cp_molar
            for x, t in zip(cp.x, cp.T_K)]
    assert np.mean(np.abs(np.array(pred) - cp.value) / cp.value) * 100 < 3.0      # ideal mixture: 7.5 %


def test_nanofluids_model_II_beats_model_I(engine, db):
    d = pd.read_csv(DATA / "nanofluid_ohanley2012.csv", comment="#")
    m2 = np.array([engine.solve({"T": r.T_K, "mass_fractions": {r.particle: r.w_table2, "water": None}}).cp_mass
                   for r in d.itertuples()])
    m1 = np.array([cp_model_I(r.phi_table2, db.get(r.particle).cp_mass("solid", r.T_K),
                              db.get("water").cp_mass("liquid", r.T_K)) for r in d.itertuples()])
    aard = lambda p: np.mean(np.abs(p - d.cp_J_kgK) / d.cp_J_kgK) * 100
    assert aard(m2) < 2.5 and aard(m1) > 10.0


def test_freezing_points_ethanol_water(db):
    d = pd.read_csv(DATA / "ethanol_water_freezing_points.csv", comment="#")
    d = d[(d.branch == "ice") & (d.w_ethanol_percent > 0) & (d.w_ethanol_percent <= 30)]
    w = db.get("water")
    dCp = w.cp_corr["liquid"].raw(w.Tm) - w.cp_corr["solid"].raw(w.Tm)
    nw, ne = (100 - d.w_ethanol_percent) / w.MW, d.w_ethanol_percent / db.get("ethanol").MW
    T = [pc.ideal_freezing_point(x, w.Tm, w.dHfus, dCp) for x in nw / (nw + ne)]
    assert np.mean(np.abs(np.array(T) - d.freezing_point_K)) < 1.0                    # K


def test_acetone_water_orientation(db):
    """γ∞ of acetone in water ≈ 6 and of water in acetone ≈ 3.6 at 298 K (corrected 2026-11-02)."""
    m = NRTL(["acetone", "water"], db, prefer="VLE")
    assert np.exp(m.ln_gamma([1e-9, 1 - 1e-9], 298.15)[0]) == pytest.approx(5.83, rel=0.01)
    assert np.exp(m.ln_gamma([1 - 1e-9, 1e-9], 298.15)[1]) == pytest.approx(3.58, rel=0.01)


def test_iapws_ice_check_values():
    from iapws_ice import props
    r = props(273.16, 611.657)                                   # IAPWS R10-06 Table 6
    assert r["cp"] == pytest.approx(0.209678431622e4, rel=1e-9)
    assert r["rho"] == pytest.approx(0.916709492200e3, rel=1e-9)
    assert r["kappa_T"] == pytest.approx(0.117793449348e-9, rel=1e-8)
