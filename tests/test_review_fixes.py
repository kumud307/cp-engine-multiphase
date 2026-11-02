"""Tests for the improvements made after the mock examiner review (pressure, electrolytes,
consistent ice-slurry model, PCM width reporting, user-supplied components)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cp_engine.models import phase_change as pc
from cp_engine.models.nrtl import NRTL
from cp_engine.system import SystemModel

DATA = Path(__file__).resolve().parent.parent / "validation" / "data"


# ---------------------------------------------------------------- pressure
@pytest.mark.parametrize("fluid", ["water", "methanol", "toluene", "n-hexane", "n-dodecane"])
def test_pressure_dependence_matches_nist(engine, fluid):
    d = pd.read_csv(DATA / f"nist_{fluid}_298K_isotherm.csv")
    pred = [engine.solve({"T": 298.15, "P": (P, "MPa"), "components": {fluid: "liquid"}}).cp_mass for P in d.P_MPa]
    assert np.max(np.abs(np.array(pred) - d.cp_J_kgK) / d.cp_J_kgK) * 100 < 1.5      # up to 100 MPa


def test_pressure_point_and_mean_consistent(engine):
    base = {"P": "300 bar", "mass_fractions": {"water": 0.6, "ethanol": 0.4}}
    r = engine.solve(dict(base, T1=290, T2=310, m=1))
    mean = [x for x in r.routes if x.tag == "A:mean"][0].cp_mass
    c = engine.cp_curve(base, np.linspace(290, 310, 401))
    assert np.trapezoid(c.total, c.T_K) / 20 == pytest.approx(mean, rel=1e-5)


def test_ideal_gas_has_no_pressure_correction(db):
    assert db.get("water").cp_pressure_correction("gas", 400.0, 5e5) == 0.0


def test_high_pressure_warning(engine):
    r = engine.solve({"T": 300, "P": "2000 bar", "components": {"water": "liquid"}})
    assert any("beyond the 100 MPa" in w for w in r.warnings)


# ---------------------------------------------------------------- electrolytes
def test_nacl_brine_independent_point(engine):
    """doi:10.1016/j.jct.2016.09.031 (not used in the fit): x(NaCl) = 0.0356, 298.15 K, Cp = 72.32 J/(mol K)."""
    r = engine.solve({"T": 298.15, "mole_fractions": {"NaCl": 0.0356, "water": 1 - 0.0356}})
    assert r.cp_molar == pytest.approx(72.32, rel=0.003)
    assert r.breakdown.excess < 0                       # dissolved NaCl lowers Cp below Neumann–Kopp


def test_brine_point_and_mean_consistent(engine):
    base = {"mass_fractions": {"water": 0.9, "NaCl": 0.1}}
    r = engine.solve(dict(base, T1=293.15, T2=303.15, m=1))
    mean = [x for x in r.routes if x.tag == "A:mean"][0].cp_mass
    c = engine.cp_curve(base, np.linspace(293.15, 303.15, 201))
    assert np.trapezoid(c.total, c.T_K) / 10 == pytest.approx(mean, rel=1e-5)


# ---------------------------------------------------------------- consistent ice slurry
def test_engine_liquidus_uses_same_nrtl_as_HE(db):
    w = db.get("water")
    dCp = w.cp_corr["liquid"].raw(w.Tm) - w.cp_corr["solid"].raw(w.Tm)
    m = NRTL(["water", "ethanol"], db)
    assert "SLE" in m.pair_sources[("water", "ethanol")]["fitted_to"]
    n = {"water": 0.8 / w.MW, "ethanol": 0.2 / db.get("ethanol").MW}
    x = n["water"] / sum(n.values())
    s = SystemModel(db, n, {})
    s.configure(220.0, 273.15)
    from scipy.optimize import brentq
    T_eng = brentq(lambda T: sum(q for a, b, st, q in s.split(T) if st == "solid") - 1e-9 * n["water"], 220, 273.14)
    T_eq = pc.liquidus_temperature(x, w.Tm, w.dHfus, lambda T: m.ln_gamma([x, 1 - x], T)[0], dCp)
    assert T_eng == pytest.approx(T_eq, abs=1e-6)


# ---------------------------------------------------------------- PCM width
def test_pcm_width_independent_quantities(engine):
    out = []
    for width in (0.5, 2.0, 5.0):
        r = engine.solve({"T": "301.3 K", "mass_fractions": {"water": 0.8, "n-octadecane": 0.2},
                          "melt_width": (width, "K")})
        out.append((r.cp_mass, r.breakdown.extras["transition_latent_J_per_kg"],
                    r.breakdown.extras["transition_mean_cp_mass"]))
        assert any("depends on the assumed width" in w for w in r.warnings)
    assert out[0][0] > 5 * out[2][0]                            # point value depends on width ...
    assert out[0][1] == pytest.approx(out[2][1])                # ... latent heat does not
    assert out[0][2] == pytest.approx(out[2][2], rel=1e-3)      # ... nor the fixed-window mean


# ---------------------------------------------------------------- user-supplied data
def test_custom_component(engine):
    r = engine.solve({"T": "25 C", "mass_fractions": {"water": 0.5, "ethylene glycol": 0.5},
                      "custom_components": {"ethylene glycol": {"MW": "62.068 g/mol", "cp": "2.43 kJ/(kg K)",
                                                                "density": "1110 kg/m3", "miscibility_class": "polar"}}})
    cp_w = engine.db.get("water").cp_mass("liquid", 298.15)
    assert r.cp_mass == pytest.approx(0.5 * cp_w + 0.5 * 2430.0, rel=1e-9)
    with pytest.raises(Exception):
        engine.db.get("ethylene glycol")                         # database itself unchanged


def test_pure_cp_override(engine):
    r = engine.solve({"T": "25 C", "components": ["water"], "pure_cp": {"water": "4.0 kJ/(kg K)"}})
    assert r.cp_mass == pytest.approx(4000.0)
    assert engine.solve({"T": "25 C", "components": ["water"]}).cp_mass == pytest.approx(4184.5, rel=1e-3)


def test_custom_component_errors(engine):
    r = engine.solve({"T": 300, "mass_fractions": {"water": 0.6, "x": 0.4}, "custom_components": {"x": {"MW": "50 g/mol"}}})
    assert not r.success and "needs at least" in r.errors[0]
    r = engine.solve({"T": 300, "components": ["water"], "custom_components": {"ethanol": {"MW": 46, "cp": 2000}}})
    assert not r.success and "clashes" in r.errors[0]
