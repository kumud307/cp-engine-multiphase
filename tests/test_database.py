"""Database loading and pure-component correlations (vs NIST reference values)."""
import pytest

from cp_engine.errors import DataError, RangeError
from cp_engine.report import Log


@pytest.mark.parametrize("name, phase, T, cp_ref, tol", [
    # NIST SRD 69 isobaric data at 0.101325 MPa (validation/data) – J/(kg K)
    ("water", "liquid", 300.0, 4180.6, 0.002),
    ("methanol", "liquid", 300.0, 2546.0, 0.003),
    ("n-hexane", "liquid", 300.0, 2280.3, 0.003),
    ("toluene", "liquid", 300.0, 1707.0, 0.003),
    ("n-dodecane", "liquid", 300.0, 2218.3, 0.004),
])
def test_liquid_cp_matches_nist(db, name, phase, T, cp_ref, tol):
    assert db.get(name).cp_mass(phase, T) == pytest.approx(cp_ref, rel=tol)


@pytest.mark.parametrize("name, cp298", [
    # NIST WebBook Cp(298.15 K) in J/(mol K) for the Shomate solids
    ("Al2O3", 79.01), ("CuO", 42.30), ("SiO2", 44.6), ("NaCl", 50.5),
])
def test_shomate_solids(db, name, cp298):
    assert db.get(name).cp("solid", 298.15) == pytest.approx(cp298, rel=0.005)


def test_ethanol_liquid_cp(db):
    # NIST WebBook: 112.3–112.4 J/(mol K) at 298.15 K
    assert db.get("ethanol").cp("liquid", 298.15) == pytest.approx(112.35, rel=0.003)


def test_water_density_and_expansion(db):
    w = db.get("water")
    assert w.density("liquid", 298.15) == pytest.approx(997.05, rel=1e-4)
    assert w.alpha("liquid", 298.15) == pytest.approx(2.57e-4, rel=0.01)
    assert abs(w.alpha("liquid", 277.13)) < 2e-6      # density maximum near 4 °C


def test_aliases_and_phase_aliases(db):
    comp, ph = db.resolve("ice")
    assert comp.key == "water" and ph == "solid"
    assert db.get("Alumina").key == "Al2O3"
    assert db.get("glycerine").key == "glycerol"


def test_unknown_component(db):
    with pytest.raises(DataError, match="Unknown component"):
        db.get("unobtainium")


def test_missing_phase_data(db):
    with pytest.raises(DataError, match="No solid-phase Cp"):
        db.get("ethanol").cp("solid", 150.0)


def test_range_warning_and_error(db):
    log = Log()
    db.get("water").cp("liquid", 265.0, log)              # 8 K below range: warn
    assert any("outside the correlation range" in w for w in log.warnings)
    with pytest.raises(RangeError):
        db.get("water").cp("liquid", 200.0, Log())         # 73 K below range: refuse


def test_clausius_clapeyron_boiling_point(db):
    w = db.get("water")
    assert w.Tb_at(101325.0) == pytest.approx(373.12)
    assert w.Tb_at(2 * 101325.0) == pytest.approx(393.4, abs=2.0)   # steam tables: 393.4 K


def test_every_component_has_mw_and_some_cp(db):
    for k, c in db.components.items():
        assert c.MW > 0 and c.cp_corr, k


def test_verify_list_not_empty(db):
    assert any("binary.nrtl" in v for v in db.verify_items())
