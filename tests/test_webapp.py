"""The web interface (webapp/) must give the same numbers as the engine called directly."""
from __future__ import annotations

import json
import sys
import threading
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "webapp"))
import server as web  # noqa: E402

from cp_engine import CpEngine  # noqa: E402

BASIS_OF = {"mass_fractions": "mass fractions", "mole_fractions": "mole fractions",
            "volume_fractions": "volume fractions", "masses": "masses (kg)", "moles": "moles (mol)"}


def form_body(inp: dict) -> dict:
    """What the page sends after an example has been filled into the form (mirrors fillForm/readForm)."""
    comps = inp.get("components") or []
    phases = comps if isinstance(comps, dict) else {}
    names, amounts, basis = list(comps), {}, "mass fractions"
    for key, label in BASIS_OF.items():
        if key in inp:
            basis, amounts = label, inp[key]
            names += [n for n in inp[key] if n not in names]
    rows = [{"name": n, "phase": phases.get(n, "auto"), "amount": "" if amounts.get(n) is None else str(amounts[n])}
            for n in names]
    skip = {"components", *BASIS_OF}
    return {"rows": rows, "basis": basis, "fields": {k: str(v) for k, v in inp.items() if k not in skip}}


@pytest.mark.parametrize("preset", web.PRESETS, ids=lambda p: p["id"])
def test_examples_match_the_engine(preset):
    direct = CpEngine().solve(preset["input"])
    out = web.solve(form_body(preset["input"]))
    json.dumps(out)                                           # the answer must be valid JSON
    assert out["success"] == direct.success
    if direct.success:
        assert out["cp_mass_J_per_kgK"] == pytest.approx(direct.cp_mass, rel=1e-12)
        assert out["report"] and "derivation" in out
    else:
        assert out["errors"] or out["missing"]


def test_options_reach_the_engine():
    body = form_body({"T": "25 C", "mole_fractions": {"ethanol": 0.3, "water": 0.7}})
    real = web.solve(body)["cp_molar_J_per_molK"]
    body["fields"]["solution_model"] = "ideal"
    ideal = web.solve(body)["cp_molar_J_per_molK"]
    assert real == pytest.approx(104.81, abs=0.01) and ideal == pytest.approx(86.47, abs=0.01)


def test_curve_shows_the_melting_peak_and_bad_ranges_are_explained():
    body = form_body({"mass_fractions": {"water": 0.9, "ethanol": 0.1}})
    out = web.curve({**body, "T1": "-15 C", "T2": "5 C", "points": 21})
    cps = [p["cp"] for p in out["points"]]
    assert len(cps) == 21 and max(cps) > 10 * min(cps)
    assert "error" in web.curve({**body, "T1": "", "T2": "5 C"})
    assert "error" in web.curve({**body, "T1": "abc", "T2": "5 C"})


def test_http_endpoints():
    srv = ThreadingHTTPServer(("127.0.0.1", 0), web.Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{srv.server_address[1]}"
    try:
        assert b"Cp Engine" in urllib.request.urlopen(base + "/").read()
        meta = json.loads(urllib.request.urlopen(base + "/api/meta").read())
        assert len(meta["components"]) >= 13 and meta["presets"]
        req = urllib.request.Request(base + "/api/solve", data=json.dumps({"input": {"T": "25 C", "components": ["water"]}}).encode(),
                                     headers={"Content-Type": "application/json"})
        assert json.loads(urllib.request.urlopen(req).read())["cp_mass_J_per_kgK"] == pytest.approx(4184.5, abs=0.1)
        bad = urllib.request.Request(base + "/api/solve", data=b"not json", headers={"Content-Type": "application/json"})
        with pytest.raises(urllib.error.HTTPError) as err:
            urllib.request.urlopen(bad)
        assert err.value.code == 400
    finally:
        srv.shutdown()
