"""Web interface for the Cp engine.

Run:  python webapp/server.py        then open http://localhost:8080
On Replit the Run button does this (see .replit in the project root).

The server uses only the Python standard library, so the only packages needed are the engine's own
(numpy, scipy, pandas). If they are missing they are installed on first start (webapp/requirements.txt).
The engine itself (cp_engine/) is imported unchanged: this folder only adds a user interface.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import threading
import time
import warnings
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))


def _import_engine():
    from cp_engine import CpEngine                      # noqa: WPS433
    from cp_engine.interactive import BASES, build_input
    return CpEngine, BASES, build_input


try:
    CpEngine, BASES, build_input = _import_engine()
except ImportError as exc:                              # fresh machine (e.g. a new Repl): install and retry
    print(f"[setup] {exc}. Installing the packages in webapp/requirements.txt (one time, about a minute)...",
          flush=True)
    base = [sys.executable, "-m", "pip", "install", "-q", "-r", str(HERE / "requirements.txt")]
    for extra in ([], ["--user"], ["--break-system-packages"]):       # hosts differ in what pip allows
        if subprocess.call(base + extra) == 0:
            break
    import importlib
    import site
    importlib.invalidate_caches()
    if site.ENABLE_USER_SITE and site.getusersitepackages() not in sys.path:
        sys.path.append(site.getusersitepackages())
    CpEngine, BASES, build_input = _import_engine()

warnings.filterwarnings("ignore")
ENGINE = CpEngine()
LOCK = threading.Lock()                                 # the engine is not written for concurrent use
MAX_BODY = 64 * 1024
MAX_CURVE_POINTS = 61

# Examples offered in the interface. `concept` names the CL 202 lecture idea behind each one.
PRESETS = [
    {"id": "water", "label": "Water, 25 °C",
     "concept": "Pure component: Cp from property data (Route A). Cv follows from Cp − Cv = Tvα²/κT.",
     "input": {"T": "25 C", "components": ["water"]}},
    {"id": "calorimeter", "label": "Calorimeter (first law)",
     "concept": "Ch 1, first law at constant pressure: ΔH = Q + Ws, so Cp = Q/(mΔT). No property data used.",
     "input": {"T1": "20 C", "T2": "30 C", "total_mass": "2 kg", "process": "constant_P", "Q": "83.7 kJ"}},
    {"id": "solution", "label": "Ethanol + water",
     "concept": "Ch 4: excess functions. Cp = Σ xᵢCpᵢ + Cpᴱ, with Cpᴱ from NRTL through Gibbs–Helmholtz. "
                "Switch the solution model to “ideal” to see the 18 % difference.",
     "input": {"T": "25 C", "mole_fractions": {"ethanol": 0.3, "water": 0.7}}},
    {"id": "nanofluid", "label": "Nanofluid (Al₂O₃ in water)",
     "concept": "Enthalpy adds over phases, so a suspension is mass-weighted, not volume-weighted.",
     "input": {"T": "25 C", "volume_fractions": {"water": 0.96, "Al2O3": 0.04}}},
    {"id": "slurry", "label": "Ice slurry, −5 °C",
     "concept": "Apparent Cp = (∂H/∂T)P. Ice melts as T rises, so the latent term (∂H/∂n)(dn/dT) dominates. "
                "Plot Cp against T from −15 to 5 °C to see the peak.",
     "input": {"T": "-5 C", "mass_fractions": {"water": 0.9, "ethanol": 0.1}},
     "curve": {"T1": "-15 C", "T2": "5 C"}},
    {"id": "pcm", "label": "PCM slurry (octadecane)",
     "concept": "Phase-change material melting at 28 °C: mean Cp over an interval that contains the melting.",
     "input": {"T1": "295 K", "T2": "308 K", "mass_fractions": {"water": 0.8, "n-octadecane": 0.2}},
     "curve": {"T1": "295 K", "T2": "308 K"}},
    {"id": "steam", "label": "Steam, 500 °C, 100 bar",
     "concept": "Equations of state: Cp = Cp(ideal gas) + residual Cpᴿ from Peng–Robinson. "
                "Change the gas model to “ideal” or “virial” to compare (NIST: 2583 J/(kg K)).",
     "input": {"T": "500 C", "P": "100 bar", "components": {"water": "gas"}}},
    {"id": "boiling", "label": "Boiling mixture, 82 °C",
     "concept": "γ–φ vapour–liquid equilibrium: bubble and dew points and a Rachford–Rice flash decide the phases.",
     "input": {"T": "82 C", "mole_fractions": {"ethanol": 0.3, "water": 0.7}},
     "curve": {"T1": "70 C", "T2": "100 C"}},
    {"id": "pressure", "label": "Water at 500 bar",
     "concept": "(∂Cp/∂P)T = −T(∂²v/∂T²)P with the Tait equation for v(T, P) (NIST: 4060 J/(kg K)).",
     "input": {"T": "25 C", "P": "500 bar", "components": ["water"]}},
    {"id": "bad", "label": "Wrong input",
     "concept": "The engine explains what is wrong instead of failing.",
     "input": {"T": "25 C", "mass_fractions": {"water": 0.7, "ethanol": 0.6}}},
]


def _clean(obj):
    """Make a result JSON-safe: numpy numbers -> float, NaN/inf -> None, tuples -> lists."""
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, (str, bool)) or obj is None:
        return obj
    try:
        f = float(obj)
    except (TypeError, ValueError):
        return str(obj)
    return f if math.isfinite(f) else None


def _input_from_form(body: dict) -> dict:
    if isinstance(body.get("input"), dict):                  # a preset or hand-written input
        return body["input"]
    rows = [(str(r.get("name", "")).strip(), str(r.get("phase") or "auto"), str(r.get("amount", "")))
            for r in body.get("rows", []) if isinstance(r, dict)]
    basis = body.get("basis") if body.get("basis") in BASES else "mass fractions"
    fields = {k: str(v) for k, v in (body.get("fields") or {}).items()}
    options = {k: fields.pop(k) for k in ("solution_model", "gas_eos") if k in fields}
    inp = build_input(rows, basis, fields)
    if inp.get("phases") and all(p == "auto" for p in inp["phases"].values()):
        del inp["phases"]                                    # keep the echoed input short
    for k, v in options.items():
        if v and v not in ("auto", "default"):
            inp[k] = v
    return inp


def solve(body: dict) -> dict:
    inp = _input_from_form(body)
    t0 = time.time()
    with LOCK:
        res = ENGINE.solve(inp)
    out = res.to_dict()
    out["report"] = str(res)
    try:
        out["derivation"] = res.derivation() if res.success else ""
    except Exception:                                        # noqa: BLE001  (never fail the whole answer for this)
        out["derivation"] = ""
    bd = res.breakdown
    if bd is not None:
        out["phase_amounts"] = bd.phase_amounts
        out["breakdown_kind"] = bd.kind
        out["breakdown_T"] = [bd.T, bd.T1, bd.T2]
        keep = ("partial_molar_cp_J_per_molK", "mixing_functions_J_per_mol", "vle", "rho_mix")
        out["extras"] = {k: v for k, v in bd.extras.items() if k in keep}
    out["input"] = inp
    out["seconds"] = round(time.time() - t0, 2)
    return _clean(out)


def curve(body: dict) -> dict:
    """Cp(T) at evenly spaced temperatures between T1 and T2 for the given composition (Route A)."""
    from cp_engine.units import to_si                         # noqa: WPS433

    def kelvin(q):
        if q in (None, ""):
            return None
        try:
            q = float(q)                                      # a plain number is taken as kelvin
        except (TypeError, ValueError):
            pass
        return float(to_si(q, "temperature", absolute_temperature=True)[0])

    inp = dict(_input_from_form(body))
    for k in ("T", "T1", "T2", "dT", "Q", "W", "dH", "dS", "Cv", "process", "total_mass", "m", "n", "gamma"):
        inp.pop(k, None)
    try:
        T1, T2 = kelvin(body.get("T1")), kelvin(body.get("T2"))
    except Exception as exc:                                  # noqa: BLE001
        return {"error": f"Could not read the temperature range: {exc}"}
    if T1 is None or T2 is None or not (math.isfinite(T1) and math.isfinite(T2)) or T1 == T2:
        return {"error": "Give two different temperatures (for example -15 C and 5 C) for the plot."}
    T1, T2 = sorted((T1, T2))
    n = max(5, min(int(body.get("points") or 41), MAX_CURVE_POINTS))
    pts, failed, first_error = [], 0, None
    with LOCK:
        for i in range(n):
            T = T1 + (T2 - T1) * i / (n - 1)
            res = ENGINE.solve({**inp, "T": T})
            if res.success and res.cp_mass is not None and math.isfinite(res.cp_mass):
                b = res.breakdown.per_mass() if res.breakdown else {}
                pts.append({"T": T, "cp": res.cp_mass, "latent": b.get("latent"), "ideal": b.get("ideal"),
                            "excess": b.get("excess")})
            else:
                failed += 1
                first_error = first_error or (res.errors[0] if res.errors else "no route to Cp")
    return _clean({"points": pts, "failed": failed, "first_error": first_error, "T1": T1, "T2": T2})


def meta() -> dict:
    db = ENGINE.db
    comps = []
    for name in db.names():
        c = db.get(name)
        comps.append({"key": name, "name": getattr(c, "name", name)})
    return {"components": comps, "bases": list(BASES), "presets": PRESETS}


class Handler(BaseHTTPRequestHandler):
    server_version = "CpEngineWeb/1.0"

    def _send(self, code: int, body: bytes, ctype: str) -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, code: int = 200) -> None:
        self._send(code, json.dumps(obj, ensure_ascii=False).encode("utf-8"), "application/json; charset=utf-8")

    def do_GET(self):                                         # noqa: N802
        path = self.path.split("?")[0]
        if path in ("/", "/index.html"):
            self._send(200, (HERE / "static" / "index.html").read_bytes(), "text/html; charset=utf-8")
        elif path == "/api/meta":
            self._json(meta())
        elif path == "/health":
            self._json({"ok": True})
        else:
            self._json({"error": "not found"}, 404)

    def do_POST(self):                                        # noqa: N802
        path = self.path.split("?")[0]
        if path not in ("/api/solve", "/api/curve"):
            return self._json({"error": "not found"}, 404)
        try:
            length = int(self.headers.get("Content-Length") or 0)
            if length > MAX_BODY:
                return self._json({"error": "request too large"}, 413)
            body = json.loads(self.rfile.read(length).decode("utf-8") or "{}")
            if not isinstance(body, dict):
                raise ValueError("the request must be a JSON object")
        except (ValueError, UnicodeDecodeError) as exc:
            return self._json({"error": f"bad request: {exc}"}, 400)
        try:
            self._json(solve(body) if path == "/api/solve" else curve(body))
        except Exception as exc:                              # noqa: BLE001  (show the reason in the page)
            self._json({"error": f"{type(exc).__name__}: {exc}"}, 500)

    def log_message(self, fmt, *args):                        # quieter console
        if "/api/" in (args[0] if args else ""):
            sys.stderr.write("%s %s\n" % (self.log_date_time_string(), fmt % args))


def main() -> None:
    port = int(os.environ.get("PORT", "8080"))
    srv = ThreadingHTTPServer(("0.0.0.0", port), Handler)
    print(f"Cp engine web interface running on port {port}  (local: http://localhost:{port})", flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
