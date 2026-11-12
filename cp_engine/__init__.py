"""cp_engine – general specific-heat (Cp) engine for multiphase fluids.

CL 202 Chemical Engineering Thermodynamics, IIT Gandhinagar – Topic 5.

Quick start
-----------
>>> from cp_engine import CpEngine
>>> engine = CpEngine()
>>> result = engine.solve({
...     "T": "25 C",
...     "components": ["water", "ethanol"],
...     "mass_fractions": {"water": 0.6, "ethanol": 0.4},
... })
>>> print(result)            # report: Cp, routes, breakdown, warnings
>>> print(result.derivation())   # chain of thermodynamic relations used
"""
import sys as _sys

# Reports contain Greek letters and symbols (α, γ, Δ, °). On Windows, output that is redirected to a file or
# a pipe uses the legacy code page (cp1252) and print() would raise UnicodeEncodeError. Use UTF-8 instead.
for _stream in (_sys.stdout, _sys.stderr):
    try:
        if _stream is not None and (_stream.encoding or "").lower().replace("-", "") != "utf8":
            _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError, OSError):   # stream replaced by a notebook or test harness
        pass

from .database import Database, load_database
from .errors import CpEngineError, DataError, InputError, PhaseError, RangeError, UnitError
from .inputs import SystemInput, describe_inputs
from .solver import CpEngine, CpResult, solve
from .system import SystemModel

__all__ = ["CpEngine", "CpResult", "solve", "SystemInput", "describe_inputs", "Database", "load_database",
           "SystemModel", "CpEngineError", "DataError", "InputError", "PhaseError", "RangeError", "UnitError"]
__version__ = "1.0.0"
