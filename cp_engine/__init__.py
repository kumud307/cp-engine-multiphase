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
from .database import Database, load_database
from .errors import CpEngineError, DataError, InputError, PhaseError, RangeError, UnitError
from .inputs import SystemInput, describe_inputs
from .solver import CpEngine, CpResult, solve
from .system import SystemModel

__all__ = ["CpEngine", "CpResult", "solve", "SystemInput", "describe_inputs", "Database", "load_database",
           "SystemModel", "CpEngineError", "DataError", "InputError", "PhaseError", "RangeError", "UnitError"]
__version__ = "1.0.0"
