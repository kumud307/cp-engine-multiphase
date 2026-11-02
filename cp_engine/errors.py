"""Custom exceptions for the Cp engine.

Every error caused by *user input* is raised as a subclass of CpEngineError
so that the top-level solver can catch it and turn it into a clear message
instead of a raw Python traceback.
"""


class CpEngineError(Exception):
    """Base class for all errors raised deliberately by the engine."""


class UnitError(CpEngineError):
    """A unit string could not be parsed or has the wrong dimension."""


class InputError(CpEngineError):
    """The user input is malformed (bad fractions, negative amounts, ...)."""


class DataError(CpEngineError):
    """Required data are missing from the component / binary database."""


class PhaseError(CpEngineError):
    """A declared phase is physically impossible at the given T and P."""


class RangeError(CpEngineError):
    """Temperature is far outside the valid range of a correlation."""
