"""A small log object that collects warnings and assumptions during a calculation.

Thermodynamic models make assumptions (ideal volume of mixing, constant Cp
over an interval, ...) and may extrapolate correlations. Instead of printing
these, every function appends them to a Log so they can be shown together in
the final result.
"""
from __future__ import annotations

from dataclasses import dataclass, field

# A correlation may be used outside its stated range by up to this many kelvin
# (with a warning). Further outside, the engine refuses (RangeError).
EXTRAPOLATION_LIMIT_K = 30.0


@dataclass
class Log:
    warnings: list[str] = field(default_factory=list)
    assumptions: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def warn(self, msg: str) -> None:
        if msg not in self.warnings:
            self.warnings.append(msg)

    def assume(self, msg: str) -> None:
        if msg not in self.assumptions:
            self.assumptions.append(msg)

    def note(self, msg: str) -> None:
        if msg not in self.notes:
            self.notes.append(msg)

    def extend(self, other: "Log") -> None:
        for m in other.warnings:
            self.warn(m)
        for m in other.assumptions:
            self.assume(m)
        for m in other.notes:
            self.note(m)


def null_log() -> Log:
    """A throw-away log for callers that do not care about messages."""
    return Log()
