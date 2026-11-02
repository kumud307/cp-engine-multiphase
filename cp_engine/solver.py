"""Forward-chaining solver and the public CpEngine API.

Algorithm (forward chaining)
----------------------------
1. Put every quantity the user supplied into a knowledge base (KB).
2. Repeat: for every rule whose inputs are all known and which has not
   fired yet, compute its outputs and add them to the KB (recording which
   rule produced them from which inputs).
   - If an output is already known and the rule is a 'check' rule, the new
     value is compared with the known one → consistency warning if they
     differ (over-specified input, e.g. Q and ΔH at constant P).
3. Stop when a full pass fires no rule.
4. Every Cp candidate in the KB is a route that reached Cp.  The chain of
   rules behind each candidate is reconstructed from the provenance record.
5. If no candidate exists, search backwards through the rules for the
   smallest sets of additional inputs that would complete a route.
"""
from __future__ import annotations

import itertools
import textwrap
from dataclasses import dataclass, field
from typing import Any

import pandas as pd

from .database import Database, load_database
from .errors import CpEngineError
from .inputs import NormalisedInput, SystemInput, normalise
from .report import Log
from .rules import BASE_SYMBOLS, CP_TAGS, SYMBOLS, Rule, build_rules, cp_symbol
from .system import CpBreakdown
from .validator import validate

DEFAULT_TOLERANCE = 0.02     # relative tolerance for consistency checks
MEASURED_ORDER = ("dH", "Cv", "gamma", "dS")


@dataclass
class Fact:
    value: Any
    source: str                        # "input" or the rule name
    parents: tuple[str, ...] = ()


@dataclass
class Step:
    rule: Rule
    inputs: dict[str, Any]
    outputs: dict[str, Any]


@dataclass
class RouteResult:
    tag: str
    label: str
    kind: str                          # "point" or "mean"
    cp_mass: float | None
    cp_molar: float | None
    steps: list[Step]
    breakdown: CpBreakdown | None = None

    def derivation(self) -> str:
        lines = [f"{self.label}"]
        for i, s in enumerate(self.steps, 1):
            ins = ", ".join(f"{k} = {_fmt(v)}" for k, v in s.inputs.items())
            outs = ", ".join(f"{k} = {_fmt(v)}" for k, v in s.outputs.items())
            lines.append(f"  {i}. [{s.rule.name}]  {s.rule.equation}")
            lines.append(f"       with {ins}" if ins else "       (no inputs)")
            lines.append(f"       → {outs}")
        return "\n".join(lines)


def _fmt(v: Any) -> str:
    if isinstance(v, float):
        return f"{v:.6g}"
    if hasattr(v, "classification"):
        return "<mixture model>"
    return str(v)


@dataclass
class CpResult:
    """Everything the engine found out. `print(result)` gives a readable report."""
    success: bool
    cp_mass: float | None = None             # J/(kg K)
    cp_molar: float | None = None            # J/(mol K)
    recommended: str | None = None
    routes: list[RouteResult] = field(default_factory=list)
    discrepancies: list[dict[str, Any]] = field(default_factory=list)
    breakdown: CpBreakdown | None = None
    classification: str | None = None
    behaviour: list[dict[str, str]] = field(default_factory=list)
    cv_mass: float | None = None
    warnings: list[str] = field(default_factory=list)
    assumptions: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    knowledge: dict[str, Any] = field(default_factory=dict)

    # ------------------------------------------------------------------
    def routes_table(self) -> pd.DataFrame:
        rows = []
        for r in self.routes:
            rows.append({"route": r.label, "kind": r.kind,
                         "cp [J/(kg K)]": r.cp_mass, "Cp [J/(mol K)]": r.cp_molar})
        return pd.DataFrame(rows)

    def breakdown_table(self) -> pd.DataFrame | None:
        if self.breakdown is None:
            return None
        pm, pn = self.breakdown.per_mass(), self.breakdown.per_mole()
        return pd.DataFrame({"J/(kg K)": pm, "J/(mol K)": pn}).rename(index={
            "ideal": "ideal (Σ n_i Cp_i)", "excess": "excess (Cp^E)",
            "latent": "latent / apparent", "total": "total"})

    def derivation(self, which: int | str | None = None) -> str:
        if not self.routes:
            return "No route reached Cp."
        if which is None:
            return "\n\n".join(r.derivation() for r in self.routes)
        if isinstance(which, int):
            return self.routes[which].derivation()
        return "\n\n".join(r.derivation() for r in self.routes if which in r.tag or which in r.label)

    def to_dict(self) -> dict[str, Any]:
        return {
            "success": self.success, "cp_mass_J_per_kgK": self.cp_mass, "cp_molar_J_per_molK": self.cp_molar,
            "recommended_route": self.recommended,
            "routes": [{"label": r.label, "kind": r.kind, "cp_mass": r.cp_mass, "cp_molar": r.cp_molar,
                        "rules": [s.rule.name for s in r.steps]} for r in self.routes],
            "discrepancies": self.discrepancies,
            "breakdown_mass": self.breakdown.per_mass() if self.breakdown else None,
            "breakdown_molar": self.breakdown.per_mole() if self.breakdown else None,
            "classification": self.classification, "cv_mass": self.cv_mass,
            "warnings": self.warnings, "assumptions": self.assumptions, "notes": self.notes,
            "errors": self.errors, "missing": self.missing,
        }

    def __str__(self) -> str:
        w = lambda s, ind="    ": textwrap.fill(s, 100, initial_indent=ind[:-2] + "- ", subsequent_indent=ind)
        out = ["=" * 72, "Cp ENGINE RESULT", "=" * 72]
        if self.success:
            out.append(f"Cp (mass basis)  = {self.cp_mass:,.2f} J/(kg K)" if self.cp_mass is not None
                       else "Cp (mass basis)  = n/a (needs composition or mass)")
            out.append(f"Cp (molar basis) = {self.cp_molar:,.4f} J/(mol K)" if self.cp_molar is not None
                       else "Cp (molar basis) = n/a (needs composition or moles)")
            out.append(f"Recommended route: {self.recommended}")
        else:
            out.append("Cp could NOT be computed.")
        if self.classification:
            out.append(f"System: {self.classification}")
        if self.routes:
            out.append("\nRoutes that reached Cp:")
            for r in self.routes:
                cm = f"{r.cp_mass:,.2f} J/(kg K)" if r.cp_mass is not None else "-"
                cn = f"{r.cp_molar:,.4f} J/(mol K)" if r.cp_molar is not None else "-"
                out.append(f"  * {r.label} [{r.kind}]: {cm} | {cn}")
        if self.discrepancies:
            out.append("\nComparison between routes:")
            for d in self.discrepancies:
                out.append(f"  * {d['route']} vs {d['reference']}: {d['relative_difference_%']:+.2f} %")
        if self.breakdown is not None:
            bd = self.breakdown
            pm = bd.per_mass()
            what = f"at T = {bd.T:.2f} K" if bd.kind == "point" else f"mean over [{bd.T1:.2f}, {bd.T2:.2f}] K"
            out.append(f"\nBreakdown (Route A, {what}), J/(kg K):")
            out.append(f"  ideal {pm['ideal']:,.2f} | excess {pm['excess']:,.2f} | "
                       f"latent/apparent {pm['latent']:,.2f} | total {pm['total']:,.2f}")
        if self.cv_mass is not None:
            out.append(f"Cv estimate (Cp − Tvα²/κ_T): {self.cv_mass:,.2f} J/(kg K)")
        for title, items in (("Errors", self.errors), ("Missing information", self.missing),
                             ("Warnings", self.warnings), ("Assumptions", self.assumptions), ("Notes", self.notes)):
            if items:
                out.append(f"\n{title}:")
                out += [w(i) for i in items]
        out.append("=" * 72)
        return "\n".join(out)

    def _repr_pretty_(self, p, cycle):      # nice display in Jupyter
        p.text(str(self))


# ===========================================================================
class Context:
    """State shared by the rules during one solve."""

    def __init__(self, db: Database, inp: NormalisedInput, kb: dict[str, Fact]):
        self.db, self.inp, self.kb, self.log = db, inp, kb, inp.log
        self.breakdowns: dict[str, CpBreakdown] = {}

    def is_known(self, s: str) -> bool:
        return s in self.kb

    def value(self, s: str) -> Any:
        return self.kb[s].value

    def known_values(self) -> dict[str, Any]:
        return {k: f.value for k, f in self.kb.items()}


class Solver:
    def __init__(self, db: Database, inp: NormalisedInput, rules: list[Rule] | None = None):
        self.db, self.inp = db, inp
        self.rules = rules or build_rules()
        self.kb: dict[str, Fact] = {k: Fact(v, "input") for k, v in inp.known.items()}
        if inp.amounts:
            self.kb["composition"] = Fact(dict(inp.amounts), "input")
        self.ctx = Context(db, inp, self.kb)
        self.steps: dict[str, Step] = {}           # output symbol -> step that produced it
        self.fired: set[str] = set()
        self.failures: list[tuple[str, str]] = []
        self.tol = float(inp.options.get("tolerance", DEFAULT_TOLERANCE))

    # ------------------------------------------------------------------
    def run(self) -> None:
        progress = True
        while progress:
            progress = False
            for rule in self.rules:
                if rule.name in self.fired or not all(i in self.kb for i in rule.inputs):
                    continue
                unknown = [o for o in rule.outputs if o not in self.kb]
                if not unknown and not rule.check:
                    continue
                self.fired.add(rule.name)
                vals = {i: self.kb[i].value for i in rule.inputs}
                try:
                    outs = rule.apply(vals, self.ctx)
                except CpEngineError as exc:
                    self.failures.append((rule.name, str(exc)))
                    continue
                except (ZeroDivisionError, ValueError, OverflowError) as exc:
                    self.failures.append((rule.name, f"numerical problem: {exc}"))
                    continue
                step = Step(rule, vals, dict(zip(rule.outputs, outs)))
                for sym, val in zip(rule.outputs, outs):
                    if sym in self.kb:
                        self._consistency(rule, sym, val)
                    else:
                        self.kb[sym] = Fact(val, rule.name, rule.inputs)
                        self.steps[sym] = step
                        progress = True

    def _consistency(self, rule: Rule, sym: str, new: Any) -> None:
        old = self.kb[sym].value
        if not isinstance(old, (int, float)) or not isinstance(new, (int, float)):
            return
        scale = max(abs(old), abs(new))
        if scale < 1e-12:
            return
        rel = abs(new - old) / scale
        desc, unit = SYMBOLS.get(sym, (sym, ""))
        msg = (f"Consistency check ({rule.name}: {rule.equation}): {desc} = {old:.6g} {unit} "
               f"(known) vs {new:.6g} {unit} (from the relation), difference {rel * 100:.2f} %.")
        if rel > self.tol:
            if rule.check_level == "warn":
                self.inp.log.warn("INCONSISTENT INPUT – " + msg)
            else:
                self.inp.log.note(msg + " The supplied value is used.")
        else:
            self.inp.log.note(msg + " Consistent.")

    # ------------------------------------------------------------------
    def chain(self, sym: str) -> list[Step]:
        """Steps (in order) that produced `sym` from user inputs."""
        seen, order = set(), []

        def visit(s: str) -> None:
            if s in seen or s not in self.steps:
                return
            seen.add(s)
            st = self.steps[s]
            for i in st.rule.inputs:
                visit(i)
            if st not in order:
                order.append(st)
        visit(sym)
        return order

    def leaves(self, sym: str) -> set[str]:
        out, stack = set(), [sym]
        while stack:
            s = stack.pop()
            f = self.kb.get(s)
            if f is None:
                continue
            if f.source == "input":
                out.add(s)
            else:
                stack.extend(f.parents)
        return out

    def route_label(self, tag: str) -> str:
        sym = cp_symbol(tag) if cp_symbol(tag) in self.kb else cp_symbol(tag, "molar")
        lv = self.leaves(sym)
        if tag == "dH":
            if "Q" in lv or "adiabatic" in lv:
                return "Route B – calorimetric (constant P: ΔH = Q + W, cp = ΔH/(mΔT))"
            if "dU" in lv:
                return "Route C – from ΔU + PΔV = ΔH, cp = ΔH/(mΔT)"
            return "Route C – mean Cp from ΔH/ΔT"
        if tag == "Cv":
            if "Q" in lv:
                return "Route B – calorimetric at constant V (Q = ΔU → cv → cp = cv + Tvα²/κ_T)"
            if "dU" in lv:
                return "Route C – from ΔU at constant V (cv → cp = cv + Tvα²/κ_T)"
            return "Route C – from Cv (cp = cv + Tvα²/κ_T)"
        if tag == "gamma":
            if "Q" in lv or "dU" in lv:
                return "Route B/C – cv from calorimetry, then cp = γ cv"
            return CP_TAGS["gamma"]
        return CP_TAGS.get(tag, tag)

    # ------------------------------------------------------------------
    def missing_suggestions(self, limit: int = 6) -> list[str]:
        """Smallest sets of extra inputs that would complete a route to Cp."""
        producers: dict[str, list[Rule]] = {}
        for r in self.rules:
            for o in r.outputs:
                producers.setdefault(o, []).append(r)

        def minimal(sets):
            sets = sorted(set(sets), key=len)
            out = []
            for s in sets:
                if not any(o <= s for o in out):
                    out.append(s)
            return out[:25]

        memo: dict[tuple[str, int], list[frozenset]] = {}

        def need(sym: str, depth: int, path: frozenset) -> list[frozenset]:
            if sym in self.kb:
                return [frozenset()]
            key = (sym, depth)
            if key in memo:
                return memo[key]
            alts = [frozenset({sym})] if sym in BASE_SYMBOLS else []
            if depth > 0:
                for r in producers.get(sym, []):
                    if r.name in self.fired or any(i in path for i in r.inputs):
                        continue
                    combos = [frozenset()]
                    for i in r.inputs:
                        sub = need(i, depth - 1, path | {sym})
                        if not sub:
                            combos = []
                            break
                        combos = minimal(a | b for a, b in itertools.product(combos, sub))
                    alts += combos
            memo[key] = minimal(alts)
            return memo[key]

        suggestions: list[tuple[int, str]] = []
        for r in self.rules:
            if not any(o.startswith("Cp_") for o in r.outputs):
                continue
            combos = [frozenset()]
            for i in r.inputs:
                sub = need(i, 5, frozenset())
                if not sub:
                    combos = []
                    break
                combos = minimal(a | b for a, b in itertools.product(combos, sub))
            tag = r.outputs[0][r.outputs[0].index("[") + 1:-1]
            route = {"A:point": "Route A (property route, Cp at T)", "A:mean": "Route A (mean Cp over T1–T2)",
                     "dH": "Route B/C (from Q, W or ΔH)", "Cv": "Route C (from Cv)",
                     "gamma": "Route C (from γ and Cv)", "dS": "Route C (from ΔS)"}[tag]
            for c in combos[:3]:
                if c:
                    items = " + ".join(SYMBOLS[s][0] for s in sorted(c))
                    suggestions.append((len(c), f"{route}: provide {items}"))
        suggestions = sorted(set(suggestions))
        one_step = [f"Rule '{r.name}' ({r.equation}) needs only: {SYMBOLS.get(miss[0], (miss[0],))[0]}"
                    for r in self.rules
                    if r.name not in self.fired and not any(o in self.kb for o in r.outputs)
                    and len(miss := [i for i in r.inputs if i not in self.kb]) == 1 and miss[0] in BASE_SYMBOLS]
        return [s for _, s in suggestions[:limit]] + one_step[:4]


# ===========================================================================
class CpEngine:
    """Public entry point.

    >>> engine = CpEngine()
    >>> res = engine.solve({"T": "25 C", "components": ["water", "ethanol"],
    ...                     "mass_fractions": {"water": 0.6, "ethanol": 0.4}})
    >>> print(res)
    """

    def __init__(self, data_dir: str | None = None):
        self.db = load_database(str(data_dir) if data_dir else None)

    def solve(self, inputs: dict | SystemInput, raise_errors: bool = False) -> CpResult:
        """Compute Cp from any valid combination of inputs (see cp_engine.inputs)."""
        try:
            inp = normalise(inputs, self.db)
            validate(inp, self.db)
        except CpEngineError as exc:
            if raise_errors:
                raise
            return CpResult(success=False, errors=[f"{type(exc).__name__}: {exc}"])
        solver = Solver(self.db, inp)
        solver.run()
        return self._assemble(solver, inp)

    # ------------------------------------------------------------------
    def system(self, inputs: dict | SystemInput) -> "SystemModel":
        """Mixture model (Route A) built from the composition part of `inputs`."""
        from .system import SystemModel
        inp = normalise(inputs, self.db)
        validate(inp, self.db)
        if not inp.amounts:
            raise CpEngineError("A composition is needed to build the mixture model.")
        return SystemModel(self.db, inp.amounts, inp.phase_spec, inp.known.get("P", 101325.0),
                           inp.options, inp.log)

    def cp_curve(self, inputs: dict | SystemInput, T) -> pd.DataFrame:
        """Route-A Cp(T) and its contributions on a temperature grid (K).

        Columns: T_K, ideal, excess, latent, total  (all J/(kg K)) and total_molar (J/(mol K)).
        """
        import numpy as np
        T = np.asarray(T, dtype=float)
        s = self.system(inputs)
        s.configure(float(T.min()), float(T.max()))
        rows = []
        for t in T:
            bd = s.cp_point(float(t))
            pm = bd.per_mass()
            rows.append({"T_K": t, **pm, "total_molar": bd.total / bd.moles})
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    def _assemble(self, s: Solver, inp: NormalisedInput) -> CpResult:
        kb, log = s.kb, inp.log
        res = CpResult(success=False)
        M = kb["M"].value if "M" in kb else None
        tags = []
        for sym in kb:
            if sym.startswith("Cp_"):
                tag = sym[sym.index("[") + 1:-1]
                if tag not in tags:
                    tags.append(tag)
        for tag in tags:
            cm = kb.get(cp_symbol(tag))
            cn = kb.get(cp_symbol(tag, "molar"))
            cp_mass = cm.value if cm else (cn.value / M if cn and M else None)
            cp_molar = cn.value if cn else (cm.value * M if cm and M else None)
            main = cp_symbol(tag) if cm else cp_symbol(tag, "molar")
            leaves = s.leaves(main)
            kind = "point" if tag == "A:point" or (tag in ("Cv", "gamma") and not ({"dU", "Q"} & leaves)) else "mean"
            res.routes.append(RouteResult(tag, s.route_label(tag), kind, cp_mass, cp_molar, s.chain(main),
                                          s.ctx.breakdowns.get(tag)))
        if "system" in kb:
            model = kb["system"].value
            res.classification = model.classification()
            res.behaviour = model.behaviour_table()
        if res.routes:
            res.success = True
            order = list(MEASURED_ORDER) + ["A:mean", "A:point"]
            best = sorted(res.routes, key=lambda r: order.index(r.tag))[0]
            res.cp_mass, res.cp_molar, res.recommended = best.cp_mass, best.cp_molar, best.label
            if best.tag.startswith("A"):
                res.recommended += " (model prediction)"
            else:
                res.recommended += " (measured data take precedence over the model prediction)"
            res.breakdown = s.ctx.breakdowns.get("A:mean") or s.ctx.breakdowns.get("A:point")
            # discrepancies: compare each measured route with the model (like with like)
            by_tag = {r.tag: r for r in res.routes}
            for r in res.routes:
                if r.tag.startswith("A"):
                    continue
                ref = by_tag.get("A:mean" if r.kind == "mean" else "A:point") or by_tag.get("A:point") or by_tag.get("A:mean")
                if ref is None and r is not best:
                    ref = best
                if ref is None or ref is r:
                    continue
                a, b = (r.cp_mass, ref.cp_mass) if r.cp_mass is not None and ref.cp_mass is not None else (r.cp_molar, ref.cp_molar)
                if a is None or b is None:
                    continue
                rel = (a - b) / b * 100
                res.discrepancies.append({"route": r.label, "reference": ref.label, "value": a,
                                          "reference_value": b, "relative_difference_%": rel})
                if abs(rel) > s.tol * 100:
                    log.warn(f"Routes disagree by {rel:+.2f} %: '{r.label}' vs '{ref.label}'. "
                             "Check the inputs or the model assumptions.")
            # Cv = Cp − Tvα²/κ_T holds at fixed phase amounts only, so skip it when phases change
            no_pc = "system" not in kb or not kb["system"].value.has_phase_change
            if (bd := s.ctx.breakdowns.get("A:point")) is not None and no_pc and {"alpha", "kappa_T", "rho"} <= set(kb):
                from .models.cv_cp import cv_from_cp
                try:
                    res.cv_mass = cv_from_cp(bd.total / bd.mass, bd.T, 1 / kb["rho"].value,
                                             kb["alpha"].value, kb["kappa_T"].value)
                except ZeroDivisionError:
                    pass
            if "A:point" in by_tag and "A:mean" in by_tag and "system" in kb and kb["system"].value.has_phase_change:
                log.note("With a phase change in the interval the point value (at the mid-point) and the "
                         "interval mean differ strongly; the mean includes all latent heat absorbed in [T1, T2].")
        else:
            res.missing = s.missing_suggestions()
            if not res.missing:
                res.missing = ["No combination of additional inputs found; check the composition and data."]
        for rule, msg in s.failures:
            (res.errors if not res.success else res.warnings).append(f"Rule '{rule}' could not be applied: {msg}")
        res.warnings += log.warnings
        res.assumptions, res.notes = list(log.assumptions), list(log.notes)
        res.knowledge = {k: f.value for k, f in kb.items() if not hasattr(f.value, "classification")}
        return res


def solve(inputs: dict | SystemInput, raise_errors: bool = False) -> CpResult:
    """Convenience function: CpEngine().solve(inputs)."""
    return CpEngine().solve(inputs, raise_errors)
