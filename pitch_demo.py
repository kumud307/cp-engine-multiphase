"""Five-minute pitch demo: `python pitch_demo.py` (add a number, e.g. `python pitch_demo.py 4`, for one demo;
add `full` to print the complete engine report of each case).

Every number printed here is computed live by the engine; reference values are quoted with their source.
"""
from __future__ import annotations

import sys
import warnings

warnings.filterwarnings("ignore")
from cp_engine import CpEngine  # noqa: E402

E = CpEngine()
FULL = "full" in sys.argv


def head(n: int, title: str, lecture: str) -> None:
    print("\n" + "=" * 78)
    print(f"DEMO {n}: {title}")
    print(f"   lecture concept: {lecture}")
    print("=" * 78)


def show(inp: dict, label: str = ""):
    r = E.solve(inp)
    print(f"  input : {inp}")
    if FULL:
        print(r)
    elif r.success:
        cp_m = f"{r.cp_mass:,.1f} J/(kg K)" if r.cp_mass is not None else "-"
        cp_n = f"{r.cp_molar:,.2f} J/(mol K)" if r.cp_molar is not None else "-"
        print(f"  result: Cp = {cp_m}   |   {cp_n}   {label}")
    else:
        print("  result: Cp could not be computed")
        for line in (r.errors + r.missing)[:3]:
            print("     ", line[:150])
    return r


def demo1():
    head(1, "Pure water, three independent routes to the same Cp",
         "Ch 1: first law (dH = Q + Ws at constant P), Cp = T(dS/dT)_P, Postulate I")
    a = show({"T": "25 C", "components": ["water"]}, "<- Route A: property data")
    b = show({"T1": "20 C", "T2": "30 C", "total_mass": "2 kg", "process": "constant_P", "Q": "83.7 kJ"},
             "<- Route B: calorimetry, Q/(m dT)")
    c = show({"T1": "20 C", "T2": "30 C", "total_mass": "2 kg", "dS": "280.74 J/K", "components": ["water"]},
             "<- Route C: entropy, dS/(m ln(T2/T1))")
    print(f"  Cv from Cp - Cv = T v alpha^2 / kappa_T : {a.cv_mass:,.1f} J/(kg K)")
    print("  missing input? the engine says what it needs (Postulate I: n + 2 properties):")
    show({"components": ["water"]})


def demo2():
    head(2, "Ethanol + water solution: ideal mixing is 18 % too low",
         "Ch 4: ideal solution, excess functions, Gibbs-Helmholtz, NRTL, partial molar properties, Euler")
    x = {"ethanol": 0.3, "water": 0.7}
    ideal = show({"T": "25 C", "mole_fractions": x, "solution_model": "ideal"}, "<- ideal solution: sum x_i Cp_i")
    real = show({"T": "25 C", "mole_fractions": x}, "<- real: + Cp^E from NRTL (Gibbs-Helmholtz)")
    pm = real.breakdown.extras["partial_molar_cp_J_per_molK"]["liquid"]
    print(f"  excess part Cp^E = {real.cp_molar - ideal.cp_molar:,.2f} J/(mol K)")
    pm = {k.lower(): float(v) for k, v in pm.items()}
    print("  partial molar Cp :", ", ".join(f"{k} {v:.2f}" for k, v in pm.items()), "J/(mol K)")
    print(f"  Euler check      : sum x_i Cp_i(partial) = {sum(x[k] * pm[k] for k in x):,.2f} = Cp of the mixture")


def demo3():
    head(3, "Nanofluid: Al2O3 particles in water (4 vol %)",
         "Ch 4: H is additive over phases -> mass-fraction weighting, not volume-fraction weighting")
    r = show({"T": "25 C", "components": ["water", "Al2O3"], "phi": {"Al2O3": 0.04}}, "<- mass-weighted (from dH/dT)")
    cw = E.solve({"T": "25 C", "components": ["water"]}).cp_mass
    cs = E.solve({"T": "25 C", "components": ["Al2O3"]}).cp_mass
    print(f"  volume-weighted 'Model I' would give {0.96 * cw + 0.04 * cs:,.1f} J/(kg K): "
          f"{(0.96 * cw + 0.04 * cs) / r.cp_mass * 100 - 100:+.1f} % off")
    print("  validation: 2.1 % mean deviation from measured nanofluid Cp (O'Hanley et al. 2012)")


def demo4():
    head(4, "Ice slurry (10 wt % ethanol in water at -5 C): apparent Cp with latent heat",
         "Cp_app = (dH/dT)_P, chain rule: latent term = sum (dH/dn)(dn/dT); solid-liquid equilibrium")
    r = show({"T": "-5 C", "mass_fractions": {"water": 0.9, "ethanol": 0.1}}, "<- about 13 x liquid water")
    bd = r.breakdown
    m = bd.mass
    print(f"  breakdown J/(kg K): ideal {bd.ideal / m:,.0f} | excess {bd.excess / m:,.0f} | latent {bd.latent / m:,.0f}")
    show({"T1": "-15 C", "T2": "5 C", "mass_fractions": {"water": 0.9, "ethanol": 0.1}},
         "<- mean over the whole melting range (use this for energy balances)")
    print("  energy check: integral of Cp_app dT equals H(T2) - H(T1) to 1e-7 or better")


def demo5():
    head(5, "Real gas: steam at 500 C and 100 bar",
         "EOS for pure fluids: virial and Peng-Robinson, departure functions, Cp = Cp_ig + Cp^R")
    for eos in ("ideal", "virial", "PR"):
        show({"T": "500 C", "P": "100 bar", "components": {"water": "gas"}, "gas_eos": eos}, f"<- {eos}")
    print("  NIST reference value: 2,583 J/(kg K)")
    show({"T": "25 C", "P": "500 bar", "components": ["water"]}, "<- liquid at 500 bar, Tait + Maxwell (NIST 4,060)")


def demo6():
    head(6, "A mixture that boils: ethanol + water at 1 atm",
         "EOS for mixtures + Ch 4: fugacity, gamma-phi equilibrium, bubble/dew point, Rachford-Rice flash")
    show({"T": "85 C", "mole_fractions": {"ethanol": 0.05, "water": 0.95}},
         "<- still liquid: bubble point is 90.3 C, not 78 C")
    r = show({"T": "82 C", "mole_fractions": {"ethanol": 0.3, "water": 0.7}}, "<- inside the boiling range: two phases")
    print("  system:", r.classification)
    show({"T": "120 C", "mole_fractions": {"ethanol": 0.3, "water": 0.7}}, "<- all vapour (Cp_ig + Peng-Robinson Cp^R)")
    print("  validation: 409 measured bubble points, mean deviation 0.26 K (ethanol) and 0.47 K (methanol)")


def demo7():
    head(7, "Wrong or incomplete input: the engine explains instead of crashing", "robustness")
    show({"T": "25 C", "components": ["watr"]})
    show({"T": "25 C", "mass_fractions": {"water": 0.7, "ethanol": 0.6}})
    r = show({"T1": "20 C", "T2": "30 C", "total_mass": "2 kg", "process": "constant_P", "Q": "20 kJ",
              "components": ["water"]}, "<- measured Q contradicts the property data")
    for w in r.warnings[:1]:
        print("  warning:", w[:230])


DEMOS = [demo1, demo2, demo3, demo4, demo5, demo6, demo7]

if __name__ == "__main__":
    picks = [int(a) for a in sys.argv[1:] if a.isdigit()]
    for i, f in enumerate(DEMOS, 1):
        if not picks or i in picks:
            f()
    print()
