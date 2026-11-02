"""User-friendly front ends.

form()    – an ipywidgets form for Jupyter / Colab
guided()  – step-by-step input() prompts (also: `python -m cp_engine`)

Both only collect what the user knows and hand a plain dict to
CpEngine.solve(); any field can be left empty.
"""
from __future__ import annotations

from .solver import CpEngine

PHASES = ["auto", "liquid", "solid", "gas", "dissolved"]
BASES = {"mass fractions": "mass_fractions", "mole fractions": "mole_fractions",
         "volume fractions": "volume_fractions", "masses (kg)": "masses", "moles (mol)": "moles"}


def _parse(text: str):
    """'25 C' -> '25 C'; '' -> None; plain numbers stay numbers."""
    text = (text or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return text


def build_input(rows, basis, fields) -> dict:
    """Assemble the solver input dict from form values.

    rows   : list of (component, phase, amount-text)
    basis  : key of BASES
    fields : dict of other inputs (text) – T, T1, T2, P, Q, W, dH, dS, Cv, total_mass, process
    """
    inp: dict = {}
    comps, phases, comp_amounts = [], {}, {}
    for name, phase, amount in rows:
        if not name:
            continue
        comps.append(name)
        phases[name] = phase
        val = _parse(amount)
        if BASES[basis] in ("masses", "moles") and isinstance(val, float):
            val = (val, "kg" if BASES[basis] == "masses" else "mol")
        comp_amounts[name] = val
    if comps:
        inp["components"] = comps
        inp["phases"] = phases
        if any(v is not None for v in comp_amounts.values()):
            inp[BASES[basis]] = comp_amounts
    for k, v in fields.items():
        val = _parse(v)
        if val is not None:
            inp[k] = val
    return inp


def form(engine: CpEngine | None = None):
    """Display an ipywidgets form; returns the widget container."""
    import ipywidgets as w
    from IPython.display import display

    engine = engine or CpEngine()
    names = [""] + engine.db.names()
    rows = [(w.Dropdown(options=names, description=f"Comp. {i + 1}", layout=w.Layout(width="260px")),
             w.Dropdown(options=PHASES, value="auto", layout=w.Layout(width="110px")),
             w.Text(placeholder="fraction / amount", layout=w.Layout(width="150px"))) for i in range(5)]
    basis = w.Dropdown(options=list(BASES), value="mass fractions", description="Composition")
    keys = [("T", "e.g. 25 C"), ("T1", "e.g. 20 C"), ("T2", "e.g. 30 C"), ("P", "e.g. 1 atm"),
            ("total_mass", "e.g. 2 kg"), ("Q", "e.g. 83.7 kJ"), ("W", "e.g. 0 J"), ("dH", "e.g. 40 kJ"),
            ("dS", "e.g. 120 J/K"), ("Cv", "e.g. 4.13 kJ/(kg K)"), ("phi", "")]
    fields = {k: w.Text(description=k, placeholder=ph, layout=w.Layout(width="260px")) for k, ph in keys if k != "phi"}
    process = w.Dropdown(options=["", "constant_P", "constant_V", "adiabatic", "adiabatic_constant_P"],
                         description="Process")
    out = w.Output()
    button = w.Button(description="Compute Cp", button_style="primary")

    def on_click(_):
        out.clear_output()
        vals = {k: f.value for k, f in fields.items()}
        if process.value:
            vals["process"] = process.value
        inp = build_input([(a.value, b.value, c.value) for a, b, c in rows], basis.value, vals)
        with out:
            print("Input passed to the engine:", inp)
            print(engine.solve(inp))

    button.on_click(on_click)
    ui = w.VBox([w.HTML("<b>Cp engine – leave unknown fields empty</b>"), basis,
                 *[w.HBox(list(r)) for r in rows],
                 w.HBox([fields["T"], fields["T1"], fields["T2"]]),
                 w.HBox([fields["P"], fields["total_mass"], process]),
                 w.HBox([fields["Q"], fields["W"], fields["dH"]]),
                 w.HBox([fields["dS"], fields["Cv"]]), button, out])
    display(ui)
    return ui


def guided(engine: CpEngine | None = None) -> None:
    """Ask for the inputs one by one with input(); Enter skips a question."""
    engine = engine or CpEngine()
    print("Cp engine – press Enter to skip anything you do not know.")
    print("Known components:", ", ".join(engine.db.names()))
    rows = []
    while True:
        name = input(f"Component {len(rows) + 1} name (Enter to finish): ").strip()
        if not name:
            break
        phase = input("  phase [auto/liquid/solid/gas/dissolved] (Enter = auto): ").strip() or "auto"
        amount = input("  mass fraction (or leave empty for the last one): ").strip()
        rows.append((name, phase, amount))
    fields = {}
    for key, prompt in (("T", "Temperature T (e.g. 25 C)"), ("T1", "Initial temperature T1"),
                        ("T2", "Final temperature T2"), ("P", "Pressure (e.g. 1 atm)"),
                        ("total_mass", "Total mass (e.g. 2 kg)"), ("Q", "Heat added Q (e.g. 10 kJ)"),
                        ("W", "Shaft/stirring work W"), ("dH", "Enthalpy change ΔH"),
                        ("dS", "Entropy change ΔS"), ("Cv", "Cv (e.g. 4.1 kJ/(kg K))"),
                        ("process", "Process [constant_P/constant_V/adiabatic/adiabatic_constant_P]")):
        fields[key] = input(f"{prompt}: ")
    inp = build_input(rows, "mass fractions", fields)
    print("\nInput:", inp)
    print(engine.solve(inp))
