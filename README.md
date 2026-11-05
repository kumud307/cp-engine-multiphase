# Cp Engine – general specific heat of multiphase fluids

CL 202 Chemical Engineering Thermodynamics, IIT Gandhinagar – Topic 5: *Specific Heats of Multiphase Fluids*.

Give the engine **any valid combination** of what you know: temperature(s), pressure, composition
(mole / mass / volume fractions, moles, masses, solid volume fraction φ), process type, heat Q,
work W, ΔH, ΔU, ΔS, Cv, γ. It uses thermodynamic relations to reach **one output: the specific
heat capacity Cp of the system** (mass and molar basis), for any number of components, including
real solutions, suspensions/nanofluids, emulsions, slurries with phase change, boiling mixtures
(vapour–liquid equilibrium) and real gases at pressure (virial and cubic equations of state).

Every result reports:
- the Cp on a mass basis, J/(kg K), and a molar basis, J/(mol K)
- the route(s) used and the chain of equations behind each
- a breakdown into ideal, excess (Cpᴱ, and the real-gas residual Cpᴿ) and latent/apparent contributions
- an expected deviation from measurements (± %), built from the project's own validation
- for liquid solutions: partial molar Cp of each component and the mixing functions ΔG, ΔH, ΔS
- for boiling mixtures: bubble and dew points, vapour fraction and compositions
- warnings and assumptions
- if Cp cannot be reached, which extra inputs would complete a route

## Installation

```bash
python -m venv .venv
.venv\Scripts\activate            # Windows   (Linux/macOS: source .venv/bin/activate)
pip install -e .[dev]             # numpy, scipy, pandas, matplotlib (+ pytest, sympy, ipywidgets, nbclient)
pytest                            # 147 tests
```

**Google Colab:** open the notebook directly from GitHub:
[colab.research.google.com/github/kumud307/cp-engine-multiphase/blob/main/notebooks/CP_Engine_Demo.ipynb](https://colab.research.google.com/github/kumud307/cp-engine-multiphase/blob/main/notebooks/CP_Engine_Demo.ipynb).
Its first cell clones this repository and installs the package. (Offline alternative: run
`python tools/make_colab_zip.py`, set `REPO_URL = ""` and upload `cp_engine_project.zip`.)

## Usage

```python
from cp_engine import CpEngine
engine = CpEngine()
res = engine.solve({...})          # any subset of inputs, see below
print(res)                         # full report
print(res.derivation())            # chain of rules (equations) used for each route
res.breakdown_table()              # ideal / excess / latent parts (pandas)
res.cp_mass, res.cp_molar          # numbers
```

Quantities are given as `(value, "unit")`, `"value unit"` strings or plain SI numbers.
Units such as `°C, K, F, bar, atm, kJ, cal, kJ/(kg K), J/g/K, g, kmol` are understood.
Sign convention: **Q > 0 = heat added to the system, W > 0 = work done on the system, ΔU = Q + W.**

### Route A – property route (state + composition)
```python
engine.solve({"T": "25 C", "components": ["water", "ethanol"],
              "mass_fractions": {"water": 0.6, "ethanol": None}})          # missing fraction is filled
engine.solve({"T": "25 C", "components": ["water", "Al2O3"], "phi": {"Al2O3": 0.04}})   # nanofluid
engine.solve({"T": "-10 C", "mass_fractions": {"water": 0.9, "ethanol": 0.1}})          # ice slurry
engine.solve({"T1": "295 K", "T2": "308 K", "mass_fractions": {"water": 0.8, "n-octadecane": 0.2}})  # PCM, interval mean
```

### Route B – calorimetric (first law)
```python
engine.solve({"T1": "20 C", "T2": "30 C", "total_mass": "2 kg", "process": "constant_P", "Q": "83.7 kJ"})
engine.solve({"dT": "2.5 K", "m": "0.5 kg", "process": "adiabatic_constant_P", "W": "5.2 kJ"})   # stirring
engine.solve({"T1": 298.15, "T2": 308.15, "components": ["water"], "m": 1,
              "process": "constant_V", "Q": "41.3 kJ"})                    # Q = ΔU → Cv → Cp
```

### Route C – derived properties
```python
engine.solve({"dT": 10, "n": "2 mol", "dH": "1.5 kJ"})                                   # ΔH/ΔT
engine.solve({"T1": 290, "T2": 320, "m": 1, "dS": "411.6 J/K", "components": ["water"]})  # ΔS/ln(T2/T1)
engine.solve({"T": 300, "Cv": "4.1302 kJ/(kg K)", "components": ["water"]})              # Cp = Cv + Tvα²/κ_T
engine.solve({"Cv": "20.8 J/(mol K)", "gamma": 1.4, "components": ["n-hexane"]})         # Cp = γ Cv
```

If several routes are possible, all are computed and compared, and over-specified energy inputs
(for example Q and ΔH at constant P) are checked for consistency.

### Pressure, electrolytes and your own data
```python
engine.solve({"T": "25 C", "P": "500 bar", "components": ["water"]})        # Tait pressure correction (4056 J/(kg K))
engine.solve({"T": "25 C", "mass_fractions": {"water": 0.965, "NaCl": 0.035}})   # brine: apparent molar Cp of NaCl(aq)
engine.solve({"T": "25 C", "mass_fractions": {"water": 0.5, "ethylene glycol": 0.5},
              "custom_components": {"ethylene glycol": {"MW": "62.07 g/mol", "cp": "2.43 kJ/(kg K)",
                                                        "density": "1110 kg/m3", "miscibility_class": "polar"}}})
engine.solve({"T": "25 C", "components": ["water"], "pure_cp": {"water": "4.18 kJ/(kg K)"}})   # override a Cp
```

### Options
| key | meaning |
|---|---|
| `phases` / `components={name: phase}` | `solid`, `liquid`, `gas`, `dissolved` or `auto` (default: inferred from T_m, T_b(P), miscibility, solubility) |
| `solution_model` | `auto` (NRTL where parameters exist), `ideal`, `nrtl` |
| `melt_width` | width of the smoothed melting step for pure components, e.g. `(2, "K")` |
| `sle_activity` | `auto` (default: NRTL γ when the parameter set was fitted to freezing points, so one liquid model gives Hᴱ and the ice equilibrium), `ideal` (γ = 1; 0.84 K on freezing points) or `nrtl` |
| `sle_dcp` | include the ΔCp (Kirchhoff) term in the solid–liquid equilibrium equation (default `True`) |
| `nrtl_parameters` | `HE` (default; NRTL refitted to measured excess enthalpies, and for water + ethanol also to freezing points) or `VLE` |
| `vle_parameters` | NRTL set for boiling of mixtures: `VLE` (default; bubble points 0.3–0.5 K) or `HE` |
| `vapour_model` | `ideal` (default, modified Raoult's law) or `lewis_randall` (virial fugacity coefficients) |
| `gas_eos` | residual Cp of gases: `PR` (default), `SRK`, `RK`, `vdW`, `virial` or `ideal` |
| `custom_components` | `{name: {MW, cp, [phase, density, alpha, kappa_T, Tm, dHfus, Tb, dHvap, miscibility_class]}}`, for one calculation |
| `pure_cp` | `{name: Cp}` or `{name: {phase: Cp}}`, overrides the database value for one calculation |
| `liquid_phases` | force the grouping of liquids into phases, e.g. `{"aqueous": ["water", "ethanol"], "organic": ["toluene"]}` |
| `tolerance` | relative tolerance for consistency checks (default 0.02) |

Interactive use: `from cp_engine.interactive import form; form()` in Jupyter, or `python -m cp_engine` in a terminal.

## Adding a new component (data only, no code)

Add an entry to `data/components.json` (every value needs a `source`; mark uncertain values `"verify": true`).
Template below: replace the `C1 … C5` placeholders with the published coefficients and check every number against its source.

```json
"ethylene glycol": {
  "name": "Ethylene glycol", "formula": "C2H6O2", "cas": "107-21-1",
  "aliases": ["meg", "glycol"], "miscibility_class": "polar",
  "MW":    {"value": 62.068, "unit": "g/mol", "source": "IUPAC", "verify": false},
  "Tm":    {"value": 260.2,  "unit": "K",     "source": "...", "verify": true},
  "dHfus": {"value": 9.96,   "unit": "kJ/mol","source": "...", "verify": true},
  "Tb":    {"value": 470.4,  "unit": "K",     "source": "...", "verify": true},
  "dHvap": {"value": 52.2,   "unit": "kJ/mol", "T": 470.4, "source": "...", "verify": true},
  "cp": {"liquid": {"type": "dippr100", "coeffs": [C1, C2, C3, C4, C5], "unit": "J/(kmol K)",
                    "Tmin": 260.2, "Tmax": 493.0, "source": "Perry's Table 2-153", "verify": true}},
  "density": {"liquid": {"type": "dippr105", "coeffs": [C1, C2, C3, C4], "Tmin": 260.2, "Tmax": 720.0,
                         "source": "Perry's Table 2-32", "verify": true}},
  "kappa_T": {"liquid": {"value": 3.6e-10, "unit": "1/Pa", "source": "...", "verify": true}}
}
```

Correlation types: `dippr100`, `dippr107` (ideal gas), `shomate`, `poly` (any unit), `constant`.
Density types: `constant`, `poly`, `dippr105`, `kell_water`, `ideal_gas`.
Binary data go into `data/binary_nrtl.json`: NRTL pairs (τᵢⱼ = aᵢⱼ + bᵢⱼ/T + eᵢⱼ ln T + fᵢⱼT, αᵢⱼ = c + d(T − 273.15)),
miscibility flags and solubility tables. Then run `python tools/make_data_to_verify.py`.

## Project layout

```
cp_engine/          units, database, inputs (normaliser), validator, rules, solver, system (Route A),
  models/           pure, ideal, nrtl, heterogeneous, phase_change, cv_cp,
                    eos (virial + cubic), vle (bubble/dew/flash), partial_molar (Gibbs–Duhem)
  uncertainty.py    expected deviation of each result from the validation table
  plotting.py, interactive.py
data/               components.json, binary_nrtl.json, model_uncertainty.json
validation/         case1 … case6 scripts, run_all.py, data/ (NIST, ThermoML, nanofluid and freezing-point data), results/
tests/              pytest suite (185 tests, incl. regression tests against literature data)
notebooks/          CP_Engine_Demo.ipynb (+ build_notebook.py that generates it)
docs/               TEAM_HANDBOOK.md (start here: demo script + defence Q&A), DATA_SOURCES.md, DATA_TO_VERIFY.md,
                    flowchart.md (professor format), flowchart_detailed.md, WORK_BREAKDOWN.md
manuscript/         main.tex (achemso, JCED two-column), references.bib
figures/            all figures (PNG + PDF, 300 dpi)
tools/              extract_thermoml.py, extract_thermoml_phase_eq.py, fetch_nist_fluid.py, fit_vapour_pressure.py,
                    make_uncertainty_table.py, digitize_ohanley2012.py, iapws_ice.py, fit_nrtl_he.py, fit_octadecane.py,
                    make_colab_zip.py, make_data_to_verify.py, apply_verified_data_2026_11_02.py
```

## Validation

`python validation/run_all.py` writes `validation/results/summary.{csv,md}` and LaTeX tables. All data
are published measurements; `docs/DATA_SOURCES.md` documents where each came from.

| Case | Data | Result |
|---|---|---|
| 1 | 5 pure liquids, NIST SRD 69 (128 points) + 15 literature point values | AARD 0.33 % (max 1.25 %); Cv → Cp route 0.02 % |
| 1 | 5 liquids at 0.1–100 MPa, NIST isotherms (48 points) | **0.39 %** with the Tait correction (1.85 % without) |
| 1 | Ideal-gas Cp of 5 components, NIST, 200–1000 K | AARD ≤ 0.46 % (max 1.25 %) |
| 2 | Hᴱ of water + ethanol / methanol (NIST ThermoML, 239 points) | RMSD 256 → 34 and 141 → 23 J/mol (fit); **46 and 48 J/mol held-out** |
| 2b | Cpᴱ of water + ethanol, source paper held out of the fit | RMSD 9.1 (VLE) → **3.7 J/(mol K)** |
| 2 | Cp of water + methanol (Tanaka et al. 2007, independent) | ideal 7.5 %, VLE-NRTL 6.5 %, **Hᴱ-NRTL 2.7 %** |
| 2 | Cp of NaCl(aq), independent point | **0.17 %** (solid-Cp estimate 2.9 %) |
| 2 | Cp of glycerol + water, 162 points | 0.52 % (ideal fallback) |
| 2 | Cp of ethanol + toluene (Vega-Maza et al. 2011), NRTL fitted to VLE + Cpᴱ | **1.3 %** (ideal mixing 8.3 %); VLE bubble T 0.6 K |
| 2b | Cpᴱ of ethanol + toluene, half the compositions held out | **3.0 J/(mol K)**; from VLE data alone 7.1 |
| 2 | γ∞ in water (independent, Straver 2005; Atik 2004) | ethanol 35 % (Cp-optimised set) / 12 % (VLE set); acetone 25 % |
| 3 | 69 nanofluid DSC points, Al₂O₃ / CuO / SiO₂ (O'Hanley et al. 2012) | **Model II = Route A 2.1 %**, Model I 18.5 % |
| 4 | Freezing points of aqueous ethanol ≤ 30 wt% (Lange 1967) | consistent NRTL 1.7 K (fit), **2.5 K held-out**; ideal + ΔCp 0.84 K |
| 5 | Route A → Q → Route B, 7 systems | agree to < 10⁻⁶ % |
| 6 | Bubble points of water + ethanol / methanol, 409 points from 24 papers (ThermoML) | **0.26 K / 0.46 K** |
| 6 | Vapour Cp at 0.5–10 MPa, 4 fluids (NIST SRD 69, 189 points) | ideal gas 12.1 % → **Peng–Robinson residual 5.3 %** |
| 6 | Vapour pressure from vdW / RK / SRK / PR, 5 fluids | 224 / 58 / 2.3 / 1.4 % |
| – | CO₂ and CO₂ + CH₄ worked examples of the CL 202 lectures | reproduced to within ±0.01 of the printed values (rounding) (tests) |

The only dataset still missing is a direct Cp measurement of water + ethanol mixtures (behind a paywall).
The template `validation/data/water_ethanol_cp_exp.csv` is read automatically once filled in.

## Main assumptions and limitations

- Ideal volume of mixing (Vᴱ = 0) for densities, α and κ_T of mixtures; interfacial enthalpy of particles is neglected.
- Cpᴱ only for pairs with NRTL parameters. Pairs without them are treated as ideal, with a warning. For water + ethanol and water + methanol the engine uses parameters refitted to measured Hᴱ (`tools/fit_nrtl_he.py`); VLE-fitted sets underpredict Cpᴱ six-fold or more. Ethanol + toluene is fitted jointly to VLE and Cpᴱ. The Cp-optimised water + ethanol set describes γ∞ less well than the VLE set (35 % vs 12 %); adding γ∞ to the fit worsened held-out Cpᴱ, so the trade-off is reported rather than hidden.
- Boiling of a miscible liquid mixture is computed (γ–φ VLE with NRTL γ and Antoine P_sat; ideal vapour by default). Three-phase boiling of immiscible liquids and solute partitioning between liquid phases are not modelled; the engine warns. Miscibility follows a measured coexistence curve where one is stored (methanol + n-hexane, UCST 314.9 K), otherwise a fixed flag.
- Every NRTL pair in use is checked for a spurious liquid–liquid split at the working temperature (the literature acetone + water set fails above ≈ 360 K; the engine warns).
- Dissolved NaCl uses an apparent molar heat capacity fitted at 298 K (warning elsewhere); other dissolved solids use the Neumann–Kopp approximation (solid-state Cp).
- Pressure: liquid and solid Cp are corrected from 1 atm with the Tait equation (validated to 100 MPa for 5 liquids; other components use a constant κ_T). Gases: ideal-gas Cp + residual Cpᴿ from the Peng–Robinson equation (5 % average error vs NIST at 0.5–10 MPa, larger within ~20 K of saturation).
- Pure-component melting is smoothed over a finite, user-adjustable width. Point values inside the step depend on it (the engine warns); the latent heat and any interval-mean Cp spanning the transition do not.
