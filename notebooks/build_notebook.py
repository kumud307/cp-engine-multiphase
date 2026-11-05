"""Generate CP_Engine_Demo.ipynb (run: python notebooks/build_notebook.py).

Keeping the notebook as code makes it easy to review and regenerate.
"""
from pathlib import Path

import nbformat as nbf

cells: list[tuple[str, str]] = []
md = lambda s: cells.append(("md", s.strip("\n")))
code = lambda s: cells.append(("code", s.strip("\n")))

# ---------------------------------------------------------------------------
md(r"""
# General Specific-Heat (Cp) Engine for Multiphase Fluids
**CL 202 – Chemical Engineering Thermodynamics, IIT Gandhinagar — Topic 5: Specific Heats of Multiphase Fluids**

The engine accepts *any valid combination* of thermodynamic information — temperature, pressure,
composition (mole / mass / volume fractions, amounts, solid volume fraction φ), heat and work,
ΔH, ΔU, ΔS, Cv, γ … — and uses thermodynamic relations (encoded as *rules*) to reach a single output:
the **specific heat capacity Cp of the (multiphase) system**, for any number of components.

| Route | Information used | Key relations |
|---|---|---|
| **A – property route** | T, P, composition | Σ xᵢCpᵢ, Cpᴱ = ∂Hᴱ/∂T (NRTL), mass-weighting over phases, latent terms |
| **B – calorimetric route** | Q, W, ΔT, mass | first law: ΔH = Q + W (const P), ΔU = Q + W (const V) |
| **C – derived routes** | ΔH, ΔU, ΔS, Cv, γ | Cp = ΔH/(mΔT), Cp = ΔS/ln(T₂/T₁), Cp = Cv + Tvα²/κ_T |

Sections: 1 Setup · 2 Database · 3 Route A · 4 Route B · 5 Route C · 6 Several routes at once ·
7 Multicomponent (≥ 4 components) · 8 System types · 9 Figures · 10 Error handling · 11 Numerics of Cpᴱ ·
12 Interactive form · 13 Validation
""")

md("## 1. Setup\nIn **Colab**: either set `REPO_URL` to your group's Git repository, or leave it empty and upload "
   "`cp_engine_project.zip` (created by `python tools/make_colab_zip.py`) when asked. Locally the cell just finds the project folder.")
code(r"""
import sys, os
from pathlib import Path

IN_COLAB = "google.colab" in sys.modules
REPO_URL = "https://github.com/kumud307/cp-engine-multiphase.git"   # public repo; leave "" to upload cp_engine_project.zip instead

if IN_COLAB:
    if REPO_URL:
        !git clone -q {REPO_URL} /content/cp_project
    else:
        from google.colab import files
        uploaded = files.upload()                       # choose cp_engine_project.zip
        !unzip -q -o cp_engine_project.zip -d /content/cp_project
    %cd /content/cp_project
    !pip install -q -e . sympy ipywidgets
    ROOT = Path("/content/cp_project")
else:
    ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "validation"))
os.chdir(ROOT)

import numpy as np, pandas as pd, matplotlib.pyplot as plt
from cp_engine import CpEngine, describe_inputs
from cp_engine import plotting
pd.set_option("display.width", 160); pd.set_option("display.max_colwidth", 90)
engine = CpEngine()
db = engine.db
print("Project root:", ROOT)
print("Components in database:", ", ".join(db.names()))
""")
code("print(describe_inputs())   # every input the engine understands")

# ---------------------------------------------------------------------------
md("## 2. Component database\nEvery value in `data/components.json` carries a source; values flagged `verify` are listed in `docs/DATA_TO_VERIFY.md`. Adding a component means adding a JSON entry — no code changes.")
code(r"""
rows = []
for k, c in db.components.items():
    T = 298.15
    phase = "liquid" if c.has_phase("liquid") else "solid"
    try:
        cp = c.cp_mass(phase, T)
        rho = c.density(phase, T)
    except Exception as e:
        cp, rho = np.nan, np.nan
    rows.append({"key": k, "name": c.name, "M / g mol⁻¹": c.MW * 1000, "phases with Cp data": ", ".join(c.cp_corr),
                 "Tm / K": c.Tm, "Tb / K": c.Tb, f"cp({phase[0]}, 298 K) / J kg⁻¹K⁻¹": round(cp, 1),
                 "ρ(298 K) / kg m⁻³": round(rho, 1)})
pd.DataFrame(rows)
""")
code(r"""
c = db.get("water")
print("Water liquid Cp correlation:", c.cp_corr["liquid"].type, c.cp_corr["liquid"].coeffs)
print("Source:", c.cp_corr["liquid"].source)
print("\nNumber of database values still flagged 'verify':", len(db.verify_items()))
""")

# ---------------------------------------------------------------------------
md("## 3. Route A — property route (state + composition only)\nWater + ethanol, 40 wt% ethanol at 25 °C. The ethanol fraction is deliberately left out (`None`) — the normaliser fills it from Σw = 1.")
code(r"""
res = engine.solve({"T": "25 C", "P": "1 atm",
                    "components": ["water", "ethanol"],
                    "mass_fractions": {"water": 0.6, "ethanol": None}})
print(res)
""")
code("print(res.derivation())      # the chain of rules used — useful for the oral defence")
code("res.breakdown_table()")
code("pd.DataFrame(res.behaviour)  # how each component was classified")

# ---------------------------------------------------------------------------
md("## 4. Route B — calorimetric route (heat and work)\nSign convention: Q > 0 heat **added to** the system; W > 0 work **done on** the system; ΔU = Q + W.")
code(r"""
# (a) constant pressure, heater supplies Q: ΔH = Q + W  ->  cp = ΔH/(m ΔT)
r = engine.solve({"T1": "20 C", "T2": "30 C", "total_mass": "2 kg", "process": "constant_P", "Q": "83.7 kJ"})
print(r.cp_mass, "J/(kg K) via", r.recommended); print(r.derivation(0))
""")
code(r"""
# (b) adiabatic stirring at constant pressure: Q = 0, ΔH = W
r = engine.solve({"dT": "2.5 K", "m": "0.5 kg", "process": "adiabatic_constant_P", "W": "5.2 kJ"})
print(r.cp_mass, "J/(kg K)"); print(r.derivation(0))
""")
code(r"""
# (c) constant volume (bomb-type cell): Q = ΔU -> cv -> cp = cv + T v α²/κ_T (needs composition for α, κ_T, ρ)
r = engine.solve({"T1": "25 C", "T2": "35 C", "components": ["water"], "total_mass": "1 kg",
                  "process": "constant_V", "Q": "41.3 kJ"})
print(r.routes_table().to_string()); print(); print(r.derivation("constant V"))
""")
code(r"""
# (d) energy per kg and no mass given: the engine chooses a 1 kg basis (Cp is intensive)
print(engine.solve({"T1": 300, "T2": 310, "process": "constant_P", "Q": "41.8 kJ/kg"}))
""")

# ---------------------------------------------------------------------------
md("## 5. Route C — derived-property routes")
code(r"""
# from ΔH over ΔT (molar basis only: no composition, so no molar mass is known)
print(engine.solve({"dT": 10, "n": "2 mol", "dH": "1.5 kJ"}))
""")
code(r"""
# from ΔS over [T1, T2] (constant P, constant Cp):  cp = ΔS / (m ln(T2/T1))
r = engine.solve({"T1": "290 K", "T2": "320 K", "m": "1 kg", "dS": "411.6 J/K", "components": ["water"]})
print(r.routes_table().to_string())
""")
code(r"""
# from Cv (NIST value for water at 300 K: cv = 4.1302 kJ/(kg K); NIST cp = 4.1806)
r = engine.solve({"T": "300 K", "Cv": "4.1302 kJ/(kg K)", "components": ["water"]})
print(r.routes_table().to_string()); print(); print(r.derivation("Cv"))
""")
code(r"""
# from γ = Cp/Cv and Cv
print(engine.solve({"Cv": "20.8 J/(mol K)", "gamma": 1.4, "components": ["n-hexane"]}).cp_molar, "J/(mol K)")
""")

# ---------------------------------------------------------------------------
md("## 6. Several routes at once: measured vs predicted\nWhen the inputs allow more than one route, all are computed and compared. Over-specified energy inputs are checked for consistency.")
code(r"""
r = engine.solve({"T1": "20 C", "T2": "30 C", "components": ["water"], "total_mass": "1 kg",
                  "process": "constant_P", "Q": "41.85 kJ", "dH": "41.85 kJ"})
print(r)
""")
code(r"""
# inconsistent over-specification: Q and ΔH disagree at constant P
r = engine.solve({"dT": 10, "m": 1, "process": "constant_P", "Q": "41.8 kJ", "dH": "50 kJ"})
print("\n".join(r.warnings))
""")

# ---------------------------------------------------------------------------
md("## 7. Multicomponent example (≥ 4 components)\nWater + methanol + ethanol + glycerol (one real liquid solution) with suspended Al₂O₃ particles.")
code(r"""
r = engine.solve({"T": "30 C",
                  "mass_fractions": {"water": 0.45, "methanol": 0.15, "ethanol": 0.15, "glycerol": 0.15, "Al2O3": 0.10}})
print(r)
""")

# ---------------------------------------------------------------------------
md("## 8. Each system type")
md("### 8.1 Ideal mixture vs real liquid solution")
code(r"""
base = {"T": "25 C", "mole_fractions": {"water": 0.7, "ethanol": 0.3}}
ideal = engine.solve(dict(base, solution_model="ideal"))
vle = engine.solve(dict(base, nrtl_parameters="VLE"))
real = engine.solve(base)                     # default: NRTL refitted to measured excess enthalpies
pd.DataFrame({"ideal mixture": ideal.breakdown.per_mole(), "NRTL, VLE parameters": vle.breakdown.per_mole(),
              "NRTL, H^E-fitted (default)": real.breakdown.per_mole()})
""")
md("VLE-regressed NRTL parameters reproduce Gᴱ but not its temperature derivatives, so they give a Cpᴱ about "
   "six times too small. The default parameter set was refitted to measured Hᴱ (NIST ThermoML data, "
   "`tools/fit_nrtl_he.py`); validation case 2 shows it predicts independent Cp data of water + methanol to 2.7 %.")
code(r"""
print({k: v for k, v in db.nrtl_pair("ethanol", "water").items() if k in ("fitted_to", "rmse_HE_J_per_mol", "rmse_HE_VLE_parameters", "n_points", "Tmin", "Tmax")})
""")
md("### 8.2 Suspension / nanofluid: Model I vs Model II\nRoute A weights phases by **mass** (thermal-equilibrium Model II). Model I (volume-weighted cp) is shown for comparison only — it is thermodynamically inconsistent because cp is per unit mass.")
code(r"""
r = engine.solve({"T": "25 C", "components": ["water", "Al2O3"], "phi": {"Al2O3": 0.04}})
print(r.classification)
print(f"Route A (= Model II): {r.cp_mass:.1f} J/(kg K)")
print(f"Model I (inconsistent): {r.breakdown.extras['model_I_cp_mass']:.1f} J/(kg K)")
print("phase volume fractions:", r.breakdown.extras["phase_volume_fractions"])
""")
md("### 8.3 Emulsion (immiscible liquids)")
code(r"""
r = engine.solve({"T": "25 C", "mass_fractions": {"water": 0.5, "toluene": 0.3, "ethanol": 0.2}})
print(r); pd.DataFrame(r.behaviour)
""")
md("Miscibility can depend on temperature. For methanol + n-hexane the engine reads the measured liquid–liquid "
   "coexistence curve (Trejo et al. 2006, UCST 314.9 K): two liquid phases at 25 °C, one solution at 50 °C. "
   "Ethanol + toluene is fully miscible; its NRTL parameters were fitted to measured VLE and Cpᴱ "
   "(1.3 % against Vega-Maza et al. 2011, against 8.3 % for ideal mixing).")
code(r"""
for T in ("25 C", "50 C"):
    print(T, "->", engine.solve({"T": T, "mole_fractions": {"methanol": 0.5, "n-hexane": 0.5}}).classification)
for model in ("ideal", "nrtl"):
    r = engine.solve({"T": "40 C", "mole_fractions": {"ethanol": 0.3, "toluene": 0.7}, "solution_model": model})
    print(f"ethanol + toluene, {model:5s}: Cp = {r.cp_molar:.1f} J/(mol K)")
""")
md("### 8.4 Slurries with phase change\n**Ice slurry** in 10 wt% aqueous ethanol: the ice fraction follows the solid–liquid equilibrium (freezing-point depression), so the apparent Cp includes the latent heat of the ice that melts as T rises.")
code(r"""
r = engine.solve({"T": "-10 C", "mass_fractions": {"water": 0.9, "ethanol": 0.1}})
print(r); r.breakdown.phase_amounts
""")
code(r"""
# interval mean: all latent heat absorbed between -15 °C and +5 °C
r = engine.solve({"T1": "-15 C", "T2": "5 C", "mass_fractions": {"water": 0.9, "ethanol": 0.1}, "total_mass": "1 kg"})
print(r.routes_table().to_string())
""")
md("**PCM slurry**: 20 wt% n-octadecane dispersed in water, melting at 301.3 K, smoothed over a user-adjustable width.")
code(r"""
for width in (1.0, 2.0, 4.0):
    r = engine.solve({"T": "301.3 K", "mass_fractions": {"water": 0.8, "n-octadecane": 0.2}, "melt_width": (width, "K")})
    print(f"melt width {width} K: apparent cp at Tm = {r.cp_mass/1000:.1f} kJ/(kg K)")
r = engine.solve({"T1": "295 K", "T2": "308 K", "mass_fractions": {"water": 0.8, "n-octadecane": 0.2}, "total_mass": 1})
print(r.routes_table().to_string())
""")
md("**Saturated salt slurry**: NaCl crystals in their saturated solution — `cp_app = Σ wᵢcpᵢ + Δh_sol · d(w_dissolved)/dT` (tabulated solubility).")
code(r"""
r = engine.solve({"T": "25 C", "mass_fractions": {"water": 0.7, "NaCl": 0.3}})
print(r)
""")

md("### 8.5 Pressure, electrolytes and your own data\n"
   "Cp correlations are for 1 atm. At other pressures the engine applies the exact relation "
   "(∂Cp/∂P)_T = −T(∂²v/∂T²)_P with v(T, P) from the Tait equation (validated against NIST to 100 MPa).")
code(r"""
for P in ("1 atm", "500 bar", "1000 bar"):
    r = engine.solve({"T": "25 C", "P": P, "components": ["water"]})
    print(f"water at {P:>8}: cp = {r.cp_mass:7.1f} J/(kg K)")
print("NIST SRD 69 (validation/data): 4060.1 J/(kg K) at 50.1 MPa, 3976.5 J/(kg K) at 100.1 MPa")
""")
md("**Brine:** dissolved NaCl is described by its apparent molar heat capacity (fitted to measured NaCl(aq) data), "
   "not by the heat capacity of solid salt.")
code(r"""
r = engine.solve({"T": "25 C", "mass_fractions": {"water": 0.965, "NaCl": 0.035}})
print(f"3.5 wt% NaCl: cp = {r.cp_mass:.1f} J/(kg K)"); r.breakdown_table()
""")
md("**Your own components and Cp values** (for one calculation; the database files are not changed):")
code(r"""
r = engine.solve({"T": "25 C", "mass_fractions": {"water": 0.5, "ethylene glycol": 0.5},
                  "custom_components": {"ethylene glycol": {"MW": "62.07 g/mol", "cp": "2.43 kJ/(kg K)",
                                                            "density": "1110 kg/m3", "miscibility_class": "polar"}}})
print(f"water + ethylene glycol (user-defined): {r.cp_mass:.1f} J/(kg K)")
r = engine.solve({"T": "25 C", "components": ["water"], "pure_cp": {"water": "4.0 kJ/(kg K)"}})
print(f"water with an overridden Cp: {r.cp_mass:.1f} J/(kg K)")
""")
md("**Point value inside a melting step:** the engine warns that it depends on the assumed melt width and reports "
   "width-independent quantities instead.")
code(r"""
r = engine.solve({"T": "301.3 K", "mass_fractions": {"water": 0.8, "n-octadecane": 0.2}})
print([w for w in r.warnings if "depends on the assumed width" in w][0])
""")

md("### 8.6 A boiling mixture (vapour–liquid equilibrium)\n"
   "A liquid mixture boils between its **bubble** and **dew** points, not at the boiling points of its components. "
   "Equal fugacities y_i φ̂_i P = x_i γ_i φ_i^sat P_i^sat (γ from NRTL, P^sat from Antoine equations fitted to NIST data, "
   "φ from the same Peng–Robinson equation as the vapour enthalpy); "
   "between the two points an isothermal flash gives the vapour fraction and both compositions, and the latent heat "
   "of what evaporates appears in the apparent Cp.")
code(r"""
for T, x in (("85 C", 0.05), ("84 C", 0.30), ("90 C", 0.50)):
    r = engine.solve({"T": T, "mole_fractions": {"ethanol": x, "water": 1 - x}})
    ex = r.breakdown.extras
    vle = list(ex.get("vle", {}).values())
    tb = f"bubble {vle[0]['T_bubble_K'] - 273.15:.1f} °C" if vle else ""
    print(f"x_ethanol = {x:.2f} at {T}: {r.classification.split(';')[0]:45s} {tb}; "
          f"vapour fraction {ex.get('vapour_mole_fraction_of_system', 0):.2f}; cp = {r.cp_mass:,.0f} J/(kg K)")
""")
code(r"""
from cp_engine.models.vle import VLE
from cp_engine.models.nrtl import NRTL
v = VLE(["ethanol", "water"], engine.db, NRTL(["ethanol", "water"], engine.db, prefer="VLE"))
for x in (0.85, 0.90, 0.95):
    print(f"x = {x}: bubble {v.bubble_T(101325, [x, 1 - x]):.2f} K, dew {v.dew_T(101325, [x, 1 - x]):.2f} K")
print("-> bubble = dew near x = 0.90, 351.3 K: the ethanol–water azeotrope")
""")

md("### 8.7 Equations of state from the CL 202 lectures\n"
   "`cp_engine/models/eos.py` implements the virial equation (Pitzer/Abbott B, Prausnitz combining rules) and the "
   "generic cubic (vdW, RK, SRK, PR). The lecture's worked CO₂ example comes out to within ±0.01 of the printed values (rounding):")
code(r"""
from cp_engine.models.eos import Critical, CubicEOS, VirialEOS
co2 = Critical("CO2", 304.13, 73.77e5, 0.224, 94.07e-6)          # constants on the lecture slide
print(f"Abbott B of CO2 at 350 K: {VirialEOS([co2]).B_mix(350, [1])[0] * 1e6:.2f} cm3/mol   (slide: -83.73)")
slide = {"vdW": (52.55, 81.13), "RK": (45.62, 61.03), "SRK": (41.99, 58.42), "PR": (41.60, 51.68)}
for kind, (ps, vl) in slide.items():
    e = CubicEOS([co2], kind); P = e.psat(280.0)
    print(f"{kind:3s}: Psat(280 K) = {P / 1e5:6.2f} bar (slide {ps}),  V_liq = {e.molar_volume(280, P, [1], 'liquid') * 1e6:6.2f} cm3/mol (slide {vl})")
""")
md("The same equations give the **residual heat capacity of real gases**, Cp = Cp^ig + Cp^R. Superheated steam at 100 bar:")
code(r"""
for eos in ("ideal", "virial", "PR"):
    r = engine.solve({"T": "500 C", "P": "100 bar", "components": {"water": "gas"}, "gas_eos": eos})
    print(f"{eos:6s}: cp = {r.cp_mass:6.0f} J/(kg K)")
print("NIST SRD 69: 2583 J/(kg K)")
""")

md("### 8.8 Partial molar properties, Gibbs–Duhem and the uncertainty of a result\n"
   "For a liquid solution the engine reports the partial molar Cp of each component (Σ x_i C̄p,i = Cp) and the "
   "mixing functions. Every result also carries an expected deviation from measurements, built from the validation.")
code(r"""
r = engine.solve({"T": "25 C", "mole_fractions": {"ethanol": 0.3, "water": 0.7}})
pc = r.breakdown.extras["partial_molar_cp_J_per_molK"]["liquid"]
print("partial molar Cp:", {k: round(float(v), 2) for k, v in pc.items()})
print(f"Σ x_i Cp_i = {0.3 * pc['Ethanol'] + 0.7 * pc['Water']:.2f}  vs  Cp = {r.cp_molar:.2f} J/(mol K)")
print({k: round(float(v), 1) for k, v in r.breakdown.extras["mixing_functions_J_per_mol"]["liquid"].items()})
from cp_engine.models import partial_molar as pmod
m = NRTL(["ethanol", "water"], engine.db)
print("Gibbs–Duhem residual Σ x_i ∂lnγ_i/∂x_1 =", pmod.gibbs_duhem_residual(lambda z: np.asarray(m.ln_gamma(list(z), 320.0)), [0.3, 0.7], 0))
u = r.uncertainty
print(f"\nexpected deviation: ± {u['relative_percent']:.1f} %")
for c in u["contributions"]:
    print(f"  {c['source']:45s} {c['percent_of_cp']:.2f} %   <- {c['validation'][:70]}")
""")

# ---------------------------------------------------------------------------
md("## 9. Figures (saved to `figures/` as PNG and PDF, 300 dpi)")
code(r"""
curve = engine.cp_curve({"mass_fractions": {"water": 0.9, "ethanol": 0.1}}, np.linspace(250, 280, 601))
plotting.plot_cp_vs_T(curve, "Ice slurry, 10 wt% ethanol(aq)", name="nb_cp_vs_T_ice_slurry"); plt.show()
""")
code(r"""
curve = engine.cp_curve({"mole_fractions": {"water": 0.7, "ethanol": 0.3}}, np.linspace(280, 340, 61))
plotting.plot_cp_vs_T(curve, "Water + ethanol, x(EtOH) = 0.3", name="nb_cp_vs_T_water_ethanol"); plt.show()
""")
code(r"""
plotting.plot_cp_vs_composition(db, "ethanol", "water", 298.15, engine, name="nb_cp_vs_x"); plt.show()
plotting.plot_cpE_vs_x(db, "ethanol", "water", [278.15, 298.15, 318.15], name="nb_cpE_vs_x"); plt.show()
""")
code(r"""
plotting.plot_nanofluid_models(db, "Al2O3", "water", 298.15, 0.10, name="nb_nanofluid"); plt.show()
plotting.plot_nanofluid_models(db, "CuO", "water", 298.15, 0.10, name="nb_nanofluid_cuo"); plt.show()
""")
code(r"""
curves = {f"melt width {w:g} K": engine.cp_curve({"mass_fractions": {"water": 0.8, "n-octadecane": 0.2},
                                                   "melt_width": w}, np.linspace(292, 310, 721)) for w in (1, 2, 4)}
plotting.plot_apparent_peak(curves, "20 wt% n-octadecane in water", name="nb_pcm_peak"); plt.show()
""")

# ---------------------------------------------------------------------------
md("## 10. Error handling: bad or missing input\nEvery problem is reported as a clear message — never a raw traceback.")
code(r"""
bad_inputs = {
    "fractions do not sum to 1": {"T": 300, "mass_fractions": {"water": 0.6, "ethanol": 0.6}},
    "negative fraction": {"T": 300, "mole_fractions": {"water": 1.2, "ethanol": -0.2}},
    "typo in an input name": {"temprature": 300, "components": ["water"]},
    "unknown component": {"T": 300, "mass_fractions": {"water": 0.5, "kerosene": 0.5}},
    "solid water at 350 K": {"T": "350 K", "components": {"water": "solid"}},
    "energy given in kg": {"T": 300, "components": ["water"], "Q": "5 kg"},
    "adiabatic but Q ≠ 0": {"process": "adiabatic", "Q": 100, "dT": 1, "m": 1},
    "far outside correlation range": {"T": "600 K", "components": {"glycerol": "liquid"}},
}
for title, inp in bad_inputs.items():
    r = engine.solve(inp)
    print(f"--- {title}:\n    " + " | ".join(r.errors)[:300])
""")
code(r"""
# under-specified input: the engine lists what would complete a route
print(engine.solve({"Q": "10 kJ", "process": "constant_P", "dT": 5}))
""")
code(r"""
print(engine.solve({"T": "298 K"}))
""")

# ---------------------------------------------------------------------------
md("## 11. Numerics of Cpᴱ: step-size convergence and symbolic cross-check\n"
   "Hᴱ = −RT² ∂(Gᴱ/RT)/∂T and Cpᴱ = ∂Hᴱ/∂T are evaluated by central differences. Truncation error ∝ h², round-off ∝ ε/h², so the optimum is h ≈ ε^¼·T ≈ 10⁻⁴T ≈ 0.03 K.")
code(r"""
from cp_engine.models.nrtl import NRTL
m = NRTL(["ethanol", "water"], db)
tab = pd.DataFrame(m.derivative_convergence([0.3, 0.7], 298.15), columns=["h / K", "H^E / J mol⁻¹", "Cp^E / J mol⁻¹K⁻¹"])
tab
""")
code(r"""
import sympy as sp
from cp_engine.database import R
T, x1 = sp.symbols("T x1", positive=True); x2 = 1 - x1
a, b = m.a, m.b
t12, t21 = a[0, 1] + b[0, 1] / T, a[1, 0] + b[1, 0] / T
G12, G21 = sp.exp(-0.3 * t12), sp.exp(-0.3 * t21)
g = x1 * x2 * (t21 * G21 / (x1 + x2 * G21) + t12 * G12 / (x2 + x1 * G12))
CpE = sp.diff(-R * T**2 * sp.diff(g, T), T)
exact = float(CpE.subs({x1: 0.3, T: 298.15}))
print(f"symbolic Cp^E = {exact:.8f}   numerical (h = 0.03 K) = {m.cpE([0.3, 0.7], 298.15):.8f}  J/(mol K)")
""")

# ---------------------------------------------------------------------------
md("## 12. Interactive form\nFill in only what you know and press **Compute Cp**. (Outside notebooks: `python -m cp_engine` asks the same questions with `input()` prompts.)")
code(r"""
from cp_engine.interactive import form
try:
    form(engine)
except Exception as e:
    print("ipywidgets not available:", e)
""")

# ---------------------------------------------------------------------------
md("## 13. Validation against literature data\n"
   "All comparisons use published measurements (provenance in `docs/DATA_SOURCES.md`): NIST SRD 69 (case 1), "
   "the NIST ThermoML archive (case 2: Hᴱ and Cp of water + alcohols, brine, glycerol, ethanol + toluene; case 2b: held-out cross-validation of every fit), O'Hanley et al. 2012 "
   "nanofluid DSC data (case 3), Lange's freezing points of aqueous ethanol (case 4), the Route A ↔ Route B "
   "round trip (case 5), and boiling mixtures, real gases and equations of state (case 6: 409 bubble points, NIST "
   "vapour isobars and saturation curves). Metric column: AARD % for relative deviations, RMSD for excess properties, "
   "|ΔT| in K for freezing and bubble points.")
code(r"""
import run_all
summary = run_all.main(save=True)
summary[["case", "system", "quantity", "metric", "N", "mean", "max", "status"]]
""")
code(r"""
from IPython.display import Image, display
for f in ["val1_toluene", "fig_pressure", "fig_HE_water_ethanol", "fig_cpE_from_HE_water_ethanol", "fig_cp_methanol_water",
          "fig_nanofluid_parity", "fig_nanofluid_models", "fig_freezing_points", "fig_cpT_ice_ethanol_10wt", "fig_pcm_peak",
          "fig_vle_water_ethanol", "fig_real_gas_steam"]:
    display(Image(filename=str(ROOT / "figures" / f"{f}.png"), width=460))
""")

# ---------------------------------------------------------------------------
nb = nbf.v4.new_notebook()
nb["cells"] = [nbf.v4.new_markdown_cell(s) if t == "md" else nbf.v4.new_code_cell(s) for t, s in cells]
nb["metadata"] = {"kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"},
                  "language_info": {"name": "python"}, "colab": {"provenance": []}}
out = Path(__file__).resolve().parent / "CP_Engine_Demo.ipynb"
nbf.write(nb, out)
print("written", out)
