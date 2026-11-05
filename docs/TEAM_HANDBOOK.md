# Team handbook: Cp engine for multiphase fluids

**CL 202 Chemical Engineering Thermodynamics, IIT Gandhinagar. Topic 5: Specific Heats of Multiphase Fluids.**
Final demo, manuscript and oral defence: **Friday 20 November 2026.**

Read this once from top to bottom, then reread the section on your own module (see
`docs/WORK_BREAKDOWN.md`). Every member must be able to explain **their module in depth** and
**the whole pipeline in outline**. Section 9 lists the questions we expect, with answers.

---

## 1. What the project does (say this in the first 30 seconds of the demo)

> The professor asked for a *general* Cp engine: the user gives **any valid combination** of what they
> know (temperature, pressure, composition in any form, heat, work, ΔH, ΔU, ΔS, Cv, γ) and the program
> uses thermodynamic relations to reach **one output: the specific heat capacity of the (multiphase)
> system**, for any number of components.

How we did it:

* Every thermodynamic relation is a small **rule** with inputs, outputs and an equation, e.g.
  `ΔH = Q + W` (constant P) or `cp = ΔH/(mΔT)`.
* A **forward-chaining solver** keeps applying any rule whose inputs are known, until nothing new can
  be derived. Each Cp it reaches is a *route*; if there are several, they are compared.
* If no route reaches Cp, the solver searches backwards and tells the user **which extra inputs would
  complete a route**.
* The **property route (Route A)** models the mixture itself. It classifies the system (ideal mixture,
  real solution, suspension, emulsion, slurry with phase change) and computes
  Cp = ideal part + excess part + latent part.

Three families of routes:

| Route | Needs | Main relations |
|---|---|---|
| A: property | T (or T1, T2), composition (and P) | Σ xᵢCpᵢ; Cpᴱ from NRTL via Gibbs–Helmholtz; mass-weighting over phases; latent terms from phase equilibrium (solid–liquid and vapour–liquid); real-gas Cpᴿ from an equation of state |
| B: calorimetric | Q and/or W, ΔT, amount | First law: ΔH = Q + W (const P); ΔU = Q + W (const V) |
| C: derived | ΔH, ΔU, ΔS, Cv or γ | cp = ΔH/(mΔT); cp = ΔS/[m ln(T2/T1)]; cp = cv + Tvα²/κ_T; cp = γ cv |

## 2. Running it

```bash
# one-time setup (Windows; on Mac/Linux use  source .venv/bin/activate)
python -m venv .venv
.venv\Scripts\activate
pip install -e .[dev]

pytest                                   # 185 tests, ~20 s
python validation/run_all.py             # all literature comparisons, figures and tables (~25 s)
python -m cp_engine                      # question-and-answer mode in the terminal
jupyter notebook notebooks/CP_Engine_Demo.ipynb
```

From Python:

```python
from cp_engine import CpEngine
e = CpEngine()
r = e.solve({"T": "25 C", "mass_fractions": {"water": 0.6, "ethanol": None}})
print(r)                 # full report
print(r.derivation())    # the chain of equations used
r.breakdown_table()      # ideal / excess / latent parts
```

**Colab:** open the notebook and run the first cell. It either clones the repo (set `REPO_URL`) or asks
you to upload `cp_engine_project.zip`, which `python tools/make_colab_zip.py` creates.

## 3. Demo script (about 10 minutes)

Run the notebook top to bottom beforehand so all outputs are visible. Then show:

1. **Section 3, Route A**, water + ethanol with one fraction left blank. Point out the filled fraction,
   the classification "real liquid solution", and the breakdown into ideal and excess parts. Run
   `print(res.derivation())` to show the rule chain.
2. **Section 4, Route B**, the same Cp from Q = 83.7 kJ heating 2 kg by 10 K. Then the adiabatic
   stirring case (Q = 0, ΔH = W).
3. **Section 6, several routes at once.** Q and ΔH are both given: the engine checks they are
   consistent, then compares the measured route with the model (+0.01 %). Then the deliberately
   inconsistent input, which produces a warning.
4. **Section 8.2, nanofluid.** Model I vs Model II on the same screen.
5. **Section 8.4, ice slurry** at −10 °C: the apparent Cp is about 13 kJ/(kg K) because ice keeps
   melting as T rises. Then the interval mean from −15 to +5 °C.
6. **Section 10, error handling.** Typo in an input name, fractions summing to 1.2, solid water at
   350 K, energy given in kg: each gives a clear message, never a traceback.
7. **Section 8.6, a boiling mixture.** Water + ethanol at 84 °C: bubble and dew points, vapour fraction
   and the latent part. Then 90 °C (all vapour) and 85 °C with 5 % ethanol (still liquid, although above
   the boiling point of pure ethanol).
8. **Section 8.7, the lecture's own numbers.** The CO₂ example of the EOS lectures (P_sat by vdW, RK, SRK,
   PR; Z of CO₂ + CH₄) comes out of our code to within ±0.01 (rounding); then steam at 100 bar (ideal gas 17 % off,
   with the Peng–Robinson residual 3 %).
9. **Section 8.8, partial molar Cp and Gibbs–Duhem**, and the uncertainty printed with every result.
10. **Section 13, validation table and figures** (see section 6 below for what to say).

## 4. Code map

```
user dict ─► inputs.py (normalise: units → SI, composition → moles)
          ─► validator.py (physical checks)
          ─► solver.py (forward chaining over rules.py)
                 ├─ helper rules (ΔT, T_mean, m ↔ n ...)
                 ├─ Route B/C rules (first law, ΔS, Cv, γ)
                 └─ Route A rules ─► system.py (classify, split phases, H(T), Cp)
                                        ├─ models/pure.py      Cp_i(T), H_i(T)
                                        ├─ models/ideal.py     Σ x_i Cp_i
                                        ├─ models/nrtl.py      G^E, H^E, Cp^E, γ
                                        ├─ models/heterogeneous.py  Model I / II
                                        ├─ models/phase_change.py   melting, SLE, solubility
                                        ├─ models/vle.py       bubble/dew points, flash (γ–φ)
                                        ├─ models/eos.py       virial + cubic EOS: Cpᴿ, φ̂, P_sat
                                        ├─ models/partial_molar.py  B̄_i, Gibbs–Duhem, ΔG/ΔH/ΔS_mix
                                        └─ models/cv_cp.py     Cp − Cv = Tvα²/κ_T
          ─► uncertainty.py (expected deviation from the validation results)
          ─► CpResult (Cp mass and molar, routes, breakdown, warnings, missing inputs)
data/components.json, data/binary_nrtl.json ─► database.py (loads, converts to SI)
```

| File | One-line purpose |
|---|---|
| `units.py` | Parses "25 °C", "4.18 kJ/(kg K)", (1, "atm"). Absolute temperatures get offsets; temperature *differences* do not |
| `database.py` | Loads the JSON, evaluates Cp and density correlations, checks validity ranges (warns outside; refuses beyond 30 K) |
| `inputs.py` | Accepts every input form, fills one missing fraction, converts specific energies, sets a calculation basis if needed |
| `validator.py` | T > 0, P > 0, fractions in [0, 1], phases physically possible, process consistent |
| `rules.py` | Every relation as a `Rule(inputs → outputs, equation)` |
| `solver.py` | Forward chaining, consistency checks, route labels, missing-input search, `CpEngine` API |
| `system.py` | Route A: behaviour of each component, liquid-phase grouping, phase split at T, enthalpy, Cp breakdown |
| `models/*.py` | The thermodynamic models (section 5) |
| `models/eos.py` | Virial (Pitzer/Abbott B, Prausnitz cross B₁₂) and generic cubic (vdW, RK, SRK, PR): Z roots, Hᴿ, Cpᴿ, φ̂ᵢ, P_sat, φ–φ bubble pressure |
| `models/vle.py` | γ–φ vapour–liquid equilibrium of liquid mixtures: bubble T, dew T, isothermal flash (Rachford–Rice) |
| `models/partial_molar.py` | Partial molar properties, Gibbs–Duhem check, ln γ from Gᴱ, mixing functions |
| `uncertainty.py` | Expected deviation of each result, built from the validation table |
| `plotting.py` | Publication figures (PNG and PDF, 300 dpi) |
| `interactive.py` | ipywidgets form and the guided `input()` mode |
| `validation/` | Five validation cases, `run_all.py`, data, results |
| `tools/` | Data extraction (NIST ThermoML: Cp, VLE, γ∞, LLE; digitiser), IAPWS ice, NRTL and octadecane fits, zip and doc generators |

## 5. The thermodynamics behind each module

### 5.1 Pure components (`database.py`, `models/pure.py`)
* Cp correlations: DIPPR eq. 100 (polynomial, liquids), DIPPR eq. 107 (ideal gas), Shomate (solids,
  NIST), IAPWS-06 fit (ice).
* **Enthalpy by Kirchhoff's law.** The reference is the liquid at T_m (H = 0). Then
  H_L(T) = ∫Cp_L dT, H_S(T) = −ΔH_fus + ∫Cp_S dT and H_G = H_L(T_b) + ΔH_vap + ∫Cp_G dT.
  The enthalpy of melting at any T is ΔH_fus + ∫(Cp_L − Cp_S)dT, so it changes with T.
* Boiling point at pressure P: Clausius–Clapeyron, 1/T_b(P) = 1/T_b − (R/ΔH_vap) ln(P/P°).
* Density and α always come from one density model (α = −(1/ρ)dρ/dT), so they are mutually consistent.

### 5.2 Ideal mixture (`models/ideal.py`)
Hᴱ = 0, so H = Σ xᵢHᵢ and Cp = Σ xᵢCpᵢ (per mole) ⇔ cp = Σ wᵢcpᵢ (per kg). The two forms are the same
statement, because wᵢ = xᵢMᵢ/M.

### 5.3 Real liquid solution (`models/nrtl.py`)
* NRTL: Gᴱ/RT = Σᵢ xᵢ (Σⱼ τⱼᵢGⱼᵢxⱼ)/(Σₖ Gₖᵢxₖ), with τᵢⱼ = aᵢⱼ + bᵢⱼ/T + eᵢⱼ ln T and Gᵢⱼ = exp(−ατᵢⱼ).
* **Gibbs–Helmholtz:** Hᴱ = −RT² ∂(Gᴱ/RT)/∂T, and **Cpᴱ = ∂Hᴱ/∂T**. So Cp = Σ xᵢCpᵢ + Cpᴱ.
* Derivatives by central differences with step h = 10⁻⁴T (≈ 0.03 K). Truncation error grows as h²
  and round-off as ε/h²; this step balances the two. `derivative_convergence()` shows the plateau,
  and a test checks against exact sympy derivatives.
* **Key finding:** NRTL parameters from VLE fitting fix Gᴱ but not its temperature derivative, so they
  predict Cpᴱ six or more times too small (about 6× for ethanol, about 9× for methanol). We refitted the parameters to measured Hᴱ
  (`tools/fit_nrtl_he.py`) and use those by default (`nrtl_parameters="HE"`). For water + ethanol the
  set is fitted **jointly** to Hᴱ (298–383 K) and to the freezing points (222–272 K), so one liquid model
  gives both the heat of mixing (→ Cpᴱ) and the ice–solution equilibrium (→ latent heat of ice slurries).
* **Stability check:** the refit is accepted only if ΔG_mix/RT stays convex in x from 220 to 390 K.
  A non-convex curve would mean a false liquid–liquid split, which is wrong for water–ethanol.
* **Held-out validation (case 2b):** every fit is repeated with part of the data left out (one isotherm,
  one paper, or one third of the freezing points) and the left-out part is predicted.
* **Ethanol + toluene** (water-free pair) is fitted jointly to isobaric VLE (68 points, Chen et al. 2012,
  activity coefficients from y P = x γ P_sat with NIST Antoine P_sat) and to measured Cpᴱ (Vega-Maza 2011).
  Mixture Cp error falls from 8.3 % (ideal) to 1.3 %. A fit to the VLE alone already predicts Cpᴱ to
  7 J/(mol K), but only the joint fit gets the Cpᴱ peak right.
* **Runtime stability warning:** every NRTL pair the engine uses is checked for a spurious liquid–liquid
  split at the working temperature. The literature acetone + water set fails above about 360 K, and the
  engine says so.
* **Infinite-dilution check (independent data):** γ∞ of ethanol in water is 35 % off with our Cp-optimised
  set and 12 % off with the VLE set. Adding γ∞ to the fit was tried: γ∞ improved, but held-out Cpᴱ and
  freezing points got worse, so we kept the Cp-optimised set (this is a Cp engine) and report the trade-off.
* Pairs without parameters are treated as ideal, with a warning.

### 5.4 Heterogeneous mixtures (`models/heterogeneous.py`, `system.py`)
* Enthalpy is extensive and each phase keeps its own properties: H = Σ m_p h_p, so
  **cp = Σ w_p cp_p (mass fractions)**.
* Nanofluid **Model II** (thermal equilibrium):
  cp = [φρ_s cp_s + (1−φ)ρ_f cp_f]/[φρ_s + (1−φ)ρ_f]. This is exactly mass weighting, since w_s = φρ_s/ρ_mix.
* **Model I**: cp = φ cp_s + (1−φ) cp_f. This weights a per-kilogram quantity by a *volume* fraction,
  which is dimensionally inconsistent unless ρ_s = ρ_f. For CuO (ρ_s/ρ_f ≈ 6) it is off by up to 63 %.
* Emulsions: immiscible liquids form separate liquid phases, and each phase is its own (real) solution.
  How a solute splits between the two phases is not modelled (the engine warns).
* **Miscibility can depend on T and x.** For methanol + n-hexane the engine reads the measured
  coexistence curve (Trejo et al. 2006, UCST 314.9 K): below the curve two liquids, above it one. If the
  T₁–T₂ window crosses the curve the engine warns, because the heat of demixing is not included.

### 5.5 Phase change and apparent Cp (`models/phase_change.py`, `system.py`)
* Total enthalpy: H(T) = Σ_phases Σᵢ nᵢᵖHᵢᵖ(T) + Σ_liquid nᵖHᴱ,ᵖ.
* By the chain rule, **dH/dT = Σ nᵢCpᵢ (ideal) + Σ nᵖCpᴱ (excess) + Σ (∂H/∂nᵢ)(dnᵢ/dT) (latent)**.
  The latent part is computed as [H(T; n(T+h)) − H(T; n(T−h))]/2h: only the phase *amounts* change,
  with properties frozen at T.
* How the phase amounts n(T) are found:
  * **Pure substance melting** (PCM, ice in pure water): a smoothed step
    f(T) = ½[1 + erf((T − T_m)/(σ√2))], σ = width/4. The peak area is exactly ΔH_fus.
  * **Solid in equilibrium with its solution** (ice in aqueous ethanol): ln(x γ) = −(ΔH_fus/R)(1/T − 1/T_m)
    + (ΔCp/R)(T_m/T − 1 − ln(T_m/T)). The ΔCp term keeps this consistent with the Kirchhoff enthalpies.
    γ comes from the **same NRTL** as Hᴱ when the parameter set was fitted to freezing points
    (`sle_activity="auto"`, default); `sle_activity="ideal"` forces γ = 1.
  * **Salt in its saturated solution**: tabulated solubility, with cp_app = Σ wᵢcpᵢ + Δh_sol·dw_dissolved/dT.
* **Electrolytes (NaCl in water):** Cp(solution) = n_w Cp_w + n_s φCp(m), with the apparent molar heat
  capacity φCp = c₀ + d₀√m + e·m (Debye–Hückel form) fitted to measured brine Cp. Dissolved salt *lowers*
  Cp (φCp < 0 at low m), which the old solid-Cp (Neumann–Kopp) estimate missed by 3 %.
* **Point vs interval near a melting point:** inside a smoothed melting step the point Cp depends on the
  assumed width (peak ∝ 1/width). The engine warns and also reports two width-independent numbers: the
  latent heat of the transition and the mean Cp over a fixed window.
* **Interval mean:** cp = [H(T2) − H(T1)]/[m(T2 − T1)]. This includes all latent heat in the interval and
  does not depend on the smoothing width.

### 5.6 Routes B and C (`rules.py`)
* Sign convention (IUPAC): Q > 0 heat **into** the system, W > 0 work **on** the system, ΔU = Q + W.
  W is shaft or stirring work; expansion work is W_pv = −∫P dV.
* Constant P: ΔH = ΔU + PΔV = Q + W_s + W_pv + PΔV = Q + W_s, since W_pv = −PΔV. Hence cp = ΔH/(mΔT).
* Constant V: ΔU = Q + W gives cv = ΔU/(mΔT), then cp = cv + Tvα²/κ_T (an exact Maxwell-relation result).
* ΔS route: at constant P with constant cp and no phase change, ΔS = m cp ln(T2/T1).
* Mixture α and κ_T come from volume-fraction averaging (ideal volume of mixing).

### 5.7 Pressure (`database.py`)
* All Cp correlations are for 1 atm. At other pressures the engine uses the exact relations
  **(∂Cp/∂P)_T = −T(∂²v/∂T²)_P** and **(∂H/∂P)_T = v − T(∂v/∂T)_P**, with v(T, P) from the **Tait equation**
  v = v₀(T)[1 − C ln(1 + ΔP/(B + P₀))], C = 0.0894, B(T) = C/κ_T(T) − P₀.
* κ_T(T) for water, methanol, hexane, toluene and dodecane comes from NIST (κ_T = (Cp/Cv)/(ρw²)).
* For an ideal gas α = 1/T makes the correction exactly zero, as it should.
* Validated against NIST isotherms 0.1–100 MPa: AARD 0.39 % (max 1.2 %), against 1.85 % (max 5.2 %)
  without the correction. Water at 500 bar: 4056 J/(kg K) vs NIST 4060.

### 5.8 User-supplied data (`custom.py`)
* `custom_components={"ethylene glycol": {"MW": "62.07 g/mol", "cp": "2.43 kJ/(kg K)", ...}}` adds a
  component for one calculation; `pure_cp={"water": "4.18 kJ/(kg K)"}` overrides a database Cp.
  The database files are never changed.

### 5.9 Solver (`solver.py`)
* Each rule fires at most once. "Check" rules also fire when their output is already known, and
  compare the two values (tolerance 2 %).
* Recommended value: measured routes (ΔH, Cv, γ, ΔS) take precedence over the model. The model
  value is still reported, along with the % difference.
* Missing inputs: backward search over the rule graph for the smallest sets of user-suppliable inputs.

### 5.10 Real gases: equations of state (`models/eos.py`)
* Gas Cp = **ideal-gas Cp** (database correlation, checked against NIST to ≤ 1.25 %) **+ residual Cpᴿ(T, P)**
  from an equation of state. Default Peng–Robinson; options `gas_eos="SRK" | "RK" | "vdW" | "virial" | "ideal"`.
* **Virial (lecture "EOS for pure fluids", slides 5–8):** Z = 1 + BP/RT, B from Pitzer/Abbott
  (BP_c/RT_c = B⁰ + ωB¹). Then Hᴿ = P(B − T dB/dT) and **Cpᴿ = −TP d²B/dT²** (derive it from
  (∂Cp/∂P)_T = −T(∂²V/∂T²)_P with V = RT/P + B). Mixtures: B_mix = ΣΣ yᵢyⱼBᵢⱼ with the Prausnitz combining
  rules for B₁₂ (lecture "EOS for mixtures", slides 7–9).
* **Generic cubic (slide 16):** P = RT/(V − b) − a(T)/[(V + εb)(V + σb)] with the vdW, RK, SRK and PR
  constants. Solve the cubic in Z (largest root = vapour, smallest = liquid; slide 17), then
  Hᴿ/RT = Z − 1 + (d ln a/d ln T − 1) q I and Cpᴿ = (∂Hᴿ/∂T)_P. Mixtures: one-fluid rules
  a = ΣΣ xᵢxⱼ√(aᵢaⱼ)(1 − kᵢⱼ), b = Σ xᵢbᵢ.
* **Critical constants and ω are data, not looked up blindly:** T_c, P_c, V_c from the critical point of
  the NIST reference EOS; ω from its definition, ω = −1 − log₁₀(P_sat/P_c) at T_r = 0.7, using NIST
  saturation data (water 0.344, methanol 0.565, n-hexane 0.300, toluene 0.266, ethanol 0.643).
* **The code reproduces the lecture's worked examples to within ±0.01 of the printed values (rounding)** (tests in
  `tests/test_eos_vle.py`): Abbott B of CO₂ at 350 K = −83.73 cm³/mol; P_sat and V_liq of CO₂ at 280 K
  from vdW/RK/SRK/PR (52.55/45.62/41.99/41.60 bar); Z of equimolar CO₂ + CH₄ at 10–60 bar (virial, PR
  k₁₂ = 0 and 0.10); B₁₂ = −47.0 cm³/mol; φ̂_CO₂ = 0.947; and the φ–φ bubble pressure: k₁₂ = 0 puts
  18.6 % CH₄ in the liquid at 6.375 MPa, k₁₂ = 0.10 gives the measured 12.5 %.
* **Validation against NIST (case 6):** vapour Cp of water, methanol, n-hexane, toluene at 0.5–10 MPa
  (189 points): ideal gas 12.1 % AARD (max 67 %), **+ PR residual 5.3 %**, + virial 5.2 %. Steam at 10 MPa,
  500 °C: 2496 vs NIST 2583 J/(kg K) (ideal gas: 2131). Largest errors are within ~20 K of saturation,
  where every cubic EOS is weak (mean-field critical behaviour, lecture slide 19).
* Vapour pressures from the EOS (equal fugacities) vs NIST, 0.55–0.95 T_c: vdW 224 %, RK 58 %,
  **SRK 2.3 %, PR 1.4 %**: Soave's idea of fitting α(T) to vapour pressures, seen on our own fluids.

### 5.11 Boiling of mixtures: vapour–liquid equilibrium (`models/vle.py`, `system.py`)
* Before this, each component "boiled" at its own pure boiling point, so water + ethanol at 90 °C came out as
  liquid water + gaseous ethanol. Now every miscible liquid mixture with a volatile component gets its
  **bubble and dew points** at P:
  below the bubble point → liquid (even above a pure boiling point); above the dew point → vapour;
  in between → **isothermal flash** at each T.
* **γ–φ equilibrium** (lecture "Property relationships for mixtures", slides 24–28): f̂ᵢᴸ = f̂ᵢⱽ gives the
  modified Raoult's law yᵢP = xᵢγᵢPᵢˢᵃᵗ (ideal vapour; option `vapour_model="lewis_randall"` adds the virial
  φᵢ). γᵢ from NRTL; Pᵢˢᵃᵗ from Antoine equations fitted to NIST saturation data (≤ 0.3 % for water).
* Bubble point: Σ Kᵢxᵢ = 1; dew point: Σ yᵢ/Kᵢ = 1; flash: Rachford–Rice Σ zᵢ(Kᵢ − 1)/[1 + V(Kᵢ − 1)] = 0,
  iterated because γ depends on the liquid composition.
* Non-volatile components (glycerol, salts, n-dodecane) stay in the liquid: K = 0, no dew point, and the
  bubble point rises (boiling-point elevation).
* Enthalpy in the two-phase region: liquid at x (with Hᴱ) + vapour at y (ideal-gas path from T_b + Hᴿ). The
  apparent Cp = dH/dT includes the latent heat of whatever evaporates as T rises; ∫cp dT = ΔH to 0.03 %.
* **Validation (case 6):** bubble points of 26 papers in the ThermoML archive, 13–101 kPa:
  water + ethanol **0.26 K** (256 points), water + methanol **0.46 K** (153 points). The azeotrope comes out
  at x = 0.90, 351.3 K. Two papers that disagree with all others by > 2 K on average are excluded (listed
  in `validation/case6_vle_real_gas.py`).
* **Two NRTL sets, on purpose.** The Hᴱ-fitted sets (best for Cp) give bubble points 4.5 K (ethanol) and
  68 K (methanol!) off; the VLE sets give 0.3–0.5 K but Cpᴱ 3–6 times too small. We tried one set fitted to
  Hᴱ + freezing points + bubble points together: Cp errors doubled (water + methanol 2.7 → 5.3 %). So caloric
  properties use the Hᴱ set and the phase split uses the VLE set (`vle_parameters`). Ethanol + toluene has
  one set fitted to both.
* Immiscible volatile liquids (water + toluene) boil together where Σ P_bubble = P (84 °C, steam
  distillation). This is detected and warned about, not modelled.

### 5.12 Partial molar properties, Gibbs–Duhem, mixing functions (`models/partial_molar.py`)
* B̄ᵢ = (∂B/∂Nᵢ)_{T,P,Nⱼ}; Euler: B = Σ Nᵢ B̄ᵢ. From a molar property: B̄ₖ = B + ∂B/∂xₖ − Σⱼ xⱼ ∂B/∂xⱼ
  (lecture slides 6–9). Both forms are coded and give the same numbers.
* γ is the partial molar excess Gibbs energy: RT ln γᵢ = ∂(NGᴱ)/∂Nᵢ. The code differentiates Gᴱ numerically
  and recovers the analytic NRTL ln γ to 10⁻⁷ (test), so the γ's are thermodynamically consistent.
* **Gibbs–Duhem** Σ xᵢ (∂ ln γᵢ/∂xⱼ)_{T,P} = 0 is checked (residual < 10⁻⁸).
* Gibbs–Helmholtz in partial molar form: H̄ᵢᴱ = −RT² ∂ ln γᵢ/∂T; C̄p,ᵢᴱ = ∂H̄ᵢᴱ/∂T.
* Every liquid solution result reports the **partial molar Cp** of each component (C̄p,ᵢ = Cp,ᵢ + C̄p,ᵢᴱ,
  with Σ xᵢC̄p,ᵢ = Cp) and the **mixing functions** ΔG_mix, ΔH_mix (= Hᴱ), ΔS_mix (`r.breakdown.extras`).

### 5.13 How accurate is a result? (`uncertainty.py`)
* Every Route-A result carries an **expected deviation from measurements** (`r.uncertainty`, also printed).
* Each ingredient's deviation in our own validation (e.g. water + methanol mixture model 2.7 %, Model II
  suspensions 2.1 %, Peng–Robinson residual 5.3 %, pure water 0.1 %) is multiplied by its share of Cp and the
  contributions are combined as √Σuₖ². `data/model_uncertainty.json` is generated from the validation table
  (`tools/make_uncertainty_table.py`), so every number traces back to a comparison with published data.
* Examples: pure water ± 0.1 %; water + methanol ± 2.4 %; water + ethanol ± 3.5 %; CuO nanofluid ± 2.1 %;
  acetone + water ± 6.5 % (pair not validated for Cp); steam at 100 bar ± 4.8 %.
* Inside a phase change the boundary location dominates (freezing points ± 2.5 K, bubble points ± 0.3 K);
  this is stated, not hidden in the number.

### 5.14 Where each CL 202 lecture concept is in the code
| Lecture concept | Where | What it does in the engine |
|---|---|---|
| System, boundary, closed/open, adiabatic/diathermal (Ch. 1) | `rules.py`, `inputs.py` (`process`) | Route B: adiabatic → Q = 0; constant P → ΔH = Q + W; constant V → ΔU = Q + W |
| Postulate I: n + 2 independent properties fix the state | `solver.py` missing-input search | Route A needs T, P and the n amounts; with fewer, the engine says which ones are missing |
| First law, dE = δQ + δW, H = U + PV (Ch. 1) | `rules.py` | ΔH = ΔU + PΔV; IUPAC sign convention (Q, W > 0 into the system) |
| Entropy, dS = δQ_rev/T | `rules.py` | cp = ΔS/[m ln(T₂/T₁)] (reversible heating at constant P) |
| Partial molar properties, Euler, Gibbs–Duhem (Ch. 4) | `models/partial_molar.py` | Partial molar Cp; Gibbs–Duhem test of the NRTL γ |
| Mixing functions, ideal solution (ΔH_mix = 0, ΔS_mix = −RΣx ln x) | `models/ideal.py`, `partial_molar.py` | Ideal fallback; ΔG/ΔH/ΔS_mix reported |
| Excess functions, activity coefficient, RT ln γᵢ = Ḡᵢᴱ | `models/nrtl.py` | Cpᴱ = ∂Hᴱ/∂T, Hᴱ = −RT² ∂(Gᴱ/RT)/∂T |
| Fugacity, Lewis–Randall rule, Henry's law (infinite dilution) | `models/vle.py`, validation case 2 | γ–φ VLE; γ∞ (Henry reference) compared with measurements |
| Virial EOS, Pitzer/Abbott B | `models/eos.py` | Cpᴿ = −TP d²B/dT²; lecture CO₂ numbers reproduced |
| Cubic EOS (vdW, RK, SRK, PR), roots, φ_L = φ_V | `models/eos.py` | Real-gas Cpᴿ; P_sat; lecture table at 280 K reproduced |
| Mixing and combining rules, k₁₂, B₁₂, φ̂ᵢ, φ–φ bubble pressure | `models/eos.py` | Gas-mixture Cpᴿ; CO₂ + CH₄ example of slides 23–26 reproduced |
| (∂Cp/∂P)_T = −T(∂²V/∂T²)_P (Maxwell relation) | `database.py` (Tait), `eos.py` (gases) | Pressure dependence of liquid and gas Cp |

## 6. Validation: what we can claim (all against real data)

"Fit" means the data were used to fit parameters; "held-out" and "independent" mean they were not.

| # | Test | Result | What it shows |
|---|---|---|---|
| 1 | 5 pure liquids vs NIST SRD 69 (128 points) | AARD 0.33 %, max 1.25 % | The pure-component data are correct |
| 1 | 5 liquids at 0.1–100 MPa vs NIST (48 points) | **0.39 %** with Tait correction (1.85 % without) | The pressure correction works |
| 1 | Water: cp from NIST cv via cp = cv + Tvα²/κ_T | AARD 0.02 % | The Cp–Cv identity and our α, κ_T are right |
| 1 | Ideal-gas Cp of 5 components vs NIST, 200–1000 K | AARD 0.10–0.46 %, max 1.25 % | The vapour-phase correlations are right |
| 2 | Water + ethanol Hᴱ (88 points) | fit 34 J/mol; **held-out 46 J/mol**; VLE parameters 256 | VLE parameters are not caloric parameters |
| 2 | Water + methanol Hᴱ (151 points) | fit 23 J/mol; **held-out 48 J/mol**; VLE 141 | Same for methanol |
| 2b | Water + ethanol Cpᴱ, source paper held out of the fit | **3.7 J/(mol K)** RMSD (VLE 9.1); peak ≈ 14.5 | Genuine prediction of Cpᴱ |
| 2 | Water + methanol Cp (independent data) | ideal 7.5 %, VLE-NRTL 6.5 %, **Hᴱ-NRTL 2.7 %** | Gibbs–Helmholtz route to Cpᴱ works |
| 2 | NaCl brine Cp, independent point | **0.17 %** (old Neumann–Kopp model 2.9 %) | Electrolyte model works |
| 2 | Glycerol + water Cp, 162 points (x ≤ 0.08) | 0.52 % with the ideal fallback | Ideal mixing is adequate for dilute glycerol |
| 2 | Ethanol + toluene Cp (NRTL fitted to VLE + Cpᴱ) | fit **1.3 %** (max 3.1 %); ideal mixing 8.3 % | Binary parameters matter even for Cp |
| 2b | Ethanol + toluene, held-out | Cpᴱ half held out **3.0 J/(mol K)**; Cpᴱ from VLE alone 7.1; VLE pressure held out RMS ln γ 0.04 | The fit predicts data it did not see |
| 2 | Ethanol + toluene VLE bubble T (68 points) | mean 0.6 K, max 2.6 K | The same parameters give the phase equilibrium |
| 2 | γ∞ in water (independent) | ethanol: Cp set 35 %, VLE set 12 %; acetone 25 % | Honest: a Cp-optimised set is not a γ∞ model |
| 3 | Nanofluids, 69 DSC points (Al₂O₃, CuO, SiO₂) | **Model II 2.1 %**, Model I 18.5 % (max 63 %) | Mass weighting is right; Model I is wrong |
| 4 | Freezing points ≤ 30 wt% ethanol | consistent NRTL: fit 1.7 K, **held-out 2.5 K**; ideal + ΔCp 0.84 K | The SLE behind ice-slurry Cp works |
| 4 | Energy conservation through melting peaks | ≤ 0.006 % | Latent heat is neither lost nor double-counted |
| 5 | Route A → Q → Route B, 7 systems | < 10⁻⁶ % | The routes are mutually consistent |
| 6 | Bubble points, water + ethanol / methanol (409 points, 24 papers) | **0.26 K / 0.46 K** (Hᴱ sets: 4.5 K / 68 K) | The engine knows when a mixture boils |
| 6 | Vapour Cp at 0.5–10 MPa, 4 fluids (189 points) | ideal gas 12.1 % → **PR residual 5.3 %** | Real-gas correction works (weakest near saturation) |
| 6 | P_sat from vdW / RK / SRK / PR, 5 fluids | 224 / 58 / **2.3 / 1.4 %** | Why α(T) was fitted to vapour pressures |
| 6 | ∫cp_app dT vs ΔH through a boiling range | ≤ 0.03 % | Latent heat of boiling counted exactly once |
| – | Lecture worked examples (CO₂, CO₂ + CH₄) | identical to within ±0.01 of the printed values (rounding) | The EOS code implements the course equations correctly |

What to say about the less good numbers:
* *Glycerol 2.4 % at 55 °C:* the literature values themselves scatter by 1.4 % at 25 °C.
* *Methanol–water at 280 K overshoots at high x:* 280 K is below the Hᴱ data range (298–378 K).
* *Held-out errors are larger than fit errors* (46 vs 34 J/mol): that is normal and is exactly why we report them.
* *The ideal + ΔCp freezing equation (0.84 K) beats the consistent NRTL (2.5 K held-out):* NRTL with 6
  parameters must describe Gᴱ from 220 to 383 K. We chose consistency (one liquid model for Hᴱ and the
  freezing equilibrium) as the default; `sle_activity="ideal"` gives the more accurate freezing point.
* *Freezing points above 40 wt% ethanol are off by several kelvin:* the liquid is deeply supercooled
  (−30 to −50 °C), far outside the range of the liquid-water Cp data. Ice slurries use 5–30 wt%.
* *Dodecane under pressure:* 0.62 % with the correction vs 0.21 % without (only 5 points up to 40 MPa);
  for the other four liquids the correction cuts the error 2–18-fold.
* *γ∞ of ethanol in water is 35 % off with our default set:* γ∞ is the x → 0 limit of Gᴱ, a region with
  almost no weight in Hᴱ. Adding γ∞ to the fit made Cpᴱ and freezing points worse (case 2b), so for a
  Cp engine we keep the Cp-optimised set; the VLE set is available (`nrtl_parameters="VLE"`).
* *One NaCl dataset was excluded* (a buffer study whose NaCl values disagree with the dedicated study by
  1.1 %); the reason is written in `tools/fit_apparent_cp.py`.

## 7. Where the data came from (short version; full detail in `docs/DATA_SOURCES.md`)

* **Downloaded:** NIST SRD 69 fluid tables; the NIST ThermoML archive (189 MB, 11,923 files), from
  which we extracted 876 binary-mixture points with the DOI of each original paper.
* **Digitised:** O'Hanley et al. 2012 (open access) gives its nanofluid data only in figures. Our script
  reads the marker positions directly from the PDF's vector drawing commands, so there is no eyeballing.
* **Implemented a standard:** IAPWS R10-06 for ice, checked to 9 digits against the official test values.
* **Corrected during verification:** acetone–water NRTL (i and j swapped), glycerol Cp (4 % high),
  ice Cp (1.4 % low), CuO range, solid densities, NaCl solubility.
* **Pressure and electrolyte data:** NIST isotherms (0.1–100 MPa) for 5 liquids; NaCl(aq) Cp and glycerol
  apparent molar Cp from the ThermoML archive (molality data are now converted too).
* **Phase-equilibrium data (ThermoML):** ethanol + toluene VLE, γ∞ in water, methanol + alkane LLE, and
  502 bubble points of water + ethanol / methanol from 26 papers.
* **NIST saturation tables and superheated-vapour isobars** for 5 fluids (`tools/fetch_nist_fluid.py`): the
  source of the Antoine equations, critical constants, ω and the real-gas validation.
* **n-octadecane (PCM):** solid Cp and liquid density fitted to measurements (Fonseca 2014, Espeau 2006).
* **Open:** low-impact flags (`docs/DATA_TO_VERIFY.md`, none can change a liquid Cp result by more than
  a few tenths of a percent); direct Cp data for water + ethanol and DSC data for PCM emulsions (paywalled).

## 8. Limitations (say these before you are asked)

1. Cpᴱ needs binary parameters. Four pairs are validated (water + ethanol, water + methanol and
   ethanol + toluene fitted; for dilute water + glycerol ideal mixing is adequate). Without parameters
   the engine falls back to ideal mixing and warns (it was 8 % for ethanol + toluene before the fit).
2. Boiling: one liquid phase + vapour is modelled (γ–φ, ideal vapour). Three-phase boiling of immiscible
   liquids (water + toluene) is detected and warned about, not modelled. Caloric properties and the phase
   split use different NRTL sets for water + alcohols (a measured trade-off, section 5.11).
3. Real gases: cubic EOS with k_ij = 0; errors up to 44 % for steam within ~20 K of saturation at 10 MPa
   (5 % on average). Acetone's critical constants are NIST averages (flagged "verify").
4. No solute partitioning between two liquid phases in emulsions; no experimental validation for emulsions
   or PCM slurries (no open DSC data found) beyond the energy-conservation checks.
5. Electrolytes: only NaCl has an apparent-molar-Cp model, at 298 K (the engine warns elsewhere).
6. Pressure correction validated to 100 MPa for 5 liquids; other components use a constant κ_T.
7. Ideal volume of mixing for ρ, α and κ_T of mixtures; interfacial energy of nanoparticles is neglected.
8. Pure-component melting is smoothed over a finite width; point values inside the step depend on it
   (the engine warns and reports the width-independent latent heat and window mean).

## 9. Questions the professor may ask, with answers

### Concept and design
1. **Why is your engine "general"?**
   Because it is not a fixed formula. It is a network of 25 relations, and whatever subset of inputs is
   given, the solver chains the relations that apply. New relations are added as new rules.
2. **What happens if I give both Q and ΔH?**
   The first-law rule computes ΔH from Q (+ W) and compares it with the given ΔH. If they differ by more
   than 2 %, it warns "INCONSISTENT INPUT". The user's ΔH is used, and the comparison is reported.
3. **What if I give too little information?**
   The backward search lists the smallest sets of extra inputs that complete a route, for example
   "Route B: provide the total mass of the system".
4. **Why report both mass and molar Cp?**
   Engineers use J/(kg K) in energy balances; thermodynamic models work in J/(mol K). They are related by
   the mean molar mass M = Σ xᵢMᵢ, which needs the composition.

### Thermodynamics
5. **Derive Cpᴱ from Gᴱ.**
   Gibbs–Helmholtz gives Hᴱ = −RT²∂(Gᴱ/RT)/∂T. Then Cpᴱ = (∂Hᴱ/∂T)_P,x = −R(2Tg′ + T²g″), with g = Gᴱ/RT.
6. **Why did the VLE NRTL parameters fail for Cpᴱ?**
   VLE fixes Gᴱ (through activity coefficients) near the boiling range. Cpᴱ is a *second* temperature
   derivative of Gᴱ/RT, which VLE data barely constrain. Measured Hᴱ at several temperatures does
   constrain it. Our refit cut the Hᴱ error from 256 to 28 J/mol and the Cpᴱ error from 9.1 to 2.3 J/(mol K).
7. **Why did you need Hᴱ at several temperatures?**
   At a single T, different (b₁₂, b₂₁) pairs fit almost equally well (a test in the suite shows this);
   the temperature dependence is what pins them down.
8. **Why check stability after fitting?**
   A freely fitted NRTL can make ΔG_mix non-convex, which predicts a fake liquid–liquid split. Our first
   6-parameter fit did exactly that at 260 K; we only accept sets that are stable from 240 to 390 K.
9. **Show that Model II equals mass weighting.**
   w_s = φρ_s/(φρ_s + (1−φ)ρ_f). Substitute into cp = w_s cp_s + (1 − w_s) cp_f; you get Model II.
10. **Why is Model I wrong?**
    cp is per kilogram. Averaging it with volume fractions is only correct if both phases have the same
    density. The CuO data (ρ ratio ≈ 6) show Model I is off by up to 63 %.
11. **What is "apparent" Cp?**
    The measured dH/dT when material changes phase as T changes. It equals the sensible Cp plus
    Δh × (rate of phase change). Near the melting point this can be 10–100 times the sensible value.
12. **How do you know how much ice there is at −10 °C in 10 % ethanol?**
    From solid–liquid equilibrium. The chemical potential of ice equals that of water in the solution:
    ln(x_w γ_w) = −(ΔH_fus/R)(1/T − 1/T_m) + ΔCp term. Solving for x_w gives the liquid composition;
    the rest of the water is ice.
13. **Why add the ΔCp term?**
    Our enthalpies use Kirchhoff's law (ΔH_fus varies with T). The equilibrium equation must use the same
    physics to be consistent, and it also improves the freezing points (1.06 → 0.84 K mean error).
14. **Why does the apparent Cp of the ice slurry drop suddenly at the freezing point?**
    Above the liquidus there is no ice left to melt, so the latent part vanishes and only the solution's Cp remains.
15. **Prove ΔH = Q + W at constant P.**
    ΔU = Q + W_s + W_pv with W_pv = −PΔV at constant P. Then ΔH = ΔU + PΔV = Q + W_s.
16. **How do you get Cp from Cv?**
    Cp − Cv = Tvα²/κ_T, exact from the Maxwell relations. For water at 25 °C the difference is about 44 J/(kg K).
    Validated against NIST: 0.12 %.
17. **When is the ΔS route valid?**
    At constant P, with Cp constant over [T1, T2] and no phase change: ΔS = mCp ln(T2/T1). The engine
    warns if the mixture model predicts a phase change in the interval.
18. **Why does the PCM peak height depend on the melt width, and is that a problem?**
    A pure substance melts at one temperature, so the ideal peak is infinite. Real samples melt over a
    range, which the width represents. The *area* under the peak (the latent heat) is exact, so the
    interval-mean Cp is width-independent.

### Numerics and software
19. **Why h = 10⁻⁴T for derivatives?**
    Central difference: truncation error ∝ h², round-off ∝ ε/h². The total error is minimised near
    h ≈ ε^(1/4)T ≈ 10⁻⁴T. We demonstrate the plateau and cross-check with sympy.
20. **How do you compute the latent part without noise?**
    By differencing H only with respect to the phase amounts: H(T; n(T±h)), with properties frozen at T.
    H is then linear in n, so there is no integration noise.
21. **How did you test the code?**
    185 pytest tests: units, limits (φ→0, φ→1, xᵢ→1, Cpᴱ→0 at the pure ends), ideal reduction,
    bad input, N = 2, 3 and 5 components, all routes, round trips, and regression tests against the
    literature data.
22. **What happens with bad input?**
    Every error is a `CpEngineError` subclass with a specific message, e.g. "Water was declared solid,
    but at 350 K it is above its melting point". The solver never shows a raw traceback.
23. **How do you add a new component?**
    Add a JSON entry to `data/components.json` (with a source for each value). No code change is needed.

### Data
24. **Where did your data come from? Did you make anything up?**
    No. NIST SRD 69 (downloaded), the NIST ThermoML archive (downloaded and parsed; every point carries the
    DOI of its paper), the O'Hanley 2012 open-access paper (digitised from vector graphics), Lange's
    Handbook (freezing points), IAPWS (ice). See `docs/DATA_SOURCES.md`.
25. **How accurate is the digitisation?**
    The markers are read from the PDF drawing commands, so precision is limited only by the plotting
    (better than 0.002 J/(g K)). The axes are calibrated from the printed tick labels. A check: the
    Model II curve in the paper starts at 4.178 J/(g K) at φ = 0, which is water at 35 °C.
26. **What errors did you find in your own data?**
    The acetone–water NRTL parameters had i and j swapped. We caught it by checking infinite-dilution
    activity coefficients. Glycerol Cp and ice Cp correlations were also replaced. See DATA_SOURCES §5.

27a. **"Why is methanol + hexane two phases at 25 °C but one at 50 °C?"**
    The pair has an upper critical solution temperature of 314.9 K (Trejo et al. 2006). The engine reads
    the measured coexistence curve: below it two liquids, above it one solution. Across the curve the
    heat of demixing is not modelled, and the engine warns.
27b. **"How did you get NRTL parameters for ethanol + toluene?"**
    From isobaric VLE (Chen et al. 2012): γᵢ = yᵢP/(xᵢP_sat,i) with NIST Antoine P_sat, fitted together with
    measured Cpᴱ. Held-out tests: half the Cpᴱ points predicted to 3.0 J/(mol K); a fit to VLE alone predicts
    Cpᴱ to 7.1. Mixture Cp: 1.3 % instead of 8.3 % for ideal mixing.
27c. **"Your ethanol–water γ∞ is 35 % off. Isn't your NRTL wrong?"**
    It is fitted for heat capacity. γ∞ is the dilute limit of Gᴱ, which barely affects Hᴱ. We tried adding
    γ∞ to the fit: γ∞ improved but held-out Cpᴱ and freezing points got worse. We report both and let the
    user switch to the VLE set for phase equilibrium.

### Questions on the course concepts (from the CL 202 lecture slides)
33. **"Where is the equation of state from my lecture in your code?"**
    `models/eos.py`: virial (Pitzer/Abbott B) and the generic cubic (vdW, RK, SRK, PR) with mixing and
    combining rules. The tests reproduce your CO₂ example (B = −83.73 cm³/mol; P_sat at 280 K 52.55 / 45.62 /
    41.99 / 41.60 bar) and the CO₂ + CH₄ example (Z table, B₁₂ = −47.0, φ̂ = 0.947, 18.6 % vs 12.5 % CH₄).
34. **"Derive the pressure dependence of the Cp of a gas."**
    (∂Cp/∂P)_T = −T(∂²V/∂T²)_P. With the virial equation V = RT/P + B(T): Cpᴿ = −TP d²B/dT². For a cubic,
    differentiate Hᴿ/RT = Z − 1 + (d ln a/d ln T − 1) q I. Ideal gas: zero. Steam at 100 bar: +17 %.
35. **"Why do vdW and RK fail for vapour pressure but SRK and PR work?"**
    Two-parameter cubics fix a and b at the critical point; Soave fitted α(T_r, ω) to vapour pressures. On our
    5 fluids: vdW 224 %, RK 58 %, SRK 2.3 %, PR 1.4 % average error. The price is liquid volume (Z_c too high),
    which is why we never use a cubic for liquid densities or liquid Cp.
36. **"Your mixture is at 90 °C and ethanol boils at 78 °C. Is it liquid?"**
    It depends on the mixture's bubble and dew points, not the pure boiling points. 5 % ethanol boils at 90.3 °C,
    so at 85 °C it is liquid; 50 % is above its dew point at 90 °C, so it is all vapour. In between, the flash
    gives x, y and the vapour fraction, and the latent heat appears in the apparent Cp.
37. **"Why γ–φ and not φ–φ for water + ethanol?"**
    Your slide 15: no constant k_ij fits strongly non-ideal polar liquids. γ–φ uses a liquid-solution model
    (NRTL) that does fit them (bubble points 0.26 K). φ–φ is in the code too, and works for CO₂ + CH₄.
38. **"You use one NRTL set for Cp and another for boiling. Isn't that inconsistent?"**
    We tried one set fitted to everything (Hᴱ + freezing points + 409 bubble points). It doubled the Cp error
    (water + methanol 2.7 → 5.3 %), because NRTL cannot describe Gᴱ(x, T) of water + alcohols that well over
    220–383 K. So each property uses the set fitted to it, and we say so. Ethanol + toluene has one set for both.
39. **"What is a partial molar property? Show it in your output."**
    B̄ᵢ = (∂B/∂Nᵢ)_{T,P,Nⱼ}. `r.breakdown.extras["partial_molar_cp_J_per_molK"]` gives C̄p of water and ethanol;
    Σ xᵢC̄p,ᵢ = Cp (Euler). We also check Gibbs–Duhem and that ln γᵢ = ∂(NGᴱ/RT)/∂Nᵢ.
40. **"What is the uncertainty of the number you just printed?"**
    It is printed: e.g. water + methanol ± 2.4 %. It is built from our validation, each ingredient's measured
    deviation × its share of Cp, combined as √Σu². Inside a phase change the boundary location (± 0.3 K for
    bubble points, ± 2.5 K for freezing points) matters more, and we say that instead of hiding it.
41. **"Which postulate tells you how many inputs Route A needs?"**
    Postulate I: a simple system's stable equilibrium state is fixed by n + 2 independent properties. Route A
    needs T, P and the n amounts; if one is missing, the engine's backward search says which.

### Questions from our own mock review (prepare these)
27. **"Change P to 500 bar. Does Cp change? Derive (∂Cp/∂P)_T."**
    From dH = Cp dT + [v − T(∂v/∂T)_P] dP and the equality of mixed second derivatives,
    (∂Cp/∂P)_T = −T(∂²v/∂T²)_P. We evaluate v(T,P) with the Tait equation. Water at 500 bar: 4056 J/(kg K)
    (NIST 4060), 3 % below the 1-atm value. For an ideal gas the correction is exactly zero.
28. **"Your Cpᴱ test uses data you fitted to. How is that validation?"**
    In the main table it is labelled a consistency check. Case 2b repeats the fit *without* that paper and
    predicts its Cpᴱ: RMSD 3.7 J/(mol K), against 9.1 for VLE parameters. All fitted sets also have
    leave-one-isotherm-out tests.
29. **"Why does the ice-slurry equilibrium use γ = 1 when the enthalpy uses NRTL?"**
    It no longer does. The water + ethanol parameters were fitted jointly to Hᴱ and freezing points, and the
    same NRTL gives γ in the equilibrium (`sle_activity="auto"`). The price: held-out freezing-point error
    2.5 K instead of 0.84 K for the simple ideal equation, which remains available as an option.
30. **"Your brine Cp was 2 % high."**
    Dissolved NaCl has a negative apparent molar heat capacity (it orders the water around the ions), so the
    solid-Cp estimate is wrong. The engine now uses φCp(m) fitted to a dedicated study and is within 0.17 %
    of an independent measurement.
31. **"Can I use a component that is not in your database?"**
    Yes: `custom_components` (MW and Cp at minimum) or `pure_cp` to override a value, for that calculation only.
32. **"The PCM Cp at 301.3 K changes 8-fold with the melt width. Which is right?"**
    Neither number is physical on its own; a pure substance has an infinite Cp at T_m. The engine warns and
    reports the width-independent latent heat (48.5 kJ per kg of slurry) and the mean Cp over a fixed window.
    For design, use T1–T2 spanning the melting range.

## 10. Glossary

| Term | Meaning |
|---|---|
| AARD | Average absolute relative deviation, (100/N) Σ \|model − exp\|/exp |
| RMSD | Root-mean-square deviation |
| Excess property (Hᴱ, Cpᴱ, Gᴱ) | Real-solution value minus ideal-solution value |
| NRTL | Non-Random Two-Liquid activity-coefficient model (Renon & Prausnitz 1968) |
| SLE | Solid–liquid equilibrium (solubility / freezing-point depression) |
| Liquidus | Temperature at which the last solid disappears on heating |
| Apparent Cp | dH/dT including latent heat of material changing phase |
| Forward chaining | Repeatedly applying any rule whose inputs are known |
| Route | One chain of rules that reaches Cp |
| ThermoML | NIST's XML format and archive of published thermophysical data |
| IAPWS | International Association for the Properties of Water and Steam |
