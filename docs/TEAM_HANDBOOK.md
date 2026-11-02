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
| A: property | T (or T1, T2), composition | Σ xᵢCpᵢ; Cpᴱ from NRTL via Gibbs–Helmholtz; mass-weighting over phases; latent terms from phase equilibrium |
| B: calorimetric | Q and/or W, ΔT, amount | First law: ΔH = Q + W (const P); ΔU = Q + W (const V) |
| C: derived | ΔH, ΔU, ΔS, Cv or γ | cp = ΔH/(mΔT); cp = ΔS/[m ln(T2/T1)]; cp = cv + Tvα²/κ_T; cp = γ cv |

## 2. Running it

```bash
# one-time setup (Windows; on Mac/Linux use  source .venv/bin/activate)
python -m venv .venv
.venv\Scripts\activate
pip install -e .[dev]

pytest                                   # 127 tests, ~3 s
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
5. **Section 8.4, ice slurry** at −10 °C: the apparent Cp is about 15 kJ/(kg K) because ice keeps
   melting as T rises. Then the interval mean from −15 to +5 °C.
6. **Section 10, error handling.** Typo in an input name, fractions summing to 1.2, solid water at
   350 K, energy given in kg: each gives a clear message, never a traceback.
7. **Section 13, validation table and figures** (see section 6 below for what to say).

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
                                        └─ models/cv_cp.py     Cp − Cv = Tvα²/κ_T
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
| `plotting.py` | Publication figures (PNG and PDF, 300 dpi) |
| `interactive.py` | ipywidgets form and the guided `input()` mode |
| `validation/` | Five validation cases, `run_all.py`, data, results |
| `tools/` | Data extraction (NIST ThermoML, digitiser), IAPWS ice, NRTL refit, zip and doc generators |

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
  (`tools/fit_nrtl_he.py`) and use those by default (`nrtl_parameters="HE"`).
* **Stability check:** the refit is accepted only if ΔG_mix/RT stays convex in x from 240 to 390 K.
  A non-convex curve would mean a false liquid–liquid split, which is wrong for water–ethanol.
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
  * **Salt in its saturated solution**: tabulated solubility, with cp_app = Σ wᵢcpᵢ + Δh_sol·dw_dissolved/dT.
* **Interval mean:** cp = [H(T2) − H(T1)]/[m(T2 − T1)]. This includes all latent heat in the interval and
  does not depend on the smoothing width.

### 5.6 Routes B and C (`rules.py`)
* Sign convention (IUPAC): Q > 0 heat **into** the system, W > 0 work **on** the system, ΔU = Q + W.
  W is shaft or stirring work; expansion work is W_pv = −∫P dV.
* Constant P: ΔH = ΔU + PΔV = Q + W_s + W_pv + PΔV = Q + W_s, since W_pv = −PΔV. Hence cp = ΔH/(mΔT).
* Constant V: ΔU = Q + W gives cv = ΔU/(mΔT), then cp = cv + Tvα²/κ_T (an exact Maxwell-relation result).
* ΔS route: at constant P with constant cp and no phase change, ΔS = m cp ln(T2/T1).
* Mixture α and κ_T come from volume-fraction averaging (ideal volume of mixing).

### 5.7 Solver (`solver.py`)
* Each rule fires at most once. "Check" rules also fire when their output is already known, and
  compare the two values (tolerance 2 %).
* Recommended value: measured routes (ΔH, Cv, γ, ΔS) take precedence over the model. The model
  value is still reported, along with the % difference.
* Missing inputs: backward search over the rule graph for the smallest sets of user-suppliable inputs.

## 6. Validation: what we can claim (all against real data)

| # | Test | Result | What it shows |
|---|---|---|---|
| 1 | 5 pure liquids vs NIST SRD 69 (128 points) | AARD 0.33 %, max 1.25 % | The pure-component data are correct |
| 1 | Water: cp from NIST cv via cp = cv + Tvα²/κ_T | AARD 0.12 % | The Cp–Cv identity and our α, κ_T are right |
| 2 | Water + ethanol Hᴱ (88 points, 298–383 K) | RMSD 256 → **28 J/mol** after refit | VLE parameters are not caloric parameters |
| 2 | Water + ethanol Cpᴱ derived from measured Hᴱ(T) | RMSD 9.1 → **2.3 J/(mol K)** | The refit captures Cpᴱ (peak ≈ 14.5 J/(mol K)) |
| 2 | Water + methanol Cp (independent data) | ideal 7.5 %, VLE-NRTL 6.5 %, **Hᴱ-NRTL 2.7 %** | Gibbs–Helmholtz route to Cpᴱ works |
| 2 | Ethanol + toluene Cp, no parameters | 8.3 % (ideal fallback) | Honest limitation: missing binary data |
| 3 | Nanofluids, 69 DSC points (Al₂O₃, CuO, SiO₂) | **Model II 2.1 %**, Model I 18.5 % (max 63 %) | Mass weighting is right; Model I is wrong |
| 4 | Freezing points, water + ethanol ≤ 30 wt% | **0.84 K** mean (ideal + ΔCp) | The SLE model behind ice-slurry Cp works |
| 4 | Energy conservation through melting peaks | 10⁻⁸ % (PCM), 0.1 % (ice slurry, grid) | Latent heat is neither lost nor double-counted |
| 5 | Route A → Q → Route B, 7 systems | 6×10⁻⁷ % | The routes are mutually consistent |

What to say about the less good numbers:
* *Glycerol 2.4 % at 55 °C:* the literature values themselves scatter by 1.4 % at 25 °C.
* *Methanol–water at 280 K overshoots at high x:* 280 K is below the fitted range (298–378 K).
* *Freezing points above 40 wt% ethanol are off by up to 19 K:* the liquid is then deeply supercooled
  (−30 to −50 °C), far outside the range of the liquid-water Cp correlation, and ΔCp is not constant there.
  Ice slurries are used at 5–30 wt%, where the error is under 1 K.
* *NRTL activity coefficients make freezing points worse:* the parameters are fitted at 298–383 K and
  extrapolated to 220–270 K, so the engine uses ideal SLE by default.

## 7. Where the data came from (short version; full detail in `docs/DATA_SOURCES.md`)

* **Downloaded:** NIST SRD 69 fluid tables; the NIST ThermoML archive (189 MB, 11,923 files), from
  which we extracted 876 binary-mixture points with the DOI of each original paper.
* **Digitised:** O'Hanley et al. 2012 (open access) gives its nanofluid data only in figures. Our script
  reads the marker positions directly from the PDF's vector drawing commands, so there is no eyeballing.
* **Implemented a standard:** IAPWS R10-06 for ice, checked to 9 digits against the official test values.
* **Corrected during verification:** acetone–water NRTL (i and j swapped), glycerol Cp (4 % high),
  ice Cp (1.4 % low), CuO range, solid densities, NaCl solubility.
* **Open:** 39 low-impact flags (`docs/DATA_TO_VERIFY.md`); direct Cp data for water + ethanol (paywalled).

## 8. Limitations (say these before you are asked)

1. Cpᴱ needs binary parameters. Without them the engine falls back to an ideal mixture and warns (8 % error for ethanol + toluene).
2. No vapour–liquid equilibrium for mixtures. Boiling of a mixture inside the interval is flagged, not computed.
3. No solute partitioning between two liquid phases in emulsions.
4. Electrolytes: dissolved NaCl uses its solid Cp (Neumann–Kopp approximation) and ideal colligative behaviour.
5. Ideal volume of mixing for ρ, α and κ_T of mixtures; interfacial energy of nanoparticles is neglected.
6. Pure-component melting is smoothed over a finite width. The *interval-mean* Cp does not depend on that width.

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
    127 pytest tests: units, limits (φ→0, φ→1, xᵢ→1, Cpᴱ→0 at the pure ends), ideal reduction,
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
