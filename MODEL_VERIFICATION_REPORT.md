# Model verification report: Cp engine for multiphase fluids

Independent audit of the repository `cp-engine-multiphase` (commit `73ffa1b`), 5 November 2026.
No file of the scientific implementation (`cp_engine/`, `data/`) was modified during the audit. All new code is in
`audit/`:

| File | Purpose |
|---|---|
| `audit/independent_checks.py` | 79 independent checks; reference values computed **without** the engine's model functions |
| `audit/test_independent.py` | pytest wrapper; known defects marked `xfail(strict=True)` |
| `audit/results.json`, `audit/results.md`, `audit/run_output.txt` | raw results of the last run |

Re-run: `python audit/independent_checks.py` and `pytest audit/test_independent.py -q`.
At the time of the audit: `pytest` → 185 passed; audit suite → 79 passed, 3 xfailed (the known defects).

> **Status after the corrections (section L):** all defects G1–G7 corrected. `pytest` → **196 passed**
> (11 new regression tests in `tests/test_audit_fixes.py`); audit checks → **77/79** at the strict criteria, with the
> two remaining differences being the predicted O(h²) and reference-truncation sizes of section E; audit suite →
> **82 passed, 0 xfailed**. Sections A–K below are the audit as it was performed (before corrections); section L
> records what was changed and the evidence afterwards.

---

## A. Executive verdict

**Mostly correct, with minor issues** (category 2 of 5), within the scope stated in section H.

Evidence:
* The enthalpy formulation C_p,app = (∂H/∂T)_P with H assembled from Kirchhoff paths, phase amounts from phase
  equilibrium, NRTL excess enthalpy and EOS departure functions was derived independently (section C) and agrees with
  the code.
* 75 of 79 independent numerical checks pass at strict tolerances (most at 10⁻⁶ % or tighter). Three of the four
  failures are numerical truncation at tolerances I set deliberately tight (sizes predicted analytically, section E).
  One is a real defect.
* Two real defects were found (section G): a wrong guard that switches off the real-gas departure term for a
  **saturated** vapour and issues a false warning (small numerical effect), and missing validation of non-finite
  inputs (`T = NaN` returns "success" with cp = 0). Neither affects the primary scope of suspensions and PCM
  emulsions at normal conditions.
* Scientific approximations that a professor can legitimately attack: the melting width of pure substances is a
  modelling choice, not physics; water + alcohols use one NRTL parameter set for H/Cp and another for the
  boiling equilibrium; the vapour is ideal in the phase equilibrium but non-ideal (Peng–Robinson) in the enthalpy;
  supercooling and hysteresis of emulsified PCMs are not modelled.

What this report does **not** establish: agreement with experiment for emulsions or PCM slurries. No open
calorimetric data for them were found, so their apparent Cp is verified only for internal thermodynamic consistency.
This is marked **UNVERIFIED (experimentally)** in section F.

---

## B. Project architecture (computational workflow)

```
user dict
  │  inputs.normalise()        units → SI; composition (x, w, φ, moles, masses) → moles of each component
  │  validator.validate()      T > 0, P > 0, fractions in [0, 1], phases physically possible
  ▼
Solver (solver.py)             forward chaining over rules.py: every rule whose inputs are known fires once
  ├─ Route B/C rules           ΔH = Q + W (const P), ΔU = Q + W (const V), cp = ΔH/(mΔT), cp = ΔS/[m ln(T2/T1)],
  │                            cp = γ cv, cp = cv + Tvα²/κ_T
  └─ Route A rule ──► SystemModel (system.py)
         configure(T1, T2)    _vle_status(): bubble/dew points of miscible volatile clusters (models/vle.py)
                              plan per component: fixed solid/liquid/gas, dissolved, melt_smooth, sle_ideal,
                              sle_table, boil_smooth, vle
                              _group_liquids(): liquid phases by miscibility (first-fit; database flags / binodal)
         split(T)             phase amounts n_i^p(T): Gaussian melt fraction, SLE iteration, tabulated solubility,
                              isothermal flash (Rachford–Rice)
         enthalpy(T)          H = Σ_p Σ_i n_i^p H_i^p(T) + Σ_liquid n^p H^E,p + electrolyte term + gas departure
         cp_frozen(T)         ideal Σ n Cp_i + excess (NRTL Cp^E, electrolyte, gas Cp^R) at fixed amounts
         cp_point(T)          + latent = [H(T; n(T+h)) − H(T; n(T−h))]/(2h), h = min(1e-3 K, width/200)
         cp_mean(T1, T2)      total = [H(T2) − H(T1)]/(T2 − T1); ideal/excess parts by QUADPACK; latent by difference
  ▼
CpResult                      routes compared; measured routes preferred; breakdown; uncertainty (uncertainty.py);
                              partial molar Cp and mixing functions (models/partial_molar.py); warnings
```
Dependencies: NumPy, SciPy (`integrate.quad` = QUADPACK, `optimize.brentq`, `least_squares`), pandas.
Data: `data/components.json` (13 components), `data/binary_nrtl.json`, `data/model_uncertainty.json`.

---

## C. Mathematical verification (independent derivation vs implementation)

Notation: n_i^p moles of i in phase p; all quantities per calculation basis (default 1 kg); P fixed.

| # | Expression | Meaning / derivation | Assumptions | Code | Status |
|---|---|---|---|---|---|
| C1 | C_p,app = (∂H/∂T)_{P, n_total} | Definition for a closed system heated at constant P; "apparent" because phase amounts change with T | Equilibrium at every T (quasi-static) | `system.py:800` `cp_point` | **Correct.** P held at `self.P` throughout; total amounts fixed; derivative w.r.t. T only |
| C2 | H = Σ_p Σ_i n_i^p H_i^p(T,P) + Σ_L n^L H^E(T, x^L) + H^R_gas | H extensive, additive over phases (no interfacial term); mixing enthalpy only inside a phase | Interfaces neglected; ideal gas mixing in vapour except departure | `system.py:667` `enthalpy` | **Correct** for the stated assumptions |
| C3 | H_L(T) = ∫_{Tm}^{T} Cp_L; H_S = −ΔH_fus + ∫_{Tm}^{T} Cp_S; H_G = H_L(Tb) + ΔH_vap + ∫_{Tb}^{T} Cp_G^ig | Kirchhoff path from a common reference (liquid at Tm) | ΔH_fus at Tm, ΔH_vap at Tb (1 atm) | `models/pure.py` `enthalpy` | **Correct.** Verified: H_L − H_S = ΔH_fus + ∫ΔCp (D2); ΔH_vap recovered exactly at Tb (G5) |
| C4 | dH/dT = Σ n_i^p Cp_i^p + Σ n^L Cp^E + Σ (∂H/∂n_i^p) dn_i^p/dT | Chain rule: sensible + excess + latent | H differentiable in T and n | `system.py:775, 808` | **Correct.** Latent part as a central difference in the phase amounts only (O(h²)); matches FD of H to 10⁻⁷ % (section F, FD) |
| C5 | f(T) = ½[1 + erf((T − Tm)/(σ√2))], σ = width/4; dn_L/dT = n f′(T) | Melted fraction of a pure dispersed component as a Gaussian CDF of melting temperatures | Width is a **modelling choice** (impurities, size, rate), not a property | `phase_change.py:43` | **Correct implementation; physical basis is an approximation.** E1, E2, E3 pass |
| C6 | latent_melt = n f′(T)[ΔH_fus + ∫_{Tm}^{T} ΔCp dT] | From C3–C5 | — | via `cp_point` | **Correct.** D1 to 6×10⁻⁵ % at the peak; area = n[ΔH_fus + ½σ² dΔCp/dT]: predicted −0.25 J, observed −0.24 J (D2) |
| C7 | ln(x_k γ_k) = −(ΔH_fus/R)(1/T − 1/Tm) + (ΔCp/R)(Tm/T − 1 − ln(Tm/T)) | Solid–liquid equilibrium, pure solid k ↔ solution; ΔCp term from Kirchhoff | Pure solid crystallises (no solid solutions) | `phase_change.py:58`, `system.py:527` | **Correct** (standard; e.g. Prausnitz et al., *Molecular Thermodynamics of Fluid-Phase Equilibria*). Liquidus reproduced to 10⁻⁸ K (validation case 4) |
| C8 | G^E/RT = Σ_i x_i Σ_j τ_ji G_ji x_j / Σ_k G_ki x_k; H^E = −RT² ∂(G^E/RT)/∂T; Cp^E = ∂H^E/∂T | NRTL (Renon & Prausnitz 1968) + Gibbs–Helmholtz | τ(T) form valid in the fitted T range | `nrtl.py:114–153` | **Correct.** ln γ from ∂(nG^E)/∂n_i agrees with analytic to 10⁻⁷; Gibbs–Duhem < 10⁻⁸ (tests/test_eos_vle.py) |
| C9 | c_p = Σ_p w_p c_p,p ⇔ Model II | Mass-weighting is C2 for non-interacting phases | No interfacial enthalpy | `heterogeneous.py:50`; Route A | **Correct.** C1 matches an independent calculation with NIST water density to 0.019 % |
| C10 | (∂Cp/∂P)_T = −T(∂²v/∂T²)_P; Tait v(T,P) | Maxwell relation; liquid departure from the 1-atm reference | Tait C = 0.0894; κ_T(T) | `database.py:266, 280` | **Correct.** 8-point Gauss–Legendre vs closed form: 2×10⁻¹² % (E4, E5) |
| C11 | y_i P = x_i γ_i P_i^sat; Σ K_i x_i = 1; Rachford–Rice | γ–φ VLE with ideal vapour | Low P; Poynting ≈ 1 | `vle.py:61–123` | **Correct.** Flash satisfies mass balance to 10⁻¹⁰ and K = y/x to 10⁻⁸; bubble T 0.26/0.46 K vs 409 measured points |
| C12 | H^R/RT = Z − 1 + (d ln a/d ln T − 1) q I; Cp^R = ∂H^R/∂T | Departure function of the generic cubic (Smith, Van Ness & Abbott) | Mean-field; k_ij = 0 | `eos.py:146–157` | **Correct.** Agrees with ∫[V − T(∂V/∂T)_P]dP to 10⁻⁵ % (G1) |
| C13 | H_gas = H_ig-path + n[H^R(T,P,y) − Σ y_i H_i^R(Tb,i, 1 atm)] | Real saturated vapour at Tb → ideal gas → heat → real gas at (T,P) | — | `system.py:686` | **Correct algebra** (G4, G5). **Defect in the guard** at line 717 (G1, section G) |
| C14 | units | n in mol, H in J, Cp in J/K; per-mass values = /m (kg) | — | throughout | **Consistent.** B1/B2 identical on both bases to 10⁻¹² |

Basis consistency: component amounts are always moles internally (`inputs.py`); mass fractions are derived with the
same molar masses; B1 (molar) and B2 (mass) agree exactly.

---

## D. Thermodynamic verification

**Enthalpy and apparent Cp.** H is built as a state function (C2, C3) and C_p,app is its total temperature derivative
at constant P and fixed total amounts. The code does not average heat capacities by any weighting that bypasses H:
for a suspension the result equals the mass-weighted rule because that is what dH/dT gives for non-interacting phases
(C9). Energy conservation ∫C_p,app dT = ΔH holds to 10⁻⁷–10⁻¹¹ through melting (D3, F4), through ice crystallisation
in a solution (F4, 255–285 K), through overlapping transitions of two PCMs (F5, 6×10⁻¹² %) and through the boiling
range of water + ethanol (validation case 6, ≤ 0.03 %).

**Sensible vs apparent heat capacity.** The breakdown keeps them separate: `ideal` and `excess` are the frozen
(fixed-amount) heat capacity; `latent` is the term caused only by the change of phase amounts (C4). Inside a
transition the latent term dominates (PCM slurry at Tm: 42.45 kJ/(kg K), of which 3.76 kJ/(kg K) is sensible; ice slurry
at 268 K: 49.71 kJ/(kg K), of which 3.84 is sensible + excess).

**Phase transitions.** Three mechanisms, each verified:
1. Pure dispersed component (PCM droplets, ice in pure water): Gaussian smoothed melting (C5, C6). The point value
   at Tm scales as 1/width (sensitivity: width 1 K → +91 %, 4 K → −46 %); interval means spanning the step are
   width-independent (≤ 8×10⁻¹¹). The engine warns and reports the width-independent quantities. **This is the
   weakest physical element:** the width is not derived from thermodynamics.
2. Component crystallising from a solution (ice slurry, PCM in a miscible oil): SLE (C7) with NRTL γ where fitted to
   freezing points. Freezing points ± 2.5 K held-out (validation case 2b).
3. Boiling of a miscible liquid mixture: γ–φ VLE (C11) between bubble and dew points.

**Latent heat.** Enters only through ∂H/∂n × dn/dT, with ∂H/∂n = H_L(T) − H_S(T) (or H_G − H_L) evaluated at the actual
temperature, i.e. ΔH_fus(T) = ΔH_fus + ∫ΔCp (Kirchhoff). It is therefore counted exactly once (no separate "add
ΔH_fus" term exists anywhere). F5b shows the expected consequence: crystallisation below Tm releases about 1 % less
than n ΔH_fus(Tm), consistent with ΔCp × (Tm − T).

**Departure function.** Used only for **gas phases** (C12, C13). It is necessary for vapours at elevated pressure
(steam at 10 MPa: ideal gas 17.5 % low; with the PR departure 3.4 % low). It is **not needed and not used** for the
primary scope (suspensions, emulsions, slurries: condensed phases), where the pressure effect is the Tait correction
(C10), itself a liquid-phase departure from the 1-atm reference. No double counting: the Tait correction is
identically zero for gases (G3); H(T, 10 MPa) − H(T, 1 atm) of steam equals n[H^R(P) − H^R(P₀)] exactly (G4); the
reference offset recovers ΔH_vap at Tb exactly (G5). Remaining issue: the VLE uses an ideal vapour (φ = 1) while the
vapour enthalpy carries the PR departure, so the phase equilibrium and the enthalpy use different vapour models.
At 1 atm H^R is of order 10² J/mol against ΔH_vap of order 4×10⁴ J/mol (≈ 0.3 %).

---

## E. Numerical method verification

### E.1 Where "Gaussian" appears in this model (and where it does not)
The audit brief assumed that Gaussian integration over a temperature distribution is central to the model. The code
uses three different things, each audited:

| Use | What is integrated / represented | Method | Check | Result |
|---|---|---|---|---|
| Gaussian **distribution** of melting temperatures | melted fraction f(T) | closed form (erf); **no numerical integration** | E1: f vs `scipy.stats.norm.cdf` | max diff 1.1×10⁻¹⁶ |
| | normalisation ∫f′dT | — | E2: QUADPACK | 1 − 1.2×10⁻¹² |
| | "95 % within ±width/2" | ±2σ | E3 | 0.95450 ✔ |
| Gauss–**Legendre** quadrature (`database._p_integral`, n = 8) | ∫_{P₀}^{P}(v − T∂v/∂T)dp and −T∫∂²v/∂T² dp (liquid pressure correction) | nodes/weights `numpy.polynomial.legendre.leggauss`; Jacobian (P − P₀)/2 included | E4: vs closed-form Tait integral, water, 298 K, 100 MPa | n = 2: 1.3×10⁻³ %; n = 4: 1.3×10⁻⁷ %; **n = 8: 2×10⁻¹² %**; n = 16: 2.7×10⁻¹² % (converged) |
| Adaptive Gauss–**Kronrod** (QUADPACK, Piessens et al. 1983) | Kirchhoff ∫Cp dT | `scipy.integrate.quad` | E7: vs closed-form antiderivatives (DIPPR-100, Shomate, DIPPR-107) | ≤ 3×10⁻¹⁴ % |
| | interval-mean ideal/excess parts | `quad` with break points at transitions | E8: vs 40-point Gauss–Legendre | 5×10⁻¹⁴ % |

Integration limits: P₀ = 101 325 Pa → P (the Cp correlations are 1-atm data, so the lower limit is the reference
pressure, not zero: correct); T₁ → T₂ for interval means; T_ref → T for Kirchhoff. Temperature derivatives inside the
pressure integral use a 0.5 K central stencil: truncation 1.9×10⁻³ % against a 0.05 K stencil (E6). Acceptable;
documented.

### E.2 Finite differences
The direct output `cp_point` was compared with [H(T + ΔT) − H(T − ΔT)]/(2ΔT) from the engine's own H(T), for seven
systems and ΔT = 1, 0.1, 0.01, 0.001 K (section F, "FD" rows). Smooth regions agree to 10⁻⁸–10⁻⁹ % at ΔT ≤ 0.01 K, with
second-order convergence (error ∝ ΔT²). Inside a transition, the FD reference itself has truncation error
ΔT²H‴/(6H′). At the centre of the Gaussian peak this is ΔT²/(6σ²) = 6.7×10⁻³ % for ΔT = 0.01 K and σ = 0.5 K; the
observed difference is 6.0×10⁻³ %, and at ΔT = 0.001 K it falls to 2×10⁻⁹ %. At a bubble or dew point C_p,app is
discontinuous (a kink in H); there a central difference straddling the boundary returns an average of the two
one-sided values and is meaningless. The engine does not difference across such points (interval means use H
directly, and quadrature break points are placed at Tm, T_bubble, T_dew).

### E.3 Stability
* NRTL Cp^E uses a second derivative with h = 10⁻⁴T. The FD of H at ΔT = 0.001 K shows round-off of 10⁻⁵ % (expected
  ε/h²); the model's own derivative is stable (checked against sympy in `tests/test_models.py`).
* The latent central difference uses h = 10⁻³ K. At T = Tm − 2σ the error is 1.2×10⁻⁴ % (D1 at 300.3 K, which fails
  my deliberately tight 10⁻⁴ % criterion). Negligible.
* The interval-mean breakdown sets `latent = 0` when it is below 10⁻⁶ of the total and there is no phase change, so the
  breakdown total equals the quadrature result rather than exactly ΔH/ΔT (E9: 1.4×10⁻⁹ relative). Negligible.

---

## F. Independent validation results

All values from `audit/results.json` (last run). "Reference from" says how each reference was obtained. Criteria were
fixed before running.

| Test | Inputs | Expected | Reference from | Reference | Model | Rel. err. % | Pass |
|---|---|---|---|---|---|---|---|
| A1 | water liquid, 298.15 K | NIST value | NIST SRD 69 table | 4180.60 J/(kg K) | 4182.95 | 0.056 | ✔ (≤ 0.5) |
| A2 | same | DIPPR-100 re-implemented | raw JSON coefficients | 4184.524 | 4184.524 | 0 | ✔ |
| A3 | Al₂O₃ 400 K / SiO₂ 350 K | Shomate re-implemented | NIST Shomate form | 96.754 / 49.470 J/(mol K) | same | 0 | ✔ |
| A4 | ethanol ideal gas 500 K | DIPPR-107 re-implemented | raw JSON | 95.926 | 95.926 | 0 | ✔ |
| B1 | ethanol/water x = 0.3, ideal, 298.15 K | Σ x_i Cp_i | independent | 86.4717 J/(mol K) | 86.4717 | 0 | ✔ |
| B2 | same | Σ w_i cp_i | independent | 3271.575 J/(kg K) | 3271.575 | 0 | ✔ |
| C1 | Al₂O₃ φ = 0.04 in water | Model II with NIST ρ_water | NIST + Shomate | 3697.07 J/(kg K) | 3696.38 | 0.019 | ✔ (≤ 0.05) |
| C2 | same, ΔH 290→330 K (basis 1 m³ = 1116.6 kg when φ is given) | Σ n ∫Cp | closed-form antiderivatives | 1.651566×10⁸ J | same | 2×10⁻¹⁴ | ✔ |
| D1 | 20 wt% C18 in water, 297–305 K (6 temperatures) | sensible + n f′[ΔH_fus + ∫ΔCp] | scipy.stats.norm + closed form | e.g. 42 453.32 J/K (basis 1 kg) at 301.3 K | 42 453.29 | ≤ 1.2×10⁻⁴ | 5 ✔, 1 ✘ (numerical, E.3) |
| D2 | same | latent area | analytic | n ΔH_fus = 48 491.94 J | 48 491.70 | 5×10⁻⁴ | ✔ (deficit matches ½σ² dΔCp/dT) |
| D3 | same, 290→312 K | sharp-melting path | closed form | 131 142.3 J | same | 7×10⁻¹⁴ | ✔ |
| E1–E9 | see E.1 | | | | | | 8 ✔, 1 ✘ (E9, numerical) |
| F1 | water, EtOH, toluene, SiO₂, C18; 121 T, 255–315 K | Σ_p n_i^p = n_i | exact | 0 | 2×10⁻¹⁶ | — | ✔ |
| F2 | same | no negative amounts | exact | 0 | 0 | — | ✔ |
| F3 | same | allocation rule respected | exact | 0 violations | 0 | — | ✔ (rule itself is a limitation) |
| F4 | same, 290–312 K / water-EtOH-toluene 255–285 K | ∫cp_app dT = ΔH | quadrature | 69 736.12 / 196 213.48 J | same | 10⁻⁷ / 5×10⁻⁸ | ✔ |
| F4a | C18 slurry 268–312 K | computable | — | — | 10 701.8 J/(kg K) | — | ✔ |
| F5 | water + C18 (301.3 K) + user PCM-B (300.3 K), overlapping | path-independent ΔH | closed form | 144 510.796 J | 144 510.796 | 6×10⁻¹² | ✔ |
| F6 | same | phase mass fractions sum to 1 | exact | 1 | 1 | 0 | ✔ |
| FD | 7 systems × 4 ΔT | FD → cp_point | engine H differenced | see results.md | | ≤ 10⁻⁵ at ΔT ≤ 0.01 K except PCM peak | 27 ✔, 1 ✘ (reference truncation, E.2) |
| G1 | PR, steam 773 K 10 MPa / n-hexane 550 K 2 MPa | two routes to H^R | ∫[V − T∂V/∂T]dP | −2317.10 / −3399.44 J/mol | same | 9×10⁻⁶ / 4×10⁻⁵ | ✔ |
| G2 | same | Cp^R = ∂H^R/∂T | FD | 6.5846 / 16.4535 | 6.5847 / 16.4537 | 3×10⁻⁴ / 10⁻³ | ✔ |
| G3–G5 | see section D | | | | | | ✔ |
| G6 | ethanol/water x = 0.1, 362 K | Cp^R of saturated vapour ≠ 0 | EOS directly | 0.0560 J/K | **0** | 100 | ✘ **defect** |

**Totals: 75/79 pass.** The failures are D1 at 300.3 K, E9 and FD at the PCM peak (all numerical; sizes predicted)
and G6 (a real defect).

**Edge cases** (from `audit/results.json` → `edge_cases`):

| Case | Outcome |
|---|---|
| x → 0, x → 1 (10⁻¹²), φ = 0, φ = 1 | finite, correct limits (4182.95 = pure water; 532.6 = pure CuO) |
| T exactly Tm | value with explicit width-dependence warning |
| T = 0, P < 0, unknown component, Σw = 1.2, T1 = T2, liquid water at 300 °C (beyond data) | clean error messages |
| **T = NaN** | **"success", cp = 0.0 J/(kg K)**: defect (section G) |
| P = ∞ | error, but the message ("min() iterable argument is empty") is not informative |
| T exactly at the bubble point | value with the false "below its dew point" warning (defect G1) |

**Experimental validation** (from `validation/`, not repeated here): over 1400 published data points, AARD 0.33 % for
pure liquids, 2.1 % for nanofluids (Model II), 2.7 % for independent water + methanol Cp, 0.26/0.46 K bubble points.
**UNVERIFIED experimentally:** apparent Cp of emulsions and PCM slurries (no open DSC data found; candidate source
listed in `docs/DATA_TO_VERIFY.md`).

---

## G. Bugs and scientific errors

| ID | Severity | Location | Description | Evidence | Consequence |
|---|---|---|---|---|---|
| **G1** | Medium (correctness of a warning; small numerical effect) | `cp_engine/system.py:686–722`, `SystemModel._gas_residual`, guard `if s > 1.02` (line 717) | The test for a supersaturated vapour uses s = Σ y_i P/P_i^sat **without γ_i**. A vapour in equilibrium with a non-ideal liquid has Σ y_i P/(γ_i P_i^sat) = 1 but s > 1 when γ > 1. | water + ethanol x = 0.1: s = 1.24 at 360 K, 1.05 at 370 K; Cp^R set to 0; warning "would condense" (G6; edge case at the bubble point) | False warning shown in demos of boiling mixtures; the departure term is dropped inside the VLE region (0.056 J/K out of about 3700 J/K at 362 K) and switches on where s falls below 1.02, a small step in H |
| **G2** | Medium (robustness) | `cp_engine/validator.py:34, 36, 47, 58` | Positivity checks `k[sym] <= 0` are False for NaN, so NaN passes; ∞ is also not rejected | `T = NaN` → success, cp = 0.0; `P = ∞` → cryptic error | A silent wrong number for invalid input |
| G3 | Low (consistency) | `models/vle.py` (K with φ = 1) vs `system.py:686` (PR departure) | Vapour ideal in the phase equilibrium, non-ideal in the enthalpy | section D | ≈ 0.3 % of the latent heat of vaporisation at 1 atm; larger at high P |
| G4 | Low (design, documented) | `system.py` (`vle_parameters="VLE"` vs `nrtl_parameters="HE"`) | Phase split from the VLE-fitted NRTL set, H and Cp^E from the H^E-fitted set: Gibbs–Helmholtz link between them broken | joint fit doubled the Cp error (handbook §5.11); sensitivity: switching set changes Cp of water + EtOH by −15 % | Defensible as a measured trade-off, but must be stated |
| G5 | Low (numerical) | `system.py:862` | Breakdown `latent` zeroed when small, so breakdown total ≠ ΔH/ΔT by 10⁻⁹ | E9 | Negligible |
| G6 | Low (documentation) | `phase_change.py` module docstring | Claims "∫ latent cp dT = ΔH_fus exactly"; exact only if ΔCp is constant (true area n[ΔH_fus + ½σ² dΔCp/dT]) | D2 | Wording only |
| G7 | Low (references) | `manuscript/references.bib` | Pak1998 and Xuan2000 marked "% VERIFY" | Crossref confirms both (doi:10.1080/08916159808946559; doi:10.1016/S0017-9310(99)00369-5) | Flags can be cleared |

No sign errors, dimensional errors, indexing or broadcasting errors were found in the audited paths (A–G, FD).
Mass balances are exact (F1, flash test).

---

## H. Limitations and unsupported conditions

**Supported and verified (within the data ranges of the components):**
* N components (tested up to 5 at once); phases: one particulate solid phase per solid component, liquid phases
  formed by miscibility, at most one vapour phase.
* Component amounts that change phase with T through Gaussian melting, SLE (several crystallising components;
  approximate below a eutectic, with a warning), tabulated salt solubility, or VLE of one liquid phase.
* Overlapping transitions (F5).

**Not supported (the engine warns or refuses; it does not silently compute):**
* Partitioning of a solute between two liquid phases: a component is assigned to one liquid phase by a first-fit rule
  (F3).
* Three-phase vapour–liquid–liquid equilibrium (immiscible volatile liquids): detected, warned, not modelled.
* Solid solutions, polymorphism, glass transitions.
* **Supercooling and heating/cooling hysteresis of emulsified PCMs**: crystallisation in small droplets needs
  supercooling, which equilibrium thermodynamics cannot describe. The model is an equilibrium (heating-like) model.
  Not verified against data in this project.
* Interfacial enthalpy of nanoparticles or droplets.
* Boiling-point elevation by dissolved salts (salts are not part of the VLE cluster).
* Electrolytes other than NaCl, and NaCl away from 298 K (warned).
* Temperatures more than 30 K outside any correlation's range (refused). For example, octadecane dissolved in toluene
  below 270 K: the ideal-solubility SLE leaves some liquid octadecane there, but its liquid-Cp data start at 300 K.

**Model choices a reviewer should know:** the melting width of pure dispersed components (default 2 K); NRTL
parameter sets per property (G4); ideal vapour in the VLE (G3); k_ij = 0 in the EOS; γ = 1 for pairs without NRTL
parameters (warned; for ethanol + toluene this was 8 % before fitting).

---

## I. Professor's evaluation

*As an examiner for an undergraduate CL 202 project.*

Strengths that would be rewarded: the formulation starts from H and differentiates it, instead of averaging heat
capacities, so latent heat cannot be double counted. Energy is conserved to machine precision through every kind of
transition. Every fitted model has held-out tests, and the lecture equations (EOS examples) are reproduced in tests.
The data provenance is unusually careful.

Weaknesses that would cost marks:
1. **The two defects in section G.** Both are small, but an examiner who types NaN or opens a boiling-mixture result
   sees them immediately.
2. **The melting width.** A student must be able to say plainly that it is a regularisation of a delta function, and
   that only interval means and the latent heat are physical.
3. **No experimental check for the emulsion and PCM cases**, which are the core of the stated scope.
4. **Different parameter sets and vapour models** for equilibrium and enthalpy (G3, G4): defensible, but it must be
   volunteered, not discovered.

Mark (out of 20): **17.5 as it stands; about 18.5 after G1–G3 are corrected.** Full marks would need experimental
validation of a PCM emulsion and a model of supercooling (out of scope for the deadline).

---

## J. Viva questions (specific to this implementation)

Each entry: (1) what is tested, (2) answer, (3) equation, (4) where in the code, (5) follow-up.

1. **What exactly is your apparent heat capacity?** (1) Definition. (2) The total T-derivative of the system enthalpy at
   constant P and fixed overall composition, including the effect of phase amounts changing. (3) C_p,app = (∂H/∂T)_{P,n}.
   (4) `system.py:800` `cp_point`. (5) Why is it not a material property?
2. **Why compute Cp from H instead of a weighted average of Cp's?** (1) Understanding of latent terms. (2) A weighted
   average misses dn/dT terms and the excess enthalpy; differentiating H includes both automatically and guarantees
   ∫Cp dT = ΔH. (3) C4. (4) `enthalpy`, `cp_frozen`. (5) Show a case where they differ: ice slurry at 268 K, 49.7 kJ/(kg K)
   against a sensible ≈ 3.8.
3. **How does latent heat enter?** (1) No double counting. (2) Only via ∂H/∂n × dn/dT; ∂H/∂n = H_L − H_S at the actual T.
   (3) C6. (4) `cp_point` line 808. (5) Why is the latent released below Tm smaller than ΔH_fus(Tm)? (Kirchhoff, F5b.)
4. **Why a Gaussian for melting?** (1) Physical honesty. (2) A pure substance melts at one T (a delta function in
   C_p,app); real samples melt over a range (impurities, droplet sizes, heating rate). The Gaussian CDF is a smooth,
   normalised regularisation with σ = width/4. (3) C5. (4) `phase_change.py:43`. (5) Which outputs do not depend on the
   width? (Latent heat; interval means: sensitivity ≤ 10⁻¹⁰.)
5. **What happens if the width doubles?** (1) Sensitivity. (2) Peak height roughly halves (−46 % for 4 K vs 2 K), area
   unchanged. (3) f′ ∝ 1/σ. (4) audit sensitivity table. (5) What width would you choose for DSC data?
6. **Is your Gaussian normalised, and over what limits?** (1) Distribution handling. (2) Yes, over (−∞, ∞), analytically
   via erf; numerically 1 − 1.2×10⁻¹². No integral of the distribution is computed numerically. (3) ∫f′dT = 1.
   (4) `melt_fraction`. (5) Why not truncate at ±3σ?
7. **Where do you use Gaussian quadrature, and why?** (1) Numerics. (2) 8-point Gauss–Legendre for the pressure
   integral of the Tait correction; smooth integrand, so 8 points give 10⁻¹² relative error. QUADPACK (adaptive
   Gauss–Kronrod) for ∫Cp dT and interval means. (3) ∫_{P₀}^{P} f dp ≈ (P − P₀)/2 Σ w_k f(p_k). (4) `database.py:259`.
   (5) Where is the Jacobian? ((P − P₀)/2.)
8. **Why does the pressure integral start at 1 atm, not 0?** (1) Limits. (2) The Cp correlations are 1-atm data; the
   correction is relative to that reference. (3) Cp(T,P) = Cp(T,P₀) − T∫_{P₀}^{P}∂²v/∂T² dp. (4) `cp_pressure_correction`.
   (5) Why is the gas correction zero there?
9. **How do you know the quadrature is converged?** (1) Verification. (2) n = 2, 4, 8, 16 against the closed-form Tait
   integral: 10⁻³, 10⁻⁷, 2×10⁻¹², 3×10⁻¹² %. (4) audit E4. (5) What limits accuracy now? (The 0.5 K T-stencil, 2×10⁻³ %.)
10. **What is a departure function and why do you have one?** (1) Residual properties. (2) H − H^ig at the same T, P; needed
    only for gases at pressure; condensed phases use the Tait correction instead. (3) C12. (4) `eos.py:146`. (5) Is it
    needed for your PCM emulsion? (No.)
11. **What is your reference state for the vapour enthalpy?** (1) Consistency. (2) Liquid at Tm → liquid at Tb →
    measured ΔH_vap gives the real saturated vapour → subtract H^R(Tb, 1 atm) → ideal gas → heat → add H^R(T,P).
    (3) C13. (4) `_gas_residual`. (5) Show that ΔH_vap is recovered at Tb (G5: exact).
12. **Does the departure term double count with the Tait correction?** (1) Double counting. (2) No: Tait is zero for gases
    (G3); the steam H difference between 10 MPa and 1 atm equals n[H^R(P) − H^R(P₀)] exactly (G4). (5) And for liquids?
13. **Your vapour is ideal in the VLE but not in the enthalpy. Defend that.** (1) Consistency (G3). (2) At 1 atm the
    effect is about 0.3 % of ΔH_vap; the `lewis_randall` option adds virial φ; a fully consistent version would use the
    same EOS in both. (5) What would change at 10 bar?
14. **How do you guarantee mass conservation?** (1) Bookkeeping. (2) Every phase split assigns n_i to phases by
    construction (melt: f n and (1 − f) n; SLE: n_L and n − n_L; flash: n_V = V n y, n_L = n − n_V); checked to 2×10⁻¹⁶
    over 121 temperatures. (4) `split`. (5) Can an amount go negative? (Checked: no; `min()` guards in the flash.)
15. **Does it really handle N components and M phases?** (1) Generality. (2) N arbitrary; phases = one per solid component
    + liquid phases by miscibility + ≤ 1 vapour; compositions inside a phase vary by SLE/VLE. Not general: no
    partitioning between liquid phases, no three-phase boiling. (4) `_group_liquids`, `split`. (5) Ethanol in
    water + toluene: where does it go and why?
16. **Two transitions overlap. What happens?** (1) Overlapping transitions. (2) Each component has its own f(T) or SLE;
    contributions add; path-independent ΔH reproduced to 6×10⁻¹² % (F5). (5) What if they form a eutectic? (Warned:
    approximate below a eutectic.)
17. **Phase fraction → 0 or 1?** (1) Limits. (2) Limits recover the pure phases (x → 0: 4182.95 = pure water; φ = 1:
    pure CuO). Gaussian tails mean a tiny amount of the other phase always exists (≈ 10⁻¹⁰⁰), which is harmless.
    (5) Could tails cause out-of-range evaluations?
18. **Why does Model I fail?** (1) Basis consistency. (2) cp is per unit mass; weighting by volume fraction is
    dimensionally inconsistent unless ρ_s = ρ_f. (3) Model II = Σ w cp. (4) `heterogeneous.py:45–50`. (5) For CuO, by
    how much? (Up to 63 % vs data.)
19. **Where does Cp^E come from?** (1) Gibbs–Helmholtz. (2) H^E = −RT²∂(G^E/RT)/∂T from NRTL; Cp^E = ∂H^E/∂T,
    numerically with h = 10⁻⁴T. (4) `nrtl.py:148–157`. (5) Why are VLE parameters bad for Cp^E? (They fix G^E, not its T
    derivative: 6–9× too small.)
20. **Why two NRTL sets for water + alcohols?** (1) Consistency (G4). (2) A joint fit doubled the Cp error; each property
    uses the set fitted to it. (5) What would a consistent model need? (More flexible G^E(T, x) or an EOS with
    association.)
21. **How is the ice fraction in a slurry found?** (1) SLE. (2) Equal chemical potential of ice and water in solution,
    with the ΔCp term and NRTL γ. (3) C7. (4) `_solve_sle`. (5) Why include ΔCp? (Consistency with the Kirchhoff
    enthalpies.)
22. **How do you know cp_point is right inside a peak?** (1) Verification. (2) Analytic Gaussian expression (D1, 6×10⁻⁵ %)
    and FD of H converging as ΔT² (E.2). (5) Why does the FD fail at ΔT = 0.01 K? (Reference truncation ΔT²/6σ².)
23. **Can you difference across the bubble point?** (1) Discontinuities. (2) No: C_p,app jumps there; the engine
    integrates with break points and evaluates means from H. (5) Is the jump physical? (Yes: evaporation starts.)
24. **What is the uncertainty of a result?** (1) Error analysis. (2) √Σ(deviation_k × share_k)² from the project's
    validation deviations. Approximation: mean absolute deviations treated as standard uncertainties, independent.
    (4) `uncertainty.py`. (5) Why not for a point inside a melting peak? (The width dominates.)
25. **What would make the PCM-emulsion results wrong in practice?** (1) Physical limits. (2) Supercooling and hysteresis,
    droplet-size effects, surfactant enthalpy: none modelled. (5) How would you test it? (DSC heating and cooling scans.)
26. **What happens with T = NaN?** (1) Robustness. (2) Today: a false "success" (defect G2); fix: reject non-finite
    inputs in the validator. (5) What other inputs need finiteness checks?
27. **Which input has the largest effect on the result?** (1) Sensitivity. (2) Inside a transition: melt width (±50–90 %)
    and the temperature (−1.8 % per 0.1 K at Tm). In a solution: the choice of NRTL set (−15 %). Sensible regions:
    composition and pure Cp almost linearly (Cp_water +1 % → cp +0.67 % at 80 wt% water). (5) Which matters for an
    energy-storage design? (Interval means, which are width-independent.)
28. **Is your H continuous through melting?** (1) Smoothness. (2) Yes for Gaussian melting (C∞); at bubble and dew points
    H is continuous with a kink. (5) What would a sharp model do to an implicit time integrator?

---

## K. Final verdict

**Mostly correct, with minor issues.**

* Thermodynamic formulation: correct; independently derived; energy conserved to machine precision through every
  transition type tested.
* Numerical methods: correct and converged (Gauss–Legendre, QUADPACK, finite differences), with predicted truncation
  sizes.
* Two real defects (G1, G2) and three consistency or documentation issues (G3, G4, G6); none changes the
  primary-scope results at normal conditions by more than about 0.3 %.
* Experimental support is strong for pure liquids, solutions, nanofluids, freezing points and boiling points, and
  **absent for emulsions and PCM slurries** (UNVERIFIED experimentally).

Not claimed: validity outside the data ranges, for unsupported phase configurations (section H), or for non-equilibrium
behaviour such as supercooling.

---

## Prioritised correction plan (as proposed; see section L for what was done)

| Priority | Item | Change | Risk | Expected effect |
|---|---|---|---|---|
| 1 | G1 | In `_gas_residual`, apply the supersaturation guard only to user-declared/fixed gases, or use Σ y_i P/(γ_i P_i^sat) with γ from the coexisting liquid; skip it for vapour produced by a flash (saturated by construction) | Low | Removes the false warning; restores Cp^R in the VLE region; removes the small H step |
| 2 | G2 | In `validator.py`, reject non-finite numbers (`math.isfinite`) for T, T1, T2, dT, P, amounts and energies | Very low | NaN/∞ give clear errors |
| 3 | G3 | Use φ̂_i from the selected gas EOS in the VLE K-values (γ–φ with real vapour), or document the 0.3 % inconsistency | Medium (re-validate bubble points) | Consistent vapour model |
| 4 | G6, G7 | Correct the docstring; clear the "% VERIFY" flags on Pak1998, Xuan2000 | None | Documentation |
| 5 | G5 | Report the breakdown with the total fixed at ΔH/ΔT and the quadrature mismatch as a diagnostic | Low | Cosmetic |
| 6 | Tails | Optionally truncate the Gaussian at ±8σ (exact 0/1 outside) to avoid evaluating correlations where the phase is absent | Low | Fewer extrapolation warnings; results unchanged (< 10⁻¹⁰) |
| 7 | Scope | Add a supercooling or hysteresis option for PCM emulsions; obtain DSC data (`docs/DATA_TO_VERIFY.md`) | High (new physics) | Addresses the main experimental gap |

---

## L. Corrections applied after the audit (authorised by the project owner)

| ID | Change | Files | Evidence after the change |
|---|---|---|---|
| G1 | The supersaturation guard of the gas departure term no longer fires for vapour produced by the phase equilibrium (flash, pure boiling), which is saturated by construction, nor for a mixture already known to be above its dew point. A user-declared gas mixture is compared with its own dew point from the γ–φ equilibrium (**with** γ); a pure gas with P_sat. | `cp_engine/system.py` (`_supersaturated`, `_gas_residual`) | G6: Cp^R of the saturated vapour 0.0552 J/K = reference (was 0); no false warning at 360–370 K or at the exact bubble point; a declared water vapour at 300 K, 1 bar is still flagged |
| G2 | Every number in the input must be finite (T, T1, T2, P, energies, properties, component amounts). | `cp_engine/validator.py` | `T = NaN`, `P = ∞`, `T1 = NaN`, NaN fraction → clear `InputError` "… is not a finite number" |
| G3 | One vapour model for equilibrium and enthalpy: new default `vapour_model="eos"`, K_i = γ_i P_i^sat φ_i^sat/(φ̂_i P) with φ from the same equation of state (`gas_eos`, default PR) that gives the vapour departure enthalpy. Root bracketed with the ideal-vapour K, then refined; φ^sat cached per T (11× faster than a naive implementation). Falls back to the ideal vapour, with a note, if critical constants are missing or no vapour root exists. | `cp_engine/models/vle.py`, `system.py`, `inputs.py` | Bubble points 0.257 K (ethanol, was 0.257) and 0.471 K (methanol, was 0.465), max errors lower (3.13/1.75 K vs 3.17/1.92 K); ethanol + toluene 0.60 K (was 0.63). Energy conservation through boiling 5×10⁻⁶ / 2×10⁻⁵ % (was 0.001 / 0.025 %, because the departure term no longer switches on and off) |
| G4 | Not changed (design choice): separate NRTL sets for caloric properties and the phase split of water + alcohols, stated in the handbook (§5.11, Q38) and the manuscript. | — | Measured trade-off: a joint fit doubles the Cp error |
| G5 | Interval-mean breakdown: without phase change the quadrature remainder is absorbed into the excess part, so the breakdown total equals ΔH/ΔT exactly; the remainder is reported as `quadrature_mismatch_J_per_K`. | `cp_engine/system.py` (`cp_mean`) | E9 exact (was 1.4×10⁻⁹); Route B vs Route A now 0 for all 7 systems of validation case 5 |
| G6 | Docstring corrected: latent area = n[ΔH_fus + ½σ² dΔCp/dT], exact only for constant ΔCp. | `cp_engine/models/phase_change.py` | — |
| G7 | Bibliography: DOIs added and flags cleared for Pak1998, Xuan2000, Kell1975, Renon1968, Dymond1988 (all checked on Crossref). The three book editions (Lange, Smith–Van Ness, Prausnitz) cannot be checked online and remain flagged. | `manuscript/references.bib` | — |
| Tails | Gaussian melting step cut at ±8σ (f = 0 or 1 exactly; neglected tail mass 1.2×10⁻¹⁵), so correlations of an absent phase are not evaluated far outside their ranges. | `cp_engine/models/phase_change.py` | PCM energy conservation 10⁻⁸ (test); D1–D3 unchanged |

Not done (out of scope for the deadline, listed honestly): supercooling and hysteresis of emulsified PCMs, and
experimental validation of emulsions and PCM slurries (no open data). These remain **UNVERIFIED experimentally**.

**Revised verdict after the corrections:** *Scientifically and computationally verified within its stated scope*
(section H) for thermodynamic consistency, numerics and the systems with experimental data (pure liquids and
vapours, liquid solutions, nanofluids, freezing points, boiling points); the PCM-emulsion and slurry results are
internally verified (energy, mass, phase equilibrium) but not experimentally validated. Professor's evaluation
(section I) after corrections: about 18.5/20.

### Sources for the methods (verified on Crossref unless stated)
* Renon, H.; Prausnitz, J. M. *AIChE J.* 1968, 14, 135–144, doi:10.1002/aic.690140124 (NRTL).
* Bonacina, C.; Comini, G.; Fasano, A.; Primicerio, M. *Int. J. Heat Mass Transfer* 1973, 16, 1825–1832,
  doi:10.1016/0017-9310(73)90202-0 (apparent-heat-capacity treatment of phase change smeared over a temperature interval).
* Piessens, R.; de Doncker-Kapenga, E.; Überhuber, C. W.; Kahaner, D. K. *QUADPACK*, Springer 1983,
  doi:10.1007/978-3-642-61786-7 (the adaptive Gauss–Kronrod routines behind `scipy.integrate.quad`).
* Peng, D.-Y.; Robinson, D. B. *Ind. Eng. Chem. Fundam.* 1976, 15, 59–64, doi:10.1021/i160057a011.
* Soave, G. *Chem. Eng. Sci.* 1972, 27, 1197–1203, doi:10.1016/0009-2509(72)80096-4.
* Rachford, H. H.; Rice, J. D. *J. Pet. Technol.* 1952, 4 (10), 19, doi:10.2118/952327-G.
* Kell, G. S. *J. Chem. Eng. Data* 1975, 20, 97–105, doi:10.1021/je60064a005 (water density).
* Pak, B. C.; Cho, Y. I. *Exp. Heat Transfer* 1998, 11, 151–170, doi:10.1080/08916159808946559 (Model I);
  Xuan, Y.; Roetzel, W. *Int. J. Heat Mass Transfer* 2000, 43, 3701–3707, doi:10.1016/S0017-9310(99)00369-5 (Model II).
* Smith, Van Ness, Abbott & Swihart, *Introduction to Chemical Engineering Thermodynamics*, 8th ed. (generic cubic,
  departure functions, Pitzer/Abbott B): textbook; edition flagged "verify" in `references.bib`.
* NIST Chemistry WebBook SRD 69 (reference data used in A1, C1 and the validation cases).

Established theory: C1–C4, C7–C12 (standard chemical-engineering thermodynamics). Modelling choices of this project: the
Gaussian melting width, the NRTL parameter sets, k_ij = 0, the ideal vapour in the VLE, and the uncertainty method.
