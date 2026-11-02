# Detailed flowchart with governing equations (for the manuscript)

Sign convention (IUPAC): Q > 0 heat added to the system; W > 0 work done on the system.

## 1. Input → normalisation (`inputs.py`)
- Units → SI (`units.py`): T [K], P [Pa], energy [J], mass [kg], amount [mol].
- Composition → moles nᵢ:
  - mole fractions: nᵢ = xᵢ (basis 1 mol)
  - mass fractions: nᵢ = wᵢ / Mᵢ (basis 1 kg)
  - volume fractions: nᵢ = φᵢ ρᵢ(T) / Mᵢ (basis 1 m³, unmixed volumes)
  - solid volume fraction φ_s with fluid composition: m_s = φ_s ρ_s, m_f = (1 − Σφ_s)/v_f, v_f = Σ w_f,i / ρᵢ
- Missing fraction: x_k = 1 − Σ_{i≠k} xᵢ.
- Specific energies: Q = q·m or Q = q·n (basis 1 kg / 1 mol if no amount; Cp is intensive).

## 2. Validation (`validator.py`)
- T > 0, P > 0, 0 ≤ xᵢ ≤ 1, |Σxᵢ − 1| ≤ 10⁻³, nᵢ ≥ 0, T₂ − T₁ = ΔT, γ ≥ 1.
- Declared solid with T > T_m + 0.5 K → error; liquid below T_m or above T_b(P) → warning.
- T_b(P) from Clausius–Clapeyron: 1/T_b(P) = 1/T_b − (R/ΔH_vap) ln(P/P°).
- Correlation range: warn outside [T_min, T_max]; refuse beyond 30 K.

## 3. Forward-chaining solver (`solver.py`, `rules.py`)
Knowledge base K ← user inputs. Repeat until no rule fires: if inputs(rule) ⊆ K, add outputs(rule);
if an output is already in K, compare (|Δ|/max > 2 % → inconsistency warning).

### Helper rules
- ΔT = T₂ − T₁; T = (T₁ + T₂)/2; T₁,₂ = T ∓ ΔT/2
- M = Σ xᵢMᵢ; m = nM; n = m/M; c_v = C_v,m / M; c_v = C_v / m

### Route A – property route (`system.py`, `models/`)
- Pure components: DIPPR-100 Cp = Σ Cₖ Tᵏ; DIPPR-107 Cp = C₁ + C₂[(C₃/T)/sinh(C₃/T)]² + C₄[(C₅/T)/cosh(C₅/T)]²; Shomate Cp = A + Bt + Ct² + Dt³ + E/t², t = T/1000.
- Enthalpies (Kirchhoff): H_L(T) = ∫_{T_m}^{T} Cp_L dT; H_S(T) = −ΔH_fus + ∫_{T_m}^{T} Cp_S dT; H_G(T) = H_L(T_b) + ΔH_vap + ∫_{T_b}^{T} Cp_G dT.
- Ideal mixture: Cp = Σ xᵢ Cpᵢ (molar) ⇔ c_p = Σ wᵢ c_p,i (mass).
- Real liquid solution (NRTL): G^E/RT = Σᵢ xᵢ Σⱼ τⱼᵢGⱼᵢxⱼ / Σₖ Gₖᵢxₖ; τᵢⱼ = aᵢⱼ + bᵢⱼ/T + eᵢⱼ ln T + fᵢⱼT; Gᵢⱼ = exp(−αᵢⱼτᵢⱼ).
  - H^E = −RT² ∂(G^E/RT)/∂T; Cp^E = (∂H^E/∂T)_{P,x} = −R(2Tg′ + T²g″), central differences with h = 10⁻⁴T.
  - Cp = Σ xᵢ Cpᵢ + Cp^E.
  - Parameter choice: sets refitted to measured H^E at several temperatures (default; accepted only if ΔG_mix/RT is convex in x from 240 to 390 K) rather than VLE-regressed sets.
- Heterogeneous mixture: c_p = Σ_phases w_p c_p,p (mass weighting).
  - Model I (inconsistent): c_p = φ c_p,s + (1 − φ) c_p,f.
  - Model II (thermal equilibrium): c_p = [φρ_s c_p,s + (1 − φ)ρ_f c_p,f] / [φρ_s + (1 − φ)ρ_f] ≡ mass weighting.
- Phase change:
  - Smooth melting of a pure component: f(T) = ½[1 + erf((T − T_m)/(σ√2))], σ = width/4.
  - Solid–liquid equilibrium: ln(x_k γ_k) = −(ΔH_fus/R)(1/T − 1/T_m) + (ΔCp/R)(T_m/T − 1 − ln(T_m/T)), ΔCp = Cp,L − Cp,S at T_m (consistent with the Kirchhoff enthalpies; γ_k = 1 by default); tabulated solubility for electrolytes.
  - Total enthalpy: H(T) = Σ_p Σᵢ nᵢᵖ Hᵢᵖ(T) + Σ_liquid nᵖ H^E,p(T, xᵖ).
  - Apparent Cp: dH/dT = Σ nᵢᵖCpᵢᵖ (ideal) + Σ nᵖCp^E,p (excess) + Σ (∂H/∂nᵢᵖ)(dnᵢᵖ/dT) (latent);
    latent ≈ [H(T; n(T+h)) − H(T; n(T−h))]/(2h).
  - Interval mean: c̄_p = [H(T₂) − H(T₁)] / [m(T₂ − T₁)].
  - Saturated salt slurry: c_p,app = Σ wᵢ c_p,i + Δh_sol · dw_dissolved/dT.

### Route B – calorimetric
- Adiabatic: Q = 0.
- Constant P: ΔH = Q + W → c_p = ΔH/(mΔT).
- Constant V: ΔU = Q + W → c_v = ΔU/(mΔT) → c_p = c_v + Tvα²/κ_T.
- General: ΔU = Q + W + W_pv.

### Route C – derived properties
- c_p = ΔH/(mΔT) (initial and final P equal).
- Constant P, constant c_p, no phase change: ΔS = m c_p ln(T₂/T₁) → c_p = ΔS/[m ln(T₂/T₁)].
- c_p = c_v + Tvα²/κ_T with mixture α = Σφᵢαᵢ, κ_T = Σφᵢκ_T,i, v = Σ wᵢ/ρᵢ (V^E = 0).
- c_p = γ c_v.
- ΔH = ΔU + PΔV (constant P).

## 4. Output
- Candidates c_p[route]; recommended = first available of {ΔH, Cv, γ, ΔS (measured)} then {Route A mean, Route A point}.
- Discrepancy: (c_p,route − c_p,A)/c_p,A × 100 %, comparing like with like (point/mean).
- If no candidate: backward search over rules for minimal sets of missing base inputs.
