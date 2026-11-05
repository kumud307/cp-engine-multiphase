| test | reference | model | abs err | rel err % | criterion | result |
|---|---|---|---|---|---|---|
| A1 water liquid vs NIST | 4180.6 | 4182.95 | 2.35 | 0.0562 | ≤ 0.5 % | PASS |
| A2 water liquid: engine = DIPPR-100 evaluated independently | 4184.52 | 4184.52 | 0 | 0 | ≤ 1e-08 % | PASS |
| A3 Al2O3 solid: engine = Shomate evaluated independently | 96.7542 | 96.7542 | 0 | 0 | ≤ 1e-08 % | PASS |
| A3 SiO2 solid: engine = Shomate evaluated independently | 49.4701 | 49.4701 | 0 | 0 | ≤ 1e-08 % | PASS |
| A4 ethanol ideal gas: engine = DIPPR-107 evaluated independently | 95.9259 | 95.9259 | 0 | 0 | ≤ 1e-08 % | PASS |
| B1 ideal mixture, molar basis Σ x_i Cp_i | 86.4717 | 86.4717 | 0 | 0 | ≤ 1e-08 % | PASS |
| B2 ideal mixture, mass basis Σ w_i cp_i | 3271.58 | 3271.58 | 0 | 0 | ≤ 1e-08 % | PASS |
| B3 real solution: excess part present (information) | 86.4717 | 104.806 | 18.3 | 21.2 | informative | PASS |
| C1 nanofluid Al2O3 φ = 0.04: mass-weighted (Model II) | 3697.07 | 3696.38 | 0.693 | 0.0187 | ≤ 0.05 % | PASS |
| C2 suspension ΔH(290→330 K) = Σ n_i ∫Cp_i dT | 1.65157e+08 | 1.65157e+08 | 2.98e-08 | 1.8e-14 | ≤ 1e-06 % | PASS |
| D1 PCM slurry apparent Cp at 297.0 K | 3729.91 | 3729.91 | 3.18e-12 | 8.53e-14 | ≤ 0.0001 % | PASS |
| D1 PCM slurry apparent Cp at 300.3 K | 8965.82 | 8965.83 | 0.0105 | 0.000117 | ≤ 0.0001 % | FAIL |
| D1 PCM slurry apparent Cp at 301.3 K | 42453.3 | 42453.3 | 0.0258 | 6.08e-05 | ≤ 0.0001 % | PASS |
| D1 PCM slurry apparent Cp at 301.8 K | 27260.6 | 27260.6 | 5.6e-07 | 2.05e-09 | ≤ 0.0001 % | PASS |
| D1 PCM slurry apparent Cp at 303.0 K | 3907.86 | 3907.86 | 0.000843 | 2.16e-05 | ≤ 0.0001 % | PASS |
| D1 PCM slurry apparent Cp at 305.0 K | 3787.2 | 3787.2 | 1.01e-09 | 2.67e-11 | ≤ 0.0001 % | PASS |
| D2 latent area ∫ latent dT = n ΔH_fus (± ΔCp term) | 48491.9 | 48491.7 | 0.25 | 0.000515 | ≤ 0.1 % | PASS |
| D3 ΔH(290→312 K) through the melting step | 131142 | 131142 | 8.73e-11 | 6.66e-14 | ≤ 1e-06 % | PASS |
| E1 melt fraction = Gaussian CDF (σ = width/4) | 0 | 6.66134e-16 | 6.66e-16 | nan | ≤ 1e-12 abs | PASS |
| E2 Gaussian pdf normalised: ∫ f'(T) dT = 1 | 1 | 1 | 1.09e-14 | 1.09e-12 | ≤ 1e-10 abs | PASS |
| E3 '95 % within ±width/2' claim | 0.9545 | 0.9545 | 2.64e-07 | 2.76e-05 | ≤ 0.0005 abs | PASS |
| E4 Gauss–Legendre n = 2: ∫(v − T∂v/∂T)dp, water 0→100 MPa | 1610.65 | 1610.63 | 0.0211 | 0.00131 | ≤ 100 % | PASS |
| E4 Gauss–Legendre n = 4: ∫(v − T∂v/∂T)dp, water 0→100 MPa | 1610.65 | 1610.65 | 2.1e-06 | 1.3e-07 | ≤ 100 % | PASS |
| E4 Gauss–Legendre n = 8: ∫(v − T∂v/∂T)dp, water 0→100 MPa | 1610.65 | 1610.65 | 3.43e-11 | 2.13e-12 | ≤ 1e-06 % | PASS |
| E4 Gauss–Legendre n = 16: ∫(v − T∂v/∂T)dp, water 0→100 MPa | 1610.65 | 1610.65 | 4.34e-11 | 2.7e-12 | ≤ 1e-06 % | PASS |
| E5 Tait Cp correction (n = 8) vs closed form | -4.17608 | -4.17608 | 4.82e-10 | 1.15e-08 | ≤ 1e-06 % | PASS |
| E6 temperature stencil h = 0.5 K vs 0.05 K (truncation error) | -4.176 | -4.17608 | 7.74e-05 | 0.00185 | ≤ 0.5 % | PASS |
| E7 ∫Cp dT water liquid 273.16–373.15 K: QUADPACK vs closed form | 7550.92 | 7550.92 | 0 | 0 | ≤ 1e-09 % | PASS |
| E7 ∫Cp dT water gas 373.15–900.0 K: QUADPACK vs closed form | 19403 | 19403 | 3.64e-12 | 1.87e-14 | ≤ 1e-09 % | PASS |
| E7 ∫Cp dT Al2O3 solid 298.15–800.0 K: QUADPACK vs closed form | 53442.7 | 53442.7 | 7.28e-12 | 1.36e-14 | ≤ 1e-09 % | PASS |
| E7 ∫Cp dT ethanol gas 351.4–600.0 K: QUADPACK vs closed form | 22893.7 | 22893.7 | 7.28e-12 | 3.18e-14 | ≤ 1e-09 % | PASS |
| E8 interval-mean ideal part: QUADPACK vs 40-pt Gauss–Legendre | 88.6479 | 88.6479 | 4.26e-14 | 4.81e-14 | ≤ 1e-08 % | PASS |
| E9 interval mean total = ΔH/ΔT (no quadrature) | 104.705 | 104.705 | 0 | 0 | ≤ 1e-09 % | PASS |
| F1 component mass conservation over 255–315 K (121 T) | 0 | 2.22045e-16 | 2.22e-16 | nan | ≤ 1e-12 abs | PASS |
| F2 no negative phase amounts | 0 | 0 | 0 | nan | ≤ 0 abs | PASS |
| F3 allocation: toluene only organic, ethanol only aqueous (fixed rule) | 0 | 0 | 0 | nan | ≤ 0 abs | PASS |
| F4a PCM slurry 268–312 K (PCM solid at 268 K) is computable | 1 | 10701.8 | nan | nan | no refusal | PASS |
| F4 energy conservation ∫cp_app dT = ΔH, 290–312 K | 69736.1 | 69736.1 | 7.27e-05 | 1.04e-07 | ≤ 0.01 % | PASS |
| F6 phase mass fractions sum to 1 at 301.0 K | 1 | 1 | 0 | 0 | ≤ 1e-12 abs | PASS |
| F4 energy conservation ∫cp_app dT = ΔH, 255–285 K | 196213 | 196213 | 9.05e-05 | 4.61e-08 | ≤ 0.01 % | PASS |
| F6 phase mass fractions sum to 1 at 270.0 K | 1 | 1 | 0 | 0 | ≤ 1e-12 abs | PASS |
| F5 two PCMs melting 1 K apart: ΔH(290→312 K) = path-independent sharp-melting value | 144511 | 144511 | 8.64e-09 | 5.98e-12 | ≤ 1e-06 % | PASS |
| F5b same: ∫latent dT vs Σ n_i ΔH_fus,i (information) | 66369 | 65666.2 | 703 | 1.06 | informative | PASS |
| FD pure water 300 K: ΔT = 1.0 K | 75.356 | 75.3558 | 0.00016 | 0.000213 | ≤ 100 % | PASS |
| FD pure water 300 K: ΔT = 0.1 K | 75.3558 | 75.3558 | 1.6e-06 | 2.13e-06 | ≤ 100 % | PASS |
| FD pure water 300 K: ΔT = 0.01 K | 75.3558 | 75.3558 | 1.59e-08 | 2.11e-08 | ≤ 0.001 % | PASS |
| FD pure water 300 K: ΔT = 0.001 K | 75.3558 | 75.3558 | 1.57e-09 | 2.09e-09 | ≤ 0.001 % | PASS |
| FD water+ethanol NRTL 300 K: ΔT = 1.0 K | 104.675 | 104.674 | 0.000888 | 0.000849 | ≤ 100 % | PASS |
| FD water+ethanol NRTL 300 K: ΔT = 0.1 K | 104.674 | 104.674 | 9.79e-06 | 9.35e-06 | ≤ 100 % | PASS |
| FD water+ethanol NRTL 300 K: ΔT = 0.01 K | 104.674 | 104.674 | 1.73e-06 | 1.65e-06 | ≤ 0.001 % | PASS |
| FD water+ethanol NRTL 300 K: ΔT = 0.001 K | 104.674 | 104.674 | 1.05e-05 | 1e-05 | ≤ 0.001 % | PASS |
| FD CuO nanofluid 310 K: ΔT = 1.0 K | 4.13497e+06 | 4.13496e+06 | 6.82 | 0.000165 | ≤ 100 % | PASS |
| FD CuO nanofluid 310 K: ΔT = 0.1 K | 4.13496e+06 | 4.13496e+06 | 0.0682 | 1.65e-06 | ≤ 100 % | PASS |
| FD CuO nanofluid 310 K: ΔT = 0.01 K | 4.13496e+06 | 4.13496e+06 | 0.000681 | 1.65e-08 | ≤ 0.001 % | PASS |
| FD CuO nanofluid 310 K: ΔT = 0.001 K | 4.13496e+06 | 4.13496e+06 | 9.78e-05 | 2.37e-09 | ≤ 0.001 % | PASS |
| FD PCM slurry at Tm: ΔT = 1.0 K | 26904.6 | 42453.3 | 1.55e+04 | 57.8 | ≤ 100 % | PASS |
| FD PCM slurry at Tm: ΔT = 0.1 K | 42196.9 | 42453.3 | 256 | 0.608 | ≤ 100 % | PASS |
| FD PCM slurry at Tm: ΔT = 0.01 K | 42450.7 | 42453.3 | 2.55 | 0.00602 | ≤ 0.001 % | FAIL |
| FD PCM slurry at Tm: ΔT = 0.001 K | 42453.3 | 42453.3 | 8.85e-07 | 2.09e-09 | ≤ 0.001 % | PASS |
| FD ice slurry 10 wt% EtOH, 268 K: ΔT = 1.0 K | 51705.9 | 49711.8 | 1.99e+03 | 3.86 | ≤ 100 % | PASS |
| FD ice slurry 10 wt% EtOH, 268 K: ΔT = 0.1 K | 49731 | 49711.8 | 19.2 | 0.0386 | ≤ 100 % | PASS |
| FD ice slurry 10 wt% EtOH, 268 K: ΔT = 0.01 K | 49712 | 49711.8 | 0.19 | 0.000382 | ≤ 0.001 % | PASS |
| FD ice slurry 10 wt% EtOH, 268 K: ΔT = 0.001 K | 49711.8 | 49711.8 | 0.000181 | 3.65e-07 | ≤ 0.001 % | PASS |
| FD boiling EtOH/water x=0.3, 357 K: ΔT = 1.0 K | 3898.05 | 3750.22 | 148 | 3.79 | ≤ 100 % | PASS |
| FD boiling EtOH/water x=0.3, 357 K: ΔT = 0.1 K | 3751.63 | 3750.22 | 1.41 | 0.0376 | ≤ 100 % | PASS |
| FD boiling EtOH/water x=0.3, 357 K: ΔT = 0.01 K | 3750.24 | 3750.22 | 0.0139 | 0.000372 | ≤ 0.001 % | PASS |
| FD boiling EtOH/water x=0.3, 357 K: ΔT = 0.001 K | 3750.22 | 3750.22 | 2e-06 | 5.32e-08 | ≤ 0.001 % | PASS |
| FD steam 10 MPa 770 K: ΔT = 1.0 K | 45.0283 | 45.0283 | 1.64e-05 | 3.64e-05 | ≤ 100 % | PASS |
| FD steam 10 MPa 770 K: ΔT = 0.1 K | 45.0283 | 45.0283 | 2.24e-05 | 4.97e-05 | ≤ 100 % | PASS |
| FD steam 10 MPa 770 K: ΔT = 0.01 K | 45.0283 | 45.0283 | 2.28e-05 | 5.06e-05 | ≤ 0.001 % | PASS |
| FD steam 10 MPa 770 K: ΔT = 0.001 K | 45.0283 | 45.0283 | 2.28e-05 | 5.06e-05 | ≤ 0.001 % | PASS |
| G1 PR departure H^R water 773.15 K 10 MPa vs ∫[V − T(∂V/∂T)]dP | -2317.1 | -2317.1 | 0.000219 | 9.45e-06 | ≤ 0.05 % | PASS |
| G2 Cp^R = ∂H^R/∂T water | 6.58463 | 6.58465 | 2.21e-05 | 0.000335 | ≤ 0.1 % | PASS |
| G1 PR departure H^R n-hexane 550.0 K 2 MPa vs ∫[V − T(∂V/∂T)]dP | -3399.44 | -3399.44 | 0.00146 | 4.3e-05 | ≤ 0.05 % | PASS |
| G2 Cp^R = ∂H^R/∂T n-hexane | 16.4535 | 16.4537 | 0.000166 | 0.00101 | ≤ 0.1 % | PASS |
| G3 no double counting: Tait correction for gases is zero | 0 | 0 | 0 | nan | ≤ 0.0 abs | PASS |
| G4 H(T, 10 MPa) − H(T, 1 atm) of steam = n[H^R(P) − H^R(P0)] | -2294.76 | -2294.76 | 4.55e-12 | 1.98e-13 | ≤ 1e-08 % | PASS |
| G5 H_gas − H_liq at T_b, 1 atm = measured ΔH_vap (reference-state offset) | 40650 | 40650 | 0 | 0 | ≤ 1e-06 % | PASS |
| G6 departure term kept for a SATURATED vapour in the VLE region | 0.0552229 | 0.0552229 | 0 | 0 | ≤ 1.0 % | PASS |
