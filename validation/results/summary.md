| case | system | quantity | metric | N | mean | max | reference | status |
|---|---|---|---|---|---|---|---|---|
| 1 | water | cp (liquid, 1 atm) | AARD % | 20 | 0.0951 | 0.241 | NIST SRD 69 (Lemmon et al.) | ok |
| 1 | methanol | cp (liquid, 1 atm) | AARD % | 28 | 0.168 | 0.35 | NIST SRD 69 (Lemmon et al.) | ok |
| 1 | n-hexane | cp (liquid, 1 atm) | AARD % | 29 | 0.234 | 0.445 | NIST SRD 69 (Lemmon et al.) | ok |
| 1 | toluene | cp (liquid, 1 atm) | AARD % | 37 | 0.647 | 1.25 | NIST SRD 69 (Lemmon et al.) | ok |
| 1 | n-dodecane | cp (liquid, 1 atm) | AARD % | 14 | 0.353 | 0.794 | NIST SRD 69 (Lemmon et al.) | ok |
| 1 | acetone (liquid) | Cp (point values) | AARD % | 2 | 0.504 | 0.684 | NIST WebBook | ok |
| 1 | ethanol (liquid) | Cp (point values) | AARD % | 4 | 0.0626 | 0.15 | NIST WebBook | ok |
| 1 | glycerol (liquid) | Cp (point values) | AARD % | 6 | 1.47 | 2.37 | NIST WebBook + O'Hanley et al. 2012 | ok |
| 1 | n-octadecane (liquid) | Cp (point values) | AARD % | 2 | 0.317 | 0.634 | NIST WebBook | ok |
| 1 | n-octadecane (solid) | Cp (point values) | AARD % | 1 | 0 | 0 | NIST WebBook | ok |
| 1 | water | cp from NIST cv via cp = cv + Tvα²/κ_T | AARD % | 20 | 0.125 | 0.762 | NIST SRD 69 | ok |
| 2 | ethanol + water | H^E, NRTL (VLE parameters) | RMSD J/mol | 88 | 256 | 630 | doi:10.1016/j.fluid.2007.06.007; doi:10.1016/j.fluid.2014.05.031; doi:10.1016/j.jct.2005.06.018; doi:10.1016/j.tca.2017.05.023 | ok |
| 2 | ethanol + water | H^E, NRTL (refitted to H^E) | RMSD J/mol | 88 | 28.4 | 104 | doi:10.1016/j.fluid.2007.06.007; doi:10.1016/j.fluid.2014.05.031; doi:10.1016/j.jct.2005.06.018; doi:10.1016/j.tca.2017.05.023 | ok |
| 2 | methanol + water | H^E, NRTL (VLE parameters) | RMSD J/mol | 151 | 141 | 344 | doi:10.1016/j.jct.2005.03.012; doi:10.1016/j.tca.2006.02.028 | ok |
| 2 | methanol + water | H^E, NRTL (refitted to H^E) | RMSD J/mol | 151 | 23.2 | 98.1 | doi:10.1016/j.jct.2005.03.012; doi:10.1016/j.tca.2006.02.028 | ok |
| 2 | methanol + water | Cp, ideal mixture | AARD % | 12 | 7.47 | 12.8 | doi:10.1016/j.jct.2006.11.001 | ok |
| 2 | methanol + water | Cp, NRTL, VLE parameters | AARD % | 12 | 6.51 | 11.4 | doi:10.1016/j.jct.2006.11.001 | ok |
| 2 | methanol + water | Cp, NRTL refitted to H^E | AARD % | 12 | 2.71 | 7.54 | doi:10.1016/j.jct.2006.11.001 | ok |
| 2 | ethanol + toluene | Cp, ideal fallback (no NRTL parameters) | AARD % | 11 | 8.33 | 12.3 | doi:10.1016/j.jct.2011.06.019 | ok |
| 2 | glycerol + water | H^E, ideal fallback (H^E = 0) | RMSD J/mol | 9 | 466 | 611 | doi:10.1016/j.tca.2009.06.004 | ok |
| 2 | ethanol + water | mean Cp^E 298-323 K from measured H^E(T), NRTL, VLE parameters | RMSD J/(mol K) | 19 | 9.12 | 13.3 | doi:10.1016/j.jct.2005.06.018 | ok |
| 2 | ethanol + water | mean Cp^E 298-323 K from measured H^E(T), NRTL refitted to H^E | RMSD J/(mol K) | 19 | 2.31 | 4.1 | doi:10.1016/j.jct.2005.06.018 | ok |
| 2 | ethanol + water | Cp (mixture) | AARD % | 0 | – | – | see docs/DATA_TO_VERIFY.md | PENDING DATA |
| 3 | Al2O3 + water | cp, Model I (volume-weighted) | AARD % | 21 | 9.8 | 17.1 | O'Hanley et al. 2012, doi:10.1155/2012/181079 (digitised) | ok |
| 3 | Al2O3 + water | cp, Model II = Route A (mass-weighted) | AARD % | 21 | 2.39 | 5.28 | O'Hanley et al. 2012, doi:10.1155/2012/181079 (digitised) | ok |
| 3 | CuO + water | cp, Model I (volume-weighted) | AARD % | 24 | 34 | 63.3 | O'Hanley et al. 2012, doi:10.1155/2012/181079 (digitised) | ok |
| 3 | CuO + water | cp, Model II = Route A (mass-weighted) | AARD % | 24 | 2.41 | 4.94 | O'Hanley et al. 2012, doi:10.1155/2012/181079 (digitised) | ok |
| 3 | SiO2 + water | cp, Model I (volume-weighted) | AARD % | 24 | 10.7 | 20.3 | O'Hanley et al. 2012, doi:10.1155/2012/181079 (digitised) | ok |
| 3 | SiO2 + water | cp, Model II = Route A (mass-weighted) | AARD % | 24 | 1.42 | 3.62 | O'Hanley et al. 2012, doi:10.1155/2012/181079 (digitised) | ok |
| 3 | all nanofluids | cp, Model I | AARD % | 69 | 18.5 | 63.3 | O'Hanley et al. 2012, doi:10.1155/2012/181079 (digitised) | ok |
| 3 | all nanofluids | cp, Model II = Route A | AARD % | 69 | 2.06 | 5.28 | O'Hanley et al. 2012, doi:10.1155/2012/181079 (digitised) | ok |
| 4 | ice + water + ethanol | freezing point, ideal, ≤ 30 wt% (ice-slurry range) | |ΔT| K | 13 | 1.06 | 3.85 | Lange's Handbook (1967) via Wikipedia 'Ethanol (data page)' | ok |
| 4 | ice + water + ethanol | freezing point, ideal, all, 2.5-72 wt% | |ΔT| K | 18 | 2.46 | 7.45 | Lange's Handbook (1967) via Wikipedia 'Ethanol (data page)' | ok |
| 4 | ice + water + ethanol | freezing point, ideal + ΔCp (engine default), ≤ 30 wt% (ice-slurry range) | |ΔT| K | 13 | 0.842 | 3.08 | Lange's Handbook (1967) via Wikipedia 'Ethanol (data page)' | ok |
| 4 | ice + water + ethanol | freezing point, ideal + ΔCp (engine default), all, 2.5-72 wt% | |ΔT| K | 18 | 2.53 | 19.3 | Lange's Handbook (1967) via Wikipedia 'Ethanol (data page)' | ok |
| 4 | ice + water + ethanol | freezing point, NRTL γ (H^E fit), ≤ 30 wt% (ice-slurry range) | |ΔT| K | 13 | 3.31 | 9 | Lange's Handbook (1967) via Wikipedia 'Ethanol (data page)' | ok |
| 4 | ice + water + ethanol | freezing point, NRTL γ (H^E fit), all, 2.5-72 wt% | |ΔT| K | 17 | 5.04 | 14.3 | Lange's Handbook (1967) via Wikipedia 'Ethanol (data page)' | ok |
| 4 | ice + water + ethanol | freezing point, NRTL γ (VLE fit), ≤ 30 wt% (ice-slurry range) | |ΔT| K | 13 | 2.06 | 6.57 | Lange's Handbook (1967) via Wikipedia 'Ethanol (data page)' | ok |
| 4 | ice + water + ethanol | freezing point, NRTL γ (VLE fit), all, 2.5-72 wt% | |ΔT| K | 18 | 4.86 | 15.1 | Lange's Handbook (1967) via Wikipedia 'Ethanol (data page)' | ok |
| 4 | ice + water + ethanol | engine liquidus vs SLE equation | |ΔT| K | 13 | 1.5e-08 | 1.5e-08 | internal check | ok |
| 4 | ice_ethanol_10wt | ∫cp_app dT vs ΔH (energy conservation) | rel. dev % | 1601 | 0.109 | 0.109 | internal check | ok |
| 4 | pcm_octadecane_20wt | ∫cp_app dT vs ΔH (energy conservation) | rel. dev % | 1601 | 7.8e-09 | 7.8e-09 | internal check | ok |
| 5 | water (pure) | Route B (Q) vs Route A (ΔH_model) | rel. dev % | 1 | 6.5e-14 | 6.5e-14 | self-consistency | ok |
| 5 | water + ethanol 40 wt% | Route B (Q) vs Route A (ΔH_model) | rel. dev % | 1 | 6.2e-07 | 6.2e-07 | self-consistency | ok |
| 5 | methanol + water + ethanol + glycerol | Route B (Q) vs Route A (ΔH_model) | rel. dev % | 1 | 6.0e-07 | 6.0e-07 | self-consistency | ok |
| 5 | CuO–water nanofluid φ = 4 % | Route B (Q) vs Route A (ΔH_model) | rel. dev % | 1 | 6.7e-14 | 6.7e-14 | self-consistency | ok |
| 5 | water–toluene–ethanol emulsion + SiO2 | Route B (Q) vs Route A (ΔH_model) | rel. dev % | 1 | 8.3e-08 | 8.3e-08 | self-consistency | ok |
| 5 | ice slurry 10 wt% ethanol | Route B (Q) vs Route A (ΔH_model) | rel. dev % | 1 | 0 | 0 | self-consistency | ok |
| 5 | PCM slurry 20 wt% n-octadecane | Route B (Q) vs Route A (ΔH_model) | rel. dev % | 1 | 0 | 0 | self-consistency | ok |
