# Data sources and provenance

This file records **where every dataset in the project came from, how it was obtained, and how it
was checked**. Nothing in `data/` or `validation/data/` was invented. When a number could only be
taken from a secondary source, that is stated.

Last updated: 2 November 2026.

---

## 1. Pure-component data (`data/components.json`)

| What | Source | How it was checked |
|---|---|---|
| Liquid Cp of water, methanol, n-hexane, toluene, n-dodecane (DIPPR eq. 100) | Perry's Handbook 8th ed., Table 2-153 | Compared with the NIST SRD 69 tables we downloaded (validation case 1): AARD 0.10–0.65 % |
| Liquid Cp of ethanol and acetone | Perry's Table 2-153 | Compared with NIST WebBook measured values at 298 K: 0.06 % and 0.5 % |
| Liquid Cp of glycerol | Omel'chenko (1962) equation, as quoted on NIST WebBook (replaced the Perry's set, which was 4 % high) | NIST values and O'Hanley et al. (2012) Table 5: AARD 1.5 % |
| n-Octadecane liquid Cp | NIST WebBook (Parks 1949) | Taken directly from NIST |
| n-Octadecane solid Cp (linear in T, pre-melting rise) | Fit (`tools/fit_octadecane.py`) to Fonseca et al. 2014 (283.85–290.46 K) and Messerly et al. 1967 (298.15 K) | Max deviation 0.8 % |
| n-Octadecane liquid density (quadratic in T) | Fit to the specific volume of Espeau & Céolin 2006 (299–400 K) | Max deviation 0.25 %; κ_T consistent with Kariznovi et al. 2012 |
| n-Octadecane T_m and ΔH_fus | NIST WebBook (Messerly 1967) | Confirmed by 300.95–301.6 K and 61.34–61.65 kJ/mol in five later papers (ThermoML) |
| Ideal-gas Cp (DIPPR 107) of ethanol, methanol, acetone, n-hexane, toluene | Perry's Table 2-156 | Compared with NIST recommended ideal-gas Cp, 200–1000 K: max deviation 1.25 % |
| Vapour pressure (Antoine) of ethanol and toluene | NIST WebBook Antoine parameters | Reproduce the boiling points measured by Chen et al. 2012 within 0.7 % in pressure |
| Solid Cp of Al₂O₃, CuO, SiO₂ (quartz), NaCl (Shomate) | NIST WebBook (NIST-JANAF) | Every coefficient compared digit by digit with the NIST pages. The CuO range was corrected to 298–2000 K |
| Ice: Cp(T), ρ, α, κ_T | **IAPWS R10-06(2009)** (Feistel & Wagner 2006), primary standard | Equation implemented in `tools/iapws_ice.py`; reproduces the release's Table 6 test values to 9 digits. Replaced the earlier Fukusako correlation, which was 1.4 % low at 0 °C |
| Water density (Kell 1975 equation) | Kell, J. Chem. Eng. Data 1975 | 997.05 kg/m³ at 25 °C, matching NIST |
| Liquid density fits for methanol, hexane, toluene, dodecane | Quadratic fits (this work) to the NIST SRD 69 tables | Maximum deviation 0.06–0.18 % |
| Ethanol density (DIPPR 105) | Perry's Table 2-32 | 32 values from Lange's Handbook (3–40 °C): maximum deviation 0.16 % |
| κ_T(T) of water, methanol, hexane, toluene, dodecane | Fit of ln κ_T vs T to NIST SRD 69, κ_T = (Cp/Cv)/(ρw²) (`tools/fit_kappa_nist.py`) | Maximum deviation 0.6–1.7 %; used by the Tait pressure correction |
| Tm, ΔH_fus, Tb, ΔH_vap of the organics | NIST WebBook phase-change pages | Read off the pages on 2026-11-02 |
| Solid densities: Al₂O₃ 3987, CuO 6315, SiO₂ 2648, NaCl 2170 kg/m³ | CRC Handbook values as quoted in Wikipedia infoboxes | Corrected from the earlier 3970, 6310 and 2165 |
| NaCl solubility (0–100 °C) and ΔH_sol = 3.87 kJ/mol | CRC Handbook (mass %), as tabulated by Wikipedia | Checked the conversion: 26.28 wt% → 35.65 g per 100 g water |

Values still flagged `"verify": true` are listed in `docs/DATA_TO_VERIFY.md`. All of them have a
small effect on the main results; that file explains why for each one.

## 2. Binary data (`data/binary_nrtl.json`)

| What | Source | Notes |
|---|---|---|
| NRTL, VLE-fitted: ethanol–water, methanol–water | Aspen Plus databank (Dortmund Data Bank VLE regressions), as tabulated in the published NRTL appendix reproduced in BIT Mesra CAPE notes | Values and orientation confirmed against that table. Kept as `"fitted_to": "VLE"` |
| NRTL, VLE-fitted: acetone–water | Same table | **Corrected 2026-11-02:** i and j had been swapped. Checked with γ∞ (acetone in water) = 5.8 and γ∞ (water in acetone) = 3.6 at 298 K |
| **NRTL, water + ethanol (`"fitted_to": "HE+SLE"`)** | **This work** (`tools/fit_nrtl_he.py`), fitted jointly to measured Hᴱ (NIST ThermoML) and to the freezing points of aqueous ethanol (Lange) | Hᴱ RMSD 256 → 34 J/mol (fit), 46 J/mol held-out; freezing points 1.7 K (fit), 2.5 K held-out. Phase-stable 220–390 K. **Default** for Cp and for the ice-slurry equilibrium |
| **NRTL, water + methanol (`"fitted_to": "HE"`)** | **This work**, fitted to measured Hᴱ (NIST ThermoML) | Hᴱ RMSD 141 → 23 J/mol (fit), 48 J/mol held-out. Default for Cp |
| **Apparent molar heat capacity of NaCl(aq)**, φCp = c₀ + d₀√m + e·m at 298 K | **This work** (`tools/fit_apparent_cp.py`), fitted to Tromans et al. 2005 (m = 1–3 mol/kg) | Independent test (Gao et al. 2017, m = 2.05): 0.17 % in Cp. A third dataset (Tyunina et al. 2019) was **excluded**: its NaCl values differ from the dedicated study by 1.1 % |
| **NRTL, ethanol + toluene (`"fitted_to": "VLE+CpE"`)** | **This work**, fitted jointly to isobaric VLE (Chen et al. 2012) and mixture Cp (Vega-Maza et al. 2011) | ln γ RMS 0.042; Cpᴱ 1.5 J/(mol K) (fit), 3.0 held-out; mixture Cp 1.3 % (ideal mixing: 8.3 %) |
| Miscibility flags | CRC solvent miscibility tables | Partially miscible pairs are flagged |
| **Methanol + n-hexane liquid–liquid coexistence curve** | Trejo et al. 2006 (14 points, UCST ≈ 314.9 K) | The engine decides one or two liquid phases from T and composition |
| Methanol + n-dodecane immiscibility | Casas et al. 2004 (mutual solubilities 278–308 K) | Confirms the flag |

## 3. Experimental validation data (`validation/data/`)

| File | Content | Source and how it was obtained |
|---|---|---|
| `nist_<fluid>_1atm.csv` (and `nist_raw_*.tsv`) | Cp, Cv, ρ, speed of sound of water, methanol, n-hexane, toluene, n-dodecane at 0.101325 MPa | **Downloaded** from NIST Chemistry WebBook SRD 69 (Lemmon, McLinden & Friend), fluid.cgi, 2026-11-02 |
| `nist_<fluid>_298K_isotherm.csv` (and `nist_raw_iso_*.tsv`) | Cp and ρ of water, methanol, toluene, n-hexane, n-dodecane at 298.15 K, 0.1–100 MPa | **Downloaded** from NIST SRD 69 (isotherms), 2026-11-02; used to validate the pressure correction |
| `pure_component_points.csv` | 15 measured Cp values (ethanol, acetone, glycerol, n-octadecane) | NIST WebBook condensed-phase tables (original references listed in each row), plus O'Hanley et al. 2012 Table 5 (glycerin) |
| `thermoml_binary_data.csv` | 1090 points of Hᴱ, Cp and apparent molar Cp for binaries of database components (compositions given as mole fraction, mass fraction or molality, all converted to mole fraction) | **Extracted** by `tools/extract_thermoml.py` from the NIST TRC ThermoML Archive, version 2020-09-30 (doi:10.18434/mds2-2422, 189 MB, 11,923 files). Each row carries the DOI of the original paper |
| `nanofluid_ohanley2012.csv` | 69 DSC cp values of Al₂O₃, CuO and SiO₂ in water, 35/45/55 °C | O'Hanley, Buongiorno, McKrell & Hu, *Adv. Mech. Eng.* 2012, 181079 (open access, CC BY). The paper gives the data only in figures. Values were **digitised from the vector graphics** by `tools/digitize_ohanley2012.py`: marker centres read from the PDF drawing commands, axes calibrated from the tick labels. Concentrations come from the paper's Table 2 |
| `ethanol_water_freezing_points.csv` | 20 freezing points of water + ethanol, 0–100 wt% | Lange's Handbook of Chemistry (1967), table "Properties of aqueous ethanol solutions", as reproduced in Wikipedia "Ethanol (data page)" (raw wikitext parsed) |
| `thermoml_vle_ethanol_toluene.csv` | Isobaric T–x–y of ethanol + toluene at 101.3–201.3 kPa (68 points) | **Extracted** by `tools/extract_thermoml_phase_eq.py` (Chen, Zhong & Xu 2012) |
| `thermoml_gamma_inf.csv` | Activity coefficients at infinite dilution in water (ethanol, acetone, toluene) | **Extracted** (Straver & de Loos 2005; Atik et al. 2004) |
| `thermoml_lle.csv` | Liquid–liquid equilibrium of methanol + n-hexane and methanol + n-dodecane | **Extracted** (Trejo et al. 2006; Casas et al. 2004) |
| `nist_ideal_gas_cp.csv` | Ideal-gas Cp of 5 components, 200–1000 K | NIST WebBook gas-phase tables (TRC 1997; Chao 1986; Scott 1974; Draeger 1985) |
| `water_ethanol_cp_exp.csv` | **Empty template.** Direct Cp of water + ethanol mixtures | No open dataset was found (it is in subscription journals). Partly replaced by the Cpᴱ derived from measured Hᴱ(T) (validation case 2g) |

### Papers whose data are used (from the ThermoML archive)

| DOI | Paper | Used for |
|---|---|---|
| 10.1016/j.jct.2005.06.018 | Nagamachi & Francesconi, *J. Chem. Thermodyn.* 2006, 38, 461 | Hᴱ of water + ethanol at 298.15 and 323.15 K (NRTL fit; Cpᴱ from Hᴱ(T)) |
| 10.1016/j.fluid.2014.05.031 | Brandt, Horstmann, Steinigeweg & Gmehling, *Fluid Phase Equilib.* 2014, 376, 48 | Hᴱ of water + ethanol at 333–383 K |
| 10.1016/j.tca.2017.05.023 | Modaresifar & Kowalski, *Thermochim. Acta* 2017, 655, 137 | One Hᴱ point, water + ethanol |
| 10.1016/j.fluid.2007.06.007 | Li, Yan & Dong, *Fluid Phase Equilib.* 2007, 261, 444 | Water + ethanol Hᴱ near infinite dilution (mostly excluded by the 0.01 < x < 0.99 filter) |
| 10.1016/j.jct.2005.03.012 | Dettmann, Ernst & Wirbser, *J. Chem. Thermodyn.* 2006, 38, 56 | Hᴱ of water + methanol, 323–378 K |
| 10.1016/j.tca.2006.02.028 | Behbehani & Ghamamy, *Thermochim. Acta* 2006, 444, 71 | Hᴱ of water + methanol at 298 K |
| 10.1016/j.jct.2006.11.001 | Tanaka, Fujita & Uematsu, *J. Chem. Thermodyn.* 2007, 39, 961 | **Cp of water + methanol** at 280 and 320 K, 0.1 MPa (independent test) |
| 10.1016/j.jct.2011.06.019 | Vega-Maza et al., *J. Chem. Thermodyn.* 2011, 43, 1893 | Cp of ethanol + toluene at 313.15 K |
| 10.1016/j.tca.2009.06.004 | Zehioua, Coquelet, Soo, Richon & Meniai, *Thermochim. Acta* 2009, 495, 72 | Hᴱ of water + glycerol at 298.15 K |
| 10.1021/je0502400 | Tromans, Königsberger, May & Hefter, *J. Chem. Eng. Data* 2005, 50, 2019 | Cp of NaCl(aq), 1–3 mol/kg, 298.15 K (φCp fit) |
| 10.1016/j.jct.2016.09.031 | Gao, Li & Li, *J. Chem. Thermodyn.* 2017, 104, 201 | Cp of NaCl(aq), 2.05 mol/kg (independent test) |
| 10.1016/j.jct.2018.10.027 | Tyunina, Badelin & Mezhevoi, *J. Chem. Thermodyn.* 2019, 131, 40 | NaCl(aq) Cp reported incidentally; **excluded** as inconsistent (see `tools/fit_apparent_cp.py`) |
| 10.1016/j.jct.2006.08.009 | Blodgett, Ziemer, Brown, Niederhauser & Woolley, *J. Chem. Thermodyn.* 2007, 39, 627 | Apparent molar Cp of glycerol in water, 278–363 K (162 points; tests the ideal fallback) |
| 10.1021/je200921u | Chen, Zhong & Xu, *J. Chem. Eng. Data* 2012, 57, 155 | Isobaric VLE of ethanol + toluene (NRTL fit) |
| 10.1021/je0495942 | Straver & de Loos, *J. Chem. Eng. Data* 2005, 50, 1171 | γ∞ of ethanol and acetone in water, 298 K |
| 10.1021/je049875+ | Atik, Gruber, Krummen & Gmehling, *J. Chem. Eng. Data* 2004, 49, 1429 | γ∞ of ethanol, acetone, toluene in water, 303–343 K |
| 10.1021/je0505321 | Trejo, Yáñez & Eustaquio-Rincón, *J. Chem. Eng. Data* 2006, 51, 1070 | Methanol + n-hexane coexistence curve |
| 10.1016/j.jct.2003.11.009 | Casás, Orge, Díaz & Tojo, *J. Chem. Thermodyn.* 2004, 36, 237 | Methanol + n-dodecane LLE |
| 10.1016/j.jct.2013.11.038 | Fonseca, Gushterov & Dohrn, *J. Chem. Thermodyn.* 2014, 73, 148 | Solid Cp, T_m, ΔH_fus of n-octadecane (measured alongside the vapour pressures) |
| 10.1016/j.tca.2006.04.001 | Espeau & Céolin, *Thermochim. Acta* 2006, 445, 32 | Liquid specific volume of n-octadecane |
| 10.1016/j.fluid.2011.10.018 | Kariznovi, Nourozieh & Abedi, *Fluid Phase Equilib.* 2012, 314, 102 | Density of n-octadecane at 323 K, 2–9.5 MPa (κ_T check) |

## 4. Re-creating everything from scratch

```bash
python validation/_clean_nist.py                      # NIST raw tables -> CSV, density fits, kappa_T
python tools/extract_thermoml.py <ThermoML folder>    # needs ThermoML.v2020-09-30.tgz extracted
python tools/digitize_ohanley2012.py <pdf>            # needs the open-access PDF from hdl.handle.net/1721.1/72068
python tools/iapws_ice.py                             # verifies the IAPWS ice equation and prints the fit
python tools/fit_kappa_nist.py                        # kappa_T(T) fits for the pressure correction
python tools/extract_thermoml_phase_eq.py <folder>    # VLE, gamma_inf and LLE datasets
python tools/fit_octadecane.py                        # n-octadecane solid Cp and liquid density
python tools/fit_nrtl_he.py                           # refits NRTL (H^E, freezing points, VLE + Cp^E) and stores them
python tools/fit_apparent_cp.py                       # NaCl(aq) apparent molar heat capacity
python validation/case2b_crossvalidation.py --recompute   # held-out (cross-validation) refits
python validation/run_all.py                          # all comparisons, figures and tables
python tools/make_data_to_verify.py                   # regenerates docs/DATA_TO_VERIFY.md
```

## 5. Corrections made during verification (worth mentioning in the defence)

1. **Acetone–water NRTL parameters were stored with i and j swapped.** This was caught by comparing
   with the published table and checking infinite-dilution activity coefficients.
2. **Glycerol liquid Cp** from the remembered Perry's set was 4 % high, so it was replaced.
3. **Ice Cp** (Fukusako correlation) was 1.4 % low at 0 °C, so it was replaced by the IAPWS standard.
4. **CuO Shomate range**, solid densities and the NaCl solubility table were corrected to the reference values.
5. **VLE-fitted NRTL parameters underpredict Cpᴱ of water + alcohols six-fold or more.**
   Refitted Hᴱ parameters were added and made the default for Cp.
6. **After a mock examiner review (2 Nov 2026):**
   * pressure was ignored → Tait pressure correction, validated on NIST isotherms to 100 MPa;
   * a "validation" reused fitted data → relabelled as a consistency check; held-out tests added (case 2b);
   * the ice-slurry equilibrium used γ = 1 while the enthalpy used NRTL → joint Hᴱ + freezing-point fit,
     same NRTL for both;
   * brine Cp 2 % high (salt counted with its solid Cp) → apparent molar heat capacity of NaCl(aq);
   * one NaCl dataset found inconsistent and excluded, with the reason documented.
7. **Closing the remaining gaps (5 Nov 2026):**
   * ethanol + toluene: from ideal mixing (8.3 % error) to NRTL fitted to VLE + Cp^E (1.3 %);
   * methanol + n-hexane: "treated as immiscible" replaced by the measured coexistence curve (UCST 314.9 K);
   * n-octadecane solid Cp, liquid density and κ_T, ideal-gas Cp of 5 components, acetone + water NRTL:
     checked against measured data and the flags cleared;
   * the acetone + water literature set predicts a spurious phase split above 360 K: the engine now
     checks every NRTL pair it uses and warns;
   * adding γ∞ data to the water + ethanol fit was tried and rejected (γ∞ better, Cp^E and freezing points
     worse); the trade-off is reported in validation case 2.
