# Work breakdown – five owners, five modules

Each member owns one part of the code, must be able to explain every equation in it,
and runs the listed tests before every commit. The interfaces between parts are
small, so the parts can be worked on in parallel.

| # | Owner | Files | Thermodynamics to defend | Tests |
|---|---|---|---|---|
| 1 | Member 1 – **Data & units** | `units.py`, `database.py` (incl. Tait pressure correction), `data/*.json`, `models/pure.py`, `tools/iapws_ice.py`, `tools/fit_kappa_nist.py`, `tools/fit_vapour_pressure.py`, `docs/DATA_SOURCES.md`, `docs/DATA_TO_VERIFY.md` | DIPPR-100/107 and Shomate forms, Kirchhoff's law for enthalpy paths, Antoine equation and T_b(P), critical constants and the definition of ω, why ρ and α come from one density model | `test_units.py`, `test_database.py` |
| 2 | Member 2 – **Input layer** | `inputs.py`, `validator.py`, `interactive.py`, `errors.py`, `custom.py` | Fraction conversions (x ↔ w ↔ φ), ideal volume of mixing, phase rule checks (T vs T_m, T_b), first-law sign convention, choice of calculation basis | `test_edge_cases.py` |
| 3 | Member 3 – **Solution thermodynamics** | `models/ideal.py`, `models/nrtl.py`, `models/partial_molar.py`, `models/cv_cp.py`, `tools/fit_nrtl_he.py`, `tools/fit_apparent_cp.py` | Ideal mixing, NRTL G^E, Gibbs–Helmholtz (H^E, Cp^E), numerical derivatives and step size, VLE vs Hᴱ-fitted parameters and the phase-stability check, partial molar properties and Gibbs–Duhem, mixing and excess functions, Cp − Cv = Tvα²/κ_T | `test_models.py` (NRTL, Cv-Cp parts), `test_eos_vle.py` (partial molar part) |
| 4 | Member 4 – **Multiphase & phase change** | `system.py`, `models/heterogeneous.py`, `models/phase_change.py`, `models/vle.py`, `models/eos.py` | Classification (solution / suspension / emulsion / slurry), Model I vs Model II, smoothed melting, solid–liquid equilibrium with the ΔCp term / freezing-point depression, vapour–liquid equilibrium (γ–φ, bubble/dew points, Rachford–Rice flash), equations of state (virial, vdW/RK/SRK/PR, residual Cp, fugacity coefficients, φ–φ bubble pressure), apparent Cp = dH/dT split into ideal + excess + latent | `test_models.py` (system parts), phase-change tests in `test_edge_cases.py`, `test_eos_vle.py` |
| 5 | Member 5 – **Solver, validation & reporting** | `rules.py`, `solver.py`, `uncertainty.py`, `plotting.py`, `validation/` (incl. `case2b_crossvalidation.py`, `case6_vle_real_gas.py`), `tools/extract_thermoml*.py`, `tools/fetch_nist_fluid.py`, `tools/make_uncertainty_table.py`, `tools/digitize_ohanley2012.py`, notebook, manuscript | Forward chaining, the three routes, consistency checks, missing-input search and unused-input warnings, AARD statistics, literature comparison, how the uncertainty estimate is built | `test_solver.py`, `validation/run_all.py` |

## Suggested schedule to Fri 20 Nov 2026

Status on 5 Nov: code complete (196 tests pass; independent audit in MODEL_VERIFICATION_REPORT.md); all validation data except a direct Cp dataset for
water + ethanol and DSC data for a PCM emulsion are in place; ethanol + toluene, methanol + n-hexane and
n-octadecane gaps closed; manuscript drafted with real results.

- **Mon 2 – Fri 6 Nov:** each owner reads their module, its tests and the matching section of
  `docs/TEAM_HANDBOOK.md`, then runs `pytest` and `python validation/run_all.py`.
- **Mon 9 – Wed 11 Nov:** Member 1 works through the category-A items in `docs/DATA_TO_VERIFY.md`.
  Member 3 or 5 tries the IITGN library for a direct water + ethanol Cp dataset (optional; the case-2
  template reads it automatically).
- **Thu 12 – Fri 13 Nov:** fill in author names, e-mail and acknowledgement in `manuscript/main.tex`;
  compile on Overleaf; proofread against `validation/results/summary.md`.
- **Mon 16 – Wed 18 Nov:** mock defence. Each member explains *another* member's module, using
  section 9 of `docs/TEAM_HANDBOOK.md`.
- **Thu 19 Nov:** run the notebook in Colab from a fresh runtime; freeze the code (tag a release).
- **Fri 20 Nov:** demo (follow section 3 of the handbook), manuscript, defence.

## Oral-defence questions to rehearse (one per module; 32 with answers in `docs/TEAM_HANDBOOK.md` §9)

1. Why must the mass basis, not a volume fraction, weight specific heats? (Member 4)
2. Where does Cp^E come from, and why do VLE-fitted NRTL parameters predict it poorly? (Member 3)
3. Show that ΔH = Q + W at constant P when W is shaft work. (Member 5)
4. How does the engine know that water at −10 °C in 10 wt% ethanol is partly ice? (Member 4)
5. What happens if a mass fraction is missing or the fractions add to 1.2? (Member 2)
6. Why is Kirchhoff's law needed to get the latent heat at a temperature other than T_m? (Member 1)
