# FLOOXS Evaluation: Status and Conclusion

## Goal

Replace TSUPREM-4 (commercial TCAD, no longer available) with FLOOXS (open-source TCAD) for generating dopant profile lookup tables used by PiezoD. The lookup table covers 3 dopants (B, P, As) x 3 doses x 3 energies x 3 temperatures x 8 times x 2 oxidation states = 1296 simulations, each requiring accurate Transient Enhanced Diffusion (TED) modeling.

## Why it failed

FLOOXS failed at two levels: its built-in physics models are inaccessible, and manually implementing equivalent PDEs hit fundamental parameter mismatches.

### Level 1: Built-in models are broken

FLOOXS has built-in Tcl procs (`DopantBulk`, `DopantFermi`, `Segregation`) that are supposed to set up diffusion equations from a parameter database. These don't work in FLOOXS 2026:

- 3 C++ bugs in the expression parser and term infrastructure (`Reduce.cc`, `Generic.cc`) cause crashes and corrupted expressions at solve time.
- 3 Tcl bugs in `Dopant.tcl` produce wrong PDB keys, spurious solved variables, and unguarded auto-creation of solution fields.
- The PDB parameter loading system is broken: the old Tcl-array loader is commented out, and the new C++ PDB has no file-loading mechanism. All parameters must be set explicitly.

After patching all 6 bugs ([flooxs_bugs.md](flooxs_bugs.md)), the built-in procs work for Fermi-level diffusion only. They do not implement TED physics (no pair diffusion, no {311} clustering, no point defect transport). So the built-in model layer is a dead end.

### Level 2: Manual PDEs hit fundamental parameter mismatches

With the built-in models unusable, we wrote all PDEs and parameters by hand (13 template iterations, 3-6 coupled PDEs per dopant). This required:

- Deriving all equations from the TSUPREM-4 manual (Appendix A tables, Eq 3-30 through 3-534)
- Implementing pair diffusion, {311} clustering, IV recombination, surface recombination, and charge neutrality in FLOOXS's Alagator equation language
- Working around FLOOXS solver limitations (no `pow()`, no `min`/`max`, brace stripping)

The manual approach got close at 1000C (Xj within 4% for boron and phosphorus) but cannot match TSUPREM-4 across the full parameter space due to:

**a. FLOOXS point defect parameters diverge from TSUPREM-4 at 900C.** The intrinsic interstitial diffusivity-concentration product (D_I * C*_I) is 17x weaker in FLOOXS than TSUPREM-4 at 900C, while similar at 1000C. This means any parameter set calibrated at 1000C is wrong at 900C.

**b. TSUPREM-4 parameters are incompatible with FLOOXS's PDE model.** TSUPREM-4's Appendix A pair kinetics (Table A-6) are designed for its internal effective diffusivity formulation (D_eff = D_intrinsic * ScaleI). In FLOOXS's explicit pair PDEs, the same parameters give pair fractions of 0.3% at 900C, making dopant effectively immobile.

**c. Contradictory per-dopant requirements.** Boron needs frozen {311} clusters (Ea=3.6 eV) to compensate for the weak pair diffusivity at 900C. Arsenic needs lower Ea because frozen clusters cause excess diffusion. No single dopant-independent parameter set works.

**d. Solver stiffness at high dose and extreme temperatures.** The coupled pair reaction / IV recombination / {311} clustering equations become extremely stiff at 2e16 dose and at 900C/1100C. Multiple crash/stuck cases persist across all 13 templates despite extensive per-dopant tuning.

**e. Profile shape errors.** Even when junction depth is reasonable, high-dose cases show 40-780% errors in sheet resistance and carrier density because the simplified physics produce wrong concentration profiles.

## Decision

FLOOXS is not viable for replacing TSUPREM-4 for lookup table generation. The existing TSUPREM-4-generated lookup table (`lookupTable.mat` -> `ionImplantLookupTable.h5`) remains the production data source.

## Best results achieved

### v5 (explicit pairs): Best 1000C accuracy

| Case | Xj error | Rs error | Status |
|------|----------|----------|--------|
| B 2e14 80keV 1000C 30min | +4.2% | -7.4% | Good |
| B 2e14 20keV 900C 30min | +6.5% | -36.3% | Xj OK, Rs poor |
| B 2e16 80keV 1100C 120min | -- | -- | CRASH |
| P 2e14 80keV 1000C 30min | -0.5% | -5.9% | Good |
| P 2e14 20keV 900C 30min | -- | -- | CRASH |
| As 2e14 80keV 1000C 30min | -4.8% | -12.5% | OK |
| As 2e14 20keV 900C 30min | +96.6% | -36.9% | Poor |

3 of 9 spot-check cases crash. Arsenic 900C severely over-predicts.

### eff_v4 (effective D): Best stability for boron

| Case | Xj error | Rs error | Status |
|------|----------|----------|--------|
| B 2e14 20keV 900C 30min | -31.7% | -28.6% | Poor |
| B 2e14 80keV 1000C 30min | -8.8% | -3.9% | Fair |
| B 2e14 80keV 1100C 120min | -3.6% | -3.5% | Good |
| B 2e16 20keV 900C 30min | +30.1% | -83.0% | Poor |
| B 2e16 80keV 1000C 30min | -7.5% | -43.8% | Poor |
| B 2e16 80keV 1100C 120min | +15.3% | -4.6% | Fair |

All 6 boron cases complete (first template to achieve this), but accuracy is insufficient for a lookup table (need <10% across all metrics, have up to 83% Rs error).

## Approaches tried (chronological)

1. **Fermi baseline** (`ion_implant_fermi.tcl`): No TED physics. Xj -47% at 1000C. Expected lower bound.
2. **QSS/Fermi hybrid** (`ion_implant_ted_qss.tcl`): B -10%, P -23%, As -54%. ScaleInter instability for As.
3. **Explicit pairs, no clustering** (`ion_implant_ted.tcl`): Xj -31% to -39%. No interstitial reservoir.
4. **5pd with TSUPREM-4 Appendix A params** (`ion_implant_5pd.tcl`): Total failure. Pair fraction too low.
5. **5pd with FLOOXS params, v4-v5** (`ion_implant_5pd_v4/v5.tcl`): Best 1000C. 3/9 crash, As 900C +97%.
6. **5pd without Vac PDE** (`ion_implant_5pd_v5_noVac.tcl`): Worse than v5 in every dimension.
7. **Effective D, v1** (`ion_implant_eff.tcl`): B 1000C -2.2%, but 900C/1000C conflict. 1100C crashes.
8. **Effective D, v2** (`ion_implant_eff_v2.tcl`): IV relaxation fixes 1100C. 900C +52.5%.
9. **Effective D, v3** (`ion_implant_eff_v3.tcl`): Best low-dose B accuracy. All 2e16 stuck/crash.
10. **Effective D, v4** (`ion_implant_eff_v4.tcl`): First 6/6 boron completion. Rs errors up to 83%.
11. **Built-in model procs** (`ion_implant_builtin_fermi.tcl`): 6 bugs patched. Fermi only, no TED.
