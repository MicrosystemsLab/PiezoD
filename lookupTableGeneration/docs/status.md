# FLOOXS Evaluation: Status and Conclusion

## Goal

Replace TSUPREM-4 (commercial TCAD, no longer available) with FLOOXS (open-source TCAD) for generating dopant profile lookup tables used by PiezoD. The lookup table covers 3 dopants (B, P, As) x 3 doses x 3 energies x 3 temperatures x 8 times x 2 oxidation states = 1296 simulations, each requiring accurate Transient Enhanced Diffusion (TED) modeling.

## Approaches tried

### 1. Fermi baseline (`ion_implant_fermi.tcl`)

Fermi-level dependent diffusion with no TED physics. Works reliably but under-predicts junction depth by ~47% at 1000C because implant damage-driven interstitial enhancement is not modeled. This is the expected lower bound.

### 2. QSS/Fermi hybrid (`ion_implant_ted_qss.tcl`)

Quasi-steady-state ScaleInter enhancement with algebraic defect equations. Boron 1000C Xj within 10%, but phosphorus and arsenic accuracy is poor (-23% to -54%). Arsenic ScaleInter gradients too steep for the Newton solver at 80 keV. Superseded by explicit pair diffusion.

### 3. Explicit pair diffusion, no clustering (`ion_implant_ted.tcl`)

5-stream pair PDEs (dopant-I and dopant-V pairs solved explicitly) without {311} clustering. Under-predicts Xj by 31-39% because there is no interstitial reservoir to sustain TED.

### 4. 5-stream + {311} with TSUPREM-4 Appendix A parameters (`ion_implant_5pd.tcl`)

All parameters from TSUPREM-4 manual Appendix A (Tables A-3, A-6, A-17, A-18/A-19, A-24). Total failure: 5 of 9 cases stuck in solver, remaining 4 show B/P Xj -52% to -78%, As 900C Xj +3025%.

Root cause: Table A-6 pair kinetics (Kr=10 s^-1 constant) give pair fractions of 0.3% at 900C and 11% at 1000C in an explicit pair PDE model. TSUPREM-4's internal effective diffusivity formulation does not have this pair fraction ceiling -- it uses D_eff = D_intrinsic * ScaleI directly. The Appendix A pair parameters are calibrated for that formulation, not for explicit pair PDEs.

### 5. 5-stream + {311} with FLOOXS parameters, v4-v5 (`ion_implant_5pd_v4.tcl`, `ion_implant_5pd_v5.tcl`)

FLOOXS 2026 binding-based pair kinetics (avoids the pair fraction problem), {311} clustering with Kr Ea=3.6 eV, per-dopant tuning (alpha_311, KinkSite). v5 adds V_eq initialization for boron Vac PDE.

Best explicit pair results: B 1000C Xj +4%, P 1000C Xj -1%, As 1000C Xj -5%. But 3 of 9 cases crash from solver stiffness (B 2e16 1100C, P 900C, P 2e16 1100C), and As has large errors at 900C (+97%) and 1100C (+45%).

### 6. 5-stream without Vac PDE (`ion_implant_5pd_v5_noVac.tcl`)

v5 with Vac PDE, BV pairs, and IV recombination removed for boron. Hypothesis: Vac/BV stiffness causes v5 B 1100C crash. Result: worse in every dimension. Without IV recombination as an interstitial sink, 1100C over-diffuses by +112%. All 2e16 cases fail. The BI pair reactions alone are stiff at high dose.

### 7. Effective diffusivity, v1 (`ion_implant_eff.tcl`)

Replace explicit pair PDEs with effective diffusivities (Table A-3 D values). Pair coupling terms derived from TSUPREM-4 Eq 3-89. 3-4 PDEs instead of 4-6. Boron with FLOOXS 2026 point defects, v5 {311}, alpha=0.2. P/As unchanged from TSUPREM-4 Appendix A.

Good B at 1000C (-2.2%) but 900C/1000C conflict persists (no single alpha works). 1100C crashes from {311}/Inter stiffness. As 1000C stuck from Vac PDE stiffness.

### 8. Effective diffusivity, v2 (`ion_implant_eff_v2.tcl`)

FLOOXS pair-derived D values (Ea=3.56 vs Table A-3's 3.46) with linear IV relaxation sink. Fixed 1100C over-diffusion (+112% -> +1.5%). Improved high-dose stability (2 of 3 v5 crash cases complete). But 900C accuracy poor (+52.5%) because FLOOXS pair D with alpha=0.1 over-diffuses when {311} clusters are frozen.

### 9. Effective diffusivity, v3 (`ion_implant_eff_v3.tcl`)

v1 Table A-3 D values + linear IV relaxation. Best low-dose boron accuracy: 1100C -1.3%, 1000C -3.6%, 900C +24.1%. But all 2e16 cases stuck/crash because Table A-3 D values create stiffer equations at high concentration than FLOOXS pair D.

### 10. Effective diffusivity, v4 (`ion_implant_eff_v4.tcl`)

v2 + soft ScaleInter cap (SmaxDop=200) in boron DiffDopI only + 10um mesh + no ramp down. First template to complete all 6 boron cases. But accuracy is mediocre: 900C -31.7%, 1000C -8.8%, 1100C -3.6%. High-dose profile shapes are wrong (Rs/Nz errors 40-83% despite reasonable Xj).

### 11. Built-in model procs (`ion_implant_builtin_fermi.tcl`)

Investigated FLOOXS built-in Tcl procs (DopantBulk, DopantFermi, Segregation). Found 3 C++ bugs and 3 Tcl bugs (see [flooxs_bugs.md](flooxs_bugs.md)). After patching, the procs work for Fermi-level diffusion but do not implement TED physics. Same as the manual Fermi baseline (-47% Xj).

## Fundamental problems

### a. FLOOXS D_I * C*_I diverges from TSUPREM-4 at 900C

| Temp | FLOOXS D_I * C*_I | TSUPREM-4 D_I * C*_I | Ratio |
|------|-------------------|----------------------|-------|
| 1000C | 4.2e6 | 3.1e6 | 1.4x |
| 900C | 1.9e4 | 3.3e5 | 0.06x (17x weaker) |

At 900C, FLOOXS intrinsic pair diffusivity is 17x too weak. Frozen {311} clusters (Ea=3.6 eV, tau=244min at 900C) compensate for boron but over-compensate for arsenic.

### b. B and As have contradictory {311} Ea requirements at 900C

Boron needs Ea=3.6 (frozen clusters to sustain ScaleI). Arsenic needs lower Ea (frozen clusters cause excess AsI diffusion at 900C). No single dopant-independent parameter set works.

### c. TSUPREM-4 Appendix A pair parameters are incompatible with explicit pair PDEs

Table A-6 pair kinetics (Kr=10 s^-1) give pair fractions too low for explicit pair PDEs. These parameters are designed for TSUPREM-4's internal effective diffusivity formulation where D_eff = D_intrinsic * ScaleI always.

### d. Solver stiffness at high dose and extreme temperatures

Multiple crash/stuck cases across all templates. The pair reaction terms, IV recombination, and {311} clustering equations become extremely stiff at 2e16 dose and at 900C/1100C extremes. No amount of per-dopant tuning eliminated all crash cases.

### e. 900C/1000C accuracy conflict

For boron with {311} clustering, the alpha parameter that optimizes 1000C (-2%) gives +24% at 900C, and vice versa. This is a fundamental consequence of the 17x D_I*C*_I mismatch combined with frozen clusters at 900C.

### f. Ramp phase instability

Ramp down (1100C -> 800C) causes solver divergence for extreme cases because EqInter drops 5700x, driving ScaleInter proportionally higher. Removing ramp down sacrifices thermal budget accuracy.

### g. Profile shape errors at high dose

Even when Xj is reasonable, high-dose cases show large Rs/Nz errors (40-780%) because the effective diffusivity cap and simplified physics produce wrong profile shapes.

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

## Decision

FLOOXS is not viable for replacing TSUPREM-4 for lookup table generation. The fundamental issues (D_I*C*_I mismatch, contradictory per-dopant requirements, solver stiffness, profile shape errors) cannot be resolved by parameter tuning alone. They stem from differences in the underlying PDE formulation and point defect parameters between FLOOXS and TSUPREM-4.

The existing TSUPREM-4-generated lookup table (`lookupTable.mat` -> `ionImplantLookupTable.h5`) remains the production data source. The FLOOXS templates, comparison data, and bug documentation are preserved for reference.
