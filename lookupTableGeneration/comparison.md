# FLOOXS vs TSUPREM-4 Comparison

## How to generate comparison values

### TSUPREM-4 reference

Extract from the lookup table (`python/src/piezod/data/ionImplantLookupTable.h5`):

```python
import h5py
f = h5py.File('python/src/piezod/data/ionImplantLookupTable.h5', 'r')
# Indices: [dopant, dose, energy, temp, time, oxidation]
# dopant: B=0, P=1, As=2
# dose: 2e14=0, 2e15=1, 2e16=2
# energy: 20keV=0, 50keV=1, 80keV=2
# temp: 900C=0, 1000C=1, 1100C=2
# time: 15=0, 30=1, 45=2, 60=3, 75=4, 90=5, 105=6, 120=7
# oxidation: no oxide=0, oxide=1
# Units: Xj in meters, Rs in ohm/sq, Nz in cm^-2
Xj_um = f['Xj'][d, do, e, t, ti, ox] * 1e6
```

### FLOOXS simulation

Generate a test script from the template and run in Docker:

```bash
cd lookupTableGeneration
bash scripts/gen_test.sh <dopant> <dose> <energy> <temp> <time> templates/ion_implant_5pd_v5.tcl > simulations/test.tcl
MSYS_NO_PATHCONV=1 docker compose run --rm flooxs /work/test.tcl
```

Rs, Nz, Beta1, Beta2 are post-processed from the concentration profile output:

```bash
python scripts/compare_results.py simulations/output.txt --dopant <dopant> --dose <dose> --energy <energy> --temp <temp> --time <time>
```

## Template parameters

`ion_implant_5pd_v5.tcl` - 5-stream pair diffusion + {311} clustering:
- All dopants: explicit interstitial transport, {311} clustering, Hobler-Moroz f_pl
- B: BI + BV pairs, IV recombination, V_eq initialization
- P: PI pairs only (I-only diffuser), alpha=0.6, constant KinkSite=1e11
- As: AsI + AsV pairs, IV recombination, V_eq initialization
- Per-dopant alpha_311: 0.1 (B,As) / 0.6 (P)
- Per-dopant KinkSite: Arrhenius (B,As) / constant 1e11 (P)

## eff_v4 template: v2 + ScaleInter cap + 10um mesh + no ramp down (2026-02-24)

Template: `ion_implant_eff_v4.tcl`

Physics: v2 effective diffusivity (FLOOXS pair D + linear IV relaxation) with a soft
ScaleInter cap in boron DiffDopI only, extended 10um mesh, and no ramp down.

Changes from v2:
1. Mesh: 5um -> 10um (captures Xj > 5um for 2e16 1100C)
2. DiffDopI cap: `ScaleCapDop = ScaleInter * 200 / (ScaleInter + 200)` applied only in
   boron DiffDopI. Prevents flux divergence at extreme supersaturation (2e16 + frozen
   {311} at 900C where ScaleInter ~ 1e5). I self-diffusion, R_clust, and IV relaxation
   use real ScaleInter.
3. Ramp down removed: causes solver divergence on 10um mesh (EqInter drops 5700x
   during 1100C -> 800C cooldown). Ramp up kept for thermal budget.

### SmaxDop sweep

| SmaxDop | Case 4 (2e16 900C) | Case 5 (2e16 1000C) | Low-dose impact at ScaleInter=50 |
|---------|--------------------|--------------------|----------------------------------|
| 100 | Complete | Complete | 33% reduction |
| 200 | Complete | Complete | 20% reduction |
| 500 | STUCK | Complete | 9% reduction |
| 1000 | Complete | STUCK | 5% reduction |

SmaxDop=200 is the highest stable value that keeps all 6 cases running.

### Boron

| Case | Xj | Rs | Beta1 | Beta2 | Nz | Status |
|------|-----|-----|-------|-------|-----|--------|
| 2e14 20keV 900C 30min | -31.7% | -28.6% | -5.3% | -30.2% | +94.0% | Complete |
| 2e14 80keV 1000C 30min | -8.8% | -3.9% | +2.3% | -11.0% | +12.7% | Complete |
| 2e16 80keV 1100C 120min | +15.3% | -4.6% | +12.3% | +24.1% | +3.2% | Complete |
| 2e16 20keV 900C 30min | +30.1% | -83.0% | -32.3% | +26.5% | +779.0% | Complete |
| 2e16 80keV 1000C 30min | -7.5% | -43.8% | +3.1% | +29.5% | +92.5% | Complete |
| 2e14 80keV 1100C 120min | -3.6% | -3.5% | +3.5% | +0.5% | +6.6% | Complete |

### Analysis

All 6 boron cases complete (v2: 5/6, case 4 crashed; no previous template achieved 6/6).

**Low-dose accuracy improved over no-ramps version**: Adding ramp up back recovers thermal
budget. 1000C Xj improved from -18.9% to -8.8%; 1100C from -15.4% to -3.6%.

- 1000C low-dose: -8.8% (v2 was +10.4%). Closer to correct, slight under-diffusion.
- 1100C low-dose: -3.6% (v2 was +1.5%). Good.
- 900C low-dose: -31.7% (v2 was +52.5%). Still poor (frozen {311} + weak FLOOXS pair D).
- 1100C high-dose: +15.3% (v2 was >5um/invalid). First valid result for this case.
- 900C high-dose: +30.1% Xj (v2 crashed). Completes but over-diffuses (ScaleCapDop=200
  allows more TED than SmaxDop=100).
- 1000C high-dose: -7.5% (v2 was -5.4%). Good.

**High-dose profile shape**: 2e16 cases show poor Rs/Nz despite reasonable Xj. The cap
changes the profile shape (more box-like), concentrating dopant near the surface.

### Ramp down investigation

The ramp down was removed after testing multiple approaches on the 10um mesh:

| Approach | Case 3 (2e16 1100C) | Cases 1-2, 4-6 |
|----------|---------------------|-----------------|
| v2 (no cap, 10um) | CRASH (ramp down) | 5/5 complete |
| DiffDopI cap only (10um) | STUCK (ramp down) | 5/5 complete |
| Cap both DiffDopI + I self-diff, Smax=100 (5um) | STUCK (ramp down) | 5/5 complete |
| Two caps: SmaxDop=100, SmaxI=100000 (5um) | Complete (Xj>5um) | 5/5 complete |
| Two caps: SmaxDop=100, SmaxI=10000 (10um) | CRASH (CIc diverged) | 5/5 complete |
| Freeze defects during ramp down (10um) | STUCK | 5/5 complete |
| No ramp down (10um) | Complete (Xj=5.24) | 5/5 complete |
| No ramps (10um) | Complete (Xj=5.06) | 6/6 complete |

Root cause: during 1100C -> 800C cooldown, EqInter drops 5700x. ScaleInter rises
proportionally, creating extreme stiffness in the Inter and {311} clustering equations.
This affects v2 equally (no cap involved) — confirmed by testing v2 on 10um mesh, which
also crashes during ramp down for this case.

## eff_v3 template: v1 Table A-3 D + linear IV relaxation (2026-02-22)

Template: `ion_implant_eff_v3.tcl`

Physics: v1 eff template (Table A-3 D values, alpha=0.2) plus linear IV relaxation sink
kIV*EqVac*(Inter-EqInter) for boron. Only change from v1 is the IV sink + algebraic EqVac
term. P/As unchanged.

### Boron

| Case | Xj | Rs | Beta1 | Beta2 | Nz | Status |
|------|-----|-----|-------|-------|-----|--------|
| 2e14 20keV 900C 30min | +24.1% | -40.0% | +1.1% | +31.8% | +91.0% | Same as v1 |
| 2e14 80keV 1000C 30min | -3.6% | -5.6% | +2.8% | -5.8% | +12.5% | Good |
| 2e16 80keV 1100C 120min | -- | -- | -- | -- | -- | STUCK (same as v1) |
| 2e16 20keV 900C 30min | -- | -- | -- | -- | -- | CRASH |
| 2e16 80keV 1000C 30min | -- | -- | -- | -- | -- | STUCK |
| 2e14 80keV 1100C 120min | -1.3% | -4.6% | +3.8% | +3.7% | +6.8% | Good |

### Analysis

The IV relaxation sink fixes 1100C without hurting 900C/1000C:
- **1100C**: -1.3% Xj (v1 was STUCK, v5-noVac was +112%). The linear damping provides the
  I sink that was missing. Rs -4.6% is excellent.
- **1000C**: -3.6% (v1 was -2.2%). The IV sink at 1000C (tau=0.75s) slightly reduces TED
  by pulling Inter toward equilibrium faster. Marginal change, still good.
- **900C**: +24.1% (v1 was +24.3%). IV sink has negligible effect at 900C (tau=30s vs
  30min anneal = only 60 time-constants, but {311} dissolution tau=244min dominates).

However, 2e16 stability is worse than v2:
- 1100C 2e16: STUCK (v2 completed). Table A-3 D values create stiffer dopant equation
  at high concentration than FLOOXS pair D values.
- 1000C 2e16: STUCK (v2 completed with Xj=2.610 um).
- 900C 2e16: CRASH (same as v2).

### Conclusion

v3 is the best template for low-dose boron: 1000C -3.6%, 1100C -1.3%, 900C +24.1% (same
as v1). The IV sink fixes 1100C with no cost to other temperatures. But Table A-3 D values
are stiffer than FLOOXS pair D at high dose, so all 2e16 cases fail. v2's FLOOXS pair D
gave better high-dose stability but worse low-dose accuracy.

## eff_v2 template: FLOOXS pair D + linear IV relaxation (2026-02-22)

Template: `ion_implant_eff_v2.tcl`

Physics: Effective diffusivity with FLOOXS_2026 pair-derived D values (QSS-equivalent to v5
explicit pairs) and linear IV relaxation sink kIV*EqVac*(Inter-EqInter). Since D_V >> D_I
(235,000x), V equilibrates in microseconds; the relaxation term replaces the stiff Vac PDE
with a linear damping (tau: 30s at 900C, 0.75s at 1000C, 0.08s at 1100C). alpha_311=0.1
(no recalibration needed with correct D values). 3 PDEs: boron, Inter, CIc.
P/As unchanged from v1.

Changes from v1:
1. D values: Table A-3 (Ea=3.46) -> FLOOXS_2026 pair D (Ea=3.56, ~20x larger pre-exp)
2. alpha_311: 0.2 -> 0.1 (matches v5)
3. Inter equation: added kIV*EqVac*(Inter - EqInter) linear IV relaxation

### Boron

| Case | Xj | Rs | Beta1 | Beta2 | Nz | Status |
|------|-----|-----|-------|-------|-----|--------|
| 2e14 20keV 900C 30min | +52.5% | -44.6% | +3.5% | +71.3% | +92.2% | Poor |
| 2e14 80keV 1000C 30min | +10.4% | -9.3% | +3.8% | +6.2% | +12.1% | Fair |
| 2e16 80keV 1100C 120min | >5um | -- | -- | -- | -- | Completed, no junction (over-diffused) |
| 2e16 20keV 900C 30min | -- | -- | -- | -- | -- | CRASH |
| 2e16 80keV 1000C 30min | 2.610 um | -- | -- | -- | -- | Completed (no ref) |
| 2e14 80keV 1100C 120min | +1.5% | -5.3% | +3.9% | +6.1% | +6.7% | Good |

### Analysis

The IV relaxation sink works as designed:
- **1100C low-dose**: +1.5% Xj (was +112% in v5-noVac without IV sink). The linear damping
  provides the I sink that surface recomb + clustering alone couldn't at 1100C.
- **1000C**: +10.4% (v5 was +4.2%, v1 eff was -2.2%). The FLOOXS pair D values are larger
  than Table A-3 at low Poni, giving more diffusion when ScaleI is moderate.
- **900C**: +52.5% (v1 was +24.3% with alpha=0.2). With alpha=0.1 and larger D values,
  the frozen {311} clusters (tau=244min) maintain high ScaleI for the entire anneal,
  producing far too much TED. The alpha=0.2 recalibration in v1 was partially compensating
  for the smaller Table A-3 D values.

Stability improved at 2e16:
- 1100C 2e16: Completed (v5 crashed, v5-noVac stuck). Over-diffused past 5um mesh.
- 1000C 2e16: Completed with Xj=2.610 um (v5 and v5-noVac both crashed).
- 900C 2e16: Still crashes (boron equation stiffness at extreme concentration).

### Conclusion

The two-fix approach partially works: IV relaxation fixes 1100C (+112% -> +1.5%) and
improves high-dose stability (2 of 3 crash cases now complete). But the FLOOXS pair D
values with alpha=0.1 over-diffuse at 900C (+52.5%) and 1000C (+10.4%). The 900C/1000C
accuracy was better in v1 with Table A-3 D + alpha=0.2, suggesting the pair D values
need their own alpha recalibration or the two changes interact differently than expected.

## eff template: FLOOXS_2026 defects + v5 {311} (2026-02-21)

Template: `ion_implant_eff.tcl`

Physics: Effective diffusivity (Table A-3) with per-dopant defect calibration. Boron uses
FLOOXS_2026 point defects, v5 reformulated {311} clustering (Kr Ea=3.6, alpha=0.2), v5
Arrhenius KinkSite surface recombination. No Vac PDE for boron (3 PDEs: boron, Inter, CIc).
P/As use unchanged TSUPREM-4 Appendix A parameters (Table A-17/A-24/A-18/A-19).

### Key findings

1. **Vac PDE has no effect for B with FLOOXS_2026**: D_V is 235,000x larger than Table A-17,
   so V equilibrates in microseconds. V depletes to V=EqVac/ScaleI, neutralizing IV recomb.
   4-PDE test gave Xj=1.612 vs 3-PDE Xj=1.615 (identical). Safely dropped.

2. **alpha recalibration needed**: v5's alpha=0.1 gives +31% Xj at B 1000C because the eff
   model's D_table*ScaleI is 18% higher per ScaleI than v5's pair-mediated D_eff (Table A-3
   vs FLOOXS_2026 pair D parameters). alpha=0.2 corrects 1000C to -2.2%.

3. **900C/1000C conflict persists**: alpha=0.2 optimizes 1000C (-2.2%) but 900C is +24%.
   alpha=0.3 optimizes 900C (+3.0%) but 1000C drops to -16.7%. Same root cause as v5:
   frozen {311} clusters (tau=244min at 900C) maintain high ScaleI, and the eff model
   amplifies this without pair saturation.

4. **1100C still crashes**: Even with 3 PDEs (no Vac), B 2e16 1100C gets stuck in the ramp
   from {311}/Inter stiffness at extreme defect supersaturation.

### Boron alpha sweep

| alpha | B 900C Xj | B 1000C Xj |
|-------|-----------|------------|
| 0.1 | +94% | +31% |
| 0.2 | +24% | -2.2% |
| 0.25 | +11% | -10.5% |
| 0.3 | +3.0% | -16.7% |
| 0.4 | -8.0% | -- |
| 0.5 | -16% | -- |

### Boron (alpha=0.2)

| Case | Xj | Rs | Beta1 | Beta2 | Nz | Status |
|------|-----|-----|-------|-------|-----|--------|
| 2e14 80keV 1000C 30min | -2.2% | -5.9% | +2.8% | -5.0% | +12.4% | Good |
| 2e14 20keV 900C 30min | +24.3% | -40.0% | +1.1% | +31.9% | +90.9% | Poor (frozen {311}) |
| 2e16 80keV 1100C 120min | -- | -- | -- | -- | -- | STUCK |

### P/As regression (unchanged from original eff)

| Case | Xj | Rs | Status |
|------|-----|-----|--------|
| P 2e14 80keV 1000C 30min | +4.8% | -6.9% | Good |
| As 2e14 80keV 1000C 30min | -- | -- | STUCK |

As 1000C gets stuck from Vac PDE stiffness with Table A-17 parameters (same as before -
the As section was not changed).

## Spot-check comparison (no oxide)

### TSUPREM-4 reference

| Case                       | Xj (um) | Rs (ohm/sq) | Nz (cm-2) | Beta1  | Beta2  |
|----------------------------|---------|-------------|-----------|--------|--------|
| B 2e14 20keV 900C 30min   | 0.590   | 584.9       | 9.05e+17  | 0.7515 | 0.1162 |
| B 2e14 80keV 1000C 30min  | 1.230   | 303.5       | 1.62e+18  | 0.7719 | 0.2809 |
| B 2e16 80keV 1100C 120min | 4.600   | 6.2         | 1.89e+20  | 0.4483 | 0.6474 |
| P 2e14 20keV 900C 30min   | 0.290   | 387.7       | 1.14e+18  | 0.6278 | 0.0344 |
| P 2e14 80keV 1000C 30min  | 0.710   | 221.6       | 1.57e+18  | 0.7051 | 0.1196 |
| P 2e16 80keV 1100C 120min | 3.640   | 4.3         | 1.91e+20  | 0.4116 | 0.4624 |
| As 2e14 20keV 900C 30min  | 0.160   | 1532.0      | 2.45e+17  | 0.6556 | 0.0124 |
| As 2e14 80keV 1000C 30min | 0.520   | 289.4       | 1.59e+18  | 0.6177 | 0.0461 |
| As 2e16 80keV 1100C 120min| 1.790   | 6.0         | 1.75e+20  | 0.3598 | 0.2534 |

### FLOOXS v5: 5-stream + {311} + V_eq init (2026-02-21)

Template: `ion_implant_5pd_v5.tcl`

Physics: 5-stream pair diffusion with {311} interstitial clustering. V_eq initialization
for B Vac (avoids extreme IV stiffness vs v4). Per-dopant tuning: alpha_311=0.1 (B,As) /
0.6 (P), KsurfI=Arrhenius (B,As) / constant 1e11 (P). f_cluster=0.9.

All metrics integrated from surface to junction depth (see README for details).

#### Summary

| Dopant | 1000C | 900C | 1100C high-dose |
|--------|-------|------|-----------------|
| B | Xj +4%, Rs -7% | Xj +7%, Rs -36% | CRASH |
| P | Xj -1%, Rs -6% | CRASH | CRASH |
| As | Xj -5%, Rs -13% | Xj +97% | Xj +45%, Rs +8% |

Boron is in good shape at 1000C and 900C. Phosphorus works at 1000C only.
Arsenic works at 1000C; 900C and 1100C have large Xj errors. Three cases crash
from solver stiffness.

#### Boron

| Case | Xj | Rs | Beta1 | Beta2 | Nz | Status |
|------|-----|-----|-------|-------|-----|--------|
| 2e14 80keV 1000C 30min | +4.2% | -7.4% | +3.3% | -0.2% | +12.1% | Good |
| 2e14 20keV 900C 30min | +6.5% | -36.3% | -1.0% | +6.6% | +90.3% | Xj good, Rs/Nz poor |
| 2e16 80keV 1100C 120min | -- | -- | -- | -- | -- | CRASH |

Xj within 7% at both 900C and 1000C. The large Rs (-36%) and Nz (+90%) errors at 900C
are amplified from the modest Xj error: shifting the junction depth cutoff into a steep
part of the profile tail changes the integrated metrics substantially.

The 2e16 crash is a solver issue: Vac and boronVac PDEs diverge from extreme defect
supersaturation. V_eq initialization (the only change from v4) did not fix this.

#### Phosphorus

| Case | Xj | Rs | Beta1 | Beta2 | Nz | Status |
|------|-----|-----|-------|-------|-----|--------|
| 2e14 80keV 1000C 30min | -0.5% | -5.9% | +7.2% | +27.4% | +7.3% | Good |
| 2e14 20keV 900C 30min | -- | -- | -- | -- | -- | CRASH |
| 2e16 80keV 1100C 120min | -- | -- | -- | -- | -- | CRASH |

1000C is excellent on Xj and Rs. Beta2 +27% is a systematic offset from profile shape
differences between the explicit pair model and TSUPREM-4's effective diffusivity.

Only 1 of 3 cases runs. Both crashes are PI pair solver stiffness (P is modeled as an
I-only diffuser with alpha=0.6 and constant KinkSite=1e11, which creates stiff PI pair
PDEs at low temperature and high dose).

#### Arsenic

| Case | Xj | Rs | Beta1 | Beta2 | Nz | Status |
|------|-----|-----|-------|-------|-----|--------|
| 2e14 80keV 1000C 30min | -4.8% | -12.5% | +17.3% | +111.1% | +1.0% | Xj good, Beta2 poor |
| 2e14 20keV 900C 30min | +96.6% | -36.9% | +28.3% | +521.2% | -8.9% | Poor |
| 2e16 80keV 1100C 120min | +44.6% | +8.3% | +10.1% | -2.3% | -4.7% | Fair |

No crashes, but accuracy varies with temperature. 1000C Xj is within 5%. The 900C case
is far off (+97% Xj) because frozen {311} clusters maintain ScaleI=10 for the entire 30min
anneal, driving excess AsI pair diffusion. The 1100C high-dose case over-predicts Xj by 45%.

Beta2 has large systematic errors (+111% at 1000C, +521% at 900C) from profile shape
differences: the explicit pair model produces a different tail shape than TSUPREM-4's
effective diffusivity, and Beta2 is very sensitive to the concentration-weighted depth
distribution.

#### Root cause: FLOOXS vs TSUPREM-4 parameter divergence

FLOOXS and TSUPREM-4 use different Arrhenius parameters for D_I and C*_I. The product
D_I * C*_I controls intrinsic pair diffusivity:

| Temp | FLOOXS D_I*C*_I | TSUPREM-4 D_I*C*_I | Ratio |
|------|-----------------|---------------------|-------|
| 1000C | 4.2e6 | 3.1e6 | 1.4x |
| 900C | 1.9e4 | 3.3e5 | 0.06x (17x weaker) |

At 1000C, the two are similar and the model works well. At 900C, FLOOXS pair diffusivity
is 17x too weak. The frozen {311} clusters (Ea=3.6 eV, tau=244min at 900C) compensate for
this deficiency in B by maintaining elevated ScaleI, but over-compensate for As (which
doesn't need as much enhancement). This is why B and As have contradictory requirements
at 900C and no single dopant-independent parameter set works across temperatures.

### FLOOXS 5pd: TSUPREM-4 Appendix A parameters (2026-02-21)

Template: `ion_implant_5pd.tcl`

Physics: 5-stream pair diffusion with {311} interstitial clustering, using all parameters
from TSUPREM-4 manual Appendix A. Dopant-independent {311} (Table A-24, Kr Ea=3.017 eV),
surface recombination (Tables A-18/A-19), pair kinetics (Table A-6: Kf=4*pi*D_I*a0,
Kr=10 s^-1), and pair diffusivities (Table A-3). Point defects from Table A-17.

Motivation: fix the 17x D_I*C*_I mismatch between FLOOXS and TSUPREM-4 at 900C by using
TSUPREM-4's own point defect parameters directly. This also removes all per-dopant hacks
(alpha_311, KinkSite constants) that v5 required.

#### Summary

| Dopant | 1000C | 900C | 1100C high-dose |
|--------|-------|------|-----------------|
| B | Xj -58% | Xj -70% | STUCK |
| P | Xj -52% | Xj -78% | STUCK |
| As | STUCK | Xj +3025% | STUCK |

Dramatically worse than v5 in all cases. 5 of 9 cases stuck in solver (never finished
ramp phase). 4 completed cases show severe Xj errors: too little diffusion for B/P,
catastrophic over-diffusion for As at 900C. No cases are within 15% Xj.

#### Boron

| Case | Xj | Rs | Beta1 | Beta2 | Nz | Status |
|------|-----|-----|-------|-------|-----|--------|
| 2e14 80keV 1000C 30min | -58.3% | +18.9% | -4.1% | -48.0% | +16.0% | Poor |
| 2e14 20keV 900C 30min | -69.8% | -13.1% | -15.1% | -70.1% | +99.6% | Poor |
| 2e16 80keV 1100C 120min | -- | -- | -- | -- | -- | STUCK |

#### Phosphorus

| Case | Xj | Rs | Beta1 | Beta2 | Nz | Status |
|------|-----|-----|-------|-------|-----|--------|
| 2e14 80keV 1000C 30min | -52.0% | +16.3% | -0.2% | -34.7% | +13.7% | Poor |
| 2e14 20keV 900C 30min | -78.3% | +78.6% | -8.3% | -74.3% | -20.2% | Poor |
| 2e16 80keV 1100C 120min | -- | -- | -- | -- | -- | STUCK |

#### Arsenic

| Case | Xj | Rs | Beta1 | Beta2 | Nz | Status |
|------|-----|-----|-------|-------|-----|--------|
| 2e14 80keV 1000C 30min | -- | -- | -- | -- | -- | STUCK |
| 2e14 20keV 900C 30min | +3025% | -82.1% | +45.5% | +19019% | +6.3% | Diverged |
| 2e16 80keV 1100C 120min | -- | -- | -- | -- | -- | STUCK |

As 900C: no junction found (profile spread beyond 5um mesh). As 1000C: solver stuck
in ramp phase with arsenicVac stiffness (330K solver iterations, still in ramp).

#### Root cause: pair fraction too low with TSUPREM-4 Table A-6 kinetics

The TSUPREM-4 pair kinetics model (Table A-6) uses diffusion-limited formation
(Kf = 4*pi*D_I*a0) and constant dissolution (Kr = 10 s^-1). In an explicit pair PDE,
the equilibrium pair fraction is Pair/Sub = Kf/Kr * I = (4*pi*D_I*a0/10) * I.

| Temp | Kf (cm^3/s) | C*_I (cm^-3) | Pair frac (eq) | Pair frac (ScaleI=10) |
|------|-------------|-------------|----------------|----------------------|
| 900C | 2.0e-17 | 1.2e15 | 0.3% | 2.5% |
| 1000C | 6.9e-17 | 1.6e16 | 10.8% | 100% (saturated) |
| 1100C | 2.0e-16 | 1.4e17 | 100% (sat) | 100% (saturated) |

At 900C, only 0.3% of dopant is in pairs at equilibrium (2.5% with TED). Since dopant
only moves via pair diffusion, the effective dopant diffusivity is 40x too low. At 1000C,
TED initially saturates pairs (100%), but after {311} dissolve and ScaleI drops to 1,
only 10.8% remains in pairs -- the post-TED phase has 10x too little diffusion.

The pair model D_eff has an intrinsic ceiling: D_pair * min(Kf/Kr * I, 1). At 1000C
with ScaleI=10, D_eff_B = D_pair = 1.44e-14 cm^2/s, matching the expected intrinsic
B diffusivity (1.53e-14, ratio = 0.945). But this only holds while ScaleI >= ~10. Once
TED ends, D_eff drops to 10.8% of the correct value. The effective diffusivity model
(used by TSUPREM-4 internally) does not have this saturation: D_eff = D_intrinsic * ScaleI
at all times.

For As at 900C, the opposite problem: D_AsI has a massive pre-exponential
(DIX.0 = 1.37e7 cm^2/s) giving D_AsI = 2.3e-8 cm^2/s at 900C. Even at 2.5% pair
fraction, the effective D = 5.7e-10 cm^2/s, enough to diffuse 15um in the anneal time.

#### Conclusion

The TSUPREM-4 Appendix A pair parameters (Tables A-3, A-6) are designed for TSUPREM-4's
internal effective diffusivity formulation, not for an explicit pair PDE model. In TSUPREM-4,
the pair kinetics produce a quasi-steady-state pair fraction that multiplies the substitutional
dopant concentration directly: D_eff * grad(C * fermi). The explicit pair model used in
FLOOXS solves separate pair PDEs where the pair concentration is a small fraction of the
total dopant (governed by Kf/Kr), producing fundamentally different diffusion rates.

The v5 template with FLOOXS_2026 binding-based kinetics (Krate * Sub * I vs Krate * Pair / Binding)
avoids this issue because the binding constant controls the pair equilibrium differently, giving
pair fractions consistent with the expected effective diffusivity.

### Kr Ea sweep (2026-02-21)

Tested {311} dissolution Ea = {3.0, 3.2, 3.4, 3.6} eV for B at 900C and 1000C,
and Ea = {3.0, 3.6} for As at 900C and 1000C.

**B 2e14 80keV 1000C 30min (Xj error):**

| Ea   | tau(1000C) | Xj    | Rs    | Notes |
|------|-----------|-------|-------|-------|
| 3.0  | 4 s       | -0.1% | -6.3% | Best 1000C |
| 3.2  | 23 s      | +0.7% | -6.4% | |
| 3.4  | 2.4 min   | +1.3% | -6.6% | |
| 3.6  | 15 min    | +4.2% | -7.4% | Current |

1000C is robust to Ea: Xj varies only 4% across the full sweep because clusters dissolve
quickly at all Ea values (tau ranges from 4s to 15min, all short vs 30min anneal).

**B 2e14 20keV 900C 30min (Xj error):**

| Ea   | tau(900C) | Xj     | Rs     | Notes |
|------|----------|--------|--------|-------|
| 3.0  | 39 s     | -32.9% | -28.1% | Clusters dissolve, too little TED |
| 3.2  | 4.7 min  | -30.2% | -28.5% | |
| 3.4  | 34 min   | -25.7% | -29.6% | |
| 3.6  | 244 min  | +6.5%  | -36.3% | Clusters frozen, compensates weak pair D |

900C is highly sensitive to Ea. Only Ea=3.6 (frozen clusters for the entire anneal)
produces enough TED to match TSUPREM-4. Lower Ea values release clusters too early,
exposing the 17x-weak pair diffusivity.

**As Ea sweep (Xj error):**

| Case  | Ea=3.0 | Ea=3.6 | Notes |
|-------|--------|--------|-------|
| 1000C | -15.1% | -4.8%  | Lower Ea worsens 1000C |
| 900C  | -25.5% | +96.6% | Ea=3.0 dramatically improves 900C |

B and As have contradictory Ea requirements at 900C. B needs Ea=3.6 (frozen clusters);
As needs lower Ea (frozen clusters cause excess AsI diffusion). No single
dopant-independent Ea satisfies both.

### FLOOXS v4: 5-stream + {311} clustering (2026-02-21)

Template: `ion_implant_5pd_v4.tcl`
Parameters: Same as v5 but with V=1.0 init for B (instead of V_eq). Results identical
to v5 for all non-crashing cases (V_eq init only affects high-dose stability).

See v5 table above for results.
Per-dopant: alpha_311=0.1 (B,As) / 0.6 (P), KsurfI=Arrhenius (B,As) / constant 1e11 (P).
Beta1/Beta2 errors have ~5% systematic offset from mobility model differences (Python vs MATLAB).

### FLOOXS v3: 5-stream with effective +n model (2026-02-20)

Template: `ion_implant_ted.tcl`
Parameters: 5-stream pair diffusion (B, P, As), Hobler-Moroz f_pl initialization

| Case                       | Xj (um) | Xj error | f_pl |
|----------------------------|---------|----------|------|
| B 2e14 80keV 1000C 30min  | 0.850   | -31%     | 1.05 |
| As 2e14 80keV 1000C 30min | 0.318   | -39%     | 2.29 |

As improved from 0.268 (+1 model) to 0.318 (+n model), but still -39% vs target 0.520.
B essentially unchanged (0.846 -> 0.850) as expected (f_pl ~ 1.05).

### FLOOXS v2: QSS/Fermi hybrid (2026-02-18)

Template: `ion_implant_ted_qss.tcl`
Parameters: B QSS (uncapped ScaleInter), P fermi-only, As QSS (ScaleInter, ScaleV=1)

| Case                       | Xj (um) | Xj error |
|----------------------------|---------|----------|
| B 2e14 80keV 1000C 30min  | 1.102   | -10%     |
| P 2e14 80keV 1000C 30min  | 0.547   | -23%     |
| P 2e16 80keV 1100C 120min | 4.606   | +27%     |
| As 2e14 80keV 1000C 30min | 0.237   | -54%     |

### FLOOXS v1: QSS all dopants (2026-02-17, superseded)

Template: `ion_implant_ted_qss.tcl` (earlier version)
Parameters: B uncapped, P fI=0.2 Smax=100, As fermi-only

| Case                       | Xj (um) | Xj error |
|----------------------------|---------|----------|
| B 2e14 80keV 1000C 30min  | 1.102   | -10%     |
| P 2e14 80keV 1000C 30min  | 0.588   | -17%     |
| As 2e14 80keV 1000C 30min | 0.237   | -54%     |

P QSS ScaleInter cap sweep (2e14 80keV 1000C 30min, fI=0.2):

| Smax  | Xj (um) | Xj error |
|-------|---------|----------|
| 100   | 0.588   | -17%     |
| 200   | 0.600   | -15%     |
| 500   | 0.617   | -13%     |
| 1000  | 0.631   | -11%     |
| 5000  | 0.660   | -7%      |
| 10000 | 0.669   | -6%      |

P QSS diverged at 2e16 dose for all Smax values tested (Inter PDE stiffness).

## v5-noVac template (ion_implant_5pd_v5_noVac.tcl)

v5 with Vac PDE, BV pairs, and IV recombination removed for boron only.
Hypothesis: Vac/BV stiffness causes v5 B 1100C crash; V pathway is only ~5% of B diffusion.
P/As unchanged from v5.

### Boron results

| Case | Xj (um) | Xj ref | Xj error | Status |
|------|---------|--------|----------|--------|
| B 2e14 20keV 900C 30min | 0.624 | 0.590 | +5.8% | OK |
| B 2e14 80keV 1000C 30min | 1.354 | 1.230 | +10.1% | Marginal (v5 was +4.4%) |
| B 2e14 80keV 1100C 120min | 5.000 | 2.360 | +111.9% | FAIL (hit mesh boundary) |
| B 2e16 20keV 900C 30min | - | - | - | CRASH (boronInt/Inter stiffness) |
| B 2e16 80keV 1000C 30min | - | - | - | CRASH (boronInt/Inter stiffness) |
| B 2e16 80keV 1100C 120min | - | - | - | STUCK (boronInt pair stiffness) |

### Analysis

Removing Vac PDE is worse than v5 in every dimension:
1. **1100C over-diffusion**: Without IV recombination, interstitials have only surface recomb + clustering as sinks. At 1100C with long anneal (120min), ScaleI stays too high -> massive TED (+112%).
2. **All 2e16 cases fail**: v5 only failed 2e16 1100C. Without Vac PDE, the BI pair reaction alone becomes stiff at high dose concentrations.
3. **1000C accuracy degrades**: +10.1% vs v5's +4.4%. Missing BV pairs change the pair fraction equilibrium.

### Conclusion

The Vac PDE cannot be simply removed. IV recombination is a critical interstitial sink, especially at 1100C. The stiffness source is not just Vac - it's the pair reactions themselves at high dose. The BI pair stiffness at 2e16 was masked in v5 by the Vac/boronVac stiffness crashing first.
