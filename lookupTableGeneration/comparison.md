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
