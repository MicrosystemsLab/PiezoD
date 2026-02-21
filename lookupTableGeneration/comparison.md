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
bash scripts/gen_test.sh <dopant> <dose> <energy> <temp> <time> templates/ion_implant_ted_qss.tcl > simulations/test.tcl
MSYS_NO_PATHCONV=1 docker compose run --rm flooxs /work/test.tcl
```

The template currently outputs Xj only. Rs, Nz, beta1, beta2 require post-processing
the full doping profile through `postProcessTables.m` (not yet implemented for FLOOXS).

## Template parameters

`ion_implant_ted_qss.tcl` - QSS/Fermi hybrid model:
- B: QSS with explicit interstitial transport (2 solved variables: B + Inter)
  - D_eff = D_I(Poni) * ScaleInter + D_V(Poni)
  - ScaleInter = I/I* (uncapped), D_V floor stabilizes solver
  - Surface recombination only (no bulk I-V recombination, no BIC clustering)
- P: Fermi-level dependent diffusivity, no TED (1 solved variable)
  - D_eff = D0 + Dn*Noni + Dnn*Noni^2
  - Inter PDE too stiff at high dose (diverges at 2e16)
- As: QSS with explicit interstitial transport, ScaleV=1 (2 solved variables: As + Inter)
  - D_eff = D_I * ScaleInter + D_V(Noni)
  - ScaleInter = I/I* (uncapped), V at equilibrium

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

### FLOOXS v4: 5-stream + {311} clustering (2026-02-21)

Template: `ion_implant_5pd_v4.tcl`
Parameters: 5-stream pair diffusion + {311} interstitial clustering, corrected KinkSite, Hobler-Moroz f_pl.
Per-dopant hacks: alpha_311=0.1 (B,As) / 0.6 (P), KsurfI=Arrhenius (B,As) / constant 1e11 (P).

All metrics integrated from surface to junction depth (see README for details).

| Case | Xj | Rs | Beta1 | Beta2 | Nz | Status |
|------|-----|-----|-------|-------|-----|--------|
| B 2e14 80keV 1000C 30min | +3.8% | -7.9% | +3.1% | -0.8% | +13.3% | OK |
| P 2e14 80keV 1000C 30min | -0.5% | -5.9% | +7.2% | +27.4% | +7.3% | OK |
| As 2e14 80keV 1000C 30min | -4.8% | -12.5% | +17.3% | +111.1% | +1.0% | OK |
| B 2e14 20keV 900C 30min | +7.0% | -36.4% | -0.9% | +7.2% | +90.3% | OK |
| P 2e14 20keV 900C 30min | -- | -- | -- | -- | -- | CRASH: PI pair stiffness |
| As 2e14 20keV 900C 30min | +96.6% | -36.9% | +28.3% | +521.2% | -8.9% | OK |
| B 2e16 80keV 1100C 120min | -- | -- | -- | -- | -- | CRASH: Vac/boronVac diverge |
| P 2e16 80keV 1100C 120min | -- | -- | -- | -- | -- | CRASH: PI pair stiffness |
| As 2e16 80keV 1100C 120min | +44.6% | +8.3% | +10.1% | -2.3% | -4.7% | OK |

Calibrated cases (2e14 80keV 1000C 30min) within 5% on Xj for all dopants.
High-dose (2e16) crashes or over-predicts -- solver stiffness from extreme defect supersaturation.
900C cases over-predict: f_cluster=0.9 leaves 10% of interstitials free at initialization,
and {311} dissolution is frozen at 900C (Kr time constant ~250 min >> 30 min anneal),
so initial free I drives excessive TED. In TSUPREM-4, CL.INI.F=1.0 (all I initially clustered).
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
