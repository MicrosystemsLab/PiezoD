# Implementation Notes

Lookup table generation uses FLOOXS, an open-source TCAD process simulator (Tcl-based, C++ solver), running in Docker. Our target is to match TSUPREM-4 (commercial TCAD) reference results, particularly for Transient Enhanced Diffusion (TED). See [../comparison.md](../comparison.md) for validation data.

## Current model: 5-stream + {311} clustering (`ion_implant_5pd.tcl`)

5-stream pair diffusion (explicit dopant-defect pairs) with {311} interstitial clustering for all dopants. All physics parameters from TSUPREM-4 manual Appendix A.

Per-dopant physics:
- **Boron**: 6 solved variables (B, BI, BV, Inter, Vac, CIc)
  - B-I pairs: DIX (neutral) + DIP (positive) charge states
  - B-V pairs: DVX (neutral) + DVP (positive) charge states
  - I-V bulk recombination, surface recombination, {311} clustering
- **Phosphorus**: 4 solved variables (P, PI, Inter, CIc)
  - P-I pairs only: DIX (neutral) + DIM (negative) + DIMM (double-negative)
  - No Vac equation (I-only diffuser, avoids kIV stiffness)
- **Arsenic**: 6 solved variables (As, AsI, AsV, Inter, Vac, CIc)
  - As-I pairs: DIX (neutral) + DIM (negative)
  - As-V pairs: DVX (neutral) + DVM (negative, electron-enhanced)
  - I-V bulk recombination, surface recombination, {311} clustering

### Parameter sources

All parameters trace to TSUPREM-4 manual Appendix A with specific table references:

| Parameter | Table | PDF page | Notes |
|-----------|-------|----------|-------|
| Point defect D, C*, charges | A-17 | 752-753 | D and C* identical for I and V |
| {311} clustering Kfc, Kr | A-24 | 755 | CF=0.9398 approximated as 1.0 |
| Surface recombination Ksurf | A-18/A-19 | 753-754 | Same for all dopants |
| Pair kinetics Kf, Kr | A-6 | 748 | Kr=10 s^-1 (constant, all dopants) |
| Pair diffusivities | A-3 | 746-747 | um^2/min except where noted cm^2/s |
| Effective +n (Hobler-Moroz) | A-46 | 761 | f_pl model |
| Solid solubility | - | - | FLOOXS_2026 (TSUPREM-4 uses tabulated) |

### Pair reaction structure

Pair formation and dissolution follow TSUPREM-4 Eq 3-61 with Table A-6 kinetics:

```
React = Kf * Sub * Defect - Kr * Pair
```

- Formation: `Kf = 4*pi*D_defect*a0` (diffusion-limited capture)
- Dissolution: `Kr = 10 s^-1` (Table A-6: R.I.S = R.V.S = 10, E.I.S = E.V.S = 0)

Pair diffusivities are pre-computed at anneal temperature from Table A-3 and used as Tcl scalars in the PDE. This avoids Arrhenius expressions in the solver while maintaining correct physics at the anneal dwell temperature.

### {311} interstitial clustering

TSUPREM-4 1-moment model (Eq 3-319, manual p. 3-73):

```
dCIc/dt = Kfc * (I/I*) * (CIc + I) - Kr * CIc
```

Parameters from Table A-24: Kfc=Arr(5.207e14, 3.774), Kr=Arr(9.431e13, 3.017). Both are pre-computed at anneal temperature. CL.CF=0.9398 is approximated as 1.0 (FLOOXS lacks pow() in Alagator).

Initialization: 90% of excess I starts in clusters (f_cluster = 0.9), 10% free. This avoids the stiff initial transient where ScaleInter >> 1.

### Surface recombination

Direct KSURF from Tables A-18/A-19, same for all dopants:
- Interstitial: KSURF.0=1.4e-6, KSURF.E=-1.75
- Vacancy: KSURF.0=4.0e-11, KSURF.E=-1.75

Pre-computed at anneal temperature and applied as scalar boundary condition coefficients.

### Effective +n model (Hobler-Moroz)

TSUPREM-4's `damage` keyword activates the Hobler-Moroz effective +n model (manual Eq 3-534), which accounts for the fact that vacancies from Frenkel pairs preferentially recombine at the surface, leaving a net excess of more than 1 interstitial per implanted ion:

```
f_pl = 1 + D.PHDF * m^D.PME * E^(D.PLF * m^D.PLME)
```

| Ion | Energy (keV) | f_pl |
|-----|-------------|------|
| B   | 80          | 1.05 |
| P   | 80          | 1.35 |
| As  | 20          | 2.78 |
| As  | 80          | 2.29 |

Initialization: `Inter = (1 - f_cluster) * f_pl * dopant`, `CIc = f_cluster * f_pl * dopant`.

### Parameters NOT changed from previous model

| Parameter | Source | Reason |
|-----------|--------|--------|
| Solubility | FLOOXS Arrhenius | TSUPREM-4 Table A-14 is tabulated, not Arrhenius. Irrelevant at our doses (peak ~1e19 < Css ~1e20). |
| f_cluster | 0.9 | TSUPREM-4 CL.INI.A=true (auto-init) not reproducible in FLOOXS. 0.9 is physically reasonable. |
| ni | Law et al. | Already matches TSUPREM-4 Table A-58. |
| f_pl | Hobler-Moroz | Already matches TSUPREM-4. |

### Pre-computation strategy

Temperature-dependent parameters are pre-computed at the anneal temperature as Tcl scalars. This:
1. Avoids Arrhenius expressions in the PDE (solver stability)
2. Uses TSUPREM-4 values directly without FLOOXS pdb machinery
3. Is accurate during the dwell (dominant diffusion phase)
4. Is approximate during the 45-min ramp (minor contribution)

Dynamic (temperature-dependent via pdbDelayDouble): C*, D_defect, charge states -- used in the diffusion/equilibrium equations where temperature variation matters.

Pre-computed (scalar at anneal temp): pair diffusivities, reaction rates, surface recombination, {311} rates -- used as coefficients where anneal-temp accuracy is sufficient.

## Gap analysis: FLOOXS 5pd vs TSUPREM-4

| Physics | TSUPREM-4 ref | 5pd template | Impact |
|---------|--------------|-------------|--------|
| Pair diffusion (M, N) | Eq 3-61/62 | Yes | Core |
| I transport | Eq 3-280 | Yes | Core |
| V transport | Eq 3-281 | Yes (B, As) | I regulator |
| Bulk I-V recomb | Eq 3-286 | Yes | Dominant I sink |
| Surface recomb | Eq 3-293 | Yes (KSURF from Tables A-18/A-19) | I boundary sink |
| {311} clustering | Eq 3-319 | Yes (direct Kfc/Kr, CF=1 approx) | Critical |
| Effective +n init | Eq 3-534 | Yes | Initialization |
| Small defect clusters | Eq 3-345 | No | Moderate |
| Dopant-assisted recomb | Eq 3-65/66 | No | Minor |
| Interstitial traps | Eq 3-313 | No | Unknown |
| BIC (B only) | Eq 3-122 | No | B-specific |

## TSUPREM-4 reference model: `pd.5str` with `ddc.full`

Our TSUPREM-4 template uses `method vertical ddc.full pd.5str` with zero parameter overrides (all defaults). Key equations:

### 5-stream dopant pair equations

The pair concentrations M (dopant-I) and N (dopant-V) are solved explicitly (Eq 3-61, 3-62):

```
dM/dt = -div(J_m) + (G_m - R_m) - (G_mv - R_mv) - R_DDCM
dN/dt = -div(J_n) + (G_n - R_n) - (G_ni - R_ni) - R_DDCN
```

### Point defect equations (Eq 3-280, 3-281)

```
dI/dt = -div(D_I transport) - pair_formation - bulk_IV_recomb - {311}_clustering - BIC
dV/dt = -div(D_V transport) - pair_formation - bulk_IV_recomb - {311}_dissolution
```

### {311} clustering (Eq 3-319)

```
dC_Ic/dt = R_IcI - R_IcV
R_IcI = K_fc * (I/I*)^ISFC * (C_Ic + alpha*I)^CF    [growth, Eq 3-323]
R_IcV = K_r * C_Ic^CR                                 [dissolution, Eq 3-329]
```

With defaults (ISFC=CF=CR=1, alpha=1, K_fi=0): `dC_Ic/dt = K_fc*(I/I*)*(C_Ic+I) - K_r*C_Ic`.

## Historical notes

### Why QSS fails for As

As at 80 keV produces a narrow implant (Rp~0.05 um) giving ScaleInter_peak ~ 2.9e8 at 800C. The ScaleInter gradient is ~8.1e12 per um, too steep for the Newton solver. Without bulk I-V recombination (QSS with ScaleV=1), there is no natural I regulator. The 5-stream model with explicit V transport provides the correct physics: V depletion acts as a natural brake on I-V recombination, bounding ScaleI to moderate levels.

### FLOOXS built-in model procs

FLOOXS has built-in Tcl procs (DopantPair, DopantReact) that implement TED physics. We cannot use them because the C++ expression parser crashes at solve time with Arrhenius expressions set up through the model proc pathway. The exact same expressions work in hand-written templates. All templates use hand-written equations.

### {311} solver stability (historical)

Early attempts used the reformulated {311} equation: `dCIc/dt = Kr * (alpha * ScaleInter * (CIc + Inter) - CIc)` to avoid cancellation stiffness from two separate Arrhenius terms. With the full TSUPREM-4 parameter set (corrected CEQUIL, KSURF, pair kinetics), the direct form `Kfc*ScaleInter*(CIc+Inter) - Kr*CIc` with pre-computed rate constants works. If solver stability issues arise, CL.KFC.0 is the single tuning knob (scales clustering strength while preserving correct temperature dependence).
