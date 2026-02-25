# Effective Diffusivity Template Plan (ion_implant_eff.tcl)

## Goal

Reproduce pd.5str (`method vertical ddc.full pd.5str`) behavior using effective
diffusivities instead of explicit pair PDEs. The pair coupling terms from
Eq 3-280/281 are preserved via the Eq 3-89 reformulation.

The 5pd template (explicit pair PDEs with Table A-6 kinetics) failed because
pair fraction saturation at ~0.3% (900C) to ~11% (1000C) made the pair PDEs
numerically stiff and physically inaccurate. The effective D approach avoids
pair PDEs entirely while retaining the coupling between dopant diffusion and
point defect transport.

## pd.5str Equation Mapping

### I equation (Eq 3-280)

```
d(I + Ic + Idds)/dt = -div[D_I * I * (grad(I/I*) + (I/I*) * grad(ln(I*_i)))]
                       - (G_m - R_m)         [I-pair formation/dissolution]
                       - (G_ni - R_ni)       [V-pair kick-out to I]
                       - R_b                 [IV bulk recombination]
                       - R_t                 [interstitial traps]
                       - R_slcI              [small {311} clustering]
                       - R_llcI              [{311} large clustering]
                       - R_l                 [dislocation loops]
                       - R_dd                [DDC: dopant-defect]
                       - sum(R_DDC,...)       [DDC reactions]
```

Term-by-term status with default TSUPREM-4 parameters:

| pd.5str term | Status | Implementation |
|---|---|---|
| d(I)/dt | Active | `ddt(Inter)` |
| d(Ic)/dt | Active | Separate CIc PDE (Eq 3-319) |
| d(Idds)/dt | = 0 | DDC params all zero (Table A-16) -> Idds = 0 |
| -div[D_I\*I\*grad(I/I\*)] | Active | `DiffI*EqInter*grad(ScaleInter)` |
| (I/I\*)\*grad(ln(I\*_i)) | = 0 | T spatially uniform -> I\*_i constant in space |
| -(G_m - R_m) | Active | `DiffDopI*grad(dopant*eta)` (derivation below) |
| -(G_ni - R_ni) | = 0 | Default K_ni = 0 (Eq 3-66) |
| R_b | Active (B, As) | `kIV*(Inter*Vac - EqInter*EqVac)` (Eq 3-286) |
| R_t | = 0 | Default TRAP.CON = 0 (Table A-21) |
| R_slcI | = 0 | Default: not CL.FULL |
| R_llcI | Active | `Kfc*ScaleI*(CIc+alpha*I)^CF - Kr*CIc` (Eq 3-323) |
| R_l | = 0 | Default: dislocation loops off |
| R_dd | = 0 | DDC.F.0 = 0 (Table A-16) -> no DDC reactions |
| sum(R_DDC,...) | = 0 | DDC.F.0 = 0 (Table A-16) -> no DDC reactions |

### V equation (Eq 3-281)

| pd.5str term | Status | Implementation |
|---|---|---|
| d(V)/dt | Active (B, As) | `ddt(Vac)` |
| d(Vc)/dt | = 0 | No V clustering: default K_rv = K_fv = 0 |
| -div[D_V\*V\*grad(V/V\*)] | Active | `DiffV*EqVac*grad(ScaleVac)` |
| -(G_n - R_n) | Active (B, As) | `DiffDopV*grad(dopant*eta)` |
| -(G_mv - R_mv) | = 0 | Default K_mv = 0 (Eq 3-65) |
| R_b | Active (B, As) | Same as I equation |
| R_slcV, R_llcV | = 0 | Default K_rv = K_fv = 0 (Eq 3-326) |
| R_DDCN,... | = 0 | DDC params all zero |

### DDC verification (Table A-16)

Our TSUPREM-4 reference uses `method vertical ddc.full pd.5str`, but Table A-16
shows ALL DDC parameters default to zero:

- DDC.T.0 = 0, DDC.T.E = 0, DDC.F.0 = 0, DDC.F.E = 0
- DDCF.D.N = 0, DDCF.I.N = 0, DDCF.N.N = 0
- DDCR.N.N = 0, DDCR.I.N = 0
- IFRACM = 0, IFRACS = 0

With DDC.F.0 = 0, the BIC formation rate K_ddF (Eq 3-124) = 0. No dopant-defect
clusters form. The `ddc.full` flag enables the framework but produces zero rates
with default parameters. DDC is NOT a deviation.

## Pair coupling derivation: (G_m - R_m) -> effective D flux

Starting from the pd.5str pair PDE (Eq 3-61):

```
dM/dt = -div(J_m) + (G_m - R_m) - (G_mv - R_mv) - R_DDCM
```

Rearranging (Eq 3-89, exact algebraic identity):

```
G_m - R_m = dM/dt + div(J_m) + R_DDCM + (G_mv - R_mv)
```

With default parameters (K_mv = 0, DDC params = 0):

```
G_m - R_m = dM/dt + div(J_m)
```

Approximation: dM/dt ~ 0. Pair fraction alpha_m ~ 0.3% at 900C, ~11% at 1000C
per Table A-6 analysis. |dM/dt| ~ alpha * |dC/dt|, negligible compared to
dominant I sinks ({311} clustering and surface recombination).

Therefore:

```
G_m - R_m ~ div(J_m)
```

J_m is the I-pathway effective dopant flux. From Eq 3-31 with M/M' = I/I*
(Eq 3-54 with K_mv = 0, which reduces to ScaleI):

```
For acceptors: J_m = -(D_I_dop * ScaleI / Poni) * grad(C * Poni)
For donors:    J_m = -(D_I_dop * ScaleI / Noni) * grad(C * Noni)
```

So the coupling term in the I equation is:

```
-(G_m - R_m) ~ -div(J_m) = div(DiffDopI * grad(C * eta))
```

where DiffDopI = D_I_dop * ScaleI / eta and eta = Poni (acceptors) or Noni (donors).

Physical meaning: I is consumed where dopant pairs form and flow out (concentration
peak), and released where pairs arrive and dissolve (profile tails). Same derivation
applies to -(G_n - R_n) for the V equation with DiffDopV.

## Equations (FLOOXS Alagator form)

### Boron (acceptor, I+V pathways)

```
# Effective D components (Table A-3, Eq 3-45/47)
DiffDopI = (DIX + DIP * Poni) * ScaleInter / Poni
DiffDopV = (DVX + DVP * Poni) * ScaleVac / Poni
DiffDop  = DiffDopI + DiffDopV

# Dopant (Eq 3-30: dC/dt = -div(J_m + J_n))
ddt(boron) - DiffDop * grad(boron * Poni) = 0

# Inter (Eq 3-280)
ddt(Inter)
  - DiffI * EqInter * grad(ScaleInter)   [defect diffusion]
  - DiffDopI * grad(boron * Poni)         [-(G_m - R_m), pair coupling]
  + kIV * (Inter*Vac - EqInter*EqVac)     [R_b, IV recombination]
  + R_clust                                [R_llcI, {311} clustering]
  = 0

# Vac (Eq 3-281)
ddt(Vac)
  - DiffV * EqVac * grad(ScaleVac)        [defect diffusion]
  - DiffDopV * grad(boron * Poni)          [-(G_n - R_n), pair coupling]
  + kIV * (Inter*Vac - EqInter*EqVac)     [R_b]
  = 0

# CIc (Eq 3-319)
ddt(CIc) - R_clust = 0
```

### Phosphorus (donor, I-only)

```
DiffDopI = (DIX + DIM * Noni + DIMM * Noni * Noni) * ScaleInter / Noni
DiffDop  = DiffDopI

# Dopant
ddt(phosphorus) - DiffDop * grad(phosphorus * Noni) = 0

# Inter (no R_b since no Vac PDE)
ddt(Inter)
  - DiffI * EqInter * grad(ScaleInter)
  - DiffDopI * grad(phosphorus * Noni)
  + R_clust
  = 0

# No Vac equation (P is I-only diffuser)

# CIc
ddt(CIc) - R_clust = 0
```

### Arsenic (donor, I+V pathways)

```
DiffDopI = (DIX + DIM * Noni) * ScaleInter / Noni
DiffDopV = (DVX + DVM * Noni) * ScaleVac / Noni
DiffDop  = DiffDopI + DiffDopV

# Dopant
ddt(arsenic) - DiffDop * grad(arsenic * Noni) = 0

# Inter
ddt(Inter)
  - DiffI * EqInter * grad(ScaleInter)
  - DiffDopI * grad(arsenic * Noni)
  + kIV * (Inter*Vac - EqInter*EqVac)
  + R_clust
  = 0

# Vac
ddt(Vac)
  - DiffV * EqVac * grad(ScaleVac)
  - DiffDopV * grad(arsenic * Noni)
  + kIV * (Inter*Vac - EqInter*EqVac)
  = 0

# CIc
ddt(CIc) - R_clust = 0
```

### {311} clustering (Eq 3-319, 3-323)

```
R_llcI = Kfc * (I/I*)^CL.IFC / (I*)^(CL.ISFC-CL.IFC) * (CIc + alpha*I)^CL.CF
       - Kr * CIc^CL.CR
```

Default parameters (Table A-24):

| Parameter | Value | Notes |
|---|---|---|
| CL.KFC.0 | 5.207e14 | Forward rate prefactor |
| CL.KFC.E | 3.774 | Forward rate Ea (eV) |
| CL.KR.0 | 9.431e13 | Reverse rate prefactor |
| CL.KR.E | 3.017 | Reverse rate Ea (eV) |
| CL.CF | 0.9398 | Forward exponent (approx as 1.0) |
| CL.CR | 1.0 | Reverse exponent (exact) |
| CL.IFC | 1 | ScaleI exponent |
| CL.ISFC | 1 | EqI exponent |
| CL.KFCI | 1.0 | alpha (confirmed from Table A-24) |

With CL.IFC = CL.ISFC = 1: I^IFC / I*^ISFC = I/I* = ScaleI.
With CL.CF ~ 1 (FLOOXS lacks pow()): (CIc + alpha*I)^CF ~ (CIc + alpha*I).

FLOOXS form:
```
R_clust = Kfc * ScaleInter * (CIc + 1.0 * Inter) - Kr * CIc
```

### Charge neutrality (Eq 3-42)

```
Charge = N_D - N_A (net ionized dopant)
Noni = 0.5 * (Charge + sqrt(Charge^2 + 4*ni^2)) / ni
Poni = 0.5 * (-Charge + sqrt(Charge^2 + 4*ni^2)) / ni
```

Uses total active dopant (not Sub = C - M - N), since pair concentrations M, N
are not tracked. This is exact when M + N << C, which holds because pair fractions
are small (Table A-6 analysis).

### Dopant equation origin (Eq 3-30, 3-31/32, 3-45, 3-54/55)

The dopant equation is Eq 3-30: dC/dt = -div(J_m + J_n).

The fluxes J_m, J_n (Eq 3-31/32) use diffusivities D_m, D_n from Eq 3-45 with
charge-state components D_mk, D_nk from Table A-3 (Eq 3-47), and enhancement
factors M/M', N/N' from Eq 3-54/55.

With default K_mv = K_ni = 0 and PAIR.SAT = false:

```
M/M' = I/I* = ScaleI     (Eq 3-54 simplifies)
N/N' = V/V* = ScaleV     (Eq 3-55 simplifies)
```

The electric field term (Eq 3-41): qE/kT = -grad(ln(eta)), where eta = n/ni = Noni.

For acceptors (z_s = -1), combining gradient and field into a single term:
```
J_m = -D_m * ScaleI * grad(C * Poni) / Poni
J_n = -D_n * ScaleV * grad(C * Poni) / Poni
```

For donors (z_s = +1):
```
J_m = -D_m * ScaleI * grad(C * Noni) / Noni
J_n = -D_n * ScaleV * grad(C * Noni) / Noni
```

## Parameters (TSUPREM-4 Appendix A)

All identical to 5pd template. Pre-computed at anneal temp as Tcl scalars.

### Effective diffusivities (Table A-3)

Units in table: um^2/min. Conversion: * 1.667e-10 -> cm^2/s.
Exceptions noted where table gives cm^2/s directly.

Boron (all Ea = 3.46 eV):

| Parameter | D0 (cm^2/s) | Table A-3 value | Charge state |
|---|---|---|---|
| DIX | 0.03517 | 2.11e8 um^2/min | neutral (I path) |
| DIP | 0.6835 | 4.10e9 um^2/min | positive (I path) |
| DVX | 1.850e-3 | 1.11e7 um^2/min | neutral (V path) |
| DVP | 0.03601 | 2.16e8 um^2/min | positive (V path) |

Phosphorus (no V pathway):

| Parameter | D0 (cm^2/s) | Ea (eV) | Table A-3 value |
|---|---|---|---|
| DIX | 3.85 | 3.66 | 2.31e10 um^2/min |
| DIM | 4.44 | 4.00 | 2.664e10 um^2/min |
| DIMM | 44.2 | 4.37 | 2.652e11 um^2/min |

Arsenic:

| Parameter | D0 (cm^2/s) | Ea (eV) | Table A-3 value |
|---|---|---|---|
| DIX | 1.37e7 | 3.44 | cm^2/s in table |
| DIM | 6.20 | 4.15 | 3.72e10 um^2/min |
| DVX | 5.47e7 | 3.44 | cm^2/s in table |
| DVM | 24.83 | 4.15 | 1.49e11 um^2/min |

### Point defects (Table A-17)

| Parameter | Interstitial | Vacancy |
|---|---|---|
| CEQUIL.0 | 1.25e29 | 1.25e29 |
| CEQUIL.E | 3.26 | 3.26 |
| D.0 | 3.65e-4 | 3.65e-4 |
| D.E | 1.58 | 1.58 |
| NEU.0 | 1.0 | 1.0 |
| NEG.0 / NEG.E | 5.68 / 0.50 | 5.68 / 0.145 |
| POS.0 / POS.E | 5.68 / 0.26 | 5.68 / 0.455 |
| DNEG.0 / DNEG.E | -- | 32.47 / 0.62 |
| KB.0 | 1.0e-21 | -- |
| KB.E | -1.0 | -- |

Note: DC.0 = 1.0, DC.E = 0 for both (defect D is Fermi-independent with defaults).
D.FACTOR = 1, D.11 = D.22 = 1 (isotropic, no scaling).

### IV bulk recombination (Eq 3-286)

With KB.HIGH default, F_IV = EqI * EqV under intrinsic conditions (Eq 3-289/291).
R_b = K_b * (I*V - EqI*EqV) where K_b is diffusion-limited:

```
kIV = 4 * pi * (D_I + D_V) * a0
```

where a0 = lattice parameter = 2.714e-8 cm.

### Surface recombination (Tables A-18/A-19, oxide/silicon interface)

| Parameter | Interstitial | Vacancy |
|---|---|---|
| KSURF.0 | 1.4e-6 | 4.0e-11 |
| KSURF.E | -1.75 | -1.75 |

Negative Ea -> rate increases with temperature.
BC: F_surf = Ksurf * (defect - Eq_defect) at oxide/silicon interface.

### {311} clustering (Table A-24)

See table above in the {311} section.

### Solubility (FLOOXS_2026)

TSUPREM-4 Table A-14 uses tabulated values. We use Arrhenius fits from FLOOXS_2026:
- B: Arr(7.68e22, 0.7086)
- P: Arr(3.89e21, 0.265)
- As: Arr(2.24e22, 0.494)

These are close to the Table A-14 values at our temperatures (900-1100C).

## Variables per dopant

| | B | P | As |
|---|---|---|---|
| dopant | yes | yes | yes |
| Inter | yes | yes | yes |
| Vac | yes | no | yes |
| CIc | yes | yes | yes |
| Total PDEs | 4 | 3 | 4 |

Compared to 5pd: 6 (B), 4 (P), 6 (As).

## Approximations

1. **dM/dt ~ 0** in the pair coupling term. Pair fraction alpha_m ~ 0.3% (900C) to
   ~11% (1000C) from Table A-6. Error is proportional to pair fraction times dC/dt,
   which is small compared to the dominant I sink rates ({311}, surface recombination).

2. **CL.CF = 0.9398 ~ 1.0** in the {311} clustering term. FLOOXS Alagator lacks
   pow(). Same approximation as the 5pd template. (CIc + I)^0.9398 ~ (CIc + I).

3. **Equilibrium activation** instead of transient activation (ACT.TRAN, Eq 3-119).
   At 1000C, tau_act ~ 33s (fast vs 30 min anneal). At 900C, tau_act ~ 14 min
   (significant fraction of anneal), but ACT.MIN = 1.0*ni means boron activates
   immediately above ni. For boron, C0.INI.F = 0 (Table A-16), so all boron starts
   active regardless. Negligible effect at our doses and temperatures.

## Non-deviations (verified zero with default parameters)

These pd.5str terms are zero with default TSUPREM-4 parameters:

- **DDC/BIC** (Table A-16): All DDC parameters = 0. No dopant-defect clustering.
- **K_ni, K_mv** (Eq 3-65/66): Default = 0. No dopant-assisted IV recombination.
- **R_t** (Table A-21): TRAP.CON = 0. No interstitial traps.
- **R_slcI/R_slcV** (Table A-25): ECLUST.0 = 0. No small cluster kinetics.
- **R_llcV** (Table A-24): CL.KRV, CL.KFV not listed = 0. No V-cluster interaction.
- **R_l** (Table A-23): Dislocation loops off by default.
- **Idds**: Follows from DDC = 0.
- **Vc**: Follows from K_rv = K_fv = 0.

## Implementation steps

1. Create `lookupTableGeneration/templates/ion_implant_eff.tcl`
   - Base structure from 5pd (mesh, background, implant, solver, ni, PD params)
   - Replace pair PDEs with effective D + pair coupling in I/V equations
   - 3-4 solution variables instead of 4-6
2. Generate 9 test scripts via gen_test.sh
3. Run all 9 cases (3 dopants x 3 temps)
4. Compare with compare_results.py, update comparison.md

### Template structure

1. Header / math settings
2. Mesh (1D, 5um depth, fine near surface)
3. Background doping (1.4e15)
4. Ion implantation
5. Solver settings
6. ni (Law et al.)
7. Point defect parameters (Table A-17: D, C*, charges)
8. {311} clustering parameters (Table A-24: Kfc, Kr)
9. kIV bulk recombination
10. Per-dopant effective D (Table A-3, pre-computed at anneal temp)
11. Per-dopant solubility (FLOOXS_2026)
12. Surface recombination parameters (Tables A-18/A-19)
13. Solution definitions: dopant, Inter, Vac (B/As only), CIc
14. Fermi level: Charge, Noni, Poni
15. Initialization: Inter/CIc from f_pl, Vac from V_eq(800C)
16. Equilibrium defects: EqInter, EqVac, ScaleInter, ScaleVac
17. Solubility-limited active dopant
18. Charge neutrality
19. Per-dopant DiffDopI, DiffDopV terms
20. Dopant equation
21. Inter equation (with pair coupling)
22. Vac equation (with pair coupling, B/As only)
23. CIc equation
24. Surface recombination BCs
25. Pre-anneal output
26. Anneal sequence (ramp up / dwell / ramp down)
27. Post-anneal output (dopant, Inter, Vac, CIc)
28. Junction depth
