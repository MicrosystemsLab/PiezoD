# Diffusion Theory for Ion-Implanted Dopants in Silicon

This document describes the physics models used in the FLOOXS process simulator for dopant diffusion, and how they relate to the TSUPREM-4 reference implementation. The goal is to reproduce TSUPREM-4's Fermi-level dependent diffusion results using FLOOXS.

## Physical Picture

After ion implantation, dopant atoms sit on substitutional silicon lattice sites. They diffuse by forming mobile pairs with point defects (interstitials I and vacancies V). The diffusion rate depends on:

1. How many point defects are available (defect equilibrium)
2. How fast the dopant-defect pairs move (pair diffusivity)
3. The local Fermi level, which affects both (1) and (2)

At low doping (below the intrinsic carrier concentration ni), the Fermi level is near midgap and diffusion is "intrinsic" -- a simple constant at each temperature. At high doping (above ni), the Fermi level shifts and diffusion accelerates. For acceptors like boron, this acceleration can be 10x or more at typical implant doses.

## Key Variables

| Symbol | Definition | Units |
|--------|-----------|-------|
| ni | Intrinsic carrier concentration at temperature T | cm^-3 |
| n, p | Electron and hole concentrations | cm^-3 |
| Noni | n/ni (normalized electron concentration) | dimensionless |
| Poni | p/ni (normalized hole concentration) | dimensionless |
| Cstar | Equilibrium point defect concentration (intrinsic) | cm^-3 |
| Eq_I, Eq_V | Equilibrium defect concentration (Fermi-corrected) | cm^-3 |

At intrinsic doping: Noni = Poni = 1, n = p = ni.

For boron (acceptor) doping above ni: p >> ni, so Poni >> 1 and Noni << 1.

## Charge Neutrality

The local carrier concentrations are determined by the net ionized dopant concentration. For charge neutrality (no space charge):

```
p - n + N_D - N_A = 0
```

where N_D is the ionized donor concentration and N_A is the ionized acceptor concentration. Combined with np = ni^2:

```
p = 0.5 * (N_A - N_D + sqrt((N_A - N_D)^2 + 4*ni^2))
n = ni^2 / p
```

For boron-only implant into n-type background (phosphorus at ~1.4e15):
- Where boron >> ni: p ~ boron, Poni ~ boron/ni
- Where boron << ni: p ~ ni, Poni ~ 1

This is the approach used by TSUPREM-4's "Fermi" model and FLOOXS's `PotentialEqns` proc in non-Poisson mode.

## Dopant Diffusion Equation

### Intrinsic (Noni = Poni = 1)

Simple Fick's law with defect enhancement:

```
dC/dt = div(D_BI * (I/I*) * grad(C)) + div(D_BV * (V/V*) * grad(C))
```

where D_BI, D_BV are the boron-interstitial and boron-vacancy pair diffusivities, and I/I* is the interstitial supersaturation ratio.

### Fermi-level dependent (FLOOXS DopantFermi model)

The diffusivity becomes Fermi-level dependent:

```
D = D0 + Dp * Poni    (for acceptors like boron)
D = D0 + Dn * Noni    (for donors like phosphorus)
```

Boron parameters (from FLOOXS_2026/Params/Silicon/Boron/):

| Path | D0 | Dp |
|------|----|----|
| Boron-Interstitial | Arrhenius(0.743, 3.56) | Arrhenius(0.617, 3.56) |
| Boron-Vacancy | Arrhenius(0.186, 3.56) | Arrhenius(0.154, 3.56) |

The transport equation for an acceptor also includes the built-in field drift:

```
dC/dt = div( (D/Poni) * grad(C_active * Poni) )
```

The grad(C * Poni) term captures both concentration gradient diffusion and electric field drift. The 1/Poni prefactor normalizes the effective diffusivity.

### Example: Boron 2e14, 80 keV, 1000C

At 1000C:
- ni ~ 7e18 cm^-3
- Peak boron ~ 1e20 cm^-3
- Poni at peak ~ 1e20 / 7e18 ~ 14
- D_total = D0 + Dp*14 ~ 0.743 + 0.617*14 ~ 9.4 (in Arrhenius prefactor units)
- D_intrinsic = D0 + Dp*1 = 1.36

The diffusivity at the peak is ~7x higher than intrinsic. This enhancement falls off away from the peak as boron concentration drops below ni.

## Point Defect Equilibrium

### Multi-charge-state model

Point defects exist in multiple charge states. The equilibrium concentration depends on which states are populated:

```
Eq_I = Cstar_I * (neutral + Noni*negative + Poni*positive) / (neutral + negative + positive)
Eq_V = Cstar_V * (neutral + Noni*(negative + Noni*dnegative) + Poni*positive) / (neutral + negative + dnegative + positive)
```

Charge-state parameters (from FLOOXS_2026/Params/Silicon/):

| Defect | neutral | negative | positive | dnegative |
|--------|---------|----------|----------|-----------|
| Interstitial | 1.0 | Arrhenius(5.68, 0.48) | Arrhenius(5.68, 0.42) | -- |
| Vacancy | 1.0 | Arrhenius(5.68, 0.145) | Arrhenius(5.68, 0.455) | Arrhenius(32.47, 0.62) |

At intrinsic (Noni = Poni = 1), the numerator equals the denominator and Eq = Cstar. In p-type silicon (Poni >> 1), the positive charge state contribution increases Eq_I.

### Detailed balance transport

The defect transport equation uses the thermodynamically correct detailed balance form:

```
dI/dt = div(D_I * Eq_I * grad(I / Eq_I)) - kIV * (I*V - Eq_I*Eq_V)
```

This reduces to Fick's law when Eq is spatially uniform, but includes drift in equilibrium gradients when the Fermi level varies with position.

## Self-Consistent Coupling

The full model is a coupled system:

```
1. Charge neutrality:  Potential = f(Charge)
2. Fermi factors:      Noni = exp(Potential/Vt), Poni = exp(-Potential/Vt)
3. Defect equilibrium: Eq_I = g(Noni, Poni), Eq_V = h(Noni, Poni)
4. Defect transport:   dI/dt = div(D_I * Eq_I * grad(I/Eq_I)) - recomb
5. Dopant transport:   dC/dt = div(D(Noni,Poni) * grad(C * Poni) / Poni)
6. Charge update:      Charge = -C_active  (for acceptors)
```

Steps 1-6 are solved simultaneously by the FLOOXS nonlinear solver at each timestep. The solver iterates until all equations are self-consistent.

## FLOOXS Built-in Model Procs

FLOOXS implements the above physics in three Tcl procs that are designed to be called in sequence:

### PotentialEqns (Potential.tcl)

Creates `Noni` and `Poni` as computed fields derived from a `Potential` variable. In non-Poisson mode (charge neutrality), solves the algebraic equation:

```
Potential / Vt = log(0.5 * (Charge + sqrt(Charge^2 + 4*ni^2)) / ni)
```

This is equivalent to TSUPREM-4's Fermi model.

### DefectBulk (Defect.tcl)

Creates `Eq_I`, `Eq_V`, `Scale_I = I/Eq_I`, `Scale_V = V/Eq_V` as computed fields. Sets up the detailed balance transport equation with I-V recombination. References `Noni` and `Poni` (must already exist from PotentialEqns).

### DopantBulk (Dopant.tcl)

Dispatches to one of four models based on the `DiffModel` parameter:
- 0: Constant (intrinsic diffusivity)
- 1: Fermi (D0 + Dp*Poni for acceptors) -- default for boron
- 2: Pair (explicit dopant-defect pair concentrations)
- 3: React (non-equilibrium kinetic model)

Adds the dopant's charge contribution to the `Charge` term, closing the self-consistent loop.

### Parameter database (pdb)

All parameters are stored in FLOOXS's material database, loaded from files in `Params/Silicon/{Species}/`. The procs use `pdbDelayDouble` for temperature-dependent lookups and `pdbGetSwitch` for model selection.

## TSUPREM-4 Reference Model

The TSUPREM-4 simulation template uses:

```
method vertical ddc.full pd.5str
```

`pd.5str` is the 5-stream point defect model with charge-state dependent equilibria. The default diffusion model for boron is "Fermi", which computes:

```
D_eff = D_intrinsic + D_acceptor * (p/ni)
```

This is equivalent to FLOOXS's DopantFermi model (DiffModel=1).

## Impact on Junction Depth

Without Fermi-level coupling (Noni = Poni = 1), the boron diffusivity is purely intrinsic and the junction depth for 2e14/80keV/1000C/30min is 0.695 um. TSUPREM-4 with Fermi coupling gives 1.23 um -- a 77% deeper junction. The missing Fermi enhancement is the dominant source of discrepancy.

## References

- M.D. Deal, "Semiconductor Process Modeling" (Stanford EE410 notes)
- FLOOXS documentation: http://www.flooxs.ece.ufl.edu/
- FLOOXS_2026/TclLib/Models/ source code
- FLOOXS_2026/Params/Silicon/ parameter files
