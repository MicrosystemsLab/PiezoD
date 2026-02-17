# TED Template: Pair Diffusion Model for All Dopants

## Context

The fermi template uses equilibrium defects (no TED) and systematically under-predicts junction depths vs TSUPREM-4's `pd.5str` reference (-21% to -54%). The react template implements full pair diffusion with BIC clusters but only handles boron. We need a single template for B, P, and As with explicit defect transport and dopant-defect pair kinetics, comparable to TSUPREM-4's 5-stream pair diffusion model.

## Approach

New template `ion_implant_ted.tcl` implementing DopantReact-style pair diffusion for all dopants. Follows the same pattern as FLOOXS's built-in `DopantDefectReact` (Dopant.tcl:213-260) but written explicitly for control and debuggability.

## Per-Dopant Physics

| | Boron | Phosphorus | Arsenic |
|---|---|---|---|
| Charge | Acceptor (Poni) | Donor (Noni) | Donor (Noni) |
| Defects | Int + Vac | Int only | Int + Vac |
| Pairs | boronInt, boronVac | phosphorusInt | arsenicInt, arsenicVac |
| Clusters | BIClust (2.5 eV) | none | none |
| Solved vars | 6 | 4 | 5 |
| Solubility | Arr(7.68e22, 0.7086) | Arr(3.89e21, 0.265) | Arr(2.24e22, 0.494) |

### Pair parameters (from FLOOXS_2026/Params/Silicon/)

**Boron-I**: Bind=Arr(8e-23,-1.0), D0=Arr(0.743,3.56), Dp=Arr(0.617,3.56)
**Boron-V**: Bind=Arr(8e-23,-0.5), D0=Arr(0.186,3.56), Dp=Arr(0.154,3.56)
**P-I**: Bind=Arr(8e-23,-1.49), D0=Arr(5.6,3.71), Dn=Arr(6.38,4.05), Dnn=Arr(2.45e-2,3.23)
**As-I**: Bind=Arr(8e-23,0.0), D0=Arr(0.0666,3.45)  [Dn=0, weak binding, quasi-eq in ~0.5ns]
**As-V**: Bind=Arr(8e-23,-0.5), Dn=Arr(12.8,4.05)  [D0=0, electron-enhanced, dominant path]

### Key physics per dopant

- **Boron**: Both I and V paths. BIClust stores excess I with 2.5 eV binding for sustained TED. Acceptor (Poni coupling).
- **Phosphorus**: I-only (Fi=1.0). Strong P-I binding (-1.49 eV) means pairs naturally trap excess I (no clusters needed). K_PI*P ~ 10^4 at 800C, so pre-partition into P-I pairs at init.
- **Arsenic**: Both I and V. Dominant diffusion via As-V (Dn=12.8, electron-enhanced). As-I binding = 0 eV (always quasi-equilibrium, negligible at moderate ScaleI). At ScaleI >> 1, the small As-I path (D0=0.0666) becomes significant. No pre-partitioning needed (pairs negligible at init).

## Critical Design Decision: V Initialization

Previous attempts crashed for P/As because V = I_implant created extreme IV stiffness (kIV*I*V ~ 10^40). The physically correct initialization for the +1 model is V = V_eq (implant creates I, not V):

- V_eq(800C) = Cstar_V(800C) ~ 9e7 cm^-3 (at intrinsic, charge states cancel)
- With V_eq: kIV*I*V ~ 10^17 (manageable vs 10^40)
- Newton converges in 2-3 iterations: V decays to ~0 in implant region within first timestep
- V maintains quasi-steady-state thereafter (replenished by diffusion from bulk)

Why V must be solved explicitly (not algebraic): assuming V = V_eq over-estimates IV recombination (tau ~ 0.7s at 1000C), killing TED prematurely. In reality, V gets depleted in the implant region, slowing I loss and prolonging TED.

## Template Structure

File: `lookupTableGeneration/templates/ion_implant_ted.tcl`
Based on structure of fermi template (switch statements) + defect physics from react template.

```
1. Mesh (same as fermi: 25nm oxide, 5um depth)
2. Background doping (switch: B->n-type P bg, P/As->p-type B bg)
3. Implant (${dopant})
4. Solver settings (same as fermi)
5. ni formula (Law et al.)
6. Point defect params (Cstar, D0, charge states for I and V)
7. Lattice constant, kIV
8. Per-dopant pair params (switch: binding, diffusivities)
9. Per-dopant solubility (switch)
10. BIClust binding (boron only)
11. Surface recombination params
12. Solution definitions (switch: different var sets per dopant)
13. Fermi level: Charge, Noni, Poni
14. Initialization (I from implant, V from V_eq, per-dopant pair partitioning)
15. Equilibrium defect terms (EqInter, EqVac, ScaleInter, ScaleVac)
16. Per-dopant pair terms (switch: Active, Sub, React, DiffPair)
17. Charge neutrality (switch: using Sub not total dopant)
18. Diffusion equations (Inter, Vac, dopant, pairs, clusters)
19. Surface recombination BCs
20. Pre-anneal output
21. Anneal sequence (ramp up, dwell, ramp down)
22. Post-anneal output (dopant, I, V, pairs, clusters)
23. Junction depth (total dopant = free + pairs + clusters)
```

## Equations (DopantReact pattern)

### Pair reaction and diffusivity

For each dopant-defect pair (e.g., boronInt = B-I pair):
```
React = Krate * (DopSub * Defect - Pair / Binding)
DiffPair = (D0 + Dn*Noni + Dp*Poni ...) / (Binding * EqDefect * chg)
```
where chg = Poni (acceptors) or Noni (donors), and Krate = 4*pi*D_defect*lattice.

### Per-dopant equations

**Boron** (6 vars: boron, boronInt, boronVac, BIClust, Inter, Vac):
```
ddt(boron) + ReactBI + ReactBV + ReactClust
ddt(boronInt) - DiffBI*grad(boronInt*Poni)/Poni - ReactBI
ddt(boronVac) - DiffBV*grad(boronVac*Poni)/Poni - ReactBV
ddt(BIClust) - ReactClust
ddt(Inter) - DI*EqI*grad(ScaleI) + kIV*(I*V-EqI*EqV) + ReactBI + ReactClust
ddt(Vac) - DV*EqV*grad(ScaleV) + kIV*(I*V-EqI*EqV) + ReactBV
```

**Phosphorus** (4 vars: phosphorus, phosphorusInt, Inter, Vac):
```
ddt(phosphorus) + ReactPI
ddt(phosphorusInt) - DiffPI*grad(phosphorusInt*Noni)/Noni - ReactPI
ddt(Inter) - DI*EqI*grad(ScaleI) + kIV*(I*V-EqI*EqV) + ReactPI
ddt(Vac) - DV*EqV*grad(ScaleV) + kIV*(I*V-EqI*EqV)
```
Note: Vac has no P coupling (P is I-only). V equation provides IV recombination sink for I.

**Arsenic** (5 vars: arsenic, arsenicInt, arsenicVac, Inter, Vac):
```
ddt(arsenic) + ReactAsI + ReactAsV
ddt(arsenicInt) - DiffAsI*grad(arsenicInt*Noni)/Noni - ReactAsI
ddt(arsenicVac) - DiffAsV*grad(arsenicVac*Noni)/Noni - ReactAsV
ddt(Inter) - DI*EqI*grad(ScaleI) + kIV*(I*V-EqI*EqV) + ReactAsI
ddt(Vac) - DV*EqV*grad(ScaleV) + kIV*(I*V-EqI*EqV) + ReactAsV
```

### Surface recombination (shared)
```
Oxide_Silicon Inter: -KsurfI*(Inter(Silicon) - EqI_surf)
Oxide_Silicon Vac:   -KsurfV*(Vac(Silicon) - EqV_surf)
```
KsurfI = pi*DI*lattice*1.3e15, KsurfV = pi*DV*lattice*1e5

## Initialization

```tcl
# I from implant (+1 model)
sel z=${dopant} name=Inter

# V at thermal equilibrium (800C)
set kT_init [expr {8.617e-5 * 1073.15}]
set V_eq_init [expr {4.0515e26 * exp(-3.97 / $kT_init)}]
sel z=$V_eq_init name=Vac

# Per-dopant pair partitioning:
switch ${dopant} {
    boron {
        # BIClust equilibrium at 800C (strongest binding dominates)
        # Quadratic: x = [(2KB+1) - sqrt(4KB+1)] / (2K)
        set K_cl_init [expr {2.0e-23 * exp(2.5 / $kT_init)}]
        sel z=(2.0*${K_cl_init}*boron+1.0-sqrt(4.0*${K_cl_init}*boron+1.0))/(2.0*${K_cl_init}) name=BIClust
        sel z=boron-BIClust name=Inter
        sel z=Inter name=boron
        sel z=1.0 name=boronInt
        sel z=1.0 name=boronVac
    }
    phosphorus {
        # P-I pair equilibrium at 800C (strong binding, K*P ~ 10^4)
        # Same quadratic form as BIClust
        set K_pi_init [expr {8e-23 * exp(1.49 / $kT_init)}]
        sel z=(2.0*${K_pi_init}*phosphorus+1.0-sqrt(4.0*${K_pi_init}*phosphorus+1.0))/(2.0*${K_pi_init}) name=phosphorusInt
        sel z=phosphorus-phosphorusInt name=Inter
        sel z=Inter name=phosphorus
    }
    arsenic {
        # No pre-partitioning needed (weak binding, V_eq tiny)
        sel z=1.0 name=arsenicInt
        sel z=1.0 name=arsenicVac
    }
}
```

## Junction depth calculation

Total dopant includes free + all pairs + clusters:
```tcl
switch ${dopant} {
    boron    { sel z=boron+boronInt+boronVac+BIClust-1.4e15 }
    phosphorus { sel z=phosphorus+phosphorusInt-1.4e15 }
    arsenic  { sel z=arsenic+arsenicInt+arsenicVac-1.4e15 }
}
```

## Other file changes

### gen_test.sh: add template argument

```bash
#!/bin/bash
# Usage: bash gen_test.sh <dopant> <dose> <energy> <temp> <time> [template]
TEMPLATE="${6:-../templates/ion_implant_fermi.tcl}"
sed -e "s/\${dopant}/$1/g" ... "$TEMPLATE"
```

Default remains fermi template. Use `../templates/ion_implant_ted.tcl` for TED model.

### README.md: add TED model description

Add entry to Templates table with physics description and validation results.

## Solver Risk: P/As with Vac

P has no V-coupling through pairs (I-only diffusion). If Newton fails to converge for P because V is weakly coupled, fallback options:
1. Drop V for phosphorus (3 vars, use surface recomb as sole I sink, P-I pairing provides natural I trapping)
2. Add tiny artificial P-V coupling for Jacobian conditioning

## Verification

Run same 6 test cases as fermi model (B and P at 2 dose/energy combos, As at 2):

```bash
cd lookupTableGeneration
# Generate from TED template
bash scripts/gen_test.sh boron 2e14 80 1000 30 ../templates/ion_implant_ted.tcl > simulations/test_ted_B_2e14_80.tcl
bash scripts/gen_test.sh boron 2e15 20 1000 30 ../templates/ion_implant_ted.tcl > simulations/test_ted_B_2e15_20.tcl
bash scripts/gen_test.sh phosphorus 2e14 80 1000 30 ../templates/ion_implant_ted.tcl > simulations/test_ted_P_2e14_80.tcl
bash scripts/gen_test.sh phosphorus 2e15 50 1000 30 ../templates/ion_implant_ted.tcl > simulations/test_ted_P_2e15_50.tcl
bash scripts/gen_test.sh arsenic 2e14 80 1000 30 ../templates/ion_implant_ted.tcl > simulations/test_ted_As_2e14_80.tcl
bash scripts/gen_test.sh arsenic 2e15 50 1000 30 ../templates/ion_implant_ted.tcl > simulations/test_ted_As_2e15_50.tcl

# Run in Docker
MSYS_NO_PATHCONV=1 docker compose run --rm flooxs test_ted_B_2e14_80.tcl
```

Success criteria:
1. All 6 cases converge (no exit 139, no inf in Jacobian)
2. Boron Xj closer to TSUPREM-4 reference than fermi model (target: within 20%)
3. Phosphorus Xj closer to reference
4. Arsenic Xj closer to reference
5. If P or As fails with Vac, apply fallback and re-test

TSUPREM-4 reference values (1000C, 30min):
| Dopant | Dose/Energy | TSUPREM-4 Xj (um) | Fermi Xj (um) |
|--------|-------------|-------------------|----------------|
| B | 2e14/80keV | 1.230 | 0.652 |
| B | 2e15/20keV | 0.880 | 0.639 |
| P | 2e14/80keV | 0.710 | 0.547 |
| P | 2e15/50keV | 0.800 | 0.832 |
| As | 2e14/80keV | 0.520 | 0.237 |
| As | 2e15/50keV | 0.400 | 0.317 |
