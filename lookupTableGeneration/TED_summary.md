# FLOOXS TED Simulation Summary

## Project Context

We're generating ion implant diffusion lookup tables for the **PiezoD** piezoresistor design tool. The lookup tables map (dopant, dose, energy, anneal temperature, anneal time) to doping profiles. We use **FLOOXS**, an open-source TCAD process simulator (Tcl-based, C++ solver), running in Docker.

Our target is to match **TSUPREM-4** (commercial TCAD) reference results, specifically for **Transient Enhanced Diffusion (TED)** -- the anomalous boost in dopant diffusion caused by excess point defects (interstitials) created during ion implantation.

**Test case**: Boron, 2e15 cm^-2 dose, 20 keV energy, 1000C/30min anneal. TSUPREM-4 reference junction depth: **Xj = 0.880 um**.

## Key Question

FLOOXS has built-in model procs (DopantPair, DopantReact) that implement TED physics including boron-interstitial clusters (BIC). We cannot use them because the C++ expression parser crashes at solve time. The same `[Arrhenius]` expressions work fine in our hand-written template, so the crash is triggered by something specific to the built-in proc equation assembly -- possibly term redefinition, `nosolve` solution variables, or expression tree depth. We need to either fix/work around this crash to unlock the built-in models, or build a physically motivated TED model by hand that can match TSUPREM-4 across a range of implant conditions.

## FLOOXS Architecture

FLOOXS is a 1D/2D process simulator with:
- **Tcl scripting layer**: User writes `.tcl` scripts using FLOOXS commands (`line`, `region`, `init`, `implant`, `diffuse`, `term`, `solution`, `sel`, `pdbSetDouble`, etc.)
- **C++ PDE solver**: Solves coupled diffusion-reaction equations on a 1D mesh
- **PDB (Parameter Database)**: Stores material/dopant parameters, accessed via `pdbSetDouble`/`pdbGetDouble`/`pdbDelayDouble`
- **Built-in model procs** (Tcl): `PotentialEqns`, `DefectBulk`, `DefectBound`, `DopantBulk`, `DopantFermi`, `DopantPair`, etc. in `TclLib/Models/`
- **`pdbDelayDouble`**: Returns a *string* expression (not a number) that gets evaluated at solve time. E.g., `pdbDelayDouble Silicon Int Cstar` might return `"[Arrhenius 3.6484e27 3.7]"` which the C++ solver evaluates as `3.6484e27 * exp(-3.7 / kT)`.

Key FLOOXS commands:
- `term name = X add eqn = "expr" Material` -- defines a named algebraic expression
- `solution name = X solve !negative add` -- declares X as a solved PDE variable
- `solution name = X add const val = "expr" Material` -- declares X as an algebraic (non-solved) variable
- `pdbSetString Material Species Equation "PDE_expr"` -- sets the PDE for a species
- `diffuse time=T temp=C init=dt damp.trbdf` -- runs the time-stepping solver

## What Works: Hand-Written Template

We have a hand-written Tcl template that sets up the coupled PDE system directly.

**Equations solved** (3 coupled PDEs + algebraic terms):
1. **Boron**: `ddt(boron) - DiffDopI * grad(boronActive * ScaleInter * Poni) - DiffDopV * grad(boronActive * ScaleVac * Poni)`
2. **Interstitials (Inter)**: `ddt(Inter) - DiffI*EqInter*grad(ScaleInter) + kIV*(Inter*Vac - EqInter*EqVac)`
3. **Vacancies (Vac)**: `ddt(Vac) - DiffV*EqVac*grad(ScaleVac) + kIV*(Inter*Vac - EqInter*EqVac)`

**Algebraic terms**:
- `EqInter = Cstar_I * (fneu + Noni*fneg + Poni*fpos) / fden` (Fermi-level-dependent equilibrium)
- `EqVac = Cstar_V * (...)` (similar)
- `ScaleInter = Inter/EqInter`, `ScaleVac = Vac/EqVac`
- `Charge = Nbackground - boronActive` (net charge)
- `Noni = n/ni`, `Poni = p/ni` (normalized carrier concentrations from charge neutrality)
- `boronActive = Css * boron / (Css + boron)` (solid solubility limited activation)

**Surface boundary conditions** at Oxide/Silicon interface:
- `KsurfI * (Inter - EqInter)` and `KsurfV * (Vac - EqVac)`

This template converges and produces results. The physics issue is how to initialize the defect concentrations.

## The TED Physics Problem

After ion implantation, excess interstitials are created (roughly 1 per implanted ion, the "+1 model"). These enhance dopant diffusion. The question is: how to initialize and model these excess defects?

### Approaches Tried

#### 1. I=V Cap Approach (works but unphysical)

Initialize both Inter and Vac to `min(boron_profile, N * EqInter_neutral)` where N is a super-saturation ratio.

| Cap N | Xj (um) | vs TSUPREM-4 |
|-------|---------|-------------|
| No TED | 0.560 | -36.4% |
| 100 | 0.517 | -41.2% |
| 1,000 | 0.523 | -40.6% |
| 10,000 | 0.646 | -26.6% |
| **30,000** | **0.847** | **-3.8%** |
| 50,000 | 0.975 | +10.8% |
| 100,000 | 1.183 | +34.4% |
| 1,000,000 | 2.149 | +144.2% |

**N=30,000 matches best** but N is a fitting parameter with no physical basis. Xj grows as ~ln(N) because time-integrated supersaturation scales logarithmically. The I=V initialization is wrong because it creates excess Frenkel pairs.

#### 2. Vac=Equilibrium Approach (physical but doesn't work)

Initialize Inter to capped boron profile, but Vac to neutral equilibrium (~9.1e7 cm^-3 at 800C).

| Cap N | Xj (um) | vs TSUPREM-4 |
|-------|---------|-------------|
| 30,000 | 0.517 | -41.2% |
| 50,000 | 0.515 | -41.5% |
| 100,000 | 0.516 | -41.4% |

**Completely insensitive to N. Actually LESS diffusion than no-TED case.** The I-V recombination term `kIV*(I*V - I*V*)` rapidly depletes vacancies below equilibrium, making ScaleVac << 1. Since boron diffuses via both I and V mechanisms, the vacancy depletion hurts more than interstitial enhancement helps.

#### 3. One-Defect Model (conceptual analysis)

Track only interstitials with a first-order bulk sink: `dI/dt = -k_bulk * (I - I*)`. With physical parameters: k_bulk = kIV * V* ~ 1.6 s^-1 at 1000C, giving tau ~ 0.6 s. Far too fast -- excess interstitials vanish in seconds, no meaningful TED over a 30-min anneal.

#### 4. Simple BIC (conceptual analysis)

BI pair binding energy ~0.5-1.0 eV -> dissolves in nanoseconds at 1000C. Too fast to act as reservoir. The TSUPREM-4 approach uses **higher-order clusters** (B3I, B4I with Eb=2.5-3.5 eV) which provide slow release over minutes.

### 5. Built-in Model Procs (blocked by FLOOXS bug)

FLOOXS has built-in Tcl procs that implement the full coupled defect-dopant system including DopantPair (DiffModel=2) and DopantReact (DiffModel=3, for BIC). We tried to use these.

**What we had to work around:**

1. **`term list` command not supported**: DefectBulk checks `[term list]` for Pressure terms. Fix: Pre-created `PressureInt` and `PressureVac` terms, redefined DefectBulk to skip the check.

2. **`term name=Charge print` not supported**: DopantBulk reads current Charge equation via `term name=Charge print`. Fix: Tracked Charge equation in a Tcl global variable `_charge_eqn`.

3. **Interface equations not initialized**: DefectBound tries `pdbGetString Oxide_Silicon Vac Equation` before it's set. Fix: `pdbSetString Oxide_Silicon Int Equation ""` before calling DefectBound.

4. **DiffLimit braced list at solve time**: `[DiffLimit Silicon {Int Vac} 0.0]` was deferred to solve time where C parser can't handle Tcl braces. Fix: Called DiffLimit at setup time.

5. **Nested Arrhenius from DiffLimit**: DiffLimit returns `[Arrhenius {(4*pi*([Arrhenius 0.138 1.37]+[Arrhenius 1.18e-4 0.1])*a)} 0.0]` -- nested Arrhenius calls. Fix: Computed kIV directly using `pdbDelayDouble` values.

6. **Trapping terms reference `Int_Silicon`**: DefectBoundSide generates terms referencing side-qualified species names. Fix: Simplified to skip trapping (Ktrap=0 in our setup).

**What we couldn't work around -- the FLOOXS C++ crash:**

After all Tcl-level fixes, setup completes perfectly (all equations generated, all terms created, all solutions declared). But the **solver crashes** at the first timestep with:

```
message that should never appear, yet
I know it will probably appear someday.

I thought was an unreachable state.
I could give you advice for what to do.
But since I never believed you would get here,
I don't have any.

called for ListExpr base class
typeinfo for base class Expr!
```

This is a C++ assertion in the expression tree evaluator/differentiator. It crashes on **any** `[Arrhenius]` expression at solve time. We tried:
- Simplifying Arrhenius sums (combining prefactors with same Ea)
- Disabling PotentialEqns (replacing with direct Noni/Poni)
- Setting ActiveModel=0 (skip Solubility Arrhenius)
- Each time it crashes on the next Arrhenius expression it encounters

**Critically**: The exact same Arrhenius expressions work perfectly in our hand-written template. The crash only happens when equations are set up through the built-in model proc pathway. The difference may be in term ordering, the number of auto-created solution variables, or some subtle interaction in how the equation system is assembled.

## Current State of test_builtin.tcl

The file sources FLOOXS built-in procs and calls them with our own overrides for DefectBulk, DopantBulk, DefectBound. It successfully:
- Sets up the mesh, regions, implant
- Creates all PDB parameters
- Calls DefectInit, DefectBulk, DopantBulk, DefectBound
- All equations and terms are generated correctly

But crashes at solve time due to the C++ expression parser bug.

## Key Files

Simulation scripts:
- `lookupTableGeneration/simulations/test_builtin.tcl` -- Built-in model attempt (crashes)
- `lookupTableGeneration/simulations/test_ted_cap30000_veq.tcl` -- Hand-written, Vac=eq, N=30k (wrong TED)
- `lookupTableGeneration/simulations/test_ted_cap30000.tcl` -- I=V cap at N=30k (Xj=0.847, closest to TSUPREM-4)
- `lookupTableGeneration/simulations/test_ted_cap50000_veq.tcl` -- Hand-written, Vac=eq, N=50k
- `lookupTableGeneration/simulations/test_ted_cap100000_veq.tcl` -- Hand-written, Vac=eq, N=100k
- `lookupTableGeneration/simulations/test_convergence_no_ted.tcl` -- Fermi-only, no defects (Xj=0.560)
- `lookupTableGeneration/templates/ion_implant_fermi.tcl` -- Fermi-level dependent diffusion (DopantPair model)
- `lookupTableGeneration/templates/ion_implant_react.tcl` -- Dopant-defect pair reactions (DopantReact/BIC model)
- `lookupTableGeneration/Dockerfile` -- FLOOXS build (from FLOOXS_2026 source)

FLOOXS full source tree (available locally at `lookupTableGeneration/FLOOXS_2026/`):
- `FLOOXS_2026/TclLib/Models/Dopant.tcl` -- DopantBulk, DopantFermi, DopantPair, DopantReact procs
- `FLOOXS_2026/TclLib/Models/Defect.tcl` -- DefectBulk, DefectInit, DefectBound procs
- `FLOOXS_2026/TclLib/Models/Potential.tcl` -- PotentialEqns proc
- `FLOOXS_2026/TclLib/Name.tcl` -- pdbName, FirstMat, SecondMat helpers
- `FLOOXS_2026/Params/paramFunc` -- DiffLimit, ConcBind, SurfDiffLimit procs
- `FLOOXS_2026/Params/Silicon/Boron/` -- Default boron parameters (Interstitial, Vacancy subdirs)
- `FLOOXS_2026/src/` -- C++ source code (math/, mesh/, etc.)
- `FLOOXS_2026/src/web/assets/` -- HTML documentation (main.html, process.html, cli.html, device.html, template.html)
- `FLOOXS_2026/Test/` -- Test scripts and examples

FLOOXS wiki documentation (261 HTML pages, mirrored locally):
- `lookupTableGeneration/www.flooxs.ece.ufl.edu/index.php/` -- Full mirror of the FLOOXS wiki
- Key pages include: `Diffusion.html`, `Diffuse_Command.html`, `Alagator.html`, `Alagator_Language_Description.html`, `Convergence.html`, `Command_Reference_Library.html`, `Example_(1D).html`, `EquationCommand.html`, `Add_Solution_Variables.html`, `Defining_a_Grid.html`, `Callback_Procedures.html`, etc.

## What We Need Help With

1. **Understanding the FLOOXS C++ crash**: Why does `[Arrhenius X Y]` work in hand-written term equations but crash when the same expressions are generated by the built-in model procs? Is there a term ordering issue, a maximum expression complexity issue, or something else?

2. **Alternative TED modeling**: Is there a way to get physically meaningful TED without full BIC clusters? The one-defect model with physical kIV is too fast. The I=V cap works but is a fitting parameter. Is there a well-known effective interstitial lifetime or effective sink term that gives the right timescale?

3. **Making built-in procs work**: Are there simpler workarounds we're missing? Could we pre-evaluate all Arrhenius expressions to numeric constants at setup time (losing temperature dependence during ramps) to avoid the parser crash?
