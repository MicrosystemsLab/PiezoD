# PiezoD Tutorial

This tutorial walks through using PiezoD to calculate the properties of a single-crystal silicon piezoresistive cantilever. It covers an epitaxial design; diffusion, ion-implantation, polycrystalline, and piezoelectric variants follow the same pattern with a different `Cantilever*` class.

## Install

```
pip install piezod
```

PiezoD requires Python 3.12 or later. The `examples/` directory in the repository has runnable scripts that mirror this tutorial.

## Import

```python
from piezod import CantileverEpitaxy
```

## Define the Cantilever

All dimensions are in MKS units (meters, Hz, V, cm^-3 for concentrations).

```python
c = CantileverEpitaxy(
    freq_min=1,                      # measurement bandwidth (Hz)
    freq_max=1000,
    l=100e-6,                        # length: 100 um
    w=4e-6,                          # width: 4 um
    t=1e-6,                          # thickness: 1 um
    l_pr_ratio=0.8,                  # piezoresistor extends 80% of length
    t_pr_ratio=0.5,                  # piezoresistor is 50% of thickness
    v_bridge=2.0,                    # 2V Wheatstone bridge bias
    doping_type="phosphorus",        # n-type, <100> orientation
    dopant_concentration=1e18,       # cm^-3
)
```

Choosing `"phosphorus"` (n-type) implicitly uses the <100> orientation where the piezoresistive coefficient peaks. Choosing `"boron"` (p-type) uses the <110> direction (E = 169 GPa) with the corresponding coefficients.

Attributes can also be set after construction:

```python
c.fluid = "vacuum"
c.number_of_piezoresistors = 2
```

## Calculate Properties

### Mechanical

```python
c.stiffness()              # spring constant (N/m)
c.omega_vacuum_hz()        # vacuum resonant frequency (Hz)
f_d, Q = c.omega_damped_hz_and_Q()  # damped frequency (Hz) and Q
```

### Electrical

```python
c.sheet_resistance()       # ohms/square
c.resistance()             # ohms
c.number_of_carriers()     # carriers in the piezoresistor
```

### Noise

```python
c.integrated_noise()       # total integrated voltage noise (V)
c.johnson_integrated()     # Johnson-Nyquist contribution (V)
c.hooge_integrated()       # 1/f (Hooge) contribution (V)
c.amplifier_integrated()   # amplifier contribution (V)
c.plot_noise_spectrum()    # voltage-noise PSD vs frequency
```

### Sensitivity and Resolution

```python
c.force_sensitivity()        # V/N
c.beta()                     # piezoresistor efficiency factor
c.force_resolution()         # minimum detectable force (N)
c.displacement_resolution()  # minimum detectable displacement (m)
```

### Summary

```python
c.print_performance()
```

prints a formatted block covering geometry, frequencies, sensitivity, resolution, resistances, noise breakdown, and thermal estimates.

## Change the Operating Environment

```python
c.fluid = "water"   # "vacuum", "air", or "water"
```

This changes the hydrodynamic damping and lowers the damped resonant frequency.

## Inspect the Doping Profile

Every cantilever exposes a `doping_profile()` method that returns the depth coordinate, net active carrier concentration, and total dopant concentration:

```python
import matplotlib.pyplot as plt

z, active, total = c.doping_profile()

fig, ax = plt.subplots()
ax.semilogy(z * 1e9, total, label="Total")
ax.semilogy(z * 1e9, active, "--", label="Net active")
ax.set_xlabel("Depth (nm)")
ax.set_ylabel("Concentration (cm$^{-3}$)")
ax.legend()
plt.show()
```

## Other Cantilever Types

The same pattern applies to the other implementations:

| Class | Doping / sensing approach |
|---|---|
| `CantileverEpitaxy` | Step-function epitaxial doping |
| `CantileverDiffusion` | POCl3 / boron diffusion (kink profile) |
| `CantileverImplantation` | Ion-implant profile from lookup tables |
| `CantileverPoly` | Polycrystalline thin-film (poly-Si, Ti, Al) |
| `CantileverPiezoelectric` | AlN or PZT piezoelectric sensing |

See [`python/examples/`](../python/examples/) for a runnable script for each.

## Optimization

Design optimization is built around three pieces:

- A **goal callable** that maps a cantilever to the scalar to minimize. PiezoD ships factories that match the MATLAB tutorial's units (pN, nm, pN/sqrt(Hz), etc.):

```python
from piezod import (
    force_resolution_goal,
    displacement_resolution_goal,
    force_noise_density_goal,
    surface_stress_resolution_goal,
)
```

- Optional **parameter_constraints** that override the default state-variable bounds (e.g. clamp thickness, cap bias voltage). Keys are `min_<name>` / `max_<name>` for any state variable returned by `c.optimization_state_vars()`.

- Optional **metric_constraints** -- inequality constraints on derived quantities like power dissipation, resonant frequency, or stiffness. Each constraint is a `CantileverMetricConstraint` referencing a `CantileverMetric` enum value.

```python
from piezod import (
    CantileverMetric,
    CantileverMetricConstraint,
    optimize_performance,
)

constraints = [
    CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
    CantileverMetricConstraint(CantileverMetric.OMEGA_VACUUM_HZ, minimum=5 * c.freq_max),
    CantileverMetricConstraint(CantileverMetric.STIFFNESS, minimum=1e-3, maximum=1e1),
]

result = optimize_performance(
    c,
    force_resolution_goal(),
    parameter_constraints={"max_v_bridge": 10.0},
    metric_constraints=constraints,
    n_starts=5,
    max_iterations=10,
    random_seed=0,
)
```

`optimize_performance` runs SciPy's SLSQP (or L-BFGS-B when there are no nonlinear constraints) from random initial conditions and returns the best result once two starts agree within `convergence_tolerance` (1% by default), capped at `max_iterations` total runs. Use `optimize_performance_from_current` for a single-shot refinement of the existing design.

```python
print(f"force_resolution: {c.force_resolution() * 1e12:.1f} pN -> "
      f"{result.optimized.force_resolution() * 1e12:.2f} pN")
```

Default geometric sanity constraints (`l/w >= 2`, `w/t >= 2`, `l_pr/w_pr >= 2`, `l_pr >= 2 um`) are added automatically; pass `default_aspect_constraints=False` to opt out.

For a runnable end-to-end example with full output, see [`python/examples/optimization.py`](../python/examples/optimization.py). The example mirrors the MATLAB tutorial flow and typically improves Harley-1999-style force resolution by several hundred times.

For implant-process-only optimization (geometry fixed, only `annealing_time/temp` and `implantation_energy/dose` varied), `CantileverImplantation.optimize_doping_for_hooge_noise` is also available with a Hooge-noise-limited default objective.

## Next Steps

- Browse [`python/examples/`](../python/examples/) for runnable end-to-end scripts.
- Read the docstrings on `Cantilever`, `CantileverEpitaxy`, etc. -- every public method is documented inline.
- For low-level customisation (custom doping profiles, manual sheet-resistance calculations) see `PiezoresistorFromProfile`.
