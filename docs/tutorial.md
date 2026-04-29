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

Design optimization (resolution / noise PSD goal functions with parameter and nonlinear constraints) is implemented in the MATLAB version but has not yet been ported to Python. Track [GitHub Issues](https://github.com/MicrosystemsLab/PiezoD/issues) for progress, or use the MATLAB implementation when optimization is required.

## Next Steps

- Browse [`python/examples/`](../python/examples/) for runnable end-to-end scripts.
- Read the docstrings on `Cantilever`, `CantileverEpitaxy`, etc. -- every public method is documented inline.
- For low-level customisation (custom doping profiles, manual sheet-resistance calculations) see `PiezoresistorFromProfile`.
