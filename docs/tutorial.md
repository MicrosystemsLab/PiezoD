# PiezoD Tutorial

This tutorial walks through using PiezoD to calculate the properties of a single crystal silicon cantilever beam.

## Setup

Start MATLAB and add the PiezoD code to your path:

```matlab
addpath('PiezoD');
```

Clear any existing variables (run twice to avoid occasional optimization slowdowns):

```matlab
clear
clear
```

## Define Cantilever Geometry

Define the length, width and thickness of the cantilever beam. All dimensions are in MKS units.

```matlab
l = 100e-6;  % 100 um
w = 4e-6;    % 4 um
t = 1e-6;    % 1 um
```

## Define Piezoresistor Properties

```matlab
l_pr_ratio = 0.8;           % PR extends 80% of cantilever length
t_pr_ratio = 0.5;           % PR is 50% of beam thickness (500 nm)

doping_type = 'phosphorus'; % n-type, <100> orientation, E = 130 GPa
dopant_concentration = 1e18;

v_bridge = 2;               % 2V across Wheatstone bridge
```

Choosing `'phosphorus'` implicitly sets the cantilever to the <100> orientation where the piezoresistive coefficient is maximum. Choosing `'boron'` would use the <110> direction (E = 169 GPa) with corresponding piezoresistive coefficients.

## Define Measurement Bandwidth

```matlab
freq_min = 1;    % Hz
freq_max = 1000; % Hz
```

## Create the Cantilever Object

```matlab
c_epitaxy = cantileverEpitaxy(freq_min, freq_max, ...
    l, w, t, l_pr_ratio, v_bridge, doping_type, dopant_concentration, ...
    t_pr_ratio);
```

Examine the object:

```matlab
c_epitaxy
```

List available methods:

```matlab
methods('cantileverEpitaxy')
```

## Calculate Properties

### Mechanical Properties

```matlab
c_epitaxy.stiffness()        % Spring constant (N/m)
c_epitaxy.omega_vacuum_hz()  % Resonant frequency in vacuum (Hz)
c_epitaxy.omega_damped_hz()  % Resonant frequency in fluid (Hz)
```

### Electrical Properties

```matlab
c_epitaxy.number_of_carriers()
c_epitaxy.sheet_resistance()
c_epitaxy.resistance()

c_epitaxy.integrated_amplifier_noise()
c_epitaxy.integrated_hooge_noise()
c_epitaxy.integrated_johnson_noise()
c_epitaxy.integrated_noise()
c_epitaxy.johnson_noise_density()

c_epitaxy.plot_noise_spectrum()
```

### Sensitivity and Resolution

```matlab
c_epitaxy.force_sensitivity()       % V/N
c_epitaxy.beta()                    % Piezoresistor efficiency factor
c_epitaxy.force_resolution()        % N
c_epitaxy.displacement_resolution() % m
```

### Print Summary

```matlab
c_epitaxy.print_performance()           % Formatted output
c_epitaxy.print_performance_for_excel() % Tab-delimited for spreadsheets
```

## Change Operating Environment

```matlab
c_epitaxy.fluid = 'water';  % Options: 'vacuum', 'air', 'water'
```

## Optimization

Define constraints:

```matlab
omega_min_hz = 5*freq_max;  % Resonant freq >= 5x measurement bandwidth
max_power = 2e-3;           % Max 2 mW power dissipation

constraints = {{'min_t', 'min_w', 'max_v_bridge'}, {1e-6, 4e-6, 10}};
```

Run optimization:

```matlab
c_epitaxy = c_epitaxy.optimize_performance_from_current(max_power, omega_min_hz, constraints);
```

Check improved performance:

```matlab
c_epitaxy.print_performance()
```

A typical optimization improves force resolution from ~125 pN to ~2.5 pN.

## Next Steps

- See `sampleCode.m` for additional examples including diffusion and ion implantation
- Check `cantilever.m` comments for advanced features like thermal actuator modeling and Monte Carlo stress analysis
