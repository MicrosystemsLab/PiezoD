# PiezoD MATLAB

The primary implementation of PiezoD, written in object-oriented MATLAB.

## Requirements

- MATLAB R2008a or later
- Optimization Toolbox

## Usage

```matlab
addpath('PiezoD');

% Create an epitaxial cantilever (Harley 1999 configuration)
freq_min = 10;
freq_max = 1000;
l = 300e-6;
w = 44e-6;
t = 89e-9;
l_pr_ratio = 45/300;
v_bridge = 5;
doping_type = 'boron';
concentration = 4e19;
t_pr_ratio = 30/89;

c = cantileverEpitaxy(freq_min, freq_max, l, w, t, l_pr_ratio, ...
    v_bridge, doping_type, concentration, t_pr_ratio);
c.fluid = 'vacuum';
c.print_performance();

% Optimize
parameter_constraints = {{'min_t', 'max_t', 'min_w', 'max_v_bridge'}, ...
    {1e-6, 1e-6, 5e-6, 10}};
nonlinear_constraints = {{'omega_min_hz', 'max_power', 'min_k', 'max_k'}, ...
    {freq_max, 2e-3, 1e-3, 1e1}};
goal = cantilever.goalForceResolution;

c_opt = c.optimize_performance_from_current(parameter_constraints, ...
    nonlinear_constraints, goal);
c_opt.print_performance();
```

See `sampleCode.m` for additional examples including diffusion and ion implantation models.

## Testing

Run all tests:

```matlab
cd matlab
addpath('PiezoD');
runtests('tests')
```

Run a specific test file:

```matlab
runtests('tests/cantileverEpitaxyTest.m')
```

Run with verbose output:

```matlab
runtests('tests', 'Verbosity', 3)
```
