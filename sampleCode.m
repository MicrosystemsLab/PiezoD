% sample_code.m
% This file demonstrates some of the basic features of piezoD
% and can also be used to verify that you have properly installed
% the code on your system.
%
% Check out cantilever.m for the many additional features included
% in the code. Examples include:
% - arbitrary parameter constraints (e.g. max temperature, diffusion temp)
% - piezoelectric and thermal actuator modeling
% - initial device deflections from instrinsic stress (Monte Carlo)
% - alternative optimization goals (e.g. surface stress)
% - alternative instrumentation amplifier noise characteristics

clear all
close all
clc

% Turn on multiprocessor support (if available)
if matlabpool('size') == 0
  matlabpool
end

% Compare the model with Harley's 89 nm thick epitaxial cantilever (1999)
freq_min = 10;
freq_max = 1000;
l = 300e-6;
w = 44e-6;
t = 89e-9;
l_pr_ratio = 45/300;
v_bridge = 5;
doping_type = 'boron';
concentration_initial = 4e19;
t_pr_ratio = 30/89;
exampleCantilever = cantileverEpitaxy(freq_min, freq_max, l, w, t, l_pr_ratio, ...
	v_bridge, doping_type, concentration_initial, t_pr_ratio);
exampleCantilever.fluid = 'vacuum'; % other options: 'air', 'water', 'arbitrary'
exampleCantilever.thermal_modeling = 'approx';
exampleCantilever.number_of_piezoresistors = 1;
exampleCantilever.print_performance();
exampleCantilever.print_performance_for_excel(); % tab delimited output
[TMax, TTip] = exampleCantilever.calculateMaxAndTipTemp();

% Setup optimization constraints
parameter_constraints = {{'min_t', 'max_t', 'min_w', 'max_v_bridge'}, ...
	{1e-6, 1e-6, 5e-6, 10}};
nonlinear_constraints = {{'omega_min_hz', 'max_power', 'min_k', 'max_k'}, ...
	{freq_max, 2e-3, 1e-3, 1e1}};
goal = cantilever.goalForceResolution;

% Optimize once from the current cantilever design
% Use optimize_performance() to randomize the optimization starting point
exampleCantilever = exampleCantilever.optimize_performance_from_current(parameter_constraints, ...
	nonlinear_constraints, goal);
exampleCantilever.print_performance();
exampleCantilever.plot_noise_spectrum();

% Diffusion Example
diffusion_time = 20*60; % 20 minutes
diffusion_temp = 800 + 273; % at 800C (1073K)
exampleCantilever = cantileverDiffusion(freq_min, freq_max, ...
    l, w, t, l_pr_ratio, v_bridge, doping_type, ...
    diffusion_time, diffusion_temp);
exampleCantilever.fluid = 'vacuum';
exampleCantilever.number_of_piezoresistors = 4;
exampleCantilever.print_performance();

% Ion implantation example
annealing_time = 20*60; % 20 minutes
annealing_temp = 950 + 273; % 950C
implantation_energy = 30; % keV
implantation_dose = 5e15; % 5e15/sq. cm
annealing_type = 'inert'; % or 'oxide'
exampleCantilever = cantileverImplantation(freq_min, freq_max, ...
        l, w, t, l_pr_ratio, v_bridge, doping_type, ...
        annealing_time, annealing_temp, annealing_type, ...
        implantation_energy, implantation_dose);
exampleCantilever.fluid = 'water';
exampleCantilever.number_of_piezoresistors = 4;
exampleCantilever.print_performance();
exampleCantilever.plotTempProfile();

% Plot the ion implanted dopant concentration profile
[x, active_doping, total_doping] = c_exampleCantilever.doping_profile();
figure
plot(x, active_doping);
set(gca, 'yscale', 'log');
xlabel('Depth (um)');
ylabel('Concentration (per cc)');

