% sample_code.m
% Demonstrate the main features of the PiezoD cantilever modeling and optimization code
clear; clear; clear all; close all; clc; % Multiple clears often speed up optimization

% Turn on multiprocessor support. Substantially speeds up optimization if you have a multi-core machine.
if matlabpool('size') == 0
  matlabpool
end

% Generate an epitaxial cantilever, calculate its performance, then optimize
freq_min = 1e0;
freq_max = 1e3;
l = 100e-6;
w = 10e-6;
t = 1e-6;
l_pr_ratio = 0.1;
v_bridge = 2;
doping_type = 'phosphorus';
concentration = 1e20;
t_pr_ratio = 0.33;

c_epitaxy = cantilever_epitaxy(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, concentration, t_pr_ratio);

% Several ways to print the performance
c_epitaxy.fluid = 'vacuum'; % other options: 'air', 'water', 'arbitrary'
c_epitaxy.print_performance();
c_epitaxy.print_performance_for_excel(); % tab delimited, ready to copy and paste into a spreadsheet

% Plot the cantilever performance
c_epitaxy.plot_noise_spectrum();
c_epitaxy.plot_thermal_conductivity_profile(); 
c_epitaxy.plotTempProfile();
c_epitaxy.plotDopantBending();

% Setup optimization constraints
% HELPFUL TRICK
% If you know the final thickness you'd like, it's much faster to:
% 1) precompute the thermal conductivity for the thickness and approximate doping
% 2) add it to cantilever.k_x()
% 3) set min_t and max_t to the value you plan to use
parameter_constraints = {{'min_t', 'max_t', 'min_w', 'max_v_bridge'}, {1e-6, 1e-6, 5e-6, 10}};
nonlinear_constraints = {{'omega_min_hz', 'max_power', 'min_k', 'max_k'}, {freq_max, 2e-3, 1e-3, 1e1}};
goal = cantilever.goalForceResolution;

% Optimize once from the current cantilever design (will not consistently find the global minimum)
c_epitaxy = c_epitaxy.optimize_performance_from_current(parameter_constraints, nonlinear_constraints, goal);

% Optimize from random initial guesses until the answer converges (so that we're reasonably sure we're at the global minimum)
c_epitaxy = c_epitaxy.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
c_epitaxy.print_performance();


% Diffusion Example
diffusion_time = 20*60; % 20 minutes
diffusion_temp = 800 + 273; % at 800C (1073K)
c_diffusion = cantilever_diffusion(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, diffusion_time, diffusion_temp);
c_diffusion.fluid = 'vacuum';
c_diffusion.number_of_piezoresistors = 4;


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
c_epitaxy = cantilever_epitaxy(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, concentration_initial, t_pr_ratio);
c_epitaxy.number_of_piezoresistors = 1;
c_epitaxy.print_performance();
