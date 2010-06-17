% sample_code.m
% PiezoD examples that demonstrate the main features of the code

clear
clear all
clc

% Generate an epitaxial cantilever, calculate its performance, then optimize
freq_min = 1e0;
freq_max = 1e3;
l = 30e-6;
w = 4e-6;
t = 1e-6;
l_pr_ratio = 0.5;
v_bridge = 2;
doping_type = 'phosphorus';
concentration = 1e18;
t_pr_ratio = 0.33;

c_epitaxy = cantilever_epitaxy(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, concentration, t_pr_ratio);
c_epitaxy.print_performance();
c_epitaxy.print_performance_for_excel();
pause

parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge'}, {1e-6, 1e-6, 5}};
nonlinear_constraints = {{'omega_min_hz', 'fluid_type', 'max_power', 'min_k', 'max_k'}, ...
                         {5*freq_max, cantilever.fluidVacuum, 2e-3, .1, 10}};
goal = cantilever.goalForceResolution;
c_epitaxy = c_epitaxy.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
c_epitaxy.print_performance();
c_epitaxy.plot_noise_spectrum();
pause


% Diffusion Example
diffusion_time = 20*60; % 20 minutes
diffusion_temp = 800 + 273; % at 800C (1073K)

c_diffusion = cantilever_diffusion(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, diffusion_time, diffusion_temp);
c_diffusion.print_performance();
pause

c_diffusion = c_diffusion.optimize_performance_from_current(parameter_constraints, nonlinear_constraints, goal);
c_diffusion.print_performance();
pause



% Play around with other features
c_epitaxy.number_of_piezoresistors = 4;
nonlinear_constraints = {{'omega_min_hz', 'fluid_type', 'delta_temp'}, {5*freq_max, cantilever.fluidWater, 10}};
c_epitaxy = c_epitaxy.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
pause




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