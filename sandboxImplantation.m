clear
clear
clear all
close all
clc

% Base parameters
freq_min = 1e0;
freq_max = 1e3;
l = 800e-6;
w = 30e-6;
t = 2e-6;
l_pr_ratio = 0.3;
v_bridge = 2;
doping_type = 'boron'; %'arsenic'; %'boron'; %'phosphorus';
annealing_time = 30*60; % 15-120min
annealing_temp = 1050 + 273; % 900-1100C
annealing_type = 'inert'; % 'inert' or 'oxide'
implantation_energy = 50; % 20-80keV
implantation_dose = 1e16; % 2e14-2e16/cm2

c_implantation = cantilever_implantation(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, ...
  annealing_time, annealing_temp, annealing_type, implantation_energy, implantation_dose);
c_implantation.fluid = 'air';
c_implantation.number_of_piezoresistors = 4;

c_implantation.print_performance();

goal = cantilever.goalForceResolution;
parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge'}, {1e-6, 20e-6, 10}};
nonlinear_constraints = {{'omega_min_hz', 'max_power'}, {5e3, .1e-3}};
c_implantation = c_implantation.optimize_performance_from_current(parameter_constraints, nonlinear_constraints, goal);
