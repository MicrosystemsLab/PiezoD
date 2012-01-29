clear all
close all
clc
Constants();

if matlabpool('size') == 0
  matlabpool
end

% ===========================================================
% ================== Compute predicted resolution for NMC cantilever ==========================
% ===========================================================

% Summary:
% - Model matches well for the as fabricated design (25 pN and 4 nm RMS noise)
% - The best possible (42 kHz, 1 Hz - 10 kHz) is about 12 pN - partly due to improved possible gamma
% - The step displacement noise is higher due to the contact noise
% - Sader model predicts f0 drops from 42 kHz to 8 kHz in water, which I haven't verified
% - the 300 Hz step resonance is likely the stage but could plausibly be the cantilever, which would be bad

freq_min = 1;
freq_max = 10e3;
l = 95e-6;
w = 6e-6;
t = 320e-9;
l_pr_ratio = 53/95;
v_bridge = 2.048;
doping_type = 'phosphorus';
diffusion_time = 15*60;
diffusion_temp = 850 + 273;
c_diffusion = cantilever_diffusion(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, diffusion_time, diffusion_temp);
c_diffusion.foobar()
c_diffusion.dopantNumber()
return
c_diffusion.fluid = 'air';
c_diffusion.print_performance()
c_diffusion.fluid = 'water';
c_diffusion.print_performance()



% what is theoretically possible?
goal = cantilever.goalForceResolution;
parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge'}, {320e-9,  6e-6,    10}};
nonlinear_constraints = {{'omega_min_hz', 'max_power'}, {42e3, .1e-3}};

c_diffusion.fluid = 'vacuum';
c_diffusion = c_diffusion.optimize_performance(parameter_constraints, nonlinear_constraints, goal);


c_diffusion.fluid = 'water';
c_diffusion = c_diffusion.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
c_diffusion.fluid = 'air';
c_diffusion.print_performance()
c_diffusion.fluid = 'water';
c_diffusion.print_performance()
return

% ===========================================================
% =========== Compare Blom and Sader Q methods ==============
% ===========================================================

freq_min = 1;
freq_max = 1e3;
l = 262e-6;
w = 5e-6;
t = 1e-6;
l_pr_ratio = 35/262;
v_bridge = 3.8;
doping_type = 'boron';
concentration_initial = 4.4e19;
t_pr_ratio = 0.224;
c_epitaxy = cantilever_epitaxy(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, concentration_initial, t_pr_ratio);


c_epitaxy = c_epitaxy.updateFluid('water');
c_epitaxy.print_performance()
sqrt(c_epitaxy.thermo_PSD(1e3))

c_epitaxy = c_epitaxy.updateFluid('air');
c_epitaxy.print_performance()
sqrt(c_epitaxy.thermo_PSD(1e3))






% ===========================================================
% ================== Thermal Plots ==========================
% ===========================================================

% Compare single designs
freq_min = 1;
freq_max = 10e3;
min_thickness = 0.3e-6;
min_width = 2e-6;

l = 30e-6;
w = 4e-6;
t = 1e-6;
l_pr_ratio = 0.5;
v_bridge = 2;
doping_type = 'phosphorus';
concentration_initial = 1e18;
t_pr_ratio = 0.33;
c_epitaxy = cantilever_epitaxy(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, concentration_initial, t_pr_ratio);
c_epitaxy.number_of_piezoresistors = 4;
goal = cantilever.goalForceResolution;

low_power = 1e-4;
high_power = 1e-3;
parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge'}, {min_thickness,  min_width,    10}};
nonlinear_constraints = {{'omega_min_hz', 'fluid_type', 'max_power'}, {2*freq_max    , cantilever.fluidVacuum    , high_power}};
c_epitaxy_highpower = c_epitaxy.optimize_performance(parameter_constraints, nonlinear_constraints, goal);

parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge'}, {min_thickness,  min_width,    10}};
nonlinear_constraints = {{'omega_min_hz', 'fluid_type', 'max_power'}, {2*freq_max    , cantilever.fluidVacuum    , low_power}};
c_epitaxy_lowpower = c_epitaxy.optimize_performance(parameter_constraints, nonlinear_constraints, goal);

max_delta_temp = 10;
parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge'}, {min_thickness,  min_width,    10}};
nonlinear_constraints = {{'omega_min_hz', 'fluid_type', 'tip_temp'}, {2*freq_max, cantilever.fluidVacuum, max_delta_temp}};
c_epitaxy_temp = c_epitaxy.optimize_performance(parameter_constraints, nonlinear_constraints, goal);

parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge'}, {min_thickness,  min_width,    10}};
nonlinear_constraints = {{'omega_min_hz', 'fluid_type', 'tip_temp', 'max_power'}, {2*freq_max, cantilever.fluidVacuum, max_delta_temp, high_power}};
c_epitaxy_highpower_temp = c_epitaxy.optimize_performance(parameter_constraints, nonlinear_constraints, goal);


clc
c_epitaxy_highpower.print_performance_for_excel()
c_epitaxy_lowpower.print_performance_for_excel()
c_epitaxy_temp.print_performance_for_excel()
c_epitaxy_highpower_temp.print_performance_for_excel()


% Look at resolution for a fixed power dissipation

freq_min = 1;
freq_max = 100e3;
min_thickness = 0.3e-6;
min_width = 2e-6;

l = 30e-6;
w = 4e-6;
t = 1e-6;
l_pr_ratio = 0.5;
v_bridge = 2;
doping_type = 'phosphorus';
concentration_initial = 1e18;
t_pr_ratio = 0.33;
c_epitaxy = cantilever_epitaxy(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, concentration_initial, t_pr_ratio);
c_epitaxy.number_of_piezoresistors = 4;
goal = cantilever.goalForceResolution;

% Look at power controlled
points_per_decade = 2;
max_power_range = logspace(-6, -3, 4*points_per_decade+1);
for ii = 1:length(max_power_range)
    max_power = max_power_range(ii);

    parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge'}, {min_thickness,  min_width,    10}};
    nonlinear_constraints = {{'omega_min_hz', 'fluid_type', 'max_power'}, {2*freq_max    , cantilever.fluidVacuum    , max_power}};

    c_epitaxy = c_epitaxy.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
    mdf(ii) = c_epitaxy.force_resolution();
    temp(ii) = c_epitaxy.approxTipDeltaTemp();
end

figure
plot(max_power_range*1e3, mdf*1e12, 'LineWidth', 2);
set(gca, 'LineWidth', axisLineWidth, 'FontSize', axisFontSize);
set(gca, 'xscale', 'log');
xlabel('Cantilever Power Dissipation (mW)', 'FontSize', labelFontSize);
ylabel('MDF (pN)', 'FontSize', labelFontSize);
print([num2str(freq_max/1000) 'kHz_Power_Force'], '-depsc');

figure
plot(max_power_range*1e3, temp, 'LineWidth', 2);
set(gca, 'LineWidth', axisLineWidth, 'FontSize', axisFontSize);
set(gca, 'xscale', 'log', 'yscale', 'log');
xlabel('Cantilever Power Dissipation (mW)', 'FontSize', labelFontSize);
ylabel('Tip \DeltaT (C)', 'FontSize', labelFontSize);
print([num2str(freq_max/1000) 'kHz_Power_Temp'], '-depsc');


clear mdf temp
max_temp_range = logspace(0, 2, 3*points_per_decade+1);
for ii = 1:length(max_temp_range)
    max_delta_temp = max_temp_range(ii);

    parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge'}, {min_thickness,  min_width,    10}};
    nonlinear_constraints = {{'omega_min_hz', 'fluid_type', 'tip_temp'}, {2*freq_max, cantilever.fluidVacuum, max_delta_temp}};

    c_epitaxy = c_epitaxy.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
    mdf(ii) = c_epitaxy.force_resolution();
    temp(ii) = c_epitaxy.approxTipDeltaTemp();
    power_dissipation(ii) = c_epitaxy.power_dissipation();
end

figure
plot(max_temp_range, mdf*1e12, 'LineWidth', 2);
set(gca, 'LineWidth', axisLineWidth, 'FontSize', axisFontSize);
set(gca, 'xscale', 'log');
xlabel('Tip \DeltaT (C)', 'FontSize', labelFontSize);
ylabel('MDF (pN)', 'FontSize', labelFontSize);
print([num2str(freq_max/1000) 'kHz_Temp_Force'], '-depsc');

figure
plot(max_temp_range, power_dissipation*1e3, 'LineWidth', 2);
set(gca, 'LineWidth', axisLineWidth, 'FontSize', axisFontSize);
set(gca, 'xscale', 'log', 'yscale', 'log');
xlabel('Tip \DeltaT (C)', 'FontSize', labelFontSize);
ylabel('Power Dissipation (mW)', 'FontSize', labelFontSize);
print([num2str(freq_max/1000) 'kHz_Temp_Power'], '-depsc');




% ========================================================================
% ================== Compare thin film and c-Si ==========================
% ========================================================================

freq_min = 1e0;
freq_max = 1e4;

min_w = 2e-6;
max_power = 1e-3;
max_v_bridge = 10;
min_t = 100e-9; % for the cSi design
min_t_top = 100e-9;
min_t_mid = 20e-9;
min_t_bot = 100e-9;
parameter_constraints = {{'min_w', 'max_v_bridge', 'min_t', 'min_t_top', 'min_t_mid', 'min_t_bot'}, {min_w, max_v_bridge, min_t, min_t_top, min_t_mid, min_t_bot}};
nonlinear_constraints = {{'omega_min_hz', 'fluid_type', 'max_power'}, {2*freq_max, cantilever.fluidVacuum, max_power}};
goal = thin_film_cantilever.goalForceResolution;

l = 100e-6;
w = min_w;
t = min_t;
l_pr_ratio = 0.5;
v_bridge = 2;

t_top = 200e-9;
t_mid = 100e-9;
t_bot = 200e-9;
matl_top = 'al';
matl_mid = 'oxide';
matl_bot = 'al';
dopant_concentration = 1e19;

% Used for diffusion
doping_type = 'phosphorus';
diffusion_time = 15*60;
diffusion_temp = 875 + 273;


c_al_single   = thin_film_cantilever(freq_min, freq_max, l, w, l_pr_ratio, v_bridge, t_top, t_mid, t_bot, 'al', 'oxide', 'al', dopant_concentration);
c_ti_single   = thin_film_cantilever(freq_min, freq_max, l, w, l_pr_ratio, v_bridge, t_top, t_mid, t_bot, 'ti', 'oxide', 'ti', dopant_concentration);
c_poly_single = thin_film_cantilever(freq_min, freq_max, l, w, l_pr_ratio, v_bridge, t_top, t_mid, t_bot, 'poly', 'oxide', 'poly', dopant_concentration);
c_al_single.number_of_piezoresistors_on_cantilever = 1;
c_ti_single.number_of_piezoresistors_on_cantilever = 1;
c_poly_single.number_of_piezoresistors_on_cantilever = 1;
c_al_single.number_of_piezoresistors = 2;
c_ti_single.number_of_piezoresistors = 2;
c_poly_single.number_of_piezoresistors = 2;

c_al_double   = thin_film_cantilever(freq_min, freq_max, l, w, l_pr_ratio, v_bridge, t_top, t_mid, t_bot, 'al', 'oxide', 'al', dopant_concentration);
c_ti_double   = thin_film_cantilever(freq_min, freq_max, l, w, l_pr_ratio, v_bridge, t_top, t_mid, t_bot, 'ti', 'oxide', 'ti', dopant_concentration);
c_poly_double = thin_film_cantilever(freq_min, freq_max, l, w, l_pr_ratio, v_bridge, t_top, t_mid, t_bot, 'poly', 'oxide', 'poly', dopant_concentration);
c_al_double.number_of_piezoresistors_on_cantilever = 2;
c_ti_double.number_of_piezoresistors_on_cantilever = 2;
c_poly_double.number_of_piezoresistors_on_cantilever = 2;
c_al_double.number_of_piezoresistors = 2;
c_ti_double.number_of_piezoresistors = 2;
c_poly_double.number_of_piezoresistors = 2;

c_diffusion = cantilever_diffusion(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, diffusion_time, diffusion_temp);
c_diffusion.number_of_piezoresistors = 2; % or 4

c_poly_single = c_poly_single.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
c_al_single = c_al_single.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
c_ti_single = c_ti_single.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
c_poly_double = c_poly_double.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
c_al_double = c_al_double.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
c_ti_double = c_ti_double.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
c_diffusion = c_diffusion.optimize_performance(parameter_constraints, nonlinear_constraints, goal);

c_al_single.print_performance_for_excel()
c_ti_single.print_performance_for_excel()
c_poly_single.print_performance_for_excel()
c_al_double.print_performance_for_excel()
c_ti_double.print_performance_for_excel()
c_poly_double.print_performance_for_excel()
c_diffusion.print_performance_for_excel()

% % Constraints
% freq_min = 1e0;
% freq_max = 2e3;
% 
% k_range = [.1 10];
% 
% % Define optimization parameters
% parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge'}, ...
%                {340e-9,  2e-6,    5}};
% nonlinear_constraints = {{'omega_min_hz', 'fluid_type', 'max_power'}, ...
%                          {5*freq_max    , VACUUM    , 2e-3}};
% goal = FORCE_RESOLUTION;
% 
% % Initial values
% l = 30e-6;
% w = 4e-6;
% t = 1e-6;
% l_pr_ratio = 0.5;
% v_bridge = 2;
% doping_type = 'boron';

% % Diffusion Example
% diffusion_time = 20*60; % 20 minutes
% diffusion_temp = 800 + 273; % at 800C (1073K)
% 
% c_diffusion = cantilever_diffusion(freq_min, freq_max, ...
%     l, w, t, l_pr_ratio, v_bridge, doping_type, ...
%     diffusion_time, diffusion_temp);
% 
% c_diffusion = c_diffusion.optimize_performance_from_current(parameter_constraints, nonlinear_constraints, goal);
% c_diffusion.print_performance()
% c_diffusion.resolution_tradeoff_plot(parameter_constraints, nonlinear_constraints, k_range)



% % E341 Designs
% freq_min = .1;
% freq_max = 160e3;
% 
% k_range = [.1 100];
% min_thickness = 6e-6;
% min_width = 12e-6;
% max_power = 1e-3;
% diffusion_time = 15*60;
% diffusion_temp = 875 + 273;
% 
% parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge'}, {min_thickness,  min_width,    10}};
% 
% 
% parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge', 'min_diffusion_temp', 'max_diffusion_temp', 'min_diffusion_time', 'max_diffusion_time'}, ...
%                          {min_thickness,  min_width,    10, diffusion_temp, diffusion_temp, diffusion_time, diffusion_time}};
% 
% nonlinear_constraints = {{'omega_min_hz', 'fluid_type', 'max_power'}, {2*freq_max    , VACUUM    , max_power}};
% goal = FORCE_RESOLUTION;
% 
% l = 100e-6;
% w = min_width;
% t = min_thickness;
% l_pr_ratio = 0.5;
% v_bridge = 5;
% doping_type = 'phosphorus';
% 
% c_diffusion = cantilever_diffusion(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, diffusion_time, diffusion_temp);
% c_diffusion = c_diffusion.optimize_performance_from_current(parameter_constraints, nonlinear_constraints, goal);
% % c_diffusion = c_diffusion.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
% c_diffusion.print_performance()
% c_diffusion.print_performance_for_excel()




% % Double-check the E341 designs
% freq_min = .1;
% freq_max = 10e3;
% 
% k_range = [.00001 100];
% min_thickness = 6e-6;
% min_wirichterdth = 12e-6;
% max_power = 1e-3;
% diffusion_time = 15*60;
% diffusion_temp = 875 + 273;
% 
% l = 300e-6;
% w = min_width;
% t = min_thickness;
% l_pr_ratio = 68/300;
% v_bridge = 1;
% doping_type = 'phosphorus';
% 
% c_diffusion = cantilever_diffusion(freq_min, freq_max, l, w, t, l_pr_ratio, ...
%     v_bridge, doping_type, diffusion_time, diffusion_temp);
% c_diffusion.print_performance_for_excel()
% 
% parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge', 'min_diffusion_temp', 'max_diffusion_temp', 'min_diffusion_time', 'max_diffusion_time'}, ...
%                          {min_thickness,  min_width,    10, diffusion_temp, diffusion_temp, diffusion_time, diffusion_time}};
% nonlinear_constraints = {{'omega_min_hz', 'fluid_type', 'max_power'}, {2*freq_max    , cantilever.fluidVacuum    , max_power}};
% goal = cantilever.goalForceResolution;
% c_diffusion = c_diffusion.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
% c_diffusion.print_performance_for_excel()
% c_diffusion.plot_noise_spectrum()
% 
% c_diffusion.force_resolution()




% % Compare 2 and 4 PRs
% freq_min = 1;
% freq_max = 10e3;
% 
% min_thickness = 300e-9;
% min_width = 2e-6;
% max_power = 1e-3;
% 
% l = 300e-6;
% w = min_width;
% t = min_thickness;
% l_pr_ratio = 68/300;
% v_bridge = 1;
% doping_type = 'phosphorus';
% diffusion_time = 15*60;
% diffusion_temp = 875 + 273;
% 
% 
% c_diffusion = cantilever_diffusion(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, diffusion_time, diffusion_temp);
% 
% parameter_constraints = {{'min_t', 'min_w', 'max_v_bridge'}, {min_thickness,  min_width,    10}};
% nonlinear_constraints = {{'omega_min_hz', 'fluid_type', 'max_power'}, {2*freq_max, cantilever.fluidVacuum, max_power}};
% goal = cantilever.goalForceResolution;
% 
% c_diffusion.number_of_piezoresistors = 2;
% c_diffusion = c_diffusion.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
% c_diffusion.print_performance()
% c_diffusion.print_performance_for_excel()
% pause
% 
% c_diffusion.number_of_piezoresistors = 4;
% c_diffusion = c_diffusion.optimize_performance(parameter_constraints, nonlinear_constraints, goal);
% c_diffusion.print_performance()
% c_diffusion.print_performance_for_excel()