clear all

% Constraints
freq_min = 1e0;
freq_max = 100e3;
omega_min_hz = 2*freq_max;
fluid_type = 'vacuum';
max_power = 2e-3;

% Bounds
constraints = {{'min_t', 'min_w', 'max_v_bias'}, ...
               {340e-9,  1e-6,    5}};

% Initial values
l = 30e-6;
w = 4e-6;
t = 1e-6;
l_pr_ratio = 0.5;
v_bias = 2;
doping_type = 'phosphorus';
concentration_initial = 1e18;
t_pr_ratio = 0.1;

% Create the cantilever and print the initial performance
c_epitaxy = cantilever_epitaxy(freq_min, freq_max, ...
    l, w, t, l_pr_ratio, v_bias, doping_type, concentration_initial, ...
    t_pr_ratio);

c_epitaxy.print_performance();

% Optimize and print the improved performance
c_epitaxy = c_epitaxy.optimize_performance(max_power, omega_min_hz, fluid_type, constraints);
c_epitaxy.print_performance();

c_epitaxy.plot_noise_spectrum();
% 
% % Diffusion Example
% diffusion_time = 20*60; % 20 minutes
% diffusion_temp = 800 + 273; % at 800C (1073K)
% 
% c_diffusion = cantilever_diffusion(freq_min, freq_max, ...
%     l, w, t, l_pr_ratio, v_bias, doping_type, ...
%     diffusion_time, diffusion_temp);
% c_diffusion.print_performance();
% 
% c_diffusion = c_diffusion.optimize_performance(max_power, omega_min_hz, fluid_type, constraints);
% c_diffusion.print_performance();