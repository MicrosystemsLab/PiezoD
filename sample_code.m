clear all

% Constraints
freq_min = 1e0;
freq_max = 1e3;
omega_min_hz = 5*freq_max;
fluid_type = 'vacuum';
max_power = 2e-3;

% Bounds
constraints = {{'min_t', 'min_w', 'max_v_bridge'}, ...
               {340e-9,  1e-6,    5}};

% Initial values
l = 30e-6;
w = 4e-6;
t = 1e-6;
l_pr_ratio = 0.5;
v_bridge = 2;
doping_type = 'phosphorus';
concentration_initial = 1e18;
t_pr_ratio = 0.33;

% Create the cantilever and print the initial performance
c_epitaxy = cantilever_epitaxy(freq_min, freq_max, ...
    l, w, t, l_pr_ratio, v_bridge, doping_type, concentration_initial, ...
    t_pr_ratio);
c_epitaxy.print_performance();

% Optimize the cantilever: use three initial random starting conditions
% in order to ensure convergence
c_epitaxy = c_epitaxy.optimize_performance(max_power, omega_min_hz, fluid_type, constraints);
c_epitaxy.print_performance();

% Or, just optimize from the current value. This is significantly faster
% but < 1% of the time (generally) won't converge and will return a
% suboptimal design.
c_epitaxy = c_epitaxy.optimize_performance_from_current(max_power, omega_min_hz, fluid_type, constraints);

% Plot the noise spectrum
c_epitaxy.plot_noise_spectrum();

% Diffusion Example
diffusion_time = 20*60; % 20 minutes
diffusion_temp = 800 + 273; % at 800C (1073K)

c_diffusion = cantilever_diffusion(freq_min, freq_max, ...
    l, w, t, l_pr_ratio, v_bridge, doping_type, ...
    diffusion_time, diffusion_temp);
c_diffusion.print_performance();

c_diffusion = c_diffusion.optimize_performance_from_current(max_power, omega_min_hz, fluid_type, constraints);
c_diffusion.print_performance();