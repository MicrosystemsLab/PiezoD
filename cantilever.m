% General Notes
% * Usage
%       Create a cantilever
%           c = cantilever_diffusion(1e0, 1e3, 100e-6, 50e-6, 10e-6, 0.5,
%           0.5, 0.5, 5, 'boron', 20*60, 1000)
%
%       Calculate some useful properties
%           resistance = c.resistance();
%           sensitivity = c.force_sensitivity();
%           [omega_d, Q] = c.omega_damped_hz_quality_factor();
%
%       Design a better cantilever
%           c_optimized = c.optimize_performance(3.5e-3, 5e3, 'water',
%               30e-9, 1e-6, 5)
%
%       More code is included in sample_code.m
%
%       The code accounts for dopant type, adjusting the elastic modulus
%       and piezoresistive coefficient accordingly. Phosphorus, boron and
%       arsenic are supported, although the piezoresistive coefficient
%       calculator is based upon Harley's fit to experimental data which is
%       only for Boron.
%
% * Optimization Scaling
%       fmincon and other l-bfgs-b optimization routines don't work
%       particularly well if there's a large difference in the size of the
%       parameters it is tweaking. After each iteration it finds the
%       sensitivity of the variable to optimize with respect to a change in
%       the other variables. If dopant concentration is 1e19 per cc and
%       voltage is 5V, a unit change for each variable will have very
%       different effects on the optimization output, and bfgs won't know
%       what to do.
%       This problem is solved by scaling the variables back to O(1)
%       before putting them into the optimization code. Careful steps are
%       taken in order to calculate the optimization goal accurately and
%       rescaling the parameters at the proper time between iterations.
%
% * Extending the Code
%       Cantilever is the abstract base class that implements most of the
%       details. In order to actually do anything, you need to implement
%       a base class which implements at least these methods:
%           doping_profile
%           optimization_scaling
%           cantilever_from_state
%           optimize_performance
%           optimization_constraints
%
%       You might also want to change or override, but it's optional
%           print_performance
%           print_performance_for_excel
%           optimization_constraints
%
%       In order to demonstrate how you might do this, and to design some
%       cantilevers, we have created and included cantilever_epitaxy and
%       cantilever_diffusion.
%
%       Right now the code is oriented towards optimization force
%       resolution, but it is easy to extend it to other variables. You
%       would just need to modify your optimize_performance function.
%
% * Assumptions
%       You're using the INA103 amplifier with a gain of at least 100. This
%       determines the input referred voltage noise of the amplifier.

classdef cantilever
  
  properties
    freq_min; % Hertz
    freq_max; % Hertz
    l; % overall cantilever length
    w; % overall cantilever width (total width of both legs)
    t; % overall cantilever thickness
    l_pr_ratio; % piezoresistor length ratio
    v_bridge; % Volts
    fluid;

    number_of_piezoresistors = 2;
    rms_actuator_displacement_noise = 1e-12; % m
    alpha = 1e-6; % unitless
    amplifier = 'INA103';
    doping_type = 'boron';
    
    % (Optional) Account for tip mass loading or a reinforcement step at the base
    tip_mass = 0;    
    l_step = 0;
    t_step = 0;
    w_step = 0;
  end
  
  % Can be referred to with cantilever.variableName
  properties (Constant)
    
    
    % Material constants
    T = 300; % kelvin
    k_b = 1.38e-23; % J/K
    k_b_eV = 8.617343e-5; % eV/K
    q = 1.60218e-19; % Coulombs
    
    numFrequencyPoints = 1000;

    % Fluid properties
    rho_water = 1e3; % kg/m^3
    eta_water = 0.9e-3; % Pa-sec
    rho_air = 1.2; % kg/m^3
    eta_air = 17e-6; % Pa-sec
        
    % Thermal properties
    k_si = 130; % W/m-K
    h_air = 1000; % W/m^2-K
    h_water = 20000; % W/m^2-k

    % Mechanical material properties
    E_Si = 130e9;
    rho_Si = 2330;
    nu_Si = 0.28;
    E_Al = 70e9;
    rho_Al = 2700;
    nu_Al = 0.35;
    maxQ = 1000;
    minQ = 1e-6;

    % Define the materials
    rho_cantilever = cantilever.rho_Si;
    rho_step = cantilever.rho_Al;
    E_step = cantilever.E_Al;
    
    % The optimization goals
    goalForceResolution = 0;
    goalDisplacementResolution = 1;
    
    % Store the Sader lookup table and vectors as a constant
    kappa_lookup = [0 0.125 0.25 0.5 0.75 1 2 3 5 7 10 20];
    reynolds_lookup = [-4 -3.5 -3 -2.5 -2 -1.5 -1 -.5 0 0.5 1 1.5 2 2.5 3 3.5 4];
    
    tau_lookup_real = ...
      [3919.41 59.3906  22.4062  9.13525  5.62175  4.05204  1.93036  1.2764   0.764081 0.545683 0.381972 0.190986;
      1531.90 59.3861  22.4061  9.13525  5.62175  4.05204  1.93036  1.2764   0.764081 0.545683 0.381972 0.190986;
      613.426 59.3420  22.4052  9.13523  5.62174  4.05204  1.93036  1.2764   0.764081 0.545683 0.381972 0.190986;
      253.109 58.9094  22.3962  9.13504  5.62172  4.05203  1.93036  1.2764   0.764081 0.545683 0.381972 0.190986;
      108.429 55.2882  22.3078  9.13319  5.62153  4.05199  1.93036  1.2764   0.764081 0.545683 0.381972 0.190986;
      48.6978 40.7883  21.5187  9.11481  5.61960  4.05160  1.93035  1.2764   0.764081 0.545683 0.381972 0.190986;
      23.2075 22.7968  17.5378  8.94370  5.60057  4.04771  1.93027  1.27639  0.76408  0.545683 0.381972 0.190986;
      11.8958 11.9511  11.0719  7.89716  5.43378  4.01051  1.92942  1.27629  0.764074 0.545682 0.381972 0.190986;
      6.64352 6.64381  6.47227  5.65652  4.64017  3.74600  1.92114  1.27536  0.764012 0.545671 0.38197  0.190986;
      4.07692 4.05940  3.99256  3.72963  3.37543  3.00498  1.85532  1.26646  0.763397 0.545564 0.381953 0.190985;
      2.74983 2.73389  2.69368  2.56390  2.39884  2.22434  1.61821  1.20592  0.757637 0.54452  0.381786 0.190981;
      2.02267 2.01080  1.98331  1.90040  1.79834  1.69086  1.31175  1.04626  0.721165 0.535593 0.38018  0.190932;
      1.60630 1.59745  1.57723  1.51690  1.44276  1.36416  1.08036  0.878177 0.635443 0.496169 0.368548 0.190459;
      1.36230 1.35532  1.33934  1.29142  1.23203  1.16842  0.932812 0.759965 0.551349 0.435586 0.334279 0.186672;
      1.21727 1.21141  1.19792  1.15718  1.10624  1.05117  0.84292  0.686229 0.493924 0.387183 0.295972 0.172722;
      1.13038 1.12518  1.11316  1.07668  1.03073  0.980721 0.78879  0.641744 0.458699 0.356289 0.268907 0.154450;
      1.07814 1.07334  1.06221  1.02827  0.985314 0.938346 0.756309 0.615164 0.437743 0.337813 0.252327 0.140852];
    
    tau_lookup_imag = ...
      [27984.8   44628.5     55176.1   71754     86311.5   100062    152411    203623    305570    407436    560225    1069521;
      9816.33   14113.1     17448.2   22690.6   27294.1   31642.3   48196.5   64391.4   96629.7   128843    177159    338212;
      3482.47   4464.16     5517.72   7175.41   8631.15   10006.2   15241.1   20362.3   30557     40743.6   56022.5   106952;
      1252.66   1415.42     1745.17   2269.09   2729.42   3164.23   4819.65   6439.14   9662.97   12884.3   17715.9   33821.2;
      458.386   457.863     552.862   717.635   863.138   1000.63   1524.112  2036.23   3055.7    4074.36   5602.25   10695.2;
      171.397   160.951     177.702   227.205   273.013   316.449   481.967   643.915   966.297   1288.43   1771.59   3382.12;
      65.8679   62.2225     61.626    72.6542   86.5364   100.144   152.418   203.625   305.57    407.436   560.225   1069.52;
      26.2106   25.21       24.1432   24.7484   27.9459   31.8957   48.2199   64.3973   96.6308  128.843   177.159   338.212;
      10.8983   10.6158     10.1909   9.7009    9.91067   10.648    15.3139   20.381    30.5604   40.7448   56.0229   106.952;
      4.78389   4.69492     4.53952   4.24925   4.09701   4.09433   5.01844   6.49605   9.67379   12.8879   17.7171   33.8214;
      2.23883   2.20681     2.14583   2.0088    1.89659   1.82463   1.85993   2.17718   3.08849   4.08581   5.60598   10.6956;
      1.12164   1.10851     1.08208   1.01654   0.953355  0.901676  0.81464   0.844519  1.04394   1.32116   1.78306   3.38349;
      0.596697  0.590686    0.578118  0.545082  0.510467  0.479247  0.403803  0.383595  0.409256  0.469688  0.589749  1.07377;
      0.332285  0.329276    1.32283   0.305262  0.285953  0.26763   0.216732  0.194409  0.186218  0.195634  0.221631  0.349855;
      0.191043  0.189434    0.185931  0.176166  0.165118  0.154323  0.122124  0.105573  0.0938839 0.0925686 0.09682   0.126835;
      0.112082  0.111181    0.109199  0.103595  0.0971392 0.0907188 0.0707736 0.059728  0.0505049 0.0476557 0.0471326 0.0534759;
      0.0665172 0.0659974   0.0648471 0.0615627 0.0577366 0.0538889 0.0416384 0.0345727 0.0282418 0.025856  0.024611  0.0252877];

  end
  
  methods (Abstract)
    doping_profile(self)
    optimization_scaling(self)
    cantilever_from_state(self, x0)
    initial_conditions_random(self)
    current_state(self)
    optimization_bounds(self, parameter_constraints)
  end
  
  methods
    
    function self = cantilever(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type)
      self.freq_min = freq_min;
      self.freq_max = freq_max;
      self.l = l;
      self.w = w;
      self.t = t;
      self.l_pr_ratio = l_pr_ratio;
      self.v_bridge = v_bridge;
      self.doping_type = doping_type;
      self.fluid = 'air'; % Default value           
    end
    
    % Calculate the actual dimensions (getter functions)
    function l_pr = l_pr(self)
      l_pr = self.l * self.l_pr_ratio;
    end
    
    function w_pr = w_pr(self)
      w_pr = self.w/2;
    end
    
    % ==================================
    % ========= Pretty output ==========
    % ==================================
    
    % fprintf performance
    function print_performance(self)
      
      % Calculate intermediate quantities
      [omega_damped_hz, Q] = self.omega_damped_hz_and_Q();
      [x, doping] = self.doping_profile();
      Nz = trapz(x, doping*1e6);
      
      [TMax_approx TTip_approx] = self.approxTempRise();
      [TMax, TTip] = self.calculateMaxAndTipTemp();
      
      thermoLimit = self.thermo_integrated()/self.force_sensitivity();
      
      fprintf('Freq range: %f to %f \n', self.freq_min, self.freq_max)
      fprintf('Operating fluid: %s \n', self.fluid);
      fprintf('Cantilever L/W/T: %f %f %f \n', self.l*1e6, self.w*1e6, self.t*1e6)
      fprintf('PR L/W: %f %f %f \n', self.l_pr()*1e6, self.w_pr()*1e6)
      fprintf('PR Length Ratio: %g \n', self.l_pr_ratio)
      fprintf('Number of silicon resistors: %f \n', self.number_of_piezoresistors)
      fprintf('\n')
      fprintf('Force resolution (N): %g \n', self.force_resolution())
      fprintf('Displacement resolution (m): %g \n', self.displacement_resolution())
      fprintf('Force Sensitivity (V/N) %g \n', self.force_sensitivity())
      fprintf('Displacement Sensitivity (V/m) %g \n', self.displacement_sensitivity())
      fprintf('Beta %g \n', self.beta())
      fprintf('Thermomechanical force noise limit: %g \n', thermoLimit);
      fprintf('\n')
      fprintf('Stiffness (N/m): %g \n', self.stiffness())
      fprintf('Vacuum freq: %f \n', self.omega_vacuum_hz())
      fprintf('Damped freq: %f \n', omega_damped_hz)
      fprintf('Quality factor: %f \n', Q)
      fprintf('\n')
      fprintf('Wheatstone bridge bias voltage: %f \n', self.v_bridge)
      fprintf('Resistance: %f \n', self.resistance())
      fprintf('Sheet Resistance: %f \n', self.sheet_resistance())
      fprintf('Power dissipation (mW): %g \n', self.power_dissipation()*1e3)
      fprintf('Approx. Temp Rises (C) - Tip: %f  Max: %f\n', TTip_approx, TMax_approx)
      fprintf('F-D Temp Rises (C)     - Tip: %f  Max: %f\n', TTip, TMax)
      fprintf('\n')
      fprintf('Integrated noise (V): %g \n', self.integrated_noise())
      fprintf('Integrated johnson noise (V): %g \n', self.johnson_integrated())
      fprintf('Integrated 1/f noise (V): %g \n', self.hooge_integrated())
      fprintf('Amplifier noise (V): %g \n', self.amplifier_integrated())
      fprintf('Thermomechanical noise (V): %g \n', self.thermo_integrated())
      fprintf('\n')
      fprintf('Johnson/Hooge: %g \n', self.johnson_integrated()/self.hooge_integrated())
      fprintf('Knee frequency (Hz): %g \n', self.knee_frequency())
      fprintf('Number of Carriers: %g \n', self.number_of_carriers());
      fprintf('Nz: %g \n', Nz)
      fprintf('\n')
    end
    
    function print_performance_for_excel(self)
      % Calculate intermediate quantities
      [omega_damped_hz, Q] = self.omega_damped_hz_and_Q();
      [TMax, TTip] = self.calculateMaxAndTipTemp();
      
      variables_to_print = [self.l*1e6, self.w*1e6, self.t*1e6, ...
        self.l_pr()*1e6, self.l_pr_ratio, ...
        self.v_bridge, self.freq_min, self.freq_max, ...
        self.force_resolution(), self.displacement_resolution(), ...
        self.omega_vacuum_hz(), omega_damped_hz, ...
        self.stiffness(), Q, self.force_sensitivity(), self.beta(), ...
        self.resistance(), self.power_dissipation()*1e3, TTip, ...
        self.integrated_noise(), self.integrated_johnson_noise(), ...
        self.integrated_hooge_noise(), self.knee_frequency()];
      
      for print_index = 1:length(variables_to_print)
        fprintf('%g \t', variables_to_print(print_index));
      end
      fprintf('\n')
    end
    
    % ==================================
    % ===== Calculate resistance =======
    % ==================================
    
    % Calculate total resistance of piezoresistor. Includes effect of other resistances (gamma)
    % Units: ohms
    function resistance = resistance(self)
      number_of_squares = self.resistor_length()/self.w_pr();
      resistance = number_of_squares * self.sheet_resistance() / self.gamma();
    end
    
    % Calculate resistor length, used to calculat resistance and number of carriers
    % Units: m
    function resistor_length = resistor_length(self)
      resistor_length = 2*self.l_pr();
    end
    
    % Calculate sheet resistance. Uses abstract method self.doping_profile() - which must be defined in a subclass.
    % Units: ohms
    function Rs = sheet_resistance(self)
      [x, doping] = self.doping_profile(); % x -> m, doping -> N/cm^3
      conductivity = self.conductivity(doping); % ohm-cm
      Rs = 1/trapz(x*1e2, conductivity); % convert x to cm
    end
    
    % Calculate conductivity for a given dopant concentration. Can use vectors or single values.
    % Units: C/V-sec-cm
    function sigma = conductivity(self, dopant_concentration)
      mu = self.mobility(dopant_concentration);
      sigma = mu.*self.q.*dopant_concentration;
    end
    
    % Data from "Modeling of Carrier Mobility Against Carrier Concentration in Arsenic-,  Phosphorus-,
    % and Boron-Doped  Silicon", Masetti, Serveri and Solmi - IEEE Trans. on Electron Devices, (1983)
    % Units: cm^2/V-sec
    function mobility = mobility(self, dopant_concentration)
      
      n = dopant_concentration;
      p = dopant_concentration;
      
      switch self.doping_type
        case 'arsenic'
          % Arsenic data
          mu_0 = 52.2;
          mu_max = 1417;
          mu_1 = 43.4;
          C_r = 9.96e16;
          C_s = 3.43e20;
          alpha = 0.680;
          beta = 2.0;
          mobility = mu_0 + (mu_max - mu_0)./(1 + (n./C_r).^alpha) - mu_1./(1 + (C_s./n).^beta);
          
        case 'phosphorus'
          % Phosphorus data
          mu_0 = 68.5;
          mu_max = 1414;
          mu_1 = 56.1;
          C_r = 9.2e16;
          C_s = 3.41e20;
          alpha = 0.711;
          beta = 1.98;
          mobility = mu_0 + (mu_max - mu_0)./(1 + (n./C_r).^alpha) - mu_1./(1 + (C_s./n).^beta);
          
        case 'boron'
          % Boron data
          mu_0 = 44.9;
          mu_max = 470.5;
          mu_1 = 29.0;
          C_r = 2.23e17;
          C_s = 6.1e20;
          alpha = 0.719;
          beta = 2.00;
          p_c = 9.23e16;
          mobility = mu_0.*exp(-p_c./n) + mu_max./(1 + (p./C_r).^alpha) - mu_1./(1 + (C_s./p).^beta);
          
      end
    end
    
    % ==================================
    % ======= Calculate noise ==========
    % ==================================
    
    % The number of current carriers in the piezoresistor. Integrate the carries to the junction depth
    % and multiply by the lateral dimensions of the piezoresistor.
    % Units: unitless
    function number_of_carriers = number_of_carriers(self)
      [x, doping] = self.doping_profile(); % Units: x -> m, doping -> N/cm^3
      Nz = trapz(x, doping*1e6); % doping: N/cm^3 -> N/m^3
      number_of_carriers = Nz*self.resistor_length()*self.w_pr();
    end
    
    % 1/f voltage power spectral density for the entire Wheatstone bridge
    % Units: V^2/Hz
    function hooge_PSD = hooge_PSD(self, freq)
      hooge_PSD = self.alpha*self.v_bridge^2*self.number_of_piezoresistors./(4*self.number_of_carriers()*freq);
    end
    
    % Integrated 1/f noise density for the entire Wheatstone bridge
    % Unit: V
    function hooge_integrated = hooge_integrated(self)
      freq = logspace( log10(self.freq_min), log10(self.freq_max), cantilever.numFrequencyPoints);
      hooge_integrated = sqrt(trapz(freq, self.hooge_PSD(freq)));
    end
    
    % Johnson noise PSD from the entire Wheatstone bridge. Equal to that of a single resistor
    % assuming that all four resistors are equal.
    % Units: V^2/Hz
    function johnson_PSD = johnson_PSD(self, freq)
      R_external = 700;
      johnson_PSD = 4*self.k_b*self.T*(self.resistance()/2 + R_external/2) * ones(1, length(freq));
    end
    
    % Integrated Johnson noise
    % Unit: V
    function johnson_integrated = johnson_integrated(self)
      freq = logspace( log10(self.freq_min), log10(self.freq_max), cantilever.numFrequencyPoints);
      johnson_integrated = sqrt(trapz(freq, self.johnson_PSD(freq)));
    end
    
    % Thermomechanical noise PSD
    % Units: V^2/Hz
    function thermo_PSD = thermo_PSD(self, freq)
      [omega_damped, Q_M] = self.omega_damped_and_Q();
      thermo_PSD = (self.force_sensitivity())^2 * 2*self.stiffness()*self.k_b*self.T/(pi*self.omega_vacuum_hz()*Q_M) * ones(1, length(freq));
    end
    
    % Integrated thermomechanical noise
    % Unit: V
    function thermo_integrated = thermo_integrated(self)
      freq = logspace( log10(self.freq_min), log10(self.freq_max), cantilever.numFrequencyPoints);
      thermo_integrated = sqrt(trapz(freq, self.thermo_PSD(freq)));
    end
    
    % Accounts for the noise of the actuator that the cantilever is mounted on
    % Units: V
    function actuator_noise_integrated = actuator_noise_integrated(self)
      actuator_noise_integrated = self.rms_actuator_displacement_noise*self.stiffness()*self.force_sensitivity(); % V
    end
    
    
    % Amplifier noise PSD
    % Units: V^2/Hz
    function amplifier_PSD = amplifier_PSD(self, freq)
      switch self.amplifier
        case 'INA103'
          A_VJ = 1.2e-9; % 1.2 nV/rtHz noise floor
          A_IJ = 2e-12; % 2 pA/rtHz noise floor
          A_VF = 6e-9; % 6 nV/rtHz @ 1 Hz
          A_IF = 25e-12; % 25 pA/rtHz @ 1 Hz
        case 'AD8221'
          A_VJ = 8e-9;
          A_IJ = 40e-15;
          A_VF = 12e-9;
          A_IF = 550e-15;
        otherwise
          fprintf('ERROR: UNKNOWN AMPLIFIER')
          pause
      end
      
      R_effective = self.resistance()/2; % resistance seen by amplifier inputs
      amplifier_PSD = (A_VJ^2 + 2*(R_effective*A_IJ)^2) + (A_VF^2 + 2*(R_effective*A_IF)^2)./freq;
    end
    
    % Integrated amplifier noise
    % Units: V
    function amplifier_integrated = amplifier_integrated(self)
      freq = logspace( log10(self.freq_min), log10(self.freq_max), cantilever.numFrequencyPoints);
      amplifier_integrated = sqrt(trapz(freq, self.amplifier_PSD(freq)));
    end
    
    % Calculate the knee frequency (equating the Hooge and Johnson noise)
    % Equating 1/f noise and johnson... numPRS*alpha*V_bridge^2/(4*N*f_knee) = 4*kb*T*R
    % Leads to f_knee = alpha*V_bridge^2/(16*N*S_j^2)
    % Units: Hz
    function knee_frequency = knee_frequency(self)
      knee_frequency = self.number_of_piezoresistors*self.alpha*self.v_bridge^2/(16*self.number_of_carriers()*self.k_b*self.T*self.resistance());
    end
    
    % Integrated cantilever noise for given bandwidth
    % Units: V
    function integrated_noise = integrated_noise(self)
      freq = logspace( log10(self.freq_min), log10(self.freq_max), cantilever.numFrequencyPoints);
      integrated_actuator_noise = self.actuator_noise_integrated();
      integrated_noise = sqrt(integrated_actuator_noise^2 + trapz(freq, self.johnson_PSD(freq) + self.hooge_PSD(freq) + self.thermo_PSD(freq) + self.amplifier_PSD(freq)));
    end
    
    % Calculate the noise in V/rtHz at a given frequency
    function voltage_noise = voltage_noise(self, freq)
      voltage_noise = sqrt(self.johnson_PSD(freq) + self.hooge_PSD(freq) + self.thermo_PSD(freq) + self.amplifier_PSD(freq));
    end
    
    function plot_noise_spectrum(self)
      freq = logspace( log10(self.freq_min), log10(self.freq_max), cantilever.numFrequencyPoints);
      noise = self.voltage_noise(freq);
      
      figure
      hold all
      plot(freq, sqrt(self.johnson_PSD(freq)), 'LineWidth', 2);
      plot(freq, sqrt(self.hooge_PSD(freq)), 'LineWidth', 2);
      plot(freq, sqrt(self.thermo_PSD(freq)), 'LineWidth', 2);
      plot(freq, sqrt(self.amplifier_PSD(freq)), 'LineWidth', 2);
      plot(freq, self.voltage_noise(freq), 'LineWidth', 2);
      hold off
      set(gca, 'xscale','log', 'yscale','log');
      set(gca, 'LineWidth', 1.5, 'FontSize', 14);
      ylabel('Noise Voltage Spectral Density (V/rtHz)', 'FontSize', 16);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      legend('Johnson', 'Hooge', 'Thermo', 'Amp', 'Total')
    end
    
    function f_min_cumulative = f_min_cumulative(self)
      frequency = logspace(log10(self.freq_min), log10(self.freq_max), cantilever.numFrequencyPoints);
      noise = self.voltage_noise(frequency);
      sensitivity = self.force_sensitivity();
      force_noise_density = noise./sensitivity;
      f_min_cumulative = sqrt(cumtrapz(frequency, force_noise_density.^2));
    end
    
    
    % ==================================
    % ===== Calculate sensitivity ======
    % ==================================
    
    
    % Piezoresistance factor. Accounts for dopant concentration dependent piezoresistivity in silicon
    % Uses Richter's 2008 model from "Piezoresistance in p-type silicon revisited" for the case of T=300K
    % Could be readily generalized to account for temperature as well
    function piezoresistance_factor = piezoresistance_factor(self, dopant_concentration)
      Nb = 6e19;
      Nc = 7e20;
      richter_alpha = 0.43;
      richter_gamma = 1.6;
      piezoresistance_factor = (1 + (dopant_concentration/Nb).^richter_alpha + (dopant_concentration/Nc).^richter_gamma).^-1;
    end
    
    function max_factor = max_piezoresistance_factor(self)
      switch self.doping_type
        case 'boron'
          max_factor = 72e-11; % Pi at low concentration in 110 direction
        case 'phosphorus'
          max_factor = 103e-11; %Pi at low concentration in 100 direction
        case 'arsenic'
          max_factor = 103e-11; %Pi at low concentration in 100 direction
      end
    end
    
    % Reduction in sensitivity from piezoresistor not located just at the surface.
    % Calculated for the general case of an arbitrarily shaped doping profile. Taken from Sung-Jin Park's Hilton Head 2008 paper.
    % Units: None
    function beta = beta(self)
      [x, doping_concentration] = self.doping_profile();
      
      % x is supposed to vary from t/2 to -t/2 as it varies from the top to bottom surface
      x = (self.t/2 - x)*1e2; % x: m -> cm
      
      mu = self.mobility(doping_concentration); % cm^2/V-s
      P = self.piezoresistance_factor(doping_concentration);
      
      numerator = trapz(x, self.q.*mu.*doping_concentration.*P.*x);
      denominator = trapz(x, self.q.*mu.*doping_concentration);
      beta = 2*numerator/(self.t*1e2*denominator); % t: m -> cm
      
      % Ensure that beta doesn't become too small or negative
      beta_epsilon = 1e-6;
      beta = max(beta, beta_epsilon);
    end
    
    % Ratio of piezoresistor resistance to total resistance (< 1)
    % Assume that metal interconnects are used and that their resistance is about 10% of the total
    function gamma = gamma(self)
      gamma = 0.9;
    end
    
    % Units: V/N
    % TODO: Account for the transverse current flow at the end of the cantilever
    function force_sensitivity = force_sensitivity(self)
      force_sensitivity = 3*(self.l - self.l_pr()/2)*self.max_piezoresistance_factor()/(2*self.w*self.t^2)*self.beta()*self.gamma()*self.v_bridge;
    end
    
    function displacement_sensitivity = displacement_sensitivity(self)
      displacement_sensitivity = self.force_sensitivity()*self.stiffness();
    end
    
    
    % ====================================
    % === Calculate thermal properties ===
    % ====================================
    
    % Power dissipation (W) in the cantilever
    function power_dissipation = power_dissipation(self)
      power_dissipation = (self.v_bridge/2)^2/self.resistance();
    end
    
    % Calculate the approximate max and tip temperatures using lumped modeling
    % This works significantly better than F-D for design optimization
    % Units: K
    function [TMax TTip] = approxTempRise(self)
      
      switch self.fluid
        case 'vacuum'
          h = 0;
        case 'air'
          h = self.h_air;
        case 'water'
          h = self.h_water;
      end
      
      exposedArea = 2*self.l*(self.w + self.t);
      R_base = self.l_pr()/(2*self.w*self.t*self.k_si);
      R_convection = 1/(2*pi*h*exposedArea); % Factor of 2*pi included for better agreement with F-D results
      R_total = 1/(1/R_base + 1/R_convection);
      
      power = self.power_dissipation();
      
      TMax = power*R_base;
      TTip = power*R_total;
    end
    
    % Model the temperature profile of a self-heated PR cantilever
    % - Assumes fixed temperature at the cantilever base, and adiabatic conditions at the tip, convection to ambient
    % - There is significant uncertainty in the convection coefficient
    
    % References used
    % - Finite differences: http://reference.wolfram.com/mathematica/tutorial/NDSolvePDE.html
    % - 1D Conduction Analysis: http://people.sc.fsu.edu/~burkardt/f_src/fd1d_heat_steady/fd1d_heat_steady.html
    function [x, T] = calculateTempProfile(self)
      n_points = 800;
      dx = self.l/(n_points - 1);
      x = 0:dx:self.l;
      power = (self.v_bridge/2)^2/self.resistance();
      l_pr = self.l_pr(); % call once to improve speed
      Qgen = power/l_pr;
      perimeter = 2*(self.w + self.t);
      
      tempBase    = 273 + 25;
      tempAmbient = 273 + 25;
      
      % Choose the convection coefficient based upon the ambient fluid
      switch self.fluid
        case 'vacuum'
          h = 0;
        case 'air'
          h = self.h_air;
        case 'water'
          h = self.h_water;
      end
      
      % Build A
      A_maindiagonal = 2*self.k_si*self.w*self.t/dx^2 + h*perimeter;
      A_offdiagonal = -self.k_si*self.w*self.t/dx^2;
      A = diag(ones(n_points, 1)  *A_maindiagonal, 0) + ...
        diag(ones(n_points-1, 1)*A_offdiagonal,  1) + ...
        diag(ones(n_points-1, 1)*A_offdiagonal, -1);
      A(1, 1) = 1; % Fixed temp at base
      A(1, 2) = 0; % Remove an extra entry
      A(n_points, n_points-1:n_points) = [1 -1]; % Adiabatic at tip
      A = sparse(A); % Use a sparse matrix because most of the entries are zero - convert after building it up
      
      % Generate the RHS matrix
      rhs = ones(n_points, 1)*h*perimeter*tempAmbient;
      rhs(x <= l_pr) = Qgen + h*perimeter*tempAmbient;
      rhs(1,1) = tempBase; % Fixed temp
      rhs(end,1) = 0; % Adiabatic
      
      % Solve and then return the temp rise relative to ambient
      T = A \ rhs;
      T = T - tempAmbient;
      
    end
    
    function [TMax, TTip] = calculateMaxAndTipTemp(self)
      [x, T] = self.calculateTempProfile();
      TMax = max(T);
      TTip = T(end);
    end
    
    % ==================================
    % ====== Calculate resolution ======
    % ==================================
    
    % Units: N
    function force_resolution = force_resolution(self)
      force_resolution = self.integrated_noise()/self.force_sensitivity();
    end
    
    function displacement_resolution = displacement_resolution(self)
      displacement_resolution = self.force_resolution()/self.stiffness();
    end
    
    % ==================================
    % ======== Beam mechanics ==========
    % ==================================
    
    % Calculate elastic modulus based upon dopant type. Assume we're using the best piezoresistor orientation.
    function elastic_modulus = modulus(self)
      switch self.doping_type
        case 'boron'
          elastic_modulus = 169e9; % <110> direction
        case 'phosphorus'
          elastic_modulus = 130e9; % <100> direction
        case 'arsenic'
          elastic_modulus = 130e9; % <100> direction
      end
    end
    
    % Bending stiffness of the cantilever to a point load at the tip
    % Units: N/m
    function stiffness = stiffness(self)
      stiffness = self.modulus() * self.w * self.t^3 / (4*self.l^3);
    end
    
    function effective_mass = effective_mass(self)
      cantilever_effective_mass = 0.243 * self.rho_cantilever * self.w * self.t * self.l;
      effective_mass = cantilever_effective_mass + self.tip_mass;
    end
    
    function [rho_fluid, eta_fluid] = lookupFluidProperties(self)
      switch self.fluid
        case 'air'
          rho_fluid = self.rho_air;
          eta_fluid = self.eta_air;
        case 'water'
          rho_fluid = self.rho_water;
          eta_fluid = self.eta_water;
        otherwise
          fprintf('ERROR - Unknown fluid: %s', self.fluid);
          pause
      end
    end
    
    % Resonant frequency for undamped vibration (first mode)
    % Units: radians/sec
    function omega_vacuum = omega_vacuum(self)
      omega_vacuum = sqrt( self.stiffness() / self.effective_mass());
    end
    
    % Resonant frequency for undamped vibration (first mode)
    % Units: cycles/sec
    function omega_vacuum_hz = omega_vacuum_hz(self)
      omega_vacuum_hz = self.omega_vacuum() / (2*pi);
    end
    
    % Calculate the damped natural frequency and Q, which we know lies between zero and the natural frequency in vacuum.
    % Per Eysden and Sader (2007)
    function [omega_damped, Q] = omega_damped_and_Q(self)
      
      % If we're in vacuum, just return the vacuum frequency
      switch self.fluid
        case 'vacuum'
          omega_damped = self.omega_vacuum();
          Q = cantilever.maxQ;
          return;
      end
      
      % Inner function for solving the transcendental equation to find omega_damped
      % We're searching for a function minimum, so return the residual squared (continuous and smooth)
      function residual_squared = find_natural_frequency(omega_damped)
        hydro = self.hydrodynamic_function(omega_damped, rho_f, eta_f);
        residual = omega_damped - self.omega_vacuum()*(1 + pi * rho_f * self.w/(4 * self.rho_cantilever * self.t) .* real(hydro)).^-0.5;
        residual_squared = residual^2;
      end
      
      % Lookup fluid properties once, then calculate omega_damped and Q
      [rho_f eta_f] = self.lookupFluidProperties();
      options = optimset('TolX', 1, 'TolFun', 1e-6, 'Display', 'off'); % Calculate omega_damped to within 1 radian/sec
      omega_damped = fminbnd(@find_natural_frequency, 0, self.omega_vacuum(), options);
      hydro = self.hydrodynamic_function(omega_damped, rho_f, eta_f);
      Q = (4 * self.rho_cantilever * self.t / (pi * rho_f * self.w) + real(hydro)) / imag(hydro);

      % Sometimes our initial guess will turn up Q = NaN because it's outside the bounds of the interpolation
      % Usually this is because the cantilever is shorter than it is wide, and it will be fixed after
      % a few iterations
      if Q < cantilever.minQ || isnan(Q)
        Q = cantilever.minQ;
      elseif Q > cantilever.maxQ
        Q = cantilever.maxQ;
      end
    end
    
    function [omega_damped_hz, Q] = omega_damped_hz_and_Q(self)
      [omega_damped, Q] = self.omega_damped_and_Q();
      omega_damped_hz =  omega_damped/(2*pi);
    end
    
    % Calculate the quality factor for a given cantilever design assuming atmospheric pressure in air
    % Implemented per "Dependence of the quality factor of micromachined silicon beam resonators on pressure and geometry" by Blom (1992)
    % Is faster to calculate than Sader's
    function Q = calculateBlomQ(self)
      if strcmp(self.fluid, 'vacuum')
        Q = cantilever.maxQ;
        return;
      end
      
      [rho_f eta_f] = self.lookupFluidProperties();
      k_0 = 1.875; % first resonant mode factor
      omega_0 = self.omega_vacuum();
      R = sqrt(self.w*self.l/pi); % effective sphere radius
      delta = sqrt(2*eta_f/rho_f/omega_0); % boundary layer thickness
      
      Q = k_0^2/(12*pi*sqrt(3))*sqrt(self.rho_cantilever*self.modulus())*self.w*self.t^2/(self.l*R*(1+R/delta)*eta_f);
      Q = min(Q, cantilever.maxQ);
    end
    
    % Calculate the Reynold's number
    function reynolds = reynolds(self, omega, rho_f, eta_f)
      reynolds = (rho_f*omega*self.w^2)/eta_f;
    end
    
    % Calculate kappa for the first mode
    function kappa = kappa(self)
      % Not used currently (for speed reasons), but kappa can be calculated for higher order modes
      %       C_n(1) = fzero(@(x) 1+cos(x)*cosh(x),1);
      %       C_n(2) = fzero(@(x) 1+cos(x)*cosh(x),5);
      %       C_n(3) = fzero(@(x) 1+cos(x)*cosh(x),7);
      %       C_n(4) = fzero(@(x) 1+cos(x)*cosh(x),10);
      %       C = C_n(1); % Right now only for 1st order modes
      C = 1.8751;
      kappa = C * self.w / self.l;
    end
    
    % Calculate the hydrodynamic function for the inviscid case (Re >> 1)
    % Not generally useful, but included for completeness.
    % From "Resonant frequencies of a rectangular cantilever beam immersed in a fluid", Sader (2006)
    function hydro = hydrodynamic_function_inviscid(self)
      kappa = self.kappa();
      hydro = (1 + 0.74273*kappa + 0.14862*kappa^2)/(1 + 0.74273*kappa + 0.35004*kappa^2 + 0.058364*kappa^3);
    end
    
    % Calculate hydrodynamic function from the lookup table provided in Eysden and Sader (2007)
    function hydro = hydrodynamic_function(self, omega, rho_f, eta_f)
      
      % Calculate Re and kappa
      kappa = self.kappa();
      reynolds = self.reynolds(omega, rho_f, eta_f);
      log_reynolds = log10(reynolds);
      
      % Lookup the tau components
      tau_real = interp2(cantilever.kappa_lookup, cantilever.reynolds_lookup, cantilever.tau_lookup_real, kappa, log_reynolds, 'linear');
      tau_imag = interp2(cantilever.kappa_lookup, cantilever.reynolds_lookup, cantilever.tau_lookup_imag, kappa, log_reynolds, 'linear');
      
      hydro = complex(tau_real, tau_imag);
    end
    
    % ==================================
    % ======= Simulate response  =======
    % ==================================
    
    
    function [t, voltageNoise] = calculateSimulinkNoise(self, tMax, Fs)
      t = 0:1/Fs:tMax;

      % Calculate the voltage noise parameters
      whiteNoiseSigma = self.voltage_noise(self.freq_max);
      fCorner = self.knee_frequency();
      overSampleRatio = Fs/(2*self.freq_max);

      % Generate voltage noise that matches the calculated spectrum
      bandwidth = self.freq_max - self.freq_min;
      whiteNoise = sqrt(bandwidth)*whiteNoiseSigma*randn(size(t))*sqrt(overSampleRatio);
      pinkNoise = 2*pi*cumsum(whiteNoise)*fCorner*sqrt(whiteNoiseSigma)/overSampleRatio;
      totalNoise = whiteNoise + pinkNoise;
      voltageNoise = [t' totalNoise'];
    end
    
    function [tSim, inputForce, actualForce, sensorForce] = simulateForceStep(self, tMax, Fs, forceMagnitude, forceDelay, forceHold)
      [time, inputNoise] = self.calculateSimulinkNoise(tMax, Fs);
      
      % Generate the force
      forceSignal = zeros(length(time), 1);
      forceSignal(find(time>forceDelay,1):find(time>forceDelay+forceHold,1)) = forceMagnitude;
      inputForce = [time' forceSignal];

      % Define the parameters
      [omega_damped Q] = self.omega_damped_and_Q();
      SFpr = self.force_sensitivity();
      k = self.stiffness();      
      m = k/omega_damped^2;
      b = k/omega_damped/Q;

      % Amplifier
      Kamp = 1e3;
      Famp = 800e3; % Hz
      Tamp = 1/(2*pi*Famp);
      Sfv = SFpr*Kamp;
      
      % Filters
      Tlowpass = 1/(2*pi*self.freq_max);
      Thighpass = 1/(2*pi*self.freq_min);
      
      options = simset('FixedStep', 1/Fs, 'SrcWorkspace', 'current');
      tSim = sim('sensorSimulation', tMax, options);
      
    end

    function simulateAndPlotForceStep(self, tMax, Fs, forceMagnitude, forceDelay, forceHold)
      [tSim, inputForce, actualForce, sensorForce] = self.simulateForceStep(tMax, Fs, forceMagnitude, forceDelay, forceHold);
      
      figure
      hold all
      plot(tSim*1e6, 1e12*inputForce(:,2)); % work from inputVoltage rather than inputForce due to diff in time
      plot(tSim*1e6, 1e12*actualForce);
      plot(tSim*1e6, 1e12*sensorForce);
      hold off
      xlabel('Time (microseconds)');
      ylabel('Output (pN)');
      legend('Applied Force', 'Cantilever Response', 'Sensor Response', 'Location', 'Best')      
    end
    
    function simulateAndPlotMultipleForceSteps(self, tMax, Fs, forceMagnitude, forceDelay, forceHold, numSims)

      rms_mdf = self.force_resolution();
      pp_mdf = 3*rms_mdf;
      
      figure
      hold all

      % Do an initial simulation to find the expected trajectory
      [tSim, inputForce, actualForce, sensorForce] = self.simulateForceStep(tMax, Fs, forceMagnitude, forceDelay, forceHold);
      
      % Fill in the grey background and draw the black nominal force line
      x = 1e6*[tSim ; flipud(tSim)];
      y = 1e12*[actualForce+pp_mdf ; flipud(actualForce-pp_mdf)];
      patch(x, y, cantilever.colorLightGrey, 'EdgeColor', 'none')

      x = 1e6*[tSim ; flipud(tSim)];
      y = 1e12*[actualForce+rms_mdf ; flipud(actualForce-rms_mdf)];
      patch(x, y, cantilever.colorDarkGrey, 'EdgeColor', 'none')      
      
      % Do the simulations and plot
      for ii = 1:numSims
        [tSim, inputForce, actualForce, sensorForce] = self.simulateForceStep(tMax, Fs, forceMagnitude, forceDelay, forceHold);
        plot(1e6*tSim, 1e12*sensorForce, 'Color', cantilever.colorBlue)
      end
      plot(1e6*tSim, 1e12*actualForce, '-', 'Color', cantilever.colorBlack, 'LineWidth', cantilever.plotLineWidth+1);      
      hold off
      xlabel('Time (microseconds)');
      ylabel('Output (pN)');
    end
    
    function simulateAndPlotMultipleForceNoise(self, tMax, Fs, numSims)
      self.simulateAndPlotMultipleForceSteps(tMax, Fs, 0, 0, 0, numSims);
    end
    
    
    % ==================================
    % ========= Optimization  ==========
    % ==================================
    
    % Calculate force resolution from the cantilever state variable vector
    % Units: pN
    function force_resolution = optimize_force_resolution(self, x0)
      self = self.cantilever_from_state(x0);
      force_resolution = self.force_resolution()*1e12;
    end
    
    % Calculate displacement resolution from the cantilever state variable vector
    % Units: nm
    function displacement_resolution = optimize_displacement_resolution(self, x0)
      self = self.cantilever_from_state(x0);
      displacement_resolution = self.displacement_resolution()*1e9;
    end
    
    
    % Nonlinear optimization constraints. For a feasible design, all constraints are negative.
    function [C, Ceq] = optimization_constraints(self, x0, nonlinear_constraints)
      
      c_new = self.cantilever_from_state(x0);
      
      % Read out the constraints as key-value pairs, e.g. {{'omega_min_hz', 'min_k'}, {1000, 10}}
      if ~isempty(nonlinear_constraints)
        keys = nonlinear_constraints{1};
        values = nonlinear_constraints{2};
        for ii = 1:length(keys)
          eval([keys{ii} '=' num2str(values{ii}) ';']);
        end
      end
      
      % Force resolution must always be positive
      % We start with this single element vector and then append any additional constraints that the user has provided.
      C(1) = -c_new.force_resolution();
      
      % Resonant frequency
      if exist('omega_min_hz', 'var')
        switch self.fluid
          case 'vacuum'
            freq_constraint = omega_min_hz - c_new.omega_vacuum_hz();
          otherwise
            [omega_damped_hz, Q] = c_new.omega_damped_hz_and_Q();
            freq_constraint = omega_min_hz - omega_damped_hz;
        end
        C = [C freq_constraint];
      end
      
      % Power dissipation
      if exist('max_power', 'var')
        power_constraint = c_new.power_dissipation() - max_power;
        C = [C power_constraint];
      end
      
      % Temp constraints
      if exist('tip_temp', 'var')
        [TMax TTip] = c_new.approxTempRise();
        temp_constraint = TTip - tip_temp;
        C = [C temp_constraint];
      end
      
      if exist('max_temp', 'var')
        [TMax TTip] = c_new.approxTempRise();
        temp_constraint = TMax - max_temp;
        C = [C temp_constraint];
      end
      
      % Min and maximum cantilever stiffness
      if exist('min_k', 'var')
        min_k_constraint = min_k - c_new.stiffness();
        C = [C min_k_constraint];
      end
      
      if exist('max_k', 'var')
        max_k_constraint = c_new.stiffness() - max_k;
        C = [C max_k_constraint];
      end
      
      
      % Cantilevera and PR aspect ratios
      if exist('min_l_w_ratio', 'var')
        length_width_ratio = min_l_w_ratio - c_new.l/c_new.w;
        C = [C length_width_ratio];
      end
      
      if exist('min_w_t_ratio', 'var')
        width_thickness_ratio = min_w_t_ratio - c_new.w/c_new.t;
        C = [C width_thickness_ratio];
      end
      
      if exist('min_pr_l_w_ratio', 'var')
        pr_length_width_ratio = min_pr_l_w_ratio - c_new.l_pr()/c_new.w_pr();
        C = [C pr_length_width_ratio];
      end
      
      
      % Now for equality based constraints
      Ceq = [];
      
      % Fix the stiffness
      if exist('fixed_k', 'var')
        fixed_k_constraint = c_new.stiffness() - fixed_k;
        Ceq = [Ceq fixed_k_constraint];
      end
      
      % Fix the resonant frequency
      if exist('fixed_f0', 'var')
        switch self.fluid
          case 'vacuum'
            fixed_f0_constraint = omega_min_hz - c_new.omega_vacuum_hz();
          otherwise
            [omega_damped_hz, Q] = c_new.omega_damped_hz_and_Q();
            fixed_f0_constraint = omega_min_hz - omega_damped_hz;
        end
        Ceq = [Ceq fixed_f0_constraint];
      end
    end
    
    % The optimization isn't convex so isn't guaranteed to converge. In practice it converges about 95% of the time
    % depending on the initial guess and constraint set. For this reason, it is best to start from a random initial
    % seed and perform the optimization and checking to make sure that it converges repeatedly.
    function optimized_cantilever = optimize_performance(self, parameter_constraints, nonlinear_constraints, goal)
      
      percent_match = 0.001; % 0.1 percent
      randomize_starting_conditions = 1;
      
      converged = 0;
      ii = 1;
      resolution = [];
      while ~converged
        % Optimize another cantilever
        [c{ii}, exitflag] = self.optimize_performance_once(parameter_constraints, nonlinear_constraints, goal, randomize_starting_conditions);
        
        % If the optimization terminated abnormally (e.g. constraints not satisfied), skip to the next iteration
        if ~(exitflag == 1 || exitflag == 2)
          c{ii}
          continue
        end
        
        % Record the resolution for the latest cantilever
        if goal == cantilever.goalForceResolution
          resolution(ii) = c{ii}.force_resolution();
        elseif goal == cantilever.goalDisplacementResolution
          resolution(ii) = c{ii}.displacement_resolution();
        end
        
        % If we have more than one result, consider stopping
        if ii > 1
          % Sort from smallest to largest, check if the two smallest values agree
          [resolution, sortIndex] = sort(resolution);
          fprintf('Resolutions so far: %s\n', mat2str(resolution, 3))
          resultsAgree = abs(1 - resolution(1)/resolution(2)) < percent_match;
          
          % If the results agree, then stop the loop. Otherwise, continue
          if resultsAgree
            fprintf('CONVERGED. Two best values: %s\n', mat2str(resolution(1:2), 3))
            optimized_cantilever = c{sortIndex(1)};
            converged = 1;
          else
            fprintf('NOT CONVERGED. Two best values: %s\n', mat2str(resolution(1:2), 3))
          end
        end
        
        % After a few tries, we'll just use the best result we came across
        if ii > 10
          [resolution, sortIndex] = sort(resolution);
          optimized_cantilever = c{sortIndex(1)};
          converged = 1;
        end
        
        % Increment
        ii = ii + 1;
      end
    end
    
    % Optimize, but don't randomize starting point
    function optimized_cantilever = optimize_performance_from_current(self, parameter_constraints, nonlinear_constraints, goal)
      randomize_starting_conditions = 0;
      [optimized_cantilever, exitflag] = self.optimize_performance_once(parameter_constraints, nonlinear_constraints, goal, randomize_starting_conditions);
    end
    
    function [optimized_cantilever, exitflag] = optimize_performance_once(self, parameter_constraints, nonlinear_constraints, goal, randomize_starting_conditions)
      
      scaling = self.optimization_scaling();
      
      % If random_flag = 1, start from random conditions. Otherwise
      % start from the current cantilever state vector
      if randomize_starting_conditions == 1
        problem.x0 = scaling.*self.initial_conditions_random(parameter_constraints);
      else
        problem.x0 = scaling.*self.current_state();
      end
      
%       % To ensure that we don't have a bad initial point
%       check_lever = self.cantilever_from_state(problem.x0);
%       check_lever.print_performance()
      
      if goal == cantilever.goalForceResolution
        problem.objective = @self.optimize_force_resolution;
      elseif goal == cantilever.goalDisplacementResolution
        problem.objective = @self.optimize_displacement_resolution;
      end
      
      [lb ub] = self.optimization_bounds(parameter_constraints);
      problem.lb = scaling.*lb;
      problem.ub = scaling.*ub;
      
      problem.options.TolFun = 1e-12;
      problem.options.TolCon = 1e-12;
      problem.options.TolX = 1e-12;
      
      problem.options.MaxFunEvals = 3000;
      problem.options.MaxIter = 2000;
      problem.options.Display = 'iter';
      problem.options.UseParallel = 'always'; % For multicore processors
      
      problem.options.Algorithm = 'Interior-point';
      problem.solver = 'fmincon';
      
      problem.nonlcon = @(x) self.optimization_constraints(x, nonlinear_constraints);
      
      [x, fval, exitflag] = fmincon(problem);
      optimized_cantilever = self.cantilever_from_state(x);
    end
  end
end
