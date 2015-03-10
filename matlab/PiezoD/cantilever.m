% Abstract cantilever base class
% Used by subclasses to implement common cantilever features
% The code is formatted to be < 85 columns wide for printing reproduction
%
% Subclasses must implement the following methods at a minimum:
%  doping_profile
%  doping_optimization_scaling
%  doping_cantilever_from_state
%  doping_current_state
%  doping_initial_conditions_random
%  doping_optimization_bounds
classdef cantilever
    properties
        % Basic pararameters that are absolutely required
        freq_min; % Minimum measured frequency, Hz
        freq_max; % Maximum measured frequency, Hz
        l; % overall cantilever/sensor length, m
        w; % overall cantilever width, m
        t; % overall cantilever thickness, m
        l_pr_ratio; % piezoresistor length ratio
        v_bridge; % bridge bias, V
        doping_type; % Dopant species: 'phosphorus', 'boron', or 'arsenic'
        
        % Optional parameters for general operation
        fluid; % 'vacuum', 'air', 'water', or 'arbitrary'. Default = 'air'
        number_of_piezoresistors; % # of PRs in the circuit. Default = 2.
        rms_actuator_displacement_noise = 1e-12; % Displacement noise from mounting, m
        default_alpha = 1e-5; % Hooge noise parameter, unitless
        amplifier; % Amplifier used: 'INA103' or 'AD8221'. Default = 'INA103'.
        T; % Ambient temp, K. Default = 300.
        T_ref; % Reference temp for calculating thermal expansion, K. Default = 300
        tip_mass; % Tip mass loading (e.g. bead or sharp tip). Default = 0.
        R_base; % Thermal resistance at the cantilever base, K/W. Default = 10e3
        R_contact; % Excess resistance per contact, default = 170 ohm
        
        % Include temp dependence in calculating electrical resistivity
        % and thermal conductivity of silicon. 'yes' or 'no'. Default = 'no'
        temperature_dependent_properties;
        
        % Degree of accuracy to use in thermal modeling (e.g. see johnson_PSD())
        % Options: 'none', 'approx' or 'exact'. Default = 'none'
        thermal_modeling;
        
        % Optional parameters for modeling arbitrary fluids (damping and heat transfer)
        rho_arb; % Density of fluid
        eta_arb; % Viscosity of fluid
        k_arb; % Thermal conductivity of fluid
        h_method; % How to calculate h_eff: 'fixed' or 'calculate'
        
        % Effective convection coefficient of fluid, W/m^2-K
        % Used in case of h_method = 'fixed' and fluid = 'arbitrary'
        h_arb;
        
        % Optional parameters for modeling a stiffener/actuator at the
        % base of the device. Mechanical/thermal impact is included in the code.
        
        % Type of the metal used in the actuator/step
        % Options: 'aluminum', 'titanium' or 'molybdenum'
        % Default = 'aluminum'
        metal_type;
        
        % Options = 'none', 'step', 'thermal', 'piezoelectric'.
        % Default = 'None'
        % 'step' = Si/SiO2/metal (t/t_oxide/t_a)
        % 'thermal' = Si/SiO2/metal (t/t_oxide/t_a)
        % 'piezoelectric' = Si/SiO2/AlN/metal/AlN/metal
        % (t/t_oxide/t_a_seed/t_electrode_bottom/t_a/t_electrode_top)
        cantilever_type;
        
        l_a; % Overall length of the step/actuator, m
        w_a; % Overall width of the step/actuator, m
        w_a_active; % Active width of the step/actuator, m
        t_a; % Thickness of the step/actuator, m
        v_actuator; % Bias voltage for either thermal or PE actuation, V
        R_heater; % Resistance of the thermal actuator heater, ohm
        t_oxide; % Oxide layer thickness, m
        t_a_seed; % Thickness of PE seed layer, m
        t_electrode_bottom; % Thickness of the bottom PE metal electrode, m
        t_electrode_top; % Thickness of the top PE metal electrode, m
        d31_manual; % Allows manual specification of the PE d31 coefficient, pm/V
        
        % Gap at the end of the actuator before the start of the sensor, m
        % Typically on the order of 20 microns to allow for contacts, etc.
        l_a_gap;
        
        % Probe deflection due to residual stress is modeled in the code
        % Options for modeling: 'nominal' or 'random'
        % The latter option performs a Monte Carlo simulation to predict
        % the tip range from the input stress range (e.g. see sigma_si_range)
        film_stress;
    end
    
    % Constants - can be referred to with cantilever.variableName
    properties (Constant)
        
        % Physical constants
        k_b = 1.38e-23; % J/K
        k_b_eV = 8.617343e-5; % eV/K
        q = 1.60218e-19; % Coulombs
        h_bar = 1.055e-34; % J-sec
        
        % Define the number of points to use in discretized calculations
        numFrequencyPoints = 1000; % For noise spectra plotting
        
        % Number of points along the cantilever for calculating
        % deflections and temperature profiles
        numXPoints = 800;
        
        % Number of points through the cantilever depth for calculating the
        % dopant profile, sheet resistance, number of carriers, thermal conductivity
        numZPoints = 200;
        
        % Number of Monte Carlo iterations for calculating the tip
        % deflection distribution
        numRandomStressIterations = 10;
        
        % Number of optimization iterations to perform without
        % convergence (to within 0.1%) before accepting
        % the best answer so far
        numOptimizationIterations = 20;
        
        
        % The width between the two cantilever legs. Afffects the PR
        % resistance and sensitivity, but only if l_pr is small
        air_gap_width = 2e-6;
        
        % Standard fluid properties
        k_water = 0.610; % W/m-K
        rho_water = 996.6; % kg/m^3
        eta_water = 7.98e-4; % Pa-sec
        h_water = 49218; % W/m^2-k
        
        k_air = 0.0262; % W/m-K
        rho_air = 1.164; % kg/m^3
        eta_air = 17e-6; % Pa-sec
        h_air = 2098; % W/m^2-K
        
        % For vacuum, use small but finite values for numerical stability
        k_vacuum = 1e-6; % W/m-K
        rho_vacuum = 1e-6; % kg/m^3
        eta_vacuum = 1e-6; % Pa-sec
        h_vacuum = 1e-6; % W/m^2-K
        
        % Thermal conductivities, W/m-K
        % AlN: "Process-dependent thin-film thermal conductivities for thermal CMOS MEMS"
        % Al: "MEMS test structure for measuring thermal conductivity of thin films"
        k_si = 148;
        k_al = 200;
        k_sio2 = 1.4;
        k_ti = 21.9;
        k_aln = 60;
        k_mo = 138;
        
        % Coefficients of thermal expansion, 1/K
        alpha_al = 23.1e-6;
        alpha_sio2 = 0.5e-6;
        alpha_si = 2.6e-6;
        alpha_ti = 8.6e-6;
        alpha_aln = 4.5e-6;
        alpha_mo = 4.8e-6;
        
        % Silicon temperature coefficient of resistance (TCR), 1/K
        % Corresponds to a peak doping of about 1e20/cc. Used to predict
        % resistance change with self-heating in the advanced thermal models
        TCR = 1372e-6;
        
        % Intrinsic stress (Pa)
        % For modeling tip deflection from film stress
        % Doping stress is treated separately elsewhere in the code
        % film_stress == 'nominal' uses the average
        % film_stress == 'random' uses a normal distribution
        sigma_si_range = 1e6*[0 0];
        sigma_sio2_range = 1e6*[-200 -300];
        sigma_al_range = 1e6*[160 200];
        sigma_aln_range = 1e6*[-300 300];
        sigma_ti_range = 1e6*[-25 25];
        sigma_mo_range = 1e6*[-25 25];
        
        % Specific heat (J/kg-K) for calculating thermal time constants
        Cv_si = 700;
        Cv_al = 910;
        Cv_sio2 = 700;
        
        % Mechanical material properties (Pa)
        nu_Si = 0.28; % <100> direction
        nu_Al = 0.3;
        nu_Ti = 0.3;
        nu_AlN = 0.24;
        nu_SiO2 = 0.17;
        nu_Mo = 0.3;
        
        % Assume plane strain conditions
        % (i.e. cantilever much wider than it is thick)
        E_si = 130e9/(1 - cantilever.nu_Si^2);
        E_al = 70e9/(1 - cantilever.nu_Al^2);
        E_ti = 90e9/(1 - cantilever.nu_Ti^2);
        E_aln = 320e9/(1 - cantilever.nu_AlN^2); % or 345 GPa
        E_sio2 = 75e9/(1 - cantilever.nu_SiO2^2);
        E_mo = 329e9/(1 - cantilever.nu_Mo^2); % or 270 GPa
        
        % Densities (kg/m^3)
        rho_si = 2330;
        rho_al = 2700;
        rho_ti = 4506;
        rho_aln = 3260;
        rho_sio2 = 2200;
        rho_mo = 10280;
        
        % Transverse piezoelectric coefficient of AlN, pm/V or pC/N
        % d31 varies with thickness, so interpolate from literature values
        % Values are from papers written by the Piazza and Roukes groups
        d31_t = 1e-9*[50 100 500 3000]; % nm
        d31_aln = 1e-12*[1.9 2.3 2.5 2.6]; % pm/V
        
        % Minimum and maximum quality factors for the cantilever
        % Q is calculated from fluid damping. Results at atmospheric pressure
        %	are accurate, but other damping mechanisms that become important
        % in vacuum (e.g. thermoelastic damping) are not included
        maxQ = 5e3;
        minQ = 1e-6;
        
        % Define the possible optimization goals.
        % Can be referred to as cantilever.goalForceResolution, etc.
        goalForceResolution = 0;
        goalDisplacementResolution = 1;
        goalForceNoiseDensity = 2;
        goalSurfaceStress = 3;
        
        % Lookup table for calculating resonant frequency and quality factor in liquid
        % Source: "Oscillations of cylinders..." by Brumley, Wilcox and Sader (2010)
        % Values are stored as constants for calculation speed
        A_lookup = [0 1/50 1/20 1/10 1/5 1/2 1 2 5 10 20 50 1000]; % t/w ratio
        Beta_lookup = [-3 -2.5 -2 -1.5 -1 -.5 0 0.5 1 1.5 2 2.5 3 100]; % log(Re)
        
        % Hydrodynamic function
        % Includes fixed bottom row from published erratum
        gamma_lookup_real = ...
            [212.184 213.310 214.977 217.701 222.978 237.780 260.256 ...
            207.210  169.667   154.616   145.909    139.855    134.720;
            91.6984 92.2467 93.0601 94.3924 96.9808 104.295 115.542 ...
            88.9011  70.8173   63.7655   59.7404    56.9653    54.6258;
            41.6417 41.9209 42.3363 43.0185 44.3487 48.1380 54.0391 ...
            39.8564  30.6996   27.2460   25.3060    23.9817    22.8730;
            20.1196 20.2683 20.4907 20.8572 21.5753 23.6370 26.8847 ...
            18.8235  13.9212   12.1457   11.1673    10.5072    9.95883;
            10.4849 10.5677 10.6926 10.8998 11.3080 12.4883 14.3601 ...
            9.43536  6.64606   5.68511   5.16801    4.82411    4.54093;
            5.96655 6.01467 6.08871 6.21279 6.45897 7.17328 8.30052 ...
            5.04739  3.35215   2.80394   2.51794    2.33126    2.17927;
            3.73387 3.76344 3.81063 3.89099 4.05154 4.51368 5.22220 ...
            2.89030  1.78322   1.45306   1.28807    1.18327    1.09943;
            2.56548 2.58563 2.61959 2.67832 2.79515 3.11907 3.58531 ...
            1.77617  0.994540  0.783333  0.684003   0.623512   0.576619;
            1.91834 1.93509 1.96437 2.01450 2.11058 2.35665 2.68270 ...
            1.17779  0.580514  0.435349  0.372208   0.336075   0.309503;
            1.54554 1.56285 1.59247 1.64069 1.72687 1.92785 2.17551 ...
            0.848104 0.357549  0.249659  0.206674   0.184001   0.168601;
            1.32633 1.34658 1.37882 1.42757 1.50844 1.68437 1.88862 ...
            0.663505 0.235193  0.148772  0.117201   0.102069   0.0928779;
            1.19577 1.2202  1.2555  1.3051  1.3833  1.5459  1.7259  ...
            0.55939  0.16703   0.093131  0.068128   0.057273   0.0515648;
            1.11746 1.1465  1.1843  1.2346  1.3117  1.4670  1.6336  ...
            0.50051  0.12874   0.062098  0.040918   0.032516   0.0287745;
            1       1.04551 1.08816 1.14064 1.21703 1.36368 1.51317 ...
            0.423881 0.0792129 0.0222121 0.00619303 0.00113212 0        ];
        
        gamma_lookup_imag = ...
            [1018.72  1021.37  1025.29  1031.66  1043.88  1077.39  1126.32 ...
            1008.65  915.159  874.583  850.149  832.704  817.599;
            374.276  375.392  377.040  379.721  384.873  399.079  420.012  ...
            370.057  331.318  314.778  304.899  297.884  291.835;
            140.659  141.144  141.862  143.031  145.284  151.534  160.848  ...
            138.825  122.228  115.278  111.167  108.266  105.776;
            54.4049  54.6253  54.9508  55.4818  56.5079  59.3754  63.7087  ...
            53.5749  46.1812  43.1534  41.3825  40.1420  39.0836;
            21.8269  21.9314  22.0855  22.3371  22.8247  24.2002  26.3169  ...
            21.4324  17.9905  16.6153  15.8210  15.2692  14.8012;
            9.16870  9.22024  9.29587  9.41936  9.65973  10.3480  11.4345  ...
            8.96804  7.28929  6.63516  6.26219  6.00523  5.78862;
            4.07467  4.10043  4.13779  4.19895  4.31957  4.67605  5.25977  ...
            3.95920  3.10274  2.77671  2.59298  2.46733  2.36186;
            1.93366  1.94552  1.96256  1.99130  2.05107  2.24127  2.56535  ...
            1.85252  1.39790  1.22868  1.13429  1.07013  1.01639;
            0.981710 0.985312 0.990956 1.00255  1.03157  1.13634  1.31768  ...
            0.915797 0.666095 0.575374 0.525354 0.491568 0.463359;
            0.527773 0.526433 0.526077 0.529479 0.543868 0.602276 0.703142 ...
            0.474037 0.333253 0.283225 0.256021 0.237799 0.222666;
            0.296143 0.291987 0.289093 0.289338 0.296683 0.328687 0.384789 ...
            0.253907 0.173548 0.145302 0.130165 0.120135 0.111868;
            0.171115 0.16564  0.16234  0.16171  0.16525  0.18260  0.21384  ...
            0.13910  0.093151 0.076988 0.068405 0.062790 0.0582134;
            0.100688 0.095021 0.092307 0.091476 0.093044 0.10247  0.11987  ...
            0.077266 0.051022 0.041760 0.036840 0.033652 0.0310905;
            0        0        0        0        0        0        0        ...
            0        0        0        0        0        0];
    end
    
    % Abstract methods that MUST be defined in subclasses
    methods (Abstract)
        doping_profile(self)
        doping_optimization_scaling(self)
        doping_cantilever_from_state(self)
        doping_current_state(self)
        doping_initial_conditions_random(self)
        doping_optimization_bounds(self, parameter_constraints)
        
        % Abstract functions for handling ion implantation
        Nz(self)
        alpha(self)
        sheet_resistance(self)
    end
    
    methods
        function self = cantilever(freq_min, freq_max, l, w, t, ...
                l_pr_ratio, v_bridge, doping_type)
            % Define the cantilever using the required parameters
            self.freq_min = freq_min;
            self.freq_max = freq_max;
            self.l = l;
            self.w = w;
            self.t = t;
            self.l_pr_ratio = l_pr_ratio;
            self.v_bridge = v_bridge;
            self.doping_type = doping_type;
            
            % Set default values for optional parameters
            % These can be overridden once the cantilever is created
            self.fluid = 'air';
            self.h_method = 'fixed';
            self.metal_type = 'aluminum';
            self.film_stress = 'nominal';
            self.cantilever_type = 'none';
            self.l_a = 0;
            self.t_a = 0;
            self.w_a = 0;
            self.w_a_active = 0;
            self.d31_manual = 0;
            self.l_a_gap = 0;
            self.tip_mass = 0;
            self.t_oxide = 100e-9;
            self.t_electrode_bottom = 50e-9;
            self.t_electrode_top = 50e-9;
            self.t_a_seed = 20e-9;
            self.number_of_piezoresistors = 2;
            self.amplifier = 'INA103';
            self.T = 273.15 + 23; % Assume operation around 23C
            self.T_ref = 273.15 + 23;
            self.R_base = 10e3; % Assume R_base is about 10e3 K/W by default
            self.R_contact = 170; % Assume 340 ohm total excess resistance
            self.temperature_dependent_properties = 'no';
            self.thermal_modeling = 'none';
        end
        
        % Helpful getter function to calculate the absolute piezoresistor length
        % Units: m
        function l_pr = l_pr(self)
            l_pr = self.l * self.l_pr_ratio;
        end
        
        % Getter function for piezoresistor width
        % Units: m
        function w_pr = w_pr(self)
            w_pr = self.w/2;
        end
        
        % Determine the ion implantation table index from the dopant type
        function dopantNumber = dopantNumber(self)
            switch self.doping_type
                case 'boron'
                    dopantNumber = 1;
                case 'phosphorus'
                    dopantNumber = 2;
                case 'arsenic'
                    dopantNumber = 3;
                otherwise
                    fprintf('ERROR: Unknown dopant type: %s\n', self.doping_type);
                    pause
            end
        end
        
        % Check if the cantilever is self-consistent
        function check_valid_cantilever(self)
            validity_checks(1) = 1;
            validity_checks(2) = ~(strcmp(self.cantilever_type, 'none') ...
                && (self.l_a > 0));
            [valid, failed_index] = min(validity_checks);
            if ~valid
                fprintf('ERROR: Invalid cantilever - failed check #%d\n', failed_index);
                pause
            end
        end
        
        % Print the cantilever performance to the Matlab prompt
        function print_performance(self)
            self.check_valid_cantilever();
            
            [omega_damped_hz, Q] = self.omega_damped_hz_and_Q();
            [x, active_doping, total_doping] = self.doping_profile();
            [TMax_approx TTip_approx] = self.approxTempRise();
            [TMax, TTip] = self.calculateMaxAndTipTemp();
            thermoLimit = self.thermo_integrated()/self.force_sensitivity();
            
            fprintf('=======================\n')
            fprintf('Freq range: %f to %f \n', self.freq_min, self.freq_max)
            fprintf('Operating fluid: %s \n', self.fluid);
            fprintf('Cantilever L/W/T: %f %f %f \n', ...
                self.l*1e6, self.w*1e6, self.t*1e6)
            fprintf('PR L/W: %f %f \n', self.l_pr()*1e6, self.w_pr()*1e6)
            fprintf('PR Length Ratio: %g \n', self.l_pr_ratio)
            fprintf('\n')
            fprintf('Force resolution (N): %g \n', self.force_resolution())
            fprintf('Force noise at 1 kHz (fN): %g \n', ...
                self.force_noise_density(1e3))
            fprintf('Displacement resolution (m): %g \n', ...
                self.displacement_resolution())
            fprintf('Force sensitivity (V/N) %g \n', self.force_sensitivity())
            fprintf('Displacement sensitivity (V/m) %g \n', ...
                self.displacement_sensitivity())
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
            fprintf('Approx. Temp Rises (C) - Tip: %f  Max: %f\n', ...
                TTip_approx, TMax_approx)
            fprintf('F-D Temp Rises (C)     - Tip: %f  Max: %f\n', ...
                TTip, TMax)
            fprintf('\n')
            fprintf('Integrated noise (V): %g \n', self.integrated_noise())
            fprintf('Integrated johnson noise (V): %g \n', self.johnson_integrated())
            fprintf('Integrated 1/f noise (V): %g \n', self.hooge_integrated())
            fprintf('Amplifier noise (V): %g \n', self.amplifier_integrated())
            fprintf('Thermomechanical noise (V): %g \n', self.thermo_integrated())
            fprintf('\n')
            fprintf('Johnson/Hooge: %g \n', ...
                self.johnson_integrated()/self.hooge_integrated())
            fprintf('Knee frequency (Hz): %g \n', self.knee_frequency())
            fprintf('Number of Carriers: %g \n', self.number_of_carriers());
            fprintf('Nz: %g \n', self.Nz())
            fprintf('\n')
            fprintf('Number of silicon resistors: %f \n', ...
                self.number_of_piezoresistors)
            fprintf('Si Thermal Conductivity (W/m-K): %f \n', self.k_base())
            fprintf('E (GPa): %f \n', self.modulus()*1e-9)
            fprintf('Alpha: %g \n', self.alpha)
            
            switch self.cantilever_type
                case 'none'
                    % Do nothing special
                case 'step'
                    fprintf('=======================\n')
                    fprintf('Step at base (um): %f thick x %f long \n', ...
                        1e6*self.t_a, 1e6*self.l_a)
                case 'thermal'
                    fprintf('=======================\n')
                    [tau, freq] = self.heaterTimeConstant();
                    fprintf('Actuator l/W/T: %f %f %f \n', ...
                        1e6*self.l_a, 1e6*self.w_a, 1e6*self.t_a)
                    fprintf('Neutral axis (um): %f \n', 1e6*self.actuatorNeutralAxis())
                    fprintf('Actuator Voltage (): %f \n', self.v_actuator)
                    fprintf('Heater resistance (kOhm): %f \n', 1e-3*self.R_heater)
                    fprintf('Actuator Power (mW): %f \n', 1e3*self.heaterPower())
                    fprintf('Tip Deflection (nm): %f \n', 1e9*self.tipDeflection())
                    fprintf('Time Constant (microseconds): %f \n', tau*1e6)
                    fprintf('-3dB frequency (kHz): %f \n', freq*1e-3)
                case 'piezoelectric'
                    fprintf('Actuator l/W/T: %f %f %f \n', ...
                        1e6*self.l_a, 1e6*self.w_a, 1e6*self.t_a)
                    fprintf('Neutral axis (um): %f \n', 1e6*self.actuatorNeutralAxis())
                    fprintf('Actuator Voltage (): %f \n', self.v_actuator)
                    fprintf('Tip Deflection (nm): %f \n', 1e9*self.tipDeflection())
            end
            fprintf('=======================\n')
        end
        
        % Tab delimited output (for Excel, etc). Optionally provide a
        % file id (fid) to print the output there. print_performance_for_excel
        % is overridden in subclasses to provide detailed doping output
        function print_performance_for_excel(self, varargin)
            varargin = varargin{1};
            optargin = size(varargin, 2);
            
            % Determine where to print the output
            if optargin == 1
                fid = varargin{1};
            elseif optargin == 0
                fid = 1; % Print to the stdout
            else
                fprintf('ERROR: Extra optional arguments')
            end
            
            % Calculate intermediate quantities
            [omega_damped_hz, Q] = self.omega_damped_hz_and_Q();
            [TMax, TTip] = self.calculateMaxAndTipTemp();
            [TMax_approx TTip_approx] = self.approxTempRise();
            thermoLimit = self.thermo_integrated()/self.force_sensitivity();
            
            variables_to_print = [self.freq_min, self.freq_max*1e-3, ...
                1e6*self.l 1e6*self.w 1e9*self.t 1e6*self.l_pr() ...
                1e6*self.l_a 1e6*self.w_a 1e9*self.t_a, ...
                self.force_resolution()*1e12, thermoLimit*1e12, ...
                self.displacement_resolution()*1e9, self.omega_vacuum_hz()*1e-3, ...
                omega_damped_hz*1e-3, Q, self.stiffness()*1e3, self.v_bridge, ...
                self.resistance()*1e-3, self.sheet_resistance(), ...
                self.power_dissipation()*1e6, TMax, TTip, TMax_approx, TTip_approx, ...
                self.number_of_piezoresistors, self.effective_mass(), self.tip_mass, ...
                self.force_sensitivity(), self.beta(), self.gamma(), ...
                self.integrated_noise()*1e6, self.johnson_integrated()*1e6, ...
                self.hooge_integrated()*1e6, self.amplifier_integrated()*1e6, ...
                self.thermo_integrated()*1e6, self.knee_frequency(), ...
                self.number_of_carriers(), 1e6*self.calculateDopantTipDeflection(), ...
                self.k_x(), self.modulus()*1e-9];
            
            fprintf(fid, '%s \t', self.doping_type);
            fprintf(fid, '%s\t', self.fluid);
            fprintf(fid, '%s\t', self.cantilever_type);
            for print_index = 1:length(variables_to_print)
                fprintf(fid, '%.4g \t', variables_to_print(print_index));
            end
        end
        
        % Calculate total resistance of piezoresistor.
        % Units: ohms
        function resistance = resistance(self)
            resistance = self.number_of_squares()*self.sheet_resistance() + ...
                2*self.R_contact;
        end
        
        % Calculate the number of resistor squares
        % 1) Longitudinal region (2*l_pr/w_pr)
        % 2) Transverse at end (air_gap_width/2*w_pr)
        % 3) the connecting corners
        % Units: squares
        function number_of_squares = number_of_squares(self)
            number_of_squares = 2*self.l_pr()/self.w_pr() + ...
                self.air_gap_width/(2*self.w_pr()) + 2;
        end
        
        % Calculate R-R0/R0 for the cantilever including the TCR
        % Used for calculating h_eff from experimental results
        % Units: -
        function dR = dR_with_temp_rise(self)
            dR = self.TCR*(self.averagePRTemp() - self.T);
        end
        
        function dR = approx_dR_with_temp_rise(self)
            dR = self.TCR*(self.approxPRTemp() - self.T);
        end
        
        % Calculate conductivity for a given dopant concentration.
        % Use vectors rather than looping for speed.
        % Units: C/V-sec-cm
        function sigma = conductivity(self, dopant_concentration)
            switch self.temperature_dependent_properties
                case 'yes'
                    [mu, sigma] = self.mobility(dopant_concentration, self.approxPRTemp());
                case 'no'
                    [mu, sigma] = self.mobility(dopant_concentration, self.T);
            end
        end
        
        % Calculate the temp dependent carrier density, mobility, conductivity
        % Current: "Electron and Hole Mobility ...", Reggiani et al. (2002)
        % Previously: "Modeling of Carrier Mobility ...", Masetti et al. (1983)
        function [mu, sigma] = mobility(self, dopantConc, T)
            Tnorm = T/300;
            Eg = 1.170-(4.730e-4*T.^2)./(T+636);
            ni = (2.4e31.*T.^3.*exp(-Eg./(self.k_b_eV.*T))).^.5;
            
            switch self.doping_type
                case 'boron'
                    mumax=470.5;
                    c=0.0;
                    gamma=2.16;
                    mu0d=90.0*Tnorm^-1.3;
                    mu0a=44.0*Tnorm^-0.7;
                    mu1d=28.2*Tnorm^-2.0;
                    mu1a=28.2*Tnorm^-0.8;
                    Cr1=1.3e18*Tnorm^2.2;
                    Cr2=2.45e17*Tnorm^3.1;
                    Cs1=1.1e18*Tnorm^6.2;
                    Cs2=6.1e20;
                    alpha1=.77;
                    alpha2=.719;
                    ni=(2.4e31*T^3*exp(-1*Eg/8.617e-5/T))^.5;
                    p=(dopantConc./2)+((dopantConc./2).^2+ni^2).^.5;
                    n=ni*ni./p;
                    ND=n;
                    NA=p;
                case 'phosphorus'
                    mumax  = 1441;
                    c      = 0.07;
                    gamma  = 2.45;
                    mu0d   = 62.2.*Tnorm.^-0.7;
                    mu0a   = 132.*Tnorm.^-1.3;
                    mu1d   = 48.6.*Tnorm.^-0.7;
                    mu1a   = 73.5.*Tnorm.^-1.25;
                    Cr1    = 8.5e16.*Tnorm.^3.65;
                    Cr2    = 1.22e17.*Tnorm.^2.65;
                    Cs1    = 4e20;
                    Cs2    = 7e20;
                    alpha1 = .68;
                    alpha2 = .72;
                    n      = (dopantConc./2)+((dopantConc./2).^2 + ni.^2).^.5;
                    p      = ni.^2./n;
                    ND     = n;
                    NA     = p;
                case 'arsenic'
                    mumax  = 1441;
                    c      = 0.07;
                    gamma  = 2.45;
                    mu0d   = 55.*Tnorm.^-0.6;
                    mu0a   = 132.*Tnorm.^-1.3;
                    mu1d   = 42.4.*Tnorm.^-0.5;
                    mu1a   = 73.5.*Tnorm.^-1.25;
                    Cr1    = 8.9e16.*Tnorm.^3.65;
                    Cr2    = 1.22e17.*Tnorm.^2.65;
                    Cs1    = 2.9e20;
                    Cs2    = 7e20;
                    alpha1 = .68;
                    alpha2 = .72;
                    n      = (dopantConc./2)+((dopantConc./2).^2 + ni.^2).^.5;
                    p      = ni.^2./n;
                    ND     = n;
                    NA     = p;
            end
            
            mu0 = (mu0d.*ND+mu0a.*NA)./(ND+NA);
            mu1 = (mu1d.*ND+mu1a.*NA)./(ND+NA);
            muL = mumax.*Tnorm.^(-gamma+c.*Tnorm);
            second = (muL-mu0)./(1+(ND./Cr1).^alpha1 + (NA./Cr2).^alpha2);
            third = mu1./(1+(ND/Cs1+NA/Cs2).^-2);
            
            mu = mu0+second-third;
            sigma = self.q.*(mu.*n + mu.*p);
        end
        
        % Calculate Rsheet(x) assuming temperature dependent cantilever properties
        % Units: ohm/square
        function [Rsheet_x] = RSheetProfile(self, x, T_x)
            [z, active_doping, total_doping] = self.doping_profile();
            % Units: z -> m, doping -> N/cm^3
            n_z_x = total_doping'*ones(1, self.numXPoints);
            
            % Generate a numZPoints x numXPoints matrix
            if length(T_x) == 1
                T_z_x = ones(self.numZPoints, 1)*T_x'*ones(1, self.numXPoints);
            else
                T_z_x = ones(self.numZPoints, 1)*T_x';
            end
            
            [mu_z_x, sigma_z_x] = self.mobility(n_z_x, T_z_x + self.T);
            Rsheet_x = 1./trapz(z*1e2, sigma_z_x); % Convert z from m to cm
        end
        
        % The number of current carriers in the piezoresistor
        % The cantilever subclass implements Nz()
        % Units: -
        function number_of_carriers = number_of_carriers(self)
            resistor_area = self.w_pr()*(2*self.l_pr() + self.w + self.air_gap_width);
            number_of_carriers = self.Nz()*resistor_area;
        end
        
        % 1/f voltage power spectral density for the entire Wheatstone bridge
        % Units: V^2/Hz
        function hooge_PSD = hooge_PSD(self, freq)
            hooge_PSD = self.alpha()*self.v_bridge^2* ...
                self.number_of_piezoresistors./(4*self.number_of_carriers()*freq);
        end
        
        % Integrated 1/f noise density for the entire Wheatstone bridge
        % Unit: V
        function hooge_integrated = hooge_integrated(self)
            hooge_integrated = sqrt(self.alpha()*self.v_bridge^2* ...
                self.number_of_piezoresistors./(4*self.number_of_carriers())* ...
                log(self.freq_max/self.freq_min));
        end
        
        % Johnson noise PSD from the entire Wheatstone bridge
        % Units: V^2/Hz
        function johnson_PSD = johnson_PSD(self, freq)
            resistance = self.resistance(); % resistance() includes contacts
            switch self.thermal_modeling
                case 'none'
                    TPR = self.T; % the ambient temperature
                case 'approx'
                    TPR = self.approxPRTemp();
                case 'exact'
                    TPR = self.averagePRTemp();
            end
            resistance = resistance*(1 + TPR*self.TCR);
            
            if self.number_of_piezoresistors == 4
                R_external = resistance;
            else
                % If we're using 2 piezoresistors, assume ideal 1 kOhm resistors
                R_external = 1e3;
            end
            
            johnson_PSD = 4*self.k_b*TPR*(resistance/2 + ...
                R_external/2) * ones(1, length(freq));
        end
        
        % Integrated Johnson noise
        % Unit: V
        function johnson_integrated = johnson_integrated(self)
            resistance = self.resistance();
            TPR = self.approxPRTemp();
            resistance = resistance*(1 + TPR*self.TCR);
            
            if self.number_of_piezoresistors == 4
                R_external = resistance;
            else
                R_external = 700;
            end
            johnson_integrated = sqrt(4*self.k_b*TPR*(resistance/2 + ...
                R_external/2)*(self.freq_max - self.freq_min));
        end
        
        % Thermomechanical noise PSD
        % Units: V^2/Hz
        function thermo_PSD = thermo_PSD(self, freq)
            [omega_damped_hz, Q_M] = self.omega_damped_hz_and_Q();
            switch self.thermal_modeling
                case 'none'
                    TPR = self.T; % the ambient temperature
                case 'approx'
                    TPR = self.approxPRTemp();
                case 'exact'
                    TPR = self.averagePRTemp();
            end
            thermo_PSD = (self.force_sensitivity())^2*4*self.stiffness()* ...
                self.k_b*TPR/(2*pi*omega_damped_hz*Q_M) * ones(1, length(freq));
        end
        
        % Integrated thermomechanical noise
        % Unit: V
        function thermo_integrated = thermo_integrated(self)
            [omega_damped_hz, Q_M] = self.omega_damped_hz_and_Q();
            switch self.thermal_modeling
                case 'none'
                    TPR = self.T; % the ambient temperature
                case 'approx'
                    TPR = self.approxPRTemp();
                case 'exact'
                    TPR = self.averagePRTemp();
            end
            thermo_integrated = sqrt((self.force_sensitivity())^2* ...
                4*self.stiffness()*self.k_b*TPR/(2*pi*omega_damped_hz*Q_M) ...
                *(self.freq_max - self.freq_min));
        end
        
        % Accounts for displacement noise (vibrations, etc)
        % Units: V
        function actuator_noise_integrated = actuator_noise_integrated(self)
            actuator_noise_integrated = self.rms_actuator_displacement_noise* ...
                self.stiffness()*self.force_sensitivity(); % V
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
            
            R = self.resistance()/2; % Resistance seen at the amp inputs
            switch self.thermal_modeling
                case 'none'
                    TPR = self.T; % the ambient temperature
                case 'approx'
                    TPR = self.approxPRTemp();
                case 'exact'
                    TPR = self.averagePRTemp();
            end
            R = R*(1 + TPR*self.TCR);
            
            amplifier_PSD = (A_VJ^2 + 2*(R*A_IJ)^2) + (A_VF^2 + 2*(R*A_IF)^2)./freq;
        end
        
        % Integrated amplifier noise
        % Units: V
        function amplifier_integrated = amplifier_integrated(self)
            
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
            
            R = self.resistance()/2; % Resistance seen at the amp inputs
            switch self.thermal_modeling
                case 'none'
                    TPR = self.T; % the ambient temperature
                case 'approx'
                    TPR = self.approxPRTemp();
                case 'exact'
                    TPR = self.averagePRTemp();
            end
            R = R*(1 + TPR*self.TCR);
            
            amplifier_integrated = sqrt(A_VJ^2*(self.freq_max - self.freq_min) ...
                + A_VF^2*log(self.freq_max/self.freq_min) + ...
                2*(R*A_IJ)^2*(self.freq_max - self.freq_min) + ...
                2*(R*A_IF)^2*log(self.freq_max/self.freq_min));
        end
        
        % Calculate the knee frequency (Hooge PSD == Johnson PSD)
        % Units: Hz
        function knee_frequency = knee_frequency(self)
            knee_frequency = self.number_of_piezoresistors*self.alpha()* ...
                self.v_bridge^2/(16*self.number_of_carriers()*self.k_b* ...
                self.approxPRTemp()*self.resistance());
        end
        
        % Integrated cantilever noise for given bandwidth
        % Pull the calculations into this function for speed (i.e. don't
        % calculate self.resistance() five separate times
        % Units: V
        function integrated_noise = integrated_noise(self)
            [omega_damped_hz, Q_M] = self.omega_damped_hz_and_Q();
            resistance = self.resistance()/self.gamma();
            force_sensitivity = self.force_sensitivity();
            spring_constant = self.stiffness();
            switch self.thermal_modeling
                case 'none'
                    TPR = self.T; % the ambient temperature
                case 'approx'
                    TPR = self.approxPRTemp();
                case 'exact'
                    TPR = self.averagePRTemp();
            end
            resistance = resistance*(1 + TPR*self.TCR);
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
            
            actuator_noise_integrated = max(0, self.rms_actuator_displacement_noise* ...
                spring_constant*force_sensitivity); % V
            johnson_integrated = sqrt(4*self.k_b*TPR*resistance* ...
                (self.freq_max - self.freq_min));
            hooge_integrated = sqrt(self.alpha()*self.v_bridge^2* ...
                self.number_of_piezoresistors./(4*self.number_of_carriers()) ...
                *log(self.freq_max/self.freq_min));
            thermo_integrated = sqrt(force_sensitivity^2*4*spring_constant* ...
                self.k_b*TPR/(2*pi*omega_damped_hz*Q_M)*(self.freq_max - self.freq_min));
            amplifier_integrated = sqrt(A_VJ^2*(self.freq_max - self.freq_min) + ...
                A_VF^2*log(self.freq_max/self.freq_min) + 2*(resistance/2*A_IJ)^2* ...
                (self.freq_max - self.freq_min) + 2*(resistance/2*A_IF)^2* ...
                log(self.freq_max/self.freq_min));
            integrated_noise = sqrt(actuator_noise_integrated^2 + ...
                johnson_integrated^2 + hooge_integrated^2 + thermo_integrated^2 ...
                + amplifier_integrated^2);
        end
        
        % Calculate the noise in V/rtHz at a given frequency
        function voltage_noise = voltage_noise(self, freq)
            voltage_noise = sqrt(self.johnson_PSD(freq) + self.hooge_PSD(freq) + ...
                self.thermo_PSD(freq) + self.amplifier_PSD(freq));
        end
        
        function f_min_cumulative = f_min_cumulative(self)
            frequency = logspace(log10(self.freq_min), log10(self.freq_max), ...
                cantilever.numFrequencyPoints);
            noise = self.voltage_noise(frequency);
            sensitivity = self.force_sensitivity();
            force_noise_density = noise./sensitivity;
            f_min_cumulative = sqrt(cumtrapz(frequency, force_noise_density.^2));
        end
        
        
        % Piezoresistance factor
        % Accounts for dopant concentration dependent piezoresistivity in silicon
        % Source: "Piezoresistance in p-type silicon revisited", Richter et al.
        function P = piezoresistance_factor(self, dopant_concentration)
            Nb = 6e19;
            Nc = 7e20;
            richter_alpha = 0.43;
            richter_gamma = 1.6;
            richter_beta = 0.1;
            richter_eta = 3;
            richter_theta = 0.9;
            
            T0 = 300;
            switch self.thermal_modeling
                case 'none'
                    average_PR_temp = self.T; % the ambient temperature
                case 'approx'
                    average_PR_temp = self.approxPRTemp();
                case 'exact'
                    average_PR_temp = self.averagePRTemp();
            end
            Theta = average_PR_temp/T0;
            
            P = Theta^-richter_theta*(1 + Theta^-richter_beta* ...
                (dopant_concentration/Nb).^richter_alpha + Theta^-richter_eta* ...
                (dopant_concentration/Nc).^richter_gamma).^-1;
        end
        
        % Low concentration longitudinal piezoresistance coefficient
        % Units: 1/Pa
        function max_factor = max_piezoresistance_factor(self)
            switch self.doping_type
                case 'boron'
                    max_factor = 72e-11; % 110 direction
                case 'phosphorus'
                    max_factor = 103e-11; % 100 direction
                case 'arsenic'
                    max_factor = 103e-11; % 100 direction
            end
        end
        
        % Calculate the sensitivity factor (beta*)
        % Accounts for finite piezoresistor thickness
        % Units: None
        function beta = beta(self)
            [z, active_doping, total_doping] = self.doping_profile();
            
            % Shift axis so that z vaies from [-t/2, t/2] and m -> cm
            z = (self.t/2 - z)*1e2;
            
            switch self.temperature_dependent_properties
                case 'yes'
                    [mu, sigma] = self.mobility(active_doping, self.approxPRTemp());
                case 'no'
                    [mu, sigma] = self.mobility(active_doping, self.T);
            end
            P = self.piezoresistance_factor(active_doping);
            numerator = trapz(z, sigma.*P.*z);
            denominator = trapz(z, sigma);
            beta = 2*numerator/(self.t*1e2*denominator); % t: m -> cm
            
            % For optimization, ensure that beta doesn't become negative
            beta = max(beta, 1e-6);
        end
        
        % Ratio of piezoresistor resistance to total resistance (< 1)
        function gamma = gamma(self)
            R = self.resistance();
            gamma = R/(R + 2*self.R_contact);
        end
        
        % Calculate the force sensitivity (V/N) for the 1/4-active Wheatstone bridge
        % Includes the transverse portion at the end of the piezoresistive loop
        % Units: V/N
        function force_sensitivity = force_sensitivity(self)
            
            % For speed: precompute these parameters
            betaStar = self.beta();
            Rs = self.sheet_resistance();
            wheatstone_bridge_sensitivity = self.v_bridge/4;
            piMax = self.max_piezoresistance_factor();
            gamma = self.gamma();
            R = self.resistance();
            l_pr = self.l_pr();
            w_pr = self.w_pr();
            
            % Pick the relative transverse PR coefficient for the doping type
            switch self.doping_type
                case 'boron'
                    transverse_factor = -1;
                case {'phosphorus', 'arsenic'}
                    transverse_factor = -.5;
            end
            
            % Average length from the base to the piezoresistor centroid
            longitudinal_l_avg = self.l - l_pr/2;
            transverse_l_avg = self.l - l_pr;
            
            % Stress prefactor
            stress_prefactor = 6*piMax/(self.w*self.t^2)*betaStar*gamma;
            
            % Calculate the longitudinal and transverse resistances
            % Assume that the transverse width is 2x the PR width
            R_longitudinal = 2*Rs*l_pr/w_pr;
            R_transverse = Rs*self.air_gap_width/(2*w_pr);
            
            % Calculate deltaR values
            longitudinal_deltaR = stress_prefactor*longitudinal_l_avg*R_longitudinal;
            transverse_deltaR = stress_prefactor*transverse_factor* ...
                transverse_l_avg*R_transverse;
            deltaR_R = (longitudinal_deltaR + transverse_deltaR)/R;
            
            force_sensitivity = deltaR_R*wheatstone_bridge_sensitivity;
            
            % Compare with the sensitivity excluding the transverse effect
            % (i.e. the simple formula that we usually include in papers)
            force_sensitivity_no_transverse_effect = 3*(self.l - l_pr/2)* ...
                piMax/(2*self.w*self.t^2)*betaStar*gamma*self.v_bridge;
        end
        
        % Calculate the input referred surface stress sensitivity
        % Units: V/Pa
        function stress_sensitivity = surface_stress_sensitivity(self)
            
            % Pick the relative transverse PR coefficient for the doping type
            switch self.doping_type
                case 'boron'
                    longitudinal_factor = 1;
                    transverse_factor = -1;
                case {'phosphorus', 'arsenic'}
                    longitudinal_factor = 1;
                    transverse_factor = -.5;
            end
            
            % The longitudinal and transverse stress is equal everywhere
            sensitivity_factor = abs(longitudinal_factor + transverse_factor);
            stress_sensitivity = 9*sensitivity_factor* ...
                self.max_piezoresistance_factor()*self.beta()*self.gamma()* ...
                self.v_bridge/(16*self.t);
        end
        
        % Calculate the input referred displacement sensitivity
        % Units: V/m
        function displacement_sensitivity = displacement_sensitivity(self)
            displacement_sensitivity = self.force_sensitivity()*self.stiffness();
        end
        
        % Power dissipation in the cantilever
        % Units: W
        function power_dissipation = power_dissipation(self)
            power_dissipation = (self.v_bridge/2)^2/self.resistance();
        end
        
        % Tip and maximum temperatures from the F-D model
        % Units: K
        function [TMax, TTip] = calculateMaxAndTipTemp(self)
            [tmp, Q, temp] = self.calculateTempProfile();
            TMax = max(temp);
            TTip = temp(end);
        end
        
        % Calculate the approximate PR temperature
        % Units: K
        function TPR = approxPRTemp(self)
            [TMax TTip] = self.approxTempRise();
            TPR = self.T + TMax/2;
        end
        
        % Calculate the exact average PR temperature
        % Units: K
        function TPR = averagePRTemp(self)
            switch self.temperature_dependent_properties
                case 'yes'
                    [x, Q, temp] = self.calculateTempProfileTempDependent();
                case 'no'
                    [x, Q, temp] = self.calculateTempProfile();
            end
            pr_indices = intersect(find(x >= self.l_a), ...
                find(x <= (self.l_a + self.l_pr())));
            TPR = self.T + mean(temp(pr_indices));
        end
        
        % Calculate the temperature at the base of the PR
        % Useful for the designing combined sensors/actuators
        % Units: K
        function TBase = tempRiseAtPRBase(self)
            [x, Q, temp] = self.calculateTempProfile();
            base_index = find(x >= self.l_a, 1);
            TBase = temp(base_index);
        end
        
        % Calculate the maximum piezoresistor temperature
        % Units: K
        function TPR = maxPRTemp(self)
            [x, Q, temp] = self.calculateTempProfile();
            pr_indices = intersect(find(x >= self.l_a), ...
                find(x <= (self.l_a + self.l_pr())));
            TPR = max(temp(pr_indices));
        end
        
        % Calculate the average temperature increase of the actuator
        % Use for calculating A_XK (nm/K) for the thermal actuators
        % Units: K
        function TActuator = averageActuatorDeltaTemp(self)
            [x, Q, temp] = self.calculateTempProfile();
            actuator_indices = x <= (self.l_a);
            TActuator = mean(temp(actuator_indices));
        end
        
        % Calculate the temp change of the PR in response to thermal actuation
        % Units: K
        function PRTempDelta = thermalCrosstalk(self)
            temp_hot = self.averagePRTemp();
            v_actuator_temp = self.v_actuator;
            self.v_actuator = 0;
            temp_cold = self.averagePRTemp();
            self.v_actuator = v_actuator_temp; %#ok<*MCHV2>
            PRTempDelta = temp_hot - temp_cold;
        end
        
        % Calculate the approx max and tip temperatures using lumped modeling
        % Useful for quickly approximating the important temperatures
        % Units: K
        function [TMax TTip] = approxTempRiseAnalytical(self)
            k_c = self.k_base();
            R_conduction_pr  = self.l_pr()/(2*self.w*self.t*k_c);
            TMax = self.power_dissipation()*R_total;
            [l_healing, tmp] = thermalHealingLengths(self);
            TTip = TMax_analytical*exp(-(self.l - 2/3*self.l_pr())/l_healing);
        end
        
        % Approximate temperature modeling via a circuit model
        % Much faster than the F-D model and slightly more accurate than the
        % lumped parameter model
        % Units: K
        function [TMax TTip] = approxTempRise(self)
            h = self.lookupHeff();
            k_c = self.k_base();
            l_pr = self.l_pr();
            W = self.power_dissipation();
            
            % Model the system as current sources (PR or heater) and resistors
            [E_metal, rho_metal, k_metal, alpha_metal] = self.lookup_metal_properties();
            
            switch self.cantilever_type
                case 'none'
                    switch self.fluid
                        case 'vacuum'
                            R_conduction_pr  = l_pr/(2*self.w*self.t*k_c) + self.R_base;
                            TMax = W*R_conduction_pr;
                            TTip = TMax;
                        otherwise
                            R_conduction_pr  = l_pr/(2*self.w*self.t*k_c) + self.R_base;
                            R_convection_pr = 1/(2*h*l_pr*(self.w + self.t));
                            R_conduction_tip  = (self.l - l_pr)/(2*self.w*self.t*k_c);
                            R_convection_tip = 1/(2*h*(self.l-l_pr)*(self.w + self.t));
                            R_total = 1/(1/R_conduction_pr + 1/R_convection_pr + ...
                                1/(R_conduction_tip + R_convection_tip));
                            
                            TMax = W*R_total;
                            TTip = W*R_total/(R_conduction_tip + ...
                                R_convection_tip)*R_convection_tip;
                    end
                case 'step'
                    R_conduction_pr  = l_pr/(2*self.w*self.t*k_c) + ...
                        self.l_a/(self.w_a*(self.t*k_c + self.t_a*k_metal)) + ...
                        self.R_base;
                    R_convection_pr = 1/(2*h*(l_pr+self.l_a)*(self.w + self.t));
                    R_conduction_tip  = (self.l - l_pr)/(2*self.w*self.t*k_c);
                    R_convection_tip = 1/(2*h*(self.l-l_pr)*(self.w + self.t));
                    R_total = 1/(1/R_conduction_pr + 1/R_convection_pr + ...
                        1/(R_conduction_tip + R_convection_tip));
                    
                    TMax = W*R_total;
                    TTip = W*R_total/(R_conduction_tip + ...
                        R_convection_tip)*R_convection_tip;
                case 'piezoelectric'
                    R_conduction_pr  = l_pr/(2*self.w*self.t*k_c) + ...
                        self.l_a/(self.w_a*(k_c*self.t + self.k_aln*(self.t_a + ...
                        self.t_a_seed) + k_metal*(self.t_electrode_bottom + ...
                        self.t_electrode_top))) + self.R_base;
                    R_convection_pr = 1/(2*h*l_pr*(self.w + self.t));
                    R_conduction_tip  = (self.l - l_pr)/(2*self.w*self.t*k_c);
                    R_convection_tip = 1/(2*h*(self.l-l_pr)*(self.w + self.t));
                    R_total = 1/(1/R_conduction_pr + 1/R_convection_pr + ...
                        1/(R_conduction_tip + R_convection_tip));
                    
                    TMax = W*R_total;
                    TTip = W*R_total/(R_conduction_tip + ...
                        R_convection_tip)*R_convection_tip;
                case 'thermal'
                    R_conduction_pr  = l_pr/(2*self.w*self.t*k_c) + ...
                        self.l_a/(self.w_a*(self.t*k_c + self.t_a*k_metal));
                    R_convection_pr = 1/(2*h*l_pr*(self.w + self.t));
                    R_conduction_tip  = (self.l - l_pr)/(2*self.w*self.t*k_c);
                    R_convection_tip = 1/(2*h*(self.l-l_pr)*(self.w + self.t));
                    R_conduction_heater = self.l_a/(2*self.w_a*(self.t*k_c + ...
                        self.t_a*k_metal));
                    R_convection_heater = 1/(2*h*self.l_a*(self.w_a + self.t_a));
                    R_total = 1/(1/(R_conduction_pr + 1/(1/R_convection_heater + ...
                        1/R_conduction_heater)) + 1/R_convection_pr + ...
                        1/(R_conduction_tip + R_convection_tip));
                    
                    T_heater = W/(1/R_convection_heater + ...
                        1/R_conduction_heater);
                    TMaxDivider = 1/(1/R_convection_pr + 1/(R_conduction_tip + ...
                        R_convection_tip)) / (R_conduction_pr + 1/(1/R_convection_pr + ...
                        1/(R_conduction_tip + R_convection_tip)));
                    TTipDivider = R_convection_tip/(R_convection_tip + R_conduction_tip);
                    TMax = T_heater*TMaxDivider + W*R_total;
                    TTip = T_heater*TMaxDivider*TTipDivider + ...
                        W*R_total/(R_conduction_tip + ...
                        R_convection_tip)*R_convection_tip;
            end
        end
        
        % Calculate the approx thermal healing lengths for the cantilever
        % Units: m
        function [l_healing_cantilever, l_healing_step] = thermalHealingLengths(self)
            A = self.w*self.t;
            P = 2*(self.w + self.t);
            h = self.lookupHeff();
            
            k_c = self.k_base();
            l_healing_cantilever = sqrt(k_c*A/h/P);
            
            [E_metal, rho_metal, k_metal, alpha_metal] = self.lookup_metal_properties();
            l_healing_step = 0;
            switch self.cantilever_type
                case {'step', 'thermal'}
                    l_healing_step = sqrt(self.w_a*(k_c*self.t + ...
                        self.k_sio2*self.t_oxide + k_metal*self.t_a)/ ...
                        (2*(self.w_a+self.t+self.t_a)*h));
                case 'piezoelectric'
                    l_healing_step = sqrt(self.w_a*(k_c*self.t + ...
                        self.k_sio2*self.t_oxide + self.k_aln*(self.t_a + self.t_a_seed) + ...
                        k_metal*(self.t_electrode_bottom + self.t_electrode_top)) / ...
                        (2*(self.w_a+self.t+self.t_a)*h));
            end
        end
        
        % Model the temp profile from Joule heating via finite differences
        % Assumes convection to ambient, adiabatic tip, and R_base to the
        % silicon die which is clamped at the ambient temperature
        function [x, Q, T_increase] = calculateTempProfile(self, varargin)
            
            % There are several ways to call calculateTempProfile()
            % No arguments: temperature independent solution
            % One argument (legacy): temp-dependent thermal conductivity
            % Two arguments: temp-dependent sheet resistance and conductivity
            switch length(varargin)
                case 0
                    k_c = self.k_base();
                    Rsheet = self.sheet_resistance();
                    k_x = ones(1, self.numXPoints)*k_c;
                    Rsheet_x = ones(1, self.numXPoints)*Rsheet;
                case 1
                    k_x = varargin{1};
                    Rsheet = self.sheet_resistance();
                    Rsheet_x = ones(1, self.numXPoints)*Rsheet;
                case 2
                    k_x = varargin{1};
                    Rsheet_x = varargin{2};
            end
            
            % Discretize the length of the cantilever
            n_points = self.numXPoints;
            totalLength = self.l + self.l_a;
            dx = totalLength/(n_points - 1);
            x = 0:dx:totalLength;
            
            % Determine the step and PR indices
            step_indices = find(x <= self.l_a);
            actuator_indices = find(x <= (self.l_a - self.l_a_gap));
            cantilever_indices = find(x > self.l_a);
            pr_indices = intersect(cantilever_indices, ...
                find(x < (self.l_a + self.l_pr())));
            
            % Calculate Qgen_x differently depending on our temperature range
            switch self.temperature_dependent_properties
                % Calculate Qgen(x) considering temperature dependent sheet resistance
                % This method does not converge quickly for design optimization,
                % so is best used for modeling
                case 'yes'
                    index_range = x <= self.l_pr();
                    R_x = 2*Rsheet_x(index_range)/self.w_pr();
                    R_calc = 2*trapz(x(index_range), Rsheet_x(index_range)/self.w_pr());
                    I_calc = (self.v_bridge/2)/R_calc;
                    Qgen_x = I_calc^2*R_x;
                    
                    % Assume power/length is constant along the piezoresistor length
                    % This method works well for modeling near the ambient
                    % temperature and for design optimization
                    % Note: calculate Qgen_x based upon the actual lengths
                    % to avoid discretization errors (line 2 here is very important)
                case 'no'
                    power = (self.v_bridge/2)^2/self.resistance();
                    Qgen_x = power/(x(pr_indices(end)) - x(pr_indices(1)))*ones(length(x),1);
            end
            
            % Setup other variables
            tempAmbient = self.T;
            h = self.lookupHeff();
            [E_metal, rho_metal, k_metal, alpha_metal] = self.lookup_metal_properties();
            K = self.w.*k_x'.*self.t.*ones(n_points, 1);
            perimeter = 2*(self.w + self.t)*ones(n_points, 1);
            Q = zeros(n_points, 1);
            Q(pr_indices) = Qgen_x(pr_indices);
            
            % Build K (area*k_c) and P
            switch self.cantilever_type
                case 'none'
                case 'step'
                    K(step_indices) = self.w_a.*(k_x(step_indices)*self.t ...
                        + k_metal*self.t_a);
                    perimeter(step_indices) = 2*(self.w_a + self.t_a);
                case 'thermal'
                    Qheater = self.heaterPower()/(x(actuator_indices(end)) - ...
                        x(actuator_indices(1)));
                    Q(actuator_indices) = Qheater;
                    K(step_indices) = self.w_a*(k_x(step_indices)*self.t + ...
                        self.t_oxide*self.k_sio2 + k_metal*self.t_a);
                    perimeter(step_indices) = 2*(self.w_a + self.t_a);
                case 'piezoelectric'
                    K(step_indices) = self.w_a*(k_x(step_indices)*self.t + ...
                        self.k_aln*(self.t_a + self.t_a_seed) + self.t_oxide*self.k_sio2 + ...
                        k_metal*(self.t_electrode_top + self.t_electrode_bottom));
                    perimeter(step_indices) = 2*(self.w_a + self.t_a);
            end
            
            % Build A and RHS matrices
            A = zeros(n_points, n_points);
            rhs = zeros(n_points, 1);
            for ii = 2:n_points-1
                A(ii, ii-1) = -K(ii-1)/dx^2;
                A(ii, ii)   = (K(ii-1) + K(ii+1))/dx^2 + h*perimeter(ii);
                A(ii, ii+1) = -K(ii+1)/dx^2;
                rhs(ii, 1) = Q(ii) + h*perimeter(ii)*tempAmbient;
            end
            A(n_points, n_points-1:n_points) = [1 -1]; % Adiabatic at tip
            A = sparse(A); % Leads to a significant speed improvement
            
            % Properly handle R_base
            A(1, 1) =  -1-self.R_base*K(1)/dx;
            A(1, 2) =     self.R_base*K(1)/dx;
            
            rhs(1, 1) = -self.T;
            T_absolute = A \ rhs;
            T_increase = T_absolute - tempAmbient;
        end
        
        % Model the temperature profile of a tip-heated cantilever
        % via finite differences. Intended for modeling cantilever heating
        % during laser doppler vibrometer or AFM measurements.
        % For optical heating applications the user needs to account for
        % reflections and limited silicon absorption.
        % Assumptions:
        % 1) Piezoresistor is turned off
        % 2) Temperature-independent material properties
        % 3) Beam intensity is constant over the beam diameter
        % 4) Accounts for the spot size being smaller/larger than the beam
        % 5) Spot size = spot diameter
        % Assumes temperature-independent properties for simplicity
        % Input parameters:
        % inputPower (W): power input to the cantilever
        % spotSize (m): length back from the tip over which the power is spread
        %
        % Note: most of the code is reused from above
        function [x, Q, T_increase] = calculateTempProfileTipHeat(self, ...
                inputPower, spotSize)
            
            k_c = self.k_base();
            k_x = ones(1, self.numXPoints)*k_c;
            n_points = self.numXPoints;
            totalLength = self.l + self.l_a;
            dx = totalLength/(n_points - 1);
            x = 0:dx:totalLength;
            step_indices = find(x <= self.l_a);
            actuator_indices = find(x <= (self.l_a - self.l_a_gap));
            cantilever_indices = find(x > self.l_a);
            pr_indices = intersect(cantilever_indices, ...
                find(x < (self.l_a + self.l_pr())));
            
            % Calculate the input power at the tip.
            % Accounts for some fraction of the beam missing the cantilever
            % if the cantilever width is comparable to the beam diameter
            % Assumes the center is aimed at the tip
            Q = zeros(n_points, 1);
            fractionAbsorbed = min((spotSize/2)*self.w/(pi*(spotSize/2)^2), 1);
            Q(x > (self.l - spotSize)) = fractionAbsorbed*inputPower/spotSize;
            
            tempAmbient = self.T;
            h = self.lookupHeff();
            [E_metal, rho_metal, k_metal, alpha_metal] = self.lookup_metal_properties();
            K = self.w.*k_x'.*self.t.*ones(n_points, 1);
            perimeter = 2*(self.w + self.t)*ones(n_points, 1);
            
            switch self.cantilever_type
                case 'none'
                case 'step'
                    K(step_indices) = self.w_a.*(k_x(step_indices)*self.t + ...
                        k_metal*self.t_a);
                    perimeter(step_indices) = 2*(self.w_a + self.t_a);
                case 'thermal'
                    Qheater = self.heaterPower()/(x(actuator_indices(end)) - ...
                        x(actuator_indices(1)));
                    Q(actuator_indices) = Qheater;
                    K(step_indices) = self.w_a_active*(k_x(step_indices)*self.t + ...
                        self.t_oxide*self.k_sio2 + k_metal*self.t_a);
                    perimeter(step_indices) = 2*(self.w_a + self.t_a);
                case 'piezoelectric'
                    K(step_indices) = self.w_a*(k_x(step_indices)*self.t + ...
                        self.k_aln*(self.t_a + self.t_a_seed) + self.t_oxide*self.k_sio2 + ...
                        k_metal*(self.t_electrode_top + self.t_electrode_bottom));
                    perimeter(step_indices) = 2*(self.w_a + self.t_a);
            end
            
            A = zeros(n_points, n_points);
            rhs = zeros(n_points, 1);
            for ii = 2:n_points-1
                A(ii, ii-1) = -K(ii-1)/dx^2;
                A(ii, ii)   = (K(ii-1) + K(ii+1))/dx^2 + h*perimeter(ii);
                A(ii, ii+1) = -K(ii+1)/dx^2;
                rhs(ii, 1) = Q(ii) + h*perimeter(ii)*tempAmbient;
            end
            A(n_points, n_points-1:n_points) = [1 -1];
            A = sparse(A);
            A(1, 1) =  -1-self.R_base*K(1)/dx;
            A(1, 2) =     self.R_base*K(1)/dx;
            rhs(1, 1) = -self.T;
            T_absolute = A \ rhs;
            T_increase = T_absolute - tempAmbient;
        end
        
        
        % Calculate the temp profile accounting for k(T) and R_s(T)
        % If temperature_dependent_properties is set to 'yes' then this
        % function should be called correctly. The temperature dependent option
        % is much slower than its temperature independent equivalent.
        function [x, Q, T_final] = calculateTempProfileTempDependent(self)
            
            % Iterative calculate the temperature distribution and k(x)
            % Initialize k(x) solver with our initial guess
            % This approach is required because the system is nonlinear
            [x, Q, T_initial] = self.calculateTempProfile();
            T_current = T_initial;
            residual = 1;
            
            while residual > 1e-6
                absoluteTemp = T_current + self.T;
                
                % Calculate k_si and Rsheet from the current temperature guess
                [x, z, k_z_x] = self.thermal_conductivity_profile(absoluteTemp);
                k_x = mean(k_z_x, 1); % Average k for a cross-section
                Rsheet_x = self.RSheetProfile(x, T_current);
                
                % Scale the resistance so that it matches our nominal value
                % This choice is geared towards modeling experimental results
                R_nominal = self.resistance();
                x_resistor = x(x <= self.l_pr());
                Rsheet_resistor = Rsheet_x(x <= self.l_pr());
                R_calc = 2*trapz(x_resistor, Rsheet_resistor/self.w_pr()) + ...
                    3.4*Rsheet_resistor(end) + 2*self.R_contact;
                Rsheet_x = Rsheet_x * (R_nominal/R_calc); % Scale
                
                % Update the temperature guess using k_x, Rsheet_x
                [x, Q, T_new] = self.calculateTempProfile(k_x, Rsheet_x);
                
                residual = sum((T_new - T_current).^2);
                T_current = T_new;
            end
            T_final = T_current;
        end
        
        % Calculate the approximate thermal conductivity of the cantilever
        % Units: W/m-K
        function k_base = k_base(self)
            k_x = self.k_x();
            k_base = k_x(1);
        end
        
        % Calculate k(x) as a function of cantilever thickness
        % Units: W/m-K
        function k_effective = k_x(self) %#ok<*MANU>
            switch self.thermal_modeling
                case {'none', 'approx'}
                    % Model from Asheghi (1997)
                    t0 = 120e-9;
                    k_effective = self.k_si*self.t./(t0 + self.t);
                case 'exact'
                    [x, z, k] = self.thermal_conductivity_profile(self.T);
                    k_effective = mean(k);
            end
        end
        
        % Calculate k(x,z) using a temp/thickness dependent model
        % Source: "Modeling and Data for Thermal Conductivity ...", Liu et al, 2006.
        % Units: W/m-K
        function [x, z, k] = thermal_conductivity_profile(self, T)
            % Phonon Group velocities (m/sec)
            nu_T =  5860;
            nu_TU = 2000;
            nu_L = 8480;
            
            % Density of state coefficients
            C_T = self.k_b/(2*pi^2*nu_T)*(self.k_b/self.h_bar)^3;
            C_TU = self.k_b/(2*pi^2*nu_TU)*(self.k_b/self.h_bar)^3;
            C_L = self.k_b/(2*pi^2*nu_L)*(self.k_b/self.h_bar)^3;
            
            % Debye tempratures
            theta_1 = 180; %K
            theta_2 = 210; %K
            theta_3 = 570; %K
            
            % Phonon specific heat/volume
            c_TU = 7.5e5; % J/m^3-K
            c_L = 3.2e5; % J/m^3-K
            
            % Other Constants
            B_T = 9.3e-13; % 1/K^3
            B_TU = 5.5e-18; % s
            B_L = 2e-24; % s/K
            beta_TU = B_TU*(self.k_b/self.h_bar)^2; % 1/s-K^2
            beta_L = B_L*(self.k_b/self.h_bar)^2; % 1/s-K^3
            beta_T = B_T*(self.k_b/self.h_bar); % 1/s-K^4
            
            % Phonon frequencies
            omega_2 = self.k_b/self.h_bar*theta_2;
            omega_3 = self.k_b/self.h_bar*theta_3;
            
            % Crystal lattice parameters
            Gamma = 2.1e-4; % computed from isotopes
            nu_s = (1/3*(nu_L^-1 + 2*nu_T^-1))^-1; % Speed of sound
            
            R_si = 5.4e-10; %m - radius of the silicon host atom
            V_si = (R_si/2)^3; % m^3 - silicon lattice volume
            M_si = 4.59e-26; % kg - mass of silicon host atom
            A_isotopes = Gamma*V_si/(4*pi*nu_s^3);
            
            % Other constants
            Q = 4.2; % Constant scattering matrix parameter
            gamma = 1.6; % Gruneisen parameter
            psi = 3.5; % Modification factor to account for inconsistency of the Si
            
            % Phonon mean free paths in nearly pure bulk silicon
            function lambda_TU = lambda_TU(T)
                lambda_TU = C_TU*(theta_2^2-theta_1^2)./(2*nu_TU*c_TU*beta_TU.*T);
            end
            
            function lambda_L = lambda_L(T)
                lambda_L = C_L*theta_3./(nu_L*c_L*beta_L*T.^3);
            end
            
            % Phonon mean free paths including thermal vibrations
            function lambda_TU_HT = lambda_TU_HT(T)
                lambda_TU_HT = (1./lambda_TU(T) + ...
                    1/(nu_TU*(A_isotopes*omega_2^4)^-1)).^-1;
            end
            
            function lambda_L_HT = lambda_L_HT(T)
                lambda_L_HT = (1./lambda_L(T) + ...
                    1/(nu_L*(A_isotopes*omega_3^4)^-1)).^-1;
            end
            
            % Contribution to thermal conductivity from phonon modes
            function k_D_layer_TU_HT = k_D_layer_TU_HT(T, d, n, dM, dR)
                A = 1./lambda_TU_HT(T) + ...
                    1./(nu_TU*((n*V_si^2./(4*pi*nu_s^3).*(dM/M_si)^2 + ...
                    2*n*V_si^2*Q^2*gamma^2./(pi*nu_s^3)*(dR/R_si)^2)*omega_2^4).^-1);
                k_D_layer_TU_HT = 2/3.*c_TU.*nu_TU.*A.^-1.*(1+3./(8*d./(A.^-1*psi))).^-1;
            end
            
            function k_D_layer_L_HT = k_D_layer_L_HT(T, d, n, dM, dR)
                A = 1./lambda_L_HT(T) + ...
                    1./(nu_L*((n*V_si^2/(4*pi*nu_s^3).*(dM/M_si)^2 + ...
                    2.*n*V_si^2*Q^2*gamma^2./(pi*nu_s^3)*(dR/R_si)^2)*omega_3^4).^-1);
                k_D_layer_L_HT = 1/3.*c_L.*nu_L.*A.^-1.*(1+3./(8*d./(A.^-1*psi))).^-1;
            end
            
            function k_D_layer_HT = k_D_layer_HT(T, d, n, dM, dR)
                n = n*1e6; % Convert from 1/cc to 1/m^3
                k_D_layer_HT = k_D_layer_L_HT(T, d, n, dM, dR) + ...
                    k_D_layer_TU_HT(T, d, n, dM, dR);
            end
            
            % Doping constants
            switch self.doping_type
                case 'arsenic'
                    fprintf('ERROR: Arsenic dM and dR values not included');
                    pause
                case 'boron'
                    dM = 18*1.673e-27;
                    dR = 0.175e-10;
                case 'phosphorus'
                    dM = 3*1.673e-27;
                    dR = 0.175e-10;
            end
            
            % Vector calculation of thermal conductivity
            % Roughly 10x faster than iterative calculation!
            [z, active_doping, total_doping] = self.doping_profile();
            x = linspace(0, self.l, self.numXPoints);
            
            % let dim(T) and dim(total_doping) be numZPoints x numXPoints
            n_matrix = total_doping'*ones(1, self.numXPoints);
            if length(T) == 1
                T_matrix = ones(self.numZPoints, 1)*T'*ones(1, self.numXPoints);
            else
                T_matrix = ones(self.numZPoints, 1)*T';
            end
            
            % TODO: Refactor to remove all of those unnecessary function calls
            k = k_D_layer_HT(T_matrix, self.t, n_matrix, dM, dR);
        end
        
        % ====== Calculate resolution ======
        function force_resolution = force_resolution(self)
            force_resolution = self.integrated_noise()/self.force_sensitivity();
        end
        
        function force_resolution = thermomechanical_force_noise(self)
            force_resolution = self.thermo_integrated()/self.force_sensitivity();
        end
        
        function displacement_resolution = displacement_resolution(self)
            displacement_resolution = self.force_resolution()/self.stiffness();
        end
        
        function force_noise_density = force_noise_density(self, freq)
            force_noise_density = self.voltage_noise(freq)/self.force_sensitivity();
        end
        
        function force_noise_density = resonant_force_noise_density(self)
            [omega_damped_hz, Q] = self.omega_damped_hz_and_Q(); %#ok<*NASGU>
            freq = omega_damped_hz;
            force_noise_density = self.force_noise_density(freq)*1e15;
        end
        
        function surface_stress_resolution = surface_stress_resolution(self)
            surface_stress_resolution = self.integrated_noise()/ ...
                self.surface_stress_sensitivity();
        end
        
        
        
        % ====== Multilayer beam mechanics and actuation ======
        
        % Calculate the cantilever curvature per unit applied moment
        function Cm = calculateActuatorNormalizedCurvature(self)
            Zm = self.actuatorNeutralAxis();
            [z, E, A, I] = self.lookupActuatorMechanics();
            Z_offset = z - Zm;
            Cm = 1./sum(E.*(I + A.*Z_offset.^2));
        end
        
        % Calculate the thickness of silicon required to have the same EI
        % as the actuator/step at the base
        function t_equivalent = calculateEquivalentThickness(self)
            function residual_squared = findEIResidual(t_equivalent_guess)
                EI_calculated = self.modulus()*self.w_a*t_equivalent_guess^3/12;
                residual_squared = (EI_effective/EI_calculated - 1)^2;
            end
            
            EI_effective = 1/self.calculateActuatorNormalizedCurvature();
            options = optimset('TolX', 1e-12, 'TolFun', 1e-12, 'Display', 'off');
            t_equivalent = fminbnd(@findEIResidual, self.t/100, self.t*100, options);
        end
        
        function Zm = actuatorNeutralAxis(self)
            [z, E, A, I] = self.lookupActuatorMechanics();
            Zm = sum(z.*E.*A)/sum(E.*A);
        end
        
        function d31 = d31(self)
            if self.d31_manual == 0
                d31 = interp1(self.d31_t, self.d31_aln, self.t_a, 'spline');
            else
                d31 = self.d31_manual;
            end
        end
        
        function [x, deflection] = calculateDeflection(self)
            n_points = self.numXPoints;
            totalLength = self.l + self.l_a;
            dx = totalLength/(n_points - 1);
            x = 0:dx:totalLength;
            
            M = 0; % external moment is zero
            P = 0; % external load is zero
            
            [z, E, A, I] = self.lookupActuatorMechanics();
            stress = self.calculateActuatorStress();
            
            % Calculate the curvature and neutral axis
            % The curvature may vary as a function of position (e.g. thermal)
            % so calculate the deflection by integrating the curvature twice
            C = zeros(length(x), 1);
            
            % Calculate the curvature, C along the cantilever length
            for ii = 1:length(x)
                C(ii) = -((M - sum(z.*A.*stress(ii,:)))*sum(E.*A) + ...
                    (P + sum(A.*stress(ii,:)))*sum(E.*z.*A))/ ...
                    (sum(E.*A)*sum(E.*(I+A.*z.^2)) - sum(z.*E.*A)^2);
            end
            
            % No curvature beyond the end of the actuator
            % i.e. no dopant bending included in this calculation
            C(x > self.l_a) = 0;
            
            % Calculate the deflection from the curvature by integrating
            theta = cumtrapz(x, C);
            deflection = cumtrapz(x, tan(theta));
        end
        
        % Calculate the stress as a function along the actuator length
        function stress = calculateActuatorStress(self)
            [z, E, A, I] = self.lookupActuatorMechanics(); %#ok<*ASGLU>
            [E_metal, rho_metal, k_metal, alpha_metal] = self.lookup_metal_properties();
            [x_temp, Q, temp] = self.calculateTempProfile();
            temp = temp + self.T; % Use absolute temperature for these calculations
            
            intrinsic_stress = ones(self.numXPoints, 1)*self.film_intrinsic_stress();
            intrinsic_stress(x_temp > self.l_a, :) = 0;
            
            switch self.cantilever_type
                case {'step', 'thermal'}
                    cte = [self.alpha_si self.alpha_sio2 alpha_metal];
                    thermal_stress = (temp-self.T_ref).*ones(self.numXPoints,1)*(cte.*E);
                    thermal_stress(x_temp > self.l_a, 2:3) = 0;
                    piezoelectric_stress = zeros(self.numXPoints, 3);
                    
                case 'piezoelectric'
                    cte = [self.alpha_si self.alpha_sio2 self.alpha_aln ...
                        alpha_metal self.alpha_aln alpha_metal];
                    thermal_stress = (temp-self.T_ref).*ones(self.numXPoints,1)*(cte.*E);
                    thermal_stress(x_temp > self.l_a, 2:6) = 0;
                    
                    % Calculate the piezoelectric stress
                    % The seed electric field will vary depending on
                    E_field = [0 0 0 0 self.v_actuator/self.t_a 0]';
                    d31 =     [0 0 0 0 self.d31()               0]';
                    
                    piezoelectric_stress = ones(self.numXPoints,1)*E*d31*E_field';
                    piezoelectric_stress(x_temp > self.l_a - self.l_a_gap, :) = 0;
            end
            
            % Scale the stress by the active width of the device
            stress = intrinsic_stress + thermal_stress + ...
                self.w_a_active/self.w_a*piezoelectric_stress;
        end
        
        function [mu_z sigma_z] = tip_deflection_distribution(self)
            v_actuator_temp = self.v_actuator;
            self.v_actuator = 0;
            self.film_stress = 'random';
            for ii = 1:cantilever.numRandomStressIterations
                [x_ii, z_ii] = self.calculateDeflection();
                
                x = x_ii;
                z(:, ii) = z_ii;
                
                z_tip(ii) = self.tipDeflection();
            end
            self.v_actuator = v_actuator_temp;
            self.film_stress = 'nominal';
            
            mu_z = mean(z_tip);
            sigma_z = std(z_tip);
        end
        
        function plot_tip_deflection_distribution(self)
            v_actuator_temp = self.v_actuator;
            self.v_actuator = 0;
            self.film_stress = 'random';
            for ii = 1:cantilever.numRandomStressIterations
                [x_ii, z_ii] = self.calculateDeflection();
                
                x = x_ii;
                z(:, ii) = z_ii;
            end
            self.v_actuator = v_actuator_temp;
            self.film_stress = 'nominal';
            
            figure
            plot(x_ii*1e6, z*1e6, 'Color', [11 132 199]./255);
            xlabel('Distance from base (\mum)');
            ylabel('Deflection (\mum)');
            box off;
            
        end
        
        function sigma_i = film_intrinsic_stress(self)
            % film_stress = 'nominal' => use average stress
            % film_stress = 'random' => use random, normally distributed value
            switch self.film_stress
                case 'nominal'
                    sigma = 0;
                case 'random'
                    sigma = 0.25; % 1/4 because 1/2 from diff and 1/2 from range = 2 STDevs
            end
            
            switch self.cantilever_type
                case {'step', 'thermal'}
                    sigma_i(1) = mean(self.sigma_si_range) + ...
                        sigma*randn*diff(self.sigma_si_range);
                    sigma_i(2) = mean(self.sigma_sio2_range) + ...
                        sigma*randn*diff(self.sigma_sio2_range);
                    sigma_i(3) = mean(self.sigma_al_range) + ...
                        sigma*randn*diff(self.sigma_al_range);
                case 'piezoelectric'
                    sigma_i(1) = mean(self.sigma_si_range) + ...
                        sigma*randn*diff(self.sigma_si_range);
                    sigma_i(2) = mean(self.sigma_sio2_range) + ...
                        sigma*randn*diff(self.sigma_sio2_range);
                    sigma_i(3) = mean(self.sigma_aln_range) + ...
                        sigma*randn*diff(self.sigma_aln_range);
                    sigma_i(5) = mean(self.sigma_aln_range) + ...
                        sigma*randn*diff(self.sigma_aln_range);
                    switch self.metal_type
                        case 'titanium'
                            sigma_i(4) = mean(self.sigma_ti_range) + ...
                                sigma*randn*diff(self.sigma_ti_range);
                            sigma_i(6) = mean(self.sigma_ti_range) + ...
                                sigma*randn*diff(self.sigma_ti_range);
                        case 'molybdenum'
                            sigma_i(4) = mean(self.sigma_mo_range) + ...
                                sigma*randn*diff(self.sigma_mo_range);
                            sigma_i(6) = mean(self.sigma_mo_range) + ...
                                sigma*randn*diff(self.sigma_mo_range);
                    end
            end
        end
        
        % The thermal actuator power dissipation
        % Units: W
        function power = heaterPower(self)
            power = self.v_actuator^2/self.R_heater;
        end
        
        % Calculate the time constant and -3 dB freq of the thermal actuator
        % Note that in practice this is roughly 4x too idealistic
        % Units: s and Hz
        function [tau, freq] = heaterTimeConstant(self)
            k_c = self.k_base();
            [E_metal, rho_metal, k_metal, alpha_metal] = self.lookup_metal_properties();
            h = self.lookupHeff();
            
            % Calculate the thermal resistances
            R_conduction = (self.l_a-self.l_a_gap)/(2*self.w_a*(self.t*k_c + ...
                self.t_a*k_metal + self.t_oxide*self.k_sio2));
            R_convection = 1/(2*h*self.l_a*(self.w_a + self.t_a));
            R = 1/(1/(R_conduction + self.R_base) + 1/R_convection); % K/W
            C = self.l_a*self.w_a*(self.t_a*self.rho_al*self.Cv_al + ...
                self.t*self.rho_si*self.Cv_si + ...
                self.t_oxide*self.rho_sio2*self.Cv_sio2); % J/K
            
            tau = R*C;
            freq = 1/(2*pi*tau);
        end
        
        % Calculate the rise time using the Q-dependent formula
        % Source: "Control Systems Engineering", Nise
        function t_rise = mechanicalRiseTime(self)
            [omega_d, Q] = self.omega_damped_and_Q();
            zeta = 1/(2*Q);
            t_rise = 1/omega_d*(1.76*zeta^3-0.417*zeta^2+1.039*zeta+1);
        end
        
        % Calculate dz_tip/dT_ambient -> temperature stability of the cantilever
        % Units: m/K
        function dzdT = tipTempDeflectionVsAmbientTemp(self)
            temp_delta_values = -5:5;
            
            T_original= self.T;
            for ii = 1:length(temp_delta_values)
                self.T = T_original + temp_delta_values(ii);
                tipDeflection(ii) = self.tipDeflection();
            end
            self.T = T_original;
            
            p = polyfit(temp_delta_values, tipDeflection, 1);
            dzdT = p(1);
        end
        
        % Calculate the electrical current in the heater
        % Units: A
        function current = heaterCurrent(self)
            current = self.v_actuator/self.R_heater;
        end
        
        % Calculate the tip deflection including all of the various effects
        % Units: m
        function z_tip = tipDeflection(self)
            [x, z] = self.calculateDeflection();
            z_tip = z(end);
        end
        
        % Calculate the tip deflection with the actuator on vs. off
        % This approach cancels any static deflections (e.g. dopant bending)
        % Units: m
        function z_tip = actuatorTipDeflectionRange(self)
            v_actuator_temp = self.v_actuator;
            self.v_actuator = 0;
            [x, z_off] = self.calculateDeflection();
            z_initial = z_off(end);
            
            self.v_actuator = v_actuator_temp;
            [x, z_on] = self.calculateDeflection();
            z_final = z_on(end);
            
            z_tip = abs(z_final - z_initial);
        end
        
        % Calculate dopant strain induced bending
        % The calculation method handles arbitrary strain gradients
        % Source: "Residual strain and resultant postrelease deflection ...", 1999
        function [x, C, theta, deflection] = calculateDopantBending(self)
            
            % Source: "Effects of phosphorus on stress...", JMM, 1999.
            delta = 4.5e-24;
            
            n_points = self.numXPoints;
            dx = self.l/(n_points - 1);
            x = 0:dx:self.l;
            
            [depth, active_doping, total_doping] = self.doping_profile();
            E = self.modulus();
            width = self.w; % But it doesn't matter
            thickness = self.t;
            
            % Shift coordinates so that strain is measured from the bottom
            y = linspace(0, thickness, length(depth))';
            dopant_concentration = flipud(active_doping);
            
            strain = dopant_concentration * delta;
            stress = strain*E;
            
            ybar = trapz(y, y.*stress')./trapz(y, stress');
            h1 = thickness - ybar;
            h2 = ybar;
            
            i1 = find(y > ybar);
            i2 = find(y <= ybar);
            y1 = y(i1);
            y2 = y(i2);
            stress1 = stress(i1);
            stress2 = stress(i2);
            
            sigma1 = 1/h1*trapz(y1, stress1);
            sigma2 = 1/h2*trapz(y2, stress2);
            
            ybar1 = trapz(y1, y1.*stress1')/trapz(y1, stress1');
            ybar2 = trapz(y2, y2.*stress2')/trapz(y2, stress2');
            
            I1 = width*h1^3/12;
            I2 = width*h2^3/12;
            I = width*thickness^3/12;
            
            sigmai = (sigma2 - sigma1)/(2 - h1*width*(ybar1 - ybar)^2/I1 + ...
                h2*width*(ybar - ybar2)^2/I2);
            
            M1 = sigmai*h1*width*(ybar1 - ybar);
            M2 = sigmai*h2*width*(ybar - ybar2);
            
            % Calculate the radius of curvature, assuming it is small
            R = (E*I)/(M1 + M2);
            C = zeros(length(x), 1);
            C_pr = 1/R;
            C(x <= self.l_pr()) = C_pr;
            
            theta = cumsum(C.*dx);
            deflection = cumsum(theta.*dx);
        end
        
        % Calculate just the tip deflection from dopant bending
        % Units: m
        function tip_deflection = calculateDopantTipDeflection(self)
            [x, C, theta, deflection] = self.calculateDopantBending();
            tip_deflection = max(deflection);
        end
        
        % ====== Handy plotting functions ======
        function plotTempProfile(self)
            switch self.temperature_dependent_properties
                case 'yes'
                    [x, Q, temp] = self.calculateTempProfileTempDependent();
                case 'no'
                    [x, Q, temp] = self.calculateTempProfile();
            end
            
            figure
            hold all
            plot(1e6*x, temp);
            plot(1e6*x, Q, '--');
            hold off
            xlabel('X (um)');
            ylabel('Temp Rise (K)');
            legend('Temperature', 'Power/Unit Length');
        end
        
        function plotTempProfileActuatorComparison(self)
            
            [x, Q, temp] = self.calculateTempProfile();
            
            figure
            subplot(2,1,1)
            plot(1e6*x, Q, '--');
            box off;
            
            ylabel('Power (\muW/\mum)');
            
            subplot(2,1,2)
            hold all
            plot(1e6*x, temp);
            
            v_actuator_temp = self.v_actuator;
            self.v_actuator = 0;
            [x, Q, temp] = self.calculateTempProfile();
            self.v_actuator = v_actuator_temp;
            
            plot(1e6*x, temp);
            box off;
            
            hold off
            xlabel('Distance from base (\mum)');
            ylabel('\DeltaT (K)');
        end
        
        function plot_thermal_conductivity_profile(self)
            [x, z, k] = self.thermal_conductivity_profile(self.T);
            
            figure
            plot(z*1e9, k);
            xlabel('Distance from Surface (nm)');
            ylabel('k (W/m-K)');
            ylim([0 150]);
            xlim([0 self.t*1e9]);
        end
        
        function plotDopantBending(self)
            [x, C, theta, deflection] = self.calculateDopantBending();
            
            figure
            plot(1e6*x, 1e9*deflection);
            xlabel('Distance from Base (um)');
            ylabel('Cantilever Deflection (nm)');
            box off;
        end
        
        function plotDeflectionAndTemp(self)
            
            v_actuator_temp = self.v_actuator;
            self.v_actuator = 0;
            [x, deflectionInitial] = self.calculateDeflection();
            
            self.v_actuator = v_actuator_temp; %#ok<*MCHV2>
            [x, deflectionFinal] = self.calculateDeflection();
            [x, Q, temp] = self.calculateTempProfile();
            deflection = deflectionFinal - deflectionInitial;
            
            figure
            subplot(2,1,1);
            plot(1e6*x, temp);
            ylabel('Temp Rise (K)');
            box off;
            
            subplot(2,1,2);
            plot(1e6*x, 1e6*deflection);
            xlabel('Distance from Base (um)');
            ylabel('Deflection (\mum)');
            box off;
        end
        
        function plot_noise_spectrum(self)
            freq = logspace( log10(self.freq_min), log10(self.freq_max), ...
                cantilever.numFrequencyPoints);
            
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
        
        % ======= Handy lookup functions =======
        % Calculate elastic modulus based upon dopant type
        % Assume we're using the best piezoresistor orientation.
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
        
        % Lookup the appropriate fluid physical properties
        function [rho_f, eta_f] = lookupFluidProperties(self)
            switch self.fluid
                case 'air'
                    rho_f = self.rho_air;
                    eta_f = self.eta_air;
                case 'water'
                    rho_f = self.rho_water;
                    eta_f = self.eta_water;
                case 'arbitrary'
                    rho_f = self.rho_arb;
                    eta_f = self.eta_arb;
                otherwise
                    fprintf('ERROR - Unknown fluid: %s', self.fluid);
                    pause
            end
        end
        
        % Lookup effective convection coefficient
        % Or try calculating it using a shape factor (not recommended)
        % Units: W/m^2-K
        function h_eff = lookupHeff(self)
            switch self.h_method
                case 'fixed'
                    switch self.fluid
                        case 'vacuum'
                            h_eff = self.h_vacuum;
                        case 'air'
                            h_eff = self.h_air;
                        case 'water'
                            h_eff = self.h_water;
                        case 'arbitrary'
                            h_eff = self.h_arb;
                    end
                case 'calculate'
                    switch self.fluid
                        case 'vacuum'
                            k_f = self.k_vacuum;
                        case 'air'
                            k_f = self.k_air;
                        case 'water'
                            k_f = self.k_water;
                        case 'arbitrary'
                            k_f = self.k_arb;
                    end
                    % Calculate h using the shape factor
                    A = self.w*self.t;
                    P = 2*(self.w + self.t);
                    
                    % Assume a thin rectangle to an infinite plane
                    z = 450e-6; % cantilever die thickness, m
                    S = 2*pi*self.l/log(2*pi*z/self.w);
                    
                    h_eff = k_f*S/(P*self.l); % calculate h from the Hu equation
            end
        end
        
        % Simple lookup function for the film stack properties
        function [E, rho, k, alpha] = lookup_metal_properties(self)
            switch self.metal_type
                case 'aluminum'
                    E = self.E_al;
                    rho = self.rho_al;
                    k = self.k_al;
                    alpha = self.alpha_al;
                case 'titanium'
                    E = self.E_ti;
                    rho = self.rho_ti;
                    k = self.k_ti;
                    alpha = self.alpha_ti;
                case 'molybdenum'
                    E = self.E_mo;
                    rho = self.rho_mo;
                    k = self.k_mo;
                    alpha = self.alpha_mo;
                otherwise
                    fprintf('ERROR: Unknown metal_type');
                    pause
            end
        end
        
        % Treat the actuator as encompassing the entire device width here
        % and scale the actuation later (to account for the metal traces
        % running alongside
        function [z_layers, E_layers, A, I] = lookupActuatorMechanics(self)
            [E_metal, rho_metal, k_metal, alpha_metal] = self.lookup_metal_properties();
            switch self.cantilever_type
                case 'none';
                    error('ERROR in lookupActuatorMechanics() - cantilever_type = none');
                case {'step', 'thermal'};
                    t_layers = [self.t self.t_oxide self.t_a];
                    w_layers = [self.w_a self.w_a self.w_a];
                    E_layers = [self.modulus() self.E_sio2 E_metal];
                case 'piezoelectric';
                    t_layers = [self.t         self.t_oxide self.t_a_seed ...
                        self.t_electrode_bottom self.t_a   self.t_electrode_top];
                    w_layers = [self.w_a       self.w_a     self.w_a ...
                        self.w_a 								self.w_a 		self.w_a];
                    E_layers = [self.modulus() self.E_sio2  self.E_aln ...
                        E_metal 								self.E_aln 		E_metal];
            end
            
            z_layers = zeros(1, length(t_layers));
            for ii = 1:length(t_layers)
                % z(1) = t(1)/2, z(2) = t(1) + t(2)/2
                z_layers(ii) = sum(t_layers) - sum(t_layers(ii:end)) + t_layers(ii)/2;
            end
            A = w_layers.*t_layers;
            I = (w_layers.*t_layers.^3)/12;
        end
        
        % ======== Beam mechanics, resonant frequency and damping ==========
        
        % Calculate the cantilever spring constant
        % Units: N/m
        function stiffness = stiffness(self)
            k_tip = self.modulus() * self.w * self.t^3 / (4*self.l^3);
            
            switch self.cantilever_type
                case 'none'
                    stiffness = k_tip;
                    
                    % For devices with an actuator/step at the base...
                    % Apply a test force and calculate the resulting deflection
                    % Treating the system as two springs in series is not valid because:
                    % 1) the force is applied beyond the tip of the thick base,
                    %			so the moment it sees is larger than would be expected
                    % 2) the angle at the end of the thick base is integrated
                    %			along the remaining length of the cantilever
                otherwise
                    F = 1e-9; % Test force (N)
                    EI_base = 1/self.calculateActuatorNormalizedCurvature();
                    EI_tip = self.modulus()*self.w*self.t^3/12;
                    
                    % Deflection at the tip of the base
                    z_base = F*self.l_a^2*(3*(self.l+self.l_a)-self.l_a)/(6*EI_base);
                    
                    % Angle at the tip of the base
                    theta_base = F*(2*(self.l_a+self.l)-self.l_a)*self.l_a/(2*EI_base);
                    
                    % Deflection at the very tip
                    z_tip = z_base + tan(theta_base)*self.l + F*self.l^2*2*self.l/(6*EI_tip);
                    stiffness = F/z_tip;
            end
        end
        
        % Effective mass of the cantilever beam
        % Does not include the effective mass of the base/actuator section
        % f0 is calculated using Rayleight-Ritz there
        % Units: kg
        function effective_mass = effective_mass(self)
            cantilever_effective_mass = 0.243 * self.rho_si * self.w * self.t * self.l;
            effective_mass = cantilever_effective_mass + self.tip_mass;
        end
        
        % Resonant frequency for undamped vibration (first mode).
        % For a simple cantilever, use Bernoulli beam theory.
        % For step/actuator designs, use Rayleigh-Ritz (i.e. U_kin = U_strain)
        % The R-R results agree with 2D and 3D FEA results to <3% when
        % (1 < t_a/t < 6, 1/5 < l_a/l < 5)
        % Units: radians/sec
        function omega_vacuum = omega_vacuum(self)
            omega_bernoulli = sqrt( self.stiffness() / self.effective_mass());
            
            switch self.cantilever_type
                case 'none'
                    omega_vacuum = omega_bernoulli;
                otherwise
                    options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'Display', 'off');
                    options = optimset('TolX', 1e-12, 'TolFun', 1e-12, 'Display', 'off');
                    omega_vacuum = fminbnd(@findEnergyResidual, 1, ...
                        100*omega_bernoulli, options);
            end
            
            % Function for iteratively finding the resonant frequency
            function residual_squared = findEnergyResidual(omega_guess)
                [U_elastic, U_kinetic] = calculateEnergies(omega_guess);
                residual_squared = (U_elastic/U_kinetic - 1)^2;
            end
            
            % Function for computing the elastic and kinetic energy of the beam
            function [U_elastic, U_kinetic] = calculateEnergies(omega)
                
                % Discretize the length of the cantilever
                totalLength = self.l + self.l_a;
                dx = totalLength/(self.numXPoints - 1);
                x = 0:dx:totalLength;
                base_indices = find(x <= self.l_a);
                tip_indices = find(x > self.l_a);
                
                % Define the multilayer mechanics
                EI_base = 1/self.calculateActuatorNormalizedCurvature();
                EI_tip = self.modulus()*self.w*self.t^3/12;
                
                % Empirical correction factors that give better agreement with FEA results
                % Account for cases where t_a >> t, w_a >> w, l_a >> l
                EI_base = (self.t/self.t_a)^.25*EI_base;
                EI_base = (self.w/self.w_a)^.1*EI_base;
                EI_base = (self.l/self.l_a)^.1*EI_base;
                
                % Generate an approximate cantilever deflection profile assuming a
                % point load force at the tip of the beam. Stitch together the two
                % sections (the moment is constant despite the EI discontinuity)
                tip_deflection = 1e-6; % Apply a test force
                F = tip_deflection*3*EI_tip/self.l^3;
                moment = F*(totalLength-x);
                deflection(base_indices) = -F*x(base_indices).^2.* ...
                    (3*totalLength-x(base_indices))/(6*EI_base);
                
                % x-coordinate from the end of l_a
                x_relative = x(tip_indices) - x(max(base_indices));
                
                % Continue with the slope that is at the end of the base section
                if max(base_indices) > 1
                    tip_slope = (deflection(max(base_indices)) - ...
                        deflection(max(base_indices)-1))/dx;
                    % Catch any cases wher l_a = 0 (shouldn't normally happen)
                else
                    tip_slope = 0;
                end
                deflection(tip_indices) = deflection(max(base_indices)) - ...
                    F*x_relative.^2.*(3*self.l-x_relative)/(6*EI_tip) + tip_slope*x_relative;
                
                [E_metal, rho_metal, k_metal, alpha_metal] = ...
                    self.lookup_metal_properties();
                dm_tip = self.w*self.t*self.rho_si;
                switch self.cantilever_type
                    case {'step', 'thermal'}
                        dm_base = self.w_a*(self.t*self.rho_si + ...
                            self.t_oxide*self.rho_sio2 + self.t_a*rho_metal);
                    case 'piezoelectric'
                        dm_base = self.w_a*(self.t*self.rho_si + ...
                            self.t_oxide*self.rho_sio2 + ...
                            self.rho_aln*(self.t_a + self.t_a_seed) + ...
                            rho_metal*(self.t_electrode_bottom + self.t_electrode_top));
                end
                
                % Piecewise kinetic and elastic energies
                Udx_elastic(base_indices) = ...
                    .5*moment(base_indices).^2*dx/EI_base;
                Udx_kinetic(base_indices) = ...
                    .5*(omega*deflection(base_indices)).^2*dx*dm_base;
                Udx_elastic(tip_indices) = ...
                    .5*moment(tip_indices).^2*dx/EI_tip;
                Udx_kinetic(tip_indices) = ...
                    .5*(omega*deflection(tip_indices)).^2*dx*dm_tip;
                
                U_elastic = trapz(x, Udx_elastic);
                U_kinetic = trapz(x, Udx_kinetic);
            end
        end
        
        % Resonant frequency for undamped vibration (first mode)
        % Units: cycles/sec
        function omega_vacuum_hz = omega_vacuum_hz(self)
            omega_vacuum_hz = self.omega_vacuum() / (2*pi);
        end
        
        % Calculate the damped natural frequency and Q (via Brumley/Sader)
        function [omega_damped, Q] = omega_damped_and_Q(self)
            
            % If we're in vacuum, just return the vacuum frequency
            switch self.fluid
                case 'vacuum'
                    omega_damped = self.omega_vacuum();
                    Q = cantilever.maxQ;
                    return;
            end
            
            % Inner function for solving the transcendental eqn to find
            % omega_damped. We're searching for a function minimum, so return
            % the residual squared (continuous and smooth)
            function residual_squared = find_natural_frequency(omega_damped)
                hydro = self.hydrodynamic_function(omega_damped, rho_f, eta_f);
                residual = omega_damped - omega_vacuum*(1 + pi * rho_f * ...
                    self.w/(4 * self.rho_si * self.t) .* real(hydro)).^-0.5;
                residual_squared = residual^2;
            end
            
            % Lookup fluid properties once, then calculate omega_damped and Q
            [rho_f eta_f] = self.lookupFluidProperties();
            omega_vacuum = self.omega_vacuum();
            
            options = optimset('TolX', 1e-12, 'TolFun', 1e-12, 'Display', 'off');
            omega_damped = fminbnd(@find_natural_frequency, 10, omega_vacuum, options);
            
            hydro = self.hydrodynamic_function(omega_damped, rho_f, eta_f);
            Q = (4*self.rho_si*self.t/(pi*rho_f*self.w)+real(hydro))/imag(hydro);
            
            % Catch cases where Q is undefined, too large or too small
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
        
        % Calculate the Reynold's number - note that 'a = w/2' in the Brumley paper
        function reynolds = reynolds(self, omega, rho_f, eta_f)
            reynolds = (rho_f*omega*(self.w/2)^2)/eta_f;
        end
        
        % Calculate hydrodynamic function from the lookup table
        function hydro = hydrodynamic_function(self, omega, rho_f, eta_f)
            
            % A = aspect ratio, Beta = Reynold's number
            A = self.t/self.w;
            Beta = self.reynolds(omega, rho_f, eta_f);
            log_Beta = log10(Beta);
            
            % If the Reynolds number gets too small, the calculation gets upset
            % This only happens during random solution generation usually
            if log_Beta < min(self.Beta_lookup)
                log_Beta = min(self.Beta_lookup);
            end
            
            gamma_real = interp2(cantilever.A_lookup, cantilever.Beta_lookup, ...
                cantilever.gamma_lookup_real, A, log_Beta, 'spline');
            gamma_imag = interp2(cantilever.A_lookup, cantilever.Beta_lookup, ...
                cantilever.gamma_lookup_imag, A, log_Beta, 'spline');
            hydro = complex(gamma_real, gamma_imag);
        end
        
        % ========= Optimization ==========
        
        % Calculate force resolution
        % Units: pN
        function force_resolution = optimize_force_resolution(self, x0)
            self = self.cantilever_from_state(x0);
            force_resolution = self.force_resolution()*1e12;
        end
        
        % Calculate displacement resolution
        % Units: nm
        function displacement_resolution = optimize_displacement_resolution(self, x0)
            self = self.cantilever_from_state(x0);
            displacement_resolution = self.displacement_resolution()*1e9;
        end
        
        % Calculate the force noise density on resonance
        % Units: pN/rtHz
        function force_noise_density = optimize_resonant_force_noise_density(self, x0)
            self = self.cantilever_from_state(x0);
            force_noise_density = self.resonant_force_noise_density()*1e12;
        end
        
        % Calculate the surface stress resolution
        % Units: Pa
        function ss_resolution = optimize_surface_stress_resolution(self, x0)
            self = self.cantilever_from_state(x0);
            ss_resolution = self.surface_stress_resolution()*1e6;
        end
        
        % Used by optimization to bring all state varibles to O(10)
        % If we didn't do that, we'd have O(1e-9) and O(1e20) variables
        function scaling = optimization_scaling(self)
            l_scale = 1e5;
            w_scale = 1e7;
            t_scale = 1e8;
            l_pr_ratio_scale = 1e2;
            v_bridge_scale = 1e1;
            scaling = [l_scale, w_scale, t_scale, l_pr_ratio_scale, ...
                v_bridge_scale, self.doping_optimization_scaling()];
            
            % Actuator specific code
            switch self.cantilever_type
                case 'step'
                    % Do nothing special
                case 'thermal'
                    l_a_scale = 1e6;
                    w_a_scale = 1e6;
                    t_a_scale = 1e9;
                    v_actuator_scale = 1;
                    R_heater_scale = 1e-3;
                    scaling = [scaling l_a_scale w_a_scale t_a_scale ...
                        v_actuator_scale R_heater_scale];
                case 'piezoelectric'
                    l_a_scale = 1e6;
                    w_a_scale = 1e6;
                    t_a_scale = 1e9;
                    v_actuator_scale = 1;
                    scaling = [scaling l_a_scale w_a_scale t_a_scale v_actuator_scale];
            end
        end
        
        % Update the changed optimization parameters
        % All optimization takes place for the same object (i.e. we update 'self')
        % so that things like 'fluid' are maintained
        function self = cantilever_from_state(self, x0)
            scaling = self.optimization_scaling();
            x0 = x0./scaling;
            
            self.l = x0(1);
            self.w = x0(2);
            self.t = x0(3);
            self.l_pr_ratio = x0(4);
            self.v_bridge = x0(5);
            self = self.doping_cantilever_from_state(x0);
            
            % Actuator specific code
            switch self.cantilever_type
                case 'step'
                    % Do nothing special
                case 'thermal'
                    self.l_a = x0(8);
                    self.w_a = x0(9);
                    self.t_a = x0(10);
                    self.v_actuator = x0(11);
                    self.R_heater = x0(12);
                case 'piezoelectric'
                    self.l_a = x0(8);
                    self.w_a = x0(9);
                    self.t_a = x0(10);
                    self.v_actuator = x0(11);
            end
        end
        
        % Return state vector for the current state
        function x = current_state(self)
            x(1) = self.l;
            x(2) = self.w;
            x(3) = self.t;
            x(4) = self.l_pr_ratio;
            x(5) = self.v_bridge;
            x = [x self.doping_current_state()];
            
            % Actuator specific code
            switch self.cantilever_type
                case 'step'
                    % Do nothing special
                case 'thermal'
                    x(8) = self.l_a;
                    x(9) = self.w_a;
                    x(10) = self.t_a;
                    x(11) = self.v_actuator;
                    x(12) = self.R_heater;
                case 'piezoelectric'
                    x(8) = self.l_a;
                    x(9) = self.w_a;
                    x(10) = self.t_a;
                    x(11) = self.v_actuator;
            end
        end
        
        % Set the minimum and maximum bounds for the cantilever state variables.
        % Bounds are written in terms of the initialization variables.
        % Secondary constraints (e.g. power dissipation, resonant frequency)
        % are applied in optimization_constraints()
        function [lb ub] = optimization_bounds(self, parameter_constraints)
            min_l = 10e-6;
            max_l = 3e-3;
            
            min_w = 2e-6;
            max_w = 100e-6;
            
            min_t = 1e-6;
            max_t = 100e-6;
            
            min_l_pr_ratio = 0.01;
            max_l_pr_ratio = 0.99;
            
            min_v_bridge = 0.1;
            max_v_bridge = 10;
            
            [doping_lb doping_ub] = ...
                self.doping_optimization_bounds(parameter_constraints);
            
            % Actuator specific code
            actuator_lb = [];
            actuator_ub = [];
            switch self.cantilever_type
                case {'step', 'none'}
                    % Override the default values if any were provided
                    % constraints = {{'min_l', 'max_v_bridge'}, {5e-6, 10}}
                    if ~isempty(parameter_constraints)
                        keys = parameter_constraints{1};
                        values = parameter_constraints{2};
                        for ii = 1:length(keys)
                            eval([keys{ii} '=' num2str(values{ii}) ';']);
                        end
                    end
                case 'thermal'
                    min_l_a = 5e-6;
                    max_l_a = 200e-6;
                    
                    min_w_a = 2e-6;
                    max_w_a = 50e-6;
                    
                    min_t_a = 200e-9;
                    max_t_a = 3e-6;
                    
                    min_v_actuator = .1;
                    max_v_actuator = 10;
                    
                    min_R_heater = 200;
                    max_R_heater = 5e3;
                    
                    % Override the default values if any were provided
                    % constraints = {{'min_l', 'max_v_bridge'}, {5e-6, 10}}
                    if ~isempty(parameter_constraints)
                        keys = parameter_constraints{1};
                        values = parameter_constraints{2};
                        for ii = 1:length(keys)
                            eval([keys{ii} '=' num2str(values{ii}) ';']);
                        end
                    end
                    
                    actuator_lb = [min_l_a min_w_a min_t_a min_v_actuator min_R_heater];
                    actuator_ub = [max_l_a max_w_a max_t_a max_v_actuator max_R_heater];
                case 'piezoelectric'
                    min_l_a = 5e-6;
                    max_l_a = 200e-6;
                    
                    min_w_a = 2e-6;
                    max_w_a = 30e-6;
                    
                    min_t_a = 200e-9;
                    max_t_a = 3e-6;
                    
                    min_v_actuator = .1;
                    max_v_actuator = 10;
                    
                    % Override the default values if any were provided
                    % constraints is a set of key value pairs, e.g.
                    % constraints = {{'min_l', 'max_v_bridge'}, {5e-6, 10}}
                    if ~isempty(parameter_constraints)
                        keys = parameter_constraints{1};
                        values = parameter_constraints{2};
                        for ii = 1:length(keys)
                            eval([keys{ii} '=' num2str(values{ii}) ';']);
                        end
                    end
                    
                    actuator_lb = [min_l_a min_w_a min_t_a min_v_actuator];
                    actuator_ub = [max_l_a max_w_a max_t_a max_v_actuator];
            end
            
            lb = [min_l, min_w, min_t, min_l_pr_ratio, ...
                min_v_bridge, doping_lb, actuator_lb];
            ub = [max_l, max_w, max_t, max_l_pr_ratio, ...
                max_v_bridge, doping_ub, actuator_ub];
        end
        
        % Generate a random cantilever design (x0) within our param bounds
        function x0 = initial_conditions_random(self, parameter_constraints)
            [lb, ub] = self.optimization_bounds(parameter_constraints);
            
            l_min = lb(1);
            l_max = ub(1);
            w_min = lb(2);
            w_max = ub(2);
            t_min = lb(3);
            t_max = ub(3);
            l_pr_ratio_min = lb(4);
            l_pr_ratio_max = ub(4);
            V_b_min = lb(5);
            V_b_max = ub(5);
            
            % Generate the random values
            l_random = l_min + rand*(l_max - l_min);
            w_random = w_min + rand*(w_max - w_min);
            t_random = t_min + rand*(t_max - t_min);
            l_pr_ratio_random = l_pr_ratio_min + rand*(l_pr_ratio_max - l_pr_ratio_min);
            v_bridge_random = V_b_min + rand*(V_b_max - V_b_min);
            
            x0_doping = self.doping_initial_conditions_random(parameter_constraints);
            
            % Actuator specific code
            x0_actuator = [];
            switch self.cantilever_type
                case 'thermal'
                    l_a_min = lb(8);
                    l_a_max = ub(8);
                    w_a_min = lb(9);
                    w_a_max = ub(9);
                    t_a_min = lb(10);
                    t_a_max = ub(10);
                    v_actuator_min = lb(11);
                    v_actuator_max = ub(11);
                    R_heater_min = lb(12);
                    R_heater_max = ub(12);
                    
                    l_a_random = l_a_min + rand*(l_a_max - l_a_min);
                    w_a_random = w_a_min + rand*(w_a_max - w_a_min);
                    t_a_random = t_a_min + rand*(t_a_max - t_a_min);
                    v_actuator_random = v_actuator_min + ...
                        rand*(v_actuator_max - v_actuator_min);
                    R_heater_random = R_heater_min + rand*(R_heater_max - R_heater_min);
                    x0_actuator = [l_a_random w_a_random t_a_random ...
                        v_actuator_random R_heater_random];
                case 'piezoelectric'
                    l_a_min = lb(8);
                    l_a_max = ub(8);
                    w_a_min = lb(9);
                    w_a_max = ub(9);
                    t_a_min = lb(10);
                    t_a_max = ub(10);
                    v_actuator_min = lb(11);
                    v_actuator_max = ub(11);
                    
                    l_a_random = l_a_min + rand*(l_a_max - l_a_min);
                    w_a_random = w_a_min + rand*(w_a_max - w_a_min);
                    t_a_random = t_a_min + rand*(t_a_max - t_a_min);
                    v_actuator_random = v_actuator_min + ...
                        rand*(v_actuator_max - v_actuator_min);
                    x0_actuator = [l_a_random w_a_random t_a_random v_actuator_random];
            end
            
            x0 = [l_random, w_random, t_random, l_pr_ratio_random, ...
                v_bridge_random, x0_doping, x0_actuator];
        end
        
        % Nonlinear optimization constraints.
        % For a feasible design, all constraints are negative.
        function [C, Ceq] = optimization_constraints(self, x0, nonlinear_constraints)
            c_new = self.cantilever_from_state(x0);
            
            % Default aspect ratios that can be overriden
            min_w_t_ratio = 2;
            min_l_w_ratio = 2;
            min_pr_l_w_ratio = 2;
            min_pr_l = 2e-6;
            
            % Read out the constraints as key-value pairs
            % , e.g. {{'omega_min_hz', 'min_k'}, {1000, 10}}
            if ~isempty(nonlinear_constraints)
                keys = nonlinear_constraints{1};
                values = nonlinear_constraints{2};
                for ii = 1:length(keys)
                    eval([keys{ii} '=' num2str(values{ii}) ';']);
                end
            end
            
            % We start with this single element vector and then append
            % any additional constraints that the user has provided.
            % If a C(ii) > 0 then constraint ii is invalid
            
            % Force resolution must always be positive
            % Use a linear function so that it's smooth/continuous
            resolution = c_new.force_resolution();
            C(1) = -resolution*1e18;
            
            % Resonant frequency
            if exist('omega_min_hz', 'var')
                switch self.fluid
                    case 'vacuum'
                        freq_constraint = omega_min_hz - c_new.omega_vacuum_hz();
                    otherwise
                        [omega_damped_hz, tmp] = c_new.omega_damped_hz_and_Q();
                        freq_constraint = omega_min_hz - omega_damped_hz;
                end
                C = [C freq_constraint];
            end
            
            % Power dissipation
            if exist('max_power', 'var')
                power_constraint = c_new.power_dissipation() - max_power;
                C = [C power_constraint];
            end
            
            % Temp constraints -- approximate lumped model
            if exist('tip_temp', 'var') || exist('max_temp', 'var')
                [TMax, TTip] = c_new.approxTempRise();
            end
            if exist('tip_temp', 'var')
                temp_constraint = TTip - tip_temp;
                C = [C temp_constraint];
            end
            if exist('max_temp', 'var')
                temp_constraint = TMax - max_temp;
                C = [C temp_constraint];
            end
            
            % Temp constraints -- 1D finite differences
            % Slower but good for refining designs
            if exist('tip_temp_exact', 'var') || exist('max_temp_exact', 'var')
                [TMax, TTip] = c_new.calculateMaxAndTipTemp();
            end
            if exist('tip_temp_exact', 'var')
                temp_constraint = TTip - tip_temp_exact;
                C = [C temp_constraint];
            end
            if exist('max_temp_exact', 'var')
                temp_constraint = TMax - max_temp_exact;
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
            
            if exist('max_v_actuator', 'var')
                max_v_actuator_constraint = c_new.v_actuator - max_v_actuator;
                C = [C max_v_actuator_constraint];
            end
            
            if exist('min_tip_deflection', 'var')
                min_tip_deflection_constraint = min_tip_deflection - ...
                    c_new.tipDeflection();
                C = [C min_tip_deflection_constraint];
            end
            
            % Aspect ratio constraints. Default ratios can be changed.
            length_width_ratio = min_l_w_ratio - c_new.l/c_new.w;
            C = [C length_width_ratio];
            
            width_thickness_ratio = min_w_t_ratio - c_new.w/c_new.t;
            C = [C width_thickness_ratio];
            
            pr_length_width_ratio = min_pr_l_w_ratio - c_new.l_pr()/c_new.w_pr();
            C = [C pr_length_width_ratio];
            
            pr_length_constraint = (min_pr_l - c_new.l_pr())*1e6;
            C = [C pr_length_constraint];
            
            
            % Now for equality based constraints
            Ceq = [];
            
            % Fix the stiffness
            if exist('fixed_k', 'var')
                fixed_k_constraint = c_new.stiffness() - fixed_k;
                Ceq = [Ceq fixed_k_constraint];
            end
            
            if exist('fixed_v_bridge', 'var')
                fixed_v_bridge_constraint = c_new.v_bridge - fixed_v_bridge;
                Ceq = [Ceq fixed_v_bridge_constraint];
            end
            
            % Fix the resonant frequency
            if exist('fixed_f0', 'var')
                switch self.fluid
                    case 'vacuum'
                        fixed_f0_constraint = fixed_f0 - c_new.omega_vacuum_hz();
                    otherwise
                        [omega_damped_hz, tmp] = c_new.omega_damped_hz_and_Q();
                        fixed_f0_constraint = fixed_f0 - omega_damped_hz;
                end
                Ceq = [Ceq fixed_f0_constraint];
            end
            
            % Useful lines for debugging constraint failures
            % (e.g. you made l_max too small for it to hit the desired f_0)
            % C
            % Ceq
            % sprintf('Active Index: %d  -- Value: %g\n', find(C==max(C)), max(C))
        end
        
        % The optimization isn't convex so isn't guaranteed to converge.
        % On average it converges 95-99% of the time.
        % For this reason, we generally optimize from random initial guesses
        % and keep going until several have converged to the same value
        function optimized_cantilever = optimize_performance(self, ...
                parameter_constraints, nonlinear_constraints, goal)
            
            percent_match = 0.01; % 1 percent match required between results
            randomize_starting_conditions = 1;
            converged = 0;
            ii = 1;
            resolution = [];
            while ~converged
                % Optimize another cantilever
                [c{ii}, exitflag] = ...
                    self.optimize_performance_once(parameter_constraints, ...
                    nonlinear_constraints, goal, randomize_starting_conditions);
                
                % If the optimization terminated abnormally, skip to the next iteration
                if ~(exitflag == 1 || exitflag == 2)
                    continue
                end
                
                % Record the resolution for the latest cantilever
                if goal == cantilever.goalForceResolution
                    resolution(ii) = c{ii}.force_resolution();
                elseif goal == cantilever.goalDisplacementResolution
                    resolution(ii) = c{ii}.displacement_resolution();
                elseif goal == cantilever.goalForceNoiseDensity
                    resolution(ii) = c{ii}.resonant_force_noise_density();
                elseif goal == cantilever.goalSurfaceStress
                    resolution(ii) = c{ii}.surface_stress_resolution();
                end
                
                % If we have more than one result, consider stopping
                if ii > 1
                    
                    % Sort from smallest to largest, check if the two smallest values agree
                    [temp_resolution, sortIndex] = sort(resolution);
                    fprintf('Resolutions so far: %s\n', mat2str(temp_resolution, 3))
                    resultsAgree = abs(1 - temp_resolution(1)/temp_resolution(2)) ...
                        < percent_match;
                    
                    % If the results agree, then stop the loop. Otherwise, continue
                    if resultsAgree
                        fprintf('CONVERGED. Two best values: %s\n', ...
                            mat2str(temp_resolution(1:2), 3))
                        optimized_cantilever = c{sortIndex(1)};
                        converged = 1;
                    else
                        fprintf('NOT CONVERGED. Two best values: %s\n', ...
                            mat2str(temp_resolution(1:2), 3))
                    end
                end
                
                if ii > cantilever.numOptimizationIterations
                    [temp_resolution, sortIndex] = sort(resolution);
                    optimized_cantilever = c{sortIndex(1)};
                    converged = 1;
                end
                
                ii = ii + 1; % Increment the loop counter
            end
        end
        
        % Optimize, but don't randomize starting point
        function optimized_cantilever = optimize_performance_from_current(self, ...
                parameter_constraints, nonlinear_constraints, goal)
            
            randomize_starting_conditions = 0;
            [optimized_cantilever, tmp] = ...
                self.optimize_performance_once(parameter_constraints, ...
                nonlinear_constraints, goal, randomize_starting_conditions);
        end
        
        function [optimized_cantilever, exitflag] = optimize_performance_once(self, ...
                parameter_constraints, nonlinear_constraints, ...
                goal, randomize_starting_conditions)
            
            scaling = self.optimization_scaling();
            self.check_valid_cantilever();
            
            % If random_flag = 1, start from random conditions
            % Else, start from the current cantilever state vector
            if randomize_starting_conditions == 1
                problem.x0 = ...
                    scaling.*self.initial_conditions_random(parameter_constraints);
            else
                problem.x0 = scaling.*self.current_state();
            end
            
            if goal == cantilever.goalForceResolution
                problem.objective = @self.optimize_force_resolution;
            elseif goal == cantilever.goalDisplacementResolution
                problem.objective = @self.optimize_displacement_resolution;
            elseif goal == cantilever.goalForceNoiseDensity
                problem.objective = @self.optimize_resonant_force_noise_density;
            elseif goal == cantilever.goalSurfaceStress
                problem.objective = @self.optimize_surface_stress_resolution;
            end
            
            [lb ub] = self.optimization_bounds(parameter_constraints);
            problem.lb = scaling.*lb;
            problem.ub = scaling.*ub;
            
            problem.options.TolFun = 1e-9;
            problem.options.TolCon = 1e-9;
            problem.options.TolX = 1e-9;
            
            problem.options.MaxFunEvals = 2000;
            problem.options.MaxIter = 2000;
            problem.options.Display = 'iter-detailed'; % iter, final, none, off
            problem.options.UseParallel = 'always'; % For multicore processors
            problem.options.TypicalX = ones(1, length(scaling));
            problem.options.InitBarrierParam = 10;
            
            % interior-point does not do well if the design is already
            % close to the optimal state (i.e. warm start). sqp is
            % supposed to do better in this condition. try experimenting
            % with both, particular for implanted devices.
            problem.options.Algorithm = 'sqp'; % sqp, interior-point
            problem.solver = 'fmincon';
            
            problem.nonlcon = @(x) self.optimization_constraints(x, ...
                nonlinear_constraints);
            
            % These errors come up occasionally for ion implanted devices
            warning off MATLAB:singularMatrix
            warning off MATLAB:nearlySingularMatrix
            
            [x, tmp, exitflag] = fmincon(problem);
            optimized_cantilever = self.cantilever_from_state(x);
        end
    end
end
