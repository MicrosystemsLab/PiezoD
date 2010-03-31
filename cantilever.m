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

    % Note that w is the total cantilever width where the piezoresistor is.
    % If the piezoresistor has two legs, then w_pr_ratio can only have a
    % maximum value of 0.5.
    properties
        l; % overall cantilever length
        w; % overall cantilever width (total width of both legs)
        t; % overall cantilever thickness

        l_pr_ratio; % piezoresistor length ratio

        v_bridge; % Volts
        number_of_piezoresistors = 2;
        
        rho_cantilever = 2330;
        rho_fluid = 1e3;
        eta_fluid = 1e-3;

        freq_min; % Hertz
        freq_max; % Hertz

        alpha = 1e-5; % unitless
        T = 300; % kelvin
        k_b = 1.38e-23; % J/K
        k_b_eV = 8.617343e-5; % eV/K
        q = 1.60218e-19; % Coulombs

        doping_type = 'boron'; % default value = boron
    end
    
    % Can be referred to with cantilever.variableName
    properties (Constant)
        colorBlue = [11/255 132/255 199/255];
        colorOrange = [255/255 102/255 0/255];
        colorGreen = [0/255 153/255 0/255];
        markerSize = 8;
        colorBlack = [0 0 0];   
        
        goalForceResolution = 0;
        goalDisplacementResolution = 1;
        
        fluidVacuum = 0;
        fluidWater = 1;
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

        function self = cantilever(freq_min, freq_max, l, w, t, l_pr_ratio, ...
                v_bridge, doping_type)

            self.freq_min = freq_min;
            self.freq_max = freq_max;

            self.l = l;
            self.w = w;
            self.t = t;

            self.l_pr_ratio = l_pr_ratio;

            self.v_bridge = v_bridge;

            self.freq_max = freq_max;
            self.freq_min = freq_min;

            self.doping_type = doping_type;
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
            fprintf('Cantilever L/W/T: %f %f %f \n', self.l*1e6, self.w*1e6, self.t*1e6)
            fprintf('PR L/W: %f %f %f \n', self.l_pr()*1e6, self.w_pr()*1e6)
            fprintf('PR Length Ratio: %g \n', self.l_pr_ratio)
            fprintf('Wheatstone bridge bias voltage: %f \n', self.v_bridge)
            fprintf('Number of silicon resistors: %f \n', self.number_of_piezoresistors)
            fprintf('\n')
            fprintf('Resistance: %f \n', self.resistance())
            fprintf('Power dissipation (mW): %g \n', self.power_dissipation()*1e3)
            fprintf('\n')
            fprintf('Force resolution (N): %g \n', self.force_resolution())
            fprintf('Displacement resolution (m): %g \n', self.displacement_resolution())
            fprintf('Sensitivity (V/N) %g \n', self.force_sensitivity())
            fprintf('Beta %g \n', self.beta())
            fprintf('\n')
            fprintf('Integrated noise (V): %g \n', self.integrated_noise())
            fprintf('Integrated johnson noise (V): %g \n', self.integrated_johnson_noise())
            fprintf('Integrated 1/f noise (V): %g \n', self.integrated_hooge_noise())
            fprintf('Amplifier noise (V): %g \n', self.integrated_amplifier_noise())
            fprintf('Knee frequency (Hz): %g \n', self.knee_frequency())
            fprintf('Johnson/Hooge: %g \n', self.integrated_johnson_noise()/self.integrated_hooge_noise())
            fprintf('\n')
            fprintf('Sheet Resistance: %g \n', self.sheet_resistance())

            [x, doping] = self.doping_profile();
            Nz = trapz(x, doping*1e6);
            fprintf('Nz: %g \n', Nz)

            fprintf('\n')
            fprintf('Number of Carriers: %g \n', self.number_of_carriers());
            fprintf('Stiffness (N/m): %g \n', self.stiffness())
            fprintf('Damped freq: %f \n', self.omega_damped_hz())
            fprintf('Quality factor: %f \n', self.quality_factor())
            fprintf('Vacuum freq: %f \n', self.omega_vacuum_hz())
            fprintf('Freq range: %f to %f \n', self.freq_min, self.freq_max)
            fprintf('\n')
        end

        function print_performance_for_excel(self)
            variables_to_print = [self.l*1e6, self.w*1e6, self.t*1e6, ...
                self.l_pr()*1e6, self.l_pr_ratio, ...
                self.v_bridge, self.freq_min, self.freq_max, ...
                self.force_resolution(), self.displacement_resolution(), ...
                self.omega_vacuum_hz(), self.omega_damped_hz(), ...
                self.stiffness(), self.quality_factor(), self.force_sensitivity(), self.beta(), ...
                self.resistance(), self.power_dissipation()*1e3, ...
                self.integrated_noise(), self.integrated_johnson_noise(), ...
                self.integrated_hooge_noise(), self.knee_frequency()];

            for print_index = 1:length(variables_to_print)
                fprintf('%g \t', variables_to_print(print_index));
            end
        end

        % ==================================
        % ===== Calculate resistance =======
        % ==================================

        % Calculate total resistance of piezoresistor. Includes effect of
        % other resistances (gamma)
        % Units: ohms
        function resistance = resistance(self)
            number_of_squares = self.resistor_length()/self.w_pr();
            resistance = number_of_squares * self.sheet_resistance() / self.gamma();
        end

        % Calculate resistor length, used to calculat resistance
        % and number of carriers
        % Units: m
        function resistor_length = resistor_length(self)
            resistor_length = 2*self.l_pr();
        end

        % Calculate sheet resistance
        % Depends upon abstract method self.doping_profile()
        function Rs = sheet_resistance(self)
            [x, doping] = self.doping_profile(); % x -> m, doping -> N/cm^3
            conductivity = self.conductivity(doping); % ohm-cm
            Rs = 1/trapz(x*1e2, conductivity); % convert x to cm
        end

        function sigma = conductivity(self, dopant_concentration)
            mu = self.mobility(dopant_concentration); % cm^2/V-sec
            sigma = mu.*self.q.*dopant_concentration; % C/V-sec-cm
        end

        % Data from "Modeling of Carrier Mobility Against Carrier Concentration
        % in Arsenic-,  Phosphorus-, and Boron-Doped  Silicon"
        % Masetti, Serveri and Solmi - IEEE Trans. on Electron Devices, (1983)
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

        % The number of current carriers in the piezoresistor
        % Calculated by integrating the number of carriers to the junction
        % depth and multiplaying by the planform dimensions of the piezoresistor
        % Units: unitless
        function number_of_carriers = number_of_carriers(self)
            [x, doping] = self.doping_profile(); % Units: x -> m, doping -> N/cm^3
            Nz = trapz(x, doping*1e6); % doping: N/cm^3 -> N/m^3
            number_of_carriers = Nz*self.resistor_length()*self.w_pr();
        end

        % 1/f noise density from the entire Wheatstone bridge
        % Assuming two piezoresistors (with 1/f noise) and two w/o 1/f
        % noise, so the 1/f noise density is sqrt(2) larger
        % Units: V
        function hooge_noise_density = hooge_noise_density(self)
            hooge_noise_density = sqrt(self.number_of_piezoresistors)/2*sqrt(self.alpha*self.v_bridge^2/self.number_of_carriers());
        end

        % Integrated 1/f noise density for the entire Wheatstone bridge
        % Unit: V
        function integrated_hooge_noise = integrated_hooge_noise(self)
            integrated_hooge_noise = self.hooge_noise_density() * sqrt(log(self.freq_max / self.freq_min));
        end

        % Johnson noise density from the entire Wheatstone bridge
        % Assuming four resistors each equal to R, so the Johnson noise
        % density is the same as for a single resistor R
        % Units: V/rtHz
        function johnson_noise_density = johnson_noise_density(self)
            johnson_noise_density = sqrt(4 * self.k_b * self.T * self.resistance());
        end

        % Integrated Johnson noise from the entire Wheatstone bridge
        % Unit: V
        function integrated_johnson_noise = integrated_johnson_noise(self)
            integrated_johnson_noise = self.johnson_noise_density()*sqrt(self.freq_max - self.freq_min);
        end

        % Assume INA103 type performance
        function integrated_amplifier_noise = integrated_amplifier_noise(self)
            
            % Additional sqrt(2) prefactor because there are two amplifier inputs
            Vamp_J = 1.2e-9*sqrt(2)*sqrt(self.freq_max - self.freq_min); % 1.8 nV/rtHz
            Vamp_H = 4e-9 *sqrt(2)*sqrt(log(self.freq_max/self.freq_min)); % 10 nV/rtHz @ 1Hz, 3.3 (sqrt of 10) nV/rtHz @ 10Hz
            Iamp_J = 2e-12 *sqrt(2)*sqrt(self.freq_max - self.freq_min);  % 2pA/rtHz
            Iamp_H = 25e-12*sqrt(2)*sqrt(log(self.freq_max/self.freq_min)); % 25 pA/rtHz @ 1Hz
            
            R_effective = self.resistance()/2; % resistance seen by amplifier inputs
            
            integrated_amplifier_noise = sqrt(Vamp_J^2 + Vamp_H^2 + (R_effective*Iamp_J)^2 + (R_effective*Iamp_H)^2);
        end

        % Calculate the knee frequency
        % Equating 1/f noise and johnson... sqrt(alpha*V_bridge^2/(2*N*f_knee)) = sqrt(4*kb*T*R)
        % Leads to f_knee = alpha*V_bridge^2/(2*N*S_j^2)
        function knee_frequency = knee_frequency(self)
            knee_frequency = self.alpha*self.v_bridge^2/(2*self.number_of_carriers()*self.johnson_noise_density()^2);
        end

        % Integrated cantilever noise for given bandwidth
        % Source are uncorrelated so they add in RMS
        % Units: V = sqrt(V^2 + V^2 + ..)
        function integrated_noise = integrated_noise(self)
            integrated_noise = sqrt(self.integrated_johnson_noise()^2 + self.integrated_hooge_noise()^2 + self.integrated_amplifier_noise()^2);
        end

        % Calculate the noise in V/rtHz at a given frequency
        function voltage_noise = voltage_noise(self, frequency)
            hooge_noise = self.hooge_noise_density()./sqrt(frequency); % V/rtHz
            johnson_noise = self.johnson_noise_density(); % V/rtHz
            voltage_noise = sqrt(hooge_noise.^2 + johnson_noise^2 );
        end

        function plot_noise_spectrum(self)
            frequency = logspace(log10(self.freq_min), log10(self.freq_max), 1e4);
            noise = self.voltage_noise(frequency);
            axis equal;

            plot(frequency, noise, 'LineWidth', 2);
            set(gca, 'xscale','log', 'yscale','log');
            set(gca, 'LineWidth', 1.5, 'FontSize', 14);
            ylabel('Noise Voltage Spectral Density (V/rtHz)', 'FontSize', 16);
            xlabel('Frequency (Hz)', 'FontSize', 16);
            
            
        end
        
        function f_min_cumulative = f_min_cumulative(self)
            frequency = logspace(log10(self.freq_min), log10(self.freq_max), 1e4);
            noise = self.voltage_noise(frequency)
            sensitivity = self.cPR.force_sensitivity();
            force_noise_density = noise./sensitivity;
            f_min_cumulative = sqrt(cumtrapz(frequency, force_noise_density.^2));
        end


        % ==================================
        % ===== Calculate sensitivity ======
        % ==================================


        % Piezoresistance factor
        % Following Jonah Harley's calculation based upon experimental data
        function piezoresistance_factor = piezoresistance_factor(self, dopant_concentration)
            b = 1.53e22;
            a = 0.2014;
            piezoresistance_factor = log10((b./dopant_concentration) .^ a);
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

        % Reduction in sensitivity from piezoresistor not located just at the
        % surface. Calculated for the general case of an arbitrarily shaped
        % doping profile. Taken from Sung-Jin Park's Hilton Head 2008 paper.
        % Units: None
        function beta = beta(self)
            [x, doping_concentration] = self.doping_profile();

            % x is supposed to vary from t/2 to -t/2
            x = (self.t/2 - x)*1e2; % x: m -> cm

            mu = self.mobility(doping_concentration); % cm^2/V-s
            P = self.piezoresistance_factor(doping_concentration);

            numerator = trapz(x, self.q.*mu.*doping_concentration.*P.*x);
            denominator = trapz(x, self.q.*mu.*doping_concentration);
            beta = 2*numerator/(self.t*1e2*denominator); % t: m -> cm
        end

        % Ratio of piezoresistor resistance to total resistance (< 1)
        % Assume that metal interconnects are used and that their resistance is about 10% of the total
        function gamma = gamma(self)
            gamma = 0.9;
        end

        % Units: V/N
        function force_sensitivity = force_sensitivity(self)
            v_bias = self.v_bridge/2;
            longitudinal_sensitivity = 3*(self.l - self.l_pr()/2)*self.max_piezoresistance_factor()/(2*self.w*self.t^2)*self.beta()*self.gamma()*v_bias;

            % TODO: Reduce the sensitivity due to transverse current flow at the tip of the cantilever.
            % This requires that the model is extended to include the gap between the legs
            % Otherwise, if we're assuming an infinitesimal gap, then there is no transverse reduction
            % And when the transverse length has a small resistance relative to the longitudinal resistance, the effect
            % will be negligible.
            
            force_sensitivity = longitudinal_sensitivity;
        end


        % ==================================
        % == Calculate power dissipation ===
        % ==================================

        % Power dissipation (W) in the cantilever
        function power_dissipation = power_dissipation(self)
            power_dissipation = self.v_bridge^2 / (4*self.resistance());
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

        % Calculate elastic modulus based upon dopant type, assuming that
        % we're orienting the piezoresistor in the direction of maximum
        % sensitivity
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

        % Resonant frequency for undamped vibration (first mode)
        % Units: radians/sec
        function omega_vacuum = omega_vacuum(self)
            m_effective = 0.243 * self.rho_cantilever * self.w * self.t * self.l;
            omega_vacuum = sqrt( self.stiffness() / m_effective);
        end

        % Resonant frequency for undamped vibration (first mode)
        % Units: cycles/sec
        function omega_vacuum_hz = omega_vacuum_hz(self)
            omega_vacuum_hz = self.omega_vacuum() / (2*pi);
        end

        % Calculate the damped natural frequency, which we know lies
        % between zero and the natural frequency in vacuum. We use fminbnd
        % to find the frequency at which the omega_damped =
        % omega_vacuum*...
        function omega_damped = omega_damped(self)
            function residual = find_natural_frequency(omega_damped)
                hydro = self.hydrodynamic_function(omega_damped);
                residual = abs(omega_damped - self.omega_vacuum()*(1 + pi * self.rho_fluid * self.w/(4 * self.rho_cantilever * self.t) .* real(hydro)).^-0.5);
            end

            omega_damped = fminbnd(@find_natural_frequency, 0, self.omega_vacuum());
        end

        function omega_damped_hz = omega_damped_hz(self)
            omega_damped_hz = self.omega_damped() / (2*pi);
        end


        % Note that this is a combination of omega_damped_hz() and
        % quality_factor(), but is 2x faster than evaluating both
        % separately.
        %
        % Care must be taken to be sure that it's kept up to date (if
        % quality_factor() changes).
        function [omega_damped_hz, quality_factor] = omega_damped_hz_quality_factor(self)
            omega_damped = self.omega_damped();
            omega_damped_hz = omega_damped / (2*pi);

            hydro = self.hydrodynamic_function(omega_damped);
            quality_factor = (4 * self.rho_cantilever * self.t / (pi * self.rho_fluid * self.w) + real(hydro)) / imag(hydro);
        end

        % Holding the cantilever width and thickness constant, figure out
        % the maximum length to meet a given natural frequency requirement
        function max_length = max_length_for_given_natural_frequency(self)

            function residual = find_max_length(length)
                self.l = length;

                residual = abs(self.omega_damped_hz() - self.freq_max);
            end

            options = optimset('TolX', 1e-9);
            min_estimated_length = 0;
            max_estimated_length = 1e-3; % If you're doing very low frequency stuff, you might need to increase this
            max_length = fminbnd(@find_max_length, min_estimated_length, max_estimated_length, options);
        end

        % Calculate the quality factor
        % TODO: Only valid for Q >> 1, should be improved from Sader papers in
        % the future
        function quality_factor = quality_factor(self)
            hydro = self.hydrodynamic_function(self.omega_damped());
            quality_factor = (4 * self.rho_cantilever * self.t / (pi * self.rho_fluid * self.w) + real(hydro)) / imag(hydro);
        end

        function reynolds = reynolds(self, omega)
            reynolds = (self.rho_fluid* omega * self.w^2) / self.eta_fluid;
        end

        % Calculate kappa for the first mode (change 1 to other values if
        % interested in higher order modes)
        function kappa = kappa(self)
            C_n(1) = fzero(@(x) 1+cos(x)*cosh(x),1);
            C_n(2) = fzero(@(x) 1+cos(x)*cosh(x),5);
            C_n(3) = fzero(@(x) 1+cos(x)*cosh(x),7);
            C_n(4) = fzero(@(x) 1+cos(x)*cosh(x),10);
            C = C_n(1); % Right now only for 1st order modes
            kappa = C * self.w / self.l;
        end

        % For the inviscid case (Re >> 1)
        % From "Resonant frequencies of a rectangular cantilever beam
        % immersed in a fluid", Sader (2006)
        function hydro = hydrodynamic_function_inviscid(self)
            kappa = self.kappa();
            hydro = (1 + 0.74273*kappa + 0.14862*kappa^2)/(1 + 0.74273*kappa + 0.35004*kappa^2 + 0.058364*kappa^3);
        end

        % Perform force and displacement resolution optimizations with the
        % spring constant fixed.
        function resolution_tradeoff_plot(self, parameter_constraints, nonlinear_constraints, k_range)

            n_k = 10;
            goal = 0; % Irrelevant for this calculation
            k_values = logspace( log10(min(k_range)), log10(max(k_range)), n_k);
            
            for ii = 1:length(k_values)
                k = k_values(ii);

                % Add the fixed k constraint
                nonlinear_constraints_k = { cat(2, nonlinear_constraints{1}, 'fixed_k'), cat(2, nonlinear_constraints{2}, k)};
                
                self = self.optimize_performance_from_current(parameter_constraints, nonlinear_constraints_k, goal);
                
                force_resolution(ii) = self.force_resolution();
                displacement_resolution(ii) = self.displacement_resolution();
            end
           
            figure
            hold on
            plot(k_values, force_resolution, '-o' ,'LineWidth', 2, ...
                'MarkerSize', cantilever.markerSize, 'Color', cantilever.colorBlue, ...
                'MarkerFaceColor', cantilever.colorBlue, 'MarkerEdgeColor', cantilever.colorBlack);
            plot(k_values, displacement_resolution, '-o' ,'LineWidth', 2, ...
                'MarkerSize', cantilever.markerSize, 'Color', cantilever.colorOrange, ...
                'MarkerFaceColor', cantilever.colorOrange, 'MarkerEdgeColor', cantilever.colorBlack);
            hold off
            set(gca, 'yscale', 'log', 'xscale', 'log', 'LineWidth', 1.5, 'FontSize', 14);
            
            xlabel('Spring Constant (N/m)', 'FontSize', 16);
            ylabel('Resolution (N or m)', 'FontSize', 16);
            legend('Force resolution', 'Displacement Resolution');
            
            print('-dpng','-r300', 'ResolutionTradeoff')
            
        end

        % Calculate hydrodynamic function
        function hydro = hydrodynamic_function(self, omega)
            % Calculate Re and kappa
            kappa = self.kappa();
            reynolds = self.reynolds(omega);
            log_reynolds = log10(reynolds);

            % Built our interpolation table
            kappa_lookup = [0 0.125 0.25 0.5 0.75 1 2 3 5 7 10 20];
            reynolds_lookup = [-4 -3.5 -3 -2.5 -2 -1.5 -1 -.5 0 0.5 1 1.5 2 2.5 3 3.5 4];

            % Hydrodynamic function lookup table (ref here)
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
                0.332285  0.329276    0.32283   0.305262  0.285953  0.26763   0.216732  0.194409  0.186218  0.195634  0.221631  0.349855;
                0.191043  0.189434    0.185931  0.176166  0.165118  0.154323  0.122124  0.105573  0.0938839 0.0925686 0.09682   0.126835;
                0.112082  0.111181    0.109199  0.103595  0.0971392 0.0907188 0.0707736 0.059728  0.0505049 0.0476557 0.0471326 0.0534759;
                0.0665172 0.0659974   0.0648471 0.0615627 0.0577366 0.0538889 0.0416384 0.0345727 0.0282418 0.025856  0.024611  0.0252877];

            tau_real = interp2(kappa_lookup, reynolds_lookup, tau_lookup_real, kappa, log_reynolds,'cubic');
            tau_imag = interp2(kappa_lookup, reynolds_lookup, tau_lookup_imag, kappa, log_reynolds,'cubic');

            hydro = complex(tau_real, tau_imag);

        end

        % ==================================
        % ========= Optimization  ==========
        % ==================================

        % Calculate force resolution (goal) from state variable vector
        function force_resolution = optimize_force_resolution(self, x0)
            self = self.cantilever_from_state(x0);
            force_resolution = self.force_resolution()*1e12;
        end
        
        function displacement_resolution = optimize_displacement_resolution(self, x0)
            self = self.cantilever_from_state(x0);
            displacement_resolution = self.displacement_resolution()*1e9;
        end


        % Nonlinear optimization constraints
        % All constraint components (e.g. C(1)) are required to be negative
        function [C, Ceq] = optimization_constraints(self, x0, nonlinear_constraints)

            c_new = self.cantilever_from_state(x0);

            % Read out the constraints as key-value pairs, e.g.
            % {{'omega_min_hz', 'min_k'}, {1000, 10}}
            if ~isempty(nonlinear_constraints)
                keys = nonlinear_constraints{1};
                values = nonlinear_constraints{2};
                for ii = 1:length(keys)
                    eval([keys{ii} '=' num2str(values{ii}) ';']);
                end
            end

            % Force resolution must always be positive
            % We start with this single element vector and then append any
            % additional constraints that the user has provided.
            C(1) = -c_new.force_resolution();

            % Natural frequency
            if exist('omega_min_hz')
                switch fluid_type
                    case cantilever.fluidVacuum
                        freq_constraint = omega_min_hz - c_new.omega_vacuum_hz();
                    case cantilever.fluidWater
                        freq_constraint = omega_min_hz - c_new.omega_damped_hz();
                end
                C = [C freq_constraint];
            end

            % Power dissipation
            if exist('max_power')
                power_constraint = c_new.power_dissipation() - max_power;
                C = [C power_constraint];
            end

            if exist('min_k')
                min_k_constraint = min_k - c_new.stiffness();
                C = [C min_k_constraint];
            end

            if exist('max_k')
                max_k_constraint = c_new.stiffness() - max_k;
                C = [C max_k_constraint];
            end

            % Silicon: L/W > 10, W/T > 3
            length_width_ratio = 5 - c_new.l/c_new.w;
            width_thickness_ratio = 2 - c_new.w/c_new.t;
            C = [C length_width_ratio width_thickness_ratio];

            % PR: L/W > 5
            pr_length_width_ratio = 2 - c_new.l_pr()/c_new.w_pr();
            C = [C pr_length_width_ratio];
            
            % Now for equality based constraints
            Ceq = [];

            if exist('fixed_k')
                fixed_k_constraint = c_new.stiffness() - fixed_k;
                Ceq = [Ceq fixed_k_constraint];
            end

            % Misc constraints
            %             C(4) = 3e-6 - c_new.w*c_new.w_gap_ratio; % Air gap needs to be at least 3 microns wide

            % Optional equality constraints
            %             Ceq(1) = omega_min_hz - c_new.omega_vacuum_hz(); % for lock-in
            %             Ceq(1) = (t - t*t_pr_ratio) - .3e-6; % Make sure that the cantilever substrate starts out a X microns thick
            %             Ceq(2) = doping - 4e19;
            %             Ceq(3) = t*t_pr_ratio - .3e-6;
        end

        % This cantilever problem isn't guaranteed to converge, and in
        % practice it fails to converge about 1% of the time for random
        % initial conditions. For this reason, it is best to start from a
        % random initial seed and perform the optimization and checking to
        % make sure that it converges repeatedly.
        function optimized_cantilever = optimize_performance(self, parameter_constraints, nonlinear_constraints, goal)

            n = 3; % the number of trials we want to try
            percent_match = 0.01;
            randomize_starting_conditions = 1;
            
            for ii = 1:n
                c{ii} = self.optimize_performance_once(parameter_constraints, nonlinear_constraints, goal, randomize_starting_conditions);
                
                if goal == cantilever.goalForceResolution
                    resolution(ii) = c{ii}.force_resolution();
                elseif goal == cantilever.goalDisplacementResolution
                    resolution(ii) = c{ii}.displacement_resolution();
                end
            end

            best_index = find(resolution == min(resolution));
            optimized_cantilever = c{best_index};

            % Output the results
            min_resolution = min(resolution);
            max_resolution = max(resolution);

            if (1 - min_resolution/max_resolution) > percent_match
                fprintf(['Optimization did not converge at least once. Values = ' num2str(resolution) '\n'])
            end
            fprintf(['Optimization converged. Values = ' num2str(resolution) '\n'])
        end

        % Optimize, but don't randomize starting point
        function optimized_cantilever = optimize_performance_from_current(self, parameter_constraints, nonlinear_constraints, goal)
            randomize_starting_conditions = 0;
            optimized_cantilever = self.optimize_performance_once(parameter_constraints, nonlinear_constraints, goal, randomize_starting_conditions);
        end

        function optimized_cantilever = optimize_performance_once(self, parameter_constraints, nonlinear_constraints, goal, randomize_starting_conditions)

            scaling = self.optimization_scaling();
           
            if goal == cantilever.goalForceResolution
                problem.objective = @self.optimize_force_resolution;
            elseif goal == cantilever.goalDisplacementResolution
                problem.objective = @self.optimize_displacement_resolution;
            end

            % If random_flag = 1, start from random conditions. Otherwise
            % start from the current cantilever state vector
            if randomize_starting_conditions == 1
                problem.x0 = scaling.*self.initial_conditions_random(parameter_constraints);
            else
                problem.x0 = scaling.*self.current_state();
            end

            [lb ub] = self.optimization_bounds(parameter_constraints);
            problem.lb = scaling.*lb;
            problem.ub = scaling.*ub;

            problem.options.TolFun = 1e-12;
            problem.options.TolCon = 1e-12;
            problem.options.TolX = 1e-12;

            problem.options.MaxFunEvals = 5000;
            problem.options.MaxIter = 1000;
            problem.options.Display = 'iter';

            problem.options.UseParallel = 'always'; % For multicore processors
            
            problem.options.Algorithm = 'Interior-point';
            problem.solver = 'fmincon';
            
            problem.nonlcon = @(x) self.optimization_constraints(x, nonlinear_constraints);

            x = fmincon(problem);
            optimized_cantilever = self.cantilever_from_state(x);
        end
    end
end
