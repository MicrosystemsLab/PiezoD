classdef cantilever_poly

    % Note that w is the total cantilever width where the piezoresistor is.
    % If the piezoresistor has two legs, then w_pr_ratio can only have a
    % maximum value of 0.5.
    properties
        l; % overall cantilever length
        w; % overall cantilever width (total width of both legs)
        
        t_top;
        t_mid;
        t_bot;

        matl_top;
        matl_mid;
        matl_bot;
        
        l_pr_ratio; % piezoresistor length ratio
        v_bridge; % Volts
        dopant_concentration;
        
        number_of_piezoresistors_on_cantilever = 2;
        number_of_piezoresistors = 2;
        
        rms_actuator_displacement_noise = 1e-9; % m
        
        rho_si = 2330;
        rho_oxide = 2634;
        rho_poly = 2330;
        rho_al = 2700;
        rho_ti = 4506;
        
        E_si = 169e9;
        E_poly = 150e9;
        E_oxide = 70e9;
        E_al = 70e9;
        E_ti = 85e9;
        
        freq_min; % Hertz
        freq_max; % Hertz

        alpha = 1e-3; % roughly independent of material type for polycrystalline films
    end
    
    % Can be referred to with cantilever.variableName
    properties (Constant)
        T = 300; % kelvin
        k_b = 1.38e-23; % J/K
        k_b_eV = 8.617343e-5; % eV/K
        q = 1.60218e-19; % Coulombs
        numFrequencyPoints = 1e3;
        
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

    methods

        function self = thin_film_cantilever(freq_min, freq_max, l, w, l_pr_ratio, ...
                v_bridge, t_top, t_mid, t_bot, matl_top, matl_mid, matl_bot, dopant_concentration)

            self.freq_min = freq_min;
            self.freq_max = freq_max;

            self.l = l;
            self.w = w;
            self.l_pr_ratio = l_pr_ratio;
            self.v_bridge = v_bridge;
            
            self.t_top = t_top;
            self.t_mid = t_mid;
            self.t_bot = t_bot;
            
            % Enforce symmetry if we're using a double PR configuration
            if self.number_of_piezoresistors_on_cantilever == 2
                self.t_bot = self.t_top;
            end
            
            self.matl_top = matl_top;
            self.matl_mid = matl_mid;
            self.matl_bot = matl_bot;
            
            self.dopant_concentration = dopant_concentration; % only used for poly
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
            fprintf('Cantilever L/W: %f %f\n', self.l*1e6, self.w*1e6)
            fprintf('Layers - Bottom:%s Mid:%s Top:%s\n', self.matl_bot, self.matl_mid, self.matl_top)
            fprintf('Thicnkesses (nm) - Bottom:%f Mid:%f Top:%f\n', self.t_bot*1e9, self.t_mid*1e9, self.t_top*1e9)
            fprintf('PR Length Ratio: %g \n', self.l_pr_ratio)
            fprintf('Wheatstone bridge bias voltage: %f \n', self.v_bridge)
            fprintf('Number of piezoresistors: %f \n', self.number_of_piezoresistors)
            fprintf('\n')
            fprintf('Resistance per resistor: %f \n', self.resistance())
            fprintf('Total cantilever power dissipation (mW): %g \n', self.power_dissipation()*1e3)
            fprintf('\n')
            fprintf('Force resolution (N): %g \n', self.force_resolution())
            fprintf('Displacement resolution (m): %g \n', self.displacement_resolution())
            fprintf('Sensitivity (V/N) %g \n', self.force_sensitivity())
            fprintf('\n')
            fprintf('Integrated noise (V): %g \n', self.integrated_noise())
            fprintf('Integrated johnson noise (V): %g \n', self.johnson_integrated())
            fprintf('Integrated 1/f noise (V): %g \n', self.hooge_integrated())
            fprintf('Amplifier noise (V): %g \n', self.amplifier_integrated())
            fprintf('Knee frequency (Hz): %g \n', self.knee_frequency())
            fprintf('Johnson/Hooge: %g \n', self.johnson_integrated()/self.hooge_integrated())
            fprintf('\n')
            fprintf('Sheet Resistance: %g \n', self.sheet_resistance())
            fprintf('\n')
            fprintf('Number of Carriers: %g \n', self.number_of_carriers());
            fprintf('Stiffness (N/m): %g \n', self.stiffness())
            fprintf('Vacuum freq (Hz): %f \n', self.resonant_frequency())
            fprintf('Freq range: %f to %f \n', self.freq_min, self.freq_max)
            fprintf('\n')
        end

        function print_performance_for_excel(self)
            fprintf('%s\t%s\t%s\t', self.matl_top, self.matl_mid, self.matl_bot)
            variables_to_print = [self.freq_min, self.freq_max*1e-3, ...
                self.l*1e6, self.w*1e6, self.l_pr()*1e6, self.dopant_concentration, ...
                self.t_bot*1e9, self.t_mid*1e9, self.t_top*1e9, ...
                self.v_bridge, self.resistance()*1e-3, self.sheet_resistance(), self.power_dissipation()*1e3, ...
                self.number_of_piezoresistors, self.number_of_piezoresistors_on_cantilever, ...
                self.force_resolution()*1e12, self.displacement_resolution()*1e9, ...
                self.resonant_frequency()*1e-3, ...
                self.force_sensitivity(), self.stiffness()*1e3, ...
                self.integrated_noise()*1e6, self.johnson_integrated()*1e6, ...
                self.hooge_integrated()*1e6, self.amplifier_integrated()*1e6, ...
                self.knee_frequency(), self.number_of_carriers()];            
            
            for print_index = 1:length(variables_to_print)
               fprintf('%4g\t', variables_to_print(print_index)); 
            end
            fprintf('\n');
        end
        
        % ==========================================================================
        % ============================ Beam Mechanics ==============================
        % ==========================================================================
        
        function Zm = neutral_axis(self)
            z = [self.t_bot/2 self.t_bot+self.t_mid/2 self.t_bot+self.t_mid+self.t_top/2];

            matls = {self.matl_bot self.matl_mid self.matl_top};
            for ii = 1:length(matls)
                switch matls{ii}
                    case 'si'
                        E(ii) = self.E_si;
                    case 'poly'
                        E(ii) = self.E_poly;
                    case 'ti'
                        E(ii) = self.E_ti;
                    case 'al'
                        E(ii) = self.E_al;
                    case 'oxide'
                        E(ii) = self.E_oxide;
                end
            end
            
            A = self.w.*[self.t_bot self.t_mid self.t_top];
            Zm = sum(z.*E.*A)/sum(E.*A);
        end
        
        function Cm = normalized_curvature(self)
            Zm = self.neutral_axis();
            z = [self.t_bot/2 self.t_bot+self.t_mid/2 self.t_bot+self.t_mid+self.t_top/2];
            Z = z - Zm;
            matls = {self.matl_bot self.matl_mid self.matl_top};
            for ii = 1:length(matls)
                switch matls{ii}
                    case 'si'
                        E(ii) = self.E_si;
                    case 'poly'
                        E(ii) = self.E_poly;
                    case 'ti'
                        E(ii) = self.E_ti;
                    case 'al'
                        E(ii) = self.E_al;
                    case 'oxide'
                        E(ii) = self.E_oxide;
                end
            end

            A = self.w.*[self.t_bot self.t_mid self.t_top];
            I = self.w.*[self.t_bot^3/12 self.t_mid^3/12 self.t_top^3/12];
            Cm = 1./sum(E.*(I + A.*Z.^2));
        end
        
        % Spring constant (N/m)
        function k = stiffness(self)
            EIeffective = 1/self.normalized_curvature();
            k = 3*EIeffective/self.l^3;
        end
        
        function freq = resonant_frequency(self)
            matls = {self.matl_bot self.matl_mid self.matl_top};
            for ii = 1:length(matls)
                switch matls{ii}
                    case 'si'
                        density(ii) = self.rho_si;
                    case 'poly'
                        density(ii) = self.rho_poly;
                    case 'ti'
                        density(ii) = self.rho_ti;
                    case 'al'
                        density(ii) = self.rho_al;
                    case 'oxide'
                        density(ii) = self.rho_oxide;
                end
            end
            
            mEff = 0.243*self.l*self.w*sum([self.t_bot self.t_mid self.t_top].*density);
            freq = 1/(2*pi)*sqrt(self.stiffness() / mEff);
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

        % Calculate resistor length, used to calculat resistance
        % and number of carriers
        % Units: m
        function resistor_length = resistor_length(self)
            resistor_length = 2*self.l_pr();
        end

        % Calculate sheet resistance
        % Depends upon abstract method self.doping_profile()
        function Rs = sheet_resistance(self)

            % Assume that the piezoresistor will be symmetric, so just calculate sheet resistance/power for the top
            switch self.matl_top
                case 'poly'
                    mobility = self.mobility(self.dopant_concentration);
                    resistivity = 1/(self.q*mobility*self.dopant_concentration); % ohm-cm
                case 'ti'
                    resistivity = 42e-6; % ohm-cm
                case 'al'
                    resistivity = 2.7e-6; % ohm-cm
            end
            Rs = resistivity/(self.t_top*1e2); % ohms            
        end

        % Data from "Modeling of Carrier Mobility Against Carrier Concentration
        % in Arsenic-,  Phosphorus-, and Boron-Doped  Silicon"
        % Masetti, Serveri and Solmi - IEEE Trans. on Electron Devices, (1983)
        % Units: cm^2/V-sec
        function mobility = mobility(self, dopant_concentration)

            n = dopant_concentration;
            p = dopant_concentration;

            doping_type = 'phosphorus'; % only P for now; could add another parameter later
            
            switch doping_type
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
            
            mobility = mobility*0.5; % assume that the mobility is 50% of the single-crystal value -> totally made up
        end

        % ==================================
        % ======= Calculate noise ==========
        % ==================================

        % The number of current carriers in the piezoresistor
        % Calculated by integrating the number of carriers to the junction
        % depth and multiplaying by the planform dimensions of the piezoresistor
        % Units: unitless
        function number_of_carriers = number_of_carriers(self)
            
            switch self.matl_top
                case 'poly'
                    Nz = self.t_top*self.dopant_concentration*1e6; % units = 1/m^3
                case 'al'
                    Nz = self.t_top*2e21*1e6;
                case 'ti'
                    Nz = self.t_top*2e21*1e6;
            end
            
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
            johnson_PSD = 4*self.k_b*self.T*self.resistance() * ones(1, length(freq));
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
            Q_M = 100; % quality factor
            thermo_PSD = (self.force_sensitivity())^2 * 2*self.stiffness()*self.k_b*self.T/(pi*self.resonant_frequency()*Q_M) * ones(1, length(freq));
        end

        % Integrated Johnson noise
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
            % INA103
            A_VJ = 1.2e-9; % 1.2 nV/rtHz noise floor
            A_IJ = 2e-12; % 2 pA/rtHz noise floor
            A_VF = 6e-9; % 6 nV/rtHz @ 1 Hz
            A_IF = 25e-12; % 25 pA/rtHz @ 1 Hz
            
%             % AD8221
%             A_VJ = 8e-9; % 1.2 nV/rtHz noise floor
%             A_IJ = 40e-15; % 2 pA/rtHz noise floor
%             A_VF = 12e-9; % 6 nV/rtHz @ 1 Hz
%             A_IF = 550e-15; % 25 pA/rtHz @ 1 Hz
            
            
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

        % Ratio of piezoresistor resistance to total resistance (< 1)
        % Assume that metal interconnects are used and that their resistance is about 10% of the total
        function gamma = gamma(self)
            gamma = 0.9;
        end

        function piezoCoefficient = piezoCoefficient(self)
            
            switch self.matl_top
                case 'poly'
                    doping_type = 'phosphorus'; % only phos for now
                    switch doping_type
                        case 'boron'
                            max_factor = 72e-11; % Pi at low concentration in 110 direction
                        case 'phosphorus'
                            max_factor = 103e-11; %Pi at low concentration in 100 direction
                        case 'arsenic'
                            max_factor = 103e-11; %Pi at low concentration in 100 direction
                    end
                    max_factor = 0.6*max_factor; % polycrystalline correction
                    
                    % Richter's model (T=300K), "Piezoresistance in p-type silicon revisited"
                    Nb = 6e19;
                    Nc = 7e20;
                    richter_alpha = 0.43;
                    richter_gamma = 1.6;
                    piezoresistance_factor = (1 + (self.dopant_concentration/Nb).^richter_alpha ...
                        + (self.dopant_concentration/Nc).^richter_gamma).^-1;

                case 'al'
                    max_factor = (1+2*.35)/self.E_al;
                    piezoresistance_factor = 1;
                case 'ti'
                    max_factor = (1+2*.35)/self.E_ti;
                    piezoresistance_factor = 1;
            end
            
            piezoCoefficient = max_factor * piezoresistance_factor;
            
        end
        
        % Units: V/N
        function force_sensitivity = force_sensitivity(self)
            
            v_bias = self.v_bridge/2;
            Zm = self.neutral_axis();
            Cm = self.normalized_curvature();
            Ztop = self.t_bot+self.t_mid+self.t_top/2;

            switch self.matl_top
                case 'poly'
                    E = self.E_poly;
                case 'al'
                    E = self.E_al;
                case 'ti'
                    E = self.E_ti;
            end

            longitudinal_sensitivity = self.number_of_piezoresistors_on_cantilever*self.piezoCoefficient()*E*(Ztop - Zm)*Cm*(self.l - self.l_pr()/2)*v_bias/4*self.gamma();

            force_sensitivity = abs(longitudinal_sensitivity);
        end


        % ==================================
        % == Calculate power dissipation ===
        % ==================================

        % Power dissipation (W) in the cantilever
        function power_dissipation = power_dissipation(self)
            power_dissipation = self.number_of_piezoresistors_on_cantilever * self.v_bridge^2 / (4*self.resistance());
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

        % ==========================================================================
        % =================== Optimization Helper Functions ========================
        % ==========================================================================
        
        % Used by optimization to bring all state varibles to O(1)
        function scaling = optimization_scaling(self)
            l_scale = 1e6;
            w_scale = 1e6;
            t_scale = 1e9;
            l_pr_ratio_scale = 1;
            v_bridge_scale = 1;
            dopant_concentration_scale = 1e-19;
            
            scaling = [l_scale w_scale l_pr_ratio_scale v_bridge_scale ...
                t_scale t_scale t_scale dopant_concentration_scale];
        end
        
        function new_cantilever = cantilever_from_state(self, x0)
            scaling = self.optimization_scaling();
            x0 = x0 ./ scaling;
            
            self.l                    = x0(1);
            self.w                    = x0(2);
            self.l_pr_ratio           = x0(3); 
            self.v_bridge             = x0(4);
            self.t_top                = x0(5);
            self.t_mid                = x0(6);
            self.t_bot                = x0(7);
            self.dopant_concentration = x0(8);
            
            new_cantilever = self;
        end
        
        % Return state vector for the current state
        function x = current_state(self)
            x(1) = self.l;
            x(2) = self.w;
            x(3) = self.l_pr_ratio;
            x(4) = self.v_bridge;
            x(5) = self.t_top;
            x(6) = self.t_mid;
            x(7) = self.t_bot;
            x(8) = self.dopant_concentration;
        end

        function [lb ub] = optimization_bounds(self, parameter_constraints)
            
            min_l = 20e-6;
            max_l = 1e-3;
            
            min_w = 5e-6;
            max_w = 100e-6;
            
            min_l_pr_ratio = 0.01;
            max_l_pr_ratio = 1;
            
            min_v_bridge = 0.1;
            max_v_bridge = 10;

            min_t_top = 50e-9;
            max_t_top = 2e-6;

            min_t_mid = 10e-9;
            max_t_mid = 2e-6;
            
            min_t_bot = 50e-9;
            max_t_bot = 2e-6;
            
            min_dopant_concentration = 1e17;
            max_dopant_concentration = 1e20;

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
            
            lb = [min_l, min_w, min_l_pr_ratio, min_v_bridge, min_t_top, min_t_mid, min_t_bot, min_dopant_concentration];
            ub = [max_l, max_w, max_l_pr_ratio, max_v_bridge, max_t_top, max_t_mid, max_t_bot, max_dopant_concentration];
        end
        
        function x0 = initial_conditions_random(self, parameter_constraints)
            [lb, ub] = self.optimization_bounds(parameter_constraints);
            
            min_l = lb(1);
            max_l = ub(1);

            min_w = lb(2);
            max_w = ub(2);
            
            min_l_pr_ratio = lb(3);
            max_l_pr_ratio = ub(3);
            
            min_v_bridge = lb(4);
            max_v_bridge = ub(4);
            
            min_t_top = lb(5);
            max_t_top = ub(5);
            
            min_t_mid = lb(6);
            max_t_mid = ub(6);
            
            min_t_bot = lb(7);
            max_t_bot = ub(7);
            
            min_dopant_concentration = lb(8);
            max_dopant_concentration = ub(8);
            
            % Generate the random values
            l = min_l + rand*(max_l - min_l);
            w = min_w + rand*(max_w - min_w);
            l_pr_ratio = min_l_pr_ratio + rand*(max_l_pr_ratio - min_l_pr_ratio);
            v_bridge = min_v_bridge + rand*(max_v_bridge - min_v_bridge);
            t_top = min_t_top + rand*(max_t_top - min_t_top);
            t_mid = min_t_mid + rand*(max_t_mid - min_t_mid);
            t_bot = min_t_bot + rand*(max_t_bot - min_t_bot);
            dopant_concentration = 10^(log10(min_dopant_concentration) + ...
                rand*(log10(max_dopant_concentration) - log10(min_dopant_concentration)));

            x0 = [l, w, l_pr_ratio, v_bridge, t_top, t_mid, t_bot, dopant_concentration];
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

            if exist('omega_min_hz', 'var')
                freq_constraint = omega_min_hz - c_new.resonant_frequency();
                C = [C freq_constraint];
            end
            
            if exist('max_power', 'var')
                power_constraint = c_new.power_dissipation() - max_power;
                C = [C power_constraint];
            end

            
            if exist('min_k', 'var')
                min_k_constraint = min_k - c_new.stiffness();
                C = [C min_k_constraint];
            end
            
            if exist('max_k', 'var')
                max_k_constraint = c_new.stiffness() - max_k;
                C = [C max_k_constraint];
            end
            
            % Silicon: L/W > 10, W/T > 3
            length_width_ratio = 5 - c_new.l/c_new.w;
            width_thickness_ratio = 5 - c_new.w/(c_new.t_top + c_new.t_mid + c_new.t_bot);
            
            C = [C length_width_ratio width_thickness_ratio];
            
            % Now for equality based constraints
            Ceq = [];
            
            if exist('fixed_k', 'var')
                fixed_k_constraint = c_new.stiffness() - fixed_k;
                Ceq = [Ceq fixed_k_constraint];
            end
        end

        % This cantilever problem isn't guaranteed to converge, and in
        % practice it fails to converge about 1% of the time for random
        % initial conditions. For this reason, it is best to start from a
        % random initial seed and perform the optimization and checking to
        % make sure that it converges repeatedly.
        function optimized_cantilever = optimize_performance(self, parameter_constraints, nonlinear_constraints, goal)

            percent_match = 0.02; % 2 percent
            randomize_starting_conditions = 1;

            converged = 0;
            ii = 1;
            while ~converged
                % Optimize another cantilever
                [c{ii}, exitflag] = self.optimize_performance_once(parameter_constraints, nonlinear_constraints, goal, randomize_starting_conditions);
                
                % If the optimization terminated abnormally (e.g. constraints not satisfied), skip to the next iteration
                if ~(exitflag == 1 || exitflag == 2)
                    continue
                end
                exitflags(ii) = exitflag;
                
                % Find the resolutions for the cantilevers made so far
                for jj = 1:ii
                    if goal == cantilever.goalForceResolution
                        resolution(jj) = c{jj}.force_resolution();
                    elseif goal == cantilever.goalDisplacementResolution
                        resolution(jj) = c{jj}.displacement_resolution();
                    end
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
                
                % Increment, saving the last calculated value - otherwise we might get stuck
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

            problem.options.TolFun = 1e-6;
            problem.options.TolCon = 1e-6;
            problem.options.TolX = 1e-6;

            problem.options.MaxFunEvals = 2000;
            problem.options.MaxIter = 1000;
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
