classdef cantilever_epitaxy < cantilever
    properties
        dopant_concentration
        t_pr_ratio
    end

    methods
        function self = cantilever_epitaxy(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, dopant_concentration, t_pr_ratio)
            self = self@cantilever(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type); % Call superclass constructor
            self.dopant_concentration = dopant_concentration;
            self.t_pr_ratio = t_pr_ratio;
        end
        
        function print_performance(self)
            print_performance@cantilever(self); % print the base stuff
            fprintf('Dopant concentration (per cc): %g \n', self.dopant_concentration);
            fprintf('PR Thickness Ratio: %g \n', self.t_pr_ratio);
            fprintf('PR Thickness: %g \n', self.junction_depth());
            fprintf('\n')
        end
        
        function x_j = junction_depth(self)
            x_j = self.t * self.t_pr_ratio;
        end
        
        
        function print_performance_for_excel(self)
            fprintf('%s\t', self.doping_type);
            variables_to_print = [self.freq_min, self.freq_max*1e-3, ...
                self.l*1e6, self.w*1e6, self.t*1e6, ...
                self.l_pr()*1e6, self.junction_depth()*1e9, ...
                self.l_pr_ratio, self.t_pr_ratio, ...
                self.v_bridge, self.resistance()*1e-3, self.sheet_resistance(), self.power_dissipation()*1e3, self.approxTipDeltaTemp(), ...
                self.number_of_piezoresistors, ...
                self.force_resolution()*1e12, self.displacement_resolution()*1e9, ...
                self.omega_vacuum_hz()*1e-3, self.omega_damped_hz()*1e-3, ...
                self.force_sensitivity(), self.beta(), self.stiffness()*1e3, ...
                self.integrated_noise()*1e6, self.johnson_integrated()*1e6, ...
                self.hooge_integrated()*1e6, self.amplifier_integrated()*1e6, ...
                self.knee_frequency(), self.number_of_carriers()];            
           
            for print_index = 1:length(variables_to_print)
               fprintf('%4g\t', variables_to_print(print_index)); 
            end
            fprintf('\n');
        end
                    
        % ==================================
        % ======== Doping Profile  =========
        % ==================================        
        
        % Return the diffusion profile down to the junction depth
        function [x, doping] = doping_profile(self)
            n_points = 5e2; % # of points of doping profile
            x = linspace(0, self.junction_depth(), n_points);
            doping = ones(1, n_points)*self.dopant_concentration;
        end

        % ==================================
        % ========= Optimization  ==========
        % ==================================        

        % Used by optimization to bring all state varibles to O(1)
        function scaling = optimization_scaling(self)
            l_scale = 1e6;
            w_scale = 1e6;
            t_scale = 1e9;
            l_pr_ratio_scale = 10;
            v_bridge_scale = 1;
            concentration_scale = 1e-18;
            t_pr_ratio_scale = 10;
            scaling = [l_scale, w_scale, t_scale, l_pr_ratio_scale, v_bridge_scale, concentration_scale, t_pr_ratio_scale];
        end
        
        function self = cantilever_from_state(self, x0)
            scaling = self.optimization_scaling();
            x0 = x0 ./ scaling;

            l = x0(1);
            w = x0(2);
            t = x0(3);
            l_pr_ratio = x0(4);
            v_bridge = x0(5);
            dopant_concentration = x0(6);
            t_pr_ratio = x0(7);
            
            self.l = l;
            self.w = w;
            self.t = t;
            self.l_pr_ratio = l_pr_ratio;
            self.v_bridge = v_bridge;
            self.dopant_concentration = dopant_concentration;
            self.t_pr_ratio = t_pr_ratio;
        end

        % Return state vector for the current state
        function x = current_state(self)
            x(1) = self.l;
            x(2) = self.w;
            x(3) = self.t;
            x(4) = self.l_pr_ratio;
            x(5) = self.v_bridge;
            x(6) = self.dopant_concentration;
            x(7) = self.t_pr_ratio;
        end
        
        % Set the minimum and maximum bounds for the cantilever state
        % variables. Bounds are written in terms of the initialization
        % variables. Secondary constraints (e.g. power dissipation,
        % piezoresistor thickness rather than ratio, resonant frequency)
        % are applied in optimization_constraints()
        function [lb ub] = optimization_bounds(self, parameter_constraints)
           
            min_l = 10e-6;
            max_l = 10e-3;
            
            min_w = 1e-6;
            max_w = 500e-6;
            
            min_t = 1e-6;
            max_t = 100e-6;
            
            min_l_pr_ratio = 0.01;
            max_l_pr_ratio = 0.99;
            
            min_v_bridge = 0.1;
            max_v_bridge = 10;
            
            min_dopant_concentration = 1e17;
            
            % Use solid solubility at 800C
            switch self.doping_type
                case 'boron'
                    max_dopant_concentration = 4.4e19;
                case 'phosphorus'
                    max_dopant_concentration = 2.9e20;
                case 'arsenic'
                    max_dopant_concentration = 2.3e19;
            end
            
            min_t_pr_ratio = 0.01;
            max_t_pr_ratio = 0.99;
            
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
            
            lb = [min_l, min_w, min_t, min_l_pr_ratio, min_v_bridge, min_dopant_concentration, min_t_pr_ratio];
            ub = [max_l, max_w, max_t, max_l_pr_ratio, max_v_bridge, max_dopant_concentration, max_t_pr_ratio];
        end
        
        function x0 = initial_conditions_random(self, parameter_constraints)
            [lb, ub] = self.optimization_bounds(parameter_constraints);
            
            % Random generation bounds. We use the conditions from
            % optimization_bounds so that we don't randomly generate
            % something outside of the allowable bounds.
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

            n_min = lb(6);
            n_max = ub(6);

            t_pr_ratio_min = lb(7);
            t_pr_ratio_max = ub(7);
            
            % Generate the random values
            l = l_min + rand*(l_max - l_min);
            w = w_min + rand*(w_max - w_min);
            t = t_min + rand*(t_max - t_min);
            l_pr_ratio = l_pr_ratio_min + rand*(l_pr_ratio_max - l_pr_ratio_min);
            v_bridge = V_b_min + rand*(V_b_max - V_b_min);
            dopant_concentration = 10^(log10(n_min) + rand*(log10(n_max) - log10(n_min))); % logarithmically distributed
            t_pr_ratio = t_pr_ratio_min + rand*(t_pr_ratio_max - t_pr_ratio_min);

            x0 = [l, w, t, l_pr_ratio, v_bridge, dopant_concentration, t_pr_ratio];
        end        
    end
end