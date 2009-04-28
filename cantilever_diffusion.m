classdef cantilever_diffusion < cantilever
    properties
        diffusion_time
        diffusion_temp
    end
    
    methods
        % Call superclass constructor
        function self = cantilever_diffusion(freq_min, freq_max, ...
                l, w, t, l_pr_ratio, v_bridge, doping_type, ...
                diffusion_time, diffusion_temp)

            self = self@cantilever(freq_min, freq_max, l, w, t, l_pr_ratio, ...
                v_bridge, doping_type);
                
            self.diffusion_time = diffusion_time;
            self.diffusion_temp = diffusion_temp;
        end
        
        
        function print_performance(self)
            print_performance@cantilever(self); % print the base stuff
            fprintf('Diffusion time (mins), temp (C): %f %f \n', self.diffusion_time/60, self.diffusion_temp-273);
            fprintf('Junction depth (nm): %f', self.junction_depth()*1e9);
            fprintf('\n')
        end

        function print_performance_for_excel(self)
            fprintf('%s \t', self.doping_type);
            variables_to_print = [self.l*1e6, self.w*1e6, self.t*1e6, ...
                self.l_pr_ratio, self.v_bridge ...
                self.diffusion_time, self.diffusion_temp, ...
                self.freq_min, self.freq_max, ...
                self.force_resolution(), self.displacement_resolution(), ...
                self.omega_vacuum_hz(), self.omega_damped_hz(), ...
                self.force_sensitivity(), self.beta(), self.junction_depth()*1e9,  ...
                self.stiffness(), self.quality_factor(), ...
                self.resistance(), self.power_dissipation()*1e3, ...
                self.integrated_noise(), self.integrated_johnson_noise(), ...
                self.integrated_hooge_noise(), self.knee_frequency()];            
           
            for print_index = 1:length(variables_to_print)
               fprintf('%g \t', variables_to_print(print_index)); 
            end
            fprintf('\n');
        end

        % ==================================
        % ======== Doping Profile  =========
        % ==================================        
          
        % Calculate the diffusion profile for a constant surface source
        % diffusion using self.diffusion_time and self.diffusion_temp as
        % as well as self.doping_type
        function [x, doping] = doping_profile(self)
            N_background = 1e14; % N/cm^3
            N_surface = 1e20; % N/cm^3
            n_points = 1000; % # of points of doping profile
            
            switch self.doping_type
                case 'arsenic'
                    D_0 = 9.17; % cm^2/sec
                    E_a = 3.99; % eV

                    % Simple diffusion model
                    D = D_0*exp(-E_a/(self.k_b_eV*self.diffusion_temp));
                    diffusion_length = sqrt(D*self.diffusion_time)*1e-2; % cm -> m

                    junction_depth = 2*diffusion_length*erfcinv(N_background/N_surface);
                    x = linspace(0, junction_depth, n_points);

                    doping = N_surface*erfc(x/(2*diffusion_length));
                case 'boron'
                    D_0 = 1.0;
                    E_a = 3.5;

                    % Simple diffusion model
                    D = D_0*exp(-E_a/(self.k_b_eV*self.diffusion_temp));
                    diffusion_length = sqrt(D*self.diffusion_time)*1e-2; % cm -> m

                    junction_depth = 2*diffusion_length*erfcinv(N_background/N_surface);
                    x = linspace(0, junction_depth, n_points);

                    doping = N_surface*erfc(x/(2*diffusion_length));
                case 'phosphorus'
                    k_b_eV = 8.617343e-5;
                    x = linspace(0, self.t*1e2, n_points); % m -> cm
                    Cs = 2.1e20;

                    T = self.diffusion_temp;
                    t = self.diffusion_time;

                    alpha = .18*exp(-1.75/k_b_eV/T);
                    Da = 200*exp(-3.77/k_b_eV/T);
                    Db = 1.2*exp(-3/k_b_eV/T);
                    Cb = 3.5*exp(-0.9/k_b_eV/T)*1e23;

                    x0 = alpha*t;
                    kappa = Cb/Cs;

                    F1 = erfc((x+alpha*t)/(2*sqrt(Da*t))) + erfc((x-3*alpha*t)/(2*sqrt(Da*t)));
                    F2 = erfc((x+alpha*t)/(2*sqrt(Db*t))) + erfc((x-3*alpha*t)/(2*sqrt(Db*t)));

                    Ca = (1-kappa)/2*Cs*exp(-alpha/(2*Da)*(x-alpha*t)).*F1;
                    Cb = kappa/2*Cs*exp(-alpha/(2*Db)*(x-alpha*t)).*F2;

                    C = Ca + Cb;
                    C(find(x <= x0)) = Cs;
                    C(find(C < 1e15)) = 1e15;
                    doping = C;
                    x = x*1e-2; % cm -> m
            end

        end         
        
        function x_j = junction_depth(self)
            [x, doping] = self.doping_profile();
            x_j = max(x);
        end   
        
        % ==================================
        % ========= Optimization  ==========
        % ==================================        
        
        % Used by optimization to bring all state varibles to O(1)
        function scaling = optimization_scaling(self)
            l_scale = 1e6;
            w_scale = 1e6;
            t_scale = 1e9;
            l_pr_ratio_scale = 1;
            v_bridge_scale = 1;
            diffusion_time_scale = 1e-3;
            diffusion_temp_scale = 1e-3;
            
            scaling = [l_scale ...
                       w_scale ...
                       t_scale ...
                       l_pr_ratio_scale ...
                       v_bridge_scale ...
                       diffusion_time_scale ...
                       diffusion_temp_scale];
        end
        
        function new_cantilever = cantilever_from_state(self, x0)
            scaling = self.optimization_scaling();
            x0 = x0 ./ scaling;

            l = x0(1);
            w = x0(2);
            t = x0(3);

            l_pr_ratio = x0(4);
            
            v_bridge = x0(5);
            
            diffusion_time = x0(6);
            diffusion_temp = x0(7);
            
            new_cantilever = cantilever_diffusion(self.freq_min, self.freq_max, ...
                l, w, t, l_pr_ratio, v_bridge, self.doping_type, ...
                diffusion_time, diffusion_temp);
        end
        
        % Return state vector for the current state
        function x = current_state(self)
            x(1) = self.l;
            x(2) = self.w;
            x(3) = self.t;
            x(4) = self.l_pr_ratio;
            x(5) = self.v_bridge;
            x(6) = self.diffusion_time;
            x(7) = self.diffusion_temp;
        end
        
        % Set the minimum and maximum bounds for the cantilever state
        % variables. Bounds are written in terms of the initialization
        % variables. Secondary constraints (e.g. power dissipation,
        % piezoresistor thickness rather than ratio, resonant frequency)
        % are applied in optimization_constraints()
        function [lb ub] = optimization_bounds(self, constraints)
            min_l = 1e-6;
            max_l = 10e-3;
            
            min_w = 2e-6;
            max_w = 10e-3;
            
            min_t = 1e-6;
            max_t = 10e-3;
            
            min_l_pr_ratio = 0.01;
            max_l_pr_ratio = 1;
            
            min_v_bridge = 0.1;
            max_v_bridge = 5;
            
            min_diffusion_time = 10*60; % seconds
            max_diffusion_time = 15*60;
            
            min_diffusion_temp = 273+800; % K
            max_diffusion_temp = 273+900;
            
            % Override the default values if any were provided
            % constraints is a set of key value pairs, e.g.
            % constraints = {{'min_l', 'max_v_bridge'}, {5e-6, 10}}
            if length(constraints) > 0
                keys = constraints{1};
                values = constraints{2};
                for ii = 1:length(keys)
                    eval([keys{ii} '=' num2str(values{ii})]);
                end
            end
            
            lb = [min_l, min_w, min_t, min_l_pr_ratio, min_v_bridge, ...
                min_diffusion_time, min_diffusion_temp];
            ub = [max_l, max_w, max_t, max_l_pr_ratio, max_v_bridge, ...
                max_diffusion_time, max_diffusion_temp];
        end
        
        function x0 = initial_conditions_random(self, constraints)
            [lb, ub] = self.optimization_bounds(constraints);
            
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

            v_bridge_min = lb(5);
            v_bridge_max = ub(5);

            diffusion_time_min = lb(6);
            diffusion_time_max = ub(6);

            diffusion_temp_min = lb(7);
            diffusion_temp_max = ub(7);

            % Generate the values
            l = l_min + rand*(l_max - l_min);
            w = w_min + rand*(w_max - w_min);
            t = t_min + rand*(t_max - t_min);

            l_pr_ratio = l_pr_ratio_min + rand*(l_pr_ratio_max - l_pr_ratio_min);

            v_bridge = v_bridge_min + rand*(v_bridge_max - v_bridge_min);

            diffusion_time = diffusion_time_min + rand*(diffusion_time_max - diffusion_time_min);
            diffusion_temp = diffusion_temp_min + rand*(diffusion_temp_max - diffusion_temp_max);

            x0 = [l, w, t, l_pr_ratio, v_bridge, diffusion_time, diffusion_temp];
        end     
    end
end