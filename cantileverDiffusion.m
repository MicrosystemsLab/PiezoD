% Model a diffused cantilever (particularly phosphorus/POCl3 diffusion)
classdef cantileverDiffusion < cantilever
  properties
    diffusion_time
    diffusion_temp
  end
  
  methods
    % Call superclass constructor
    function self = cantileverDiffusion(freq_min, freq_max, l, w, t, ...
			l_pr_ratio, v_bridge, doping_type, diffusion_time, diffusion_temp)
      
      self = self@cantilever(freq_min, freq_max, l, w, t, l_pr_ratio, ...
				v_bridge, doping_type);
      
      self.diffusion_time = diffusion_time;
      self.diffusion_temp = diffusion_temp;
    end
    
    function print_performance(self)
      print_performance@cantilever(self); % print the base class stuff
      fprintf('Diffusion time (mins), temp (C): %f %f \n', ...
				self.diffusion_time/60, self.diffusion_temp-273);
      fprintf('Junction depth (nm): %f \n', self.junction_depth()*1e9);
      fprintf('=======================\n\n')
    end
    
    function print_performance_for_excel(self, varargin)
      
      % Call the superclass method first
      print_performance_for_excel@cantilever(self, varargin);
      
      % varargin gets another set of {} from the cantilever subclasses
      optargin = size(varargin, 2);
      if optargin == 1
        fid = varargin{1};
      elseif optargin == 0
        fid = 1; % Print to the stdout
      else
        fprintf('ERROR: Extra optional arguments')
      end

      % Then print out our additional variables
      variables_to_print = [self.diffusion_time/60, self.diffusion_temp-273];
      for print_index = 1:length(variables_to_print)
        fprintf(fid, '%4g\t', variables_to_print(print_index));
      end
      fprintf(fid, '\n');
    end
    
    % Calculate the diffusion profile for a constant surface source
    % diffusion using self.diffusion_time and self.diffusion_temp as
    % as well as self.doping_type
    function [x, active_doping, total_doping] = doping_profile(self)
      N_background = 1e15; % N/cm^3
      N_surface = 1e20; % N/cm^3
      n_points = self.numZPoints; % # of points of doping profile
      
      switch self.doping_type
        case 'arsenic'
          D_0 = 9.17; % cm^2/sec
          E_a = 3.99; % eV
          
          % Simple diffusion model
          D = D_0*exp(-E_a/(self.k_b_eV*self.diffusion_temp));
          diffusion_length = sqrt(D*self.diffusion_time)*1e-2; % cm -> m
          
          junction_depth = 2*diffusion_length*erfcinv(N_background/N_surface);
          x = linspace(0, self.t, n_points);
          
          active_doping = N_surface*erfc(x/(2*diffusion_length));
          total_doping = active_doping;
          
        case 'boron'
          D_0 = 1.0;
          E_a = 3.5;
          
          % Simple diffusion model
          D = D_0*exp(-E_a/(self.k_b_eV*self.diffusion_temp));
          diffusion_length = sqrt(D*self.diffusion_time)*1e-2; % cm -> m
          
          junction_depth = 2*diffusion_length*erfcinv(N_background/N_surface);
          x = linspace(0, self.t, n_points);
          
          active_doping = N_surface*erfc(x/(2*diffusion_length));
          total_doping = active_doping;
          
        % A model for phosphorus diffusion which accounts for the
        % kink. Takes the solid solubility to be 2.1e20/cc, which
        % is accurate for ~850C diffusion temp.
        case 'phosphorus'
          k_b_eV = 8.617343e-5;
          x = linspace(0, self.t*1e2, n_points); % m -> cm
          
          % Current (7/18/10)
          if self.diffusion_temp < 780+273
            time_offset = 0;
          else
            time_offset = 2;
          end
          
          % Alpha activation energy varies quadratically with temp
          T_values = 273+[775 800 850 950];
          Ea_values = [1.755 1.73 1.71 1.66];
          p = polyfit(T_values, Ea_values, 2);
          Ea_alpha = polyval(p, self.diffusion_temp);
          alpha = .18*exp(-Ea_alpha./(k_b_eV*self.diffusion_temp));
          
          Da = 100*exp(-3.77/k_b_eV/self.diffusion_temp);
          Db = 2.3*exp(-1.95./(k_b_eV*self.diffusion_temp))*1e-5;
          Cb = 3*exp(-0.88./(k_b_eV*self.diffusion_temp))*1e23;
          
					% Add some offset time to account for the oxidation/purge steps
          t = self.diffusion_time + time_offset*60;
          
          % Temperature dependent electrically active concentration
					% Source: Solmi and Nobili papers (various)
          surface_concentration_total = 2.5e23* ...
						exp(-0.62./(k_b_eV*self.diffusion_temp));
          surface_concentration_active = 1.1e22* ...
						exp(-0.37./(k_b_eV*self.diffusion_temp));
          
          % Calculate the profile
          x0 = alpha*t;
          kappa = Cb/surface_concentration_total;
          
          F1 = erfc((x+alpha*t)/(2*sqrt(Da*t))) + ...
						erfc((x-3*alpha*t)/(2*sqrt(Da*t)));
          F2 = erfc((x+alpha*t)/(2*sqrt(Db*t))) + ...
						erfc((x-3*alpha*t)/(2*sqrt(Db*t)));
          
          Ca = (1-kappa)/2*surface_concentration_total* ...
						exp(-alpha/(2*Da)*(x-alpha*t)).*F1;
          Cb = kappa/2*surface_concentration_total* ...
						exp(-alpha/(2*Db)*(x-alpha*t)).*F2;
          
          C_total = Ca + Cb;
          C_total(find(C_total < 1e15)) = 1e15;
          C_total(x <= x0) = surface_concentration_total;
          total_doping = C_total;
          
          C_active = C_total;
          C_active(C_active >= surface_concentration_active) = ...
						surface_concentration_active;
          active_doping = C_active;    
          
          x = x*1e-2; % cm -> m
      end
    end
    
    % Effective carrier density, N/m^2
    function Nz = Nz(self)
      [z, N_active, N_total] = self.doping_profile();
      [mu, sigma] = self.mobility(N_active, self.T);
%       Nz_total = trapz(z, N_active*1e6) % doping: N/cm^3 -> N/m^3
      Nz = 1e4*trapz(z*1e2, N_active.*mu)^2/trapz(z*1e2, N_active.*mu.^2);    
    end

    function plot_doping_profile(self)
      [x, active_doping, total_doping] = self.doping_profile();
      figure
      hold all
      plot(x, total_doping)
      plot(x, active_doping)
      hold off
      set(gca, 'yscale', 'log');
      xlim([0 self.t])
    end
    
    function x_j = junction_depth(self)
      [x, active_doping, total_doping] = self.doping_profile();
      x_j = x(find(active_doping == 1e15, 1));
    end
    
    function alpha = alpha(self)
      alpha = self.default_alpha; % use the alpha from the superclass
    end
    
    % ==================================
    % ========= Optimization  ==========
    % ==================================
    
    function scaling = doping_optimization_scaling(self)
      diffusion_time_scale = 1e-3;
      diffusion_temp_scale = 1e-3;
      scaling = [diffusion_time_scale diffusion_temp_scale];
    end
    
    function self = doping_cantilever_from_state(self, x0)
      self.diffusion_time = x0(6);
      self.diffusion_temp = x0(7);
    end
    
    function x = doping_current_state(self)
      x(1) = self.diffusion_time;
      x(2) = self.diffusion_temp;
    end
    
    % Set the minimum and maximum bounds for the cantilever state
    % variables. Bounds are written in terms of the initialization
    % variables. Secondary constraints (e.g. power dissipation,
    % piezoresistor thickness rather than ratio, resonant frequency)
    % are applied in optimization_constraints()
    function [lb ub] = doping_optimization_bounds(self, parameter_constraints)
      min_diffusion_time = 5*60; % seconds
      max_diffusion_time = 90*60;
      
      min_diffusion_temp = 273+800; % K
      max_diffusion_temp = 273+1000;
      
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
      
      lb = [min_diffusion_time, min_diffusion_temp];
      ub = [max_diffusion_time, max_diffusion_temp];
    end
    
    function x0 = doping_initial_conditions_random(self, parameter_constraints)
      [lb, ub] = self.doping_optimization_bounds(parameter_constraints);
      
      diffusion_time_min = lb(1);
      diffusion_time_max = ub(1);
      
      diffusion_temp_min = lb(2);
      diffusion_temp_max = ub(2);
      
      
      % Generate the random values
      diffusion_time_random = diffusion_time_min + ...
				rand*(diffusion_time_max - diffusion_time_min);
      diffusion_temp_random = diffusion_temp_min + ...
				rand*(diffusion_temp_max - diffusion_temp_min);
      
      x0 = [diffusion_time_random, diffusion_temp_random];
    end
  end
end