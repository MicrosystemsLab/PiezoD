classdef cantilever_diffusion < cantilever
  properties
    diffusion_time
    diffusion_temp
  end
  
  methods
    % Call superclass constructor
    function self = cantilever_diffusion(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, ...
        diffusion_time, diffusion_temp)
      
      self = self@cantilever(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type);
      
      self.diffusion_time = diffusion_time;
      self.diffusion_temp = diffusion_temp;
    end
    
    
    function print_performance(self)
      print_performance@cantilever(self); % print the base stuff
      fprintf('Diffusion time (mins), temp (C): %f %f \n', self.diffusion_time/60, self.diffusion_temp-273);
      fprintf('Junction depth (nm): %f', self.junction_depth()*1e9);
      fprintf('\n')
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
    
    % ==================================
    % ======== Doping Profile  =========
    % ==================================
    
    % Calculate the diffusion profile for a constant surface source
    % diffusion using self.diffusion_time and self.diffusion_temp as
    % as well as self.doping_type
    function [x, doping] = doping_profile(self)
      N_background = 1e14; % N/cm^3
      N_surface = 1e20; % N/cm^3
      n_points = 10e2; % # of points of doping profile
      
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
          
          % A model for phosphorus diffusion which accounts for the
          % kink. Takes the solid solubility to be 2.1e20/cc, which
          % is accurate for ~850C diffusion temp.
        case 'phosphorus'
          k_b_eV = 8.617343e-5;
          x = linspace(0, self.t*1e2, n_points); % m -> cm
          T = self.diffusion_temp;
          t = self.diffusion_time;
          
          % Temperature dependent solid solubility for phosphorus
          % from the Trumbore data included in the TSuprem manual
          SS_temp = [650 700 800 900 1000 1100] + 273; %K
          SS_phos = [1.2e20 1.2e20 2.9e20 6e20 1e21 1.2e21]; % N/cc
          SS_phosfit = polyfit(SS_temp, SS_phos, 2);
          
          % % For debuggin the fit
          % figure
          % hold all
          % plot(SS_temp, SS_phos, 'o');
          % plot(SS_temp, polyval(SS_phosfit, SS_temp));
          % hold off
          % pause
          
          Cs = polyval(SS_phosfit, T)*.5; % factor of 1/2 to fit experimentally observed data
          
          
          %                     alpha = .18*exp(-1.75/k_b_eV/T);
          %                     Da = 200*exp(-3.77/k_b_eV/T);
          %                     Db = 1.2*exp(-3/k_b_eV/T);
          %                     Cb = 3.5*exp(-0.9/k_b_eV/T)*1e23;
          
          alpha = 6e2*exp(-2.45/k_b_eV/T);
          Da = 100*exp(-3.77/k_b_eV/T);
          Db = .8*exp(-3/k_b_eV/T);
          Cb = 1.5*exp(-0.9/k_b_eV/T)*1e23;
          
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
      x_j = x(find(doping == 1e15, 1));
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
      min_diffusion_time = 10*60; % seconds
      max_diffusion_time = 60*60;
      
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
      diffusion_time_random = diffusion_time_min + rand*(diffusion_time_max - diffusion_time_min);
      diffusion_temp_random = diffusion_temp_min + rand*(diffusion_temp_max - diffusion_temp_max);
      
      x0 = [diffusion_time_random, diffusion_temp_random];
    end
  end
end