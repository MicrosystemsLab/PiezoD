classdef cantilever_diffusion_actuator < cantilever_diffusion
  properties
    
    % Actuator/step dimensions
    l_a;
    w_a;
    t_a;
    actuator_type; % thermal or piezoelectric
    
    % Thermal specific
    v_heater;
    R_heater;
    
    % Piezoelectric specific
    t_electrode;
    v_bias;    
  end
  
  properties (Constant)
    
  end
  
  methods
    function self = cantilever_diffusion_actuator(freq_min, freq_max, ...
        l, w, t, l_pr_ratio, v_bridge, doping_type, diffusion_time, diffusion_temp, ...
        l_a, w_a, t_a, v_heater, R_heater)
      
      self = self@cantilever_diffusion(freq_min, freq_max, l, w, t, l_pr_ratio, v_bridge, doping_type, diffusion_time, diffusion_temp);
      
      self.l_a = l_a;
      self.w_a = w_a;
      self.t_a = t_a;
      self.v_heater = v_heater;
      self.R_heater = R_heater;
    end
    
    function print_performance(self)
    end
    
    function stiffness = stiffness(self)
      
      Cm = 
      k_tip = self.modulus() * self.w * self.t^3 / (4*self.l^3);
      k_base = 3/Cm/self.l_step^3;
      
      % If there is an actuator/reinforcement step at the base, model as two springs in series
      if self.l_step > 0 && self.w_step > 0 && self.t_step > 0

        % Calculate the stiffness of the actuator/step at the base as a bimorph
        z = [self.t/2 self.t+self.t_step/2];
        E = [self.modulus() self.E_step];
        A = [self.w*t_c self.w_step*self.t_step];
        I = 1/12*[self.w*self.t self.w_step*self.t_step].^3;
        z_n = sum(z.*E.*A)/sum(E.*A);
        Cm = 1/sum(E.*(I + A.*(z - z_n).^2));
        
        stiffness = 1/(1/k_base + 1/k_tip);
      % Otherwise, just 
      else
        stiffness = k_tip;
      end
      
      
    end
    
    % ==================================
    % ========= Multilayer mechanics  ==========
    % ==================================
    
    function Zm = neutral_axis(self)
      [z, E, A, I] = self.calculateMechanicsParameters();
      Zm = sum(z.*E.*A)/sum(E.*A);
    end
    
    function Cm = normalized_curvature(self)
      Zm = self.neutral_axis();
      [z, E, A, I] = self.calculateMechanicsParameters();
      Z_offset = z - Zm;
      Cm = 1./sum(E.*(I + A.*Z_offset.^2));
    end
    
    function [z, E, A, I] = calculateMechanicsParameters(self)
      switch self.base_type
        case 'none'
        case 'thermal'
          t = [self.t         self.t_oxide self.t_actuator];
          w = [self.w_base    self.w_c     self.w_actuator];
          E = [self.modulus() self.E_oxide self.E_al];
        case 'piezoelectric'
          
      end
      
      
      for ii = 1:length(t)
        z(ii) = sum(t) - sum(t(ii:end)) + t(ii)/2; % z(1) = t(1)/2, z(2) = t(1) + t(2)/2
      end
      A = w.*t;
      I = (w.*t.^3)/12;
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

      l_a_scale = 1e6;
      w_a_scale = 1e6;
      t_a_scale = 1e9;
      v_heater_scale = 1;
      R_heater_scale = 1e-3;
      
      scaling = [l_scale ...
        w_scale ...
        t_scale ...
        l_pr_ratio_scale ...
        v_bridge_scale ...
        diffusion_time_scale ...
        diffusion_temp_scale ...
        l_a_scale ...
        w_a_scale ...
        t_a_scale ...
        v_heater_scale ...
        R_heater_scale
        ];
    end
    
    function self = cantilever_from_state(self, x0)
      scaling = self.optimization_scaling();
      x0 = x0 ./ scaling;
      
      l = x0(1);
      w = x0(2);
      t = x0(3);
      l_pr_ratio = x0(4);
      v_bridge = x0(5);
      diffusion_time = x0(6);
      diffusion_temp = x0(7);
      
      self.l = l;
      self.w = w;
      self.t = t;
      self.l_pr_ratio = l_pr_ratio;
      self.v_bridge = v_bridge;
      self.diffusion_time = diffusion_time;
      self.diffusion_temp = diffusion_temp;
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
    function [lb ub] = optimization_bounds(self, parameter_constraints)
      min_l = 1e-6;
      max_l = 1e-3;
      
      min_w = 2e-6;
      max_w = 100e-6;
      
      min_t = 1e-6;
      max_t = 50e-6;
      
      min_l_pr_ratio = 0.01;
      max_l_pr_ratio = 1;
      
      min_v_bridge = 0.1;
      max_v_bridge = 10;
      
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
      
      lb = [min_l, min_w, min_t, min_l_pr_ratio, min_v_bridge, ...
        min_diffusion_time, min_diffusion_temp];
      ub = [max_l, max_w, max_t, max_l_pr_ratio, max_v_bridge, ...
        max_diffusion_time, max_diffusion_temp];
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