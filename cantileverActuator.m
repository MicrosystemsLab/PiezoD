classdef cantileverActuator
  properties
    l_c;
    w_c;
    t_c;
    
    l_a;
    w_a;
    t_a;    
  end
  
  properties (Constant)
    numXPoints = 500;
    k_si = 130; % W/m-K
    h_air = 1000; % W/m^2-K
    h_water = 20000; % W/m^2-k
    
    E_si = 130e9; % Assume <100> orientation
    E_oxide = 70e9;
    E_al = 70e9;
    E_aln = 396e9;
    E_ti = 90e9;
    d31_aln = 2e-12;
  end
  
  methods (Abstract)
    calculateActuatorStress(self)
    calculateMechanicsParameters(self)
  end
  
  methods
    function self = cantileverActuator(l_c, w_c, t_c, l_a, w_a, t_a)
      self.l_c = l_c;
      self.w_c = w_c;
      self.t_c = t_c;
      self.l_a = l_a;
      self.w_a = w_a;
      self.t_a = t_a;
    end
    
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
    
    function [x, deflection] = calculateDeflection(self)
      x = linspace(0, self.l_c, self.numXPoints);
      dx = self.l_c/(self.numXPoints - 1);
      
      M = 0; % external moment is zero
      P = 0; % external load is zero
      
      [z, E, A, I] = self.calculateMechanicsParameters();
      stress = self.calculateActuatorStress();

      % Calculate the curvature and neutral axis
      % The curvature may vary as a function of position (especially for thermal actuation), so calculate the
      % deflection by calculating the angle (from the curvature) and numerically integrating.
      C = zeros(length(x), 1);
      Zn = self.t_c/2*ones(length(x), 1); % At centroid by default
      
      % Calculate the curvature, C, and the neutral axis, Zn, along the cantilever length
      for ii = 1:length(x)
        if x(ii) <= self.l_a
          C(ii) = -((M - sum(z.*A.*stress(ii,:)))*sum(E.*A) + (P + sum(A.*stress(ii,:)))*sum(E.*z.*A))/ ...
            (sum(E.*A)*sum(E.*(I+A.*z.^2)) - sum(z.*E.*A)^2);
          Zn(ii) = -((M - sum(z.*A.*stress(ii,:)))*sum(z.*E.*A) + (P + sum(A.*stress(ii,:)))*sum(E.*(I + A.*z.^2)))/ ...
            ((M - sum(z.*A.*stress(ii,:)))*sum(E.*A) + (P + sum(A.*stress(ii,:)))*sum(E.*z.*A));
        end
      end
      
      theta = cumsum(C.*dx);
      deflection = cumsum(theta.*dx);
    end
    
    function z_tip = tipDeflection(self)
      [x, z] = self.calculateDeflection();
      z_tip = max(z);
    end
  end
end