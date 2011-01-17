% Notes

% M-Lint Settings
%#ok<*PROP>


classdef cantilever_piezoelectric
  
  properties
    % Constants
    n_plot_points = 1e3; % # of points to use for calculating and plotting resolutions
    
    freq_min;
    freq_max;
    l_si;
    w_si;
    t_si;
    l_pe;
    w_pe;
    t_pe;
    material;
    r_shunt;
    
    resistivity;
    permittivity;
    d31;
    fluid;

    
    % Store optional flags in a key->value map
    % amplifierNoise -> if true, include the effect of amplifier noise
    % thermomechanicalNoise -> if true, include thermomechanical noise
    keyAmplifierNoise = 'amplifierNoise';
    keyThermomechanicalNoise = 'thermomechanicalNoise';
    flags = containers.Map();
    
    bondpad_size = 200e-6;
    C_amp = 0.2e-12; % input capacitance of amplifier
    
    t_bottom_electrode = 50e-9;
    t_top_electrode = 50e-9;
    rho_electrode = 4506;
    E_electrode = 90e9;
    
    rhoPE;
    EPE;
    
    rhoSi = 2330; % kg/m^3
    ESi = 169e9; % Pa
  end
  
  properties (Constant)
    epsilon_0 = 8.854e-12; % F/m
    k_b = 1.38e-23; % joules/kelvin
    T = 300; % kelvin
    q = 1.6e-19; % C/electron
    
    rho_water = 1e3; % kg/m^3
    eta_water = 0.9e-3; % Pa-sec
    rho_air = 1.2; % kg/m^3
    eta_air = 17e-6; % Pa-sec
    maxQ = 1000;
    minQ = 1e-6;
    
    goalVoltageForce = 0;
    goalVoltageDisplacement = 1;
    goalChargeForce = 2;
    goalChargeDisplacement = 3;
    
    fluidVacuum = 0;
    fluidWater = 1;
    
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
  
  
  methods
    function self = PECantilever(freq_min, freq_max, l_si, w_si, t_si, t_pe, material, r_shunt)
      self.freq_min = freq_min;
      self.freq_max = freq_max;
      self.l_si = l_si;
      self.w_si = w_si;
      self.t_si = t_si;
      
      self.l_pe = l_si;
      self.w_pe = w_si;
      self.t_pe = t_pe;
      
      self = self.updateMaterial(material);
      self.r_shunt = r_shunt;
      
      self.fluid = 'air'; % Default value
      
      self.defaultFlags();
    end
    
    function changeFlag(self, key, value)
      self.flags(key) = value;
    end
    
    % Set the defalt flags
    function defaultFlags(self)
      self.flags(self.keyAmplifierNoise) = 1;
      self.flags(self.keyThermomechanicalNoise) = 1;
    end
    
    function self = updateMaterial(self, material)
      % Based upon a slide from Justin Black (in my Quals folder)
      self.material = material;
      switch material
        case 'AlN'
          self.permittivity = 10*self.epsilon_0; % F/m
          self.resistivity = 1e12; % ohm-cm
          self.d31 = 2e-12; % C/N
          self.EPE = 396e9; % Pa
          self.rhoPE = 3260; % kg/m^3
        case 'PZT'
          self.permittivity = 900*self.epsilon_0; % F/m
          self.resistivity = 1e8; % ohm-cm
          self.d31 = 70e-12; % C/N
          self.EPE = 55e9; % Pa
          self.rhoPE = 7550; % kg/m^3
      end
    end
    
    
    % ==========================================================================
    % ============================ Misc Functions ==============================
    % ==========================================================================
    
    function print_performance(self)
      [omega_damped_hz, Q] = self.omega_damped_hz_and_Q();
      fprintf('Cantilever L/W/T: %.2f %.2f %.2f \n', self.l_si*1e6, self.w_si*1e6, self.t_si*1e6)
      fprintf('PE L/W/T: %.2f %.2f %.2f \n', self.l_pe*1e6, self.w_pe*1e6, self.t_pe*1e6)
      fprintf('PE Thickness Ratio: %.2f \n', self.t_pe/self.t_si)
      fprintf('Shunt Resistance: %.2g \n', self.r_shunt)
      fprintf('\n')
      fprintf('f0       (kHz): %.3g \n', self.resonant_frequency()*1e-3);
      fprintf('f_damped (kHz): %.3g \n', omega_damped_hz*1e-3);
      fprintf('Q: %.3g \n', Q);
      fprintf('k (N/m): %.3g \n', self.stiffness());
      fprintf('\n')
      fprintf('Resistance PE/Half circuit: %.3g %.3g \n', self.Rpe(), self.R_half())
      fprintf('Capacitance PE/Half circuit: %.3g %.3g \n', self.Cpe(), self.C_half())
      fprintf('\n')
      fprintf('Voltage Mode:\n')
      fprintf('Force resolution (N): %.3g \n', self.Fminv())
      fprintf('Force sensitivity (V/N): %.3g \n', self.v_f_sensitivity(1e3)); % calculated @ 1kHz
      fprintf('Displacement resolution (m): %.3g \n', self.Xminv())
      fprintf('Integrated noise (V): %.3g \n', self.Vn_total())
      fprintf('\n')
      fprintf('Charge Mode:\n')
      fprintf('Force resolution (N): %.3g \n', self.Fminq())
      fprintf('Force sensitivity (C/N): %.3g \n', self.q_f_sensitivity(1e3)); % calculated @ 1kHz
      fprintf('Displacement resolution (m): %.3g \n', self.Xminq())
      fprintf('Integrated noise (C): %.3g \n', self.Qn_total())
      fprintf('\n')
      fprintf('Thermomechanical force noise limit: %g \n', self.thermomechanicalForceLimit());
      fprintf('\n\n\n')
    end
    
    % ==========================================================================
    % ============================ Beam Mechanics ==============================
    % ==========================================================================
    
    function [z, E, A, I] = calculateMechanicsParameters(self)
      t = [self.t_si self.t_bottom_electrode self.t_pe self.t_top_electrode];
      
      for ii = 1:length(t)
        z(ii) = sum(t) - sum(t(ii:end)) + t(ii)/2; % z(1) = t(1)/2, z(2) = t(1) + t(2)/2
      end
      E = [self.ESi self.E_electrode self.EPE self.E_electrode];
      A = self.w_si*t;
      I = (self.w_si*t.^3)/12;
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
    
    % Spring constant (N/m)
    function k = stiffness(self)
      EIeffective = 1/self.normalized_curvature();
      k = 3*EIeffective/self.l_pe^3;
    end
    
    function mEff = effectiveMass(self)
      mEff = 0.243*self.w_si*self.l_si*(self.t_si*self.rhoSi + self.t_bottom_electrode*self.rho_electrode + self.t_pe*self.rhoPE + self.t_top_electrode*self.rho_electrode);
    end
    
    function freq = resonant_frequency(self)
      freq = 1/(2*pi)*sqrt(self.stiffness() / self.effectiveMass());
    end
    
    % The average stress normalized per unit force
    function sigmaN = normalized_stress(self)
      sigmaN = self.EPE*self.normalized_curvature()*self.l_pe/2*(self.neutral_axis() - (self.t_si + self.t_pe/2));
    end
    
    % Strain profile per unit force applied to the cantilever tip
    function [z, strain] = normalized_strain_profile(self)
      z = linspace(0, self.t_si + self.t_pe, self.n_plot_points);
      Cm = self.normalized_curvature();
      Zm = self.neutral_axis();
      strain = Cm*self.l_pe*(Zm - z);
    end
    
    
    % ==========================================================================
    % ============================= Sensitivity ================================
    % ==========================================================================
    
    function Qf = q_f_sensitivity_intrinsic(self)
      Zm = self.neutral_axis();
      [z, E, A, I] = self.calculateMechanicsParameters();
      zPE = z(3);
      Cm = self.normalized_curvature();
      Qf = abs(.5*self.d31*self.EPE*(Zm - zPE)*Cm*self.w_pe*self.l_pe^2);
    end
    
    % Note: this works out to be about the same as the intrinsic sensitivity and actually should depend on the force
    % magnitude
    function Qf = q_f_sensitivity(self, freq)
      Qf = ones(1,length(freq))*self.q_f_sensitivity_intrinsic();
    end
    
    function Qx = q_x_sensitivity(self, freq)
      Qx = self.q_f_sensitivity(freq)*self.stiffness();
    end
    
    function Vf = v_f_sensitivity(self, freq)
      Vf = 2*(2*pi*freq*self.q_f_sensitivity_intrinsic())*self.R_half()./(1+2*pi*freq*self.R_half()*self.C_half());
    end
    
    function Vx = v_x_sensitivity(self, freq)
      Vx = self.v_f_sensitivity(freq)*self.stiffness();
    end
    
    % Current sensitivity in units of amps/newton (If = dQ/dt = Q*2*pi*f)
    function If = current_sensitivity(self, freq)
      If = 2*pi*self.q_f_sensitivity(freq).*freq;
    end
    
    % ==========================================================================
    % ======================== R, C calculations ===============================
    % ==========================================================================
    
    function Rpe = Rpe(self)
      area = self.l_pe*self.w_pe;
      Rpe = self.resistivity*1e-2*self.t_pe/area;
    end
    
    function Cpe = Cpe(self)
      area = self.l_pe*self.w_pe;
      Cpe = self.permittivity*area/self.t_pe;
    end
    
    % Rpe/2 || Rshunt (half circuit)
    function R_half = R_half(self)
      R_half = self.Rpe()*self.r_shunt/(self.Rpe() + 2*self.r_shunt);
    end
    
    function C_half = C_half(self)
      C_amp = 0.2e-12;
      C_half = 2*(self.Cpe() + C_amp);
    end
    
    % Rpe || Cpe
    function Z = Z(self, freq)
      Z = 1./(1/self.Rpe() + 2*pi*freq*1i*self.Cpe());
    end
    
    % Half circuit Z, Rpe/2 || 2Cpe
    function Z_half = Z_half(self, freq)
      Z_half = 1./(1/self.R_half() + 2*pi*freq*1i*self.C_half());
    end
    
    
    % ==========================================================================
    % =============================== Noise ====================================
    % ==========================================================================
    
    % Thermomechanical voltage noise
    function Vth = Vth(self, freq)
      [omega_damped, Q_M] = self.omega_damped_and_Q();
      Fth = sqrt(2*self.k_b*self.stiffness()*self.T/(Q_M*pi*self.resonant_frequency()));
      Vth = Fth*self.v_f_sensitivity(freq);
    end
    
    function thermoLimit = thermomechanicalForceLimit(self)
      [omega_damped, Q_M] = self.omega_damped_and_Q();
      thermoLimit = sqrt(self.freq_max - self.freq_min)*sqrt(2*self.k_b*self.stiffness()*self.T/(Q_M*pi*self.resonant_frequency()));
    end
    
    % Amplifier voltage noise. Assumes 'voltage-mode' sensing
    % Considers both the current and voltage noise of the amplifier
    function Vamp = Vamp(self, freq)
      
      % INA116 Amplifier Noise Parameters
      Vamp_J = 28e-9;
      Vamp_H = 300e-9;
      Iamp_J = .1e-15;
      Iamp_H = 1e-15;
      
      Vamp_VJ = Vamp_J;
      Vamp_IJ = Iamp_J*real(self.Z_half(freq));
      Vamp_VH = Vamp_H./sqrt(freq);
      Vamp_IH = Iamp_H.*real(self.Z_half(freq))./sqrt(freq);
      
      % Amplifier current noise power is multiplied by two because it's for the half circuit impedance
      Vamp = sqrt(Vamp_VJ.^2 + Vamp_IJ.^2 + Vamp_VH.^2 + Vamp_IH.^2);
    end
    
    % Includes Johnson, thermomechanical and amplifier
    % Units: V/rtHz
    function Vn = Vn(self, freq)
      Vj = sqrt(4*self.k_b*self.T*self.R_half())./(1 + 2*pi*freq*self.R_half()*self.C_half());
      
      if self.flags(self.keyThermomechanicalNoise)
        Vth = self.Vth(freq);
      else
        Vth = 0;
      end
      
      if self.flags(self.keyAmplifierNoise)
        Vamp = self.Vamp(freq);
      else
        Vamp = 0;
      end
      
      Vn = sqrt(Vj.^2 + Vth.^2 + Vamp.^2);
    end
    
    function Qamp = Qamp(self, freq)
      
      % INA116 Noise Parameters
      Vamp_J = 28e-9; % 28 nV/rtHz
      Vamp_H = 300e-9; % 300 nV/rtHz @ 1Hz
      Iamp_J = .1e-15; % .1 fA/rtHz
      Iamp_H = 1e-15; % 1 fA/rtHz (estimated from ratio of Vh to Vj
      
      % Charge noise induced by the amplifier voltage noise - independent of frequency
      Qamp_VJ = Vamp_J*self.Cpe();
      Qamp_IJ = Iamp_J./(2*pi*freq); % integrated amplifier current noise
      Qamp_VH = Vamp_H*self.Cpe()./sqrt(freq);
      Qamp_IH = Iamp_H./sqrt(freq)./(2*pi*freq);
      
      % Integrated amplifier current noise
      Qamp = sqrt(Qamp_VJ.^2 + Qamp_IJ.^2 + Qamp_VH.^2 + Qamp_IH.^2);
      
      
      
      %             % Debugging plots
      %
      %             % Plot the calculated amplifier noise
      %             figure
      %             plot(freq, sqrt(Vamp_J.^2 + (Vamp_H./sqrt(freq)).^2));
      %             set(gca,'xscale','log','yscale','log');
      %             pause
      %
      %             % Plot plot the calculated charge noise
      %             figure
      %             hold all
      %             plot(freq, Qamp_VJ.*ones(1, length(freq)))
      %             plot(freq, Qamp_VH)
      %             plot(freq, Qamp_IJ)
      %             plot(freq, Qamp_IH)
      %             hold off
      %             set(gca,'xscale','log','yscale','log');
      %             legend('VJ','VH','IJ','IH');
      %             pause
      
      
    end
    
    
    function Qn = Qn(self, freq)
      Qj = self.Cpe()*sqrt(4*self.k_b*self.T*self.Rpe())./(1 + 2*pi*freq*self.Rpe()*self.Cpe());
      Qth = self.Vth(freq)*self.Cpe();
      Qamp = self.Qamp(freq);
      Qn = sqrt(Qj.^2 + Qth.^2 + Qamp.^2);
      
      
      %             % Debugging plot
      %             figure
      %             hold all
      %             plot(freq, Qj);
      %             plot(freq, Qth);
      %             plot(freq, Qamp);
      %             set(gca, 'yscale','log','xscale','log');
      %             legend('Johnson','Thermal','Amp');
      %             pause
    end
    
    function Vn_total = Vn_total(self)
      freq = self.freq_for_plot();
      Vn = self.Vn(freq);
      Vn_total = sqrt(trapz(freq, Vn.^2));
    end
    
    function Qn_total = Qn_total(self)
      freq = self.freq_for_plot();
      Qn = self.Qn(freq);
      Qn_total = sqrt(trapz(freq, Qn.^2));
    end
    
    
    function Vn_cumulative = Vn_cumulative(self)
      freq = self.freq_for_plot();
      Vn = self.Vn(freq);
      Vn_cumulative = sqrt(cumtrapz(freq, Vn.^2));
    end
    
    function Qn_cumulative = Qn_cumulative(self)
      freq = self.freq_for_plot();
      Qn = self.Qn(freq);
      Qn_cumulative = sqrt(cumtrapz(freq, Qn.^2));
    end
    
    
    % ==========================================================================
    % ============================= Resolution =================================
    % ==========================================================================
    
    
    function Sfminv = Sfminv(self, freq)
      Sfminv = self.Vn(freq)./self.v_f_sensitivity(freq);
    end
    
    function Sxminv = Sxminv(self, freq)
      Sxminv = self.Vn(freq)./self.v_x_sensitivity(freq);
    end
    
    function Sfminq = Sfminq(self, freq)
      Sfminq = self.Qn(freq)./self.q_f_sensitivity(freq);
    end
    
    function Sxminq = Sxminq(self, freq)
      Sxminq = self.Qn(freq)./self.q_x_sensitivity(freq);
    end
    
    
    function Fminv = Fminv(self)
      freq = self.freq_for_plot();
      Sfminv = self.Sfminv(freq);
      Fminv = sqrt(trapz(freq, Sfminv.^2));
    end
    
    function Xminv = Xminv(self)
      freq = self.freq_for_plot();
      Sxminv = self.Sxminv(freq);
      Xminv = sqrt(trapz(freq, Sxminv.^2));
    end
    
    function Fminq = Fminq(self)
      freq = self.freq_for_plot();
      Sfminq = self.Sfminq(freq);
      Fminq = sqrt(trapz(freq, Sfminq.^2));
    end
    
    function Xminq = Xminq(self)
      freq = self.freq_for_plot();
      Sxminq = self.Sxminq(freq);
      Xminq = sqrt(trapz(freq, Sxminq.^2));
    end
    
    
    
    
    function Fminv_cumulative = Fminv_cumulative(self)
      freq = self.freq_for_plot();
      Sfminv = self.Sfminv(freq);
      Fminv_cumulative = sqrt(cumtrapz(freq, Sfminv.^2));
    end
    
    function Xminv_cumulative = Xminv_cumulative(self)
      freq = self.freq_for_plot();
      Sxminv = self.Sxminv(freq);
      Xminv_cumulative = sqrt(cumtrapz(freq, Sxminv.^2));
    end
    
    function Fminq_cumulative = Fminq_cumulative(self)
      freq = self.freq_for_plot();
      Sfminq = self.Sfminq(freq);
      Fminq_cumulative = sqrt(cumtrapz(freq, Sfminq.^2));
    end
    
    function Xminq_cumulative = Xminq_cumulative(self)
      freq = self.freq_for_plot();
      Sxminq = self.Sxminq(freq);
      Xminq_cumulative = sqrt(cumtrapz(freq, Sxminq.^2));
    end
    
    % ==========================================================================
    % ========================== Plotting functions ============================
    % ==========================================================================
    
    function freq = freq_for_plot(self)
      freq = logspace(log10(self.freq_min), log10(self.freq_max), self.n_plot_points);
    end
    
    function plotStrain(self, force)
      [z, normalizedStrain] = self.normalized_strain_profile();
      
      strain = force*normalizedStrain;
      
      figure
      plot(z, strain, 'LineWidth', 2);
      xlabel('Z (m)', 'FontSize', 16);
      ylabel('Strain', 'FontSize', 16);
      
      stress = [self.ESi*strain(z <= self.t_si) self.EPE*strain(z >= self.t_si)];
      
      figure
      plot(z, stress, 'LineWidth', 2);
      xlabel('Z (m)', 'FontSize', 16);
      ylabel('Stress (Pa)', 'FontSize', 16);
    end
    
    function plotZ(self)
      freq = self.freq_for_plot();
      
      Z = self.Z(freq);
      
      hold all
      plot(freq, abs(Z), 'LineWidth', 2);
      plot(freq, real(Z), 'LineWidth', 2);
      hold off
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('Z (Ohms)', 'FontSize', 16);
      set(gca, 'LineWidth', 2, 'FontSize', 14);
      legend('|Z|', 'Re(Z)');
      
      set(gca, 'xscale', 'log', 'yscale', 'log');
    end
    
    function p = plotVn(self, lineStyle, color)
      freq = self.freq_for_plot();
      Vn = self.Vn(freq);
      
      p = plot(freq, Vn, 'LineWidth', 2, 'Linestyle', lineStyle, 'Color', color);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('Voltage Noise PSD (V-Hz^{-0.5})', 'FontSize', 16);
      set(gca, 'LineWidth', 2, 'FontSize', 14);
      set(gca, 'xscale', 'log', 'yscale', 'log');
      xlim([min(freq) max(freq)])
    end
    
    function p = plotQn(self, lineStyle, color)
      freq = self.freq_for_plot();
      Qn = self.Qn(freq);
      
      p = plot(freq, Qn, 'LineWidth', 2, 'Linestyle', lineStyle, 'Color', color);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('Charge Noise PSD (C-Hz^{-0.5})', 'FontSize', 16);
      set(gca, 'LineWidth', 2, 'FontSize', 14);
      set(gca, 'xscale', 'log', 'yscale', 'log');
      xlim([min(freq) max(freq)])
    end
    
    function plotSfv(self)
      freq = self.freq_for_plot();
      Sfv = self.v_f_sensitivity(freq);
      
      plot(freq, Sfv, 'LineWidth', 2, 'Linestyle', lineStyle);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('Voltage Sensitivity (V/N)', 'FontSize', 16);
      set(gca, 'xscale', 'log', 'yscale', 'log');
      set(gca, 'LineWidth', 2, 'FontSize', 14);
    end
    
    function plotSxv(self, lineStyle)
      freq = self.freq_for_plot();
      Sxv = self.v_x_sensitivity(freq);
      
      plot(freq, Sxv, 'LineWidth', 2, 'Linestyle', lineStyle);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('Displacement Sensitivity (V/m)', 'FontSize', 16);
      set(gca, 'xscale', 'log', 'yscale', 'log');
      set(gca, 'LineWidth', 2, 'FontSize', 14);
    end
    
    function plotSfq(self, lineStyle)
      freq = self.freq_for_plot();
      Sfq = self.q_f_sensitivity(freq);
      
      plot(freq, Sfq, 'LineWidth', 2, 'Linestyle', lineStyle);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('Charge Sensitivity (C/N)', 'FontSize', 16);
      set(gca, 'xscale', 'log', 'yscale', 'log');
      set(gca, 'LineWidth', 2, 'FontSize', 14);
    end
    
    function plotSxq(self, lineStyle)
      freq = self.freq_for_plot();
      Sxq = self.q_x_sensitivity(freq);
      
      plot(freq, Sxq, 'LineWidth', 2, 'Linestyle', lineStyle);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('Charge Sensitivity (C/m)', 'FontSize', 16);
      set(gca, 'xscale', 'log', 'yscale', 'log');
      set(gca, 'LineWidth', 2, 'FontSize', 14);
    end
    
    
    function plotSfminV(self, lineStyle)
      freq = self.freq_for_plot();
      SfminV = self.Sfminv(freq);
      
      plot(freq, SfminV, 'LineWidth', 2, 'Linestyle', lineStyle);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('Force Noise (N/rtHz)', 'FontSize', 16);
      set(gca, 'xscale', 'log', 'yscale' ,'log');
      set(gca, 'LineWidth', 2, 'FontSize', 14);
    end
    
    function plotSxminV(self, lineStyle)
      freq = self.freq_for_plot();
      SxminV = self.Sxminv(freq);
      
      plot(freq, SxminV, 'LineWidth', 2, 'Linestyle', lineStyle);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('Displacement Noise (m/rtHz)', 'FontSize', 16);
      set(gca, 'xscale', 'log', 'yscale' ,'log');
      set(gca, 'LineWidth', 2, 'FontSize', 14);
    end
    
    
    function plotSfminQ(self, lineStyle)
      freq = self.freq_for_plot();
      Sfminq = self.Sfminq(freq);
      
      plot(freq, Sfminq, 'LineWidth', 2, 'Linestyle', lineStyle);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('Force Noise (N/rtHz)', 'FontSize', 16);
      set(gca, 'xscale', 'log', 'yscale' ,'log');
      set(gca, 'LineWidth', 2, 'FontSize', 14);
    end
    
    function plotSxminQ(self, lineStyle)
      freq = self.freq_for_plot();
      Sxminq = self.Sxminq(freq);
      
      plot(freq, Sxminq, 'LineWidth', 2, 'Linestyle', lineStyle);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('Displacement Noise (m/rtHz)', 'FontSize', 16);
      set(gca, 'xscale', 'log', 'yscale' ,'log');
      set(gca, 'LineWidth', 2, 'FontSize', 14);
    end
    
    
    function plotFminV(self, lineStyle)
      freq = self.freq_for_plot();
      Fmin = self.Fminv_cumulative();
      
      plot(freq, Fmin, 'LineWidth', 2, 'Linestyle', lineStyle);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('MDF (N)', 'FontSize', 16);
      set(gca, 'xscale', 'log', 'yscale' ,'log');
      set(gca, 'LineWidth', 2, 'FontSize', 14);
    end
    
    function plotXminV(self, lineStyle)
      freq = self.freq_for_plot();
      Xmin = self.Xminv_cumulative();
      
      plot(freq, Xmin, 'LineWidth', 2, 'Linestyle', lineStyle);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('MDD (m)', 'FontSize', 16);
      set(gca, 'xscale', 'log', 'yscale' ,'log');
      set(gca, 'LineWidth', 2, 'FontSize', 14);
    end
    
    
    function plotFminQ(self, lineStyle)
      freq = self.freq_for_plot();
      Fmin = self.Fminq_cumulative();
      
      plot(freq, Fmin, 'LineWidth', 2, 'Linestyle', lineStyle);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('MDF (N)', 'FontSize', 16);
      set(gca, 'xscale', 'log', 'yscale' ,'log');
      set(gca, 'LineWidth', 2, 'FontSize', 14);
    end
    
    function plotXminQ(self, lineStyle)
      freq = self.freq_for_plot();
      Xmin = self.Xminq_cumulative();
      
      plot(freq, Xmin, 'LineWidth', 2, 'Linestyle', lineStyle);
      xlabel('Frequency (Hz)', 'FontSize', 16);
      ylabel('MDD (m)', 'FontSize', 16);
      set(gca, 'xscale', 'log', 'yscale' ,'log');
      set(gca, 'LineWidth', 2, 'FontSize', 14);
    end
    
    % =================== Damping calculations ========================
    
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
    
    % Calculate the quality factor for a given cantilever design assuming atmospheric pressure in air
    % Implemented per "Dependence of the quality factor of micromachined silicon beam resonators
    % on pressure and geometry" by Blom, 1992
    function Q = calculateQ(self)
      k_0 = 1.875; % first resonant mode factor
      [rho_f eta_f] = self.lookupFluidProperties();
      
      omega_0 = 2*pi*self.resonant_frequency();
      R = sqrt(self.w_si*self.l_si/pi); % effective sphere radius
      delta = sqrt(2*eta_f/rho_f/omega_0); % boundary layer thickness
      
      Q = k_0^2/(12*pi*sqrt(3))*sqrt(self.rhoSi*self.ESi)*self.w_si*self.t_si^2/(self.l_si*R*(1+R/delta)*eta_f);
    end
    
    
    % Calculate the damped natural frequency and Q, which we know lies between zero and the natural frequency in vacuum.
    % Per Eysden and Sader (2007)
    function [omega_damped, Q] = omega_damped_and_Q(self)
      
      % If we're in vacuum, just return the vacuum frequency
      if strcmp(self.fluid, 'vacuum')
        omega_damped = 2*pi*self.resonant_frequency();
        Q = PECantilever.maxQ;
        return;
      end
      
      % Inner function for solving the transcendental equation to find omega_damped
      % We're searching for a function minimum, so return the residual squared (continuous and smooth)
      function residual_squared = find_natural_frequency(omega_damped)
        hydro = self.hydrodynamic_function(omega_damped, rho_f, eta_f);
        residual = omega_damped - 2*pi*self.resonant_frequency()*(1 + pi * rho_f * self.w_si/(4 * self.rhoSi * self.t_si) .* real(hydro)).^-0.5;
        residual_squared = residual^2;
      end
      
      % Lookup fluid properties once, then calculate omega_damped and Q
      [rho_f eta_f] = self.lookupFluidProperties();
      options = optimset('TolX', 1, 'TolFun', 1e-6, 'Display', 'off'); % Calculate omega_damped to within 1 radian/sec
      omega_damped = fminbnd(@find_natural_frequency, 0, 2*pi*self.resonant_frequency(), options);
      hydro = self.hydrodynamic_function(omega_damped, rho_f, eta_f);
      Q = (4 * self.rhoSi * self.t_si / (pi * rho_f * self.w_si) + real(hydro)) / imag(hydro);
      
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
    
    % Calculate the Reynold's number
    function reynolds = reynolds(self, omega, rho_f, eta_f)
      reynolds = (rho_f*omega*self.w_si^2)/eta_f;
    end
    
    % Calculate kappa for the first mode
    function kappa = kappa(self)
      C = 1.8751;
      kappa = C * self.w_si / self.l_si;
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
    
    
    
    % ==========================================================================
    % =================== Optimization Helper Functions ========================
    % ==========================================================================
    
    % Used by optimization to bring all state varibles to O(1)
    function scaling = optimization_scaling(self)
      l_scale = 1e6;
      w_scale = 1e6;
      t_scale = 1e9;
      r_shunt_scale = 1e-9;
      scaling = [l_scale w_scale t_scale t_scale r_shunt_scale];
    end
    
    function new_cantilever = cantilever_from_state(self, x0)
      scaling = self.optimization_scaling();
      x0 = x0 ./ scaling;
      
      l_si = x0(1);
      w_si = x0(2);
      t_si = x0(3);
      t_pe = x0(4);
      r_shunt = x0(5);
      
      self.l_si = l_si;
      self.w_si = w_si;
      self.t_si = t_si;
      self.l_pe = l_si;
      self.w_pe = w_si;
      self.t_pe = t_pe;
      self.r_shunt = r_shunt;
      
      new_cantilever = self;
    end
    
    % Return state vector for the current state
    function x = current_state(self)
      x(1) = self.l_si;
      x(2) = self.w_si;
      x(3) = self.t_si;
      x(4) = self.t_pe;
      x(5) = self.r_shunt;
    end
    
    % ==========================================================================
    % ================ Optimization Goals and Constraints ======================
    % ==========================================================================
    
    function force_resolution = voltage_force_resolution(self, x0)
      c_new = self.cantilever_from_state(x0);
      force_resolution = c_new.Fminv()*1e12;
    end
    
    function displacement_resolution = voltage_displacement_resolution(self, x0)
      c_new = self.cantilever_from_state(x0);
      displacement_resolution = c_new.Xminv()*1e9;
    end
    
    function force_resolution = charge_force_resolution(self, x0)
      c_new = self.cantilever_from_state(x0);
      force_resolution = c_new.Fminq()*1e12;
    end
    
    function displacement_resolution = charge_displacement_resolution(self, x0)
      c_new = self.cantilever_from_state(x0);
      displacement_resolution = c_new.Xminq()*1e9;
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
      C(1) = -c_new.Fminv();
      
      % Natural frequency
      if exist('omega_min_hz', 'var')
        freq_constraint = omega_min_hz - c_new.resonant_frequency();
      end
      C = [C freq_constraint];
      
      if exist('min_k', 'var')
        min_k_constraint = min_k - c_new.stiffness();
        C = [C min_k_constraint];
      end
      
      if exist('max_k', 'var')
        max_k_constraint = c_new.stiffness() - max_k;
        C = [C max_k_constraint];
      end
      
      if exist('min_l_w_ratio', 'var')
        length_width_ratio = min_l_w_ratio - c_new.l_si/c_new.w_si;
        C = [C length_width_ratio];
      end
      
      if exist('min_w_t_ratio', 'var')
        width_thickness_ratio = min_w_t_ratio - c_new.w_si/c_new.t_si;
        C = [C width_thickness_ratio];
      end
      
      % Now for equality based constraints
      Ceq = [];
      
      if exist('fixed_k', 'var')
        fixed_k_constraint = c_new.stiffness() - fixed_k;
        Ceq = [Ceq fixed_k_constraint];
      end
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
      max_t = 10e-6;
      
      min_t_pe = 200e-9;
      max_t_pe = 1e-6;
      
      min_r_shunt = 1e12;
      max_r_shunt = 1e12;
      
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
      
      lb = [min_l, min_w, min_t, min_t_pe min_r_shunt];
      ub = [max_l, max_w, max_t, max_t_pe max_r_shunt];
    end
    
    function x0 = initial_conditions_random(self, parameter_constraints)
      [lb, ub] = self.optimization_bounds(parameter_constraints);
      
      min_l = lb(1);
      max_l = ub(1);
      
      min_w = lb(2);
      max_w = ub(2);
      
      min_t = lb(3);
      max_t = ub(3);
      
      min_t_pe = lb(4);
      max_t_pe = ub(4);
      
      min_r_shunt = lb(5);
      max_r_shunt = ub(5);
      
      % Generate the random values
      l = min_l + rand*(max_l - min_l);
      w = min_w + rand*(max_w - min_w);
      t = min_t + rand*(max_t - min_t);
      t_pe = min_t_pe + rand*(max_t_pe - min_t_pe);
      r_shunt = 10^(log10(min_r_shunt) + rand*(log10(max_r_shunt) - log10(min_r_shunt))); % logarithmically distributed
      
      x0 = [l, w, t, t_pe, r_shunt];
    end
    
    
    
    
    % ==========================================================================
    % ==================== Core Optimization Functions =========================
    % ==========================================================================
    
    
    % This cantilever problem isn't guaranteed to converge, and in
    % practice it fails to converge about 1% of the time for random
    % initial conditions. For this reason, it is best to start from a
    % random initial seed and perform the optimization and checking to
    % make sure that it converges repeatedly.
    function optimized_cantilever = optimize_performance(self, parameter_constraints, nonlinear_constraints, goal)
      
      percent_match = 0.001;
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
          if goal == self.goalVoltageForce
            resolution(jj) = c{jj}.Fminv();
          elseif goal == self.goalVoltageDisplacement
            resolution(jj) = c{jj}.Xminv();
          elseif goal == self.goalChargeForce
            resolution(jj) = c{jj}.Fminq();
          elseif goal == self.goalChargeDisplacement
            resolution(jj) = c{jj}.Xminq();
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
    function [optimized_cantilever, exitflag] = optimize_performance_from_current(self, parameter_constraints, nonlinear_constraints, goal)
      randomize_starting_conditions = 0;
      [optimized_cantilever, exitflag] = self.optimize_performance_once(parameter_constraints, nonlinear_constraints, goal, randomize_starting_conditions);
    end
    
    function [optimized_cantilever, exitflag] = optimize_performance_once(self, parameter_constraints, nonlinear_constraints, goal, random_flag)
      
      % Load the scaling vector
      scaling = self.optimization_scaling();
      
      % If random_flag = 1, start from random conditions. Otherwise
      % start from the current cantilever state vector
      if random_flag == 1
        problem.x0 = scaling.*self.initial_conditions_random(parameter_constraints);
      else
        problem.x0 = scaling.*self.current_state();
      end

      % To ensure that we don't have a bad initial point
      check_lever = self.cantilever_from_state(problem.x0);
      check_lever.print_performance()
      
      if goal == self.goalVoltageForce
        problem.objective = @self.voltage_force_resolution;
      elseif goal == self.goalVoltageDisplacement
        problem.objective = @self.voltage_displacement_resolution;
      elseif goal == self.goalChargeForce
        problem.objective = @self.charge_force_resolution;
      elseif goal == self.goalChargeDisplacement
        problem.objective = @self.charge_displacement_resolution;
      end
      
      [lb ub] = self.optimization_bounds(parameter_constraints);
      problem.lb = scaling.*lb;
      problem.ub = scaling.*ub;
      
      problem.options.TolFun = 1e-10;
      problem.options.TolCon = 1e-10;
      problem.options.TolX = 1e-10;
      
      problem.options.MaxFunEvals = 3000;
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
