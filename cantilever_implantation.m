classdef cantilever_implantation < cantilever
  properties
    implantation_energy; % 20 - 80 keV
    implantation_dose; % 2e14 - 2e16 per sq cm
    annealing_temp; % 900 - 1100C
    annealing_time; % 15 - 120min
    annealing_type; % 'inert' or 'oxide' -> inert only or regrow 1500A passivation oxide
    lookupTableData
  end
  
  methods
    % Call superclass constructor
    function self = cantilever_implantation(freq_min, freq_max, ...
        l, w, t, l_pr_ratio, v_bridge, doping_type, ...
        annealing_time, annealing_temp, annealing_type, ...
        implantation_energy, implantation_dose)
      
      self = self@cantilever(freq_min, freq_max, l, w, t, l_pr_ratio, ...
        v_bridge, doping_type);
      
      self.implantation_energy = implantation_energy;
      self.implantation_dose = implantation_dose;
      self.annealing_type = annealing_type;
      self.annealing_temp = annealing_temp;
      self.annealing_time = annealing_time;

      % Load the lookup table
      try
        data = open('lookupTable.mat');
        self.lookupTableData = data;
      catch
        fprintf('Dopant type: %s\n', self.doping_type);
        fprintf('Anneal type: %s\n', self.annealing_type);
        fprintf('Not available\n');
        return;
      end
    end
    
    
    function print_performance(self)
      print_performance@cantilever(self); % print the base stuff
      fprintf('Implantation energy (keV), dose (cm-2): %.1f %1g \n', self.implantation_energy, self.implantation_dose);
      fprintf('Annealing time (mins), temp (C): %.1f %.1f \n', self.annealing_time/60, self.annealing_temp-273);
      fprintf('Diffusion length (cm): %1g \n', self.diffusion_length());
      fprintf('Junction depth (um): %.3f \n', self.junction_depth()*1e6);
      fprintf('Alpha: %3g \n', self.alpha());
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
    function annealNumber = annealNumber(self)
      switch self.annealing_type
        case 'inert'
          annealNumber = 1;
        case 'oxide'
          annealNumber = 2;
        otherwise
          fprintf('Unknown anneal condition: %s\n', self.annealing_type);
          pause
      end
    end
    
    % Lookup the concentration profile from the lookup table
    function [x, active_doping, total_doping] = doping_profile(self)
      x = self.lookupTableData.z; % 10nm spacing from 0 to 5um

      n = interpn(x, self.lookupTableData.ImplantDopants, self.lookupTableData.ImplantDoses, ...
        self.lookupTableData.ImplantEnergies, self.lookupTableData.AnnealTemps, ...
        self.lookupTableData.AnnealTimes, self.lookupTableData.AnnealOxidation, ...
        self.lookupTableData.n, ...
        x, self.dopantNumber(), self.implantation_dose, self.implantation_energy, ...
        self.annealing_temp-273, self.annealing_time/60, self.annealNumber(), 'linear');
      
      % Remove data beyond the device thickness
      % We don't need to remove data beyond the junction because beta/Rs/Nz
      % are already just calculated to the junction in the lookup table
      n(x > self.t) = [];
      x(x > self.t) = [];
      
      % Active = total unless the doping is higher than the solid solubility limit
      % which is generally not the case in the ion implantation data
      active_doping = n;
      total_doping = n; 
    end
    
    % Junction depth. Units = m
    function Xj = junction_depth(self)
      Xj = interpn(self.lookupTableData.ImplantDopants, self.lookupTableData.ImplantDoses, ...
        self.lookupTableData.ImplantEnergies, self.lookupTableData.AnnealTemps, ...
        self.lookupTableData.AnnealTimes, self.lookupTableData.AnnealOxidation, ...
        self.lookupTableData.Xj, ...
        self.dopantNumber(), self.implantation_dose, self.implantation_energy, self.annealing_temp-273, ...
        self.annealing_time/60, self.annealNumber(), 'linear');
    end
    
    function Rs = sheet_resistance(self)
      Rs = interpn(self.lookupTableData.ImplantDopants, self.lookupTableData.ImplantDoses, ...
        self.lookupTableData.ImplantEnergies, self.lookupTableData.AnnealTemps, ...
        self.lookupTableData.AnnealTimes, self.lookupTableData.AnnealOxidation, ...
        self.lookupTableData.Rs, ...
        self.dopantNumber(), self.implantation_dose, self.implantation_energy, self.annealing_temp-273, ...
        self.annealing_time/60, self.annealNumber(), 'linear');
    end
    
    % Nz is an abstract function specifically for this class
    % In the other classes (epitaxy, diffusion) we calculate Nz from the profile
    % Note that Nz != Nz_total due to current crowding
    function Nz = Nz(self)
      Nz = interpn(self.lookupTableData.ImplantDopants, self.lookupTableData.ImplantDoses, ...
        self.lookupTableData.ImplantEnergies, self.lookupTableData.AnnealTemps, ...
        self.lookupTableData.AnnealTimes, self.lookupTableData.AnnealOxidation, ...
        self.lookupTableData.Nz, ...
        self.dopantNumber(), self.implantation_dose, self.implantation_energy, self.annealing_temp-273, ...
        self.annealing_time/60, self.annealNumber(), 'linear');
    end
    
    function beta = beta(self)
      Beta1 = interpn(self.lookupTableData.ImplantDopants, self.lookupTableData.ImplantDoses, ...
        self.lookupTableData.ImplantEnergies, self.lookupTableData.AnnealTemps, ...
        self.lookupTableData.AnnealTimes, self.lookupTableData.AnnealOxidation, ...
        self.lookupTableData.Beta1, ...
        self.dopantNumber(), self.implantation_dose, self.implantation_energy, self.annealing_temp-273, ...
        self.annealing_time/60, self.annealNumber(), 'linear');
      
      Beta2 = interpn(self.lookupTableData.ImplantDopants, self.lookupTableData.ImplantDoses, ...
        self.lookupTableData.ImplantEnergies, self.lookupTableData.AnnealTemps, ...
        self.lookupTableData.AnnealTimes, self.lookupTableData.AnnealOxidation, ...
        self.lookupTableData.Beta2, ...
        self.dopantNumber(), self.implantation_dose, self.implantation_energy, self.annealing_temp-273, ...
        self.annealing_time/60, self.annealNumber(), 'linear');
      
      beta = Beta1 - 2/(self.t*1e6)*Beta2; % Beta2 is in units of microns, so convert t
      
%       self.lookupTableData.ImplantDopants
%       self.dopantNumber()
% 
%       self.lookupTableData.ImplantDoses
%       self.implantation_dose
% 
%       self.lookupTableData.ImplantEnergies
%       self.implantation_energy
% 
%       self.lookupTableData.AnnealTemps
%       self.annealing_temp-273
% 
%       self.lookupTableData.AnnealTimes
%       self.annealing_time/60
% 
%       self.lookupTableData.AnnealOxidation
%       self.annealNumber()
    end
    
    % Note that this is only an approximate diffusion length
    % TODO: Is the data for the dopant or silicon diffusion length?
    % TODO: Does it hold for dopants besides boron?
    function diffusion_length = diffusion_length(self)
      switch self.doping_type
        case 'arsenic'
          D0 = 22.9; %cm2/s
          Ea = 4.1; %eV
        case 'boron'
          D0 = 0.76; %cm2/s
          Ea = 3.46; %eV
        case 'phosphorus'
          D0 = 3.85; %cm2/s
          Ea = 3.66; %eV
      end
      diffusivity = D0*exp(-1*Ea/self.k_b_eV/self.annealing_temp);
      diffusion_length = (diffusivity*self.annealing_time)^0.5; %cm
    end
    
    % Calculate alpha from sqrt(Dt) using the data I compiled for the book
    function alpha = alpha(self)
      alpha = 2.469e-10*self.diffusion_length()^-0.598;
    end
    
    % ==================================
    % ========= Optimization  ==========
    % ==================================
    
    function scaling = doping_optimization_scaling(self)
      annealing_time_scale = 1e-3;
      annealing_temp_scale = 1e-3;
      implantation_energy_scale = 1;
      implantation_dose_scale = 1e-14;
      scaling = [annealing_time_scale, annealing_temp_scale, ...
        implantation_energy_scale, implantation_dose_scale];
    end
    
    function self = doping_cantilever_from_state(self, x0)
      self.annealing_time = x0(6);
      self.annealing_temp = x0(7);
      self.implantation_energy = x0(8);
      self.implantation_dose = x0(9);
    end
    
    function x = doping_current_state(self)
      x(1) = self.annealing_time;
      x(2) = self.annealing_temp;
      x(3) = self.implantation_energy;
      x(4) = self.implantation_dose;
    end
    
    % Set the minimum and maximum bounds for the cantilever state
    % variables. Bounds are written in terms of the initialization
    % variables. Secondary constraints (e.g. power dissipation,
    % piezoresistor thickness rather than ratio, resonant frequency)
    % are applied in optimization_constraints()
    function [lb ub] = doping_optimization_bounds(self, parameter_constraints)
      min_annealing_time = 15*60; % seconds
      max_annealing_time = 120*60;
      
      min_annealing_temp = 273+900; % K
      max_annealing_temp = 273+1100;
      
      min_implantation_energy = 20; %keV
      max_implantation_energy = 80;
      
      min_implantation_dose = 2e14; %cm-2
      max_implantation_dose = 2e16;
      
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
      
      lb = [min_annealing_time, min_annealing_temp min_implantation_energy min_implantation_dose];
      ub = [max_annealing_time, max_annealing_temp max_implantation_energy max_implantation_dose];
    end
    
    function x0 = doping_initial_conditions_random(self, parameter_constraints)
      [lb, ub] = self.doping_optimization_bounds(parameter_constraints);
      
      annealing_time_min = lb(1);
      annealing_time_max = ub(1);
      
      annealing_temp_min = lb(2);
      annealing_temp_max = ub(2);
      
      implantation_energy_min = lb(3);
      implantation_energy_max = ub(3);
      
      implantation_dose_min = lb(4);
      implantation_dose_max = ub(4);
      
      % Generate the random values
      annealing_time_random = annealing_time_min + rand*(annealing_time_max - annealing_time_min);
      annealing_temp_random = annealing_temp_min + rand*(annealing_temp_max - annealing_temp_min);
      implantation_energy_random = implantation_energy_min + rand*(implantation_energy_max - implantation_energy_min);
      implantation_dose_random = implantation_dose_min + rand*(implantation_dose_max - implantation_dose_min);
      
      x0 = [annealing_time_random, annealing_temp_random, implantation_energy_random, implantation_dose_random];
    end
  end
end