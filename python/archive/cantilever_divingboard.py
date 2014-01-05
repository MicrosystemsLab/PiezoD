from numpy import *
from scipy import *
from scipy.linalg import *
from pylab import *
# import matplotlib.axes3d as p3

from scikits.openopt import NLP

# Note: All units are MKS unless otherwise noted
# Piezoresistor dimensions are specified in percentages of total dimensions
# Ex: t = 10e-6, t_pr = 0.5 => 5 micron thick piezoresistor
class cantilever_divingboard:
    def __init__(self, 
                freq_min = 1e2, freq_max = 100e3,
                l = 100e-6, w = 20e-6, t = 10e-6,
                l_pr_ratio = 0.5, w_pr_ratio = 0.5, w_gap_ratio = 0.5, t_pr_ratio = 0.5,
                V_bias = 5., N = 4e19,
                E = 169.e9, rho_c = 2330., rho_l = 1000., eta = 1.e-3, beta = 1.):
        
        # Cantilever geometry
        self.l = l
        self.w = w
        self.t = t
        
        # Piezo geometry
        self.l_pr_ratio = l_pr_ratio
        self.w_pr_ratio = w_pr_ratio
        self.w_gap_ratio = w_gap_ratio
        self.t_pr_ratio = t_pr_ratio
        
        # Mechanical properties
        self.E = E
        
        # Fluid properties
        self.rho_c = rho_c
        self.rho_l = rho_l
        self.eta = eta
        
        # Electrical properties
        self.V_bias = V_bias
        self.N = N
        self.freq_max = freq_max
        self.freq_min = freq_min
        self.beta_measured = beta # used for passing in values of beta
        
        # Material constants
        self.alpha = 1.e-5 # questionable
        self.T = 300 # Temperature (K)
        self.k_b = 1.38e-23 # Boltzmann's constant
        self.q = 1.60218e-19 # electron charge (C)
        self.rho_aluminum = 2.7e-6 # Resistivity of aluminum (ohm-m)
        
    
    # Calculate the actual piezo dimensions
    def l_pr(self):
        return self.l*self.l_pr_ratio
        
    def w_pr(self):
        return self.w*self.w_pr_ratio
        
    def w_gap(self):
        return self.w*self.w_gap_ratio
        
    def t_pr(self):
        return self.t*self.t_pr_ratio
    
    
    
    # Calculate the resistance (ohms)
    # Units: V-sec/C
    def resistance(self):
        resistivity_in_meters = self.resistivity()*1e-2
        return resistivity_in_meters*self.resistor_length()/self.resistor_cross_section()
    
    def resistor_cross_section(self):
        return self.w_pr()*self.t_pr()
    
    # Print out the key cantilever properties
    def print_performance(self):
        print 'Cantilever L/W/T: %f %f %f' % (self.l*1e6, self.w*1e6, self.t*1e6)
        print 'Piezo L/W/G/T: %f %f %f %f' % (self.l_pr()*1e6, self.w_pr()*1e6, self.w_gap()*1e6, self.t_pr()*1e6)
        print 'Piezo Ratios L/W/G/T: %g %g %g %g' % (self.l_pr_ratio, self.w_pr_ratio, self.w_gap_ratio, self.t_pr_ratio)
        print 'Bias voltage: %f' % self.V_bias
        print 'Doping Concentration: %g' % self.N
        print ""
        print 'Resistance: %f' % self.resistance()
        print 'Power dissipation (mW): %g' % (self.power_dissipation()*1e3)
        print ""
        print 'Force resolution: %g' % self.force_resolution()
        print 'Sensitivity (V/N) %g' % self.force_sensitivity()
        print 'Piezo factor(m^2/N): %g' % self.piezoresistance_factor()
        print 'Beta %g' % self.beta()
        print ""
        print 'Integrated noise (V): %g' % self.integrated_noise()
        print 'Integrated johnson noise (V): %g' % self.integrated_johnson_noise()
        print 'Integrated 1/f noise (V): %g' % self.integrated_hooge_noise()
        print 'Knee frequency (Hz): %g' % self.knee_frequency()
        print ""
        print 'Stiffness (N/m): %g' % self.stiffness()
        print 'Damped freq: %f' % self.omega_damped_hz()
        print 'Vacuum freq: %f' % self.omega_vacuum_hz()
        print 'Freq range: %f to %f' % (self.freq_min, self.freq_max)        
        print ""
    
    
    # Mobility of majority carriers in silicon for As, P or B.
    # Units: cm^2/V-sec
    def mobility(self):
        # Arsenic data
        umin_As = 52.2
        umax_As = 1417
        Nr_As = 9.68e16
        alpha_As = 0.68

        # Phosphorous data
        umin_P = 68.5
        umax_P = 1414
        Nr_P = 9.20e16
        alpha_P = 0.711

        # Boron data
        umin_B = 44.9
        umax_B = 470.5
        Nr_B = 2.23e17
        alpha_B = 0.719

        # Boron
        umin = umin_B
        umax = umax_B
        Nr = Nr_B
        alpha = alpha_B

        mob = umin + (umax-umin)/(1+(self.N/Nr)**alpha)
        return mob
        
    # Calculate the resistivity
    # Units: 1/(1/cm^3 * cm^2/V-sec * C) = cm-V-sec/C
    def resistivity(self):
        return 1/(self.N*self.mobility()*self.q)
    
    
    # Calculate the effective length of the resistor
    # Units: m
    def resistor_length(self):
        return 2*self.l_pr() + self.w_gap()
    

    # Piezoresistance factor, which decreases with doping level
    def piezoresistance_factor(self):
        b = 1.53e22
        a = 0.2014
        P = log10((b/self.N) ** a)
        P_base = 72.e-11 # Pi at low concentration in 110 from Smith
        return P*P_base
    

    # 1/f noise density
    # Units: V^2
    def hooge_noise_density(self):
        return self.alpha*self.V_bias**2/self.number_of_carriers()

    # Johnson noise density
    # Units: V^2/Hz
    def johnson_noise_density(self):
        return 4*self.k_b*self.T*self.resistance()

    # V
    def integrated_johnson_noise(self):
        return sqrt(self.johnson_noise_density()*(self.freq_max - self.freq_min))

    # V
    def integrated_hooge_noise(self):
        return sqrt(self.alpha * self.V_bias**2 / self.number_of_carriers() * log(self.freq_max / self.freq_min))
    
    # Assume noise of 2nV/rtHz
    def integrated_amplifier_noise(self):
        return 2e-9 * sqrt(self.freq_max - self.freq_min)
        
    # Johnson noise voltage
    # Units = V/rtHz
    def johnson_noise_voltage(self):
        return sqrt(self.johnson_noise_density())

    # Power dissipation (W)
    def power_dissipation(self):
        current = self.V_bias / self.resistance()
        return current**2 * self.resistance()

    # Calculate the noise in V/rtHz at a given frequency
    def voltage_noise(self, frequency):
        hooge_noise = sqrt(self.hooge_noise_density()/frequency) # V/rtHz
        johnson_noise = sqrt(self.johnson_noise_density()) # V/rtHz
        return sqrt(hooge_noise**2 + johnson_noise**2)
        
    def plot_noise_spectrum(self):
        frequency = linspace(self.freq_min, self.freq_max, 1000)
        noise = self.voltage_noise(frequency)
        loglog(frequency, noise)

    # Calculate the knee frequency
    # Equating 1/f noise and johnson leads to: f = alpha*V_bias^2/N/Johnson
    def knee_frequency(self):
        return self.alpha*self.V_bias**2/(self.number_of_carriers()*self.johnson_noise_density())
        

    # Integrated cantilever noise for given bandwidth
    # Source are uncorrelated so they add in RMS
    # Units: V = sqrt(V^2 + V^2/Hz*Hz)
    def integrated_noise(self):
        return sqrt(self.integrated_johnson_noise()**2 + 
                    self.integrated_hooge_noise()**2 +
                    self.integrated_amplifier_noise()**2)

    # Reduction in sensitivity from piezoresistor not located just at the surface
    def beta(self):
        if self.beta < 1:
            return self.beta_measured
        else:
            return (1 - self.t_pr_ratio)
        
    #Check to make sure that all of the values are kosher
    def check_if_valid_solution(self):
        if not(self.l > 0 and self.w > 0 and self.t > 0).all():
            raise NameError, 'Negative dimensions'
        elif not(self.l > self.l_pr() and self.w > self.w_pr() and self.t > self.t_pr()):
            raise NameError, 'The piezoresistor is larger than the cantilever'
        elif self.resistance() < 0:
            raise NameError, 'The resistance of the piezoresistor is negative'
        elif self.force_resolution() < 0:
            raise NameError, 'The force resolution is negative'
        else:
            return
    
    # Ratio of piezoresistor resistance to total resistance (< 1)
    def gamma(self):
        return 1
        
    # Units: V/N
    def force_sensitivity(self):
        return 3*(self.l - self.l_pr()/2)*self.piezoresistance_factor()/(2*self.w*self.t**2)*self.beta()*self.gamma()*self.V_bias
        
    # 
    # Units: N
    def force_resolution(self):
        return self.integrated_noise()/self.force_sensitivity()

    # The number of current carriers in the piezoresistor
    # Used for calculating 1/f noise
    # Units: unitless
    def number_of_carriers(self):
        volume = self.resistor_length()*self.resistor_cross_section() # cubic meters
        volume_cc = volume*(1e2)**3
        return volume_cc*self.N

    
    # Bending stiffness of the cantilever to a point load at the tip
    # Units: N/m
    def stiffness(self):
        return self.E*self.w*self.t**3/(4*self.l**3)
    

    # Resonant frequency for undamped vibration (first mode)
    # Units: radians/sec
    def omega_vacuum(self):
        m_effective = 0.243*self.rho_c*self.w*self.t*self.l
        omega_radians = sqrt( self.stiffness() / m_effective)
        return omega_radians
    

    # Resonant frequency for undamped vibration (first mode)
    # Units: cycles/sec
    def omega_vacuum_hz(self):
        return self.omega_vacuum() / (2*pi)
    

    # Calculate the damped natural frequency (assuming Q > 1) from Sader paper
    # Units: radians/sec
    def omega_damped(self):
        def find_freq(omega):
            tau_r = real( self.hydrodynamic_function(omega) )
            residual = omega - self.omega_vacuum() * (1. + pi * self.rho_l * self.w /(4 * self.rho_c * self.t) * tau_r) ** -0.5
            return residual
        
        omega_radians = optimize.fsolve( find_freq, self.omega_vacuum() )
        return omega_radians
    

    # Calculate the damped natural frequency (assuming Q > 1) from Sader paper
    # Units: cycles/sec
    def omega_damped_hz(self):
        return self.omega_damped() / (2*pi)
    

    # Reynold's number for flow around the cantilever for a given frequency (in rad/sec)
    # Units: Dimensionless
    def reynolds_number(self, omega):
        return (self.rho_l* omega * self.w**2) / self.eta
        
    

    # Normalized bending mode number. The shape parameter (C) is hte solution to 1+cos(x)*cosh(x) = 0
    # Currently just used with the first vibrational mode, but higher order modes could
    # be solved by finding higher zeros of the hyperbolic function.
    def kappa(self):
        def beam_modes(x):
            return 1. + cos(x)*cosh(x)
        
        C = optimize.fsolve(beam_modes, 1.)
        return C * self.w / self.l
    
    
    
    # Find hydrodynamic function for a cantilever with a given reynold's number
    # and normalized mode number (kappa). Data from Sader 2007 paper.
    # The data is based upon the log10 of the reynold's number.
    def hydrodynamic_function(self, omega):
        log_reynolds = log10( self.reynolds_number(omega) )
        kappa = self.kappa()
        
        kappa_lookup = array([0, 0.125, 0.25, 0.5, 0.75, 1, 2, 3, 5, 7, 10, 20])
        reynolds_lookup = array([-4, -3.5, -3, -2.5, -2, -1.5, -1, -.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4])
        
        tau_lookup_real = array([
        [3919.41, 59.3906, 22.4062, 9.13525, 5.62175,  4.05204,  1.93036,  1.2764,   0.764081, 0.545683, 0.381972, 0.190986],
        [1531.90, 59.3861, 22.4061, 9.13525, 5.62175,  4.05204,  1.93036,  1.2764,   0.764081, 0.545683, 0.381972, 0.190986],
        [613.426, 59.3420, 22.4052, 9.13523, 5.62174,  4.05204,  1.93036,  1.2764,   0.764081, 0.545683, 0.381972, 0.190986],
        [253.109, 58.9094, 22.3962, 9.13504, 5.62172,  4.05203,  1.93036,  1.2764,   0.764081, 0.545683, 0.381972, 0.190986],
        [108.429, 55.2882, 22.3078, 9.13319, 5.62153,  4.05199,  1.93036,  1.2764,   0.764081, 0.545683, 0.381972, 0.190986],
        [48.6978, 40.7883, 21.5187, 9.11481, 5.61960,  4.05160,  1.93035,  1.2764,   0.764081, 0.545683, 0.381972, 0.190986],
        [23.2075, 22.7968, 17.5378, 8.94370, 5.60057,  4.04771,  1.93027,  1.27639,  0.76408,  0.545683, 0.381972, 0.190986],
        [11.8958, 11.9511, 11.0719, 7.89716, 5.43378,  4.01051,  1.92942,  1.27629,  0.764074, 0.545682, 0.381972, 0.190986],
        [6.64352, 6.64381, 6.47227, 5.65652, 4.64017,  3.74600,  1.92114,  1.27536,  0.764012, 0.545671, 0.38197,  0.190986],
        [4.07692, 4.05940, 3.99256, 3.72963, 3.37543,  3.00498,  1.85532,  1.26646,  0.763397, 0.545564, 0.381953, 0.190985],
        [2.74983, 2.73389, 2.69368, 2.56390, 2.39884,  2.22434,  1.61821,  1.20592,  0.757637, 0.54452,  0.381786, 0.190981],
        [2.02267, 2.01080, 1.98331, 1.90040, 1.79834,  1.69086,  1.31175,  1.04626,  0.721165, 0.535593, 0.38018,  0.190932],
        [1.60630, 1.59745, 1.57723, 1.51690, 1.44276,  1.36416,  1.08036,  0.878177, 0.635443, 0.496169, 0.368548, 0.190459],
        [1.36230, 1.35532, 1.33934, 1.29142, 1.23203,  1.16842,  0.932812, 0.759965, 0.551349, 0.435586, 0.334279, 0.186672],
        [1.21727, 1.21141, 1.19792, 1.15718, 1.10624,  1.05117,  0.84292,  0.686229, 0.493924, 0.387183, 0.295972, 0.172722],
        [1.13038, 1.12518, 1.11316, 1.07668, 1.03073,  0.980721, 0.78879,  0.641744, 0.458699, 0.356289, 0.268907, 0.154450],
        [1.07814, 1.07334, 1.06221, 1.02827, 0.985314, 0.938346, 0.756309, 0.615164, 0.437743, 0.337813, 0.252327, 0.140852]])
        
        tau_lookup_imag = array([
        [27984.8,   44628.5,   55176.1,   71754,     86311.5,   100062,    152411,    203623,    305570,    407436,    560225,    1069521],
        [9816.33,   14113.1,   17448.2,   22690.6,   27294.1,   31642.3,   48196.5,   64391.4,   96629.7,   128843,    177159,    338212],
        [3482.47,   4464.16,   5517.72,   7175.41,   8631.15,   10006.2,   15241.1,   20362.3,   30557,     40743.6,   56022.5,   106952],
        [1252.66,   1415.42,   1745.17,   2269.09,   2729.42,   3164.23,   4819.65,   6439.14,   9662.97,   12884.3,   17715.9,   33821.2],
        [458.386,   457.863,   552.862,   717.635,   863.138,   1000.63,   1524.112,  2036.23,   3055.7,    4074.36,   5602.25,   10695.2],
        [171.397,   160.951,   177.702,   227.205,   273.013,   316.449,   481.967,   643.915,   966.297,   1288.43,   1771.59,   3382.12],
        [65.8679,   62.2225,   61.626,    72.6542,   86.5364,   100.144,   152.418,   203.625,   305.57,    407.436,   560.225,   1069.52],
        [26.2106,   25.21,     24.1432,   24.7484,   27.9459,   31.8957,   48.2199,   64.3973,   96.6308,  128.843,    177.159,   338.212],
        [10.8983,   10.6158,   10.1909,   9.7009,    9.91067,   10.648,    15.3139,   20.381,    30.5604,   40.7448,   56.0229,   106.952],
        [4.78389,   4.69492,   4.53952,   4.24925,   4.09701,   4.09433,   5.01844,   6.49605,   9.67379,   12.8879,   17.7171,   33.8214],
        [2.23883,   2.20681,   2.14583,   2.0088,    1.89659,   1.82463,   1.85993,   2.17718,   3.08849,   4.08581,   5.60598,   10.6956],
        [1.12164,   1.10851,   1.08208,   1.01654,   0.953355,  0.901676,  0.81464,   0.844519,  1.04394,   1.32116,   1.78306,   3.38349],
        [0.596697,  0.590686,  0.578118,  0.545082,  0.510467,  0.479247,  0.403803,  0.383595,  0.409256,  0.469688,  0.589749,  1.07377],
        [0.332285,  0.329276,  0.32283,   0.305262,  0.285953,  0.26763,   0.216732,  0.194409,  0.186218,  0.195634,  0.221631,  0.349855],
        [0.191043,  0.189434,  0.185931,  0.176166,  0.165118,  0.154323,  0.122124,  0.105573,  0.0938839, 0.0925686, 0.09682,   0.126835],
        [0.112082,  0.111181,  0.109199,  0.103595,  0.0971392, 0.0907188, 0.0707736, 0.059728,  0.0505049, 0.0476557, 0.0471326, 0.0534759],
        [0.0665172, 0.0659974, 0.0648471, 0.0615627, 0.0577366, 0.0538889, 0.0416384, 0.0345727, 0.0282418, 0.025856,  0.024611,  0.0252877]])
    
        # Interpolate for our tau values using a bivariate spline (order kx/ky)
        x,y = meshgrid(kappa_lookup, reynolds_lookup)
        tau_real_grid = interpolate.RectBivariateSpline(kappa_lookup, reynolds_lookup, transpose(tau_lookup_real), kx=3, ky=3)
        tau_real = tau_real_grid(kappa, log_reynolds)
        tau_imag_grid = interpolate.RectBivariateSpline(kappa_lookup, reynolds_lookup, transpose(tau_lookup_imag), kx=3, ky=3)
        tau_imag = tau_imag_grid(kappa, log_reynolds)
        
        return complex(tau_real, tau_imag)
    
    # Quality factor
    def Q(self):
        hydro = self.hydrodynamic_function( self.omega_damped() )
        return (4*self.rho_c*self.t/(pi*self.rho_l*self.w) + hydro.real)/hydro.imag
        
        
        
    # ===== Optimization ======
    
    def optimize_cantilever_new(self, omega_min, min_dimensions, max_voltage, max_power, min_N, max_N, fluid_type):
        scale_factor = array([1e6, 1e6, 1e9, 1e6, 1e6, 1e6, 1e9, 1, 1e-18])
        def optimize_cantilever(x, omega_min, min_dimensions, max_voltage, max_power, min_N, max_N, fluid_type):
            x = x.copy()/scale_factor            
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x)
            # print c.force_resolution()        
            return c.force_resolution()
        
        def enforce_positive_resolution(x, omega_min, min_dimensions, max_voltage, max_power, min_N, max_N, fluid_type):
            x = x.copy()/scale_factor
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x)
            return c.force_resolution()
            
        def enforce_minimum_dimensions(x, omega_min, min_dimensions, max_voltage, max_power, min_N, max_N, fluid_type):
            x = x.copy()/scale_factor
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x)
            (min_length, min_width, min_thickness,
             min_pr_length, min_pr_width, min_pr_thickness) = min_dimensions

            return min(c.l - min_length, c.w - min_width, c.t - min_thickness)
        
        def enforce_positive_dimensions(x, omega_min, min_dimensions, max_voltage, max_power, min_N, max_N, fluid_type):
            x = x.copy()/scale_factor
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x)
            
            return min(c.l, c.w, c.t, c.l_pr(), c.w_pr(), c.w_gap(), c.t_pr())
        
        def enforce_maximum_piezo_dimensions(x, omega_min, min_dimensions, max_voltage, max_power, min_N, max_N, fluid_type):
            x = x.copy()/scale_factor
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x)
            return min(c.l - c.l_pr(), c.w - c.w_pr(), c.t - c.t_pr())
        
        def enforce_minimum_piezo_dimensions(x, omega_min, min_dimensions, max_voltage, max_power, min_N, max_N, fluid_type):
            x = x.copy()/scale_factor
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x)            
            (min_length, min_width, min_thickness,
             min_pr_length, min_pr_width, min_pr_thickness) = min_dimensions
            
            return min(c.l_pr() - min_pr_length, c.w_pr() - min_pr_width, c.t_pr() - min_pr_thickness)

        def enforce_minimum_natural_frequency(x, omega_min, min_dimensions, max_voltage, max_power, min_N, max_N, fluid_type):
            x = x.copy()/scale_factor
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x)
                        
            if fluid_type == 'vacuum':
                return c.omega_vacuum_hz() - omega_min
            elif fluid_type == 'water':
                return c.omega_damped_hz() - omega_min
            else:
                return c.omega_vacuum_hz()  
        
        def enforce_N_range(x, omega_min, min_dimensions, max_voltage, max_power, min_N, max_N, fluid_type):
            x = x.copy()/scale_factor
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x)
            if c.N > max_N:
                return max_N - c.N
            elif c.N < min_N:
                return c.N - min_N
            else:
                return c.N
        
        def enforce_voltage_range(x, omega_min, min_dimensions, max_voltage, max_power, min_N, max_N, fluid_type):
            x = x.copy()/scale_factor
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x)
            if c.V_bias > max_voltage:
                return max_voltage - c.V_bias
            else:
                return c.V_bias
                
        def enforce_max_power(x, omega_min, min_dimensions, max_voltage, max_power, min_N, max_N, fluid_type):
            x = x.copy()/scale_factor
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x)
            return max_power - c.power_dissipation()
        
        initial_guess = (self.l, self.w, self.t, 
                            self.l_pr_ratio, self.w_pr_ratio, self.w_gap_ratio, self.t_pr_ratio, 
                            self.V_bias, self.N)
        initial_guess *= scale_factor
        x = optimize.fmin_cobyla( func = optimize_cantilever, 
                                  x0 = initial_guess,
                                  args = (omega_min, min_dimensions, max_voltage, max_power, min_N, max_N, fluid_type),
                                  cons = (enforce_positive_resolution, 
                                          enforce_positive_dimensions,
                                          enforce_minimum_dimensions,
                                          enforce_maximum_piezo_dimensions,
                                          enforce_minimum_piezo_dimensions,
                                          enforce_minimum_natural_frequency,
                                          enforce_N_range,
                                          enforce_voltage_range,
                                          enforce_max_power),
                                  rhobeg = 1e-9, rhoend = 1e-12, maxfun = 10000, iprint = 1)
        x *= 1/scale_factor
        (self.l, self.w, self.t, 
              self.l_pr_ratio, self.w_pr_ratio, self.w_gap_ratio, self.t_pr_ratio, 
              self.V_bias, self.N) = x
        
    
    def optimize_cantilever(self, constraints, omega_min, fluid_type):
        """
        Optimize the force resolution of a piezoresistive cantilever over a given
        frequency band with the following assumptions:

        * the piezoresistive layer is fabricated by epitaxy i.e. it is a step function
        * no higher temperature processes (> 800C) are applied so sqrt(Dt) is small

        Keyword arguments:
        initial_guess -- a tuple of the initial cantilever state
        constraints -- a list of the form (min, max) for each state parameter.
        freq_min -- the minimum limit of the frequency range to integrate the signal over
        freq_max -- the maximum limit of the frequency range to integrate the signal over
        omega_min -- the minimum **damped** natural frequency

        Notes:
        minimum frequency -- approximate as the inverse of the measurement time
                             e.g. 10 seconds = 0.1Hz
        maximum frequency -- choose based upon the fastest signal that you want to
                             be able to capture. Remember to actually sample the
                             system at 2x the maximum frequency specified here
        frequency range -- noise is integrated over the entire frequency range
        constraints -- constraints are specified as upper and lower limits. If there
                        is no limit on one side you can say 'None'. Constraints can
                        be tricky with the piezoresistor dimensions, because if for
                        example the max limit of the piezoresistor thickness is
                        larger than the min limit of the cantilever thickness, you
                        could have a problem.
        general optimization behavior -- thin, narrow cantilevers result in the best
        """
        scale_factor = array([1e6, 1e6, 1e9, 1e-1, 1e-1, 1e-1, 1e-1, 1., 1e-18])

        def optimize_force_resolution(x0):
            x0 = x0.copy()/scale_factor
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x0)

            # Apply a penalty if the natural frequency is too low
            # penalty = 0
            # if fluid_type == 'vacuum':
            #     penalty += exp(omega_min - c.omega_vacuum_hz())
            # elif fluid_type == 'water':
            #     penalty += exp(omega_min - c.omega_damped_hz())
                
            # Apply a penalty if the piezoresistor is larger than the cantilever
            # if c.invalid_dimensions():
            #     penalty += c.l_pr - c.l
            #     penalty += c.w_pr - c.w
            #     penalty += c.t_pr - c.t
            #     penalty *= penalty*1e-3

            # print x0
            return c.force_resolution()

        # Scale our initial guess and constraints
        initial_guess = (self.l, self.w, self.t, 
                            self.l_pr_ratio, self.w_pr_ratio, self.w_gap_ratio, self.t_pr_ratio, 
                            self.V_bias, self.N)
        initial_guess *= scale_factor
        new_constraints = []
        for ii in range(0, len(constraints)):
            first = constraints[ii][0]*scale_factor[ii]
            second = constraints[ii][1]*scale_factor[ii]
            new_constraints.append((first, second))

        # Find the optimal solution and descale the result
        x, f, d = optimize.fmin_l_bfgs_b(func = optimize_force_resolution,
                                    x0 = initial_guess,
                                    bounds = new_constraints,
                                    approx_grad=1,
                                    m = 20,
                                    pgtol = 1e-8,
                                    factr = 1e6,
                                    epsilon = 1e-9)
        x *= 1/scale_factor
        
        if d['warnflag'] > 0:
            print d
                        
        # Check that the solution is valid and return
        # c = cantilever_divingboard(self.freq_min, self.freq_max, x)
        # c.check_if_valid_solution() # TODO
        
        (self.l, self.w, self.t, 
            self.l_pr_ratio, self.w_pr_ratio, self.w_gap_ratio, self.t_pr_ratio, 
            self.V_bias, self.N) = x
 
 
 
 
    # openopt version    
    def optimize_openopt(self, omega_min, lower_bounds, upper_bounds, max_power, fluid_type):

        def optimize_force_resolution(x):
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x)
            return c.force_resolution()*1e12
        
        # We always want power dissipation to be less than max_power
        def enforce_max_power(x):
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x)
            return c.power_dissipation() - max_power

        # We always want the natural frequency to be greater than the minimum allowed
        def enforce_minimum_natural_frequency(x):
            c = cantilever_divingboard(self.freq_min, self.freq_max, *x)

            if fluid_type == 'vacuum':
                return omega_min - c.omega_vacuum_hz()
            elif fluid_type == 'water':
                return omega_min - c.omega_damped_hz()
            else:
                return -c.omega_vacuum_hz()  
        
        # Our initial guess is based upon the existing cantilever dimensions
        initial_guess = (self.l, self.w, self.t, 
                            self.l_pr_ratio, self.w_pr_ratio, self.w_gap_ratio, self.t_pr_ratio, 
                            self.V_bias, self.N)
                        
        # Setup the solver
        p = NLP(optimize_force_resolution, initial_guess)
        p.lb = lower_bounds
        p.ub = upper_bounds
        p.scale=(1e-6, 1e-6, 1e-9, 1, 1, 1, 1, 1e2, 1e19)
                 
        # c(x) <= 0 constraints (use these if the solver supports NL constraints)
        # p.c = (enforce_minimum_natural_frequency,
        #        enforce_max_power)
        
        
        
        # Solving parameters
        p.ftol = 1e-12 # one of stop criteria, default 1e-6
        p.xtol = 1e-12 # one of stop criteria, default 1e-6
        p.gradtol = 1e-7
        p.contol = 1e-7
        
        # Set limits of how long to solve for
        p.maxIter = 1e4
        p.maxFunEvals = 1e8
        p.maxTime = 1e4
        
        # p.plot = 1
        # p.debug=1
        
        r = p.solve('algencan')
        
        (self.l, self.w, self.t, 
            self.l_pr_ratio, self.w_pr_ratio, self.w_gap_ratio, self.t_pr_ratio, 
            self.V_bias, self.N) = r.xf
        