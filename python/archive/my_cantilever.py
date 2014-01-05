from numpy import *
from scipy import *
from scipy.linalg import *
from pylab import *
import matplotlib.axes3d as p3

class cantilever:
    def __init__(self, l, w, t, E = 169.e9, rho_c = 2330., rho_l = 1000., eta = 1.e-3):

        # Cantilever geometry
        self.l = l
        self.w = w
        self.t = t
        
        self.w_leg
        self.w_gap
        
        self.t_piezoresistor
        self.l_piezoresistor
        
        self.l_piezoelectric
        self.w_piezoelectric

        # Mechanical properties
        self.E = E

        # Fluid properties
        self.rho_c = rho_c
        self.rho_l = rho_l
        self.eta = eta

        # Electrical properties
        self.V_bias = V_bias
        self.alpha = 1.e-6
        
        self.T = 300
        self.k_b = 1.38e-23

        e = 1.6e-19 # electron charge (C)
        mu = 480 # hole mobility (cm^2/V-s)
        
        self.N = 4e19 # Dopant concentration (N/cm^3)
        self.rho_nominal = 1/(self.N*e*mu)*1e-2 # Resistivity (ohm-m)
        self.rho_aluminum = 2.7e-6 # Resistivity of aluminum (ohm-m)

        self.freq_max = 100e3
        self.freq_min = 10

    # 
    def resistance(self):
        

    # Piezoresistance factor, which decreases with doping level
    def piezoresistance_factor(self):
        b = 1.53e22
        a = 0.2014
        P = log10((b/self.N) ** a)
        P_base = 72.e-11 # Pi at low concentration in 110 from Smith
        return P*P_base
        
    # Bending stiffness of the cantilever to a point load at the tip
    def stiffness(self):
        return self.E*self.w*self.t**3/(4*self.l**3)
    
    # Resonant frequency for undamped vibration (first mode)
    def omega_vacuum(self):
        m_effective = 0.243*self.rho_c*self.w*self.t*self.l
        omega_radians = sqrt( self.stiffness() / m_effective)
        return omega_radians

    def omega_vacuum_hz(self):
        return self.omega_vacuum() / (2*pi)

    # Calculate the damped natural frequency (assuming Q > 1) from Sader paper
    def omega_damped(self):
        def find_freq(omega):
            tau_r = real( self.hydrodynamic_function(omega) )
            residual = omega - self.omega_vacuum() * (1. + pi * self.rho_l * self.w /(4 * self.rho_c * self.t) * tau_r) ** -0.5
            return residual

        omega_radians = optimize.fsolve( find_freq, self.omega_vacuum() )
        return omega_radians

    def omega_damped_hz(self):
        return self.omega_damped() / (2*pi)
        
    # Reynold's number for flow around the cantilever for a given frequency (in rad/sec)
    def reynolds_number(self, omega):
        return (self.rho_l* omega * self.w**2) / self.eta

        
    # Normalized bending mode number. The first resonant frequency corresponds to the first
    # solution of the hyperbolic function 1+cos(x)*cosh(x) = 0
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