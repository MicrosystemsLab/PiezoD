import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt
import abc

class Cantilever:
    __metaclass__ = abc.ABCMeta

    # Physical constants
    k_b = 1.38e-23                  # J/K
    k_b_eV = 8.617343e-5            # eV/K
    q = 1.60218e-19                 # Coulombs
    h_bar = 1.055e-34               # J-sec

    # Define the number of points to use in discretized calculations
    numFrequencyPoints = 1000       # For noise spectra plotting
    numXPoints = 800                # Points along cantilever length for deflection / temperature
    numZPoints = 200                # Points along cantilever depth for electrical / thermal calculations
    numRandomStressIterations = 10  # Number of Monte Carlo iterations for initial tip deflection calculations
    numOptimizationIterations = 20  # Max number of optimization attempts before giving up

    # Standard fluid properties
    k_water = 0.610                 # W/m-K
    rho_water = 996.6               # kg/m^3
    eta_water = 7.98e-4             # Pa-sec
    h_water = 49218                 # W/m^2-k

    k_air = 0.0262                  # W/m-K
    rho_air = 1.164                 # kg/m^3
    eta_air = 17e-6                 # Pa-sec
    h_air = 2098                    # W/m^2-K

    # For vacuum, use small but finite values for numerical stability
    k_vacuum = 1e-6                 # W/m-K
    rho_vacuum = 1e-6               # kg/m^3
    eta_vacuum = 1e-6               # Pa-sec
    h_vacuum = 1e-6                 # W/m^2-K

    # Thermal conductivities (W/m-K)
    # AlN: "Process-dependent thin-film thermal conductivities for thermal CMOS MEMS"
    # Al: "MEMS test structure for measuring thermal conductivity of thin films"
    k_si = 148
    k_al = 200
    k_sio2 = 1.4
    k_ti = 21.9
    k_aln = 60
    k_mo = 138

    # Coefficients of thermal expansion (1/K)
    alpha_al = 23.1e-6
    alpha_sio2 = 0.5e-6
    alpha_si = 2.6e-6
    alpha_ti = 8.6e-6
    alpha_aln = 4.5e-6
    alpha_mo = 4.8e-6

    # Silicon temperature coefficient of resistance (1/K)
    # Corresponds to a peak doping of about 1e20/cc.
    # Used to predict resistance change with self-heating in the advanced thermal models
    TCR = 1372e-6

    # Intrinsic stress (Pa)
    # For modeling tip deflection from film stress. Doping stress is treated separately elsewhere in the code
    # film_stress == 'nominal' uses the average
    # film_stress == 'random' uses a normal distribution
    sigma_si_range = 1e6*np.array([0, 0])
    sigma_sio2_range = 1e6*np.array([-200, -300])
    sigma_al_range = 1e6*np.array([160, 200])
    sigma_aln_range = 1e6*np.array([-300, 300])
    sigma_ti_range = 1e6*np.array([-25, 25])
    sigma_mo_range = 1e6*np.array([-25, 25])

    # Specific heat (J/kg-K) for calculating thermal time constants
    Cv_si = 700
    Cv_al = 910
    Cv_sio2 = 700

    # Poisson ratio (-)
    # Assume <100> direction for silicon, other materials are isotropic in-plane)
    nu_Si = 0.28
    nu_Al = 0.3
    nu_Ti = 0.3
    nu_AlN = 0.24
    nu_SiO2 = 0.17
    nu_Mo = 0.3

    # Elastic modulus (Pa)
    # Assume plane strain conditions (i.e. cantilever much wider than it is thick)
    E_si = 130e9/(1 - nu_Si**2)
    E_al = 70e9/(1 - nu_Al**2)
    E_ti = 90e9/(1 - nu_Ti**2)
    E_aln = 320e9/(1 - nu_AlN**2)
    E_sio2 = 75e9/(1 - nu_SiO2**2)
    E_mo = 329e9/(1 - nu_Mo**2)

    # Densities (kg/m^3)
    rho_si = 2330
    rho_al = 2700
    rho_ti = 4506
    rho_aln = 3260
    rho_sio2 = 2200
    rho_mo = 10280

    # Transverse piezoelectric coefficient of AlN (pm/V == pC/N)
    # d31 varies with thickness, so interpolate from literature values
    # Values are from papers written by the Piazza and Roukes groups
    d31_t = 1e-9*np.array([50, 100, 500, 3000])         # m
    d31_aln = -1e-12*np.array([1.9, 2.3, 2.5, 2.6])     # pm/V

    # Minimum and maximum physically realistic quality factors
    # Only viscous damping is modeled - TED, anchor loss, Akhiezer, etc are ignore
    maxQ = 5e3
    minQ = 1e-6

    # Define the possible optimization goals.
    goalForceResolution = 0
    goalDisplacementResolution = 1
    goalForceNoiseDensity = 2
    goalSurfaceStress = 3

    # Lookup table for calculating resonant frequency and quality factor in liquid
    # Source: "Oscillations of cylinders..." by Brumley, Wilcox and Sader (2010)
    # A = t/w ratio, beta = log(Re)
    A_lookup = np.array([0, 1/50, 1/20, 1/10, 1/5, 1/2, 1, 2, 5, 10, 20, 50, 1000])
    Beta_lookup = np.array([-3, -2.5, -2, -1.5, -1, -.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 100])

    # Hydrodynamic function
    # Includes fixed bottom row from published erratum
    gamma_lookup_real = np.array([
        [212.184, 213.310, 214.977, 217.701, 222.978, 237.780, 260.256, 207.210,  169.667,   154.616,   145.909,    139.855,    134.720],
        [91.6984, 92.2467, 93.0601, 94.3924, 96.9808, 104.295, 115.542, 88.9011,  70.8173,   63.7655,   59.7404,    56.9653,    54.6258],
        [41.6417, 41.9209, 42.3363, 43.0185, 44.3487, 48.1380, 54.0391, 39.8564,  30.6996,   27.2460,   25.3060,    23.9817,    22.8730],
        [20.1196, 20.2683, 20.4907, 20.8572, 21.5753, 23.6370, 26.8847, 18.8235,  13.9212,   12.1457,   11.1673,    10.5072,    9.95883],
        [10.4849, 10.5677, 10.6926, 10.8998, 11.3080, 12.4883, 14.3601, 9.43536,  6.64606,   5.68511,   5.16801,    4.82411,    4.54093],
        [5.96655, 6.01467, 6.08871, 6.21279, 6.45897, 7.17328, 8.30052, 5.04739,  3.35215,   2.80394,   2.51794,    2.33126,    2.17927],
        [3.73387, 3.76344, 3.81063, 3.89099, 4.05154, 4.51368, 5.22220, 2.89030,  1.78322,   1.45306,   1.28807,    1.18327,    1.09943],
        [2.56548, 2.58563, 2.61959, 2.67832, 2.79515, 3.11907, 3.58531, 1.77617,  0.994540,  0.783333,  0.684003,   0.623512,   0.576619],
        [1.91834, 1.93509, 1.96437, 2.01450, 2.11058, 2.35665, 2.68270, 1.17779,  0.580514,  0.435349,  0.372208,   0.336075,   0.309503],
        [1.54554, 1.56285, 1.59247, 1.64069, 1.72687, 1.92785, 2.17551, 0.848104, 0.357549,  0.249659,  0.206674,   0.184001,   0.168601],
        [1.32633, 1.34658, 1.37882, 1.42757, 1.50844, 1.68437, 1.88862, 0.663505, 0.235193,  0.148772,  0.117201,   0.102069,   0.0928779],
        [1.19577, 1.2202,  1.2555,  1.3051,  1.3833,  1.5459,  1.7259,  0.55939,  0.16703,   0.093131,  0.068128,   0.057273,   0.0515648],
        [1.11746, 1.1465,  1.1843,  1.2346,  1.3117,  1.4670,  1.6336,  0.50051,  0.12874,   0.062098,  0.040918,   0.032516,   0.0287745],
        [1,       1.04551, 1.08816, 1.14064, 1.21703, 1.36368, 1.51317, 0.423881, 0.0792129, 0.0222121, 0.00619303, 0.00113212, 0]
    ])

    gamma_lookup_imag = np.array([
        [1018.72,  1021.37,  1025.29,  1031.66,  1043.88,  1077.39,  1126.32,  1008.65,  915.159,  874.583,  850.149,  832.704,  817.599],
        [374.276,  375.392,  377.040,  379.721,  384.873,  399.079,  420.012,  370.057,  331.318,  314.778,  304.899,  297.884,  291.835],
        [140.659,  141.144,  141.862,  143.031,  145.284,  151.534,  160.848,  138.825,  122.228,  115.278,  111.167,  108.266,  105.776],
        [54.4049,  54.6253,  54.9508,  55.4818,  56.5079,  59.3754,  63.7087,  53.5749,  46.1812,  43.1534,  41.3825,  40.1420,  39.0836],
        [21.8269,  21.9314,  22.0855,  22.3371,  22.8247,  24.2002,  26.3169,  21.4324,  17.9905,  16.6153,  15.8210,  15.2692,  14.8012],
        [9.16870,  9.22024,  9.29587,  9.41936,  9.65973,  10.3480,  11.4345,  8.96804,  7.28929,  6.63516,  6.26219,  6.00523,  5.78862],
        [4.07467,  4.10043,  4.13779,  4.19895,  4.31957,  4.67605,  5.25977,  3.95920,  3.10274,  2.77671,  2.59298,  2.46733,  2.36186],
        [1.93366,  1.94552,  1.96256,  1.99130,  2.05107,  2.24127,  2.56535,  1.85252,  1.39790,  1.22868,  1.13429,  1.07013,  1.01639],
        [0.981710, 0.985312, 0.990956, 1.00255,  1.03157,  1.13634,  1.31768,  0.915797, 0.666095, 0.575374, 0.525354, 0.491568, 0.463359],
        [0.527773, 0.526433, 0.526077, 0.529479, 0.543868, 0.602276, 0.703142, 0.474037, 0.333253, 0.283225, 0.256021, 0.237799, 0.222666],
        [0.296143, 0.291987, 0.289093, 0.289338, 0.296683, 0.328687, 0.384789, 0.253907, 0.173548, 0.145302, 0.130165, 0.120135, 0.111868],
        [0.171115, 0.16564,  0.16234,  0.16171,  0.16525,  0.18260,  0.21384,  0.13910,  0.093151, 0.076988, 0.068405, 0.062790, 0.0582134],
        [0.100688, 0.095021, 0.092307, 0.091476, 0.093044, 0.10247,  0.11987,  0.077266, 0.051022, 0.041760, 0.036840, 0.033652, 0.0310905],
        [0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0]
    ])

    # Abstract methods
    @abc.abstractmethod
    def doping_profile(self):
        """Method documentation"""
        return

    @abc.abstractmethod
    def doping_optimization_scaling(self):
        return

    @abc.abstractmethod
    def doping_cantilever_from_state(self):
        return

    @abc.abstractmethod
    def doping_current_state(self):
        return

    @abc.abstractmethod
    def doping_initial_conditions_random(self):
        return

    @abc.abstractmethod
    def doping_optimization_bounds(self, parameter_constraints):
        return

    @abc.abstractmethod
    def Nz(self):
        return

    @abc.abstractmethod
    def alpha(self):
        return

    @abc.abstractmethod
    def sheet_resistance(self):
        return

    # TODO
    def __init__(self):
        # Initialize with reasonable defaults
        self.doping_type = 'phosphorus'
        self.freq_min = 1
        self.freq_max = 1e3
        self.l = 100e-6
        self.w = 10e-6
        self.t = 1e-6
        self.l_pr_ratio = 0.3
        self.v_bridge = 1
        self.air_gap_width = 2e-6            # The width between the two cantilever legs. Significant if l_pr is small

        self.fluid = 'air'
        self.h_method = 'fixed'
        self.metal_type = 'aluminum'
        self.film_stress = 'nominal'

        self.number_of_piezoresistors = 2
        self.amplifier = 'INA103'
        self.R_contact = 100
        self.tip_mass = 0

        self.cantilever_type = 'none'
        self.l_a = 0
        self.t_a = 0
        self.w_a = 0
        self.w_a_active = 0
        self.d31_manual = 0
        self.l_a_gap = 0
        self.t_oxide = 100e-9
        self.t_electrode_bottom = 50e-9
        self.t_electrode_top = 50e-9
        self.t_a_seed = 20e-9

        # Use simple thermal models by default
        self.T = 273.15 + 23
        self.T_ref = 273.15 + 23
        self.temperature_dependent_properties = 'no'
        self.thermal_modeling = 'none'
        self.R_base = 10e3                              # Assume R_base is 10e3 K/W by default

    def l_pr(self):
        return (self.l * self.l_pr_ratio)

    def w_pr(self):
        return (self.w/2)