import math
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate, interpolate, optimize


@dataclass
class GapConfig:
    """Configuration for gap geometry (piezoresistor or beam).

    The gap defines a reduction in effective width near the cantilever root.
    This generalizes the "air gap" concept for U-shaped piezoresistors and
    can also be used for cantilever beam geometry (e.g., tuning-fork designs).

    The effective width at position x is:
        - For 0 <= x <= gap_extent: w_nominal - gap_width
        - For x > gap_extent: w_nominal (full width)

    Attributes:
        gap_width: Absolute gap width reduction (m). If zero, gap_fraction is used.
        gap_fraction: Gap as fraction of nominal width (-). Only used if gap_width is 0.
        gap_extent: Length from root where gap applies (m). If None, uses full
            reference length (e.g., piezoresistor length or cantilever length).
    """

    gap_width: float = 0.0  # Absolute gap width (m)
    gap_fraction: float = 0.0  # Fractional gap (0-1)
    gap_extent: float | None = None  # Length from root where gap applies (m)

    def get_gap_width(self, nominal_width: float) -> float:
        """Calculate the actual gap width.

        Args:
            nominal_width: The nominal (full) width (m)

        Returns:
            The gap width reduction (m)
        """
        if self.gap_width > 0:
            return self.gap_width
        return self.gap_fraction * nominal_width

    def effective_width(self, x: float, reference_length: float, nominal_width: float) -> float:
        """Calculate effective width at position x.

        The gap reduction applies from the root (x=0) up to gap_extent.
        Beyond gap_extent, the full nominal width is used.

        Args:
            x: Position along cantilever from root (m)
            reference_length: Reference length for gap extent default (m)
            nominal_width: Nominal (unreduced) width (m)

        Returns:
            Effective width at position x (m)
        """
        extent = self.gap_extent if self.gap_extent is not None else reference_length
        gap_w = self.get_gap_width(nominal_width)

        if x <= extent:
            return nominal_width - gap_w
        return nominal_width

    def effective_width_array(self, x: np.ndarray, reference_length: float, nominal_width: float) -> np.ndarray:
        """Calculate effective width at each position in array x.

        Args:
            x: Array of positions along cantilever from root (m)
            reference_length: Reference length for gap extent default (m)
            nominal_width: Nominal (unreduced) width (m)

        Returns:
            Array of effective widths at each position (m)
        """
        extent = self.gap_extent if self.gap_extent is not None else reference_length
        gap_w = self.get_gap_width(nominal_width)

        w_eff = np.where(x <= extent, nominal_width - gap_w, nominal_width)
        return w_eff


class Cantilever:
    # Physical constants
    k_b = 1.38e-23  # J/K
    k_b_eV = 8.617343e-5  # eV/K
    q = 1.60218e-19  # Coulombs
    h_bar = 1.055e-34  # J-sec

    dopantOptions = {"boron": 1, "phosphorus": 2, "arsenic": 3}

    # Define the number of points to use in discretized calculations
    numFrequencyPoints = 1000  # For noise spectra plotting
    numXPoints = 800  # Points along cantilever length for deflection / temperature
    numZPoints = 200  # Points along cantilever depth for electrical / thermal calculations
    numRandomStressIterations = 10  # Number of Monte Carlo iterations for initial tip deflection calculations
    numOptimizationIterations = 20  # Max number of optimization attempts before giving up

    # Standard fluid properties
    k_water = 0.610  # W/m-K
    rho_water = 996.6  # kg/m^3
    eta_water = 7.98e-4  # Pa-sec
    h_water = 49218  # W/m^2-k

    k_air = 0.0262  # W/m-K
    rho_air = 1.164  # kg/m^3
    eta_air = 17e-6  # Pa-sec
    h_air = 2098  # W/m^2-K

    # For vacuum, use small but finite values for numerical stability
    k_vacuum = 1e-6  # W/m-K
    rho_vacuum = 1e-6  # kg/m^3
    eta_vacuum = 1e-6  # Pa-sec
    h_vacuum = 1e-6  # W/m^2-K

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

    # Default Hooge noise parameter (unitless)
    # Used by epitaxy and diffusion subclasses
    default_alpha = 1e-5

    # Intrinsic stress (Pa)
    # For modeling tip deflection from film stress. Doping stress is treated separately elsewhere in the code
    # film_stress == 'nominal' uses the average
    # film_stress == 'random' uses a normal distribution
    sigma_si_range = 1e6 * np.array([0, 0])
    sigma_sio2_range = 1e6 * np.array([-200, -300])
    sigma_al_range = 1e6 * np.array([160, 200])
    sigma_aln_range = 1e6 * np.array([-300, 300])
    sigma_ti_range = 1e6 * np.array([-25, 25])
    sigma_mo_range = 1e6 * np.array([-25, 25])

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
    E_si = 130e9 / (1 - nu_Si**2)
    E_al = 70e9 / (1 - nu_Al**2)
    E_ti = 90e9 / (1 - nu_Ti**2)
    E_aln = 320e9 / (1 - nu_AlN**2)
    E_sio2 = 75e9 / (1 - nu_SiO2**2)
    E_mo = 329e9 / (1 - nu_Mo**2)

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
    d31_t = 1e-9 * np.array([50, 100, 500, 3000])  # m
    d31_aln = -1e-12 * np.array([1.9, 2.3, 2.5, 2.6])  # pm/V

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
    A_lookup = np.array([0, 1 / 50, 1 / 20, 1 / 10, 1 / 5, 1 / 2, 1, 2, 5, 10, 20, 50, 1000])
    Beta_lookup = np.array([-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 100])

    # Hydrodynamic function
    # Includes fixed bottom row from published erratum
    gamma_lookup_real = np.array(
        [
            [
                212.184,
                213.310,
                214.977,
                217.701,
                222.978,
                237.780,
                260.256,
                207.210,
                169.667,
                154.616,
                145.909,
                139.855,
                134.720,
            ],
            [
                91.6984,
                92.2467,
                93.0601,
                94.3924,
                96.9808,
                104.295,
                115.542,
                88.9011,
                70.8173,
                63.7655,
                59.7404,
                56.9653,
                54.6258,
            ],
            [
                41.6417,
                41.9209,
                42.3363,
                43.0185,
                44.3487,
                48.1380,
                54.0391,
                39.8564,
                30.6996,
                27.2460,
                25.3060,
                23.9817,
                22.8730,
            ],
            [
                20.1196,
                20.2683,
                20.4907,
                20.8572,
                21.5753,
                23.6370,
                26.8847,
                18.8235,
                13.9212,
                12.1457,
                11.1673,
                10.5072,
                9.95883,
            ],
            [
                10.4849,
                10.5677,
                10.6926,
                10.8998,
                11.3080,
                12.4883,
                14.3601,
                9.43536,
                6.64606,
                5.68511,
                5.16801,
                4.82411,
                4.54093,
            ],
            [
                5.96655,
                6.01467,
                6.08871,
                6.21279,
                6.45897,
                7.17328,
                8.30052,
                5.04739,
                3.35215,
                2.80394,
                2.51794,
                2.33126,
                2.17927,
            ],
            [
                3.73387,
                3.76344,
                3.81063,
                3.89099,
                4.05154,
                4.51368,
                5.22220,
                2.89030,
                1.78322,
                1.45306,
                1.28807,
                1.18327,
                1.09943,
            ],
            [
                2.56548,
                2.58563,
                2.61959,
                2.67832,
                2.79515,
                3.11907,
                3.58531,
                1.77617,
                0.994540,
                0.783333,
                0.684003,
                0.623512,
                0.576619,
            ],
            [
                1.91834,
                1.93509,
                1.96437,
                2.01450,
                2.11058,
                2.35665,
                2.68270,
                1.17779,
                0.580514,
                0.435349,
                0.372208,
                0.336075,
                0.309503,
            ],
            [
                1.54554,
                1.56285,
                1.59247,
                1.64069,
                1.72687,
                1.92785,
                2.17551,
                0.848104,
                0.357549,
                0.249659,
                0.206674,
                0.184001,
                0.168601,
            ],
            [
                1.32633,
                1.34658,
                1.37882,
                1.42757,
                1.50844,
                1.68437,
                1.88862,
                0.663505,
                0.235193,
                0.148772,
                0.117201,
                0.102069,
                0.0928779,
            ],
            [
                1.19577,
                1.2202,
                1.2555,
                1.3051,
                1.3833,
                1.5459,
                1.7259,
                0.55939,
                0.16703,
                0.093131,
                0.068128,
                0.057273,
                0.0515648,
            ],
            [
                1.11746,
                1.1465,
                1.1843,
                1.2346,
                1.3117,
                1.4670,
                1.6336,
                0.50051,
                0.12874,
                0.062098,
                0.040918,
                0.032516,
                0.0287745,
            ],
            [
                1,
                1.04551,
                1.08816,
                1.14064,
                1.21703,
                1.36368,
                1.51317,
                0.423881,
                0.0792129,
                0.0222121,
                0.00619303,
                0.00113212,
                0,
            ],
        ]
    )

    gamma_lookup_imag = np.array(
        [
            [
                1018.72,
                1021.37,
                1025.29,
                1031.66,
                1043.88,
                1077.39,
                1126.32,
                1008.65,
                915.159,
                874.583,
                850.149,
                832.704,
                817.599,
            ],
            [
                374.276,
                375.392,
                377.040,
                379.721,
                384.873,
                399.079,
                420.012,
                370.057,
                331.318,
                314.778,
                304.899,
                297.884,
                291.835,
            ],
            [
                140.659,
                141.144,
                141.862,
                143.031,
                145.284,
                151.534,
                160.848,
                138.825,
                122.228,
                115.278,
                111.167,
                108.266,
                105.776,
            ],
            [
                54.4049,
                54.6253,
                54.9508,
                55.4818,
                56.5079,
                59.3754,
                63.7087,
                53.5749,
                46.1812,
                43.1534,
                41.3825,
                40.1420,
                39.0836,
            ],
            [
                21.8269,
                21.9314,
                22.0855,
                22.3371,
                22.8247,
                24.2002,
                26.3169,
                21.4324,
                17.9905,
                16.6153,
                15.8210,
                15.2692,
                14.8012,
            ],
            [
                9.16870,
                9.22024,
                9.29587,
                9.41936,
                9.65973,
                10.3480,
                11.4345,
                8.96804,
                7.28929,
                6.63516,
                6.26219,
                6.00523,
                5.78862,
            ],
            [
                4.07467,
                4.10043,
                4.13779,
                4.19895,
                4.31957,
                4.67605,
                5.25977,
                3.95920,
                3.10274,
                2.77671,
                2.59298,
                2.46733,
                2.36186,
            ],
            [
                1.93366,
                1.94552,
                1.96256,
                1.99130,
                2.05107,
                2.24127,
                2.56535,
                1.85252,
                1.39790,
                1.22868,
                1.13429,
                1.07013,
                1.01639,
            ],
            [
                0.981710,
                0.985312,
                0.990956,
                1.00255,
                1.03157,
                1.13634,
                1.31768,
                0.915797,
                0.666095,
                0.575374,
                0.525354,
                0.491568,
                0.463359,
            ],
            [
                0.527773,
                0.526433,
                0.526077,
                0.529479,
                0.543868,
                0.602276,
                0.703142,
                0.474037,
                0.333253,
                0.283225,
                0.256021,
                0.237799,
                0.222666,
            ],
            [
                0.296143,
                0.291987,
                0.289093,
                0.289338,
                0.296683,
                0.328687,
                0.384789,
                0.253907,
                0.173548,
                0.145302,
                0.130165,
                0.120135,
                0.111868,
            ],
            [
                0.171115,
                0.16564,
                0.16234,
                0.16171,
                0.16525,
                0.18260,
                0.21384,
                0.13910,
                0.093151,
                0.076988,
                0.068405,
                0.062790,
                0.0582134,
            ],
            [
                0.100688,
                0.095021,
                0.092307,
                0.091476,
                0.093044,
                0.10247,
                0.11987,
                0.077266,
                0.051022,
                0.041760,
                0.036840,
                0.033652,
                0.0310905,
            ],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ]
    )

    # Abstract methods
    def doping_profile(self):
        raise NotImplementedError("Implement doping_profile() in your Cantilever subclass")

    def doping_optimization_scaling(self):
        raise NotImplementedError("Implement doping_optimization_scaling() in your Cantilever subclass")

    def doping_cantilever_from_state(self, x0):
        raise NotImplementedError("Implement doping_cantilever_from_state() in your Cantilever subclass")

    def doping_current_state(self):
        raise NotImplementedError("Implement doping_current_state() in your Cantilever subclass")

    def doping_initial_conditions_random(self):
        raise NotImplementedError("Implement doping_initial_conditions_random() in your Cantilever subclass")

    def doping_optimization_bounds(self, parameter_constraints):
        raise NotImplementedError("Implement doping_optimization_bounds() in your Cantilever subclass")

    def Nz(self):
        raise NotImplementedError("Implement Nz() in your Cantilever subclass")

    def alpha(self):
        raise NotImplementedError("Implement alpha() in your Cantilever subclass")

    def sheet_resistance(self):
        raise NotImplementedError("Implement sheet_resistance() in your Cantilever subclass")

    # Method for iteratively finding the equivalent thickness
    def findEIResidual(self, t_equivalent_guess):
        EI_calculated = self.modulus() * self.w_a * t_equivalent_guess**3 / 12
        EI_effective = 1 / self.calculateActuatorNormalizedCurvature()
        residual = (EI_effective / EI_calculated - 1) ** 2
        return residual

    # Function for iteratively finding the resonant frequency
    def findEnergyResidual(self, omega_guess):
        U_elastic, U_kinetic = self.calculateEnergies(omega_guess)
        return (U_elastic / U_kinetic - 1) ** 2

    # Function for computing the elastic and kinetic energy of the beam
    def calculateEnergies(self, omega):
        # Discretize the length of the cantilever
        totalLength = self.l + self.l_a
        dx = totalLength / (self.numXPoints - 1)
        x = np.arange(0, totalLength + 1, dx)
        base_indices = np.nonzero(x <= self.l_a)
        tip_indices = np.nonzero(x > self.l_a)

        deflection = np.zeros((self.numXPoints, 1))
        Udx_elastic = np.zeros((self.numXPoints, 1))
        Udx_kinetic = np.zeros((self.numXPoints, 1))

        # Calculate effective beam width profile (for beam gap handling)
        # x coordinates in tip region are relative to start of l_a
        x_tip_relative = x[tip_indices] - self.l_a if self.l_a > 0 else x[tip_indices]
        w_eff_tip = self.effective_beam_width_array(x_tip_relative.flatten())

        # Define the multilayer mechanics
        EI_base = 1 / self.calculateActuatorNormalizedCurvature()
        # Use average effective width for the simple EI_tip estimate
        # (the detailed calculation below uses spatially varying EI)
        w_avg_tip = np.mean(w_eff_tip) if len(w_eff_tip) > 0 else self.w
        EI_tip_avg = self.modulus() * w_avg_tip * self.t**3 / 12

        # Empirical correction factors that give better agreement with FEA results
        # Account for cases where t_a >> t, w_a >> w, l_a >> l
        EI_base *= (self.t / self.t_a) ** 0.25
        EI_base *= (self.w / self.w_a) ** 0.1
        EI_base *= (self.l / self.l_a) ** 0.1

        # Generate an approximate cantilever deflection profile assuming a
        # point load force at the tip of the beam. Stitch together the two
        # sections (the moment is constant despite the EI discontinuity)
        tip_deflection = 1e-6  # Apply a test force
        F = tip_deflection * 3 * EI_tip_avg / self.l**3
        moment = F * (totalLength - x)
        deflection[base_indices] = -F * x[base_indices] ** 2 * (3 * totalLength - x[base_indices]) / (6 * EI_base)

        # x-coordinate from the end of l_a
        x_relative = x[tip_indices] - x[max(base_indices)]

        # Continue with the slope that is at the end of the base section
        if max(base_indices) > 1:
            tip_slope = (deflection[max(base_indices)] - deflection[max(base_indices) - 1]) / dx
        else:
            tip_slope = 0

        deflection[tip_indices] = (
            deflection[max(base_indices)]
            - F * x_relative**2 * (3 * self.l - x_relative) / (6 * EI_tip_avg)
            + tip_slope * x_relative
        )

        E_metal, rho_metal, k_metal, alpha_metal = self.lookup_metal_properties()

        # Spatially varying EI and mass for tip region (beam gap support)
        EI_tip = self.modulus() * w_eff_tip * self.t**3 / 12
        dm_tip = w_eff_tip * self.t * self.rho_si

        if self.cantilever_type in ("step", "thermal"):
            dm_base = self.w_a * (self.t * self.rho_si + self.t_oxide * self.rho_sio2 + self.t_a * rho_metal)
        else:
            dm_base = self.w_a * (
                self.t * self.rho_si
                + self.t_oxide * self.rho_sio2
                + self.rho_aln * (self.t_a + self.t_a_seed)
                + rho_metal * (self.t_electrode_bottom + self.t_electrode_top)
            )

        # Piecewise kinetic and elastic energies
        Udx_elastic[base_indices] = 0.5 * moment[base_indices] ** 2 * dx / EI_base
        Udx_kinetic[base_indices] = 0.5 * (omega * deflection[base_indices]) ** 2 * dx * dm_base

        # Tip region uses spatially varying EI and dm (for beam gap)
        tip_idx_flat = tip_indices[0]
        Udx_elastic[tip_idx_flat] = (0.5 * moment[tip_idx_flat] ** 2 * dx / EI_tip).reshape(-1, 1)
        Udx_kinetic[tip_idx_flat] = (0.5 * (omega * deflection[tip_idx_flat].flatten()) ** 2 * dx * dm_tip).reshape(
            -1, 1
        )

        U_elastic = np.trapezoid(x, Udx_elastic)
        U_kinetic = np.trapezoid(x, Udx_kinetic)
        return U_elastic, U_kinetic

    # Methods
    def __init__(self):
        # Initialize with reasonable defaults
        self.doping_type = "phosphorus"
        self.freq_min = 1
        self.freq_max = 1e3
        self.l = 100e-6
        self.w = 10e-6
        self.t = 1e-6
        self.l_pr_ratio = 0.3
        self.v_bridge = 1

        # Gap configuration for U-shaped piezoresistor
        # Default: 2 um gap width applied over entire piezoresistor length (root region)
        self._gap_config = GapConfig(gap_width=2e-6)

        # Beam gap configuration for cantilever geometry
        # Used in Rayleigh-Ritz frequency calculation when the beam has a slot/gap
        # at the root (e.g., tuning fork designs). Default: no beam gap.
        self._beam_gap_config: GapConfig | None = None

        self.fluid = "air"
        self.rho_arb = 1.0  # Arbitrary fluid density (kg/m^3)
        self.eta_arb = 1e-3  # Arbitrary fluid viscosity (Pa-s)
        self.k_arb = 0.1  # Arbitrary fluid thermal conductivity (W/m-K)
        self.h_arb = 100  # Arbitrary fluid heat transfer coefficient (W/m^2-K)
        self.h_method = "fixed"
        self.metal_type = "aluminum"
        self.film_stress = "nominal"

        self.number_of_piezoresistors = 2
        self.amplifier = "INA103"
        self.R_contact = 100
        self.tip_mass = 0
        self.rms_actuator_displacement_noise = 1e-12  # Displacement noise from mounting (m)

        self.cantilever_type = "none"
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
        self.temperature_dependent_properties = "no"
        self.thermal_modeling = "none"
        self.R_base = 10e3
        self.R_heater = 1e3

    def l_pr(self):
        return self.l * self.l_pr_ratio

    def w_pr(self):
        return self.w / 2

    @property
    def gap_config(self) -> GapConfig:
        """Get the gap configuration for the piezoresistor."""
        return self._gap_config

    @gap_config.setter
    def gap_config(self, config: GapConfig) -> None:
        """Set the gap configuration for the piezoresistor."""
        self._gap_config = config

    @property
    def air_gap_width(self) -> float:
        """Get the air gap width (backward compatibility).

        This property provides backward compatibility with the old air_gap_width
        attribute. For new code, use gap_config instead.

        Returns:
            The gap width in meters.
        """
        return self._gap_config.get_gap_width(self.w_pr())

    @air_gap_width.setter
    def air_gap_width(self, width: float) -> None:
        """Set the air gap width (backward compatibility).

        This property provides backward compatibility with the old air_gap_width
        attribute. For new code, use gap_config instead.

        Args:
            width: The gap width in meters.
        """
        self._gap_config = GapConfig(gap_width=width)

    def effective_pr_width_at(self, x: float) -> float:
        """Calculate effective piezoresistor width at position x.

        The gap reduces the conducting width in the root region of the
        piezoresistor. This method returns the effective width at any
        position along the cantilever.

        Args:
            x: Position along cantilever from root (m)

        Returns:
            Effective piezoresistor width at position x (m)
        """
        return self._gap_config.effective_width(x, self.l_pr(), self.w_pr())

    def gap_extent(self) -> float:
        """Get the extent of the piezoresistor gap region from the root.

        Returns:
            Length from root where gap applies (m). Defaults to l_pr if not set.
        """
        if self._gap_config.gap_extent is not None:
            return self._gap_config.gap_extent
        return self.l_pr()

    @property
    def beam_gap_config(self) -> GapConfig | None:
        """Get the beam gap configuration for Rayleigh-Ritz calculations.

        The beam gap reduces the effective cantilever width in the root region,
        affecting stiffness and resonant frequency calculations. This is used
        for tuning-fork or slotted cantilever designs.

        Returns:
            GapConfig for beam geometry, or None if no beam gap is configured.
        """
        return self._beam_gap_config

    @beam_gap_config.setter
    def beam_gap_config(self, config: GapConfig | None) -> None:
        """Set the beam gap configuration for Rayleigh-Ritz calculations.

        Args:
            config: GapConfig for beam geometry, or None to disable beam gap.
        """
        self._beam_gap_config = config

    def effective_beam_width_at(self, x: float) -> float:
        """Calculate effective cantilever beam width at position x.

        If a beam gap is configured, the width is reduced in the root region.
        This affects stiffness, mass, and resonant frequency calculations.

        Args:
            x: Position along cantilever from root (m)

        Returns:
            Effective beam width at position x (m)
        """
        if self._beam_gap_config is None:
            return self.w
        return self._beam_gap_config.effective_width(x, self.l, self.w)

    def effective_beam_width_array(self, x: np.ndarray) -> np.ndarray:
        """Calculate effective cantilever beam width at each position in array.

        If a beam gap is configured, the width is reduced in the root region.
        This affects stiffness, mass, and resonant frequency calculations.

        Args:
            x: Array of positions along cantilever from root (m)

        Returns:
            Array of effective beam widths at each position (m)
        """
        if self._beam_gap_config is None:
            return np.full_like(x, self.w)
        return self._beam_gap_config.effective_width_array(x, self.l, self.w)

    def beam_gap_extent(self) -> float:
        """Get the extent of the beam gap region from the root.

        Returns:
            Length from root where beam gap applies (m). Returns 0 if no gap.
        """
        if self._beam_gap_config is None:
            return 0.0
        if self._beam_gap_config.gap_extent is not None:
            return self._beam_gap_config.gap_extent
        return self.l

    # Determine the ion implantation table index from the dopant type
    # TODO Catch exceptions
    def dopantNumber(self):
        return Cantilever.dopantOptions[self.doping_type]

    # Check if the cantilever is self-consistent
    def check_valid_cantilever(self):
        if self.cantilever_type == "none" and self.l_a > 0:
            raise RuntimeError("Cantilever is not valid: type = 'none' and l_a > 0")
        return

    def print_performance(self) -> None:
        """Print cantilever performance summary."""
        self.check_valid_cantilever()

        omega_damped_hz, Q = self.omega_damped_hz_and_Q()
        x, active_doping, total_doping = self.doping_profile()
        TMax_approx, TTip_approx = self.approxTempRise()
        TMax, TTip = self.calculateMaxAndTipTemp()
        thermoLimit = self.thermo_integrated() / self.force_sensitivity()

        print("=======================")
        print(f"Freq range: {self.freq_min} to {self.freq_max}")
        print(f"Operating fluid: {self.fluid}")
        print(f"Cantilever L/W/T: {self.l * 1e6} {self.w * 1e6} {self.t * 1e6}")
        print(f"PR L/W: {self.l_pr() * 1e6} {self.w_pr() * 1e6}")
        print(f"PR Length Ratio: {self.l_pr_ratio:g}")

        print(f"Force resolution (N): {self.force_resolution():g}")
        print(f"Force noise at 1 kHz (fN): {self.force_noise_density(1e3):g}")
        print(f"Displacement resolution (m): {self.displacement_resolution():g}")
        print(f"Force sensitivity (V/N): {self.force_sensitivity():g}")
        print(f"Displacement sensitivity (V/m): {self.displacement_sensitivity():g}")
        print(f"Beta: {self.beta():g}")
        print(f"Thermomechanical force noise limit: {thermoLimit:g}")

        print(f"Stiffness (N/m): {self.stiffness():g}")
        print(f"Vacuum freq: {self.omega_vacuum_hz()}")
        print(f"Damped freq: {omega_damped_hz}")
        print(f"Quality factor: {Q}")

        print(f"Wheatstone bridge bias voltage: {self.v_bridge}")
        print(f"Resistance: {self.resistance()}")
        print(f"Sheet Resistance: {self.sheet_resistance()}")
        print(f"Power dissipation (mW): {self.power_dissipation() * 1e3:g}")
        print(f"Approx. Temp Rises (C) - Tip: {TTip_approx}  Max: {TMax_approx}")
        print(f"F-D Temp Rises (C)     - Tip: {TTip}  Max: {TMax}")

        print(f"Integrated noise (V): {self.integrated_noise():g}")
        print(f"Integrated johnson noise (V): {self.johnson_integrated():g}")
        print(f"Integrated 1/f noise (V): {self.hooge_integrated():g}")
        print(f"Amplifier noise (V): {self.amplifier_integrated():g}")
        print(f"Thermomechanical noise (V): {self.thermo_integrated():g}")

        print(f"Johnson/Hooge: {self.johnson_integrated() / self.hooge_integrated():g}")
        print(f"Knee frequency (Hz): {self.knee_frequency():g}")
        print(f"Number of Carriers: {self.number_of_carriers():g}")
        print(f"Nz: {self.Nz():g}")

        print(f"Number of silicon resistors: {self.number_of_piezoresistors}")
        print(f"Si Thermal Conductivity (W/m-K): {self.k_base()}")
        print(f"E (GPa): {self.modulus() * 1e-9}")
        print(f"Alpha: {self.alpha:g}")

        if self.cantilever_type == "step":
            print("=======================")
            print(f"Step at base (um): {1e6 * self.t_a} thick x {1e6 * self.l_a} long")
        elif self.cantilever_type == "thermal":
            print("=======================")
            tau, freq = self.heaterTimeConstant()
            print(f"Actuator l/W/T: {1e6 * self.l_a} {1e6 * self.w_a} {1e6 * self.t_a}")
            print(f"Neutral axis (um): {1e6 * self.actuatorNeutralAxis()}")
            print(f"Actuator Voltage: {self.v_actuator}")
            print(f"Heater resistance (kOhm): {1e-3 * self.R_heater}")
            print(f"Actuator Power (mW): {1e3 * self.heaterPower()}")
            print(f"Tip Deflection (nm): {1e9 * self.tipDeflection()}")
            print(f"Time Constant (microseconds): {tau * 1e6}")
            print(f"-3dB frequency (kHz): {freq * 1e-3}")
        elif self.cantilever_type == "piezoelectric":
            print(f"Actuator l/W/T: {1e6 * self.l_a} {1e6 * self.w_a} {1e6 * self.t_a}")
            print(f"Neutral axis (um): {1e6 * self.actuatorNeutralAxis()}")
            print(f"Actuator Voltage: {self.v_actuator}")
            print(f"Tip Deflection (nm): {1e9 * self.tipDeflection()}")
        print("=======================")

    def print_performance_for_excel(self):
        # TODO
        return

    # Calculate total resistance of piezoresistor (ohms)
    def resistance(self):
        return self.number_of_squares() * self.sheet_resistance() + 2 * self.R_contact

    # Calculate the number of resistor squares (-)
    # 1) Longitudinal region (2*l_pr/w_pr)
    # 2) Transverse at end (air_gap_width/2*w_pr)
    # 3) the connecting corners
    def number_of_squares(self):
        return 2 * self.l_pr() / self.w_pr() + self.air_gap_width / (2 * self.w_pr()) + 2

    # Calculate R-R0/R0 for the cantilever including the TCR (-)
    # Used for calculating h_eff from experimental results
    # TODO Refactor
    def dR_with_temp_rise(self):
        return self.TCR * (self.averagePRTemp() - self.T)

    def approx_dR_with_temp_rise(self):
        return self.TCR * (self.approxPRTemp() - self.T)

    # Calculate conductivity for a given dopant concentration (C/V-sec-cm)
    def conductivity(self, dopant_concentration):
        TPR = self.piezoresistor_temp()
        mu, sigma = self.mobility(dopant_concentration, TPR)
        return sigma

    # Calculate the temp dependent carrier density, mobility, conductivity (C/V-sec-cm)
    # Current: "Electron and Hole Mobility ...", Reggiani et al. (2002)
    # Previously: "Modeling of Carrier Mobility ...", Masetti et al. (1983)
    def mobility(self, dopantConc, T):
        Tnorm = T / 300
        Eg = 1.170 - (4.730e-4 * T**2) / (T + 636)
        ni = (2.4e31 * T**3 * np.exp(-Eg / (self.k_b_eV * T))) ** 0.5

        if self.doping_type == "boron":
            mumax = 470.5
            c = 0.0
            gamma = 2.16
            mu0d = 90.0 * Tnorm**-1.3
            mu0a = 44.0 * Tnorm**-0.7
            mu1d = 28.2 * Tnorm**-2.0
            mu1a = 28.2 * Tnorm**-0.8
            Cr1 = 1.3e18 * Tnorm**2.2
            Cr2 = 2.45e17 * Tnorm**3.1
            Cs1 = 1.1e18 * Tnorm**6.2
            Cs2 = 6.1e20
            alpha1 = 0.77
            alpha2 = 0.719
            ni = (2.4e31 * T**3 * np.exp(-1 * Eg / 8.617e-5 / T)) ** 0.5
            p = (dopantConc / 2) + ((dopantConc / 2) ** 2 + ni**2) ** 0.5
            n = ni * ni / p
            ND = n
            NA = p
        elif self.doping_type == "phosphorus":
            mumax = 1441
            c = 0.07
            gamma = 2.45
            mu0d = 62.2 * Tnorm**-0.7
            mu0a = 132 * Tnorm**-1.3
            mu1d = 48.6 * Tnorm**-0.7
            mu1a = 73.5 * Tnorm**-1.25
            Cr1 = 8.5e16 * Tnorm**3.65
            Cr2 = 1.22e17 * Tnorm**2.65
            Cs1 = 4e20
            Cs2 = 7e20
            alpha1 = 0.68
            alpha2 = 0.72
            n = (dopantConc / 2) + ((dopantConc / 2) ** 2 + ni**2) ** 0.5
            p = ni**2 / n
            ND = n
            NA = p
        elif self.doping_type == "arsenic":
            mumax = 1441
            c = 0.07
            gamma = 2.45
            mu0d = 55 * Tnorm**-0.6
            mu0a = 132 * Tnorm**-1.3
            mu1d = 42.4 * Tnorm**-0.5
            mu1a = 73.5 * Tnorm**-1.25
            Cr1 = 8.9e16 * Tnorm**3.65
            Cr2 = 1.22e17 * Tnorm**2.65
            Cs1 = 2.9e20
            Cs2 = 7e20
            alpha1 = 0.68
            alpha2 = 0.72
            n = (dopantConc / 2) + ((dopantConc / 2) ** 2 + ni**2) ** 0.5
            p = ni**2 / n
            ND = n
            NA = p
        else:
            raise RuntimeError("Unknown dopant type!")

        mu0 = (mu0d * ND + mu0a * NA) / (ND + NA)
        mu1 = (mu1d * ND + mu1a * NA) / (ND + NA)
        muL = mumax * Tnorm ** (-gamma + c * Tnorm)
        second = (muL - mu0) / (1 + (ND / Cr1) ** alpha1 + (NA / Cr2) ** alpha2)
        third = mu1 / (1 + (ND / Cs1 + NA / Cs2) ** -2)

        mu = mu0 + second - third
        sigma = self.q * (mu * n + mu * p)
        return (mu, sigma)

    ## Calculate Rsheet(x) assuming temperature dependent cantilever properties (ohm/sq)
    def RSheetProfile(self, x, T_x):
        z, active_doping, total_doping = self.doping_profile()  # Units: z -> m, doping -> N/cm^3
        n_z_x = np.transpose(total_doping) * np.ones((1, self.numXPoints))

        # Generate a numZPoints x numXPoints matrix
        if len(T_x) == 1:
            T_z_x = np.dot(np.ones((self.numZPoints, 1)), np.transpose(T_x) * np.ones((1, self.numXPoints)))
        else:
            T_z_x = np.dot(np.ones((self.numZPoints, 1)), np.transpose(T_x))

        # TODO
        mu_z_x, sigma_z_x = self.mobility(n_z_x, T_z_x + self.T)
        Rsheet_x = 1 / np.trapezoid(z * 1e2, sigma_z_x)  # Convert z from m to cm
        return Rsheet_x

    # The number of current carriers in the piezoresistor (-)
    def number_of_carriers(self):
        resistor_area = self.w_pr() * (2 * self.l_pr() + self.w + self.air_gap_width)
        return self.Nz() * resistor_area

    # 1/f voltage power spectral density for the entire Wheatstone bridge (V^2/Hz)
    def hooge_PSD(self, freq):
        return self.alpha() * self.v_bridge**2 * self.number_of_piezoresistors / (4 * self.number_of_carriers() * freq)

    # Integrated 1/f noise density for the entire Wheatstone bridge (V)
    def hooge_integrated(self):
        return math.sqrt(
            self.alpha()
            * self.v_bridge**2
            * self.number_of_piezoresistors
            / (4 * self.number_of_carriers())
            * math.log(self.freq_max / self.freq_min)
        )

    # Johnson noise PSD from the entire Wheatstone bridge (V^2/Hz)
    def johnson_PSD(self, freq):
        resistance = self.resistance()  # resistance() includes contacts
        TPR = self.piezoresistor_temp()
        resistance *= 1 + TPR * self.TCR

        # If using 2 piezoresistors, assume ideal 1 kOhm external resistors
        R_external = resistance if self.number_of_piezoresistors == 4 else 1e3

        return 4 * self.k_b * TPR * (resistance / 2 + R_external / 2) * np.ones((1, freq.size))

    # Integrated Johnson noise
    # Unit: V
    def johnson_integrated(self):
        resistance = self.resistance()
        TPR = self.approxPRTemp()
        resistance *= 1 + TPR * self.TCR

        R_external = resistance if self.number_of_piezoresistors == 4 else 700
        return math.sqrt(4 * self.k_b * TPR * (resistance / 2 + R_external / 2) * (self.freq_max - self.freq_min))

    def piezoresistor_temp(self):
        if self.thermal_modeling == "approx":
            TPR = self.approxPRTemp()
        elif self.thermal_modeling == "exact":
            TPR = self.averagePRTemp()
        else:
            TPR = self.T  # the ambient temperature
        return TPR

    # Thermomechanical noise PSD
    # Units: V^2/Hz
    def thermo_PSD(self, freq):
        omega_damped_hz, Q_M = self.omega_damped_hz_and_Q()
        TPR = self.piezoresistor_temp()
        return (
            self.force_sensitivity() ** 2
            * 4
            * self.stiffness()
            * Cantilever.k_b
            * TPR
            / (2 * math.pi * omega_damped_hz * Q_M)
            * np.ones((1, freq.size))
        )

    # Integrated thermomechanical noise
    # Unit: V
    def thermo_integrated(self):
        [omega_damped_hz, Q_M] = self.omega_damped_hz_and_Q()
        TPR = self.piezoresistor_temp()

        return math.sqrt(
            (self.force_sensitivity()) ** 2
            * 4
            * self.stiffness()
            * Cantilever.k_b
            * TPR
            / (2 * math.pi * omega_damped_hz * Q_M)
            * (self.freq_max - self.freq_min)
        )

    # Accounts for displacement noise, vibrations, etc (V)
    def actuator_noise_integrated(self):
        return self.rms_actuator_displacement_noise * self.stiffness() * self.force_sensitivity()

    # J coefficients specify the noise floor (V/rtHz and A/rtHz)
    # F coefficients specify the 1/f noise at 1 Hz (V/rtHz and A/rtHz)
    def amplifier_noise_coefficients(self):
        if self.amplifier == "INA103":
            A_VJ = 1.2e-9
            A_IJ = 2e-12
            A_VF = 6e-9
            A_IF = 25e-12
        elif self.amplifier == "AD8221":
            A_VJ = 8e-9
            A_IJ = 40e-15
            A_VF = 12e-9
            A_IF = 550e-15
        else:
            raise RuntimeError("Unknown amplifier!")
        return A_VJ, A_IJ, A_VF, A_IF

    # Amplifier noise PSD
    # Units: V^2/Hz
    def amplifier_PSD(self, freq):
        A_VJ, A_IJ, A_VF, A_IF = self.amplifier_noise_coefficients()
        R = self.resistance() / 2  # Resistance seen at the amp inputs
        TPR = self.piezoresistor_temp()
        R *= 1 + TPR * self.TCR
        return (A_VJ**2 + 2 * (R * A_IJ) ** 2) + (A_VF**2 + 2 * (R * A_IF) ** 2) / freq

    # Integrated amplifier noise
    # Units: V
    def amplifier_integrated(self):
        A_VJ, A_IJ, A_VF, A_IF = self.amplifier_noise_coefficients()
        R = self.resistance() / 2  # Resistance seen at the amp inputs
        TPR = self.piezoresistor_temp()
        R *= 1 + TPR * self.TCR

        return math.sqrt(
            A_VJ**2 * (self.freq_max - self.freq_min)
            + A_VF**2 * math.log(self.freq_max / self.freq_min)
            + 2 * (R * A_IJ) ** 2 * (self.freq_max - self.freq_min)
            + 2 * (R * A_IF) ** 2 * math.log(self.freq_max / self.freq_min)
        )

    # Calculate the 1/f corner frequency (Hz)
    def knee_frequency(self):
        TPR = self.piezoresistor_temp()
        return (
            self.number_of_piezoresistors
            * self.alpha()
            * self.v_bridge**2
            / (16 * self.number_of_carriers() * Cantilever.k_b * TPR * self.resistance())
        )

    # Integrated cantilever noise for given bandwidth
    # Pull the calculations into this function for speed (i.e. don't
    # calculate self.resistance() five separate times
    # Units: V
    def integrated_noise(self):
        [omega_damped_hz, Q_M] = self.omega_damped_hz_and_Q()
        resistance = self.resistance() / self.gamma()
        force_sensitivity = self.force_sensitivity()
        spring_constant = self.stiffness()
        TPR = self.piezoresistor_temp()
        resistance *= 1 + TPR * self.TCR

        A_VJ, A_IJ, A_VF, A_IF = self.amplifier_noise_coefficients()

        actuator_noise_integrated = max(0, self.rms_actuator_displacement_noise * spring_constant * force_sensitivity)
        johnson_integrated = math.sqrt(4 * Cantilever.k_b * TPR * resistance * (self.freq_max - self.freq_min))

        hooge_integrated = math.sqrt(
            self.alpha()
            * self.v_bridge**2
            * self.number_of_piezoresistors
            / (4 * self.number_of_carriers())
            * math.log(self.freq_max / self.freq_min)
        )

        thermo_integrated = math.sqrt(
            force_sensitivity**2
            * 4
            * spring_constant
            * Cantilever.k_b
            * TPR
            / (2 * math.pi * omega_damped_hz * Q_M)
            * (self.freq_max - self.freq_min)
        )

        amplifier_integrated = math.sqrt(
            A_VJ**2 * (self.freq_max - self.freq_min)
            + A_VF**2 * math.log(self.freq_max / self.freq_min)
            + 2 * (resistance / 2 * A_IJ) ** 2 * (self.freq_max - self.freq_min)
            + 2 * (resistance / 2 * A_IF) ** 2 * math.log(self.freq_max / self.freq_min)
        )

        return math.sqrt(
            actuator_noise_integrated**2
            + johnson_integrated**2
            + hooge_integrated**2
            + thermo_integrated**2
            + amplifier_integrated**2
        )

    # Calculate the noise at a given frequency (V/rtHz)
    def voltage_noise(self, freq):
        return math.sqrt(
            self.johnson_PSD(freq) + self.hooge_PSD(freq) + self.thermo_PSD(freq) + self.amplifier_PSD(freq)
        )

    def f_min_cumulative(self):
        frequency = np.logspace(math.log10(self.freq_min), math.log10(self.freq_max), Cantilever.numFrequencyPoints)
        noise = self.voltage_noise(frequency)
        sensitivity = self.force_sensitivity()
        force_noise_density = noise / sensitivity
        return np.sqrt(integrate.cumulative_trapezoid(force_noise_density**2, frequency, initial=0))

    # Piezoresistance factor
    # Accounts for dopant concentration dependent piezoresistivity in silicon
    # Source: "Piezoresistance in p-type silicon revisited", Richter et al.
    def piezoresistance_factor(self, dopant_concentration):
        Nb = 6e19
        Nc = 7e20
        richter_alpha = 0.43
        richter_gamma = 1.6
        richter_beta = 0.1
        richter_eta = 3
        richter_theta = 0.9

        T0 = 300
        average_PR_temp = self.piezoresistor_temp()
        Theta = average_PR_temp / T0

        return (
            Theta**-richter_theta
            * (
                1
                + Theta**-richter_beta * (dopant_concentration / Nb) ** richter_alpha
                + Theta**-richter_eta * (dopant_concentration / Nc) ** richter_gamma
            )
            ** -1
        )

    # Low concentration longitudinal piezoresistance coefficient (1/Pa)
    def max_piezoresistance_factor(self):
        if self.doping_type == "boron":
            max_factor = 72e-11  # 110 direction
        elif self.doping_type == "phosphorus":
            max_factor = 103e-11  # 100 direction
        else:
            max_factor = 103e-11  # 100 direction
        return max_factor

    # Calculate the sensitivity factor (beta*)
    # Accounts for finite piezoresistor thickness
    # Units: None
    def beta(self):
        z, active_doping, total_doping = self.doping_profile()

        # Shift axis so that z vaies from [-t/2, t/2] and m -> cm
        z = (self.t / 2 - z) * 1e2

        TPR = self.piezoresistor_temp()
        mu, sigma = self.mobility(active_doping, TPR)

        P = self.piezoresistance_factor(active_doping)
        numerator = np.trapezoid(z, sigma * P * z)
        denominator = np.trapezoid(z, sigma)
        beta = 2 * numerator / (self.t * 1e2 * denominator)  # t: m -> cm
        return max(beta, 1e-6)  # For optimization, ensure that beta doesn't become negative

    # Ratio of piezoresistor resistance to total resistance (< 1)
    def gamma(self):
        R = self.resistance()
        return R / (R + 2 * self.R_contact)

    # Calculate the force sensitivity (V/N) for the 1/4-active Wheatstone bridge
    # Includes the transverse portion at the end of the piezoresistive loop
    def force_sensitivity(self):
        # For speed: precompute these parameters
        betaStar = self.beta()
        Rs = self.sheet_resistance()
        wheatstone_bridge_sensitivity = self.v_bridge / 4
        piMax = self.max_piezoresistance_factor()
        gamma = self.gamma()
        R = self.resistance()
        l_pr = self.l_pr()
        w_pr = self.w_pr()

        # Pick the relative transverse PR coefficient for the doping type
        transverse_factor = -1 if self.doping_type == "boron" else -0.5

        # Average length from the base to the piezoresistor centroid
        longitudinal_l_avg = self.l - l_pr / 2
        transverse_l_avg = self.l - l_pr

        # Stress prefactor
        stress_prefactor = 6 * piMax / (self.w * self.t**2) * betaStar * gamma

        # Calculate the longitudinal and transverse resistances
        # Assume that the transverse width is 2x the PR width
        R_longitudinal = 2 * Rs * l_pr / w_pr
        R_transverse = Rs * self.air_gap_width / (2 * w_pr)

        # Calculate deltaR values
        longitudinal_deltaR = stress_prefactor * longitudinal_l_avg * R_longitudinal
        transverse_deltaR = stress_prefactor * transverse_factor * transverse_l_avg * R_transverse
        deltaR_R = (longitudinal_deltaR + transverse_deltaR) / R

        return deltaR_R * wheatstone_bridge_sensitivity

    # Calculate the input referred surface stress sensitivity
    # Units: V/Pa
    def surface_stress_sensitivity(self):
        # Pick the relative transverse PR coefficient for the doping type
        if self.doping_type == "boron":
            longitudinal_factor = 1
            transverse_factor = -1
        else:
            longitudinal_factor = 1
            transverse_factor = -0.5

        # The longitudinal and transverse stress is equal everywhere
        sensitivity_factor = abs(longitudinal_factor + transverse_factor)
        return (
            9
            * sensitivity_factor
            * self.max_piezoresistance_factor()
            * self.beta()
            * self.gamma()
            * self.v_bridge
            / (16 * self.t)
        )

    # Calculate the input referred displacement sensitivity (V/m)
    def displacement_sensitivity(self):
        return self.force_sensitivity() * self.stiffness()

    # Power dissipation in the cantilever (W)
    def power_dissipation(self):
        return (self.v_bridge / 2) ** 2 / self.resistance()

    # Tip and maximum temperatures from the F-D model (K)
    def calculateMaxAndTipTemp(self):
        tmp, Q, temp = self.calculateTempProfile()
        TMax = max(temp)
        TTip = temp[-1]
        return TMax, TTip

    # Calculate the approximate PR temperature (K)
    def approxPRTemp(self):
        TMax, TTip = self.approxTempRise()
        return self.T + TMax / 2

    # Calculate the exact average PR temperature (K)
    def averagePRTemp(self):
        if self.temperature_dependent_properties == "yes":
            x, Q, temp = self.calculateTempProfileTempDependent()
        else:
            x, Q, temp = self.calculateTempProfile()

        pr_indices = np.intersect1d(np.nonzero(x >= self.l_a), np.nonzero(x <= (self.l_a + self.l_pr())))
        return self.T + np.mean(temp(pr_indices))

    # Calculate the temperature at the base of the PR (K)
    # Useful for the designing combined sensors/actuators
    def tempRiseAtPRBase(self):
        [x, Q, temp] = self.calculateTempProfile()
        base_index = np.nonzero(x >= self.l_a)[0][0]
        return temp(base_index)

    # Calculate the maximum piezoresistor temperature (K)
    def maxPRTemp(self):
        [x, Q, temp] = self.calculateTempProfile()
        pr_indices = np.intersect1d(np.nonzero(x >= self.l_a), np.nonzero(x <= (self.l_a + self.l_pr())))
        return max(temp(pr_indices))

    # Calculate the average temperature increase of the actuator (K)
    # Use for calculating A_XK (nm/K) for the thermal actuators
    def averageActuatorDeltaTemp(self):
        [x, Q, temp] = self.calculateTempProfile()
        actuator_indices = x <= (self.l_a)
        return np.mean(temp(actuator_indices))

    # Calculate the temp change of the PR in response to thermal actuation (K)
    def thermalCrosstalk(self):
        temp_hot = self.averagePRTemp()
        v_actuator_temp = self.v_actuator
        self.v_actuator = 0
        temp_cold = self.averagePRTemp()
        self.v_actuator = v_actuator_temp
        return temp_hot - temp_cold

    # Calculate the approx max and tip temperatures using lumped modeling (K)
    # Useful for quickly approximating the important temperatures
    def approxTempRiseAnalytical(self):
        k_c = self.k_base()
        R_conduction_pr = self.l_pr() / (2 * self.w * self.t * k_c)
        TMax = self.power_dissipation() * R_conduction_pr
        l_healing, tmp = self.thermalHealingLengths()
        TTip = TMax * np.exp(-(self.l - 2 / 3 * self.l_pr()) / l_healing)
        return TMax, TTip

    # Approximate temperature modeling via a circuit model
    # Much faster than the F-D model and slightly more accurate than the
    # lumped parameter model
    # Units: K
    def approxTempRise(self):
        h = self.lookupHeff()
        k_c = self.k_base()
        l_pr = self.l_pr()
        W = self.power_dissipation()

        # Model the system as current sources (PR or heater) and resistors
        [E_metal, rho_metal, k_metal, alpha_metal] = self.lookup_metal_properties()

        if self.cantilever_type == "none":
            if self.fluid == "vacuum":
                R_conduction_pr = l_pr / (2 * self.w * self.t * k_c) + self.R_base
                TMax = W * R_conduction_pr
                TTip = TMax
            else:
                R_conduction_pr = l_pr / (2 * self.w * self.t * k_c) + self.R_base
                R_convection_pr = 1 / (2 * h * l_pr * (self.w + self.t))
                R_conduction_tip = (self.l - l_pr) / (2 * self.w * self.t * k_c)
                R_convection_tip = 1 / (2 * h * (self.l - l_pr) * (self.w + self.t))
                R_total = 1 / (1 / R_conduction_pr + 1 / R_convection_pr + 1 / (R_conduction_tip + R_convection_tip))
                TMax = W * R_total
                TTip = W * R_total / (R_conduction_tip + R_convection_tip) * R_convection_tip

        elif self.cantilever_type == "step":
            R_conduction_pr = (
                l_pr / (2 * self.w * self.t * k_c)
                + self.l_a / (self.w_a * (self.t * k_c + self.t_a * k_metal))
                + self.R_base
            )
            R_convection_pr = 1 / (2 * h * (l_pr + self.l_a) * (self.w + self.t))
            R_conduction_tip = (self.l - l_pr) / (2 * self.w * self.t * k_c)
            R_convection_tip = 1 / (2 * h * (self.l - l_pr) * (self.w + self.t))
            R_total = 1 / (1 / R_conduction_pr + 1 / R_convection_pr + 1 / (R_conduction_tip + R_convection_tip))
            TMax = W * R_total
            TTip = W * R_total / (R_conduction_tip + R_convection_tip) * R_convection_tip

        elif self.cantilever_type == "piezoelectric":
            R_conduction_pr = (
                l_pr / (2 * self.w * self.t * k_c)
                + self.l_a
                / (
                    self.w_a
                    * (
                        k_c * self.t
                        + self.k_aln * (self.t_a + self.t_a_seed)
                        + k_metal * (self.t_electrode_bottom + self.t_electrode_top)
                    )
                )
                + self.R_base
            )
            R_convection_pr = 1 / (2 * h * l_pr * (self.w + self.t))
            R_conduction_tip = (self.l - l_pr) / (2 * self.w * self.t * k_c)
            R_convection_tip = 1 / (2 * h * (self.l - l_pr) * (self.w + self.t))
            R_total = 1 / (1 / R_conduction_pr + 1 / R_convection_pr + 1 / (R_conduction_tip + R_convection_tip))
            TMax = W * R_total
            TTip = W * R_total / (R_conduction_tip + R_convection_tip) * R_convection_tip

        elif self.cantilever_type == "thermal":
            R_conduction_pr = l_pr / (2 * self.w * self.t * k_c) + self.l_a / (
                self.w_a * (self.t * k_c + self.t_a * k_metal)
            )
            R_convection_pr = 1 / (2 * h * l_pr * (self.w + self.t))
            R_conduction_tip = (self.l - l_pr) / (2 * self.w * self.t * k_c)
            R_convection_tip = 1 / (2 * h * (self.l - l_pr) * (self.w + self.t))
            R_conduction_heater = self.l_a / (2 * self.w_a * (self.t * k_c + self.t_a * k_metal))
            R_convection_heater = 1 / (2 * h * self.l_a * (self.w_a + self.t_a))
            R_total = 1 / (
                1 / (R_conduction_pr + 1 / (1 / R_convection_heater + 1 / R_conduction_heater))
                + 1 / R_convection_pr
                + 1 / (R_conduction_tip + R_convection_tip)
            )

            T_heater = W / (1 / R_convection_heater + 1 / R_conduction_heater)
            TMaxDivider = (
                1
                / (1 / R_convection_pr + 1 / (R_conduction_tip + R_convection_tip))
                / (R_conduction_pr + 1 / (1 / R_convection_pr + 1 / (R_conduction_tip + R_convection_tip)))
            )
            TTipDivider = R_convection_tip / (R_convection_tip + R_conduction_tip)
            TMax = T_heater * TMaxDivider + W * R_total
            TTip = (
                T_heater * TMaxDivider * TTipDivider
                + W * R_total / (R_conduction_tip + R_convection_tip) * R_convection_tip
            )
        return TMax, TTip

    # Calculate the approx thermal healing lengths for the cantilever (m)
    def thermalHealingLengths(self):
        A = self.w * self.t
        P = 2 * (self.w + self.t)
        h = self.lookupHeff()

        k_c = self.k_base()
        l_healing_cantilever = math.sqrt(k_c * A / h / P)

        E_metal, rho_metal, k_metal, alpha_metal = self.lookup_metal_properties()
        l_healing_step = 0
        if self.cantilever_type == "step" or self.cantilever_type == "thermal":
            l_healing_step = math.sqrt(
                self.w_a
                * (k_c * self.t + self.k_sio2 * self.t_oxide + k_metal * self.t_a)
                / (2 * (self.w_a + self.t + self.t_a) * h)
            )

        elif self.cantilever_type == "piezoelectric":
            l_healing_step = math.sqrt(
                self.w_a
                * (
                    k_c * self.t
                    + self.k_sio2 * self.t_oxide
                    + self.k_aln * (self.t_a + self.t_a_seed)
                    + k_metal * (self.t_electrode_bottom + self.t_electrode_top)
                )
                / (2 * (self.w_a + self.t + self.t_a) * h)
            )
        return l_healing_cantilever, l_healing_step

    # Model the temp profile from Joule heating via finite differences
    # Assumes convection to ambient, adiabatic tip, and R_base to the
    # silicon die which is clamped at the ambient temperature
    def calculateTempProfile(self, *arg):
        # There are several ways to call calculateTempProfile()
        # No arguments: temperature independent solution
        # One argument (legacy): temp-dependent thermal conductivity
        # Two arguments: temp-dependent sheet resistance and conductivity
        if len(arg) == 1:
            k_x = arg[0]
            Rsheet = self.sheet_resistance()
            Rsheet_x = np.ones((1, self.numXPoints)) * Rsheet
        elif len(arg) == 2:
            k_x = arg[0]
            Rsheet_x = arg[1]
        else:
            k_c = self.k_base()
            Rsheet = self.sheet_resistance()
            k_x = np.ones((1, self.numXPoints)) * k_c
            Rsheet_x = np.ones((1, self.numXPoints)) * Rsheet

        # Discretize the length of the cantilever
        n_points = self.numXPoints
        totalLength = self.l + self.l_a
        dx = totalLength / (n_points - 1)
        x = np.arange(0, totalLength + 1, dx)

        # Determine the step and PR indices
        step_indices = np.nonzero(x <= self.l_a)
        actuator_indices = np.nonzero(x <= (self.l_a - self.l_a_gap))
        cantilever_indices = np.nonzero(x > self.l_a)
        pr_indices = np.intersect1d(cantilever_indices, np.nonzero(x < (self.l_a + self.l_pr())))

        # Calculate Qgen_x differently depending on our temperature range
        if self.temperature_dependent_properties == "yes":
            # Calculate Qgen(x) considering temperature dependent sheet resistance
            # This method does not converge quickly for design optimization,
            # so is best used for modeling
            index_range = np.nonzero(x <= self.l_pr())
            R_x = 2 * Rsheet_x[index_range] / self.w_pr()
            R_calc = 2 * np.trapezoid(x[index_range], Rsheet_x(index_range) / self.w_pr())
            I_calc = (self.v_bridge / 2) / R_calc
            Qgen_x = I_calc**2 * R_x
        else:
            # Assume power/length is constant along the piezoresistor length
            # This method works well for modeling near the ambient
            # temperature and for design optimization
            # Note: calculate Qgen_x based upon the actual lengths
            # to avoid discretization errors (line 2 here is very important)
            power = (self.v_bridge / 2) ** 2 / self.resistance()
            Qgen_x = power / (x[pr_indices[-1]] - x[pr_indices[1]]) * np.ones((x.size, 1))

        # Setup other variables
        tempAmbient = self.T
        h = self.lookupHeff()
        E_metal, rho_metal, k_metal, alpha_metal = self.lookup_metal_properties()
        K = self.w * np.transpose(k_x) * self.t * np.ones(((n_points, 1)))
        perimeter = 2 * (self.w + self.t) * np.ones((n_points, 1))
        Q = np.zeros((n_points, 1))
        Q[pr_indices] = Qgen_x[pr_indices]

        # Build K (area*k_c) and P
        if self.cantilever_type == "step":
            K[step_indices] = self.w_a * (k_x[step_indices] * self.t + k_metal * self.t_a)
            perimeter[step_indices] = 2 * (self.w_a + self.t_a)
        elif self.cantilever_type == "thermal":
            Qheater = self.heaterPower() / (x[actuator_indices[-1]] - x[actuator_indices[0]])
            Q[actuator_indices] = Qheater
            K[step_indices] = self.w_a * (k_x[step_indices] * self.t + self.t_oxide * self.k_sio2 + k_metal * self.t_a)
            perimeter[step_indices] = 2 * (self.w_a + self.t_a)
        elif self.cantilever_type == "piezoelectric":
            K[step_indices] = self.w_a * (
                k_x[step_indices] * self.t
                + self.k_aln * (self.t_a + self.t_a_seed)
                + self.t_oxide * self.k_sio2
                + k_metal * (self.t_electrode_top + self.t_electrode_bottom)
            )
            perimeter[step_indices] = 2 * (self.w_a + self.t_a)

        # Build A and RHS matrices
        A = np.zeros((n_points, n_points))
        rhs = np.zeros((n_points, 1))
        for ii in range(2, n_points):
            A[ii, ii - 1] = -K[ii - 1] / dx**2
            A[ii, ii] = (K[ii - 1] + K[ii + 1]) / dx**2 + h * perimeter[ii]
            A[ii, ii + 1] = -K[ii + 1] / dx**2
            rhs[ii, 1] = Q[ii] + h * perimeter[ii] * tempAmbient

        A[n_points, n_points - 1 : n_points] = [1 - 1]  # Adiabatic at tip
        # TODO A = sparse(A) # Leads to a significant speed improvement

        # Properly handle R_base
        A[0, 0] = -1 - self.R_base * K[0] / dx
        A[0, 1] = self.R_base * K[0] / dx

        rhs[0, 0] = -self.T
        T_absolute = np.linalg.solve(A, rhs)
        T_increase = T_absolute - tempAmbient

        return x, Q, T_increase

    # Model the temperature profile of a tip-heated cantilever
    # via finite differences. Intended for modeling cantilever heating
    # during laser doppler vibrometer or AFM measurements.
    # For optical heating applications the user needs to account for
    # reflections and limited silicon absorption.
    # Assumptions:
    # 1) Piezoresistor is turned off
    # 2) Temperature-independent material properties
    # 3) Beam intensity is constant over the beam diameter
    # 4) Accounts for the spot size being smaller/larger than the beam
    # 5) Spot size = spot diameter
    # Assumes temperature-independent properties for simplicity
    # Input parameters:
    # inputPower (W): power input to the cantilever
    # spotSize (m): length back from the tip over which the power is spread
    #
    # Note: most of the code is reused from above
    # TODO Refactor this method and the one above
    def calculateTempProfileTipHeat(self, inputPower, spotSize):
        k_c = self.k_base()
        k_x = np.ones((1, self.numXPoints)) * k_c
        n_points = self.numXPoints
        totalLength = self.l + self.l_a
        dx = totalLength / (n_points - 1)
        x = np.arange(0, totalLength + 1, dx)
        step_indices = np.nonzero(x <= self.l_a)
        actuator_indices = np.nonzero(x <= (self.l_a - self.l_a_gap))
        cantilever_indices = np.nonzero(x > self.l_a)
        _pr_indices = np.intersect1d(cantilever_indices, np.nonzero(x < (self.l_a + self.l_pr())))

        # Calculate the input power at the tip.
        # Accounts for some fraction of the beam missing the cantilever
        # if the cantilever width is comparable to the beam diameter
        # Assumes the center is aimed at the tip
        Q = np.zeros((n_points, 1))
        fractionAbsorbed = min((spotSize / 2) * self.w / (math.pi * (spotSize / 2) ** 2), 1)
        Q[x > (self.l - spotSize)] = fractionAbsorbed * inputPower / spotSize

        tempAmbient = self.T
        h = self.lookupHeff()
        E_metal, rho_metal, k_metal, alpha_metal = self.lookup_metal_properties()
        K = self.w * np.transpose(k_x) * self.t * np.ones((n_points, 1))
        perimeter = 2 * (self.w + self.t) * np.ones((n_points, 1))

        if self.cantilever_type == "step":
            K[step_indices] = self.w_a * (k_x[step_indices] * self.t + k_metal * self.t_a)
            perimeter[step_indices] = 2 * (self.w_a + self.t_a)
        elif self.cantilever_type == "thermal":
            Qheater = self.heaterPower() / (x[actuator_indices[-1]] - x[actuator_indices[0]])
            Q[actuator_indices] = Qheater
            K[step_indices] = self.w_a_active * (
                k_x[step_indices] * self.t + self.t_oxide * self.k_sio2 + k_metal * self.t_a
            )
            perimeter[step_indices] = 2 * (self.w_a + self.t_a)
        elif self.cantilever_type == "piezoelectric":
            K[step_indices] = self.w_a * (
                k_x[step_indices] * self.t
                + self.k_aln * (self.t_a + self.t_a_seed)
                + self.t_oxide * self.k_sio2
                + k_metal * (self.t_electrode_top + self.t_electrode_bottom)
            )
            perimeter[step_indices] = 2 * (self.w_a + self.t_a)

        A = np.zeros((n_points, n_points))
        rhs = np.zeros((n_points, 1))
        for ii in range(2, n_points):
            A[ii, ii - 1] = -K[ii - 1] / dx**2
            A[ii, ii] = (K[ii - 1] + K[ii + 1]) / dx**2 + h * perimeter[ii]
            A[ii, ii + 1] = -K[ii + 1] / dx**2
            rhs[ii, 1] = Q[ii] + h * perimeter[ii] * tempAmbient
        A[n_points, n_points - 1 : n_points] = [1 - 1]
        # A = sparse(A)
        A[0, 0] = -1 - self.R_base * K[0] / dx
        A[0, 1] = self.R_base * K[0] / dx
        rhs[0, 0] = -self.T
        T_absolute = np.linalg.solve(A, rhs)
        T_increase = T_absolute - tempAmbient

        return x, Q, T_increase

    # Calculate the temp profile accounting for k(T) and R_s(T)
    # If temperature_dependent_properties is set to 'yes' then this
    # function should be called correctly. The temperature dependent option
    # is much slower than its temperature independent equivalent.
    def calculateTempProfileTempDependent(self):
        # Iterative calculate the temperature distribution and k(x)
        # Initialize k(x) solver with our initial guess
        # This approach is required because the system is nonlinear
        x, Q, T_initial = self.calculateTempProfile()
        T_current = T_initial
        residual = 1

        while residual > 1e-6:
            absoluteTemp = T_current + self.T

            # Calculate k_si and Rsheet from the current temperature guess
            x, z, k_z_x = self.thermal_conductivity_profile(absoluteTemp)
            k_x = np.mean(k_z_x, 1)  # Average k for a cross-section
            Rsheet_x = self.RSheetProfile(x, T_current)
            pr_indices = np.nonzero(x <= self.l_pr())

            # Scale the resistance so that it matches our nominal value
            # This choice is geared towards modeling experimental results
            R_nominal = self.resistance()
            x_resistor = x[pr_indices]
            Rsheet_resistor = Rsheet_x[pr_indices]
            R_calc = (
                2 * np.trapezoid(x_resistor, Rsheet_resistor / self.w_pr())
                + 3.4 * Rsheet_resistor[-1]
                + 2 * self.R_contact
            )
            Rsheet_x = Rsheet_x * (R_nominal / R_calc)

            # Update the temperature guess using k_x, Rsheet_x
            [x, Q, T_new] = self.calculateTempProfile(k_x, Rsheet_x)

            residual = sum((T_new - T_current) ** 2)
            T_current = T_new
        T_final = T_current
        return x, Q, T_final

    # Calculate the approximate thermal conductivity of the cantilever (W/m-K)
    def k_base(self):
        k_x = self.k_x()
        return k_x[0]

    # Calculate k(x) as a function of cantilever thickness
    # None/approx = Asheghi (1997) model
    # Units: W/m-K
    def k_x(self):
        if self.thermal_modeling == "none" or self.thermal_modeling == "approx":
            t0 = 120e-9
            k_effective = self.k_si * self.t / (t0 + self.t)
        else:
            [x, z, k] = self.thermal_conductivity_profile(self.T)
            k_effective = np.mean(k)
        return k_effective

    # Calculate k(x,z) using a temp/thickness dependent model
    # Source: "Modeling and Data for Thermal Conductivity ...", Liu et al, 2006.
    # Units: W/m-K
    def thermal_conductivity_profile(self, T):
        # Phonon Group velocities (m/sec)
        nu_T = 5860
        nu_TU = 2000
        nu_L = 8480

        # Density of state coefficients
        _C_T = Cantilever.k_b / (2 * math.pi**2 * nu_T) * (Cantilever.k_b / self.h_bar) ** 3
        C_TU = Cantilever.k_b / (2 * math.pi**2 * nu_TU) * (Cantilever.k_b / self.h_bar) ** 3
        C_L = Cantilever.k_b / (2 * math.pi**2 * nu_L) * (Cantilever.k_b / self.h_bar) ** 3

        # Debye temperatures (K)
        theta_1 = 180
        theta_2 = 210
        theta_3 = 570

        # Phonon specific heat/volume (J/m^3-K)
        c_TU = 7.5e5
        c_L = 3.2e5

        # Other Constants
        B_T = 9.3e-13  # 1/K^3
        B_TU = 5.5e-18  # s
        B_L = 2e-24  # s/K
        beta_TU = B_TU * (Cantilever.k_b / self.h_bar) ** 2  # 1/s-K^2
        beta_L = B_L * (Cantilever.k_b / self.h_bar) ** 2  # 1/s-K^3
        _beta_T = B_T * (Cantilever.k_b / self.h_bar)  # 1/s-K^4 (unused, for reference)

        # Phonon frequencies
        omega_2 = Cantilever.k_b / self.h_bar * theta_2
        omega_3 = Cantilever.k_b / self.h_bar * theta_3

        # Crystal lattice parameters
        Gamma = 2.1e-4  # computed from isotopes
        nu_s = (1 / 3 * (nu_L**-1 + 2 * nu_T**-1)) ** -1  # Speed of sound

        R_si = 5.4e-10  # m - radius of the silicon host atom
        V_si = (R_si / 2) ** 3  # m^3 - silicon lattice volume
        M_si = 4.59e-26  # kg - mass of silicon host atom
        A_isotopes = Gamma * V_si / (4 * math.pi * nu_s**3)

        # Other constants
        Q = 4.2  # Constant scattering matrix parameter
        gamma = 1.6  # Gruneisen parameter
        psi = 3.5  # Modification factor to account for inconsistency of the Si

        # Doping constants
        if self.doping_type == "boron":
            dM = 18 * 1.673e-27
            dR = 0.175e-10
        elif self.doping_type == "phosphorus":
            dM = 3 * 1.673e-27
            dR = 0.175e-10
        else:
            raise RuntimeError("Arsenic dM and dR values not defined!")

        # Vector calculation of thermal conductivity
        # Roughly 10x faster than iterative calculation!
        z, active_doping, total_doping = self.doping_profile()
        x = np.linspace(0, self.l, self.numXPoints)

        # let dim(T) and dim(total_doping) be numZPoints x numXPoints
        n_matrix = 1e6 * np.transpose(total_doping) * np.ones((1, self.numXPoints))  # convert from 1/cc to 1/m^3
        if T.size == 1:
            T_matrix = np.ones((self.numZPoints, 1)) * np.transpose(T) * np.ones((1, self.numXPoints))
        else:
            T_matrix = np.ones((self.numZPoints, 1)) * np.transpose(T)

        # Calculate thermal conductivity
        lambda_TU = C_TU * (theta_2**2 - theta_1**2) / (2 * nu_TU * c_TU * beta_TU * T_matrix)
        lambda_TU_HT = (1 / lambda_TU + 1 / (nu_TU * (A_isotopes * omega_2**4) ** -1)) ** -1
        A = 1 / lambda_TU_HT + 1 / (
            nu_TU
            * (
                (
                    n_matrix * V_si**2 / (4 * math.pi * nu_s**3) * (dM / M_si) ** 2
                    + 2 * n_matrix * V_si**2 * Q**2 * gamma**2 / (math.pi * nu_s**3) * (dR / R_si) ** 2
                )
                * omega_2**4
            )
            ** -1
        )
        k_D_layer_TU_HT = 2 / 3 * c_TU * nu_TU * A**-1 * (1 + 3 / (8 * self.t / (A**-1 * psi))) ** -1

        lambda_L = C_L * theta_3 / (nu_L * c_L * beta_L * T_matrix**3)
        lambda_L_HT = (1 / lambda_L + 1 / (nu_L * (A_isotopes * omega_3**4) ** -1)) ** -1
        A = 1 / lambda_L_HT + 1 / (
            nu_L
            * (
                (
                    n_matrix * V_si**2 / (4 * math.pi * nu_s**3) * (dM / M_si) ** 2
                    + 2 * n_matrix * V_si**2 * Q**2 * gamma**2 / (math.pi * nu_s**3) * (dR / R_si) ** 2
                )
                * omega_3**4
            )
            ** -1
        )
        k_D_layer_L_HT = 1 / 3 * c_L * nu_L * A**-1 * (1 + 3 / (8 * self.t / (A**-1 * psi))) ** -1

        k = k_D_layer_L_HT + k_D_layer_TU_HT

        return x, z, k

    # ====== Calculate resolution ======
    def force_resolution(self):
        return self.integrated_noise() / self.force_sensitivity()

    def thermomechanical_force_noise(self):
        return self.thermo_integrated() / self.force_sensitivity()

    def displacement_resolution(self):
        return self.force_resolution() / self.stiffness()

    def force_noise_density(self, freq):
        return self.voltage_noise(freq) / self.force_sensitivity()

    def resonant_force_noise_density(self):
        [omega_damped_hz, Q] = self.omega_damped_hz_and_Q()  ##ok<*NASGU>
        freq = omega_damped_hz
        return self.force_noise_density(freq) * 1e15

    def surface_stress_resolution(self):
        return self.integrated_noise() / self.surface_stress_sensitivity()

    # ====== Multilayer beam mechanics and actuation ======

    # Calculate the cantilever curvature per unit applied moment
    def calculateActuatorNormalizedCurvature(self):
        Zm = self.actuatorNeutralAxis()
        [z, E, A, I] = self.lookupActuatorMechanics()
        Z_offset = z - Zm
        return 1 / sum(E * (I + A * Z_offset**2))

    # Calculate the thickness of silicon required to have the same EI
    # as the actuator/step at the base
    def calculateEquivalentThickness(self):
        t_equivalent = optimize.fminbound(self.findEIResidual, self.t / 100, self.t * 100, "xtol", 1e-12)
        return t_equivalent

    def actuatorNeutralAxis(self):
        [z, E, A, I] = self.lookupActuatorMechanics()
        return sum(z * E * A) / sum(E * A)

    def d31(self):
        if self.d31_manual == 0:
            d31 = interpolate.interp1d(self.d31_t, self.d31_aln, self.t_a, "spline")
        else:
            d31 = self.d31_manual
        return d31

    def calculateDeflection(self):
        n_points = self.numXPoints
        totalLength = self.l + self.l_a
        dx = totalLength / (n_points - 1)
        x = np.arange(0, totalLength + 1, dx)

        M = 0  # external moment is zero
        P = 0  # external load is zero

        z, E, A, I = self.lookupActuatorMechanics()
        stress = self.calculateActuatorStress()

        # Calculate the curvature and neutral axis
        # The curvature may vary as a function of position (e.g. thermal)
        # so calculate the deflection by integrating the curvature twice
        C = np.zeros((x.size, 1))

        # Calculate the curvature, C along the cantilever length
        for ii in range(1, x.size + 1):
            C[ii] = -(
                (M - np.sum(z * A * stress[ii, :])) * np.sum(E * A) + (P + sum(A * stress[ii, :])) * np.sum(E * z * A)
            ) / (np.sum(E * A) * np.sum(E * (I + A * z**2)) - np.sum(z * E * A) ** 2)

        # No curvature beyond the end of the actuator
        # i.e. no dopant bending included in this calculation
        C[x > self.l_a] = 0

        # Calculate the deflection from the curvature by integrating
        theta = integrate.cumulative_trapezoid(C, x, initial=0)
        deflection = integrate.cumulative_trapezoid(np.tan(theta), x, initial=0)
        return x, deflection

    # Calculate the stress as a function along the actuator length
    def calculateActuatorStress(self):
        z, E, A, I = self.lookupActuatorMechanics()
        E_metal, rho_metal, k_metal, alpha_metal = self.lookup_metal_properties()
        x_temp, Q, temp = self.calculateTempProfile()
        temp = temp + self.T  # Use absolute temperature for these calculations

        intrinsic_stress = np.ones((self.numXPoints, 1)) * self.film_intrinsic_stress()
        intrinsic_stress[x_temp > self.l_a, :] = 0

        if self.cantilever_type == "step" or self.cantilever_type == "thermal":
            cte = np.array((self.alpha_si, self.alpha_sio2, alpha_metal))
            thermal_stress = (temp - self.T_ref) * np.ones((self.numXPoints, 1)) * (cte * E)
            thermal_stress[x_temp > self.l_a, 2:3] = 0
            piezoelectric_stress = np.zeros((self.numXPoints, 3))
        elif self.cantilever_type == "piezoelectric":
            cte = np.array((self.alpha_si, self.alpha_sio2, self.alpha_aln, alpha_metal, self.alpha_aln, alpha_metal))
            thermal_stress = (temp - self.T_ref) * np.ones((self.numXPoints, 1)) * (cte * E)
            thermal_stress[x_temp > self.l_a, 2:6] = 0

            # Calculate the piezoelectric stress
            # The seed electric field will vary depending on
            E_field = np.array((0, 0, 0, 0, self.v_actuator / self.t_a, 0))
            d31 = np.array((0, 0, 0, 0, self.d31(), 0)).transpose()

            piezoelectric_stress = np.ones((self.numXPoints, 1)) * E * d31 * E_field
            piezoelectric_stress[x_temp > self.l_a - self.l_a_gap, :] = 0
        else:
            raise RuntimeError("Unknown cantilever_type")

        # Scale the stress by the active width of the device
        return intrinsic_stress + thermal_stress + self.w_a_active / self.w_a * piezoelectric_stress

    def tip_deflection_distribution(self):
        v_actuator_temporary = self.v_actuator
        self.v_actuator = 0
        self.film_stress = "random"
        z = np.zeros((Cantilever.numRandomStressIterations, 1))
        z_tip = np.zeros((Cantilever.numRandomStressIterations, 1))

        for ii in range(Cantilever.numRandomStressIterations):
            _, z_ii = self.calculateDeflection()
            z[:, ii] = z_ii
            z_tip[ii] = self.tipDeflection()

        self.v_actuator = v_actuator_temporary
        self.film_stress = "nominal"

        mu_z = np.mean(z_tip)
        sigma_z = np.std(z_tip)
        return mu_z, sigma_z

    def plot_tip_deflection_distribution(self):
        v_actuator_temporary = self.v_actuator
        self.v_actuator = 0
        self.film_stress = "random"
        z = np.zeros((Cantilever.numRandomStressIterations, 1))

        for ii in range(Cantilever.numRandomStressIterations):
            x_vals, z_ii = self.calculateDeflection()
            z[:, ii] = z_ii
        self.v_actuator = v_actuator_temporary
        self.film_stress = "nominal"

        plt.figure()
        plt.plot(x_vals * 1e6, z * 1e6)
        plt.xlabel(r"Distance from base ($\mu$m)")
        plt.ylabel(r"Deflection ($\mu$m)")
        plt.box(False)

    def film_intrinsic_stress(self):
        # film_stress = 'nominal' => use average stress
        # film_stress = 'random' => use random, normally distributed value (sigma = 1/4 so that user input = 2 sigma)
        sigma = 0
        if self.film_stress == "random":
            sigma = 1 / 4

        if self.cantilever_type == "step" or self.cantilever_type == "thermal":
            sigma_i = np.zeros((3, 1))
            sigma_i[0] = np.mean(self.sigma_si_range) + sigma * np.random.randn() * np.diff(self.sigma_si_range)
            sigma_i[1] = np.mean(self.sigma_sio2_range) + sigma * np.random.randn() * np.diff(self.sigma_sio2_range)
            sigma_i[2] = np.mean(self.sigma_al_range) + sigma * np.random.randn() * np.diff(self.sigma_al_range)
        elif self.cantilever_type == "piezoelectric":
            sigma_i = np.zeros((5, 1))
            sigma_i[0] = np.mean(self.sigma_si_range) + sigma * np.random.randn() * np.diff(self.sigma_si_range)
            sigma_i[1] = np.mean(self.sigma_sio2_range) + sigma * np.random.randn() * np.diff(self.sigma_sio2_range)
            sigma_i[2] = np.mean(self.sigma_aln_range) + sigma * np.random.randn() * np.diff(self.sigma_aln_range)
            sigma_i[4] = np.mean(self.sigma_aln_range) + sigma * np.random.randn() * np.diff(self.sigma_aln_range)
            if self.metal_type == "titanium":
                sigma_i[3] = np.mean(self.sigma_ti_range) + sigma * np.random.randn() * np.diff(self.sigma_ti_range)
                sigma_i[5] = np.mean(self.sigma_ti_range) + sigma * np.random.randn() * np.diff(self.sigma_ti_range)
            elif self.metal_type == "molybdenum":
                sigma_i[3] = np.mean(self.sigma_mo_range) + sigma * np.random.randn() * np.diff(self.sigma_mo_range)
                sigma_i[5] = np.mean(self.sigma_mo_range) + sigma * np.random.randn() * np.diff(self.sigma_mo_range)
        return sigma_i

    # The thermal actuator power dissipation (W)
    def heaterPower(self):
        return self.v_actuator**2 / self.R_heater

    # Calculate the time constant and -3 dB freq of the thermal actuator
    # Note that in practice this is roughly 4x too idealistic
    # Units: s and Hz
    def heaterTimeConstant(self):
        k_c = self.k_base()
        E_metal, rho_metal, k_metal, alpha_metal = self.lookup_metal_properties()
        h = self.lookupHeff()

        # Calculate the thermal resistances
        R_conduction = (self.l_a - self.l_a_gap) / (
            2 * self.w_a * (self.t * k_c + self.t_a * k_metal + self.t_oxide * self.k_sio2)
        )
        R_convection = 1 / (2 * h * self.l_a * (self.w_a + self.t_a))
        R = 1 / (1 / (R_conduction + self.R_base) + 1 / R_convection)  # K/W
        C = (
            self.l_a
            * self.w_a
            * (
                self.t_a * self.rho_al * self.Cv_al
                + self.t * self.rho_si * self.Cv_si
                + self.t_oxide * self.rho_sio2 * self.Cv_sio2
            )
        )  # J/K

        tau = R * C
        freq = 1 / (2 * math.pi * tau)
        return tau, freq

    # Calculate the rise time using the Q-dependent formula
    # Source: "Control Systems Engineering", Nise
    def mechanicalRiseTime(self):
        omega_d, Q = self.omega_damped_and_Q()
        zeta = 1 / (2 * Q)
        return 1 / omega_d * (1.76 * zeta**3 - 0.417 * zeta**2 + 1.039 * zeta + 1)

    # Calculate dz_tip/dT_ambient -> temperature stability of the cantilever (m/K)
    def tipTempDeflectionVsAmbientTemp(self):
        temp_delta_values = np.arange(-5, 6)

        tipDeflection = np.zeros((temp_delta_values.size, 1))
        T_original = self.T
        for ii in range(temp_delta_values.size):
            self.T = T_original + temp_delta_values[ii]
            tipDeflection[ii] = self.tipDeflection()
        self.T = T_original

        p = np.polyfit(temp_delta_values, tipDeflection, 1)
        return p[0]

    # Calculate the electrical current in the heater (A)
    def heaterCurrent(self):
        return self.v_actuator / self.R_heater

    # Calculate the tip deflection including all of the various effects (m)
    def tipDeflection(self):
        [x, z] = self.calculateDeflection()
        return z[-1]

    # Calculate the tip deflection with the actuator on vs. off (m)
    # This approach cancels any static deflections (e.g. dopant bending)
    def actuatorTipDeflectionRange(self):
        v_actuator_temp = self.v_actuator
        self.v_actuator = 0
        x, z_off = self.calculateDeflection()
        z_initial = z_off[-1]

        self.v_actuator = v_actuator_temp
        x, z_on = self.calculateDeflection()
        z_final = z_on[-1]

        return np.abs(z_final - z_initial)

    # Calculate dopant strain induced bending
    # The calculation method handles arbitrary strain gradients
    # Source: "Residual strain and resultant postrelease deflection ...", 1999
    def calculateDopantBending(self):
        # Source: "Effects of phosphorus on stress...", JMM, 1999.
        delta = 4.5e-24

        n_points = self.numXPoints
        dx = self.l / (n_points - 1)
        x = np.arange(0, self.l + 1, dx)

        [depth, active_doping, total_doping] = self.doping_profile()
        E = self.modulus()
        width = self.w  # But it doesn't matter
        thickness = self.t

        # Shift coordinates so that strain is measured from the bottom
        y = np.linspace(0, thickness, depth.size).transpose()
        dopant_concentration = np.flipud(active_doping)

        strain = dopant_concentration * delta
        stress = strain * E

        ybar = np.trapezoid(y, y * stress.transpose()) / np.trapezoid(y, stress.transpose())
        h1 = thickness - ybar
        h2 = ybar

        i1 = np.nonzero(y > ybar)
        i2 = np.nonzero(y <= ybar)
        y1 = y(i1)
        y2 = y(i2)
        stress1 = stress[i1]
        stress2 = stress[i2]

        sigma1 = 1 / h1 * np.trapezoid(y1, stress1)
        sigma2 = 1 / h2 * np.trapezoid(y2, stress2)

        ybar1 = np.trapezoid(y1, y1 * stress1.transpose()) / np.trapezoid(y1, stress1.transpose())
        ybar2 = np.trapezoid(y2, y2 * stress2.transpose()) / np.trapezoid(y2, stress2.transpose())

        I1 = width * h1**3 / 12
        I2 = width * h2**3 / 12
        I = width * thickness**3 / 12

        sigmai = (sigma2 - sigma1) / (2 - h1 * width * (ybar1 - ybar) ** 2 / I1 + h2 * width * (ybar - ybar2) ** 2 / I2)

        M1 = sigmai * h1 * width * (ybar1 - ybar)
        M2 = sigmai * h2 * width * (ybar - ybar2)

        # Calculate the radius of curvature, assuming it is small
        R = (E * I) / (M1 + M2)
        C = np.zeros((x.size, 1))
        C_pr = 1 / R
        C[x <= self.l_pr()] = C_pr

        theta = np.cumsum(C * dx)
        deflection = np.cumsum(theta * dx)

        return x, C, theta, deflection

    # Calculate just the tip deflection from dopant bending
    # Units: m
    def calculateDopantTipDeflection(self):
        x, C, theta, deflection = self.calculateDopantBending()
        tipDeflection = max(deflection)
        return tipDeflection

    # ====== Handy plotting functions ======
    def plotTempProfile(self):
        if self.temperature_dependent_properties == "yes":
            x, Q, temp = self.calculateTempProfileTempDependent()
        else:
            x, Q, temp = self.calculateTempProfile()

        plt.figure()
        plt.plot(1e6 * x, temp)
        plt.plot(1e6 * x, Q, "--")
        plt.xlabel("X (um)")
        plt.ylabel("Temp Rise (K)")
        plt.legend("Temperature", "Power/Unit Length")

    def plotTempProfileActuatorComparison(self):
        x, Q, temp = self.calculateTempProfile()

        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(1e6 * x, Q, "--")
        plt.ylabel(r"Power ($\mu$W/$\mu$m)")

        plt.subplot(2, 1, 2)
        plt.plot(1e6 * x, temp)

        v_actuator_temp = self.v_actuator
        self.v_actuator = 0
        [x, Q, temp] = self.calculateTempProfile()
        self.v_actuator = v_actuator_temp

        plt.plot(1e6 * x, temp)
        plt.xlabel(r"Distance from base ($\mu$m)")
        plt.ylabel(r"$\Delta$T (K)")
        plt.box(False)

    def plot_thermal_conductivity_profile(self):
        x, z, k = self.thermal_conductivity_profile(self.T)

        plt.figure()
        plt.plot(z * 1e9, k)
        plt.xlabel("Distance from Surface (nm)")
        plt.ylabel("k (W/m-K)")
        plt.gca().set_ylim([0, 150])
        plt.gca().set_xlim([0, self.t * 1e9])

    def plotDopantBending(self):
        x, C, theta, deflection = self.calculateDopantBending()

        plt.figure()
        plt.plot(1e6 * x, 1e9 * deflection)
        plt.xlabel("Distance from Base (um)")
        plt.ylabel("Cantilever Deflection (nm)")

    def plotDeflectionAndTemp(self):
        v_actuator_temp = self.v_actuator
        self.v_actuator = 0
        [x, deflectionInitial] = self.calculateDeflection()

        self.v_actuator = v_actuator_temp
        x, deflectionFinal = self.calculateDeflection()
        x, Q, temp = self.calculateTempProfile()
        deflection = deflectionFinal - deflectionInitial

        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(1e6 * x, temp)
        plt.ylabel("Temp Rise (K)")

        plt.subplot(2, 1, 2)
        plt.plot(1e6 * x, 1e6 * deflection)
        plt.xlabel("Distance from Base (um)")
        plt.ylabel(r"Deflection ($\mu$m)")

    def plot_noise_spectrum(self):
        freq = np.logspace(math.log10(self.freq_min), math.log10(self.freq_max), Cantilever.numFrequencyPoints)

        plt.figure()
        plt.plot(freq, math.sqrt(self.johnson_PSD(freq)))
        plt.plot(freq, math.sqrt(self.hooge_PSD(freq)))
        plt.plot(freq, math.sqrt(self.thermo_PSD(freq)))
        plt.plot(freq, math.sqrt(self.amplifier_PSD(freq)))
        plt.plot(freq, self.voltage_noise(freq))
        plt.gca().set_xscale("log")
        plt.gca().set_yscale("log")
        plt.ylabel("Noise Voltage Spectral Density (V/rtHz)")
        plt.xlabel("Frequency (Hz)")

    # ======= Handy lookup functions =======
    # Calculate elastic modulus based upon dopant type
    # Assume we're using the best piezoresistor orientation.
    def modulus(self):
        if self.doping_type == "boron":
            return 169e9
        elif self.doping_type in ("phosphorus", "arsenic"):
            return 130e9
        else:
            raise RuntimeError("Unknown dopant type!")

    # Lookup the appropriate fluid physical properties
    def lookupFluidProperties(self):
        if self.fluid == "air":
            return self.rho_air, self.eta_air
        elif self.fluid == "water":
            return self.rho_water, self.eta_water
        elif self.fluid == "arbitrary":
            return self.rho_arb, self.eta_arb
        else:
            raise RuntimeError("Unknown fluid type!")

    # Lookup effective convection coefficient
    # Or try calculating it using a shape factor (not recommended)
    # Units: W/m^2-K
    def lookupHeff(self):
        if self.h_method == "fixed":
            if self.fluid == "vacuum":
                return self.h_vacuum
            elif self.fluid == "air":
                return self.h_air
            elif self.fluid == "water":
                return self.h_water
            elif self.fluid == "arbitrary":
                return self.h_arb
            else:
                raise RuntimeError("Unknown fluid!")
        elif self.h_method == "calculate":
            if self.fluid == "vacuum":
                k_f = self.k_vacuum
            elif self.fluid == "air":
                k_f = self.k_air
            elif self.fluid == "water":
                k_f = self.k_water
            elif self.fluid == "arbitrary":
                k_f = self.k_arb
            else:
                raise RuntimeError("Unknown fluid!")

            # Calculate h using the shape factor
            P = 2 * (self.w + self.t)

            # Assume a thin rectangle to an infinite plane
            z = 450e-6  # cantilever die thickness, m
            S = 2 * math.pi * self.l / math.log(2 * math.pi * z / self.w)

            return k_f * S / (P * self.l)  # calculate h from the Hu equation
        else:
            raise RuntimeError("Unknown h_method!")

    # Simple lookup function for the film stack properties
    def lookup_metal_properties(self):
        if self.metal_type == "aluminum":
            E = self.E_al
            rho = self.rho_al
            k = self.k_al
            alpha = self.alpha_al
        elif self.metal_type == "titanium":
            E = self.E_ti
            rho = self.rho_ti
            k = self.k_ti
            alpha = self.alpha_ti
        elif self.metal_type == "molybdenum":
            E = self.E_mo
            rho = self.rho_mo
            k = self.k_mo
            alpha = self.alpha_mo
        else:
            raise RuntimeError("Unknown metal_type!")
        return E, rho, k, alpha

    # Treat the actuator as encompassing the entire device width here
    # and scale the actuation later (to account for the metal traces
    # running alongside
    def lookupActuatorMechanics(self):
        E_metal, rho_metal, k_metal, alpha_metal = self.lookup_metal_properties()

        if self.cantilever_type == "step" or self.cantilever_type == "thermal":
            t_layers = np.array((self.t, self.t_oxide, self.t_a))
            w_layers = np.array((self.w_a, self.w_a, self.w_a))
            E_layers = np.array((self.modulus(), self.E_sio2, E_metal))
        elif self.cantilever_type == "piezoelectric":
            t_layers = np.array(
                (self.t, self.t_oxide, self.t_a_seed, self.t_electrode_bottom, self.t_a, self.t_electrode_top)
            )
            w_layers = np.array((self.w_a, self.w_a, self.w_a, self.w_a, self.w_a, self.w_a))
            E_layers = np.array((self.modulus(), self.E_sio2, self.E_aln, E_metal, self.E_aln, E_metal))
        else:
            raise RuntimeError("Unknown cantilever_type")

        z_layers = np.zeros((1, t_layers.size))
        for ii in range(1, t_layers.size + 1):
            # z(1) = t(1)/2, z(2) = t(1) + t(2)/2
            z_layers[ii] = np.sum(t_layers) - np.sum(t_layers[ii:]) + t_layers[ii] / 2
        A = w_layers * t_layers
        I = (w_layers * t_layers**3) / 12
        return z_layers, E_layers, A, I

    # ======== Beam mechanics, resonant frequency and damping ==========

    # Calculate the cantilever spring constant
    # Units: N/m
    def stiffness(self):
        k_tip = self.modulus() * self.w * self.t**3 / (4 * self.l**3)

        # For devices with an actuator/step at the base, appply a test force and calculate k = F/x
        # Treating the system as two springs in series is not valid because:
        # 1) the thick base sees a larger moment than in the springs-in-series approach
        # 2) the angle at the end of the thick base is integrated along the remaining length of the cantilever
        if self.cantilever_type == "none" and self._beam_gap_config is None:
            stiffness = k_tip
        elif self.cantilever_type == "none" and self._beam_gap_config is not None:
            # For beam with gap, calculate stiffness numerically using spatially varying EI
            stiffness = self._calculate_stiffness_with_beam_gap()
        else:
            F = 1e-9  # Test force (N)
            EI_base = 1 / self.calculateActuatorNormalizedCurvature()
            EI_tip = self.modulus() * self.w * self.t**3 / 12

            # Deflection at the tip of the base
            z_base = F * self.l_a**2 * (3 * (self.l + self.l_a) - self.l_a) / (6 * EI_base)

            # Angle at the tip of the base
            theta_base = F * (2 * (self.l_a + self.l) - self.l_a) * self.l_a / (2 * EI_base)

            # Deflection at the very tip
            z_tip = z_base + np.tan(theta_base) * self.l + F * self.l**2 * 2 * self.l / (6 * EI_tip)
            stiffness = F / z_tip
        return stiffness

    def _calculate_stiffness_with_beam_gap(self) -> float:
        """Calculate stiffness for a cantilever with beam gap using numerical integration.

        For a beam with spatially varying width, we integrate 1/EI(x) along the length
        to find the tip deflection under a unit load.

        Returns:
            Spring constant (N/m)
        """
        # Discretize the cantilever length
        n_points = self.numXPoints
        x = np.linspace(0, self.l, n_points)
        dx = x[1] - x[0]

        # Get effective width at each position
        w_eff = self.effective_beam_width_array(x)

        # Calculate EI at each position
        EI = self.modulus() * w_eff * self.t**3 / 12

        # For a point load F at the tip, M(x) = F * (L - x)
        # Deflection: z'' = M/EI = F*(L-x)/EI(x)
        # Integrate twice to get deflection at tip
        F = 1.0  # Unit force

        # First integration: slope
        M_over_EI = F * (self.l - x) / EI
        theta = np.zeros(n_points)
        for i in range(1, n_points):
            theta[i] = theta[i - 1] + M_over_EI[i] * dx

        # Second integration: deflection
        z = np.zeros(n_points)
        for i in range(1, n_points):
            z[i] = z[i - 1] + theta[i] * dx

        # Stiffness = F / z_tip
        z_tip = z[-1]
        if z_tip > 0:
            return F / z_tip
        # Fallback to uniform beam formula
        return self.modulus() * self.w * self.t**3 / (4 * self.l**3)

    # Effective mass of the cantilever beam
    # Does not include the effective mass of the base/actuator section
    # f0 is calculated using Rayleight-Ritz there
    # Units: kg
    def effective_mass(self):
        if self._beam_gap_config is None:
            cantilever_effective_mass = 0.243 * self.rho_si * self.w * self.t * self.l
        else:
            # For beam with gap, calculate effective mass accounting for varying width
            cantilever_effective_mass = self._calculate_effective_mass_with_beam_gap()
        effective_mass = cantilever_effective_mass + self.tip_mass
        return effective_mass

    def _calculate_effective_mass_with_beam_gap(self) -> float:
        """Calculate effective mass for a cantilever with beam gap.

        The effective mass accounts for the mode shape distribution.
        For a simple cantilever, m_eff = 0.243 * m_total.
        With varying width, we integrate the mass distribution weighted by mode shape.

        Returns:
            Effective mass (kg)
        """
        # Discretize the cantilever length
        n_points = self.numXPoints
        x = np.linspace(0, self.l, n_points)

        # Get effective width at each position
        w_eff = self.effective_beam_width_array(x)

        # Mass per unit length at each position
        dm_dx = w_eff * self.t * self.rho_si

        # First-mode shape for a cantilever beam (approximate)
        # phi(x) ~ (x/L)^2 * (3 - x/L) for a point-loaded beam
        # Normalized so phi(L) = 1
        xi = x / self.l
        phi = 1.5 * xi**2 - 0.5 * xi**3

        # Effective mass = integral of dm_dx * phi^2
        m_eff = np.trapezoid(dm_dx * phi**2, x)

        return m_eff

    # Resonant frequency for undamped vibration (first mode).
    # For a simple cantilever, use Bernoulli beam theory.
    # For step/actuator designs or beams with gap, use Rayleigh-Ritz (i.e. U_kin = U_strain)
    # The R-R results agree with 2D and 3D FEA results to <3# when
    # (1 < t_a/t < 6, 1/5 < l_a/l < 5)
    # Units: radians/sec
    def omega_vacuum(self):
        omega_bernoulli = math.sqrt(self.stiffness() / self.effective_mass())

        if self.cantilever_type == "none" and self._beam_gap_config is None:
            omega_vacuum = omega_bernoulli
        elif self.cantilever_type == "none" and self._beam_gap_config is not None:
            # For beam with gap but no actuator, use Rayleigh-Ritz
            omega_vacuum = self._omega_vacuum_with_beam_gap()
        else:
            omega_vacuum = optimize.fminbound(self.findEnergyResidual, 1, 100 * omega_bernoulli, "xtol", 1e-6)
        return omega_vacuum

    def _omega_vacuum_with_beam_gap(self) -> float:
        """Calculate vacuum resonant frequency for beam with gap using Rayleigh-Ritz.

        For a beam with spatially varying width, we use the Rayleigh-Ritz method
        to find the fundamental frequency by equating kinetic and potential energy.

        Returns:
            Resonant frequency (rad/s)
        """
        # Discretize the cantilever length
        n_points = self.numXPoints
        x = np.linspace(0, self.l, n_points)

        # Get effective width at each position
        w_eff = self.effective_beam_width_array(x)

        # EI(x) and dm/dx(x)
        EI = self.modulus() * w_eff * self.t**3 / 12
        dm_dx = w_eff * self.t * self.rho_si

        # Trial mode shape: static deflection under tip load
        # phi''(x) = M(x)/EI(x) = (L-x)/EI(x), integrated twice
        # We use a simplified polynomial that satisfies BCs: phi(0)=0, phi'(0)=0
        xi = x / self.l
        phi = 1.5 * xi**2 - 0.5 * xi**3  # Static deflection shape

        # Second derivative of mode shape (curvature)
        # d^2(phi)/dx^2 = d^2/dx^2 [1.5*(x/L)^2 - 0.5*(x/L)^3]
        #               = 3/L^2 - 3x/L^3
        phi_ddot = 3 / self.l**2 - 3 * x / self.l**3

        # Potential energy: U = 0.5 * integral(EI * phi''^2 dx)
        U = 0.5 * np.trapezoid(EI * phi_ddot**2, x)

        # Kinetic energy coefficient: T = 0.5 * omega^2 * integral(dm_dx * phi^2 dx)
        T_coeff = np.trapezoid(dm_dx * phi**2, x)

        # omega^2 = 2*U / T_coeff
        # Fallback to simple formula if T_coeff is zero or negative
        omega = math.sqrt(2 * U / T_coeff) if T_coeff > 0 else math.sqrt(self.stiffness() / self.effective_mass())

        return omega

    # Resonant frequency for undamped vibration (first mode)
    # Units: cycles/sec
    def omega_vacuum_hz(self):
        omega_vacuum_hz = self.omega_vacuum() / (2 * math.pi)
        return omega_vacuum_hz

    # Calculate the damped natural frequency and Q (via Brumley/Sader)
    def omega_damped_and_Q(self):
        # If we're in vacuum, just return the vacuum frequency
        if self.fluid == "vacuum":
            omega_damped = self.omega_vacuum()
            Q = Cantilever.maxQ
            return omega_damped, Q

        # Inner function for solving the transcendental eqn to find
        # omega_damped. We're searching for a function minimum, so return
        # the residual squared (continuous and smooth)
        def find_natural_frequency(omega_damped):
            hydro = self.hydrodynamic_function(omega_damped, rho_f, eta_f)
            residual = (
                omega_damped
                - omega_vacuum * (1 + math.pi * rho_f * self.w / (4 * self.rho_si * self.t) * np.real(hydro)) ** -0.5
            )
            residual_squared = residual**2
            return residual_squared

        # Lookup fluid properties once, then calculate omega_damped and Q
        rho_f, eta_f = self.lookupFluidProperties()
        omega_vacuum = self.omega_vacuum()

        omega_damped = optimize.fminbound(find_natural_frequency, 10, omega_vacuum, xtol=1e-12)

        hydro = self.hydrodynamic_function(omega_damped, rho_f, eta_f)
        Q = (4 * self.rho_si * self.t / (math.pi * rho_f * self.w) + np.real(hydro)) / np.imag(hydro)

        # Catch cases where Q is undefined, too large or too small
        if Cantilever.minQ > Q or np.isnan(Q):
            Q = Cantilever.minQ
        elif Cantilever.maxQ < Q:
            Q = Cantilever.maxQ
        return omega_damped, Q

    def omega_damped_hz_and_Q(self):
        omega_damped, Q = self.omega_damped_and_Q()
        omega_damped_hz = omega_damped / (2 * math.pi)
        return omega_damped_hz, Q

    # Calculate the Reynold's number - note that 'a = w/2' in the Brumley paper
    def reynolds(self, omega, rho_f, eta_f):
        reynolds = (rho_f * omega * (self.w / 2) ** 2) / eta_f
        return reynolds

    # Calculate hydrodynamic function from the lookup table
    def hydrodynamic_function(self, omega, rho_f, eta_f):
        # A = aspect ratio, Beta = Reynold's number
        A = self.t / self.w
        Beta = self.reynolds(omega, rho_f, eta_f)
        log_Beta = math.log10(Beta)

        # If the Reynolds number gets too small, the calculation gets upset
        log_Beta = max(log_Beta, min(self.Beta_lookup))

        # Interpolate from lookup tables using RectBivariateSpline
        # gamma arrays are shaped (len(Beta_lookup), len(A_lookup))
        interp_real = interpolate.RectBivariateSpline(
            Cantilever.Beta_lookup, Cantilever.A_lookup, Cantilever.gamma_lookup_real
        )
        interp_imag = interpolate.RectBivariateSpline(
            Cantilever.Beta_lookup, Cantilever.A_lookup, Cantilever.gamma_lookup_imag
        )
        gamma_real = interp_real(log_Beta, A)[0, 0]
        gamma_imag = interp_imag(log_Beta, A)[0, 0]
        hydro = complex(gamma_real, gamma_imag)
        return hydro

    # # ========= Optimization ==========
    #
    # # Calculate force resolution (pN)
    # def optimize_force_resolution(self, x0):
    #     self = self.cantilever_from_state(x0)
    #     force_resolution = self.force_resolution()*1e12
    #     return force_resolution
    #
    # # Calculate displacement resolution (nm)
    # def optimize_displacement_resolution(self, x0):
    #     self = self.cantilever_from_state(x0)
    #     displacement_resolution = self.displacement_resolution()*1e9
    #     return displacement_resolution
    #
    # # Calculate the force noise density on resonance (pN/rtHz)
    # def optimize_resonant_force_noise_density(self, x0):
    #     self = self.cantilever_from_state(x0)
    #     force_noise_density = self.resonant_force_noise_density()*1e12
    #     return force_noise_density
    #
    # # Calculate the surface stress resolution (Pa)
    # def optimize_surface_stress_resolution(self, x0):
    #     self = self.cantilever_from_state(x0)
    #     ss_resolution = self.surface_stress_resolution()*1e6
    #     return ss_resolution
    #
    # # Used by optimization to bring all state varibles to O(10)
    # # If we didn't do that, we'd have O(1e-9) and O(1e20) variables
    # def optimization_scaling(self):
    #     l_scale = 1e5
    #     w_scale = 1e7
    #     t_scale = 1e8
    #     l_pr_ratio_scale = 1e2
    #     v_bridge_scale = 1e1
    #     scaling = np.array((l_scale, w_scale, t_scale, l_pr_ratio_scale,
    #                         v_bridge_scale, self.doping_optimization_scaling()))
    #
    #     # Actuator specific code
    #     if self.cantilever_type == 'thermal':
    #         l_a_scale = 1e6
    #         w_a_scale = 1e6
    #         t_a_scale = 1e9
    #         v_actuator_scale = 1
    #         R_heater_scale = 1e-3
    #         scaling = np.concatenate((scaling, (l_a_scale, w_a_scale, t_a_scale, v_actuator_scale, R_heater_scale)))
    #     elif self.cantilever_type == 'piezoelectric':
    #         l_a_scale = 1e6
    #         w_a_scale = 1e6
    #         t_a_scale = 1e9
    #         v_actuator_scale = 1
    #         scaling = np.concatenate((scaling, (l_a_scale, w_a_scale, t_a_scale, v_actuator_scale)))
    #     return scaling
    #
    # # Update the changed optimization parameters
    # # All optimization takes place for the same object (i.e. we update 'self')
    # # so that things like 'fluid' are maintained
    # def cantilever_from_state(self, x0):
    #     scaling = self.optimization_scaling()
    #     x0 = x0/scaling
    #
    #     self.l = x0[0]
    #     self.w = x0[1]
    #     self.t = x0[2]
    #     self.l_pr_ratio = x0[3]
    #     self.v_bridge = x0[4]
    #     self = self.doping_cantilever_from_state(x0)
    #
    #     # Actuator specific code
    #     if self.cantilever_type == 'thermal':
    #         self.l_a = x0[7]
    #         self.w_a = x0[8]
    #         self.t_a = x0[9]
    #         self.v_actuator = x0[10]
    #         self.R_heater = x0[11]
    #     elif self.cantilever_type == 'piezoelectric':
    #         self.l_a = x0[7]
    #         self.w_a = x0[8]
    #         self.t_a = x0[9]
    #         self.v_actuator = x0[10]
    #     return self
    #
    # # Return state vector for the current state
    # def current_state(self):
    #     if self.cantilever_type == 'thermal':
    #         x = np.zeros((11,1))
    #     elif self.cantilever_type == 'piezoelectric':
    #         x = np.zeros((10,1))
    #     else:
    #         x = np.zeros((6,1))
    #
    #     x[0] = self.l
    #     x[1] = self.w
    #     x[2] = self.t
    #     x[3] = self.l_pr_ratio
    #     x[4] = self.v_bridge
    #     x = np.concatenate((x, self.doping_current_state()))
    #
    #     # Actuator specific code
    #     if self.cantilever_type == 'thermal':
    #         x[7] = self.l_a
    #         x[8] = self.w_a
    #         x[9] = self.t_a
    #         x[10] = self.v_actuator
    #         x[11] = self.R_heater
    #     elif self.cantilever_type == 'piezoelectric':
    #         x[7] = self.l_a
    #         x[8] = self.w_a
    #         x[9] = self.t_a
    #         x[10] = self.v_actuator
    #     return x
    #
    # # Set the minimum and maximum bounds for the cantilever state variables.
    # # Bounds are written in terms of the initialization variables.
    # # Secondary constraints (e.g. power dissipation, resonant frequency)
    # # are applied in optimization_constraints()
    # def optimization_bounds(self, parameter_constraints):
    #     min_l = 10e-6
    #     max_l = 3e-3
    #
    #     min_w = 2e-6
    #     max_w = 100e-6
    #
    #     min_t = 1e-6
    #     max_t = 100e-6
    #
    #     min_l_pr_ratio = 0.01
    #     max_l_pr_ratio = 0.99
    #
    #     min_v_bridge = 0.1
    #     max_v_bridge = 10
    #
    #     doping_lb, doping_ub = self.doping_optimization_bounds(parameter_constraints)
    #
    #     # Actuator specific code
    #     # TODO
    #     actuator_lb = []
    #     actuator_ub = []
    #     if self.cantilever_type == 'step' or 'none':
    #         1
    #         # TODO Use a dictionary to allow for arbitrary parameters to be overriden
    #         # Return lb/ub with dictionaries to fix this (or assemble arrays from dictionaries)
    #         # Override the default values if any were provided
    #         # constraints = {{'min_l', 'max_v_bridge'}, {5e-6, 10}}
    #         #if ~isempty(parameter_constraints)
    #         #    keys = parameter_constraints{1}
    #         #    values = parameter_constraints{2}
    #         #    for ii in range(1, keys.size + 1):
    #         #        eval([keys{ii} '=' num2str(values{ii}) ''])
    #         #    end
    #         #end
    #     elif self.cantilever_type == 'thermal':
    #         min_l_a = 5e-6
    #         max_l_a = 200e-6
    #
    #         min_w_a = 2e-6
    #         max_w_a = 50e-6
    #
    #         min_t_a = 200e-9
    #         max_t_a = 3e-6
    #
    #         min_v_actuator = .1
    #         max_v_actuator = 10
    #
    #         min_R_heater = 200
    #         max_R_heater = 5e3
    #
    #         # TODO Fix me
    #         # Override the default values if any were provided
    #         # constraints = {{'min_l', 'max_v_bridge'}, {5e-6, 10}}
    #         #if ~isempty(parameter_constraints)
    #         #    keys = parameter_constraints{1}
    #         #    values = parameter_constraints{2}
    #         #    for ii in range(1, keys.size + 1):
    #         #        eval([keys{ii} '=' num2str(values{ii}) ''])
    #         #    end
    #         #end
    #
    #         actuator_lb = [min_l_a, min_w_a, min_t_a, min_v_actuator, min_R_heater]
    #         actuator_ub = [max_l_a, max_w_a, max_t_a, max_v_actuator, max_R_heater]
    #     elif self.cantilever_type == 'piezoelectric':
    #         min_l_a = 5e-6
    #         max_l_a = 200e-6
    #
    #         min_w_a = 2e-6
    #         max_w_a = 30e-6
    #
    #         min_t_a = 200e-9
    #         max_t_a = 3e-6
    #
    #         min_v_actuator = .1
    #         max_v_actuator = 10
    #
    #         # TODO Fix me
    #         # Override the default values if any were provided
    #         # constraints is a set of key value pairs, e.g.
    #         # constraints = {{'min_l', 'max_v_bridge'}, {5e-6, 10}}
    #         # if ~isempty(parameter_constraints)
    #         #     keys = parameter_constraints{1}
    #         #     values = parameter_constraints{2}
    #         #     for ii in range(1, keys.size + 1):
    #         #         eval([keys{ii} '=' num2str(values{ii}) ''])
    #         #     end
    #         # end
    #
    #         actuator_lb = [min_l_a, min_w_a, min_t_a, min_v_actuator]
    #         actuator_ub = [max_l_a, max_w_a, max_t_a, max_v_actuator]
    #
    #     lb = [min_l, min_w, min_t, min_l_pr_ratio, min_v_bridge, doping_lb, actuator_lb]
    #     ub = [max_l, max_w, max_t, max_l_pr_ratio, max_v_bridge, doping_ub, actuator_ub]
    #
    #     return lb, ub
    #
    # # Generate a random cantilever design (x0) within our param bounds
    # def initial_conditions_random(self, parameter_constraints):
    #     lb, ub = self.optimization_bounds(parameter_constraints)
    #
    #     l_min = lb[0]
    #     l_max = ub[0]
    #     w_min = lb[1]
    #     w_max = ub[1]
    #     t_min = lb[2]
    #     t_max = ub[2]
    #     l_pr_ratio_min = lb[3]
    #     l_pr_ratio_max = ub[3]
    #     V_b_min = lb[4]
    #     V_b_max = ub[4]
    #
    #     # Generate the random values
    #     l_random = l_min + np.random.rand()*(l_max - l_min)
    #     w_random = w_min + np.random.rand()*(w_max - w_min)
    #     t_random = t_min + np.random.rand()*(t_max - t_min)
    #     l_pr_ratio_random = l_pr_ratio_min + np.random.rand()*(l_pr_ratio_max - l_pr_ratio_min)
    #     v_bridge_random = V_b_min + np.random.rand()*(V_b_max - V_b_min)
    #
    #     x0_doping = self.doping_initial_conditions_random(parameter_constraints)
    #
    #     # Actuator specific code
    #     x0_actuator = []
    #     if self.cantilever_type == 'thermal':
    #         l_a_min = lb[7]
    #         l_a_max = ub[7]
    #         w_a_min = lb[8]
    #         w_a_max = ub[8]
    #         t_a_min = lb[9]
    #         t_a_max = ub[9]
    #         v_actuator_min = lb[10]
    #         v_actuator_max = ub[10]
    #         R_heater_min = lb[11]
    #         R_heater_max = ub[11]
    #
    #         l_a_random = l_a_min + np.random.rand()*(l_a_max - l_a_min)
    #         w_a_random = w_a_min + np.random.rand()*(w_a_max - w_a_min)
    #         t_a_random = t_a_min + np.random.rand()*(t_a_max - t_a_min)
    #         v_actuator_random = v_actuator_min + np.random.rand()*(v_actuator_max - v_actuator_min)
    #         R_heater_random = R_heater_min + np.random.rand()*(R_heater_max - R_heater_min)
    #         x0_actuator = [l_a_random, w_a_random, t_a_random, v_actuator_random, R_heater_random]
    #     elif self.cantilever_type == 'piezoelectric':
    #         l_a_min = lb[7]
    #         l_a_max = ub[7]
    #         w_a_min = lb[8]
    #         w_a_max = ub[8]
    #         t_a_min = lb[9]
    #         t_a_max = ub[9]
    #         v_actuator_min = lb[10]
    #         v_actuator_max = ub[10]
    #
    #         l_a_random = l_a_min + np.random.rand()*(l_a_max - l_a_min)
    #         w_a_random = w_a_min + np.random.rand()*(w_a_max - w_a_min)
    #         t_a_random = t_a_min + np.random.rand()*(t_a_max - t_a_min)
    #         v_actuator_random = v_actuator_min + np.random.rand()*(v_actuator_max - v_actuator_min)
    #         x0_actuator = [l_a_random, w_a_random, t_a_random, v_actuator_random]
    #
    #     x0 = [l_random, w_random, t_random, l_pr_ratio_random, v_bridge_random, x0_doping, x0_actuator]
    #     return x0
    #
    # # Nonlinear optimization constraints.
    # # For a feasible design, all constraints are negative.
    # def optimization_constraints(self, x0, nonlinear_constraints):
    #     c_new = self.cantilever_from_state(x0)
    #
    #     # Default aspect ratios that can be overriden
    #     min_w_t_ratio = 2
    #     min_l_w_ratio = 2
    #     min_pr_l_w_ratio = 2
    #     min_pr_l = 2e-6
    #
    #     # TODO Fix passing key/values for optimization
    #     # Read out the constraints as key-value pairs
    #     # , e.g. {{'omega_min_hz', 'min_k'}, {1000, 10}}
    #     #if ~isempty(nonlinear_constraints)
    #     #     keys = nonlinear_constraints{1}
    #     #     values = nonlinear_constraints{2}
    #     #     for ii in range(1, keys.size + 1)
    #     #         eval([keys{ii} '=' num2str(values{ii}) ''])
    #     #     end
    #     # end
    #
    #     # We start with this single element vector and then append
    #     # any additional constraints that the user has provided.
    #     # If a C(ii) > 0 then constraint ii is invalid
    #     C = []
    #
    #     # Force resolution must always be positive
    #     # Use a linear function so that it's smooth/continuous
    #     resolution = c_new.force_resolution()
    #     C.append(-resolution*1e18)
    #
    #     # Resonant frequency
    #     if nonlinear_constraints.get('omega_min_hz'):
    #         omega_min_hz = nonlinear_constraints.get('omega_min_hz')
    #         if self.fluid == 'vacuum':
    #             freq_constraint = omega_min_hz - c_new.omega_vacuum_hz()
    #         else:
    #             [omega_damped_hz, tmp] = c_new.omega_damped_hz_and_Q()
    #             freq_constraint = omega_min_hz - omega_damped_hz
    #         C.append(freq_constraint)
    #
    #     # Power dissipation
    #     if nonlinear_constraints.get('max_power'):
    #         max_power = nonlinear_constraints.get('max_power')
    #         power_constraint = c_new.power_dissipation() - max_power
    #         C.append(power_constraint)
    #
    #     # Temp constraints -- approximate lumped model
    #     if nonlinear_constraints.get('tip_temp') or nonlinear_constraints.get('max_temp'):
    #         [TMax, TTip] = c_new.approxTempRise()
    #
    #     if nonlinear_constraints.get('tip_temp'):
    #         tip_temp = nonlinear_constraints.get('tip_temp')
    #         C.append(TTip - tip_temp)
    #
    #     if nonlinear_constraints.get('max_temp'):
    #         max_temp = nonlinear_constraints.get('max_temp')
    #         C.append(TMax - max_temp)
    #
    #     # Temp constraints -- 1D finite differences
    #     # Slower but good for refining designs
    #     if nonlinear_constraints.get('tip_temp_exact') or nonlinear_constraints.get('max_temp_exact'):
    #         [TMax, TTip] = c_new.calculateMaxAndTipTemp()
    #
    #     if nonlinear_constraints.get('tip_temp_exact'):
    #         tip_temp_exact = nonlinear_constraints.get('tip_temp_exact')
    #         C.append(TTip - tip_temp_exact)
    #
    #     if nonlinear_constraints.get('max_temp_exact'):
    #         max_temp_exact = nonlinear_constraints.get('max_temp_exact')
    #         C.append(TMax - max_temp_exact)
    #
    #     # Min and maximum cantilever stiffness
    #     if nonlinear_constraints.get('min_k'):
    #         min_k = nonlinear_constraints.get('min_k')
    #         C.append(min_k - c_new.stiffness())
    #
    #     if nonlinear_constraints.get('max_k'):
    #         max_k = nonlinear_constraints.get('max_k')
    #         C.append(c_new.stiffness() - max_k)
    #
    #     # Aspect ratio constraints. Default ratios can be changed.
    #     C.append(min_l_w_ratio - c_new.l/c_new.w)
    #     C.append(min_w_t_ratio - c_new.w/c_new.t)
    #     C.append(min_pr_l_w_ratio - c_new.l_pr()/c_new.w_pr())
    #     C.append(min_pr_l - c_new.l_pr())
    #
    #     # Now for equality based constraints
    #     Ceq = []
    #
    #     # Fix the stiffness
    #     if nonlinear_constraints.get('fixed_k'):
    #         fixed_k = nonlinear_constraints.get('fixed_k')
    #         Ceq.append(c_new.stiffness() - fixed_k)
    #
    #     if nonlinear_constraints.get('fixed_v_bridge'):
    #         fixed_v_bridge = nonlinear_constraints.get('fixed_v_bridge')
    #         Ceq.append(c_new.v_bridge - fixed_v_bridge)
    #
    #     # Fix the resonant frequency
    #     if nonlinear_constraints.get('fixed_f0'):
    #         fixed_f0 = nonlinear_constraints.get('fixed_f0')
    #         if self.fluid == 'vacuum':
    #             fixed_f0_constraint = fixed_f0 - c_new.omega_vacuum_hz()
    #         else:
    #             [omega_damped_hz, tmp] = c_new.omega_damped_hz_and_Q()
    #             fixed_f0_constraint = fixed_f0 - omega_damped_hz
    #         Ceq.append(fixed_f0_constraint)
    #
    #     # Useful lines for debugging constraint failures (e.g. you made l_max too small for it to hit the desired f_0)
    #     # print C, Ceq,
    #     # print 'Active Index: %d -- Value: %g' % find(C==max(C)), max(C)
    #
    #     return C, Ceq
    #
    # # The optimization isn't convex so isn't guaranteed to converge.
    # # On average it converges 95-99# of the time.
    # # For this reason, we generally optimize from random initial guesses
    # # and keep going until several have converged to the same value
    # def optimized_cantilever = optimize_performance(self, parameter_constraints, nonlinear_constraints, goal):
    #
    #     percent_match = 0.01 # 1 percent match required between results
    #     randomize_starting_conditions = 1
    #     converged = 0
    #     ii = 1
    #     resolution = []
    #     c = []
    #     while ~converged:
    #         # Optimize another cantilever
    #         c_current, exitflag = self.optimize_performance_once(parameter_constraints,
    #                                                 nonlinear_constraints, goal, randomize_starting_conditions)
    #         c.append(c_current)
    #
    #         # If the optimization terminated abnormally, skip to the next iteration
    #         #if ~(exitflag == 1 or exitflag == 2):
    #         #    continue
    #
    #         # Record the resolution for the latest cantilever
    #         goal = []
    #         # if goal == Cantilever.goalForceResolution:
    #         #     resolution[ii] = c_current.force_resolution()
    #         # elif goal == Cantilever.goalDisplacementResolution:
    #         #     resolution[ii] = c_current.displacement_resolution()
    #         # elif goal == Cantilever.goalForceNoiseDensity:
    #         #     resolution[ii] = c_current.resonant_force_noise_density()
    #         # elif goal == Cantilever.goalSurfaceStress:
    #         #     resolution[ii] = c_current.surface_stress_resolution()
    #
    #         # If we have more than one result, consider stopping
    #         if ii > 1:
    #             # Sort from smallest to largest, check if the two smallest values agree
    #             [temp_resolution, sortIndex] = np.sort(resolution)
    #             print 'Resolutions so far: %s\n' % temp_resolution
    #             resultsAgree = abs(1 - temp_resolution(1)/temp_resolution(2)) < percent_match
    #
    #             # If the results agree, then stop the loop. Otherwise, continue
    #             if resultsAgree:
    #                 print 'CONVERGED. Two best values: %s' % temp_resolution
    #                 # optimized_cantilever = c{sortIndex(1)}
    #                 converged = 1
    #             else:
    #                 print 'NOT CONVERGED. Two best values: %s' % temp_resolution
    #
    #         if ii > Cantilever.numOptimizationIterations:
    #             temp_resolution, sortIndex = np.sort(resolution)
    #             # optimized_cantilever = c{sortIndex(1)}
    #             converged = 1
    #
    #         ii = ii + 1 # Increment the loop counter
    #     return optimized_cantilever
    #
    # # Optimize, but don't randomize starting point
    # def optimized_cantilever = optimize_performance_from_current(self,
    #                                  parameter_constraints, nonlinear_constraints, goal):
    #     randomize_starting_conditions = 0
    #     [optimized_cantilever, tmp] = self.optimize_performance_once(parameter_constraints,
    #         nonlinear_constraints, goal, randomize_starting_conditions)
    #
    # def optimized_cantilever, exitflag = optimize_performance_once(self, parameter_constraints,
    #                                                     nonlinear_constraints, goal, randomize_starting_conditions):
    #
    #     scaling = self.optimization_scaling()
    #     self.check_valid_cantilever()
    #
    #     # If random_flag = 1, start from random conditions
    #     # Else, start from the current cantilever state vector
    #     if randomize_starting_conditions == 1:
    #         problem.x0 = scaling*self.initial_conditions_random(parameter_constraints)
    #     else:
    #         problem.x0 = scaling*self.current_state()
    #
    #     #if goal == cantilever.goalForceResolution:
    #     #    problem.objective = self.optimize_force_resolution
    #     #elif goal == cantilever.goalDisplacementResolution:
    #     #    problem.objective = self.optimize_displacement_resolution
    #     #elif goal == cantilever.goalForceNoiseDensity:
    #     #    problem.objective = self.optimize_resonant_force_noise_density
    #     #elif goal == cantilever.goalSurfaceStress:
    #     #    problem.objective = self.optimize_surface_stress_resolution
    #
    #     # [lb ub] = self.optimization_bounds(parameter_constraints)
    #     # problem.lb = scaling*lb
    #     # problem.ub = scaling*ub
    #     #
    #     # problem.options.TolFun = 1e-9
    #     # problem.options.TolCon = 1e-9
    #     # problem.options.TolX = 1e-9
    #     #
    #     # problem.options.MaxFunEvals = 2000
    #     # problem.options.MaxIter = 2000
    #     # problem.options.Display = 'iter-detailed' # iter, final, none, off
    #     # problem.options.UseParallel = 'always' # For multicore processors
    #     # problem.options.TypicalX = np.ones((1, scaling.size))
    #     # problem.options.InitBarrierParam = 10
    #
    #     # interior-point does not do well if the design is already
    #     # close to the optimal state (i.e. warm start). sqp is
    #     # supposed to do better in this condition. try experimenting
    #     # with both, particular for implanted devices.
    #     # problem.options.Algorithm = 'sqp' # sqp, interior-point
    #     # problem.solver = 'fmincon'
    #     #
    #     # problem.nonlcon = @(x) self.optimization_constraints(x, ...
    #     #     nonlinear_constraints)
    #     #
    #     # # These errors come up occasionally for ion implanted devices
    #     # warning off MATLAB:singularMatrix
    #     # warning off MATLAB:nearlySingularMatrix
    #     #
    #     # [x, tmp, exitflag] = fmincon(problem)
    #     # optimized_cantilever = self.cantilever_from_state(x)
