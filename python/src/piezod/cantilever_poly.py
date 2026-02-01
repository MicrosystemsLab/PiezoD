"""Polycrystalline thin-film cantilever with multi-layer beam mechanics.

This module models polycrystalline thin-film cantilevers with 3 layers (top, mid, bottom)
using different materials (Si, poly-Si, Ti, Al, oxide). It calculates beam mechanics,
sheet resistance for poly/Ti/Al piezoresistors, and noise characteristics.

This is a standalone class that does NOT inherit from the base Cantilever class,
as it uses fundamentally different multi-layer beam mechanics.
"""

from enum import Enum

import numpy as np
from numpy.typing import NDArray


class Material(Enum):
    """Material types for cantilever layers."""

    SI = "si"
    POLY = "poly"
    TI = "ti"
    AL = "al"
    OXIDE = "oxide"


class CantileverPoly:
    """Polycrystalline thin-film cantilever with multi-layer beam mechanics.

    Models a cantilever with three layers (top, mid, bottom) that can be made of
    different materials. Supports piezoresistors made of polycrystalline silicon,
    titanium, or aluminum.

    Attributes:
        l: Cantilever length (m)
        w: Cantilever width (m)
        t_top: Top layer thickness (m)
        t_mid: Middle layer thickness (m)
        t_bot: Bottom layer thickness (m)
        matl_top: Top layer material
        matl_mid: Middle layer material
        matl_bot: Bottom layer material
        l_pr_ratio: Piezoresistor length to cantilever length ratio
        v_bridge: Wheatstone bridge voltage (V)
        dopant_concentration: Dopant concentration for poly-Si (cm^-3)
    """

    # Physical constants
    k_b = 1.38e-23  # J/K
    k_b_eV = 8.617343e-5  # eV/K
    q = 1.60218e-19  # Coulombs
    T = 300  # Temperature (K)

    # Material densities (kg/m^3)
    rho_si = 2330
    rho_poly = 2330
    rho_oxide = 2634
    rho_al = 2700
    rho_ti = 4506

    # Elastic moduli (Pa)
    E_si = 169e9
    E_poly = 150e9
    E_oxide = 70e9
    E_al = 70e9
    E_ti = 85e9

    # Hooge parameter for polycrystalline materials (higher than single-crystal)
    alpha = 1e-3

    # Number of frequency points for integration
    numFrequencyPoints = 1000

    def __init__(
        self,
        freq_min: float = 1.0,
        freq_max: float = 1e3,
        l: float = 100e-6,
        w: float = 10e-6,
        t_top: float = 100e-9,
        t_mid: float = 500e-9,
        t_bot: float = 100e-9,
        matl_top: Material = Material.POLY,
        matl_mid: Material = Material.SI,
        matl_bot: Material = Material.POLY,
        l_pr_ratio: float = 0.5,
        v_bridge: float = 1.0,
        dopant_concentration: float = 1e19,
        number_of_piezoresistors: int = 2,
        number_of_piezoresistors_on_cantilever: int = 2,
    ) -> None:
        """Initialize a polycrystalline cantilever.

        Args:
            freq_min: Minimum frequency for analysis (Hz)
            freq_max: Maximum frequency for analysis (Hz)
            l: Cantilever length (m)
            w: Cantilever width (m)
            t_top: Top layer thickness (m)
            t_mid: Middle layer thickness (m)
            t_bot: Bottom layer thickness (m)
            matl_top: Top layer material
            matl_mid: Middle layer material
            matl_bot: Bottom layer material
            l_pr_ratio: Piezoresistor length / cantilever length ratio
            v_bridge: Bridge voltage (V)
            dopant_concentration: Dopant concentration for poly-Si (cm^-3)
            number_of_piezoresistors: Number of piezoresistors in bridge
            number_of_piezoresistors_on_cantilever: Number of piezoresistors on cantilever
        """
        self.freq_min = freq_min
        self.freq_max = freq_max
        self.l = l
        self.w = w
        self.t_top = t_top
        self.t_mid = t_mid
        self.t_bot = t_bot
        self.matl_top = matl_top
        self.matl_mid = matl_mid
        self.matl_bot = matl_bot
        self.l_pr_ratio = l_pr_ratio
        self.v_bridge = v_bridge
        self.dopant_concentration = dopant_concentration
        self.number_of_piezoresistors = number_of_piezoresistors
        self.number_of_piezoresistors_on_cantilever = number_of_piezoresistors_on_cantilever

        # Enforce symmetry for double piezoresistor configuration
        if self.number_of_piezoresistors_on_cantilever == 2:
            self.t_bot = self.t_top

    @property
    def l_pr(self) -> float:
        """Piezoresistor length (m)."""
        return self.l * self.l_pr_ratio

    @property
    def w_pr(self) -> float:
        """Piezoresistor width (m)."""
        return self.w / 2

    @property
    def t_total(self) -> float:
        """Total cantilever thickness (m)."""
        return self.t_top + self.t_mid + self.t_bot

    def _get_elastic_modulus(self, material: Material) -> float:
        """Get elastic modulus for a material."""
        modulus_map = {
            Material.SI: self.E_si,
            Material.POLY: self.E_poly,
            Material.TI: self.E_ti,
            Material.AL: self.E_al,
            Material.OXIDE: self.E_oxide,
        }
        return modulus_map[material]

    def _get_density(self, material: Material) -> float:
        """Get density for a material."""
        density_map = {
            Material.SI: self.rho_si,
            Material.POLY: self.rho_poly,
            Material.TI: self.rho_ti,
            Material.AL: self.rho_al,
            Material.OXIDE: self.rho_oxide,
        }
        return density_map[material]

    def neutral_axis(self) -> float:
        """Calculate the neutral axis position from the bottom surface.

        Returns:
            Neutral axis position (m)
        """
        # Centroid positions of each layer from bottom
        z = np.array(
            [
                self.t_bot / 2,
                self.t_bot + self.t_mid / 2,
                self.t_bot + self.t_mid + self.t_top / 2,
            ]
        )

        # Elastic moduli
        E = np.array(
            [
                self._get_elastic_modulus(self.matl_bot),
                self._get_elastic_modulus(self.matl_mid),
                self._get_elastic_modulus(self.matl_top),
            ]
        )

        # Cross-sectional areas
        A = self.w * np.array([self.t_bot, self.t_mid, self.t_top])

        # Neutral axis from weighted average
        return float(np.sum(z * E * A) / np.sum(E * A))

    def normalized_curvature(self) -> float:
        """Calculate normalized curvature (1/EI_effective).

        Returns:
            Normalized curvature (1/N-m^2)
        """
        Zm = self.neutral_axis()

        # Distance from neutral axis to each layer centroid
        z = np.array(
            [
                self.t_bot / 2,
                self.t_bot + self.t_mid / 2,
                self.t_bot + self.t_mid + self.t_top / 2,
            ]
        )
        Z = z - Zm

        # Elastic moduli
        E = np.array(
            [
                self._get_elastic_modulus(self.matl_bot),
                self._get_elastic_modulus(self.matl_mid),
                self._get_elastic_modulus(self.matl_top),
            ]
        )

        # Cross-sectional areas and second moments of area
        A = self.w * np.array([self.t_bot, self.t_mid, self.t_top])
        I = self.w * np.array(
            [
                self.t_bot**3 / 12,
                self.t_mid**3 / 12,
                self.t_top**3 / 12,
            ]
        )

        # EI_effective using parallel axis theorem
        return float(1.0 / np.sum(E * (I + A * Z**2)))

    def stiffness(self) -> float:
        """Calculate spring constant of the cantilever.

        Returns:
            Stiffness (N/m)
        """
        EI_effective = 1.0 / self.normalized_curvature()
        return 3 * EI_effective / self.l**3

    def resonant_frequency(self) -> float:
        """Calculate vacuum resonant frequency.

        Returns:
            Resonant frequency (Hz)
        """
        densities = np.array(
            [
                self._get_density(self.matl_bot),
                self._get_density(self.matl_mid),
                self._get_density(self.matl_top),
            ]
        )
        thicknesses = np.array([self.t_bot, self.t_mid, self.t_top])

        # Effective mass (0.243 is the first mode coefficient)
        m_eff = 0.243 * self.l * self.w * np.sum(thicknesses * densities)
        return np.sqrt(self.stiffness() / m_eff) / (2 * np.pi)

    def mobility(self, dopant_concentration: float) -> float:
        """Calculate carrier mobility using Masetti model with poly correction.

        Based on Masetti et al. 1983, with 50% reduction for polycrystalline silicon.

        Args:
            dopant_concentration: Dopant concentration (cm^-3)

        Returns:
            Carrier mobility (cm^2/V-s)
        """
        n = dopant_concentration

        # Phosphorus parameters (Masetti 1983)
        mu_0 = 68.5
        mu_max = 1414
        mu_1 = 56.1
        C_r = 9.2e16
        C_s = 3.41e20
        alpha = 0.711
        beta = 1.98

        mobility = mu_0 + (mu_max - mu_0) / (1 + (n / C_r) ** alpha) - mu_1 / (1 + (C_s / n) ** beta)

        # 50% reduction for polycrystalline silicon
        return mobility * 0.5

    def sheet_resistance(self) -> float:
        """Calculate sheet resistance of the piezoresistor layer.

        Returns:
            Sheet resistance (ohms/square)
        """
        if self.matl_top == Material.POLY:
            mobility = self.mobility(self.dopant_concentration)
            resistivity = 1.0 / (self.q * mobility * self.dopant_concentration)
        elif self.matl_top == Material.TI:
            resistivity = 42e-6  # ohm-cm
        elif self.matl_top == Material.AL:
            resistivity = 2.7e-6  # ohm-cm
        else:
            # Default to poly-Si calculation
            mobility = self.mobility(self.dopant_concentration)
            resistivity = 1.0 / (self.q * mobility * self.dopant_concentration)

        return resistivity / (self.t_top * 1e2)  # Convert thickness to cm

    def resistor_length(self) -> float:
        """Calculate total resistor length.

        Returns:
            Resistor length (m)
        """
        return 2 * self.l_pr

    def gamma(self) -> float:
        """Ratio of piezoresistor resistance to total resistance.

        Accounts for metal interconnect resistance (~10% of total).

        Returns:
            Resistance ratio (dimensionless)
        """
        return 0.9

    def resistance(self) -> float:
        """Calculate total resistance of the piezoresistor.

        Returns:
            Resistance (ohms)
        """
        number_of_squares = self.resistor_length() / self.w_pr
        return number_of_squares * self.sheet_resistance() / self.gamma()

    def number_of_carriers(self) -> float:
        """Calculate number of charge carriers in the piezoresistor.

        Returns:
            Number of carriers (dimensionless)
        """
        if self.matl_top == Material.POLY:
            Nz = self.t_top * self.dopant_concentration * 1e6  # 1/m^2
        else:
            # For Ti and Al, use approximate carrier density
            Nz = self.t_top * 2e21 * 1e6  # 1/m^2

        return Nz * self.resistor_length() * self.w_pr

    def hooge_PSD(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate 1/f noise power spectral density.

        Args:
            freq: Frequency array (Hz)

        Returns:
            Voltage PSD (V^2/Hz)
        """
        return self.alpha * self.v_bridge**2 * self.number_of_piezoresistors / (4 * self.number_of_carriers() * freq)

    def hooge_integrated(self) -> float:
        """Calculate integrated 1/f noise.

        Returns:
            Integrated noise voltage (V)
        """
        freq = np.logspace(np.log10(self.freq_min), np.log10(self.freq_max), self.numFrequencyPoints)
        return float(np.sqrt(np.trapezoid(self.hooge_PSD(freq), freq)))

    def johnson_PSD(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate Johnson noise power spectral density.

        Args:
            freq: Frequency array (Hz)

        Returns:
            Voltage PSD (V^2/Hz)
        """
        return 4 * self.k_b * self.T * self.resistance() * np.ones_like(freq)

    def johnson_integrated(self) -> float:
        """Calculate integrated Johnson noise.

        Returns:
            Integrated noise voltage (V)
        """
        freq = np.logspace(np.log10(self.freq_min), np.log10(self.freq_max), self.numFrequencyPoints)
        return float(np.sqrt(np.trapezoid(self.johnson_PSD(freq), freq)))

    def thermo_PSD(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate thermomechanical noise power spectral density.

        Args:
            freq: Frequency array (Hz)

        Returns:
            Voltage PSD (V^2/Hz)
        """
        Q_M = 100  # Quality factor
        return (
            self.force_sensitivity() ** 2
            * 2
            * self.stiffness()
            * self.k_b
            * self.T
            / (np.pi * self.resonant_frequency() * Q_M)
            * np.ones_like(freq)
        )

    def thermo_integrated(self) -> float:
        """Calculate integrated thermomechanical noise.

        Returns:
            Integrated noise voltage (V)
        """
        freq = np.logspace(np.log10(self.freq_min), np.log10(self.freq_max), self.numFrequencyPoints)
        return float(np.sqrt(np.trapezoid(self.thermo_PSD(freq), freq)))

    def amplifier_PSD(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate amplifier noise power spectral density.

        Uses INA103 amplifier parameters.

        Args:
            freq: Frequency array (Hz)

        Returns:
            Voltage PSD (V^2/Hz)
        """
        # INA103 parameters
        A_VJ = 1.2e-9  # V/sqrt(Hz) noise floor
        A_IJ = 2e-12  # A/sqrt(Hz) noise floor
        A_VF = 6e-9  # V/sqrt(Hz) @ 1 Hz
        A_IF = 25e-12  # A/sqrt(Hz) @ 1 Hz

        R_effective = self.resistance() / 2
        return (A_VJ**2 + 2 * (R_effective * A_IJ) ** 2) + (A_VF**2 + 2 * (R_effective * A_IF) ** 2) / freq

    def amplifier_integrated(self) -> float:
        """Calculate integrated amplifier noise.

        Returns:
            Integrated noise voltage (V)
        """
        freq = np.logspace(np.log10(self.freq_min), np.log10(self.freq_max), self.numFrequencyPoints)
        return float(np.sqrt(np.trapezoid(self.amplifier_PSD(freq), freq)))

    def actuator_noise_integrated(self) -> float:
        """Calculate integrated actuator/mounting noise.

        Returns:
            Integrated noise voltage (V)
        """
        rms_displacement_noise = 1e-9  # m
        return rms_displacement_noise * self.stiffness() * self.force_sensitivity()

    def knee_frequency(self) -> float:
        """Calculate knee frequency where Johnson = Hooge noise.

        Returns:
            Knee frequency (Hz)
        """
        return (
            self.number_of_piezoresistors
            * self.alpha
            * self.v_bridge**2
            / (16 * self.number_of_carriers() * self.k_b * self.T * self.resistance())
        )

    def integrated_noise(self) -> float:
        """Calculate total integrated noise voltage.

        Returns:
            Integrated noise voltage (V)
        """
        freq = np.logspace(np.log10(self.freq_min), np.log10(self.freq_max), self.numFrequencyPoints)
        actuator = self.actuator_noise_integrated()
        total_psd = self.johnson_PSD(freq) + self.hooge_PSD(freq) + self.thermo_PSD(freq) + self.amplifier_PSD(freq)
        return float(np.sqrt(actuator**2 + np.trapezoid(total_psd, freq)))

    def piezo_coefficient(self) -> float:
        """Calculate piezoresistive coefficient.

        Uses Richter's model for poly-Si, geometric factor for metals.

        Returns:
            Piezoresistive coefficient (1/Pa)
        """
        if self.matl_top == Material.POLY:
            # Pi at low concentration in <100> direction for n-type
            max_factor = 103e-11
            # Polycrystalline correction (60%)
            max_factor *= 0.6

            # Richter's model (T=300K)
            Nb = 6e19
            Nc = 7e20
            richter_alpha = 0.43
            richter_gamma = 1.6
            piezoresistance_factor = (
                1
                + (self.dopant_concentration / Nb) ** richter_alpha
                + (self.dopant_concentration / Nc) ** richter_gamma
            ) ** (-1)
        elif self.matl_top == Material.AL:
            poisson = 0.35
            max_factor = (1 + 2 * poisson) / self.E_al
            piezoresistance_factor = 1.0
        elif self.matl_top == Material.TI:
            poisson = 0.35
            max_factor = (1 + 2 * poisson) / self.E_ti
            piezoresistance_factor = 1.0
        else:
            # Default to poly-Si
            max_factor = 103e-11 * 0.6
            piezoresistance_factor = 1.0

        return max_factor * piezoresistance_factor

    def force_sensitivity(self) -> float:
        """Calculate force sensitivity.

        Returns:
            Force sensitivity (V/N)
        """
        v_bias = self.v_bridge / 2
        Zm = self.neutral_axis()
        Cm = self.normalized_curvature()
        Z_top = self.t_bot + self.t_mid + self.t_top / 2

        E = self._get_elastic_modulus(self.matl_top)

        longitudinal_sensitivity = (
            self.number_of_piezoresistors_on_cantilever
            * self.piezo_coefficient()
            * E
            * (Z_top - Zm)
            * Cm
            * (self.l - self.l_pr / 2)
            * v_bias
            / 4
            * self.gamma()
        )

        return abs(longitudinal_sensitivity)

    def power_dissipation(self) -> float:
        """Calculate power dissipation in the cantilever.

        Returns:
            Power dissipation (W)
        """
        return self.number_of_piezoresistors_on_cantilever * self.v_bridge**2 / (4 * self.resistance())

    def force_resolution(self) -> float:
        """Calculate minimum detectable force.

        Returns:
            Force resolution (N)
        """
        return self.integrated_noise() / self.force_sensitivity()

    def displacement_resolution(self) -> float:
        """Calculate minimum detectable displacement.

        Returns:
            Displacement resolution (m)
        """
        return self.force_resolution() / self.stiffness()

    def print_performance(self) -> None:
        """Print cantilever performance summary."""
        print(f"Cantilever L/W: {self.l * 1e6:.1f} {self.w * 1e6:.1f} um")
        print(f"Layers - Bottom:{self.matl_bot.value} Mid:{self.matl_mid.value} Top:{self.matl_top.value}")
        print(f"Thicknesses (nm) - Bottom:{self.t_bot * 1e9:.1f} Mid:{self.t_mid * 1e9:.1f} Top:{self.t_top * 1e9:.1f}")
        print(f"PR Length Ratio: {self.l_pr_ratio}")
        print(f"Bridge voltage: {self.v_bridge} V")
        print(f"Number of piezoresistors: {self.number_of_piezoresistors}")
        print()
        print(f"Resistance: {self.resistance():.1f} ohms")
        print(f"Power dissipation: {self.power_dissipation() * 1e3:.3f} mW")
        print()
        print(f"Force resolution: {self.force_resolution():.3g} N")
        print(f"Displacement resolution: {self.displacement_resolution():.3g} m")
        print(f"Sensitivity: {self.force_sensitivity():.3g} V/N")
        print()
        print(f"Integrated noise: {self.integrated_noise():.3g} V")
        print(f"Johnson noise: {self.johnson_integrated():.3g} V")
        print(f"1/f noise: {self.hooge_integrated():.3g} V")
        print(f"Amplifier noise: {self.amplifier_integrated():.3g} V")
        print(f"Knee frequency: {self.knee_frequency():.1f} Hz")
        print()
        print(f"Sheet resistance: {self.sheet_resistance():.1f} ohms/sq")
        print(f"Number of carriers: {self.number_of_carriers():.3g}")
        print(f"Stiffness: {self.stiffness():.3g} N/m")
        print(f"Resonant frequency: {self.resonant_frequency():.1f} Hz")
