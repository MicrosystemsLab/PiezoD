"""Piezoelectric cantilever with AlN or PZT materials.

This module models a piezoelectric cantilever sensor with multi-layer beam mechanics,
supporting both voltage-mode and charge-mode sensing. The Sader hydrodynamic damping
model is used for calculating Q factor in fluid environments.
"""

from enum import Enum
from typing import Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize_scalar


class PiezoMaterial(Enum):
    """Piezoelectric material options."""

    ALN = "AlN"
    PZT = "PZT"


class FluidType(Enum):
    """Fluid environment options."""

    VACUUM = "vacuum"
    AIR = "air"
    WATER = "water"


class CantileverPiezoelectric:
    """Piezoelectric cantilever with multi-layer beam mechanics.

    Models a cantilever with silicon substrate, electrodes, and piezoelectric layer.
    Supports both AlN and PZT piezoelectric materials, and vacuum/air/water fluids.

    Attributes:
        freq_min: Minimum frequency for analysis (Hz)
        freq_max: Maximum frequency for analysis (Hz)
        l_si: Silicon cantilever length (m)
        w_si: Silicon cantilever width (m)
        t_si: Silicon thickness (m)
        t_pe: Piezoelectric layer thickness (m)
        material: Piezoelectric material type (AlN or PZT)
        r_shunt: Shunt resistance (ohms)
        fluid: Fluid environment type
    """

    # Physical constants
    EPSILON_0 = 8.854e-12  # F/m
    K_B = 1.38e-23  # J/K
    T = 300  # K
    Q_E = 1.6e-19  # C

    # Fluid properties
    RHO_WATER = 1e3  # kg/m^3
    ETA_WATER = 0.9e-3  # Pa-s
    RHO_AIR = 1.2  # kg/m^3
    ETA_AIR = 17e-6  # Pa-s

    # Q factor limits
    MAX_Q = 1000.0
    MIN_Q = 1e-6

    # Number of plot points
    N_PLOT_POINTS = 1000

    # Silicon properties
    RHO_SI = 2330  # kg/m^3
    E_SI = 169e9  # Pa

    # Electrode properties (titanium)
    T_BOTTOM_ELECTRODE = 50e-9  # m
    T_TOP_ELECTRODE = 50e-9  # m
    RHO_ELECTRODE = 4506  # kg/m^3
    E_ELECTRODE = 90e9  # Pa

    # Amplifier properties
    C_AMP = 0.2e-12  # F

    # Material properties lookup
    MATERIAL_PROPERTIES = {
        PiezoMaterial.ALN: {
            "permittivity": 10 * EPSILON_0,  # F/m
            "resistivity": 1e12,  # ohm-cm
            "d31": 2e-12,  # C/N
            "E": 396e9,  # Pa
            "rho": 3260,  # kg/m^3
        },
        PiezoMaterial.PZT: {
            "permittivity": 900 * EPSILON_0,  # F/m
            "resistivity": 1e8,  # ohm-cm
            "d31": 70e-12,  # C/N
            "E": 55e9,  # Pa
            "rho": 7550,  # kg/m^3
        },
    }

    # Sader lookup table for hydrodynamic function
    # kappa = C * w / l, where C = 1.8751 for first mode
    KAPPA_LOOKUP = np.array([0, 0.125, 0.25, 0.5, 0.75, 1, 2, 3, 5, 7, 10, 20])
    REYNOLDS_LOOKUP = np.array([-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4])

    # fmt: off
    TAU_LOOKUP_REAL = np.array([
        [3919.41, 59.3906, 22.4062, 9.13525, 5.62175, 4.05204, 1.93036, 1.2764, 0.764081, 0.545683, 0.381972, 0.190986],
        [1531.90, 59.3861, 22.4061, 9.13525, 5.62175, 4.05204, 1.93036, 1.2764, 0.764081, 0.545683, 0.381972, 0.190986],
        [613.426, 59.3420, 22.4052, 9.13523, 5.62174, 4.05204, 1.93036, 1.2764, 0.764081, 0.545683, 0.381972, 0.190986],
        [253.109, 58.9094, 22.3962, 9.13504, 5.62172, 4.05203, 1.93036, 1.2764, 0.764081, 0.545683, 0.381972, 0.190986],
        [108.429, 55.2882, 22.3078, 9.13319, 5.62153, 4.05199, 1.93036, 1.2764, 0.764081, 0.545683, 0.381972, 0.190986],
        [48.6978, 40.7883, 21.5187, 9.11481, 5.61960, 4.05160, 1.93035, 1.2764, 0.764081, 0.545683, 0.381972, 0.190986],
        [23.2075, 22.7968, 17.5378, 8.94370, 5.60057, 4.04771, 1.93027, 1.27639, 0.76408, 0.545683, 0.381972, 0.190986],
        [11.8958, 11.9511, 11.0719, 7.89716, 5.43378, 4.01051, 1.92942, 1.27629, 0.764074, 0.545682, 0.381972, 0.190986],
        [6.64352, 6.64381, 6.47227, 5.65652, 4.64017, 3.74600, 1.92114, 1.27536, 0.764012, 0.545671, 0.38197, 0.190986],
        [4.07692, 4.05940, 3.99256, 3.72963, 3.37543, 3.00498, 1.85532, 1.26646, 0.763397, 0.545564, 0.381953, 0.190985],
        [2.74983, 2.73389, 2.69368, 2.56390, 2.39884, 2.22434, 1.61821, 1.20592, 0.757637, 0.54452, 0.381786, 0.190981],
        [2.02267, 2.01080, 1.98331, 1.90040, 1.79834, 1.69086, 1.31175, 1.04626, 0.721165, 0.535593, 0.38018, 0.190932],
        [1.60630, 1.59745, 1.57723, 1.51690, 1.44276, 1.36416, 1.08036, 0.878177, 0.635443, 0.496169, 0.368548, 0.190459],
        [1.36230, 1.35532, 1.33934, 1.29142, 1.23203, 1.16842, 0.932812, 0.759965, 0.551349, 0.435586, 0.334279, 0.186672],
        [1.21727, 1.21141, 1.19792, 1.15718, 1.10624, 1.05117, 0.84292, 0.686229, 0.493924, 0.387183, 0.295972, 0.172722],
        [1.13038, 1.12518, 1.11316, 1.07668, 1.03073, 0.980721, 0.78879, 0.641744, 0.458699, 0.356289, 0.268907, 0.154450],
        [1.07814, 1.07334, 1.06221, 1.02827, 0.985314, 0.938346, 0.756309, 0.615164, 0.437743, 0.337813, 0.252327, 0.140852],
    ])

    TAU_LOOKUP_IMAG = np.array([
        [27984.8, 44628.5, 55176.1, 71754, 86311.5, 100062, 152411, 203623, 305570, 407436, 560225, 1069521],
        [9816.33, 14113.1, 17448.2, 22690.6, 27294.1, 31642.3, 48196.5, 64391.4, 96629.7, 128843, 177159, 338212],
        [3482.47, 4464.16, 5517.72, 7175.41, 8631.15, 10006.2, 15241.1, 20362.3, 30557, 40743.6, 56022.5, 106952],
        [1252.66, 1415.42, 1745.17, 2269.09, 2729.42, 3164.23, 4819.65, 6439.14, 9662.97, 12884.3, 17715.9, 33821.2],
        [458.386, 457.863, 552.862, 717.635, 863.138, 1000.63, 1524.11, 2036.23, 3055.7, 4074.36, 5602.25, 10695.2],
        [171.397, 160.951, 177.702, 227.205, 273.013, 316.449, 481.967, 643.915, 966.297, 1288.43, 1771.59, 3382.12],
        [65.8679, 62.2225, 61.626, 72.6542, 86.5364, 100.144, 152.418, 203.625, 305.57, 407.436, 560.225, 1069.52],
        [26.2106, 25.21, 24.1432, 24.7484, 27.9459, 31.8957, 48.2199, 64.3973, 96.6308, 128.843, 177.159, 338.212],
        [10.8983, 10.6158, 10.1909, 9.7009, 9.91067, 10.648, 15.3139, 20.381, 30.5604, 40.7448, 56.0229, 106.952],
        [4.78389, 4.69492, 4.53952, 4.24925, 4.09701, 4.09433, 5.01844, 6.49605, 9.67379, 12.8879, 17.7171, 33.8214],
        [2.23883, 2.20681, 2.14583, 2.0088, 1.89659, 1.82463, 1.85993, 2.17718, 3.08849, 4.08581, 5.60598, 10.6956],
        [1.12164, 1.10851, 1.08208, 1.01654, 0.953355, 0.901676, 0.81464, 0.844519, 1.04394, 1.32116, 1.78306, 3.38349],
        [0.596697, 0.590686, 0.578118, 0.545082, 0.510467, 0.479247, 0.403803, 0.383595, 0.409256, 0.469688, 0.589749, 1.07377],
        [0.332285, 0.329276, 0.32283, 0.305262, 0.285953, 0.26763, 0.216732, 0.194409, 0.186218, 0.195634, 0.221631, 0.349855],
        [0.191043, 0.189434, 0.185931, 0.176166, 0.165118, 0.154323, 0.122124, 0.105573, 0.0938839, 0.0925686, 0.09682, 0.126835],
        [0.112082, 0.111181, 0.109199, 0.103595, 0.0971392, 0.0907188, 0.0707736, 0.059728, 0.0505049, 0.0476557, 0.0471326, 0.0534759],
        [0.0665172, 0.0659974, 0.0648471, 0.0615627, 0.0577366, 0.0538889, 0.0416384, 0.0345727, 0.0282418, 0.025856, 0.024611, 0.0252877],
    ])
    # fmt: on

    def __init__(
        self,
        freq_min: float = 1.0,
        freq_max: float = 1e3,
        l_si: float = 100e-6,
        w_si: float = 10e-6,
        t_si: float = 1e-6,
        t_pe: float = 500e-9,
        material: PiezoMaterial = PiezoMaterial.ALN,
        r_shunt: float = 1e12,
        fluid: FluidType = FluidType.AIR,
        amplifier_noise: bool = True,
        thermomechanical_noise: bool = True,
    ) -> None:
        """Initialize piezoelectric cantilever.

        Args:
            freq_min: Minimum frequency for analysis (Hz)
            freq_max: Maximum frequency for analysis (Hz)
            l_si: Silicon cantilever length (m)
            w_si: Silicon cantilever width (m)
            t_si: Silicon thickness (m)
            t_pe: Piezoelectric layer thickness (m)
            material: Piezoelectric material type
            r_shunt: Shunt resistance (ohms)
            fluid: Fluid environment type
            amplifier_noise: Include amplifier noise in calculations
            thermomechanical_noise: Include thermomechanical noise in calculations
        """
        self.freq_min = freq_min
        self.freq_max = freq_max
        self.l_si = l_si
        self.w_si = w_si
        self.t_si = t_si
        self.l_pe = l_si
        self.w_pe = w_si
        self.t_pe = t_pe
        self.material = material
        self.r_shunt = r_shunt
        self.fluid = fluid
        self.amplifier_noise = amplifier_noise
        self.thermomechanical_noise = thermomechanical_noise

        # Set material properties
        props = self.MATERIAL_PROPERTIES[material]
        self.permittivity = props["permittivity"]
        self.resistivity = props["resistivity"]
        self.d31 = props["d31"]
        self.E_pe = props["E"]
        self.rho_pe = props["rho"]

        # Create interpolators for hydrodynamic lookup
        self._tau_real_interp = RegularGridInterpolator(
            (self.REYNOLDS_LOOKUP, self.KAPPA_LOOKUP),
            self.TAU_LOOKUP_REAL,
            method="linear",
            bounds_error=False,
            fill_value=None,
        )
        self._tau_imag_interp = RegularGridInterpolator(
            (self.REYNOLDS_LOOKUP, self.KAPPA_LOOKUP),
            self.TAU_LOOKUP_IMAG,
            method="linear",
            bounds_error=False,
            fill_value=None,
        )

    # ==================== Beam Mechanics ====================

    def _calculate_mechanics_parameters(
        self,
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """Calculate layer-wise mechanical parameters.

        Returns:
            Tuple of (z_centers, E_layers, A_layers, I_layers)
        """
        t = np.array([self.t_si, self.T_BOTTOM_ELECTRODE, self.t_pe, self.T_TOP_ELECTRODE])
        z = np.zeros(4)
        for i in range(4):
            z[i] = np.sum(t) - np.sum(t[i:]) + t[i] / 2
        E = np.array([self.E_SI, self.E_ELECTRODE, self.E_pe, self.E_ELECTRODE])
        A = self.w_si * t
        I = (self.w_si * t**3) / 12
        return z, E, A, I

    @property
    def neutral_axis(self) -> float:
        """Calculate neutral axis position (m)."""
        z, E, A, _ = self._calculate_mechanics_parameters()
        return float(np.sum(z * E * A) / np.sum(E * A))

    @property
    def normalized_curvature(self) -> float:
        """Calculate normalized curvature (1/N-m)."""
        zm = self.neutral_axis
        z, E, A, I = self._calculate_mechanics_parameters()
        z_offset = z - zm
        return float(1.0 / np.sum(E * (I + A * z_offset**2)))

    @property
    def stiffness(self) -> float:
        """Calculate spring constant (N/m)."""
        EI_effective = 1 / self.normalized_curvature
        return 3 * EI_effective / self.l_pe**3

    @property
    def effective_mass(self) -> float:
        """Calculate effective mass (kg)."""
        return (
            0.243
            * self.w_si
            * self.l_si
            * (
                self.t_si * self.RHO_SI
                + self.T_BOTTOM_ELECTRODE * self.RHO_ELECTRODE
                + self.t_pe * self.rho_pe
                + self.T_TOP_ELECTRODE * self.RHO_ELECTRODE
            )
        )

    @property
    def resonant_frequency(self) -> float:
        """Calculate resonant frequency in vacuum (Hz)."""
        return 1 / (2 * np.pi) * np.sqrt(self.stiffness / self.effective_mass)

    # ==================== R, C Calculations ====================

    @property
    def Rpe(self) -> float:
        """Calculate piezoelectric layer resistance (ohms)."""
        area = self.l_pe * self.w_pe
        return self.resistivity * 1e-2 * self.t_pe / area

    @property
    def Cpe(self) -> float:
        """Calculate piezoelectric layer capacitance (F)."""
        area = self.l_pe * self.w_pe
        return self.permittivity * area / self.t_pe

    @property
    def R_half(self) -> float:
        """Calculate half-circuit resistance: Rpe/2 || Rshunt (ohms)."""
        return self.Rpe * self.r_shunt / (self.Rpe + 2 * self.r_shunt)

    @property
    def C_half(self) -> float:
        """Calculate half-circuit capacitance (F)."""
        return 2 * (self.Cpe + self.C_AMP)

    def Z(self, freq: NDArray[np.float64]) -> np.ndarray:
        """Calculate impedance Rpe || Cpe (ohms)."""
        return 1.0 / (1 / self.Rpe + 2j * np.pi * freq * self.Cpe)

    def Z_half(self, freq: NDArray[np.float64]) -> np.ndarray:
        """Calculate half-circuit impedance (ohms)."""
        return 1.0 / (1 / self.R_half + 2j * np.pi * freq * self.C_half)

    # ==================== Sensitivity ====================

    @property
    def q_f_sensitivity_intrinsic(self) -> float:
        """Calculate intrinsic charge sensitivity (C/N)."""
        zm = self.neutral_axis
        z, _, _, _ = self._calculate_mechanics_parameters()
        z_pe = z[2]
        cm = self.normalized_curvature
        return abs(0.5 * self.d31 * self.E_pe * (zm - z_pe) * cm * self.w_pe * self.l_pe**2)

    def q_f_sensitivity(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate charge force sensitivity (C/N)."""
        return np.ones_like(freq) * self.q_f_sensitivity_intrinsic

    def q_x_sensitivity(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate charge displacement sensitivity (C/m)."""
        return self.q_f_sensitivity(freq) * self.stiffness

    def v_f_sensitivity(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate voltage force sensitivity (V/N)."""
        omega = 2 * np.pi * freq
        return 2 * omega * self.q_f_sensitivity_intrinsic * self.R_half / (1 + omega * self.R_half * self.C_half)

    def v_x_sensitivity(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate voltage displacement sensitivity (V/m)."""
        return self.v_f_sensitivity(freq) * self.stiffness

    # ==================== Damping ====================

    def _lookup_fluid_properties(self) -> Tuple[float, float]:
        """Get fluid density and viscosity."""
        if self.fluid == FluidType.AIR:
            return self.RHO_AIR, self.ETA_AIR
        elif self.fluid == FluidType.WATER:
            return self.RHO_WATER, self.ETA_WATER
        else:
            return 1e-6, 1e-6

    @property
    def kappa(self) -> float:
        """Calculate kappa for first mode."""
        C = 1.8751
        return C * self.w_si / self.l_si

    def _reynolds(self, omega: float, rho_f: float, eta_f: float) -> float:
        """Calculate Reynolds number."""
        return rho_f * omega * self.w_si**2 / eta_f

    def _hydrodynamic_function(self, omega: float, rho_f: float, eta_f: float) -> complex:
        """Calculate hydrodynamic function from lookup table."""
        kappa = self.kappa
        reynolds = self._reynolds(omega, rho_f, eta_f)
        log_reynolds = np.log10(reynolds)

        # Clamp to table bounds
        log_reynolds = np.clip(log_reynolds, self.REYNOLDS_LOOKUP[0], self.REYNOLDS_LOOKUP[-1])
        kappa_clamped = np.clip(kappa, self.KAPPA_LOOKUP[0], self.KAPPA_LOOKUP[-1])

        tau_real = float(self._tau_real_interp((log_reynolds, kappa_clamped)))
        tau_imag = float(self._tau_imag_interp((log_reynolds, kappa_clamped)))

        return complex(tau_real, tau_imag)

    def omega_damped_and_Q(self) -> Tuple[float, float]:
        """Calculate damped natural frequency (rad/s) and Q factor."""
        if self.fluid == FluidType.VACUUM:
            return 2 * np.pi * self.resonant_frequency, self.MAX_Q

        rho_f, eta_f = self._lookup_fluid_properties()
        omega_vac = 2 * np.pi * self.resonant_frequency

        def residual(omega_d: float) -> float:
            if omega_d <= 0:
                return 1e10
            hydro = self._hydrodynamic_function(omega_d, rho_f, eta_f)
            target = omega_vac * (1 + np.pi * rho_f * self.w_si / (4 * self.RHO_SI * self.t_si) * hydro.real) ** -0.5
            return (omega_d - target) ** 2

        result = minimize_scalar(residual, bounds=(0, omega_vac), method="bounded")
        omega_damped = result.x

        hydro = self._hydrodynamic_function(omega_damped, rho_f, eta_f)
        Q = (4 * self.RHO_SI * self.t_si / (np.pi * rho_f * self.w_si) + hydro.real) / hydro.imag

        # Clamp Q to valid range
        if np.isnan(Q) or Q < self.MIN_Q:
            Q = self.MIN_Q
        elif Q > self.MAX_Q:
            Q = self.MAX_Q

        return omega_damped, Q

    # ==================== Noise ====================

    def Vth(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate thermomechanical voltage noise (V/sqrt(Hz))."""
        omega_damped, Q_M = self.omega_damped_and_Q()
        Fth = np.sqrt(2 * self.K_B * self.stiffness * self.T / (Q_M * np.pi * self.resonant_frequency))
        return Fth * self.v_f_sensitivity(freq)

    def Vamp(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate amplifier voltage noise (V/sqrt(Hz))."""
        # INA116 parameters
        Vamp_J = 28e-9
        Vamp_H = 300e-9
        Iamp_J = 0.1e-15
        Iamp_H = 1e-15

        Z_half_real = np.real(self.Z_half(freq))

        Vamp_VJ = Vamp_J
        Vamp_IJ = Iamp_J * Z_half_real
        Vamp_VH = Vamp_H / np.sqrt(freq)
        Vamp_IH = Iamp_H * Z_half_real / np.sqrt(freq)

        return np.sqrt(Vamp_VJ**2 + Vamp_IJ**2 + Vamp_VH**2 + Vamp_IH**2)

    def Vn(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate total voltage noise (V/sqrt(Hz))."""
        omega = 2 * np.pi * freq
        Vj = np.sqrt(4 * self.K_B * self.T * self.R_half) / (1 + omega * self.R_half * self.C_half)

        Vth = self.Vth(freq) if self.thermomechanical_noise else np.zeros_like(freq)
        Vamp = self.Vamp(freq) if self.amplifier_noise else np.zeros_like(freq)

        return np.sqrt(Vj**2 + Vth**2 + Vamp**2)

    def Qamp(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate amplifier charge noise (C/sqrt(Hz))."""
        # INA116 parameters
        Vamp_J = 28e-9
        Vamp_H = 300e-9
        Iamp_J = 0.1e-15
        Iamp_H = 1e-15

        omega = 2 * np.pi * freq

        Qamp_VJ = Vamp_J * self.Cpe
        Qamp_IJ = Iamp_J / omega
        Qamp_VH = Vamp_H * self.Cpe / np.sqrt(freq)
        Qamp_IH = Iamp_H / np.sqrt(freq) / omega

        return np.sqrt(Qamp_VJ**2 + Qamp_IJ**2 + Qamp_VH**2 + Qamp_IH**2)

    def Qn(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate total charge noise (C/sqrt(Hz))."""
        omega = 2 * np.pi * freq
        Qj = self.Cpe * np.sqrt(4 * self.K_B * self.T * self.Rpe) / (1 + omega * self.Rpe * self.Cpe)
        Qth = self.Vth(freq) * self.Cpe if self.thermomechanical_noise else np.zeros_like(freq)
        Qamp = self.Qamp(freq) if self.amplifier_noise else np.zeros_like(freq)

        return np.sqrt(Qj**2 + Qth**2 + Qamp**2)

    # ==================== Resolution ====================

    def freq_for_plot(self) -> NDArray[np.float64]:
        """Generate frequency array for calculations."""
        return np.logspace(np.log10(self.freq_min), np.log10(self.freq_max), self.N_PLOT_POINTS)

    def Sfminv(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Force noise density in voltage mode (N/sqrt(Hz))."""
        return self.Vn(freq) / self.v_f_sensitivity(freq)

    def Sxminv(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Displacement noise density in voltage mode (m/sqrt(Hz))."""
        return self.Vn(freq) / self.v_x_sensitivity(freq)

    def Sfminq(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Force noise density in charge mode (N/sqrt(Hz))."""
        return self.Qn(freq) / self.q_f_sensitivity(freq)

    def Sxminq(self, freq: NDArray[np.float64]) -> NDArray[np.float64]:
        """Displacement noise density in charge mode (m/sqrt(Hz))."""
        return self.Qn(freq) / self.q_x_sensitivity(freq)

    def Fminv(self) -> float:
        """Minimum detectable force in voltage mode (N)."""
        freq = self.freq_for_plot()
        Sfminv = self.Sfminv(freq)
        return float(np.sqrt(np.trapezoid(Sfminv**2, freq)))

    def Xminv(self) -> float:
        """Minimum detectable displacement in voltage mode (m)."""
        freq = self.freq_for_plot()
        Sxminv = self.Sxminv(freq)
        return float(np.sqrt(np.trapezoid(Sxminv**2, freq)))

    def Fminq(self) -> float:
        """Minimum detectable force in charge mode (N)."""
        freq = self.freq_for_plot()
        Sfminq = self.Sfminq(freq)
        return float(np.sqrt(np.trapezoid(Sfminq**2, freq)))

    def Xminq(self) -> float:
        """Minimum detectable displacement in charge mode (m)."""
        freq = self.freq_for_plot()
        Sxminq = self.Sxminq(freq)
        return float(np.sqrt(np.trapezoid(Sxminq**2, freq)))

    def thermomechanical_force_limit(self) -> float:
        """Calculate thermomechanical force noise limit (N)."""
        omega_damped, Q_M = self.omega_damped_and_Q()
        return np.sqrt(self.freq_max - self.freq_min) * np.sqrt(
            2 * self.K_B * self.stiffness * self.T / (Q_M * np.pi * self.resonant_frequency)
        )

    def Vn_total(self) -> float:
        """Total integrated voltage noise (V)."""
        freq = self.freq_for_plot()
        Vn = self.Vn(freq)
        return float(np.sqrt(np.trapezoid(Vn**2, freq)))

    def Qn_total(self) -> float:
        """Total integrated charge noise (C)."""
        freq = self.freq_for_plot()
        Qn = self.Qn(freq)
        return float(np.sqrt(np.trapezoid(Qn**2, freq)))
