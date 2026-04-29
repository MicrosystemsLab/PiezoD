"""Cantilever with diffused doping profile.

This module models piezoresistive cantilevers with diffused dopant profiles,
particularly for POCl3 phosphorus diffusion, arsenic, and boron diffusion.
The diffusion profiles are calculated using complementary error functions
based on Fick's laws of diffusion.

References:
    - Fair, R.B. (1981). "Concentration profiles of diffused dopants in silicon."
      In F.F.Y. Wang (Ed.), Impurity Doping Processes in Silicon (pp. 315-442).
    - Tsai, M. Y. (1983). "Modeling of phosphorus diffusion from POCl3 source."
      Journal of the Electrochemical Society, 130(10), 2095-2099.
"""

import numpy as np
import numpy.typing as npt
from scipy.special import erfc

from piezod.cantilever import Cantilever
from piezod.optimization.state import StateVar


class CantileverDiffusion(Cantilever):
    """Model a diffused cantilever (particularly phosphorus/POCl3 diffusion).

    This class extends the base Cantilever class to model cantilevers with
    diffused piezoresistor doping profiles. It supports arsenic, boron, and
    phosphorus diffusion with temperature-dependent diffusion coefficients.

    Attributes:
        diffusion_time: Diffusion time in seconds.
        diffusion_temp: Diffusion temperature in Kelvin.
    """

    def __init__(
        self,
        freq_min: float,
        freq_max: float,
        l: float,
        w: float,
        t: float,
        l_pr_ratio: float,
        v_bridge: float,
        doping_type: str,
        diffusion_time: float,
        diffusion_temp: float,
    ) -> None:
        """Initialize a CantileverDiffusion instance.

        Args:
            freq_min: Minimum frequency of interest (Hz).
            freq_max: Maximum frequency of interest (Hz).
            l: Cantilever length (m).
            w: Cantilever width (m).
            t: Cantilever thickness (m).
            l_pr_ratio: Piezoresistor length to cantilever length ratio (dimensionless).
            v_bridge: Wheatstone bridge bias voltage (V).
            doping_type: Dopant type: "arsenic", "boron", or "phosphorus".
            diffusion_time: Diffusion time (s).
            diffusion_temp: Diffusion temperature (K).
        """
        super().__init__()

        # Set basic cantilever parameters
        self.freq_min = freq_min
        self.freq_max = freq_max
        self.l = l
        self.w = w
        self.t = t
        self.l_pr_ratio = l_pr_ratio
        self.v_bridge = v_bridge
        self.doping_type = doping_type

        # Set diffusion-specific parameters
        self.diffusion_time = diffusion_time
        self.diffusion_temp = diffusion_temp

    def doping_profile(self) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Calculate the diffusion profile for a constant surface source.

        The doping profile is calculated using complementary error function (erfc)
        solutions to Fick's diffusion equation. For phosphorus, a more sophisticated
        model accounting for concentration-dependent diffusion is used.

        Returns:
            Tuple of (x, active_doping, total_doping) where:
                x: Position array from surface into silicon (m).
                active_doping: Net active carriers of the resistor type (cm^-3) =
                    ``max(0, electrically_active_dopant - substrate_background_cm3)``.
                    Goes to zero below the junction.
                total_doping: Total concentration ``max(dopant_species,
                    substrate_background_cm3)`` (cm^-3); floors at the substrate
                    so the profile is directly plottable.

        References:
            - Arsenic/Boron: Simple erfc profile from Fick's law
            - Phosphorus: Tsai model for POCl3 diffusion (J. Electrochem. Soc., 1983)
        """
        N_surface = 1e20  # Surface concentration for arsenic/boron (1/cm^3)
        n_points = self.numZPoints
        background = self.substrate_background_cm3

        if self.doping_type == "arsenic":
            # Arsenic diffusion parameters
            D_0 = 9.17  # Pre-exponential factor (cm^2/s)
            E_a = 3.99  # Activation energy (eV)

            # Temperature-dependent diffusion coefficient
            D = D_0 * np.exp(-E_a / (self.k_b_eV * self.diffusion_temp))
            diffusion_length = np.sqrt(D * self.diffusion_time) * 1e-2  # cm -> m

            # Position array
            x = np.linspace(0, self.t, n_points)

            # Complementary error function profile
            dopant = N_surface * erfc(x / (2 * diffusion_length))
            total_doping = np.maximum(dopant, background)
            active_doping = np.maximum(0.0, dopant - background)

        elif self.doping_type == "boron":
            # Boron diffusion parameters
            D_0 = 1.0  # Pre-exponential factor (cm^2/s)
            E_a = 3.5  # Activation energy (eV)

            # Temperature-dependent diffusion coefficient
            D = D_0 * np.exp(-E_a / (self.k_b_eV * self.diffusion_temp))
            diffusion_length = np.sqrt(D * self.diffusion_time) * 1e-2  # cm -> m

            # Position array
            x = np.linspace(0, self.t, n_points)

            # Complementary error function profile
            dopant = N_surface * erfc(x / (2 * diffusion_length))
            total_doping = np.maximum(dopant, background)
            active_doping = np.maximum(0.0, dopant - background)

        elif self.doping_type == "phosphorus":
            # Phosphorus diffusion using Tsai model for POCl3
            k_b_eV = 8.617343e-5  # Boltzmann constant (eV/K)
            x = np.linspace(0, self.t * 1e2, n_points)  # m -> cm

            # Time offset for oxidation predeposition
            time_offset = 0 if self.diffusion_temp < 780 + 273 else 2  # minutes

            # Alpha activation energy varies quadratically with temperature
            T_values = np.array([273 + 775, 273 + 800, 273 + 850, 273 + 950])
            Ea_values = np.array([1.755, 1.73, 1.71, 1.66])
            p = np.polyfit(T_values, Ea_values, 2)
            Ea_alpha = np.polyval(p, self.diffusion_temp)
            alpha = 0.18 * np.exp(-Ea_alpha / (k_b_eV * self.diffusion_temp))

            # Two-component diffusion model parameters
            Da = 100 * np.exp(-3.77 / k_b_eV / self.diffusion_temp)  # Fast diffusion (cm^2/s)
            Db = 2.3 * np.exp(-1.95 / (k_b_eV * self.diffusion_temp)) * 1e-5  # Slow diffusion (cm^2/s)
            Cb = 3 * np.exp(-0.88 / (k_b_eV * self.diffusion_temp)) * 1e23  # Slow component concentration (1/cm^3)

            # Effective diffusion time including offset
            t = self.diffusion_time + time_offset * 60  # Convert minutes to seconds

            # Surface concentrations (temperature-dependent)
            surface_concentration_total = 2.5e23 * np.exp(-0.62 / (k_b_eV * self.diffusion_temp))
            surface_concentration_active = 1.1e22 * np.exp(-0.37 / (k_b_eV * self.diffusion_temp))

            # Oxide layer thickness parameter
            x0 = alpha * t
            kappa = Cb / surface_concentration_total

            # Two-component erfc model
            F1 = erfc((x + alpha * t) / (2 * np.sqrt(Da * t))) + erfc((x - 3 * alpha * t) / (2 * np.sqrt(Da * t)))
            F2 = erfc((x + alpha * t) / (2 * np.sqrt(Db * t))) + erfc((x - 3 * alpha * t) / (2 * np.sqrt(Db * t)))

            Ca = (1 - kappa) / 2 * surface_concentration_total * np.exp(-alpha / (2 * Da) * (x - alpha * t)) * F1
            Cb_profile = kappa / 2 * surface_concentration_total * np.exp(-alpha / (2 * Db) * (x - alpha * t)) * F2

            # Dopant species (no floor) and active concentration limited by
            # solid solubility.
            dopant = Ca + Cb_profile
            dopant[x <= x0] = surface_concentration_total
            total_doping = np.maximum(dopant, background)

            C_active = np.minimum(dopant, surface_concentration_active)
            active_doping = np.maximum(0.0, C_active - background)

            # Convert back to meters
            x = x * 1e-2  # cm -> m

        else:
            raise ValueError(f"Unknown doping type: {self.doping_type}")

        return x, active_doping, total_doping

    def sheet_resistance(self) -> float:
        """Calculate sheet resistance of the doped layer.

        The sheet resistance is calculated by integrating the conductivity
        profile through the thickness of the doped layer.

        Returns:
            Sheet resistance (ohms/square).
        """
        x, active_doping, total_doping = self.doping_profile()
        conductivity = self.conductivity(active_doping)
        # Integrate conductivity using trapezoidal rule (x in m, convert to cm)
        Rs = 1 / np.trapezoid(conductivity, x * 1e2)
        return Rs

    def Nz(self) -> float:
        """Calculate effective carrier density.

        The effective carrier density is a mobility-weighted average that
        accounts for the non-uniform doping profile.

        Returns:
            Effective carrier density (1/m^2).

        References:
            Weighted average: Nz = (integral(N*mu*dz))^2 / integral(N*mu^2*dz)
        """
        z, N_active, N_total = self.doping_profile()
        mu, sigma = self.mobility(N_active, self.T)

        # Calculate mobility-weighted effective density
        # Convert z from m to cm for consistency
        z_cm = z * 1e2

        numerator = np.trapezoid(N_active * mu, z_cm) ** 2
        denominator = np.trapezoid(N_active * mu**2, z_cm)

        # Convert from 1/cm^2 to 1/m^2
        Nz = 1e4 * numerator / denominator
        return Nz

    @property
    def junction_depth(self) -> float:
        """Junction depth where active carriers fall to zero.

        Linear interpolation of where ``active_doping`` crosses zero (i.e. where
        the dopant species drops to ``substrate_background_cm3``). Returns the
        deepest sample if ``active_doping`` never reaches zero.
        """
        x, active_doping, _ = self.doping_profile()
        below = active_doping <= 0
        if not np.any(below):
            return float(x[-1])
        first = int(np.argmax(below))
        if first == 0:
            return float(x[0])
        # Linear interpolation between the bracketing samples (active_doping
        # is positive at first-1 and zero at first).
        x0, x1 = x[first - 1], x[first]
        a0, a1 = active_doping[first - 1], active_doping[first]
        if a0 == a1:
            return float(x0)
        return float(x0 + (x1 - x0) * a0 / (a0 - a1))

    def alpha(self) -> float:
        """Return Hooge parameter for 1/f noise modeling.

        Returns:
            Hooge parameter (dimensionless), typically 1e-5.
        """
        return self.default_alpha

    # ========= Optimization Methods ==========

    def doping_optimization_scaling(self) -> list[float]:
        """Return scaling factors for optimization parameters.

        These scaling factors help normalization in optimization algorithms.

        Returns:
            List of [diffusion_time_scale, diffusion_temp_scale].
        """
        diffusion_time_scale = 1e-3
        diffusion_temp_scale = 1e-3
        return [diffusion_time_scale, diffusion_temp_scale]

    def doping_cantilever_from_state(self, x0: list[float]) -> None:
        """Update cantilever diffusion parameters from optimization state vector.

        Args:
            x0: State vector with elements [l, w, t, l_pr_ratio, v_bridge,
                diffusion_time, diffusion_temp].
        """
        self.diffusion_time = x0[5]
        self.diffusion_temp = x0[6]

    def doping_current_state(self) -> list[float]:
        """Return current diffusion parameters as a state vector.

        Returns:
            List of [diffusion_time, diffusion_temp].
        """
        return [self.diffusion_time, self.diffusion_temp]

    def doping_optimization_bounds(
        self, parameter_constraints: dict[str, float] | None = None
    ) -> tuple[list[float], list[float]]:
        """Return bounds for diffusion optimization parameters.

        Args:
            parameter_constraints: Optional dict of parameter constraints.
                Keys can include: min_diffusion_time, max_diffusion_time,
                min_diffusion_temp, max_diffusion_temp.

        Returns:
            Tuple of (lower_bounds, upper_bounds) for [diffusion_time, diffusion_temp].
        """
        # Default bounds
        min_diffusion_time = 5 * 60  # 5 minutes in seconds
        max_diffusion_time = 90 * 60  # 90 minutes in seconds
        min_diffusion_temp = 273 + 800  # K
        max_diffusion_temp = 273 + 1000  # K

        # Override with user constraints if provided
        if parameter_constraints is not None:
            min_diffusion_time = parameter_constraints.get("min_diffusion_time", min_diffusion_time)
            max_diffusion_time = parameter_constraints.get("max_diffusion_time", max_diffusion_time)
            min_diffusion_temp = parameter_constraints.get("min_diffusion_temp", min_diffusion_temp)
            max_diffusion_temp = parameter_constraints.get("max_diffusion_temp", max_diffusion_temp)

        lb = [min_diffusion_time, min_diffusion_temp]
        ub = [max_diffusion_time, max_diffusion_temp]
        return lb, ub

    def doping_initial_conditions_random(self, parameter_constraints: dict[str, float] | None = None) -> list[float]:
        """Generate random initial conditions within bounds for optimization.

        Args:
            parameter_constraints: Optional dict of parameter constraints.

        Returns:
            List of random [diffusion_time, diffusion_temp] within bounds.
        """
        lb, ub = self.doping_optimization_bounds(parameter_constraints)

        diffusion_time_min, diffusion_temp_min = lb
        diffusion_time_max, diffusion_temp_max = ub

        diffusion_time_random = diffusion_time_min + np.random.rand() * (diffusion_time_max - diffusion_time_min)
        diffusion_temp_random = diffusion_temp_min + np.random.rand() * (diffusion_temp_max - diffusion_temp_min)

        return [diffusion_time_random, diffusion_temp_random]

    def optimization_state_vars(self) -> tuple[StateVar, ...]:
        """Declarative state spec for joint geometry + diffusion optimization.

        Seven state variables: cantilever length, width, thickness,
        piezoresistor length ratio, bridge bias voltage, diffusion time,
        and diffusion temperature. Default bounds match MATLAB
        ``cantileverDiffusion``.
        """
        return (
            StateVar("l", 1e5, 10e-6, 3e-3),
            StateVar("w", 1e7, 2e-6, 100e-6),
            StateVar("t", 1e8, 1e-6, 100e-6),
            StateVar("l_pr_ratio", 1e2, 0.01, 0.99),
            StateVar("v_bridge", 1e1, 0.1, 10.0),
            StateVar("diffusion_time", 1e-3, 5 * 60, 90 * 60),
            StateVar("diffusion_temp", 1e-3, 273 + 800, 273 + 1000),
        )
