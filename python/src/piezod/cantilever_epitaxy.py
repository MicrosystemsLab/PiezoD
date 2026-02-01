"""Cantilever with epitaxial doping profile.

This module models an epitaxial cantilever where doping is achieved through
epitaxial growth, resulting in a step-function doping profile. Supports boron,
phosphorus, and arsenic dopants with negligible diffusion.
"""

from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from piezod.cantilever import Cantilever


class CantileverEpitaxy(Cantilever):
    """Cantilever with epitaxial piezoresistor doping.

    Models a cantilever with epitaxially grown piezoresistor layer.
    Assumes negligible dopant diffusion, resulting in a sharp junction
    between the epitaxial layer and substrate.

    Attributes:
        dopant_concentration: Dopant concentration in epitaxial layer (cm^-3)
        t_pr_ratio: Ratio of piezoresistor thickness to total thickness (-)
    """

    def __init__(
        self,
        freq_min: float = 1.0,
        freq_max: float = 1e3,
        l: float = 100e-6,
        w: float = 10e-6,
        t: float = 1e-6,
        l_pr_ratio: float = 0.3,
        v_bridge: float = 1.0,
        doping_type: str = "phosphorus",
        dopant_concentration: float = 1e19,
        t_pr_ratio: float = 0.3,
    ) -> None:
        """Initialize epitaxial cantilever.

        Args:
            freq_min: Minimum frequency for analysis (Hz)
            freq_max: Maximum frequency for analysis (Hz)
            l: Cantilever length (m)
            w: Cantilever width (m)
            t: Cantilever thickness (m)
            l_pr_ratio: Ratio of piezoresistor length to cantilever length (-)
            v_bridge: Bridge voltage (V)
            doping_type: Dopant type - "boron", "phosphorus", or "arsenic"
            dopant_concentration: Dopant concentration in epitaxial layer (cm^-3)
            t_pr_ratio: Ratio of piezoresistor thickness to total thickness (-)
        """
        super().__init__()

        # Set base parameters
        self.freq_min = freq_min
        self.freq_max = freq_max
        self.l = l
        self.w = w
        self.t = t
        self.l_pr_ratio = l_pr_ratio
        self.v_bridge = v_bridge
        self.doping_type = doping_type

        # Set epitaxy-specific parameters
        self.dopant_concentration = dopant_concentration
        self.t_pr_ratio = t_pr_ratio

    @property
    def junction_depth(self) -> float:
        """Calculate junction depth of epitaxial layer.

        Returns:
            Junction depth (m)
        """
        return self.t * self.t_pr_ratio

    def doping_profile(self) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """Calculate the electrically active dopant concentration profile.

        Returns a step-function profile with constant doping in the epitaxial
        layer and background concentration in the substrate.

        Returns:
            Tuple of (z, active_doping, total_doping) where:
                z: Depth coordinates from surface (m)
                active_doping: Electrically active dopant concentration (cm^-3)
                total_doping: Total dopant concentration (cm^-3)
        """
        n_points = self.numZPoints
        z = np.linspace(0, self.t, n_points)

        # Initialize to background concentration
        background_concentration = 1e15
        active_doping = np.ones(n_points) * background_concentration

        # Fill in epitaxial region (assume active = total for epitaxy)
        active_doping[z <= self.junction_depth] = self.dopant_concentration
        total_doping = active_doping.copy()

        return z, active_doping, total_doping

    def sheet_resistance(self) -> float:
        """Calculate sheet resistance of the piezoresistor.

        Uses the conductivity at the epitaxial doping concentration and
        junction depth to compute sheet resistance.

        Returns:
            Sheet resistance (ohms/square)
        """
        # Get conductivity in ohm-cm
        conductivity = self.conductivity(self.dopant_concentration)

        # Convert junction depth to cm and calculate sheet resistance
        junction_depth_cm = self.junction_depth * 1e2
        Rs = 1.0 / (junction_depth_cm * conductivity)

        return Rs

    def Nz(self) -> float:
        """Calculate integrated carrier concentration per unit area.

        For a constant concentration epitaxial layer, this is simply the
        product of junction depth and dopant concentration.

        Returns:
            Integrated carrier concentration (carriers/m^2)
        """
        return self.junction_depth * self.dopant_concentration * 1e6

    def alpha(self) -> float:
        """Calculate Hooge noise parameter.

        For epitaxial layers, uses the default Hooge parameter from the
        base class.

        Returns:
            Hooge parameter (dimensionless)
        """
        return self.default_alpha

    def doping_optimization_scaling(self) -> NDArray[np.float64]:
        """Get optimization scaling factors for doping parameters.

        Scaling factors help optimization algorithms handle parameters
        with different orders of magnitude.

        Returns:
            Array of [concentration_scale, t_pr_ratio_scale]
        """
        concentration_scale = 1e-19
        t_pr_ratio_scale = 10.0
        return np.array([concentration_scale, t_pr_ratio_scale])

    def doping_cantilever_from_state(self, x0: NDArray[np.float64]) -> None:  # type: ignore[override]
        """Reconstruct cantilever doping parameters from optimization state.

        Args:
            x0: Optimization state vector with doping parameters at indices 5-6
        """
        self.dopant_concentration = x0[5]
        self.t_pr_ratio = x0[6]

    def doping_current_state(self) -> NDArray[np.float64]:
        """Get current doping state for optimization.

        Returns:
            Array of [dopant_concentration, t_pr_ratio]
        """
        return np.array([self.dopant_concentration, self.t_pr_ratio])

    def doping_optimization_bounds(
        self, parameter_constraints: Optional[dict] = None
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Get optimization bounds for doping parameters.

        Default bounds based on solid solubility limits at ~1000°C.
        Constraints can override defaults using a dictionary.

        Args:
            parameter_constraints: Optional dict to override defaults.
                Example: {'min_dopant_concentration': 1e18, 'max_t_pr_ratio': 0.5}

        Returns:
            Tuple of (lower_bounds, upper_bounds) as numpy arrays
        """
        # Default bounds
        min_dopant_concentration = 1e17

        # Approximate solid solubilities at 1000C
        if self.doping_type == "boron":
            max_dopant_concentration = 2e20
        elif self.doping_type == "phosphorus":
            max_dopant_concentration = 4e20
        elif self.doping_type == "arsenic":
            max_dopant_concentration = 8e20
        else:
            # Fallback for unknown dopant types
            max_dopant_concentration = 1e20

        min_t_pr_ratio = 0.01
        max_t_pr_ratio = 0.99

        # Override with provided constraints
        if parameter_constraints:
            min_dopant_concentration = parameter_constraints.get("min_dopant_concentration", min_dopant_concentration)
            max_dopant_concentration = parameter_constraints.get("max_dopant_concentration", max_dopant_concentration)
            min_t_pr_ratio = parameter_constraints.get("min_t_pr_ratio", min_t_pr_ratio)
            max_t_pr_ratio = parameter_constraints.get("max_t_pr_ratio", max_t_pr_ratio)

        lb = np.array([min_dopant_concentration, min_t_pr_ratio])
        ub = np.array([max_dopant_concentration, max_t_pr_ratio])

        return lb, ub

    def doping_initial_conditions_random(self, parameter_constraints: Optional[dict] = None) -> NDArray[np.float64]:
        """Generate random initial conditions for optimization.

        Dopant concentration is sampled logarithmically (uniform in log space)
        to cover the wide range of possible values. t_pr_ratio is sampled
        uniformly.

        Args:
            parameter_constraints: Optional dict to override default bounds

        Returns:
            Array of [dopant_concentration, t_pr_ratio] sampled randomly
                within bounds
        """
        lb, ub = self.doping_optimization_bounds(parameter_constraints)

        n_min = lb[0]
        n_max = ub[0]
        t_pr_ratio_min = lb[1]
        t_pr_ratio_max = ub[1]

        # Generate random values
        # Dopant concentration: logarithmically distributed
        log_n_min = np.log10(n_min)
        log_n_max = np.log10(n_max)
        dopant_concentration_random = 10 ** (log_n_min + np.random.rand() * (log_n_max - log_n_min))

        # t_pr_ratio: uniformly distributed
        t_pr_ratio_random = t_pr_ratio_min + np.random.rand() * (t_pr_ratio_max - t_pr_ratio_min)

        x0 = np.array([dopant_concentration_random, t_pr_ratio_random])

        return x0
