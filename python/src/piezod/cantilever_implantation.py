"""Cantilever with ion-implanted doping profile.

Model an ion implanted cantilever using lookup tables from TSuprem.
B, P and As are supported with either inert or oxidizing anneal environments.

For all conditions, a 250A protection oxide layer is grown before the ion implantation.

After the implantation step:
- 'oxide' -> strip protection oxide, regrow 1500A oxide, inert anneal
- 'inert' -> leave protection oxide, inert anneal, strip oxide

Lookup tables are linearly interpolated. Splines can result in negative values
(e.g. sheet resistance). Gradient/Hessian discontinuities don't seem to pose a
problem for optimization.

The range of the simulation conditions is equal to the limits of the highest
quality data in TSuprem:
- Implantation energy: 20 - 80 keV
- Implantation dose: 2e14 - 2e16 per sq cm
- Annealing temp: 900 - 1100C (1173 - 1373K)
- Annealing time: 15 - 120 min
"""

from importlib import resources
from typing import Literal, Tuple

import h5py
import numpy as np
from scipy.interpolate import interpn

from piezod.cantilever import Cantilever

# Module-level cache for lookup table data (loaded once)
_LOOKUP_TABLE_CACHE: dict | None = None


def _load_lookup_table() -> dict:
    """Load the bundled ion implantation lookup table.

    Returns:
        Dictionary containing lookup table arrays:
            - z: Depth points (microns)
            - ImplantDopants: Dopant type indices (1=B, 2=P, 3=As)
            - ImplantDoses: Dose values (cm^-2)
            - ImplantEnergies: Energy values (keV)
            - AnnealTemps: Temperature values (C)
            - AnnealTimes: Time values (minutes)
            - AnnealOxidation: Oxidation type (1=inert, 2=oxide)
            - n: Doping concentration (cm^-3)
            - Xj: Junction depth (m)
            - Rs: Sheet resistance (Ohm/sq)
            - Nz: Effective doping (cm^-2)
            - Beta1, Beta2: Piezoresistive coefficients
    """
    global _LOOKUP_TABLE_CACHE
    if _LOOKUP_TABLE_CACHE is not None:
        return _LOOKUP_TABLE_CACHE

    with (
        resources.files("piezod.data").joinpath("ionImplantLookupTable.h5").open("rb") as f,
        h5py.File(f, "r") as hf,
    ):
        _LOOKUP_TABLE_CACHE = {key: np.array(hf[key]) for key in hf}

    return _LOOKUP_TABLE_CACHE


class CantileverImplantation(Cantilever):
    """Cantilever with ion-implanted piezoresistor doping.

    Lookup table data is automatically loaded from the bundled HDF5 file.

    Attributes:
        implantation_energy: Ion implantation energy (keV), range: 20-80
        implantation_dose: Ion implantation dose (cm^-2), range: 2e14-2e16
        annealing_temp: Annealing temperature (K), range: 1173-1373 (900-1100C)
        annealing_time: Annealing time (seconds), range: 900-7200 (15-120 min)
        annealing_type: Annealing environment ('inert' or 'oxide')
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
        doping_type: Literal["boron", "phosphorus", "arsenic"],
        annealing_time: float,
        annealing_temp: float,
        annealing_type: Literal["inert", "oxide"],
        implantation_energy: float,
        implantation_dose: float,
    ):
        """Initialize ion-implanted cantilever.

        Args:
            freq_min: Minimum target frequency (Hz)
            freq_max: Maximum target frequency (Hz)
            l: Cantilever length (m)
            w: Cantilever width (m)
            t: Cantilever thickness (m)
            l_pr_ratio: Piezoresistor length to cantilever length ratio (-)
            v_bridge: Bridge voltage (V)
            doping_type: Dopant type ('boron', 'phosphorus', or 'arsenic')
            annealing_time: Annealing time (seconds)
            annealing_temp: Annealing temperature (K)
            annealing_type: Annealing environment ('inert' or 'oxide')
            implantation_energy: Ion implantation energy (keV)
            implantation_dose: Ion implantation dose (cm^-2)
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

        # Set implantation-specific parameters
        self.implantation_energy = implantation_energy
        self.implantation_dose = implantation_dose
        self.annealing_type = annealing_type
        self.annealing_temp = annealing_temp
        self.annealing_time = annealing_time

        # Load bundled lookup table
        self._lookup_data = _load_lookup_table()

    def anneal_number(self) -> int:
        """Convert annealing type to lookup table index.

        Returns:
            1 for 'inert', 2 for 'oxide'

        Raises:
            ValueError: If annealing_type is not 'inert' or 'oxide'
        """
        if self.annealing_type == "inert":
            return 1
        elif self.annealing_type == "oxide":
            return 2
        else:
            raise ValueError(f"Unknown anneal condition: {self.annealing_type}. Must be 'inert' or 'oxide'.")

    def _interpolate_lookup(self, field_name: str) -> float:
        """Interpolate a single value from lookup table.

        This helper method handles the common interpolation pattern to avoid
        code duplication across multiple methods.

        Args:
            field_name: Name of the field to interpolate from lookup table

        Returns:
            Interpolated value at current cantilever parameters
        """
        # Convert temperatures from K to C and time from seconds to minutes
        anneal_temp_c = self.annealing_temp - 273
        anneal_time_min = self.annealing_time / 60

        result = interpn(
            (
                self._lookup_data["ImplantDopants"],
                self._lookup_data["ImplantDoses"],
                self._lookup_data["ImplantEnergies"],
                self._lookup_data["AnnealTemps"],
                self._lookup_data["AnnealTimes"],
                self._lookup_data["AnnealOxidation"],
            ),
            self._lookup_data[field_name],
            np.array(
                [
                    [
                        self.dopantNumber(),
                        self.implantation_dose,
                        self.implantation_energy,
                        anneal_temp_c,
                        anneal_time_min,
                        self.anneal_number(),
                    ]
                ]
            ),
            method="linear",
        )
        # interpn returns an array, extract the scalar value
        return float(result[0])

    def doping_profile(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Lookup the concentration profile from the lookup table.

        Returns:
            Tuple of (x, active_doping, total_doping):
                - x: Depth from surface (m)
                - active_doping: Active dopant concentration (cm^-3)
                - total_doping: Total dopant concentration (cm^-3)

        Note:
            Active = total unless the doping is higher than the solid solubility
            limit, which is generally not the case in the ion implantation data.
            Data beyond the device thickness is removed.
        """
        # Convert from microns to meters
        x = self._lookup_data["z"] * 1e-6  # 10nm spacing from 0 to 5um

        # Convert temperatures from K to C and time from seconds to minutes
        anneal_temp_c = self.annealing_temp - 273
        anneal_time_min = self.annealing_time / 60

        # For doping_profile, x is the first dimension, so we interpolate at each x point
        # Create query points: for each x, use the same dopant/dose/energy/temp/time/oxidation
        xi = np.column_stack(
            [
                x,  # First dimension is depth
                np.full(len(x), self.dopantNumber()),
                np.full(len(x), self.implantation_dose),
                np.full(len(x), self.implantation_energy),
                np.full(len(x), anneal_temp_c),
                np.full(len(x), anneal_time_min),
                np.full(len(x), self.anneal_number()),
            ]
        )

        # Interpolate doping concentration
        n = interpn(
            (
                x,
                self._lookup_data["ImplantDopants"],
                self._lookup_data["ImplantDoses"],
                self._lookup_data["ImplantEnergies"],
                self._lookup_data["AnnealTemps"],
                self._lookup_data["AnnealTimes"],
                self._lookup_data["AnnealOxidation"],
            ),
            self._lookup_data["n"],
            xi,
            method="linear",
        )

        # Remove data beyond the device thickness
        mask = x <= self.t
        x = x[mask]
        n = n[mask]

        # Active = total for ion implantation
        active_doping = n
        total_doping = n

        return x, active_doping, total_doping

    @property
    def junction_depth(self) -> float:
        """Junction depth from lookup table.

        Returns:
            Junction depth (m)
        """
        return self._interpolate_lookup("Xj")

    def sheet_resistance(self) -> float:
        """Sheet resistance from lookup table.

        Returns:
            Sheet resistance (Ohm/sq)
        """
        return self._interpolate_lookup("Rs")

    def Nz(self) -> float:
        """Effective doping concentration from lookup table.

        Note:
            Nz != Nz_total due to current crowding effects.

        Returns:
            Effective doping concentration (cm^-2)
        """
        return self._interpolate_lookup("Nz")

    def beta(self) -> float:
        """Piezoresistive coefficient from lookup table.

        Returns:
            Piezoresistive coefficient (Pa^-1)

        Note:
            Beta2 is in units of microns, so t is converted accordingly.
        """
        Beta1 = self._interpolate_lookup("Beta1")
        Beta2 = self._interpolate_lookup("Beta2")

        # Beta2 is in units of microns, so convert t
        return Beta1 - 2 / (self.t * 1e6) * Beta2

    @property
    def diffusion_length(self) -> float:
        """Calculate diffusion length from dopant diffusion.

        Based on dopant diffusion length, not Si diffusion length.
        From personal communications with L.K.J. Vandamme (Jan 2012).
        Note that most experimental data is from boron resistors.

        Returns:
            Diffusion length (cm)

        Raises:
            ValueError: If doping_type is not recognized
        """
        # Diffusion parameters for different dopants
        diffusion_params = {
            "arsenic": {"D0": 22.9, "Ea": 4.1},  # cm^2/s, eV
            "boron": {"D0": 0.76, "Ea": 3.46},  # cm^2/s, eV
            "phosphorus": {"D0": 3.85, "Ea": 3.66},  # cm^2/s, eV
        }

        if self.doping_type not in diffusion_params:
            raise ValueError(f"Unknown doping type: {self.doping_type}. Must be 'arsenic', 'boron', or 'phosphorus'.")

        params = diffusion_params[self.doping_type]
        D0 = params["D0"]
        Ea = params["Ea"]

        # Calculate diffusivity (cm^2/s)
        diffusivity = D0 * np.exp(-Ea / self.k_b_eV / self.annealing_temp)

        # Calculate diffusion length (cm)
        return np.sqrt(diffusivity * self.annealing_time)

    def alpha(self) -> float:
        """Calculate Hooge noise parameter from diffusion length.

        Calculated from sqrt(Dt) using empirical data compiled for the book.

        Returns:
            Hooge noise parameter (dimensionless)
        """
        return 2.469e-10 * self.diffusion_length**-0.598

    # ========= Optimization Methods ==========

    def doping_optimization_scaling(self) -> np.ndarray:
        """Get scaling factors for optimization variables.

        Returns:
            Array of scaling factors [time_scale, temp_scale, energy_scale, dose_scale]
        """
        return np.array([1e-2, 1e-2, 1e0, 1e-13])

    def doping_cantilever_from_state(self, x0: np.ndarray) -> "CantileverImplantation":
        """Update cantilever doping parameters from state vector.

        Args:
            x0: State vector with elements:
                [0-5]: Other cantilever parameters (not used here)
                [6]: annealing_time (seconds)
                [7]: annealing_temp (K)
                [8]: implantation_energy (keV)
                [9]: implantation_dose (cm^-2)

        Returns:
            Self (for method chaining)
        """
        self.annealing_time = x0[6]
        self.annealing_temp = x0[7]
        self.implantation_energy = x0[8]
        self.implantation_dose = x0[9]
        return self

    def doping_current_state(self) -> np.ndarray:
        """Get current doping parameters as state vector.

        Returns:
            Array of [annealing_time, annealing_temp, implantation_energy, implantation_dose]
        """
        return np.array(
            [
                self.annealing_time,
                self.annealing_temp,
                self.implantation_energy,
                self.implantation_dose,
            ]
        )

    def doping_optimization_bounds(self, parameter_constraints: dict | None = None) -> Tuple[np.ndarray, np.ndarray]:
        """Get optimization bounds for doping parameters.

        Args:
            parameter_constraints: Optional dict to override default bounds.
                Keys can be any of:
                - 'min_annealing_time', 'max_annealing_time' (seconds)
                - 'min_annealing_temp', 'max_annealing_temp' (K)
                - 'min_implantation_energy', 'max_implantation_energy' (keV)
                - 'min_implantation_dose', 'max_implantation_dose' (cm^-2)

        Returns:
            Tuple of (lower_bounds, upper_bounds) arrays
        """
        # Default bounds based on TSuprem data range
        bounds = {
            "min_annealing_time": 15 * 60,  # seconds
            "max_annealing_time": 120 * 60,
            "min_annealing_temp": 273 + 900,  # K
            "max_annealing_temp": 273 + 1100,
            "min_implantation_energy": 20,  # keV
            "max_implantation_energy": 80,
            "min_implantation_dose": 2e14,  # cm^-2
            "max_implantation_dose": 2e16,
        }

        # Override with user-provided constraints
        if parameter_constraints is not None:
            bounds.update(parameter_constraints)

        lb = np.array(
            [
                bounds["min_annealing_time"],
                bounds["min_annealing_temp"],
                bounds["min_implantation_energy"],
                bounds["min_implantation_dose"],
            ]
        )

        ub = np.array(
            [
                bounds["max_annealing_time"],
                bounds["max_annealing_temp"],
                bounds["max_implantation_energy"],
                bounds["max_implantation_dose"],
            ]
        )

        return lb, ub

    def doping_initial_conditions_random(self, parameter_constraints: dict | None = None) -> np.ndarray:
        """Generate random initial conditions within optimization bounds.

        Args:
            parameter_constraints: Optional dict to override default bounds

        Returns:
            Random state vector within bounds
        """
        lb, ub = self.doping_optimization_bounds(parameter_constraints)
        return lb + np.random.rand(4) * (ub - lb)
