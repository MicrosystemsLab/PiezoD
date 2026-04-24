"""Cantilever with ion-implanted doping profile.

Model an ion implanted cantilever using lookup tables from process simulators.
Two lookup table sources are supported:

- "tsuprem4": TSUPREM-4 lookup table (ionImplantLookupTable_tsuprem.h5)
  - Dopants: B, P, As
  - Anneal types: "inert" (1), "oxide" (2)
  - Energy: 20-80 keV, Dose: 2e14-2e16 cm^-2
  - Temperature: 900-1100C, Time: 15-120 min

- "dopedealer": DopeDealer lookup table (ionImplantLookupTable_dopedealer.h5)
  - Dopants: B, P, As
  - Anneal types: "inert" (1), "dry_o2" (2), "wet_o2" (3)
  - Energy: 10-120 keV, Dose: 1e13-8e16 cm^-2
  - Temperature: 850-1100C, Time: 15-120 min

For all conditions, a 250A protection oxide layer is grown before the ion
implantation.

After the implantation step (TSUPREM-4 table):
- 'oxide' -> strip protection oxide, regrow 1500A oxide, inert anneal
- 'inert' -> leave protection oxide, inert anneal, strip oxide

Lookup tables are linearly interpolated. Splines can result in negative values
(e.g. sheet resistance). Gradient/Hessian discontinuities don't seem to pose a
problem for optimization.
"""

from importlib import resources
from typing import Literal, Tuple

import h5py
import numpy as np
from scipy.interpolate import interpn

from piezod.cantilever import Cantilever

# Module-level cache for lookup table data, keyed by source name
_LOOKUP_TABLE_CACHE: dict[str, dict] = {}

# Valid annealing types per lookup source
_ANNEAL_TYPES: dict[str, dict[str, int]] = {
    "tsuprem4": {"inert": 1, "oxide": 2},
    "dopedealer": {"inert": 1, "dry_o2": 2, "wet_o2": 3},
}

# Default optimization bounds per lookup source
_DEFAULT_BOUNDS: dict[str, dict[str, float]] = {
    "tsuprem4": {
        "min_annealing_time": 15 * 60,
        "max_annealing_time": 120 * 60,
        "min_annealing_temp": 273 + 900,
        "max_annealing_temp": 273 + 1100,
        "min_implantation_energy": 20,
        "max_implantation_energy": 80,
        "min_implantation_dose": 2e14,
        "max_implantation_dose": 2e16,
    },
    "dopedealer": {
        "min_annealing_time": 15 * 60,
        "max_annealing_time": 120 * 60,
        "min_annealing_temp": 273 + 850,
        "max_annealing_temp": 273 + 1100,
        "min_implantation_energy": 10,
        "max_implantation_energy": 120,
        "min_implantation_dose": 1e13,
        "max_implantation_dose": 8e16,
    },
}

# HDF5 filenames per lookup source
_LOOKUP_FILENAMES: dict[str, str] = {
    "tsuprem4": "ionImplantLookupTable_tsuprem.h5",
    "dopedealer": "ionImplantLookupTable_dopedealer.h5",
}


def _load_lookup_table(source: str) -> dict:
    """Load a bundled ion implantation lookup table.

    Args:
        source: Lookup table source ("tsuprem4" or "dopedealer")

    Returns:
        Dictionary containing lookup table arrays.

    Raises:
        ValueError: If source is not recognized.
        FileNotFoundError: If the lookup table file is not bundled.
    """
    if source in _LOOKUP_TABLE_CACHE:
        return _LOOKUP_TABLE_CACHE[source]

    if source not in _LOOKUP_FILENAMES:
        raise ValueError(
            f"Unknown lookup source: {source!r}. "
            f"Must be one of: {sorted(_LOOKUP_FILENAMES.keys())}"
        )

    filename = _LOOKUP_FILENAMES[source]

    with (
        resources.files("piezod.data").joinpath(filename).open("rb") as f,
        h5py.File(f, "r") as hf,
    ):
        _LOOKUP_TABLE_CACHE[source] = {key: np.array(hf[key]) for key in hf}

    return _LOOKUP_TABLE_CACHE[source]


class CantileverImplantation(Cantilever):
    """Cantilever with ion-implanted piezoresistor doping.

    Lookup table data is automatically loaded from the bundled HDF5 file.

    Args:
        lookup_source: Which process simulator lookup table to use.
            "tsuprem4" uses the TSUPREM-4 table (default).
            "dopedealer" uses the DopeDealer table.

    Attributes:
        implantation_energy: Ion implantation energy (keV)
        implantation_dose: Ion implantation dose (cm^-2)
        annealing_temp: Annealing temperature (K)
        annealing_time: Annealing time (seconds)
        annealing_type: Annealing environment (source-dependent, see module docstring)
        lookup_source: Which lookup table is in use ("tsuprem4" or "dopedealer")
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
        annealing_type: str,
        implantation_energy: float,
        implantation_dose: float,
        lookup_source: Literal["tsuprem4", "dopedealer"] = "tsuprem4",
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
            annealing_type: Annealing environment. For tsuprem4: 'inert' or
                'oxide'. For dopedealer: 'inert', 'dry_o2', or 'wet_o2'.
            implantation_energy: Ion implantation energy (keV)
            implantation_dose: Ion implantation dose (cm^-2)
            lookup_source: Lookup table source ('tsuprem4' or 'dopedealer')
        """
        if lookup_source not in _LOOKUP_FILENAMES:
            raise ValueError(
                f"Unknown lookup_source: {lookup_source!r}. "
                f"Must be one of: {sorted(_LOOKUP_FILENAMES.keys())}"
            )

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
        self.lookup_source = lookup_source

        # Load bundled lookup table
        self._lookup_data = _load_lookup_table(lookup_source)

        # Resolve the ambient axis key (differs between tables)
        if "AnnealAmbient" in self._lookup_data:
            self._ambient_key = "AnnealAmbient"
        else:
            self._ambient_key = "AnnealOxidation"

    def anneal_number(self) -> int:
        """Convert annealing type to lookup table index.

        Returns:
            Integer index for the annealing type in the active lookup table.

        Raises:
            ValueError: If annealing_type is not valid for the active table.
        """
        anneal_map = _ANNEAL_TYPES[self.lookup_source]
        if self.annealing_type in anneal_map:
            return anneal_map[self.annealing_type]
        valid = sorted(anneal_map.keys())
        raise ValueError(
            f"Unknown anneal condition: {self.annealing_type!r} "
            f"for source {self.lookup_source!r}. Must be one of: {valid}"
        )

    def _interpolate_lookup(self, field_name: str) -> float:
        """Interpolate a single value from lookup table.

        Args:
            field_name: Name of the field to interpolate from lookup table

        Returns:
            Interpolated value at current cantilever parameters
        """
        anneal_temp_c = self.annealing_temp - 273
        anneal_time_min = self.annealing_time / 60

        result = interpn(
            (
                self._lookup_data["ImplantDopants"],
                self._lookup_data["ImplantDoses"],
                self._lookup_data["ImplantEnergies"],
                self._lookup_data["AnnealTemps"],
                self._lookup_data["AnnealTimes"],
                self._lookup_data[self._ambient_key],
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
        x = self._lookup_data["z"] * 1e-6  # Convert from microns to meters

        anneal_temp_c = self.annealing_temp - 273
        anneal_time_min = self.annealing_time / 60

        xi = np.column_stack(
            [
                x,
                np.full(len(x), self.dopantNumber()),
                np.full(len(x), self.implantation_dose),
                np.full(len(x), self.implantation_energy),
                np.full(len(x), anneal_temp_c),
                np.full(len(x), anneal_time_min),
                np.full(len(x), self.anneal_number()),
            ]
        )

        n = interpn(
            (
                x,
                self._lookup_data["ImplantDopants"],
                self._lookup_data["ImplantDoses"],
                self._lookup_data["ImplantEnergies"],
                self._lookup_data["AnnealTemps"],
                self._lookup_data["AnnealTimes"],
                self._lookup_data[self._ambient_key],
            ),
            self._lookup_data["n"],
            xi,
            method="linear",
        )

        mask = x <= self.t
        x = x[mask]
        n = n[mask]

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
        """Mobility-weighted effective carrier sheet density from lookup table.

        Harmon current-crowding formula:
            Nz = (integral mu*N dz)^2 / integral mu^2*N dz
        Nz <= Nz_total (the raw retained active dose); the difference arises
        because higher-doping regions carry a disproportionate share of
        current, reducing the effective carrier count for 1/f noise.

        Returns:
            Effective carrier sheet density (1/m^2). Consistent with
            cantilever.number_of_carriers = Nz * resistor_area[m^2].
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
        diffusion_params = {
            "arsenic": {"D0": 22.9, "Ea": 4.1},  # cm^2/s, eV
            "boron": {"D0": 0.76, "Ea": 3.46},  # cm^2/s, eV
            "phosphorus": {"D0": 3.85, "Ea": 3.66},  # cm^2/s, eV
        }

        if self.doping_type not in diffusion_params:
            raise ValueError(
                f"Unknown doping type: {self.doping_type}. "
                f"Must be 'arsenic', 'boron', or 'phosphorus'."
            )

        params = diffusion_params[self.doping_type]
        D0 = params["D0"]
        Ea = params["Ea"]

        diffusivity = D0 * np.exp(-Ea / self.k_b_eV / self.annealing_temp)

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

        Default bounds depend on the active lookup table source.

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
        bounds = dict(_DEFAULT_BOUNDS[self.lookup_source])

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
