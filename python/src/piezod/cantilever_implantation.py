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
  - Anneal types: "inert" (1), "dry_o2" (2)
  - Energy: 10-120 keV, Dose: 1e13-5e16 cm^-2
  - Temperature: 900-1100C, Time: 15-150 min

For all conditions, a 250A protection oxide layer is grown before the ion
implantation.

After the implantation step (TSUPREM-4 table):
- 'oxide' -> strip protection oxide, regrow 1500A oxide, inert anneal
- 'inert' -> leave protection oxide, inert anneal, strip oxide

Lookup tables are linearly interpolated. Splines can result in negative values
(e.g. sheet resistance). Gradient/Hessian discontinuities don't seem to pose a
problem for optimization.
"""

from __future__ import annotations

import copy
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from enum import Enum
from importlib import resources
from typing import Literal, Tuple

import h5py
import numpy as np
from scipy.interpolate import interpn
from scipy.optimize import OptimizeResult, minimize

from piezod.cantilever import Cantilever
from piezod.optimization.state import StateVar

# Module-level cache for lookup table data, keyed by source name
_LOOKUP_TABLE_CACHE: dict[str, dict] = {}

# Valid annealing types per lookup source
_ANNEAL_TYPES: dict[str, dict[str, int]] = {
    "tsuprem4": {"inert": 1, "oxide": 2},
    "dopedealer": {"inert": 1, "dry_o2": 2},
}

# Default optimization bounds per lookup source
_DEFAULT_BOUNDS: dict[str, dict[str, float]] = {
    "tsuprem4": {
        "min_annealing_time": 15 * 60,
        "max_annealing_time": 120 * 60,
        "min_annealing_temp": 273.15 + 900,
        "max_annealing_temp": 273.15 + 1100,
        "min_implantation_energy": 20,
        "max_implantation_energy": 80,
        "min_implantation_dose": 2e14,
        "max_implantation_dose": 2e16,
    },
    "dopedealer": {
        "min_annealing_time": 15 * 60,
        "max_annealing_time": 150 * 60,
        "min_annealing_temp": 273.15 + 900,
        "max_annealing_temp": 273.15 + 1100,
        "min_implantation_energy": 10,
        "max_implantation_energy": 120,
        "min_implantation_dose": 1e13,
        "max_implantation_dose": 5e16,
    },
}

# HDF5 filenames per lookup source
_LOOKUP_FILENAMES: dict[str, str] = {
    "tsuprem4": "ionImplantLookupTable_tsuprem.h5",
    "dopedealer": "ionImplantLookupTable_dopedealer.h5",
}

# The TSUPREM-4 lookup tables were generated with a uniform substrate
# background of ~1.36e15 cm^-3 (same dopant species as the implant) baked
# into every simulation, so the stored `n` is dopant + 1.36e15. The DopeDealer
# tables have no baked-in background. We subtract the constant at load time
# so both sources expose the dopant species alone, matching the
# `doping_profile()` convention. Verified empirically: tail values of `n`
# cluster at ~1.36e15 across most TSUPREM parameter combinations.
_TSUPREM_BAKED_BACKGROUND_CM3: float = 1.36e15

_DOPING_STATE_NAMES: tuple[str, ...] = (
    "annealing_time",
    "annealing_temp",
    "implantation_energy",
    "implantation_dose",
)

# Boltzmann constant in eV/K, matching Cantilever.k_b_eV. Local copy avoids an
# import cycle when this module is imported standalone.
_K_B_EV: float = 8.617343e-5

# Dopant diffusion parameters used by `diffusion_length_cm` (D0 in cm^2/s, Ea in eV).
# From personal communications with L.K.J. Vandamme (Jan 2012). Most experimental
# data is from boron resistors.
_DIFFUSION_PARAMS: dict[str, dict[str, float]] = {
    "arsenic": {"D0": 22.9, "Ea": 4.1},
    "boron": {"D0": 0.76, "Ea": 3.46},
    "phosphorus": {"D0": 3.85, "Ea": 3.66},
}


def diffusion_length_cm(doping_type: str, annealing_temp_K: float, annealing_time_s: float) -> float:
    """Dopant diffusion length sqrt(D*t) in cm.

    Args:
        doping_type: One of 'boron', 'phosphorus', 'arsenic'.
        annealing_temp_K: Annealing temperature in Kelvin.
        annealing_time_s: Annealing time in seconds.

    Returns:
        Diffusion length in cm.

    Raises:
        ValueError: If doping_type is not recognized.
    """
    if doping_type not in _DIFFUSION_PARAMS:
        raise ValueError(f"Unknown doping_type: {doping_type!r}. Must be 'boron', 'phosphorus', or 'arsenic'.")
    params = _DIFFUSION_PARAMS[doping_type]
    diffusivity = params["D0"] * np.exp(-params["Ea"] / _K_B_EV / annealing_temp_K)
    return float(np.sqrt(diffusivity * annealing_time_s))


def hooge_alpha_from_anneal(
    doping_type: str,
    annealing_temp_K: float,
    annealing_time_s: float,
) -> float:
    """Hooge 1/f noise parameter from anneal conditions.

    Empirical fit `2.469e-10 * sqrt(D*t)^-0.598` against the dopant diffusion
    length (Vandamme, personal communication; data primarily from boron PRs).

    Args:
        doping_type: One of 'boron', 'phosphorus', 'arsenic'.
        annealing_temp_K: Annealing temperature in Kelvin.
        annealing_time_s: Annealing time in seconds.

    Returns:
        Hooge alpha (dimensionless).
    """
    return 2.469e-10 * diffusion_length_cm(doping_type, annealing_temp_K, annealing_time_s) ** -0.598


@dataclass(frozen=True)
class DopingProcessMetrics:
    """Ion implantation process metrics used for optimization."""

    alpha_h: float
    nz: float
    beta1: float
    beta2_um: float
    beta: float
    sheet_resistance: float
    junction_depth_m: float
    peak_concentration_cm3: float
    retained_dose_cm2: float


class DopingMetric(str, Enum):
    """Names of process metrics that can be hard-constrained during optimization.

    Each value matches the corresponding attribute on DopingProcessMetrics so
    the metric can be looked up generically.
    """

    ALPHA_H = "alpha_h"
    NZ = "nz"
    BETA1 = "beta1"
    BETA2_UM = "beta2_um"
    BETA = "beta"
    SHEET_RESISTANCE = "sheet_resistance"
    JUNCTION_DEPTH_M = "junction_depth_m"
    PEAK_CONCENTRATION_CM3 = "peak_concentration_cm3"
    RETAINED_DOSE_CM2 = "retained_dose_cm2"


@dataclass(frozen=True)
class MetricConstraint:
    """Hard inequality constraint on a derived doping process metric.

    At least one of `minimum` or `maximum` must be provided. Bounds are in the
    units exposed by DopingProcessMetrics: sheet_resistance in ohm/sq,
    junction_depth_m in meters, beta2_um in microns, retained_dose_cm2 in
    cm^-2, peak_concentration_cm3 in cm^-3.
    """

    metric: DopingMetric
    minimum: float | None = None
    maximum: float | None = None


@dataclass(frozen=True)
class DopingOptimizationResult:
    """Result returned by doping process optimization."""

    optimized: CantileverImplantation
    metrics: DopingProcessMetrics
    objective_value: float
    state: np.ndarray
    scipy_result: OptimizeResult
    all_results: tuple[OptimizeResult, ...]
    boundary_flags: dict[str, bool]


_METRIC_CONSTRAINT_TOLERANCE = 1e-6


def _default_doping_objective(metrics: DopingProcessMetrics) -> float:
    """Evaluate the default noise and sensitivity process objective."""
    beta_magnitude = abs(metrics.beta)
    if metrics.alpha_h <= 0 or metrics.nz <= 0 or beta_magnitude <= 0:
        return float("inf")

    return float(np.log(metrics.alpha_h) - np.log(metrics.nz) - 2 * np.log(beta_magnitude))


def _doping_boundary_flags(
    state: np.ndarray,
    lower_bounds: np.ndarray,
    upper_bounds: np.ndarray,
) -> dict[str, bool]:
    """Return flags indicating whether the public state is on each bound."""
    flags: dict[str, bool] = {}
    for index, name in enumerate(_DOPING_STATE_NAMES):
        flags[f"{name}_at_min"] = bool(np.isclose(state[index], lower_bounds[index]))
        flags[f"{name}_at_max"] = bool(np.isclose(state[index], upper_bounds[index]))
    return flags


def _validate_metric_constraints(
    constraints: Sequence[MetricConstraint],
) -> tuple[MetricConstraint, ...]:
    """Validate user-supplied metric constraints.

    Each entry must be a MetricConstraint with a DopingMetric metric, at least
    one finite bound, and minimum <= maximum if both are present.

    Args:
        constraints: Sequence of MetricConstraint instances.

    Returns:
        Immutable tuple of validated constraints.

    Raises:
        ValueError: If any constraint is malformed.
    """
    validated: list[MetricConstraint] = []
    for index, constraint in enumerate(constraints):
        if not isinstance(constraint, MetricConstraint):
            raise ValueError(
                f"metric_constraints[{index}] must be a MetricConstraint, got {type(constraint).__name__}."
            )
        if not isinstance(constraint.metric, DopingMetric):
            raise ValueError(
                f"metric_constraints[{index}].metric must be a DopingMetric enum, "
                f"got {type(constraint.metric).__name__}."
            )

        minimum = constraint.minimum
        maximum = constraint.maximum
        if minimum is None and maximum is None:
            raise ValueError(
                f"metric_constraints[{index}] for {constraint.metric.value!r} must specify "
                "at least one of minimum or maximum."
            )
        if minimum is not None and not np.isfinite(minimum):
            raise ValueError(f"metric_constraints[{index}].minimum must be finite, got {minimum!r}.")
        if maximum is not None and not np.isfinite(maximum):
            raise ValueError(f"metric_constraints[{index}].maximum must be finite, got {maximum!r}.")
        if minimum is not None and maximum is not None and minimum > maximum:
            raise ValueError(
                f"metric_constraints[{index}] for {constraint.metric.value!r}: "
                f"minimum ({minimum}) must be <= maximum ({maximum})."
            )

        validated.append(constraint)

    return tuple(validated)


def _metric_value(metrics: DopingProcessMetrics, metric: DopingMetric) -> float:
    """Extract a metric value from DopingProcessMetrics by enum."""
    return float(getattr(metrics, metric.value))


def _metric_constraint_residuals(
    metrics: DopingProcessMetrics,
    constraints: Sequence[MetricConstraint],
) -> np.ndarray:
    """Build the SLSQP inequality residual vector.

    Each constraint contributes one residual per provided bound, scaled by
    max(abs(limit), 1.0) for numerical conditioning. A residual >= 0 means the
    bound is satisfied.

    Args:
        metrics: Process metrics evaluated at the candidate state.
        constraints: Validated metric constraints.

    Returns:
        Residual array (length equals total number of bounds across all constraints).
    """
    residuals: list[float] = []
    for constraint in constraints:
        value = _metric_value(metrics, constraint.metric)
        if constraint.minimum is not None:
            scale = max(abs(constraint.minimum), 1.0)
            residuals.append((value - constraint.minimum) / scale)
        if constraint.maximum is not None:
            scale = max(abs(constraint.maximum), 1.0)
            residuals.append((constraint.maximum - value) / scale)
    return np.asarray(residuals, dtype=float)


def _build_slsqp_constraints(
    constraints: Sequence[MetricConstraint],
    state_to_metrics: Callable[[np.ndarray], DopingProcessMetrics],
) -> list[dict]:
    """Build a single vector-valued SLSQP inequality constraint.

    Using a single dict means doping_process_metrics() is computed once per
    optimizer state, not once per bound, which matters because each metric
    evaluation interpolates the lookup table.

    Args:
        constraints: Validated metric constraints.
        state_to_metrics: Callable that maps optimizer state to metrics.

    Returns:
        List with one SLSQP constraint dict.
    """

    def fun(optimizer_state: np.ndarray) -> np.ndarray:
        metrics = state_to_metrics(optimizer_state)
        return _metric_constraint_residuals(metrics, constraints)

    return [{"type": "ineq", "fun": fun}]


def _constraints_satisfied(
    metrics: DopingProcessMetrics,
    constraints: Sequence[MetricConstraint],
    tol: float = _METRIC_CONSTRAINT_TOLERANCE,
) -> bool:
    """Check whether all constraint residuals are nonnegative within tol."""
    if not constraints:
        return True
    residuals = _metric_constraint_residuals(metrics, constraints)
    return bool(np.all(residuals >= -tol))


def _impute_nan_harmonic(arr: np.ndarray, max_iter: int = 500, tol: float = 1e-12) -> np.ndarray:
    """Fill NaN cells with the harmonic (discrete-Laplace) extension of valid cells.

    Iterates Gauss-Seidel "cell = mean of face neighbors" on NaN cells until
    convergence, with valid cells acting as Dirichlet boundary conditions. In
    1D this converges to exact linear interpolation between the bounding
    valid cells; in higher dim it gives the curvature-minimizing fill with
    C^0 continuity at hole boundaries (no kinks). `np.roll` provides
    reflecting boundaries; this is acceptable as long as holes do not touch
    array edges.

    Args:
        arr: Array containing NaN cells. Returned unchanged when no NaN.
        max_iter: Maximum Gauss-Seidel iterations.
        tol: Convergence tolerance on the maximum absolute update.

    Returns:
        Array with NaN cells filled. The input is not modified.

    Raises:
        ValueError: If the array is entirely NaN (no boundary data).
    """
    mask = np.isnan(arr)
    if not mask.any():
        return arr
    if mask.all():
        raise ValueError("Cannot inpaint: array is entirely NaN.")

    out = arr.astype(np.float64, copy=True)
    out[mask] = float(np.nanmean(arr))

    ndim = out.ndim
    for _ in range(max_iter):
        s = np.zeros_like(out)
        for ax in range(ndim):
            s += np.roll(out, 1, axis=ax) + np.roll(out, -1, axis=ax)
        s /= 2 * ndim
        delta = np.where(mask, s - out, 0.0)
        out = out + delta
        if float(np.max(np.abs(delta))) < tol:
            break
    return out


def _impute_lookup_holes(data: dict) -> None:
    """Inpaint NaN cells in metric arrays per (dopant, ambient) sub-cube.

    LUT generators occasionally fail at certain (dopant, dose, energy, temp,
    ambient) combinations, leaving NaN scalars in the metric arrays. We fill
    those holes with the harmonic extension of valid neighbors restricted to
    the same dopant species and same anneal ambient -- never borrowing
    physics across species or anneal types. Without this, downstream
    `interpn(method="linear")` returns NaN whenever the enclosing hypercube
    touches a hole, even via a 0-weight corner, because IEEE-754 arithmetic
    gives `0 * NaN = NaN`. Modifies `data` in place.
    """
    metric_keys = ("Beta1", "Beta2", "Nz", "Nz_total", "Rs", "Xj")
    n_dopants = len(data["ImplantDopants"])
    ambient_key = "AnnealAmbient" if "AnnealAmbient" in data else "AnnealOxidation"
    n_ambients = len(data[ambient_key])

    for key in metric_keys:
        if key not in data:
            continue
        arr = data[key]
        if not np.issubdtype(arr.dtype, np.floating) or not np.isnan(arr).any():
            continue
        arr = arr.astype(np.float64, copy=True)
        for dopant_idx in range(n_dopants):
            for ambient_idx in range(n_ambients):
                sub = arr[dopant_idx, ..., ambient_idx]
                if not np.isnan(sub).any():
                    continue
                arr[dopant_idx, ..., ambient_idx] = _impute_nan_harmonic(sub)
        data[key] = arr


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
        raise ValueError(f"Unknown lookup source: {source!r}. Must be one of: {sorted(_LOOKUP_FILENAMES.keys())}")

    filename = _LOOKUP_FILENAMES[source]

    with (
        resources.files("piezod.data").joinpath(filename).open("rb") as f,
        h5py.File(f, "r") as hf,
    ):
        data = {key: np.array(hf[key]) for key in hf}

    if source == "tsuprem4":
        data["n"] = np.maximum(data["n"].astype(np.float64) - _TSUPREM_BAKED_BACKGROUND_CM3, 0.0)

    _impute_lookup_holes(data)

    _LOOKUP_TABLE_CACHE[source] = data
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
                'oxide'. For dopedealer: 'inert' or 'dry_o2'.
            implantation_energy: Ion implantation energy (keV)
            implantation_dose: Ion implantation dose (cm^-2)
            lookup_source: Lookup table source ('tsuprem4' or 'dopedealer')
        """
        if lookup_source not in _LOOKUP_FILENAMES:
            raise ValueError(
                f"Unknown lookup_source: {lookup_source!r}. Must be one of: {sorted(_LOOKUP_FILENAMES.keys())}"
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
        anneal_temp_c = self.annealing_temp - 273.15
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
        """Lookup the doping profile from the lookup table.

        Returns:
            Tuple of (x, active_doping, total_doping):
                - x: Depth from surface (m).
                - active_doping: Net active carriers of the resistor type
                  (cm^-3) = ``max(0, dopant_species - substrate_background_cm3)``.
                  Goes to zero below the junction.
                - total_doping: Total concentration ``max(dopant_species,
                  substrate_background_cm3)`` (cm^-3); floors at the substrate
                  background so the profile is directly plottable.

        Note:
            For the TSUPREM-4 lookup, the simulator's baked-in 1.36e15 cm^-3
            substrate is removed at table-load time so the underlying dopant
            species curve matches the DopeDealer convention before the
            substrate floor is reapplied. The lookup-table-derived metrics
            (Rs, Nz, Xj, Beta1, Beta2) still embed the simulator's substrate
            assumption and are unaffected by ``substrate_background_cm3``; for
            metrics that track a custom background, pass this profile to
            :class:`PiezoresistorFromProfile`. Data beyond the device
            thickness is removed.
        """
        x = self._lookup_data["z"] * 1e-6  # Convert from microns to meters

        anneal_temp_c = self.annealing_temp - 273.15
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

        total_doping = np.maximum(n, self.substrate_background_cm3)
        active_doping = np.maximum(0.0, n - self.substrate_background_cm3)

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
        """Dopant diffusion length sqrt(D*t) in cm."""
        return diffusion_length_cm(self.doping_type, self.annealing_temp, self.annealing_time)

    def alpha(self) -> float:
        """Hooge 1/f noise parameter (dimensionless), from anneal conditions."""
        return hooge_alpha_from_anneal(self.doping_type, self.annealing_temp, self.annealing_time)

    def doping_process_metrics(self, device_thickness_m: float | None = None) -> DopingProcessMetrics:
        """Return process-only doping metrics from the implantation lookup table.

        Args:
            device_thickness_m: Optional device thickness used for the combined
                beta calculation. If omitted, the cantilever thickness is used.

        Returns:
            Process metrics for the current implant and anneal conditions.
        """
        if device_thickness_m is None:
            device_thickness_m = self.t

        device_thickness_um = device_thickness_m * 1e6
        beta1 = self._interpolate_lookup("Beta1")
        beta2_um = self._interpolate_lookup("Beta2")
        beta = beta1 - 2 * beta2_um / device_thickness_um

        depth_m, _active_doping, total_doping = self.doping_profile()
        depth_cm = depth_m * 100

        return DopingProcessMetrics(
            alpha_h=float(self.alpha()),
            nz=float(self.Nz()),
            beta1=float(beta1),
            beta2_um=float(beta2_um),
            beta=float(beta),
            sheet_resistance=float(self.sheet_resistance()),
            junction_depth_m=float(self.junction_depth),
            peak_concentration_cm3=float(np.max(total_doping)),
            retained_dose_cm2=float(np.trapezoid(total_doping, depth_cm)),
        )

    # ========= Optimization Methods ==========

    def _copy_with_doping_state(self, state: np.ndarray) -> "CantileverImplantation":
        """Return a shallow copy with process variables replaced by state."""
        candidate = copy.copy(self)
        candidate.annealing_time = float(state[0])
        candidate.annealing_temp = float(state[1])
        candidate.implantation_energy = float(state[2])
        candidate.implantation_dose = float(state[3])
        return candidate

    @staticmethod
    def _validate_doping_state(state: np.ndarray) -> np.ndarray:
        """Validate and normalize a process-only state vector."""
        state = np.asarray(state, dtype=float)
        if state.shape != (4,):
            raise ValueError(
                "Doping state must contain four values: "
                "annealing_time_s, annealing_temp_K, implantation_energy_keV, implantation_dose_cm2."
            )
        return state

    @staticmethod
    def _optimizer_bounds(
        lower_bounds: np.ndarray,
        upper_bounds: np.ndarray,
        log_dose: bool,
    ) -> tuple[list[tuple[float, float]], np.ndarray, np.ndarray]:
        """Return optimizer bounds and transformed lower and upper arrays."""
        optimizer_lower_bounds = lower_bounds.astype(float).copy()
        optimizer_upper_bounds = upper_bounds.astype(float).copy()

        if log_dose:
            if lower_bounds[3] <= 0 or upper_bounds[3] <= 0:
                raise ValueError("Dose bounds must be positive when log_dose=True.")
            optimizer_lower_bounds[3] = np.log10(lower_bounds[3])
            optimizer_upper_bounds[3] = np.log10(upper_bounds[3])

        optimizer_bounds = list(zip(optimizer_lower_bounds, optimizer_upper_bounds, strict=True))
        return optimizer_bounds, optimizer_lower_bounds, optimizer_upper_bounds

    @staticmethod
    def _to_optimizer_state(state: np.ndarray, log_dose: bool) -> np.ndarray:
        """Convert a physical process state to optimizer coordinates."""
        optimizer_state = state.astype(float).copy()
        if log_dose:
            if optimizer_state[3] <= 0:
                raise ValueError("Initial dose must be positive when log_dose=True.")
            optimizer_state[3] = np.log10(optimizer_state[3])
        return optimizer_state

    @staticmethod
    def _to_physical_state(
        optimizer_state: np.ndarray,
        log_dose: bool,
        dose_bounds: tuple[float, float] | None = None,
    ) -> np.ndarray:
        """Convert optimizer coordinates to physical process state."""
        state = np.asarray(optimizer_state, dtype=float).copy()
        if log_dose:
            dose = 10 ** state[3]
            if dose_bounds is not None:
                dose_min, dose_max = dose_bounds
                if np.isclose(dose, dose_min, rtol=1e-14, atol=0.0):
                    dose = dose_min
                elif np.isclose(dose, dose_max, rtol=1e-14, atol=0.0):
                    dose = dose_max
            state[3] = dose
        return state

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

    def optimization_state_vars(self) -> tuple[StateVar, ...]:
        """Declarative state spec for joint geometry + implant optimization.

        Nine state variables: cantilever length, width, thickness,
        piezoresistor length ratio, bridge bias voltage, annealing time,
        annealing temperature, implant energy, and implant dose
        (log-scale). Default doping bounds come from
        ``_DEFAULT_BOUNDS[self.lookup_source]`` so the bounds match
        whichever lookup table the cantilever was constructed with.
        """
        bounds = _DEFAULT_BOUNDS[self.lookup_source]
        return (
            StateVar("l", 1e5, 10e-6, 3e-3),
            StateVar("w", 1e7, 2e-6, 100e-6),
            StateVar("t", 1e8, 1e-6, 100e-6),
            StateVar("l_pr_ratio", 1e2, 0.01, 0.99),
            StateVar("v_bridge", 1e1, 0.1, 10.0),
            StateVar(
                "annealing_time",
                1e-2,
                bounds["min_annealing_time"],
                bounds["max_annealing_time"],
            ),
            StateVar(
                "annealing_temp",
                1e-2,
                bounds["min_annealing_temp"],
                bounds["max_annealing_temp"],
            ),
            StateVar(
                "implantation_energy",
                1.0,
                bounds["min_implantation_energy"],
                bounds["max_implantation_energy"],
            ),
            StateVar(
                "implantation_dose",
                1.0,
                bounds["min_implantation_dose"],
                bounds["max_implantation_dose"],
                log_scale=True,
            ),
        )

    def optimize_doping_for_hooge_noise(
        self,
        objective: Callable[[DopingProcessMetrics], float] | None = None,
        *,
        device_thickness_m: float | None = None,
        parameter_constraints: dict[str, float] | None = None,
        metric_constraints: Sequence[MetricConstraint] | None = None,
        initial_states: Sequence[np.ndarray] | None = None,
        n_random_starts: int = 0,
        method: str | None = None,
        log_dose: bool = True,
        random_seed: int | None = None,
    ) -> DopingOptimizationResult:
        """Optimize the implant and anneal process variables for Hooge-noise-limited resolution.

        The default objective is the Hooge-noise-limited force-resolution figure
        of merit, log(alpha_h) - log(Nz) - 2*log(|beta|). Pass a custom
        `objective` to optimize a different criterion against the same process
        variables.

        Mirrors the MATLAB PiezoD fmincon pattern: direct parameter bounds
        (parameter_constraints) plus optional nonlinear constraints on derived
        metrics (metric_constraints).

        Args:
            objective: Maps DopingProcessMetrics to a scalar to minimize.
                Defaults to the Hooge-noise-limited resolution figure of merit.
            device_thickness_m: Optional thickness used for the combined beta
                metric and for any beta-based metric constraints. If omitted,
                the cantilever thickness is used. Threaded into both the
                objective and the constraint metric evaluation so they agree.
            parameter_constraints: Direct bounds on the four process variables
                (annealing time/temp, implant energy/dose). Used as the
                optimizer's box bounds. Temperatures are Kelvin. See
                doping_optimization_bounds for keys.
            metric_constraints: Hard inequality constraints on derived process
                metrics from DopingProcessMetrics. Each MetricConstraint can
                set a minimum, maximum, or both. Units: sheet_resistance in
                ohm/sq, junction_depth_m in meters, beta2_um in microns,
                retained_dose_cm2 in cm^-2, peak_concentration_cm3 in cm^-3.
                When provided, optimization uses SLSQP by default; passing
                method="L-BFGS-B" raises ValueError because L-BFGS-B cannot
                enforce nonlinear constraints.
            initial_states: Physical process state vectors to use as starts.
            n_random_starts: Number of additional random physical starts.
            method: SciPy minimize method. If None, defaults to "L-BFGS-B"
                without metric_constraints and "SLSQP" with them.
            log_dose: If true, optimize log10(dose) internally.
            random_seed: Optional seed for random starts.

        Returns:
            Optimization result with the optimized cantilever copy and metrics.

        Raises:
            ValueError: If metric_constraints are malformed, or method
                "L-BFGS-B" is combined with metric_constraints.
            RuntimeError: If no start produces a finite objective, or if no
                start produces a feasible-and-successful result when
                metric_constraints are supplied.
        """
        if n_random_starts < 0:
            raise ValueError("n_random_starts must be nonnegative.")

        if objective is None:
            objective = _default_doping_objective

        if metric_constraints is None:
            validated_constraints: tuple[MetricConstraint, ...] = ()
        else:
            validated_constraints = _validate_metric_constraints(metric_constraints)

        if validated_constraints:
            if method == "L-BFGS-B":
                raise ValueError("L-BFGS-B cannot enforce nonlinear metric constraints; use SLSQP.")
            resolved_method = method if method is not None else "SLSQP"
        else:
            resolved_method = method if method is not None else "L-BFGS-B"

        lower_bounds, upper_bounds = self.doping_optimization_bounds(parameter_constraints)
        optimizer_bounds, _, _ = self._optimizer_bounds(
            lower_bounds,
            upper_bounds,
            log_dose,
        )
        dose_bounds = (float(lower_bounds[3]), float(upper_bounds[3]))

        physical_starts: list[np.ndarray] = []
        if initial_states is not None:
            physical_starts.extend(self._validate_doping_state(state).copy() for state in initial_states)

        if n_random_starts > 0:
            if random_seed is None:
                physical_starts.extend(
                    self.doping_initial_conditions_random(parameter_constraints) for _ in range(n_random_starts)
                )
            else:
                random_state = np.random.get_state()
                try:
                    np.random.seed(random_seed)
                    physical_starts.extend(
                        self.doping_initial_conditions_random(parameter_constraints) for _ in range(n_random_starts)
                    )
                finally:
                    np.random.set_state(random_state)

        if not physical_starts:
            physical_starts.append(self.doping_current_state())

        def metrics_from_optimizer_state(optimizer_state: np.ndarray) -> DopingProcessMetrics:
            physical_state = self._to_physical_state(optimizer_state, log_dose, dose_bounds)
            candidate = self._copy_with_doping_state(physical_state)
            return candidate.doping_process_metrics(device_thickness_m)

        def objective_from_optimizer_state(optimizer_state: np.ndarray) -> float:
            return float(objective(metrics_from_optimizer_state(optimizer_state)))

        slsqp_constraints = (
            _build_slsqp_constraints(validated_constraints, metrics_from_optimizer_state)
            if validated_constraints
            else ()
        )

        all_results: list[OptimizeResult] = []
        for physical_start in physical_starts:
            optimizer_start = self._to_optimizer_state(physical_start, log_dose)
            minimize_kwargs: dict = {
                "method": resolved_method,
                "bounds": optimizer_bounds,
            }
            if slsqp_constraints:
                minimize_kwargs["constraints"] = slsqp_constraints
            result = minimize(
                objective_from_optimizer_state,
                optimizer_start,
                **minimize_kwargs,
            )
            all_results.append(result)

        if validated_constraints:
            feasible_results: list[OptimizeResult] = []
            for result in all_results:
                if not result.success or not np.isfinite(result.fun):
                    continue
                candidate_metrics = metrics_from_optimizer_state(result.x)
                if _constraints_satisfied(candidate_metrics, validated_constraints):
                    feasible_results.append(result)

            if not feasible_results:
                raise RuntimeError("No feasible doping solution satisfies the metric constraints.")

            scipy_result = min(feasible_results, key=lambda result: result.fun)
        else:
            finite_results = [result for result in all_results if np.isfinite(result.fun)]
            if not finite_results:
                raise RuntimeError("Doping optimization did not produce a finite objective value.")
            scipy_result = min(finite_results, key=lambda result: result.fun)

        state = self._to_physical_state(scipy_result.x, log_dose, dose_bounds)
        optimized = self._copy_with_doping_state(state)
        metrics = optimized.doping_process_metrics(device_thickness_m)

        return DopingOptimizationResult(
            optimized=optimized,
            metrics=metrics,
            objective_value=float(scipy_result.fun),
            state=state,
            scipy_result=scipy_result,
            all_results=tuple(all_results),
            boundary_flags=_doping_boundary_flags(state, lower_bounds, upper_bounds),
        )
