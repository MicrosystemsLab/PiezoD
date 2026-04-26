"""Piezoresistor metrics computed from a user-supplied dopant profile.

Use when you have a measured (e.g. SIMS) profile and want piezod to compute
the bending efficiency factor `beta`, the uniform efficiency factor `beta1`,
sheet resistance, etc., without going through the implantation lookup table.

Profile-only metrics: the doping profile sets `beta`, `beta1`, `beta2_um`,
`sheet_resistance`, `nz`, `peak_concentration_cm3`, `retained_dose_cm2`, and
`junction_depth_m`. The Hooge alpha (`alpha_h`) requires anneal conditions
because it is an empirical fit against sqrt(D*t); pass `annealing_temp` and
`annealing_time`. Mobility mu(n) and the doping-dependent piezoresistance
reduction P(n) come from `Cantilever`'s Reggiani / Richter lookups.

The substrate is always counter-doped to the piezoresistor with
`substrate_background_cm3`; net active carriers in the resistor are
`max(0, active_cm3 - substrate_background_cm3)`. All carrier integrals
(sheet resistance, Nz, beta, beta1) use the net active concentration so they
naturally truncate at the junction.
"""

from __future__ import annotations

from typing import Literal, Tuple

import numpy as np
from numpy.typing import ArrayLike

from piezod.cantilever import Cantilever
from piezod.cantilever_implantation import DopingProcessMetrics, hooge_alpha_from_anneal
from piezod.piezoresistance import (
    CrystalOrientation,
    default_orientation,
    pi_low_doping,
    rotate_in_plane_stress,
)


class PiezoresistorFromProfile:
    """Piezoresistor model driven by a user-supplied dopant profile.

    Args:
        depth_m: 1D monotonic depth array (m). z=0 is the implant surface;
            increasing depth points into the substrate. Must have at least
            two points.
        active_cm3: Electrically active dopant species concentration at each
            depth (cm^-3), same length as depth_m. Should *not* include the
            substrate background -- supply ``substrate_background_cm3``
            separately and the class computes net active carriers internally.
        device_thickness_m: Total thickness of the device layer (m). Used
            for the bending weight in `beta` and for clipping the profile.
        doping_type: 'boron', 'phosphorus', or 'arsenic'.
        annealing_temp: Anneal temperature (K). Required for `alpha_h` /
            `doping_process_metrics`; otherwise unused.
        annealing_time: Anneal time (s). Required for `alpha_h` /
            `doping_process_metrics`; otherwise unused.
        crystal_orientation: Wafer surface and current direction. Defaults
            to <100> for n-type and <110> for boron.
        temperature: Operating temperature for mobility / P(n) lookups (K).
            Defaults to 300 K.
        substrate_background_cm3: Counter-doped substrate concentration
            (cm^-3). Defaults to 1e15. Subtracted from `active_cm3` to obtain
            net active resistor carriers; also sets where the profile crosses
            into the substrate (junction).
    """

    def __init__(
        self,
        depth_m: ArrayLike,
        active_cm3: ArrayLike,
        *,
        device_thickness_m: float,
        doping_type: Literal["boron", "phosphorus", "arsenic"],
        annealing_temp: float,
        annealing_time: float,
        crystal_orientation: CrystalOrientation | None = None,
        temperature: float = 300.0,
        substrate_background_cm3: float = 1e15,
    ):
        depth = np.asarray(depth_m, dtype=float)
        active = np.asarray(active_cm3, dtype=float)
        if depth.ndim != 1 or active.ndim != 1:
            raise ValueError("depth_m and active_cm3 must be 1D arrays.")
        if depth.shape != active.shape:
            raise ValueError(f"depth_m and active_cm3 must have the same length, got {depth.shape} and {active.shape}.")
        if depth.size < 2:
            raise ValueError("depth_m / active_cm3 must contain at least two samples.")
        if not np.all(np.diff(depth) > 0):
            raise ValueError("depth_m must be strictly increasing.")
        if device_thickness_m <= 0:
            raise ValueError(f"device_thickness_m must be positive, got {device_thickness_m!r}.")
        if doping_type not in ("boron", "phosphorus", "arsenic"):
            raise ValueError(f"Unknown doping_type: {doping_type!r}. Must be 'boron', 'phosphorus', or 'arsenic'.")
        if annealing_temp <= 0:
            raise ValueError(f"annealing_temp must be positive (Kelvin), got {annealing_temp!r}.")
        if annealing_time <= 0:
            raise ValueError(f"annealing_time must be positive (seconds), got {annealing_time!r}.")
        if substrate_background_cm3 < 0:
            raise ValueError(f"substrate_background_cm3 must be nonnegative, got {substrate_background_cm3!r}.")

        # Clip profile to the device thickness so integrals don't include
        # samples that are physically outside the device layer.
        mask = depth <= device_thickness_m
        if mask.sum() < 2:
            raise ValueError(
                "device_thickness_m clips the profile to fewer than two samples; "
                "supply a thicker device or a finer profile."
            )

        self.depth_m: np.ndarray = depth[mask]
        self.active_cm3: np.ndarray = active[mask]
        self.device_thickness_m: float = float(device_thickness_m)
        self.doping_type: str = doping_type
        self.annealing_temp: float = float(annealing_temp)
        self.annealing_time: float = float(annealing_time)
        self.temperature: float = float(temperature)
        self.substrate_background_cm3: float = float(substrate_background_cm3)
        self.crystal_orientation: CrystalOrientation = (
            crystal_orientation if crystal_orientation is not None else default_orientation(doping_type)
        )

        # Reuse Cantilever's mobility / P(n) implementations as the single source
        # of truth. We bypass cantilever-specific defaults by overriding only
        # what those two methods read.
        self._physics: Cantilever = Cantilever()
        self._physics.doping_type = self.doping_type
        self._physics.T = self.temperature
        # piezoresistor_temp() returns self.T when thermal_modeling == "none"
        # (the Cantilever default), so P(n) uses self.temperature directly.

    # ----- piezoresistance coefficients & dR/R -----

    def pi_longitudinal(self) -> float:
        """Signed low-doping longitudinal piezoresistance coefficient (1/Pa)."""
        pi_l, _ = pi_low_doping(self.doping_type, self.crystal_orientation)
        return pi_l

    def pi_transverse(self) -> float:
        """Signed low-doping transverse piezoresistance coefficient (1/Pa)."""
        _, pi_t = pi_low_doping(self.doping_type, self.crystal_orientation)
        return pi_t

    def dr_over_r(self, sigma_l: ArrayLike, sigma_t: ArrayLike) -> np.ndarray | float:
        """Bending dR/R for an in-plane biaxial stress at the resistor surface.

        See `Cantilever.dr_over_r` for the formula and the shallow-implant
        approximation note.
        """
        beta = self.beta()
        pi_l = self.pi_longitudinal()
        pi_t = self.pi_transverse()
        sigma_l_a = np.asarray(sigma_l, dtype=float)
        sigma_t_a = np.asarray(sigma_t, dtype=float)
        result = beta * (pi_l * sigma_l_a + pi_t * sigma_t_a)
        if result.ndim == 0:
            return float(result)
        return result

    def dr_over_r_from_tensor(
        self,
        sxx: ArrayLike,
        syy: ArrayLike,
        sxy: ArrayLike,
        theta_rad: float,
    ) -> np.ndarray | float:
        """Bending dR/R for an in-plane stress tensor; rotates then calls dr_over_r."""
        sigma_l, sigma_t = rotate_in_plane_stress(sxx, syy, sxy, theta_rad)
        return self.dr_over_r(sigma_l, sigma_t)

    # ----- profile-derived metrics -----

    def net_active_cm3(self) -> np.ndarray:
        """Net active resistor carriers ``max(0, active_cm3 - substrate_background_cm3)`` (cm^-3)."""
        return np.maximum(0.0, self.active_cm3 - self.substrate_background_cm3)

    def total_cm3(self) -> np.ndarray:
        """Total concentration ``max(active_cm3, substrate_background_cm3)`` (cm^-3).

        Plottable magnitude that floors at the substrate background, mirroring
        the ``total_doping`` returned by ``Cantilever.doping_profile()``.
        """
        return np.maximum(self.active_cm3, self.substrate_background_cm3)

    def _mobility_and_conductivity(self) -> Tuple[np.ndarray, np.ndarray]:
        """Vectorized mobility and conductivity from net active carriers."""
        mu, sigma = self._physics.mobility(self.net_active_cm3(), self.temperature)
        return np.asarray(mu, dtype=float), np.asarray(sigma, dtype=float)

    def _piezoresistance_factor(self) -> np.ndarray:
        """Vectorized P(n) from net active carriers."""
        return np.asarray(self._physics.piezoresistance_factor(self.net_active_cm3()), dtype=float)

    def sheet_resistance(self) -> float:
        """Sheet resistance Rs in ohm/sq, from integral of conductivity over depth."""
        _, sigma = self._mobility_and_conductivity()
        depth_cm = self.depth_m * 1e2
        integral = np.trapezoid(sigma, depth_cm)
        if integral <= 0:
            return float("inf")
        return float(1.0 / integral)

    def Nz(self) -> float:
        """Mobility-weighted effective carrier sheet density (1/m^2).

        Harmon current-crowding formula on net active carriers:
            Nz = (integral mu*N dz)^2 / integral mu^2*N dz
        """
        mu, _ = self._mobility_and_conductivity()
        depth_cm = self.depth_m * 1e2
        n_net = self.net_active_cm3()
        numerator = np.trapezoid(mu * n_net, depth_cm) ** 2
        denominator = np.trapezoid(mu**2 * n_net, depth_cm)
        if denominator <= 0:
            return 0.0
        # cm^-2 to m^-2
        return float(1e4 * numerator / denominator)

    @property
    def junction_depth(self) -> float:
        """Junction depth (m) where the dopant species crosses ``substrate_background_cm3``.

        Linearly interpolated between the bracketing samples. Returns the
        deepest sample depth if the profile never crosses the background.
        """
        diff = self.active_cm3 - self.substrate_background_cm3
        below = diff < 0
        if not below.any():
            return float(self.depth_m[-1])
        first = int(np.argmax(below))
        if first == 0:
            return float(self.depth_m[0])
        z0, z1 = self.depth_m[first - 1], self.depth_m[first]
        d0, d1 = diff[first - 1], diff[first]
        return float(z0 + (z1 - z0) * d0 / (d0 - d1))

    def beta1(self) -> float:
        """Uniform-stress efficiency factor.

            beta1 = integral(sigma * P) / integral(sigma)

        where sigma = q * mu * n is the conductivity profile and P(n) is the
        Richter doping reduction. Equals the response to membrane (mid-plane)
        stress that is constant through the thickness.
        """
        _, sigma = self._mobility_and_conductivity()
        P = self._piezoresistance_factor()
        depth_cm = self.depth_m * 1e2
        denominator = np.trapezoid(sigma, depth_cm)
        if denominator <= 0:
            return 0.0
        numerator = np.trapezoid(sigma * P, depth_cm)
        return float(numerator / denominator)

    def beta(self) -> float:
        """Bending-weighted efficiency factor for stress at z=0 (resistor surface).

            beta = integral(sigma * P * (1 - 2 z / t_dev)) / integral(sigma)

        z=0 is the resistor surface; t_dev is the device thickness. Park
        (JMEMS 2010 Part I, Eq. 18).
        """
        _, sigma = self._mobility_and_conductivity()
        P = self._piezoresistance_factor()
        depth_cm = self.depth_m * 1e2
        t_dev_cm = self.device_thickness_m * 1e2
        weight = 1.0 - 2.0 * depth_cm / t_dev_cm
        denominator = np.trapezoid(sigma, depth_cm)
        if denominator <= 0:
            return 0.0
        numerator = np.trapezoid(sigma * P * weight, depth_cm)
        return float(numerator / denominator)

    def alpha(self) -> float:
        """Hooge 1/f noise parameter from anneal conditions."""
        return hooge_alpha_from_anneal(self.doping_type, self.annealing_temp, self.annealing_time)

    def doping_process_metrics(self) -> DopingProcessMetrics:
        """Process-only metrics from the supplied profile and anneal conditions.

        Returns the same dataclass shape as
        `CantileverImplantation.doping_process_metrics()`.

        `beta2_um` is recovered from beta1 and beta via
            beta2_um = (beta1 - beta) * t_dev_um / 2
        so that `beta = beta1 - 2 * beta2_um / t_dev_um` continues to hold
        (matches the implantation-table convention).
        """
        beta1 = self.beta1()
        beta = self.beta()
        device_thickness_um = self.device_thickness_m * 1e6
        beta2_um = (beta1 - beta) * device_thickness_um / 2.0
        depth_cm = self.depth_m * 1e2

        return DopingProcessMetrics(
            alpha_h=float(self.alpha()),
            nz=float(self.Nz()),
            beta1=float(beta1),
            beta2_um=float(beta2_um),
            beta=float(beta),
            sheet_resistance=float(self.sheet_resistance()),
            junction_depth_m=float(self.junction_depth),
            peak_concentration_cm3=float(np.max(self.active_cm3)),
            retained_dose_cm2=float(np.trapezoid(self.active_cm3, depth_cm)),
        )

    # ----- plotting -----

    def plot_doping_profile(self):
        """Plot the total doping magnitude and the net active carriers.

        Mirrors :meth:`Cantilever.plot_doping_profile`. ``total_cm3()`` floors at
        the substrate background so the substrate region is visible directly.
        """
        import matplotlib.pyplot as plt

        z_nm = self.depth_m * 1e9
        n_net = self.net_active_cm3()

        plt.figure()
        plt.semilogy(z_nm, self.total_cm3(), label="Total dopant")
        plt.semilogy(z_nm, np.maximum(n_net, 1.0), label="Net active carriers")
        xj = self.junction_depth
        if 0.0 < xj < float(self.depth_m[-1]):
            plt.axvline(xj * 1e9, color="green", linestyle="--", label="Junction")
        plt.xlabel("Depth (nm)")
        plt.ylabel(r"Concentration (cm$^{-3}$)")
        plt.legend()
