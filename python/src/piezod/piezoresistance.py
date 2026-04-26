"""Crystal-orientation-dependent piezoresistance coefficients and dR/R helpers.

Single source of truth for the low-doping longitudinal and transverse
piezoresistance coefficients (pi_l, pi_t) of single-crystal silicon at 300 K
along the standard wafer / current orientations used in MEMS work.
Values follow Smith (1954) and Kanda (1982).

The doping-concentration reduction P(n) is applied separately (see
`Cantilever.piezoresistance_factor`) -- the values exposed here are the
low-doping limits, signed.

Convention: positive sigma_l / sigma_t means tension along / perpendicular
to the current direction.

    dR/R = beta * (pi_l * sigma_l + pi_t * sigma_t)

with beta the bending-weighted process efficiency factor.
"""

from __future__ import annotations

from enum import Enum
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike


class CrystalOrientation(str, Enum):
    """Wafer surface orientation and in-plane current direction.

    `WAFER_100_DIR_100`: (100) wafer, current along a <100> direction.
    `WAFER_100_DIR_110`: (100) wafer, current along a <110> direction.
    """

    WAFER_100_DIR_100 = "100_along_100"
    WAFER_100_DIR_110 = "100_along_110"


# Smith (1954) / Kanda (1982) low-doping coefficients in 1e-11 Pa^-1 (i.e. 1/TPa).
# Signs included. Keys: (doping_type, orientation).
_PI_LOW_DOPING_TPA: dict[Tuple[str, CrystalOrientation], Tuple[float, float]] = {
    ("phosphorus", CrystalOrientation.WAFER_100_DIR_100): (-102.2, 53.4),
    ("phosphorus", CrystalOrientation.WAFER_100_DIR_110): (-31.2, -17.6),
    ("arsenic", CrystalOrientation.WAFER_100_DIR_100): (-102.2, 53.4),
    ("arsenic", CrystalOrientation.WAFER_100_DIR_110): (-31.2, -17.6),
    ("boron", CrystalOrientation.WAFER_100_DIR_100): (6.6, -1.1),
    ("boron", CrystalOrientation.WAFER_100_DIR_110): (71.8, -66.3),
}


def default_orientation(doping_type: str) -> CrystalOrientation:
    """Default crystal orientation for a given dopant.

    n-type (phosphorus, arsenic) is most sensitive along <100>;
    p-type (boron) is most sensitive along <110>.
    """
    if doping_type in ("phosphorus", "arsenic"):
        return CrystalOrientation.WAFER_100_DIR_100
    if doping_type == "boron":
        return CrystalOrientation.WAFER_100_DIR_110
    raise ValueError(f"Unknown doping_type: {doping_type!r}. Must be 'boron', 'phosphorus', or 'arsenic'.")


def pi_low_doping(doping_type: str, orientation: CrystalOrientation) -> Tuple[float, float]:
    """Return signed (pi_l, pi_t) low-doping coefficients in Pa^-1.

    Args:
        doping_type: One of 'boron', 'phosphorus', 'arsenic'.
        orientation: Wafer orientation enum.

    Returns:
        (pi_longitudinal, pi_transverse) in Pa^-1, signed.

    Raises:
        ValueError: If the (doping_type, orientation) combination is unknown.
    """
    if not isinstance(orientation, CrystalOrientation):
        raise TypeError(f"orientation must be a CrystalOrientation enum, got {type(orientation).__name__}.")
    key = (doping_type, orientation)
    if key not in _PI_LOW_DOPING_TPA:
        raise ValueError(
            f"No piezoresistance coefficients for doping_type={doping_type!r}, orientation={orientation.value!r}."
        )
    pi_l_tpa, pi_t_tpa = _PI_LOW_DOPING_TPA[key]
    return pi_l_tpa * 1e-11, pi_t_tpa * 1e-11


def rotate_in_plane_stress(
    sxx: ArrayLike,
    syy: ArrayLike,
    sxy: ArrayLike,
    theta_rad: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Rotate an in-plane stress tensor into the resistor's longitudinal frame.

    Standard 2D tensor rotation by `theta_rad` from the wafer x-axis to the
    current direction.

        sigma_l =  c^2 sxx + s^2 syy + 2 c s sxy
        sigma_t =  s^2 sxx + c^2 syy - 2 c s sxy

    Args:
        sxx, syy, sxy: In-plane stress components in Pa. Scalars or arrays.
        theta_rad: Angle from the wafer x-axis to the current direction.

    Returns:
        (sigma_l, sigma_t) arrays of the same broadcast shape as the inputs.
    """
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)
    sxx_a = np.asarray(sxx, dtype=float)
    syy_a = np.asarray(syy, dtype=float)
    sxy_a = np.asarray(sxy, dtype=float)
    sigma_l = c * c * sxx_a + s * s * syy_a + 2 * c * s * sxy_a
    sigma_t = s * s * sxx_a + c * c * syy_a - 2 * c * s * sxy_a
    return sigma_l, sigma_t
