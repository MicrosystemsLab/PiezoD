"""Named objective factories matching MATLAB's ``goal*`` constants.

MATLAB exposed four goals (``goalForceResolution``,
``goalDisplacementResolution``, ``goalForceNoiseDensity``,
``goalSurfaceStress``) and the piezoelectric subclass added voltage-mode
and charge-mode variants. Each MATLAB ``optimize_*`` wrapper applied a
unit-scaling factor purely for fmincon's numerical conditioning: pN, nm,
pN/sqrt(Hz), MPa. The Python factories preserve those scaling choices so
optimizer-space objective values stay in a sensible range and so test
expectations match the MATLAB tutorial's reported numbers.

The MATLAB ``goalForceNoiseDensity`` factory is intentionally not provided
here. Its Python counterpart depends on ``Cantilever.force_noise_density``,
which currently calls ``voltage_noise`` with a scalar; ``voltage_noise``
uses ``math.sqrt`` and is unable to evaluate at a single frequency
without raising. Once that scalar/array handling is fixed in the noise
chain, a force-noise-density goal can be added in a follow-up PR.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

# All goals share the signature ``Callable[[Cantilever], float]`` and are
# typed as ``Callable[[Any], float]`` here because the cantilever
# subclasses (Epitaxy/Diffusion/Implantation, Poly, Piezoelectric) do not
# share a common base class and this avoids an import cycle with
# cantilever.py.
Goal = Callable[[Any], float]


def force_resolution_goal() -> Goal:
    """Goal: minimum detectable force in pN (force_resolution * 1e12)."""
    return lambda c: float(c.force_resolution()) * 1e12


def displacement_resolution_goal() -> Goal:
    """Goal: minimum detectable displacement in nm (displacement_resolution * 1e9)."""
    return lambda c: float(c.displacement_resolution()) * 1e9


def surface_stress_resolution_goal() -> Goal:
    """Goal: minimum detectable surface stress in N/mm (surface_stress_resolution * 1e6).

    Matches MATLAB's ``optimize_surface_stress_resolution`` scaling. The
    underlying SI value is in N/m, so 1e6 converts to N/mm.
    """
    return lambda c: float(c.surface_stress_resolution()) * 1e6


def voltage_force_resolution_goal() -> Goal:
    """Piezoelectric voltage-mode force resolution in pN (Fminv * 1e12)."""
    return lambda c: float(c.Fminv()) * 1e12


def voltage_displacement_resolution_goal() -> Goal:
    """Piezoelectric voltage-mode displacement resolution in nm (Xminv * 1e9)."""
    return lambda c: float(c.Xminv()) * 1e9


def charge_force_resolution_goal() -> Goal:
    """Piezoelectric charge-mode force resolution in pN (Fminq * 1e12)."""
    return lambda c: float(c.Fminq()) * 1e12


def charge_displacement_resolution_goal() -> Goal:
    """Piezoelectric charge-mode displacement resolution in nm (Xminq * 1e9)."""
    return lambda c: float(c.Xminq()) * 1e9
