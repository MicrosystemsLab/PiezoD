"""State-vector machinery for cantilever optimization.

Each cantilever subclass declares the variables it wants to optimize over by
returning a tuple of :class:`StateVar` from ``optimization_state_vars()``.
This module provides the helpers that convert between physical attribute
values on a cantilever, optimizer-space state vectors (log-transformed and
scaled to O(1)), and bounds derived from defaults plus user overrides.

The design mirrors the MATLAB implementation's pattern of a scaling vector
plus per-class bound/initial-condition hooks, but expressed as a single
declarative spec instead of five parallel methods.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class StateVar:
    """Declarative spec for one optimizable state variable.

    Attributes:
        name: Cantilever attribute name (e.g. ``"l"``, ``"dopant_concentration"``).
            The optimizer reads/writes this attribute via getattr/setattr.
        scale: Multiplicative factor applied in optimizer space so the scaled
            state lands near O(1). Mirrors MATLAB's ``optimization_scaling``
            vector. Example: ``"l"`` (meters, ~1e-4) uses ``scale = 1e5``.
        default_min: Lower bound in physical units, used unless overridden by
            ``parameter_constraints`` keyed ``"min_<name>"``.
        default_max: Upper bound in physical units, used unless overridden by
            ``parameter_constraints`` keyed ``"max_<name>"``.
        log_scale: If True, optimize log10(value) instead of value. Used for
            quantities that span many orders of magnitude (e.g. dopant
            concentration, implant dose). Bounds are interpreted in physical
            space; the log transform is applied internally.
    """

    name: str
    scale: float
    default_min: float
    default_max: float
    log_scale: bool = False


def physical_state_from_cantilever(
    cantilever: Any,
    state_vars: Sequence[StateVar],
) -> np.ndarray:
    """Read the current physical state vector from a cantilever instance."""
    return np.array([float(getattr(cantilever, var.name)) for var in state_vars], dtype=float)


def apply_physical_state(
    cantilever: Any,
    state_vars: Sequence[StateVar],
    physical_state: np.ndarray,
) -> None:
    """Apply a physical state vector by setting the corresponding attributes.

    After all state attributes are set, calls
    ``cantilever.optimization_post_apply()`` if the method exists. Subclasses
    use this to keep derived attributes in sync (e.g.
    ``CantileverPiezoelectric`` mirrors ``l_si``/``w_si`` to ``l_pe``/``w_pe``,
    ``CantileverPoly`` enforces top/bottom thickness symmetry).
    """
    if len(physical_state) != len(state_vars):
        raise ValueError(f"physical_state has length {len(physical_state)}, expected {len(state_vars)}.")
    for var, value in zip(state_vars, physical_state, strict=True):
        setattr(cantilever, var.name, float(value))
    post_apply = getattr(cantilever, "optimization_post_apply", None)
    if post_apply is not None:
        post_apply()


def to_optimizer_state(
    physical_state: np.ndarray,
    state_vars: Sequence[StateVar],
) -> np.ndarray:
    """Convert a physical state vector to optimizer (scaled, log-where-applicable) space."""
    physical_state = np.asarray(physical_state, dtype=float)
    if len(physical_state) != len(state_vars):
        raise ValueError(f"physical_state has length {len(physical_state)}, expected {len(state_vars)}.")
    out = np.empty_like(physical_state)
    for i, var in enumerate(state_vars):
        value = physical_state[i]
        if var.log_scale:
            if value <= 0:
                raise ValueError(f"State variable {var.name!r} has log_scale=True but value {value} is not positive.")
            value = np.log10(value)
        out[i] = value * var.scale
    return out


def to_physical_state(
    optimizer_state: np.ndarray,
    state_vars: Sequence[StateVar],
) -> np.ndarray:
    """Invert the optimizer-space transform back to physical units."""
    optimizer_state = np.asarray(optimizer_state, dtype=float)
    if len(optimizer_state) != len(state_vars):
        raise ValueError(f"optimizer_state has length {len(optimizer_state)}, expected {len(state_vars)}.")
    out = np.empty_like(optimizer_state)
    for i, var in enumerate(state_vars):
        value = optimizer_state[i] / var.scale
        if var.log_scale:
            value = 10.0**value
        out[i] = value
    return out


def physical_bounds(
    state_vars: Sequence[StateVar],
    parameter_constraints: dict[str, float] | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Build (lower, upper) bound vectors in physical units.

    ``parameter_constraints`` overrides defaults using keys ``"min_<name>"``
    and ``"max_<name>"``. Unknown keys raise :class:`ValueError` so typos
    fail loudly instead of silently being ignored.
    """
    lower = np.array([var.default_min for var in state_vars], dtype=float)
    upper = np.array([var.default_max for var in state_vars], dtype=float)

    if parameter_constraints is None:
        return lower, upper

    name_to_index = {var.name: i for i, var in enumerate(state_vars)}
    for key, value in parameter_constraints.items():
        if key.startswith("min_"):
            target = key[4:]
            if target not in name_to_index:
                raise ValueError(
                    f"parameter_constraints key {key!r} does not match any state variable. "
                    f"Known names: {sorted(name_to_index)}"
                )
            lower[name_to_index[target]] = float(value)
        elif key.startswith("max_"):
            target = key[4:]
            if target not in name_to_index:
                raise ValueError(
                    f"parameter_constraints key {key!r} does not match any state variable. "
                    f"Known names: {sorted(name_to_index)}"
                )
            upper[name_to_index[target]] = float(value)
        else:
            raise ValueError(f"parameter_constraints key {key!r} must start with 'min_' or 'max_'.")

    if np.any(lower > upper):
        bad = [(var.name, float(lower[i]), float(upper[i])) for i, var in enumerate(state_vars) if lower[i] > upper[i]]
        raise ValueError(f"parameter_constraints produced lower > upper for: {bad}")

    return lower, upper


def optimizer_bounds(
    state_vars: Sequence[StateVar],
    parameter_constraints: dict[str, float] | None,
) -> list[tuple[float, float]]:
    """Build SciPy-style bound list in optimizer (scaled, log) space."""
    lower, upper = physical_bounds(state_vars, parameter_constraints)
    lo_opt = to_optimizer_state(lower, state_vars)
    hi_opt = to_optimizer_state(upper, state_vars)
    return list(zip(lo_opt.tolist(), hi_opt.tolist(), strict=True))


def random_physical_state(
    state_vars: Sequence[StateVar],
    parameter_constraints: dict[str, float] | None,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample a uniformly-random physical state within bounds.

    For ``log_scale`` variables, samples uniformly in log10 space so that
    e.g. dopant concentration spans the full bound range rather than
    clustering near the upper end. Mirrors MATLAB's
    ``doping_initial_conditions_random`` for epitaxy.
    """
    lower, upper = physical_bounds(state_vars, parameter_constraints)
    out = np.empty(len(state_vars), dtype=float)
    for i, var in enumerate(state_vars):
        lo = lower[i]
        hi = upper[i]
        if var.log_scale:
            if lo <= 0 or hi <= 0:
                raise ValueError(
                    f"State variable {var.name!r} has log_scale=True but bounds ({lo}, {hi}) are not strictly positive."
                )
            log_lo = np.log10(lo)
            log_hi = np.log10(hi)
            out[i] = 10.0 ** (log_lo + rng.random() * (log_hi - log_lo))
        else:
            out[i] = lo + rng.random() * (hi - lo)
    return out
