"""Nonlinear constraint scaffolding for cantilever optimization.

User-supplied :class:`CantileverMetricConstraint` entries plus a small set of
always-on geometric sanity constraints are translated into the
``scipy.optimize.minimize`` SLSQP constraint format. The MATLAB
``optimization_constraints`` method dispatched on string keys
(``'omega_min_hz'``, ``'min_k'``, etc.); the Python equivalent uses an enum
so typos fail at construction time.
"""

from __future__ import annotations

import copy
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from enum import Enum
from typing import Any

import numpy as np

from piezod.optimization.state import (
    StateVar,
    apply_physical_state,
    to_physical_state,
)

_RESIDUAL_SCALE_FLOOR = 1.0
_DEFAULT_MIN_PR_LENGTH_M = 2e-6
_DEFAULT_MIN_L_W_RATIO = 2.0
_DEFAULT_MIN_W_T_RATIO = 2.0
_DEFAULT_MIN_PR_L_W_RATIO = 2.0


class CantileverMetric(str, Enum):
    """Derived metrics that can be constrained during optimization.

    Each value resolves to a method (or a method plus tuple index) on a
    :class:`Cantilever`. Constraints reference metrics by enum so
    misspellings fail at construction time rather than silently being
    skipped, the way MATLAB's ``exist('omega_min_hz', 'var')`` pattern fails.
    """

    FORCE_RESOLUTION = "force_resolution"
    DISPLACEMENT_RESOLUTION = "displacement_resolution"
    POWER_DISSIPATION = "power_dissipation"
    OMEGA_VACUUM_HZ = "omega_vacuum_hz"
    OMEGA_DAMPED_HZ = "omega_damped_hz"
    STIFFNESS = "stiffness"
    # Lumped-circuit (approx) and finite-difference (exact) temperature
    # rises. APPROX is fast and good for optimization; EXACT solves the
    # 1-D FD heat equation along the cantilever and is more accurate when
    # the actuator contributes significant heating.
    TEMP_TIP_APPROX = "temp_tip_approx"
    TEMP_MAX_APPROX = "temp_max_approx"
    TEMP_TIP_EXACT = "temp_tip_exact"
    TEMP_MAX_EXACT = "temp_max_exact"
    TIP_DEFLECTION = "tip_deflection"
    SURFACE_STRESS_RESOLUTION = "surface_stress_resolution"


@dataclass(frozen=True)
class CantileverMetricConstraint:
    """Inequality constraint on a derived cantilever metric.

    Either ``minimum`` or ``maximum`` (or both) must be set. With both, the
    constraint becomes a closed interval. To pin a metric to a specific
    value, set ``minimum == maximum``.
    """

    metric: CantileverMetric
    minimum: float | None = None
    maximum: float | None = None

    def __post_init__(self) -> None:
        if self.minimum is None and self.maximum is None:
            raise ValueError(
                f"CantileverMetricConstraint for {self.metric.value!r} must specify at least one of minimum or maximum."
            )
        if self.minimum is not None and not np.isfinite(self.minimum):
            raise ValueError(f"minimum must be finite, got {self.minimum!r}.")
        if self.maximum is not None and not np.isfinite(self.maximum):
            raise ValueError(f"maximum must be finite, got {self.maximum!r}.")
        if self.minimum is not None and self.maximum is not None and self.minimum > self.maximum:
            raise ValueError(f"minimum ({self.minimum}) must be <= maximum ({self.maximum}).")


def _call_or_attr(obj: Any, name: str) -> float:
    """Read a quantity that might be a method, property, or plain attribute."""
    value = getattr(obj, name)
    if callable(value):
        value = value()
    return float(value)


def evaluate_metric(cantilever: Any, metric: CantileverMetric) -> float:
    """Evaluate a metric on a cantilever instance, handling tuple-returning methods."""
    if metric == CantileverMetric.OMEGA_DAMPED_HZ:
        omega_damped_hz, _ = cantilever.omega_damped_hz_and_Q()
        return float(omega_damped_hz)
    if metric == CantileverMetric.TEMP_MAX_APPROX:
        t_max, _ = cantilever.approxTempRise()
        return float(t_max)
    if metric == CantileverMetric.TEMP_TIP_APPROX:
        _, t_tip = cantilever.approxTempRise()
        return float(t_tip)
    if metric == CantileverMetric.TEMP_MAX_EXACT:
        t_max, _ = cantilever.calculateMaxAndTipTemp()
        return float(t_max)
    if metric == CantileverMetric.TEMP_TIP_EXACT:
        _, t_tip = cantilever.calculateMaxAndTipTemp()
        return float(t_tip)
    if metric == CantileverMetric.TIP_DEFLECTION:
        return float(cantilever.tipDeflection())
    return _call_or_attr(cantilever, metric.value)


def _residual_scale(limit: float) -> float:
    return max(abs(limit), _RESIDUAL_SCALE_FLOOR)


def _make_state_applier(
    cantilever_template: Any,
    state_vars: Sequence[StateVar],
) -> Callable[[np.ndarray], Any]:
    """Build a closure that turns an optimizer-space state into a fresh cantilever copy.

    Uses ``copy.copy`` (shallow) so attached lookup tables / arrays are
    shared across iterations rather than re-allocated. Mirrors the pattern
    used by ``CantileverImplantation._copy_with_doping_state``.
    """

    def apply(optimizer_state: np.ndarray) -> Any:
        physical = to_physical_state(optimizer_state, state_vars)
        candidate = copy.copy(cantilever_template)
        apply_physical_state(candidate, state_vars, physical)
        return candidate

    return apply


def build_metric_constraint_funcs(
    cantilever_template: Any,
    state_vars: Sequence[StateVar],
    metric_constraints: Sequence[CantileverMetricConstraint],
) -> list[dict]:
    """Translate :class:`CantileverMetricConstraint` entries into SLSQP ``ineq`` dicts.

    Each constraint becomes one or two SLSQP entries (one for min, one for max).
    SLSQP treats inequality constraints as ``g(x) >= 0`` for feasible, so the
    residual returned is positive when the bound is satisfied.
    """
    apply_state = _make_state_applier(cantilever_template, state_vars)
    constraints: list[dict] = []

    for constraint in metric_constraints:
        metric = constraint.metric

        if constraint.minimum is not None:
            lo = constraint.minimum
            scale = _residual_scale(lo)

            def make_min_fun(metric=metric, lo=lo, scale=scale):
                def fun(x: np.ndarray) -> float:
                    candidate = apply_state(x)
                    try:
                        value = evaluate_metric(candidate, metric)
                    except (ValueError, FloatingPointError):
                        return -1.0
                    if not np.isfinite(value):
                        return -1.0
                    return (value - lo) / scale

                return fun

            constraints.append({"type": "ineq", "fun": make_min_fun()})

        if constraint.maximum is not None:
            hi = constraint.maximum
            scale = _residual_scale(hi)

            def make_max_fun(metric=metric, hi=hi, scale=scale):
                def fun(x: np.ndarray) -> float:
                    candidate = apply_state(x)
                    try:
                        value = evaluate_metric(candidate, metric)
                    except (ValueError, FloatingPointError):
                        return -1.0
                    if not np.isfinite(value):
                        return -1.0
                    return (hi - value) / scale

                return fun

            constraints.append({"type": "ineq", "fun": make_max_fun()})

    return constraints


def _has_attr(obj: Any, name: str) -> bool:
    try:
        getattr(obj, name)
    except (AttributeError, NotImplementedError):
        return False
    return True


def build_default_aspect_ratio_constraints(
    cantilever_template: Any,
    state_vars: Sequence[StateVar],
) -> list[dict]:
    """Add the always-on geometric sanity constraints from MATLAB.

    MATLAB enforces ``l/w >= 2``, ``w/t >= 2``, ``l_pr/w_pr >= 2``, and
    ``l_pr >= 2 um`` on every standard optimization. These prevent the
    optimizer from picking geometrically silly cantilevers. We add them
    only when the cantilever exposes the relevant attributes -- e.g.
    Piezoelectric uses ``l_si``/``w_si``/``t_si`` and has no piezoresistor,
    so its applicable defaults are different.
    """
    apply_state = _make_state_applier(cantilever_template, state_vars)
    constraints: list[dict] = []

    has_l = _has_attr(cantilever_template, "l")
    has_w = _has_attr(cantilever_template, "w")
    has_t = _has_attr(cantilever_template, "t")
    has_l_si = _has_attr(cantilever_template, "l_si")
    has_w_si = _has_attr(cantilever_template, "w_si")
    has_t_si = _has_attr(cantilever_template, "t_si")
    has_l_pr = _has_attr(cantilever_template, "l_pr")
    has_w_pr = _has_attr(cantilever_template, "w_pr")

    def safe_eval(x: np.ndarray, fn: Callable[[Any], float]) -> float:
        try:
            value = fn(apply_state(x))
        except (ValueError, FloatingPointError, ZeroDivisionError):
            return -1.0
        if not np.isfinite(value):
            return -1.0
        return float(value)

    if has_l and has_w:
        constraints.append(
            {
                "type": "ineq",
                "fun": lambda x: safe_eval(x, lambda c: c.l / c.w - _DEFAULT_MIN_L_W_RATIO),
            }
        )
    elif has_l_si and has_w_si:
        constraints.append(
            {
                "type": "ineq",
                "fun": lambda x: safe_eval(x, lambda c: c.l_si / c.w_si - _DEFAULT_MIN_L_W_RATIO),
            }
        )

    if has_w and has_t:
        constraints.append(
            {
                "type": "ineq",
                "fun": lambda x: safe_eval(x, lambda c: c.w / c.t - _DEFAULT_MIN_W_T_RATIO),
            }
        )
    elif has_w_si and has_t_si:
        constraints.append(
            {
                "type": "ineq",
                "fun": lambda x: safe_eval(x, lambda c: c.w_si / c.t_si - _DEFAULT_MIN_W_T_RATIO),
            }
        )

    if has_l_pr and has_w_pr:
        constraints.append(
            {
                "type": "ineq",
                "fun": lambda x: safe_eval(
                    x,
                    lambda c: _call_or_attr(c, "l_pr") / _call_or_attr(c, "w_pr") - _DEFAULT_MIN_PR_L_W_RATIO,
                ),
            }
        )
        constraints.append(
            {
                "type": "ineq",
                "fun": lambda x: safe_eval(
                    x,
                    lambda c: (_call_or_attr(c, "l_pr") - _DEFAULT_MIN_PR_LENGTH_M) * 1e6,
                ),
            }
        )

    return constraints
