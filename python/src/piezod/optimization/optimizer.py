"""Multi-start cantilever optimization built on ``scipy.optimize.minimize``.

Two entry points mirror MATLAB:

- :func:`optimize_performance_from_current` -- single shot from the
  cantilever's current state, equivalent to MATLAB
  ``optimize_performance_from_current``.
- :func:`optimize_performance` -- multi-start with random initial
  conditions and a convergence check that stops early when the two best
  results agree within ``convergence_tolerance``. Equivalent to MATLAB
  ``optimize_performance`` driven by ``cantilever.numOptimizationIterations``.

The state machinery in :mod:`piezod.optimization.state` and constraint
builders in :mod:`piezod.optimization.constraints` keep this module
focused on the optimizer dispatch and convergence logic.
"""

from __future__ import annotations

import copy
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.optimize import OptimizeResult, minimize

from piezod.optimization.constraints import (
    CantileverMetricConstraint,
    build_default_aspect_ratio_constraints,
    build_metric_constraint_funcs,
)
from piezod.optimization.state import (
    StateVar,
    apply_physical_state,
    optimizer_bounds,
    physical_state_from_cantilever,
    random_physical_state,
    to_optimizer_state,
    to_physical_state,
)

_DEFAULT_TOL = 1e-9
_DEFAULT_MAXITER = 2000
_INFEASIBLE_OBJECTIVE = 1e30


@dataclass(frozen=True)
class OptimizationResult:
    """Return value of :func:`optimize_performance` and friends.

    Attributes:
        optimized: A copy of the input cantilever with the optimal state
            applied. Original cantilever is left untouched.
        objective_value: Objective function value at ``optimized``, in the
            same units as the goal callable returned (e.g. pN for
            :func:`force_resolution_goal`).
        physical_state: Final state vector in physical units, in the order
            given by ``cantilever.optimization_state_vars()``.
        scipy_result: The ``OptimizeResult`` from the best run.
        all_results: Tuple of ``OptimizeResult`` for every multi-start.
    """

    optimized: object
    objective_value: float
    physical_state: np.ndarray
    scipy_result: OptimizeResult
    all_results: tuple[OptimizeResult, ...]


def _resolve_method(
    method: str | None,
    has_constraints: bool,
) -> str:
    if method is not None:
        if method == "L-BFGS-B" and has_constraints:
            raise ValueError("L-BFGS-B cannot enforce nonlinear constraints; use SLSQP.")
        return method
    return "SLSQP" if has_constraints else "L-BFGS-B"


def _run_one(
    cantilever_template: Any,
    state_vars: Sequence[StateVar],
    objective: Callable[[object], float],
    optimizer_start: np.ndarray,
    bounds: list[tuple[float, float]],
    constraints: list[dict],
    method: str,
    options: dict,
) -> OptimizeResult:
    def objective_from_optimizer_state(x: np.ndarray) -> float:
        physical = to_physical_state(x, state_vars)
        candidate = copy.copy(cantilever_template)
        try:
            apply_physical_state(candidate, state_vars, physical)
            value = float(objective(candidate))
        except (ValueError, FloatingPointError, ZeroDivisionError):
            return _INFEASIBLE_OBJECTIVE
        if not math.isfinite(value):
            return _INFEASIBLE_OBJECTIVE
        return value

    minimize_kwargs: dict = {
        "method": method,
        "bounds": bounds,
        "options": options,
    }
    if constraints:
        minimize_kwargs["constraints"] = constraints

    return minimize(objective_from_optimizer_state, optimizer_start, **minimize_kwargs)


def _build_optimized_cantilever(
    cantilever_template: Any,
    state_vars: Sequence[StateVar],
    optimizer_state: np.ndarray,
) -> tuple[object, np.ndarray]:
    physical = to_physical_state(optimizer_state, state_vars)
    optimized = copy.copy(cantilever_template)
    apply_physical_state(optimized, state_vars, physical)
    return optimized, physical


def optimize_performance_from_current(
    cantilever: Any,
    objective: Callable[[object], float],
    *,
    parameter_constraints: dict[str, float] | None = None,
    metric_constraints: Sequence[CantileverMetricConstraint] | None = None,
    default_aspect_constraints: bool = True,
    method: str | None = None,
    options: dict | None = None,
) -> OptimizationResult:
    """Single-shot optimization starting from the cantilever's current state.

    Mirrors MATLAB ``optimize_performance_from_current``. Useful for
    refining an existing design rather than searching the full bound box.

    Args:
        cantilever: Instance of a :class:`Cantilever` subclass that
            implements ``optimization_state_vars()``. The instance is not
            mutated; the result is a copy.
        objective: Callable mapping a cantilever to a scalar to minimize.
            Use one of the named goal factories from
            :mod:`piezod.optimization.goals` or supply your own.
        parameter_constraints: Optional dict overriding default state
            bounds. Keys are ``"min_<name>"`` / ``"max_<name>"`` for any
            ``StateVar.name`` returned by ``optimization_state_vars()``.
        metric_constraints: Sequence of :class:`CantileverMetricConstraint`
            inequality constraints on derived metrics
            (``power_dissipation``, ``stiffness``, etc.).
        default_aspect_constraints: When True, add the standard
            ``l/w >= 2``, ``w/t >= 2``, ``l_pr/w_pr >= 2``, ``l_pr >= 2 um``
            sanity constraints. Mirrors MATLAB's always-on aspect ratio
            checks. Disable for advanced cases (sub-2 aspect ratios).
        method: SciPy ``minimize`` method. Default chooses SLSQP if any
            nonlinear constraints are present, L-BFGS-B otherwise.
        options: Override SciPy ``minimize`` options. Defaults: ``ftol``,
            ``maxiter``, set to match MATLAB's ``TolFun``/``MaxIter``.

    Returns:
        :class:`OptimizationResult` with the optimized cantilever and full
        SciPy diagnostic information.
    """
    state_vars = cantilever.optimization_state_vars()
    bounds = optimizer_bounds(state_vars, parameter_constraints)
    physical_start = physical_state_from_cantilever(cantilever, state_vars)
    physical_start = _clip_to_bounds(physical_start, state_vars, parameter_constraints)
    optimizer_start = to_optimizer_state(physical_start, state_vars)

    constraints = list(metric_constraints) if metric_constraints else []
    slsqp_constraints = build_metric_constraint_funcs(cantilever, state_vars, constraints)
    if default_aspect_constraints:
        slsqp_constraints.extend(build_default_aspect_ratio_constraints(cantilever, state_vars))

    resolved_method = _resolve_method(method, has_constraints=bool(slsqp_constraints))
    resolved_options = _build_options(options)

    scipy_result = _run_one(
        cantilever_template=cantilever,
        state_vars=state_vars,
        objective=objective,
        optimizer_start=optimizer_start,
        bounds=bounds,
        constraints=slsqp_constraints,
        method=resolved_method,
        options=resolved_options,
    )

    optimized, physical = _build_optimized_cantilever(cantilever, state_vars, scipy_result.x)
    return OptimizationResult(
        optimized=optimized,
        objective_value=float(scipy_result.fun),
        physical_state=physical,
        scipy_result=scipy_result,
        all_results=(scipy_result,),
    )


def optimize_performance(
    cantilever: Any,
    objective: Callable[[object], float],
    *,
    parameter_constraints: dict[str, float] | None = None,
    metric_constraints: Sequence[CantileverMetricConstraint] | None = None,
    default_aspect_constraints: bool = True,
    n_starts: int = 5,
    max_iterations: int = 20,
    convergence_tolerance: float = 0.01,
    include_current_state: bool = False,
    random_seed: int | None = None,
    method: str | None = None,
    options: dict | None = None,
) -> OptimizationResult:
    """Multi-start optimization with early-stop convergence check.

    Mirrors MATLAB ``optimize_performance``: run optimization from random
    initial conditions until the two best objective values agree within
    ``convergence_tolerance``, capped at ``max_iterations`` runs total.

    Args:
        cantilever: Instance to optimize (not mutated).
        objective: Callable returning the scalar to minimize.
        parameter_constraints: Optional state-bound overrides.
        metric_constraints: Optional derived-metric inequality constraints.
        default_aspect_constraints: Add standard geometric sanity checks.
        n_starts: Number of random starts to run before checking convergence.
            After that point, additional starts run one at a time until
            convergence or ``max_iterations`` is reached. Mirrors MATLAB's
            iterative ``ii > 1`` check.
        max_iterations: Hard upper bound on total optimization runs.
            Equivalent to MATLAB's ``cantilever.numOptimizationIterations``.
        convergence_tolerance: Two best objective values must agree within
            this relative tolerance for early stop. ``0.01`` = 1%, the
            MATLAB default.
        include_current_state: When True, the cantilever's current state
            is used as the first start before any random ones.
        random_seed: Optional seed for reproducibility.
        method: SciPy method override.
        options: SciPy options override.

    Returns:
        :class:`OptimizationResult` whose ``optimized`` cantilever is the
        best of all starts. ``all_results`` contains every SciPy run.
    """
    if n_starts < 1:
        raise ValueError(f"n_starts must be >= 1, got {n_starts}.")
    if max_iterations < n_starts:
        raise ValueError(f"max_iterations ({max_iterations}) must be >= n_starts ({n_starts}).")

    state_vars = cantilever.optimization_state_vars()
    bounds = optimizer_bounds(state_vars, parameter_constraints)

    user_constraints = list(metric_constraints) if metric_constraints else []
    slsqp_constraints = build_metric_constraint_funcs(cantilever, state_vars, user_constraints)
    if default_aspect_constraints:
        slsqp_constraints.extend(build_default_aspect_ratio_constraints(cantilever, state_vars))

    resolved_method = _resolve_method(method, has_constraints=bool(slsqp_constraints))
    resolved_options = _build_options(options)

    rng = np.random.default_rng(random_seed)
    all_results: list[OptimizeResult] = []
    objective_values: list[float] = []

    def run_start(physical_start: np.ndarray) -> None:
        clipped = _clip_to_bounds(physical_start, state_vars, parameter_constraints)
        optimizer_start = to_optimizer_state(clipped, state_vars)
        result = _run_one(
            cantilever_template=cantilever,
            state_vars=state_vars,
            objective=objective,
            optimizer_start=optimizer_start,
            bounds=bounds,
            constraints=slsqp_constraints,
            method=resolved_method,
            options=resolved_options,
        )
        all_results.append(result)
        # SciPy SLSQP often reports success=False on constrained problems
        # ("inequality constraints incompatible") even when it found a
        # useful local minimum. Filter only on the objective being finite
        # so we don't discard genuinely good solutions. MATLAB filtered on
        # exitflag in [1, 2]; SciPy's analog is too coarse.
        if math.isfinite(result.fun) and result.fun < _INFEASIBLE_OBJECTIVE / 2:
            objective_values.append(float(result.fun))
        else:
            objective_values.append(math.inf)

    if include_current_state:
        run_start(physical_state_from_cantilever(cantilever, state_vars))

    while len(all_results) < n_starts:
        run_start(random_physical_state(state_vars, parameter_constraints, rng))

    while len(all_results) < max_iterations:
        finite = sorted(v for v in objective_values if math.isfinite(v))
        if len(finite) >= 2:
            best, second = finite[0], finite[1]
            if best > 0 and abs(1.0 - best / second) < convergence_tolerance:
                break
        run_start(random_physical_state(state_vars, parameter_constraints, rng))

    finite_indices = [i for i, v in enumerate(objective_values) if math.isfinite(v)]
    if not finite_indices:
        raise RuntimeError(
            f"All {len(all_results)} optimization starts failed to produce a finite objective. "
            "Check parameter_constraints, metric_constraints, and the objective callable."
        )

    best_index = min(finite_indices, key=lambda i: objective_values[i])
    best_result = all_results[best_index]

    optimized, physical = _build_optimized_cantilever(cantilever, state_vars, best_result.x)
    return OptimizationResult(
        optimized=optimized,
        objective_value=float(best_result.fun),
        physical_state=physical,
        scipy_result=best_result,
        all_results=tuple(all_results),
    )


def _clip_to_bounds(
    physical_state: np.ndarray,
    state_vars: Sequence[StateVar],
    parameter_constraints: dict[str, float] | None,
) -> np.ndarray:
    """Project a physical state vector inside the bound box.

    SciPy's SLSQP requires the initial guess to satisfy bounds (otherwise
    it can produce NaNs in the first step). Pure-clipping is sufficient
    here -- we don't need a feasible-with-respect-to-nonlinear-constraints
    start, just a bound-feasible one.
    """
    from piezod.optimization.state import physical_bounds

    lower, upper = physical_bounds(state_vars, parameter_constraints)
    return np.clip(physical_state, lower, upper)


def _build_options(user_options: dict | None) -> dict:
    options = {"ftol": _DEFAULT_TOL, "maxiter": _DEFAULT_MAXITER}
    if user_options:
        options.update(user_options)
    return options
