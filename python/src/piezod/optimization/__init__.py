"""Cantilever optimization framework.

Public API mirrors the MATLAB optimization toolbox:

- :func:`optimize_performance` -- multi-start with convergence check
- :func:`optimize_performance_from_current` -- single-shot from current state
- :class:`OptimizationResult` -- structured return value
- :class:`StateVar` -- declarative state-variable spec
- :class:`CantileverMetric` -- enum of constrainable metrics
- :class:`CantileverMetricConstraint` -- inequality constraint dataclass
- Goal factories: :func:`force_resolution_goal`, etc.

Example:
    >>> from piezod import CantileverEpitaxy
    >>> from piezod.optimization import (
    ...     optimize_performance,
    ...     force_resolution_goal,
    ...     CantileverMetric,
    ...     CantileverMetricConstraint,
    ... )
    >>> c = CantileverEpitaxy()  # default values
    >>> result = optimize_performance(
    ...     c,
    ...     force_resolution_goal(),
    ...     metric_constraints=[
    ...         CantileverMetricConstraint(
    ...             CantileverMetric.POWER_DISSIPATION, maximum=2e-3
    ...         ),
    ...     ],
    ...     random_seed=0,
    ... )
    >>> result.optimized.force_resolution()  # doctest: +SKIP
"""

from piezod.optimization.constraints import (
    CantileverMetric,
    CantileverMetricConstraint,
)
from piezod.optimization.goals import (
    charge_displacement_resolution_goal,
    charge_force_resolution_goal,
    displacement_resolution_goal,
    force_resolution_goal,
    surface_stress_resolution_goal,
    voltage_displacement_resolution_goal,
    voltage_force_resolution_goal,
)
from piezod.optimization.optimizer import (
    OptimizationResult,
    optimize_performance,
    optimize_performance_from_current,
)
from piezod.optimization.state import StateVar

__all__ = [
    "CantileverMetric",
    "CantileverMetricConstraint",
    "OptimizationResult",
    "StateVar",
    "charge_displacement_resolution_goal",
    "charge_force_resolution_goal",
    "displacement_resolution_goal",
    "force_resolution_goal",
    "optimize_performance",
    "optimize_performance_from_current",
    "surface_stress_resolution_goal",
    "voltage_displacement_resolution_goal",
    "voltage_force_resolution_goal",
]
