"""Tests for the generic geometry+process cantilever optimizer."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from piezod import (
    CantileverDiffusion,
    CantileverEpitaxy,
    CantileverImplantation,
    CantileverMetric,
    CantileverMetricConstraint,
    CantileverPiezoelectric,
    CantileverPoly,
    Material,
    OptimizationResult,
    PiezoMaterial,
    StateVar,
    charge_force_resolution_goal,
    displacement_resolution_goal,
    force_resolution_goal,
    optimize_performance,
    optimize_performance_from_current,
    surface_stress_resolution_goal,
    voltage_force_resolution_goal,
)
from piezod.optimization.state import (
    physical_bounds,
    physical_state_from_cantilever,
    random_physical_state,
    to_optimizer_state,
    to_physical_state,
)


@pytest.fixture(autouse=True)
def _silence_warnings():
    """SciPy and the cantilever model emit RuntimeWarnings during normal optimizer
    iteration (singular interpolation, divide-by-zero). They are harmless but
    noisy in test output."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def _epitaxy_default() -> CantileverEpitaxy:
    """Return an epitaxial cantilever inside the default optimization bound box."""
    c = CantileverEpitaxy(
        freq_min=10,
        freq_max=1000,
        l=200e-6,
        w=20e-6,
        t=2e-6,
        l_pr_ratio=0.3,
        v_bridge=2.0,
        doping_type="boron",
        dopant_concentration=1e19,
        t_pr_ratio=0.3,
    )
    c.fluid = "vacuum"
    c.number_of_piezoresistors = 4
    return c


def _diffusion_default() -> CantileverDiffusion:
    c = CantileverDiffusion(
        freq_min=10,
        freq_max=1000,
        l=200e-6,
        w=20e-6,
        t=2e-6,
        l_pr_ratio=0.3,
        v_bridge=2.0,
        doping_type="phosphorus",
        diffusion_time=30 * 60,
        diffusion_temp=273 + 850,
    )
    c.fluid = "vacuum"
    c.number_of_piezoresistors = 4
    return c


def _implantation_default() -> CantileverImplantation:
    c = CantileverImplantation(
        freq_min=10,
        freq_max=1000,
        l=200e-6,
        w=20e-6,
        t=2e-6,
        l_pr_ratio=0.3,
        v_bridge=2.0,
        doping_type="boron",
        annealing_time=30 * 60,
        annealing_temp=273.15 + 1000,
        annealing_type="inert",
        implantation_energy=40,
        implantation_dose=1e15,
    )
    c.fluid = "vacuum"
    c.number_of_piezoresistors = 4
    return c


def _poly_default() -> CantileverPoly:
    return CantileverPoly(
        freq_min=10,
        freq_max=1000,
        l=200e-6,
        w=20e-6,
        t_top=200e-9,
        t_mid=1e-6,
        t_bot=200e-9,
        matl_top=Material.POLY,
        matl_mid=Material.SI,
        matl_bot=Material.POLY,
        l_pr_ratio=0.3,
        v_bridge=2.0,
        dopant_concentration=1e19,
        number_of_piezoresistors=2,
        number_of_piezoresistors_on_cantilever=2,
    )


def _piezoelectric_default() -> CantileverPiezoelectric:
    return CantileverPiezoelectric(
        freq_min=10,
        freq_max=1000,
        l_si=200e-6,
        w_si=20e-6,
        t_si=2e-6,
        t_pe=500e-9,
        material=PiezoMaterial.ALN,
        r_shunt=1e12,
    )


# -----------------------------------------------------------------------------
# State machinery
# -----------------------------------------------------------------------------


class TestStateMachinery:
    def test_state_var_round_trip_linear(self):
        var = StateVar("l", 1e5, 10e-6, 3e-3)
        physical = np.array([1e-4])
        opt = to_optimizer_state(physical, [var])
        back = to_physical_state(opt, [var])
        assert np.allclose(back, physical)

    def test_state_var_round_trip_log(self):
        var = StateVar("dopant_concentration", 1.0, 1e17, 1e20, log_scale=True)
        physical = np.array([1e19])
        opt = to_optimizer_state(physical, [var])
        # log10(1e19) * 1.0 = 19.0
        assert np.isclose(opt[0], 19.0)
        back = to_physical_state(opt, [var])
        assert np.isclose(back[0], physical[0])

    def test_log_scale_rejects_nonpositive(self):
        var = StateVar("x", 1.0, 1e-3, 1.0, log_scale=True)
        with pytest.raises(ValueError, match="log_scale"):
            to_optimizer_state(np.array([0.0]), [var])

    def test_physical_bounds_default(self):
        c = _epitaxy_default()
        lower, upper = physical_bounds(c.optimization_state_vars(), None)
        assert lower[0] == 10e-6  # min_l
        assert upper[0] == 3e-3  # max_l

    def test_physical_bounds_override(self):
        c = _epitaxy_default()
        lower, upper = physical_bounds(
            c.optimization_state_vars(),
            {"min_l": 50e-6, "max_v_bridge": 5.0},
        )
        assert lower[0] == 50e-6
        assert upper[4] == 5.0  # v_bridge upper

    def test_physical_bounds_unknown_key(self):
        c = _epitaxy_default()
        with pytest.raises(ValueError, match="not match any state variable"):
            physical_bounds(c.optimization_state_vars(), {"min_nonexistent": 1.0})

    def test_physical_bounds_bad_prefix(self):
        c = _epitaxy_default()
        with pytest.raises(ValueError, match="must start with"):
            physical_bounds(c.optimization_state_vars(), {"length": 1.0})

    def test_physical_bounds_lower_above_upper(self):
        c = _epitaxy_default()
        with pytest.raises(ValueError, match="lower > upper"):
            physical_bounds(
                c.optimization_state_vars(),
                {"min_l": 1e-3, "max_l": 1e-4},
            )

    def test_random_physical_state_inside_bounds(self):
        c = _epitaxy_default()
        rng = np.random.default_rng(0)
        state_vars = c.optimization_state_vars()
        lower, upper = physical_bounds(state_vars, None)
        for _ in range(10):
            sample = random_physical_state(state_vars, None, rng)
            assert np.all(sample >= lower)
            assert np.all(sample <= upper)

    def test_random_physical_state_log_scale_distribution(self):
        c = _epitaxy_default()
        rng = np.random.default_rng(0)
        state_vars = c.optimization_state_vars()
        # Sample many times; the geometric mean of dopant_concentration
        # should be near 10**((17+log10(2e20))/2) ~= 10**18.65 = 4.5e18
        # for boron (max 2e20)
        samples = np.array([random_physical_state(state_vars, None, rng) for _ in range(500)])
        dopant = samples[:, 5]
        log_mean = np.mean(np.log10(dopant))
        expected_log_mean = (np.log10(1e17) + np.log10(2e20)) / 2
        assert abs(log_mean - expected_log_mean) < 0.2


# -----------------------------------------------------------------------------
# Per-subclass state-var declarations
# -----------------------------------------------------------------------------


class TestStateVarDeclarations:
    @pytest.mark.parametrize(
        "factory,expected_count,expected_first",
        [
            (_epitaxy_default, 7, "l"),
            (_diffusion_default, 7, "l"),
            (_implantation_default, 9, "l"),
            (_poly_default, 8, "l"),
            (_piezoelectric_default, 5, "l_si"),
        ],
    )
    def test_subclass_state_vars(self, factory, expected_count, expected_first):
        c = factory()
        state_vars = c.optimization_state_vars()
        assert len(state_vars) == expected_count
        assert state_vars[0].name == expected_first
        # All declared attributes exist on the cantilever
        for var in state_vars:
            assert hasattr(c, var.name), f"{type(c).__name__} missing attribute {var.name!r}"
        # All bounds finite and lower < upper (or equal for pinned)
        for var in state_vars:
            assert np.isfinite(var.default_min)
            assert np.isfinite(var.default_max)
            assert var.default_min <= var.default_max

    def test_epitaxy_concentration_bound_depends_on_dopant(self):
        c = _epitaxy_default()
        c.doping_type = "boron"
        boron_max = next(v for v in c.optimization_state_vars() if v.name == "dopant_concentration").default_max
        c.doping_type = "arsenic"
        arsenic_max = next(v for v in c.optimization_state_vars() if v.name == "dopant_concentration").default_max
        assert arsenic_max > boron_max  # 8e20 > 2e20 per MATLAB defaults

    def test_implantation_bounds_depend_on_lookup_source(self):
        c_tsuprem = CantileverImplantation(
            freq_min=10,
            freq_max=1000,
            l=200e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.3,
            v_bridge=2.0,
            doping_type="boron",
            annealing_time=30 * 60,
            annealing_temp=273.15 + 1000,
            annealing_type="inert",
            implantation_energy=40,
            implantation_dose=1e15,
            lookup_source="tsuprem4",
        )
        c_dope = CantileverImplantation(
            freq_min=10,
            freq_max=1000,
            l=200e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.3,
            v_bridge=2.0,
            doping_type="boron",
            annealing_time=30 * 60,
            annealing_temp=273.15 + 1000,
            annealing_type="inert",
            implantation_energy=40,
            implantation_dose=1e15,
            lookup_source="dopedealer",
        )
        ts = next(v for v in c_tsuprem.optimization_state_vars() if v.name == "implantation_energy")
        dd = next(v for v in c_dope.optimization_state_vars() if v.name == "implantation_energy")
        # dopedealer table covers a wider energy range
        assert dd.default_min < ts.default_min
        assert dd.default_max > ts.default_max


# -----------------------------------------------------------------------------
# Goal factories
# -----------------------------------------------------------------------------


class TestGoalFactories:
    def test_force_resolution_goal_units(self):
        c = _epitaxy_default()
        goal = force_resolution_goal()
        # MATLAB-style units: pN
        assert np.isclose(goal(c), c.force_resolution() * 1e12)

    def test_displacement_resolution_goal_units(self):
        c = _epitaxy_default()
        goal = displacement_resolution_goal()
        assert np.isclose(goal(c), c.displacement_resolution() * 1e9)

    def test_surface_stress_goal_units(self):
        c = _epitaxy_default()
        goal = surface_stress_resolution_goal()
        assert np.isclose(goal(c), c.surface_stress_resolution() * 1e6)

    def test_piezoelectric_voltage_goal(self):
        c = _piezoelectric_default()
        goal = voltage_force_resolution_goal()
        assert np.isclose(goal(c), c.Fminv() * 1e12)

    def test_piezoelectric_charge_goal(self):
        c = _piezoelectric_default()
        goal = charge_force_resolution_goal()
        assert np.isclose(goal(c), c.Fminq() * 1e12)


# -----------------------------------------------------------------------------
# Constraint validation
# -----------------------------------------------------------------------------


class TestConstraintValidation:
    def test_metric_constraint_requires_min_or_max(self):
        with pytest.raises(ValueError, match="at least one"):
            CantileverMetricConstraint(CantileverMetric.STIFFNESS)

    def test_metric_constraint_min_above_max(self):
        with pytest.raises(ValueError, match="<= maximum"):
            CantileverMetricConstraint(CantileverMetric.STIFFNESS, minimum=10.0, maximum=1.0)

    def test_metric_constraint_rejects_inf(self):
        with pytest.raises(ValueError, match="finite"):
            CantileverMetricConstraint(CantileverMetric.STIFFNESS, minimum=float("inf"))


# -----------------------------------------------------------------------------
# End-to-end optimization smoke tests
# -----------------------------------------------------------------------------


class TestOptimizeFromCurrent:
    def test_epitaxy_improves(self):
        c = _epitaxy_default()
        before = c.force_resolution()
        result = optimize_performance_from_current(
            c,
            force_resolution_goal(),
            metric_constraints=[
                CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
            ],
        )
        assert isinstance(result, OptimizationResult)
        assert result.optimized is not c
        # Original cantilever untouched
        assert np.isclose(c.force_resolution(), before)
        assert result.optimized.force_resolution() <= before

    def test_diffusion_runs(self):
        c = _diffusion_default()
        result = optimize_performance_from_current(
            c,
            force_resolution_goal(),
            metric_constraints=[
                CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
            ],
        )
        assert np.isfinite(result.objective_value)
        assert result.optimized.diffusion_time > 0

    def test_implantation_runs(self):
        c = _implantation_default()
        result = optimize_performance_from_current(
            c,
            force_resolution_goal(),
            metric_constraints=[
                CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
            ],
        )
        assert np.isfinite(result.objective_value)
        assert result.optimized.implantation_dose > 0

    def test_poly_runs(self):
        c = _poly_default()
        result = optimize_performance_from_current(
            c,
            force_resolution_goal(),
            metric_constraints=[
                CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
            ],
        )
        assert np.isfinite(result.objective_value)
        # Symmetric layout enforces t_bot == t_top via optimization_post_apply
        assert result.optimized.t_bot == result.optimized.t_top

    def test_piezoelectric_runs(self):
        c = _piezoelectric_default()
        result = optimize_performance_from_current(
            c,
            voltage_force_resolution_goal(),
            metric_constraints=[
                CantileverMetricConstraint(CantileverMetric.STIFFNESS, minimum=1e-3, maximum=1e2),
            ],
        )
        assert np.isfinite(result.objective_value)
        # l_pe / w_pe must track l_si / w_si via optimization_post_apply
        assert result.optimized.l_pe == result.optimized.l_si
        assert result.optimized.w_pe == result.optimized.w_si


# -----------------------------------------------------------------------------
# Multi-start + reproducibility
# -----------------------------------------------------------------------------


class TestOptimizePerformance:
    def test_multi_start_reproducible_with_seed(self):
        c = _epitaxy_default()
        kwargs = {
            "cantilever": c,
            "objective": force_resolution_goal(),
            "metric_constraints": [
                CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
            ],
            "n_starts": 2,
            "max_iterations": 2,
            "random_seed": 42,
        }
        result_a = optimize_performance(**kwargs)
        result_b = optimize_performance(**kwargs)
        assert np.isclose(result_a.objective_value, result_b.objective_value)
        assert np.allclose(result_a.physical_state, result_b.physical_state)

    def test_multi_start_finds_better_than_single_start(self):
        c = _epitaxy_default()
        constraints = [
            CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
            CantileverMetricConstraint(CantileverMetric.OMEGA_VACUUM_HZ, minimum=5000),
        ]
        single = optimize_performance_from_current(c, force_resolution_goal(), metric_constraints=constraints)
        multi = optimize_performance(
            c,
            force_resolution_goal(),
            metric_constraints=constraints,
            n_starts=4,
            max_iterations=8,
            random_seed=0,
        )
        # Multi-start should be at least as good as single-start (or close)
        # within a generous tolerance since SLSQP is local and may converge
        # to comparable basins.
        assert multi.objective_value <= single.objective_value * 1.5

    def test_n_starts_must_be_at_least_one(self):
        c = _epitaxy_default()
        with pytest.raises(ValueError, match="n_starts"):
            optimize_performance(c, force_resolution_goal(), n_starts=0)

    def test_max_iterations_below_n_starts_rejected(self):
        c = _epitaxy_default()
        with pytest.raises(ValueError, match="max_iterations"):
            optimize_performance(
                c,
                force_resolution_goal(),
                n_starts=5,
                max_iterations=2,
            )


# -----------------------------------------------------------------------------
# Constraint enforcement
# -----------------------------------------------------------------------------


class TestConstraintEnforcement:
    def test_power_dissipation_respected(self):
        c = _epitaxy_default()
        result = optimize_performance_from_current(
            c,
            force_resolution_goal(),
            metric_constraints=[
                CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=1e-3),
            ],
        )
        # SLSQP can violate constraints by ~1e-6 of scale; allow 1% tolerance
        assert result.optimized.power_dissipation() <= 1e-3 * 1.01

    def test_stiffness_minimum_respected(self):
        c = _epitaxy_default()
        result = optimize_performance_from_current(
            c,
            force_resolution_goal(),
            metric_constraints=[
                CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
                CantileverMetricConstraint(CantileverMetric.STIFFNESS, minimum=1e-3),
            ],
        )
        assert result.optimized.stiffness() >= 1e-3 * 0.99

    def test_aspect_ratio_default_constraints_respected(self):
        c = _epitaxy_default()
        result = optimize_performance_from_current(
            c,
            force_resolution_goal(),
            metric_constraints=[
                CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
            ],
        )
        opt = result.optimized
        assert opt.l / opt.w >= 2.0 * 0.99
        assert opt.w / opt.t >= 2.0 * 0.99
        assert opt.l_pr() / opt.w_pr() >= 2.0 * 0.99
        assert opt.l_pr() >= 2e-6 * 0.99

    def test_disabling_aspect_constraints(self):
        c = _epitaxy_default()
        result = optimize_performance_from_current(
            c,
            force_resolution_goal(),
            metric_constraints=[
                CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
            ],
            default_aspect_constraints=False,
        )
        assert np.isfinite(result.objective_value)

    def test_parameter_constraint_override(self):
        c = _epitaxy_default()
        result = optimize_performance_from_current(
            c,
            force_resolution_goal(),
            parameter_constraints={"max_v_bridge": 3.0},
            metric_constraints=[
                CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
            ],
        )
        assert result.optimized.v_bridge <= 3.0 + 1e-9


# -----------------------------------------------------------------------------
# Custom objectives
# -----------------------------------------------------------------------------


class TestCustomObjective:
    def test_lambda_objective_runs(self):
        c = _epitaxy_default()
        result = optimize_performance_from_current(
            c,
            lambda c: c.force_resolution() * 1e12 + c.power_dissipation() * 1e3,
        )
        assert np.isfinite(result.objective_value)

    def test_l_bfgs_b_with_no_constraints(self):
        c = _epitaxy_default()
        # Disable defaults so we can use unconstrained L-BFGS-B
        result = optimize_performance_from_current(
            c,
            force_resolution_goal(),
            default_aspect_constraints=False,
            method="L-BFGS-B",
        )
        assert np.isfinite(result.objective_value)

    def test_l_bfgs_b_rejected_with_metric_constraints(self):
        c = _epitaxy_default()
        with pytest.raises(ValueError, match="L-BFGS-B"):
            optimize_performance_from_current(
                c,
                force_resolution_goal(),
                metric_constraints=[
                    CantileverMetricConstraint(CantileverMetric.STIFFNESS, minimum=1e-3),
                ],
                method="L-BFGS-B",
            )


# -----------------------------------------------------------------------------
# Original cantilever immutability
# -----------------------------------------------------------------------------


def test_original_cantilever_not_mutated():
    c = _epitaxy_default()
    state_vars = c.optimization_state_vars()
    before = physical_state_from_cantilever(c, state_vars)
    optimize_performance_from_current(
        c,
        force_resolution_goal(),
        metric_constraints=[
            CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
        ],
    )
    after = physical_state_from_cantilever(c, state_vars)
    assert np.allclose(before, after)
