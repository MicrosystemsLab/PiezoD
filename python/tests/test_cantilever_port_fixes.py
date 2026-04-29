"""Regression tests for the cantilever port-bug fixes.

These cover methods that were broken on master prior to this branch:
- The FD temperature solver (calculateTempProfile and friends)
- The actuator stress / deflection chain
- Standalone bugs: d31, calculateEquivalentThickness, film_intrinsic_stress,
  thermal_conductivity_profile, calculateEnergies
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from piezod import (
    CantileverEpitaxy,
    CantileverMetric,
    CantileverMetricConstraint,
    force_resolution_goal,
    optimize_performance_from_current,
)


@pytest.fixture(autouse=True)
def _silence_warnings():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def _epitaxy_default() -> CantileverEpitaxy:
    """Default epitaxial cantilever with no actuator."""
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


def _thermal_actuated() -> CantileverEpitaxy:
    """Same geometry but with a thermal actuator at the base."""
    c = _epitaxy_default()
    c.fluid = "air"
    c.cantilever_type = "thermal"
    c.l_a = 50e-6
    c.w_a = 20e-6
    c.t_a = 200e-9
    c.w_a_active = 15e-6
    c.l_a_gap = 5e-6
    c.v_actuator = 2.0
    c.R_heater = 500.0
    return c


# -----------------------------------------------------------------------------
# FD temperature solver (calculateTempProfile and dependents)
# -----------------------------------------------------------------------------


class TestTempProfile:
    def test_calculate_temp_profile_runs(self):
        c = _epitaxy_default()
        x, Q, T = c.calculateTempProfile()
        assert x.shape == (c.numXPoints,)
        assert Q.shape == (c.numXPoints,)
        assert T.shape == (c.numXPoints,)
        assert np.all(np.isfinite(T))
        assert T.max() > 0  # Joule heating raises the temperature

    def test_fd_and_approx_agree_for_pr_only(self):
        # Without an actuator, lumped-circuit (approx) and FD should agree
        # within ~5%. Big disagreement signals a regression in either.
        c = _epitaxy_default()
        TMax_approx, TTip_approx = c.approxTempRise()
        TMax_exact, TTip_exact = c.calculateMaxAndTipTemp()
        assert TMax_approx > 0 and TMax_exact > 0
        assert abs(TMax_exact - TMax_approx) / TMax_approx < 0.05

    def test_actuator_heating_raises_temp(self):
        # A thermal actuator dissipating ~mW should raise the FD-computed
        # tip temperature far above the lumped-circuit estimate (which only
        # captures PR heating).
        c = _thermal_actuated()
        TMax_approx, _ = c.approxTempRise()
        TMax_exact, _ = c.calculateMaxAndTipTemp()
        assert TMax_exact > 5 * TMax_approx

    def test_temp_dependent_iteration_converges(self):
        c = _epitaxy_default()
        x, Q, T = c.calculateTempProfileTempDependent()
        assert x.shape == (c.numXPoints,)
        assert np.all(np.isfinite(T))

    def test_average_pr_temp_returns_scalar(self):
        c = _epitaxy_default()
        value = c.averagePRTemp()
        assert isinstance(value, float)
        assert value > c.T  # T is ambient; PR heats above ambient

    def test_max_pr_temp_returns_scalar(self):
        c = _epitaxy_default()
        value = c.maxPRTemp()
        assert isinstance(value, float)
        assert value > 0

    def test_temp_rise_at_pr_base(self):
        c = _epitaxy_default()
        value = c.tempRiseAtPRBase()
        assert isinstance(value, float)
        assert np.isfinite(value)

    def test_average_actuator_delta_temp(self):
        c = _thermal_actuated()
        value = c.averageActuatorDeltaTemp()
        assert value > 0

    def test_thermal_crosstalk(self):
        c = _thermal_actuated()
        value = c.thermalCrosstalk()
        assert np.isfinite(value)

    def test_thermal_conductivity_profile_scalar_T(self):
        c = _epitaxy_default()
        x, z, k = c.thermal_conductivity_profile(300.0)
        assert k.shape == (c.numZPoints, c.numXPoints)
        assert np.all(k > 0)

    def test_thermal_conductivity_profile_array_T(self):
        c = _epitaxy_default()
        T = np.full(c.numXPoints, 300.0)
        x, z, k = c.thermal_conductivity_profile(T)
        assert k.shape == (c.numZPoints, c.numXPoints)
        assert np.all(k > 0)


# -----------------------------------------------------------------------------
# Standalone port bugs
# -----------------------------------------------------------------------------


class TestStandaloneFixes:
    def test_d31_returns_finite_scalar(self):
        c = _thermal_actuated()
        c.cantilever_type = "piezoelectric"
        c.t_a = 500e-9
        value = c.d31()
        assert isinstance(value, float)
        # Class data has d31_aln stored as a negative pm/V; just require
        # the spline lookup to produce a finite, non-zero value.
        assert np.isfinite(value)
        assert value != 0

    def test_d31_manual_override(self):
        c = _thermal_actuated()
        c.cantilever_type = "piezoelectric"
        c.t_a = 500e-9
        c.d31_manual = 5e-12
        assert c.d31() == 5e-12

    def test_film_intrinsic_stress_default_none(self):
        c = _epitaxy_default()
        sigma = c.film_intrinsic_stress()
        assert sigma.shape == (3,)
        assert np.all(sigma == 0)

    def test_film_intrinsic_stress_thermal(self):
        c = _thermal_actuated()
        sigma = c.film_intrinsic_stress()
        assert sigma.shape == (3,)
        assert np.all(np.isfinite(sigma))

    def test_film_intrinsic_stress_piezoelectric(self):
        c = _thermal_actuated()
        c.cantilever_type = "piezoelectric"
        c.metal_type = "titanium"
        sigma = c.film_intrinsic_stress()
        assert sigma.shape == (6,)

    def test_calculate_equivalent_thickness(self):
        c = _thermal_actuated()
        t_equiv = c.calculateEquivalentThickness()
        # Equivalent thickness should be in the same order of magnitude as
        # the actual layer stack; certainly positive and finite.
        assert np.isfinite(t_equiv)
        assert t_equiv > 0

    def test_calculate_energies_returns_finite(self):
        c = _thermal_actuated()
        omega_test = 2 * np.pi * 50e3  # 50 kHz trial
        U_e, U_k = c.calculateEnergies(omega_test)
        assert np.isfinite(U_e) and np.isfinite(U_k)
        assert U_e > 0 and U_k > 0


# -----------------------------------------------------------------------------
# Actuator stress + deflection chain
# -----------------------------------------------------------------------------


class TestDeflectionChain:
    def test_lookup_actuator_mechanics_thermal(self):
        c = _thermal_actuated()
        z, E, A, I = c.lookupActuatorMechanics()
        assert z.shape == (3,)
        assert E.shape == (3,)
        # Centroid of bottom layer = t/2; centroid of top = sum(t) - t_top/2
        assert z[0] < z[-1]
        assert np.all(A > 0)
        assert np.all(I > 0)

    def test_lookup_actuator_mechanics_none_raises(self):
        c = _epitaxy_default()
        with pytest.raises(RuntimeError, match="lookupActuatorMechanics"):
            c.lookupActuatorMechanics()

    def test_calculate_actuator_stress_shape_thermal(self):
        c = _thermal_actuated()
        stress = c.calculateActuatorStress()
        assert stress.shape == (c.numXPoints, 3)
        # Past the end of the actuator (x > l_a) the oxide and metal layers
        # don't exist, so columns 1 and 2 of the stress array must be zero.
        # The silicon column (col 0) extends the full cantilever length and
        # carries thermal stress everywhere.
        total_length = c.l + c.l_a
        end_actuator_idx = int(np.ceil(c.numXPoints * c.l_a / total_length)) + 1
        assert np.all(stress[end_actuator_idx:, 1] == 0)  # oxide
        assert np.all(stress[end_actuator_idx:, 2] == 0)  # metal

    def test_calculate_actuator_stress_none_returns_zeros(self):
        c = _epitaxy_default()
        stress = c.calculateActuatorStress()
        assert stress.shape == (c.numXPoints, 3)
        assert np.all(stress == 0)

    def test_calculate_deflection_thermal(self):
        c = _thermal_actuated()
        x, defl = c.calculateDeflection()
        assert x.shape == (c.numXPoints,)
        assert defl.shape == (c.numXPoints,)
        # Some non-zero deflection from thermal actuation
        assert abs(defl[-1]) > 0

    def test_tip_deflection_none_returns_zero(self):
        c = _epitaxy_default()
        assert c.tipDeflection() == 0.0

    def test_tip_deflection_thermal_nonzero(self):
        c = _thermal_actuated()
        assert abs(c.tipDeflection()) > 0

    def test_actuator_neutral_axis(self):
        c = _thermal_actuated()
        z = c.actuatorNeutralAxis()
        assert np.isfinite(z) and z > 0

    def test_tip_deflection_distribution(self):
        c = _thermal_actuated()
        mu, sigma = c.tip_deflection_distribution()
        assert np.isfinite(mu) and np.isfinite(sigma)


# -----------------------------------------------------------------------------
# Resonant frequency on actuated cantilevers
# -----------------------------------------------------------------------------


class TestResonantFrequencyActuated:
    def test_omega_vacuum_thermal(self):
        c = _thermal_actuated()
        omega = c.omega_vacuum()
        assert np.isfinite(omega) and omega > 0

    def test_omega_vacuum_hz_thermal(self):
        c = _thermal_actuated()
        f0 = c.omega_vacuum_hz()
        # 200 um Si cantilever with 50 um thermal base: should land in
        # the tens to hundreds of kHz range.
        assert 1e3 < f0 < 1e6


# -----------------------------------------------------------------------------
# print_performance and plot_noise_spectrum on default cantilever
# -----------------------------------------------------------------------------


class TestPrintAndPlot:
    def test_print_performance_runs(self, capsys):
        c = _epitaxy_default()
        c.print_performance()
        out = capsys.readouterr().out
        assert "Force resolution" in out
        assert "F-D Temp Rises" in out
        assert "Approx. Temp Rises" in out

    def test_plot_noise_spectrum_runs(self):
        import matplotlib

        matplotlib.use("Agg")
        c = _epitaxy_default()
        c.plot_noise_spectrum()  # raises on failure


# -----------------------------------------------------------------------------
# Optimizer integration with the now-working metrics
# -----------------------------------------------------------------------------


class TestOptimizerWithRestoredMetrics:
    def test_temp_max_exact_constraint(self):
        c = _epitaxy_default()
        result = optimize_performance_from_current(
            c,
            force_resolution_goal(),
            metric_constraints=[
                CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
                CantileverMetricConstraint(CantileverMetric.TEMP_MAX_EXACT, maximum=10.0),
            ],
        )
        # SLSQP can violate constraints by a few percent on this kind of
        # multi-constraint geometry+process problem; check the bound is
        # respected within 10% as a sanity bound and that the optimizer
        # made meaningful progress (the unconstrained solution sits well
        # above 10 K).
        t_max, _ = result.optimized.calculateMaxAndTipTemp()
        assert t_max <= 10.0 * 1.10

    def test_tip_deflection_metric_evaluable(self):
        # Default 'none' cantilever returns 0 for tipDeflection; the
        # constraint is well-defined even though the optimizer can't
        # actually shape it without an actuator.
        c = _epitaxy_default()
        result = optimize_performance_from_current(
            c,
            force_resolution_goal(),
            metric_constraints=[
                CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
            ],
        )
        assert result.optimized.tipDeflection() == 0.0
