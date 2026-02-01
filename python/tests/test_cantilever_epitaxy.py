"""Tests for CantileverEpitaxy class."""

import numpy as np
import pytest

from piezod.cantilever_epitaxy import CantileverEpitaxy


class TestCantileverEpitaxyInitialization:
    """Test initialization and basic properties."""

    def test_default_initialization(self) -> None:
        """Test cantilever initializes with default values."""
        cant = CantileverEpitaxy()

        assert cant.dopant_concentration == 1e19
        assert cant.t_pr_ratio == 0.3
        assert cant.doping_type == "phosphorus"
        assert cant.freq_min == 1.0
        assert cant.freq_max == 1e3
        assert cant.t == 1e-6

    def test_custom_initialization(self) -> None:
        """Test cantilever initializes with custom values."""
        cant = CantileverEpitaxy(
            freq_min=10.0,
            freq_max=1e4,
            l=200e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.4,
            v_bridge=2.0,
            doping_type="boron",
            dopant_concentration=5e19,
            t_pr_ratio=0.5,
        )

        assert cant.freq_min == 10.0
        assert cant.freq_max == 1e4
        assert cant.l == 200e-6
        assert cant.w == 20e-6
        assert cant.t == 2e-6
        assert cant.l_pr_ratio == 0.4
        assert cant.v_bridge == 2.0
        assert cant.doping_type == "boron"
        assert cant.dopant_concentration == 5e19
        assert cant.t_pr_ratio == 0.5

    def test_inherits_from_cantilever(self) -> None:
        """Test that CantileverEpitaxy properly inherits from Cantilever."""
        from piezod.cantilever import Cantilever

        cant = CantileverEpitaxy()
        assert isinstance(cant, Cantilever)
        assert hasattr(cant, "k_b")
        assert hasattr(cant, "numZPoints")


class TestJunctionDepth:
    """Test junction depth property."""

    def test_junction_depth_calculation(self) -> None:
        """Test junction depth is correctly calculated."""
        cant = CantileverEpitaxy(t=1e-6, t_pr_ratio=0.3)
        assert cant.junction_depth == pytest.approx(0.3e-6)

    def test_junction_depth_varies_with_t(self) -> None:
        """Test junction depth changes with thickness."""
        cant = CantileverEpitaxy(t=2e-6, t_pr_ratio=0.3)
        assert cant.junction_depth == pytest.approx(0.6e-6)

    def test_junction_depth_varies_with_ratio(self) -> None:
        """Test junction depth changes with t_pr_ratio."""
        cant = CantileverEpitaxy(t=1e-6, t_pr_ratio=0.5)
        assert cant.junction_depth == pytest.approx(0.5e-6)

    def test_junction_depth_is_property(self) -> None:
        """Test that junction_depth is a property, not a method."""
        cant = CantileverEpitaxy()
        # Should be accessible without calling as a function
        depth = cant.junction_depth
        assert isinstance(depth, float)


class TestDopingProfile:
    """Test doping profile calculation."""

    def test_doping_profile_returns_correct_types(self) -> None:
        """Test doping_profile returns three numpy arrays."""
        cant = CantileverEpitaxy()
        z, active, total = cant.doping_profile()

        assert isinstance(z, np.ndarray)
        assert isinstance(active, np.ndarray)
        assert isinstance(total, np.ndarray)

    def test_doping_profile_correct_shape(self) -> None:
        """Test all arrays have length numZPoints."""
        cant = CantileverEpitaxy()
        z, active, total = cant.doping_profile()

        assert len(z) == cant.numZPoints
        assert len(active) == cant.numZPoints
        assert len(total) == cant.numZPoints

    def test_doping_profile_range(self) -> None:
        """Test z ranges from 0 to t."""
        cant = CantileverEpitaxy(t=1e-6)
        z, _, _ = cant.doping_profile()

        assert z[0] == pytest.approx(0.0)
        assert z[-1] == pytest.approx(1e-6)

    def test_doping_profile_step_function(self) -> None:
        """Test profile is a step function at junction depth."""
        cant = CantileverEpitaxy(t=1e-6, t_pr_ratio=0.3, dopant_concentration=1e19)
        z, active, total = cant.doping_profile()

        # Points below junction should have epitaxial concentration
        epitaxial_points = z <= cant.junction_depth
        assert np.all(active[epitaxial_points] == 1e19)

        # Points above junction should have background concentration
        substrate_points = z > cant.junction_depth
        if np.any(substrate_points):
            assert np.all(active[substrate_points] == 1e15)

    def test_doping_profile_active_equals_total(self) -> None:
        """Test active doping equals total doping for epitaxy."""
        cant = CantileverEpitaxy()
        _, active, total = cant.doping_profile()

        np.testing.assert_array_equal(active, total)

    def test_doping_profile_background_concentration(self) -> None:
        """Test background concentration is 1e15."""
        cant = CantileverEpitaxy(t=1e-6, t_pr_ratio=0.1)
        z, active, _ = cant.doping_profile()

        # Last point should be in substrate
        assert active[-1] == pytest.approx(1e15)


class TestSheetResistance:
    """Test sheet resistance calculation."""

    def test_sheet_resistance_returns_float(self) -> None:
        """Test sheet_resistance returns a float."""
        cant = CantileverEpitaxy()
        Rs = cant.sheet_resistance()
        assert isinstance(Rs, (float, np.floating))

    def test_sheet_resistance_positive(self) -> None:
        """Test sheet resistance is positive."""
        cant = CantileverEpitaxy()
        Rs = cant.sheet_resistance()
        assert Rs > 0

    def test_sheet_resistance_decreases_with_concentration(self) -> None:
        """Test sheet resistance decreases with higher doping."""
        cant1 = CantileverEpitaxy(dopant_concentration=1e18)
        cant2 = CantileverEpitaxy(dopant_concentration=1e19)

        Rs1 = cant1.sheet_resistance()
        Rs2 = cant2.sheet_resistance()

        assert Rs2 < Rs1

    def test_sheet_resistance_decreases_with_thickness(self) -> None:
        """Test sheet resistance decreases with thicker epitaxial layer."""
        cant1 = CantileverEpitaxy(t=1e-6, t_pr_ratio=0.2)
        cant2 = CantileverEpitaxy(t=1e-6, t_pr_ratio=0.4)

        Rs1 = cant1.sheet_resistance()
        Rs2 = cant2.sheet_resistance()

        assert Rs2 < Rs1


class TestNz:
    """Test integrated carrier concentration."""

    def test_Nz_returns_float(self) -> None:
        """Test Nz returns a float."""
        cant = CantileverEpitaxy()
        Nz = cant.Nz()
        assert isinstance(Nz, (float, np.floating))

    def test_Nz_positive(self) -> None:
        """Test Nz is positive."""
        cant = CantileverEpitaxy()
        Nz = cant.Nz()
        assert Nz > 0

    def test_Nz_calculation(self) -> None:
        """Test Nz calculation matches expected formula."""
        cant = CantileverEpitaxy(t=1e-6, t_pr_ratio=0.3, dopant_concentration=1e19)
        Nz = cant.Nz()

        # Nz = junction_depth * dopant_concentration * 1e6
        # junction_depth = 1e-6 * 0.3 = 3e-7 m
        # Nz = 3e-7 * 1e19 * 1e6 = 3e18 carriers/m^2
        expected = 0.3e-6 * 1e19 * 1e6
        assert Nz == pytest.approx(expected)

    def test_Nz_scales_with_concentration(self) -> None:
        """Test Nz scales linearly with dopant concentration."""
        cant1 = CantileverEpitaxy(dopant_concentration=1e18)
        cant2 = CantileverEpitaxy(dopant_concentration=2e18)

        Nz1 = cant1.Nz()
        Nz2 = cant2.Nz()

        assert Nz2 == pytest.approx(2 * Nz1)


class TestAlpha:
    """Test Hooge parameter."""

    def test_alpha_returns_default(self) -> None:
        """Test alpha returns the default value."""
        cant = CantileverEpitaxy()
        alpha = cant.alpha()
        assert alpha == cant.default_alpha

    def test_alpha_value(self) -> None:
        """Test alpha has expected value of 1e-5."""
        cant = CantileverEpitaxy()
        alpha = cant.alpha()
        assert alpha == pytest.approx(1e-5)


class TestOptimizationScaling:
    """Test optimization scaling factors."""

    def test_scaling_returns_array(self) -> None:
        """Test scaling returns numpy array."""
        cant = CantileverEpitaxy()
        scaling = cant.doping_optimization_scaling()
        assert isinstance(scaling, np.ndarray)

    def test_scaling_correct_length(self) -> None:
        """Test scaling array has length 2."""
        cant = CantileverEpitaxy()
        scaling = cant.doping_optimization_scaling()
        assert len(scaling) == 2

    def test_scaling_values(self) -> None:
        """Test scaling values match expected values."""
        cant = CantileverEpitaxy()
        scaling = cant.doping_optimization_scaling()

        assert scaling[0] == pytest.approx(1e-19)
        assert scaling[1] == pytest.approx(10.0)


class TestOptimizationStateMethods:
    """Test optimization state conversion methods."""

    def test_doping_current_state(self) -> None:
        """Test current state extraction."""
        cant = CantileverEpitaxy(dopant_concentration=5e19, t_pr_ratio=0.4)
        state = cant.doping_current_state()

        assert isinstance(state, np.ndarray)
        assert len(state) == 2
        assert state[0] == pytest.approx(5e19)
        assert state[1] == pytest.approx(0.4)

    def test_doping_cantilever_from_state(self) -> None:
        """Test state reconstruction."""
        cant = CantileverEpitaxy()
        # Create a state vector with doping params at indices 5-6
        x0 = np.array([0, 0, 0, 0, 0, 7e19, 0.6])

        cant.doping_cantilever_from_state(x0)

        assert cant.dopant_concentration == pytest.approx(7e19)
        assert cant.t_pr_ratio == pytest.approx(0.6)


class TestOptimizationBounds:
    """Test optimization bounds for different dopant types."""

    def test_bounds_return_tuple_of_arrays(self) -> None:
        """Test bounds returns tuple of two arrays."""
        cant = CantileverEpitaxy()
        lb, ub = cant.doping_optimization_bounds()

        assert isinstance(lb, np.ndarray)
        assert isinstance(ub, np.ndarray)
        assert len(lb) == 2
        assert len(ub) == 2

    def test_bounds_lower_less_than_upper(self) -> None:
        """Test lower bounds are less than upper bounds."""
        cant = CantileverEpitaxy()
        lb, ub = cant.doping_optimization_bounds()

        assert np.all(lb < ub)

    def test_boron_bounds(self) -> None:
        """Test boron has correct solid solubility limit."""
        cant = CantileverEpitaxy(doping_type="boron")
        lb, ub = cant.doping_optimization_bounds()

        assert lb[0] == pytest.approx(1e17)
        assert ub[0] == pytest.approx(2e20)
        assert lb[1] == pytest.approx(0.01)
        assert ub[1] == pytest.approx(0.99)

    def test_phosphorus_bounds(self) -> None:
        """Test phosphorus has correct solid solubility limit."""
        cant = CantileverEpitaxy(doping_type="phosphorus")
        lb, ub = cant.doping_optimization_bounds()

        assert lb[0] == pytest.approx(1e17)
        assert ub[0] == pytest.approx(4e20)

    def test_arsenic_bounds(self) -> None:
        """Test arsenic has correct solid solubility limit."""
        cant = CantileverEpitaxy(doping_type="arsenic")
        lb, ub = cant.doping_optimization_bounds()

        assert lb[0] == pytest.approx(1e17)
        assert ub[0] == pytest.approx(8e20)

    def test_custom_constraints_override_defaults(self) -> None:
        """Test custom constraints override default bounds."""
        cant = CantileverEpitaxy()
        constraints = {
            "min_dopant_concentration": 5e17,
            "max_dopant_concentration": 1e20,
            "min_t_pr_ratio": 0.1,
            "max_t_pr_ratio": 0.8,
        }

        lb, ub = cant.doping_optimization_bounds(constraints)

        assert lb[0] == pytest.approx(5e17)
        assert ub[0] == pytest.approx(1e20)
        assert lb[1] == pytest.approx(0.1)
        assert ub[1] == pytest.approx(0.8)

    def test_partial_constraints(self) -> None:
        """Test partial constraint override."""
        cant = CantileverEpitaxy(doping_type="boron")
        constraints = {"max_t_pr_ratio": 0.5}

        lb, ub = cant.doping_optimization_bounds(constraints)

        # Doping bounds should be defaults
        assert lb[0] == pytest.approx(1e17)
        assert ub[0] == pytest.approx(2e20)
        # t_pr_ratio upper bound should be custom
        assert lb[1] == pytest.approx(0.01)
        assert ub[1] == pytest.approx(0.5)

    def test_unknown_dopant_type(self) -> None:
        """Test fallback for unknown dopant type."""
        cant = CantileverEpitaxy()
        cant.doping_type = "unknown"

        lb, ub = cant.doping_optimization_bounds()

        # Should use fallback value
        assert ub[0] == pytest.approx(1e20)


class TestRandomInitialConditions:
    """Test random initial condition generation."""

    def test_returns_array(self) -> None:
        """Test returns numpy array of length 2."""
        cant = CantileverEpitaxy()
        x0 = cant.doping_initial_conditions_random()

        assert isinstance(x0, np.ndarray)
        assert len(x0) == 2

    def test_within_bounds(self) -> None:
        """Test random values are within bounds."""
        cant = CantileverEpitaxy()
        lb, ub = cant.doping_optimization_bounds()

        # Test multiple random samples
        for _ in range(10):
            x0 = cant.doping_initial_conditions_random()
            assert lb[0] <= x0[0] <= ub[0]
            assert lb[1] <= x0[1] <= ub[1]

    def test_within_custom_bounds(self) -> None:
        """Test random values respect custom constraints."""
        cant = CantileverEpitaxy()
        constraints = {
            "min_dopant_concentration": 5e18,
            "max_dopant_concentration": 1e19,
            "min_t_pr_ratio": 0.2,
            "max_t_pr_ratio": 0.6,
        }

        lb, ub = cant.doping_optimization_bounds(constraints)

        for _ in range(10):
            x0 = cant.doping_initial_conditions_random(constraints)
            assert lb[0] <= x0[0] <= ub[0]
            assert lb[1] <= x0[1] <= ub[1]

    def test_concentration_logarithmic_distribution(self) -> None:
        """Test dopant concentration is logarithmically distributed."""
        cant = CantileverEpitaxy()

        # Generate many samples
        samples = [cant.doping_initial_conditions_random()[0] for _ in range(100)]

        # Convert to log space
        log_samples = np.log10(samples)

        # In log space, should be approximately uniformly distributed
        # Check that we get samples across the full range
        lb, ub = cant.doping_optimization_bounds()
        log_min = np.log10(lb[0])
        log_max = np.log10(ub[0])

        # At least some samples in each quartile
        q1 = log_min + 0.25 * (log_max - log_min)
        q2 = log_min + 0.5 * (log_max - log_min)
        q3 = log_min + 0.75 * (log_max - log_min)

        assert any(log_samples < q1)
        assert any(q1 <= s < q2 for s in log_samples)
        assert any(q2 <= s < q3 for s in log_samples)
        assert any(s >= q3 for s in log_samples)

    def test_t_pr_ratio_uniform_distribution(self) -> None:
        """Test t_pr_ratio is uniformly distributed."""
        cant = CantileverEpitaxy()

        # Generate many samples
        samples = [cant.doping_initial_conditions_random()[1] for _ in range(100)]

        lb, ub = cant.doping_optimization_bounds()
        t_min = lb[1]
        t_max = ub[1]

        # Check samples span the range
        assert min(samples) < t_min + 0.2 * (t_max - t_min)
        assert max(samples) > t_max - 0.2 * (t_max - t_min)

    def test_randomness(self) -> None:
        """Test that multiple calls produce different values."""
        cant = CantileverEpitaxy()

        x0_1 = cant.doping_initial_conditions_random()
        x0_2 = cant.doping_initial_conditions_random()

        # Very unlikely to be exactly equal if truly random
        assert not np.array_equal(x0_1, x0_2)


class TestIntegrationWithBaseClass:
    """Test integration with base Cantilever class."""

    def test_can_call_conductivity(self) -> None:
        """Test can call inherited conductivity method."""
        cant = CantileverEpitaxy()
        conductivity = cant.conductivity(1e19)

        assert isinstance(conductivity, (float, np.floating))
        assert conductivity > 0

    def test_uses_correct_doping_type_in_conductivity(self) -> None:
        """Test conductivity uses the correct doping type."""
        cant_p = CantileverEpitaxy(doping_type="phosphorus")
        cant_b = CantileverEpitaxy(doping_type="boron")

        # At same concentration, different dopants have different conductivities
        cond_p = cant_p.conductivity(1e19)
        cond_b = cant_b.conductivity(1e19)

        # They should be different (though both positive)
        assert cond_p != cond_b
        assert cond_p > 0
        assert cond_b > 0
