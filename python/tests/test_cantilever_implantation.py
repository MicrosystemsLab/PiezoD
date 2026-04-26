"""Tests for CantileverImplantation class."""

import numpy as np
import pytest

from piezod.cantilever_implantation import (
    CantileverImplantation,
    DopingMetric,
    MetricConstraint,
)


class TestCantileverImplantationInitialization:
    """Test CantileverImplantation initialization."""

    @pytest.fixture
    def basic_params(self):
        """Common parameters for creating test cantilevers."""
        return {
            "freq_min": 1e3,
            "freq_max": 1e6,
            "l": 100e-6,
            "w": 20e-6,
            "t": 2e-6,
            "l_pr_ratio": 0.5,
            "v_bridge": 5.0,
            "doping_type": "boron",
            "annealing_time": 3600,  # 60 min
            "annealing_temp": 1173.15,  # 900 C
            "annealing_type": "inert",
            "implantation_energy": 50,
            "implantation_dose": 1e15,
        }

    def test_initialization_inert(self, basic_params):
        """Test initialization with inert annealing."""
        cantilever = CantileverImplantation(**basic_params)

        assert cantilever.implantation_energy == 50
        assert cantilever.implantation_dose == 1e15
        assert cantilever.annealing_type == "inert"
        assert cantilever.annealing_temp == 1173.15
        assert cantilever.annealing_time == 3600
        # Lookup table is auto-loaded
        assert cantilever._lookup_data is not None

    def test_initialization_oxide(self, basic_params):
        """Test initialization with oxide annealing."""
        basic_params["annealing_type"] = "oxide"
        cantilever = CantileverImplantation(**basic_params)

        assert cantilever.annealing_type == "oxide"

    def test_initialization_all_dopants(self, basic_params):
        """Test initialization with all dopant types."""
        for dopant in ["boron", "phosphorus", "arsenic"]:
            basic_params["doping_type"] = dopant
            cantilever = CantileverImplantation(**basic_params)
            assert cantilever.doping_type == dopant


class TestAnnealNumber:
    """Test anneal_number method."""

    @pytest.fixture
    def cantilever_inert(self):
        """Create cantilever with inert annealing."""
        return CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=3600,
            annealing_temp=1173.15,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

    @pytest.fixture
    def cantilever_oxide(self):
        """Create cantilever with oxide annealing."""
        return CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=3600,
            annealing_temp=1173.15,
            annealing_type="oxide",
            implantation_energy=50,
            implantation_dose=1e15,
        )

    def test_anneal_number_inert(self, cantilever_inert):
        """Test anneal_number returns 1 for inert."""
        assert cantilever_inert.anneal_number() == 1

    def test_anneal_number_oxide(self, cantilever_oxide):
        """Test anneal_number returns 2 for oxide."""
        assert cantilever_oxide.anneal_number() == 2

    def test_anneal_number_invalid(self, cantilever_inert):
        """Test anneal_number raises error for invalid type."""
        cantilever_inert.annealing_type = "invalid"
        with pytest.raises(ValueError, match="Unknown anneal condition"):
            cantilever_inert.anneal_number()


class TestDiffusionLength:
    """Test diffusion_length property."""

    @pytest.fixture
    def cantilever_params(self):
        """Common parameters for diffusion length tests."""
        return {
            "freq_min": 1e3,
            "freq_max": 1e6,
            "l": 100e-6,
            "w": 20e-6,
            "t": 2e-6,
            "l_pr_ratio": 0.5,
            "v_bridge": 5.0,
            "annealing_time": 3600,  # 60 min
            "annealing_temp": 1173.15,  # 900 C
            "annealing_type": "inert",
            "implantation_energy": 50,
            "implantation_dose": 1e15,
        }

    def test_diffusion_length_boron(self, cantilever_params):
        """Test diffusion length calculation for boron."""
        cantilever_params["doping_type"] = "boron"
        cantilever = CantileverImplantation(**cantilever_params)

        # Calculate expected value
        D0 = 0.76  # cm^2/s
        Ea = 3.46  # eV
        k_b_eV = 8.617343e-5  # eV/K
        T = 1173.15  # K
        t = 3600  # s
        D = D0 * np.exp(-Ea / k_b_eV / T)
        expected = np.sqrt(D * t)

        assert cantilever.diffusion_length == pytest.approx(expected, rel=1e-10)

    def test_diffusion_length_phosphorus(self, cantilever_params):
        """Test diffusion length calculation for phosphorus."""
        cantilever_params["doping_type"] = "phosphorus"
        cantilever = CantileverImplantation(**cantilever_params)

        # Calculate expected value
        D0 = 3.85  # cm^2/s
        Ea = 3.66  # eV
        k_b_eV = 8.617343e-5  # eV/K
        T = 1173.15  # K
        t = 3600  # s
        D = D0 * np.exp(-Ea / k_b_eV / T)
        expected = np.sqrt(D * t)

        assert cantilever.diffusion_length == pytest.approx(expected, rel=1e-10)

    def test_diffusion_length_arsenic(self, cantilever_params):
        """Test diffusion length calculation for arsenic."""
        cantilever_params["doping_type"] = "arsenic"
        cantilever = CantileverImplantation(**cantilever_params)

        # Calculate expected value
        D0 = 22.9  # cm^2/s
        Ea = 4.1  # eV
        k_b_eV = 8.617343e-5  # eV/K
        T = 1173.15  # K
        t = 3600  # s
        D = D0 * np.exp(-Ea / k_b_eV / T)
        expected = np.sqrt(D * t)

        assert cantilever.diffusion_length == pytest.approx(expected, rel=1e-10)

    def test_diffusion_length_temperature_dependence(self, cantilever_params):
        """Test diffusion length increases with temperature."""
        cantilever_params["doping_type"] = "boron"

        cantilever_low = CantileverImplantation(**cantilever_params)
        cantilever_params["annealing_temp"] = 1373.15  # Higher temp
        cantilever_high = CantileverImplantation(**cantilever_params)

        assert cantilever_high.diffusion_length > cantilever_low.diffusion_length

    def test_diffusion_length_time_dependence(self, cantilever_params):
        """Test diffusion length increases with time."""
        cantilever_params["doping_type"] = "boron"

        cantilever_short = CantileverImplantation(**cantilever_params)
        cantilever_params["annealing_time"] = 7200  # Double the time
        cantilever_long = CantileverImplantation(**cantilever_params)

        # Diffusion length should scale with sqrt(time)
        expected_ratio = np.sqrt(2)
        actual_ratio = cantilever_long.diffusion_length / cantilever_short.diffusion_length
        assert actual_ratio == pytest.approx(expected_ratio, rel=1e-10)

    def test_diffusion_length_invalid_dopant(self, cantilever_params):
        """Test diffusion length raises error for invalid dopant."""
        cantilever_params["doping_type"] = "boron"
        cantilever = CantileverImplantation(**cantilever_params)
        cantilever.doping_type = "invalid"

        with pytest.raises(ValueError, match="Unknown doping_type"):
            _ = cantilever.diffusion_length


class TestAlpha:
    """Test alpha (Hooge noise parameter) calculation."""

    @pytest.fixture
    def cantilever(self):
        """Create a test cantilever."""
        return CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=3600,
            annealing_temp=1173.15,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

    def test_alpha_calculation(self, cantilever):
        """Test alpha calculation from diffusion length."""
        diffusion_length = cantilever.diffusion_length
        expected_alpha = 2.469e-10 * diffusion_length**-0.598

        assert cantilever.alpha() == pytest.approx(expected_alpha, rel=1e-10)

    def test_alpha_positive(self, cantilever):
        """Test alpha is always positive."""
        assert cantilever.alpha() > 0

    def test_alpha_decreases_with_diffusion_length(self):
        """Test alpha decreases as diffusion length increases."""
        # Longer anneal time -> longer diffusion length -> smaller alpha
        cantilever_short = CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=900,  # 15 min
            annealing_temp=1173.15,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

        cantilever_long = CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=7200,  # 120 min
            annealing_temp=1173.15,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

        assert cantilever_short.alpha() > cantilever_long.alpha()


class TestOptimizationScaling:
    """Test doping_optimization_scaling method."""

    @pytest.fixture
    def cantilever(self):
        """Create a test cantilever."""
        return CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=3600,
            annealing_temp=1173.15,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

    def test_scaling_values(self, cantilever):
        """Test optimization scaling returns correct values."""
        scaling = cantilever.doping_optimization_scaling()

        assert len(scaling) == 4
        assert scaling[0] == pytest.approx(1e-2)  # time
        assert scaling[1] == pytest.approx(1e-2)  # temp
        assert scaling[2] == pytest.approx(1e0)  # energy
        assert scaling[3] == pytest.approx(1e-13)  # dose

    def test_scaling_is_array(self, cantilever):
        """Test scaling returns numpy array."""
        scaling = cantilever.doping_optimization_scaling()
        assert isinstance(scaling, np.ndarray)


class TestOptimizationBounds:
    """Test doping_optimization_bounds method."""

    @pytest.fixture
    def cantilever(self):
        """Create a test cantilever."""
        return CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=3600,
            annealing_temp=1173.15,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

    def test_default_bounds(self, cantilever):
        """Test default optimization bounds."""
        lb, ub = cantilever.doping_optimization_bounds()

        assert len(lb) == 4
        assert len(ub) == 4

        # Check time bounds (seconds)
        assert lb[0] == pytest.approx(15 * 60)
        assert ub[0] == pytest.approx(120 * 60)

        # Check temp bounds (K)
        assert lb[1] == pytest.approx(273.15 + 900)
        assert ub[1] == pytest.approx(273.15 + 1100)

        # Check energy bounds (keV)
        assert lb[2] == pytest.approx(20)
        assert ub[2] == pytest.approx(80)

        # Check dose bounds (cm^-2)
        assert lb[3] == pytest.approx(2e14)
        assert ub[3] == pytest.approx(2e16)

    def test_bounds_are_arrays(self, cantilever):
        """Test bounds return numpy arrays."""
        lb, ub = cantilever.doping_optimization_bounds()
        assert isinstance(lb, np.ndarray)
        assert isinstance(ub, np.ndarray)

    def test_lower_bounds_less_than_upper(self, cantilever):
        """Test all lower bounds are less than upper bounds."""
        lb, ub = cantilever.doping_optimization_bounds()
        assert np.all(lb < ub)

    def test_custom_bounds(self, cantilever):
        """Test custom constraint override."""
        constraints = {
            "min_annealing_time": 1800,
            "max_annealing_time": 5400,
        }

        lb, ub = cantilever.doping_optimization_bounds(constraints)

        assert lb[0] == pytest.approx(1800)
        assert ub[0] == pytest.approx(5400)
        # Others should remain default
        assert lb[1] == pytest.approx(273.15 + 900)
        assert ub[1] == pytest.approx(273.15 + 1100)

    def test_multiple_custom_bounds(self, cantilever):
        """Test multiple custom constraints."""
        constraints = {
            "min_annealing_temp": 273.15 + 950,
            "max_annealing_temp": 273.15 + 1050,
            "min_implantation_energy": 30,
            "max_implantation_energy": 70,
        }

        lb, ub = cantilever.doping_optimization_bounds(constraints)

        assert lb[1] == pytest.approx(273.15 + 950)
        assert ub[1] == pytest.approx(273.15 + 1050)
        assert lb[2] == pytest.approx(30)
        assert ub[2] == pytest.approx(70)


class TestDopingState:
    """Test doping state methods."""

    @pytest.fixture
    def cantilever(self):
        """Create a test cantilever."""
        return CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=3600,
            annealing_temp=1173.15,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

    def test_current_state(self, cantilever):
        """Test doping_current_state returns correct values."""
        state = cantilever.doping_current_state()

        assert len(state) == 4
        assert state[0] == pytest.approx(3600)
        assert state[1] == pytest.approx(1173.15)
        assert state[2] == pytest.approx(50)
        assert state[3] == pytest.approx(1e15)

    def test_current_state_is_array(self, cantilever):
        """Test current state returns numpy array."""
        state = cantilever.doping_current_state()
        assert isinstance(state, np.ndarray)

    def test_from_state(self, cantilever):
        """Test doping_cantilever_from_state updates parameters."""
        # Create state vector with 10 elements (6 base + 4 doping)
        new_state = np.array([0, 0, 0, 0, 0, 0, 7200, 1273.15, 60, 5e15])

        result = cantilever.doping_cantilever_from_state(new_state)

        assert result is cantilever  # Should return self
        assert cantilever.annealing_time == pytest.approx(7200)
        assert cantilever.annealing_temp == pytest.approx(1273.15)
        assert cantilever.implantation_energy == pytest.approx(60)
        assert cantilever.implantation_dose == pytest.approx(5e15)

    def test_state_round_trip(self, cantilever):
        """Test get/set state round trip preserves values."""
        original_state = cantilever.doping_current_state()

        # Modify the cantilever
        cantilever.annealing_time = 7200
        cantilever.annealing_temp = 1273.15

        # Restore original state
        full_state = np.concatenate([np.zeros(6), original_state])
        cantilever.doping_cantilever_from_state(full_state)

        restored_state = cantilever.doping_current_state()
        np.testing.assert_array_almost_equal(original_state, restored_state)


class TestRandomInitialConditions:
    """Test doping_initial_conditions_random method."""

    @pytest.fixture
    def cantilever(self):
        """Create a test cantilever."""
        return CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=3600,
            annealing_temp=1173.15,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

    def test_random_within_bounds(self, cantilever):
        """Test random initial conditions are within bounds."""
        lb, ub = cantilever.doping_optimization_bounds()

        for _ in range(10):  # Test multiple times due to randomness
            x0 = cantilever.doping_initial_conditions_random()
            assert np.all(x0 >= lb)
            assert np.all(x0 <= ub)

    def test_random_is_array(self, cantilever):
        """Test random initial conditions returns numpy array."""
        x0 = cantilever.doping_initial_conditions_random()
        assert isinstance(x0, np.ndarray)
        assert len(x0) == 4

    def test_random_with_custom_constraints(self, cantilever):
        """Test random initial conditions with custom bounds."""
        constraints = {
            "min_annealing_time": 1800,
            "max_annealing_time": 3600,
        }

        lb, ub = cantilever.doping_optimization_bounds(constraints)

        for _ in range(10):
            x0 = cantilever.doping_initial_conditions_random(constraints)
            assert x0[0] >= lb[0]
            assert x0[0] <= ub[0]

    def test_random_produces_different_values(self, cantilever):
        """Test random produces different values on repeated calls."""
        x1 = cantilever.doping_initial_conditions_random()
        x2 = cantilever.doping_initial_conditions_random()

        # Very unlikely to get exactly the same random values
        assert not np.allclose(x1, x2)


class TestLookupTableInterpolation:
    """Test methods that use lookup table interpolation with bundled data."""

    @pytest.fixture
    def cantilever(self):
        """Create cantilever (lookup table auto-loaded)."""
        return CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=3600,  # 60 min
            annealing_temp=1173.15,  # 900 C
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

    def test_lookup_data_loaded(self, cantilever):
        """Test lookup data is automatically loaded."""
        assert cantilever._lookup_data is not None
        assert "z" in cantilever._lookup_data
        assert "ImplantDopants" in cantilever._lookup_data
        assert "n" in cantilever._lookup_data

    def test_junction_depth(self, cantilever):
        """Test junction_depth property."""
        Xj = cantilever.junction_depth
        assert isinstance(Xj, (float, np.floating))
        assert Xj > 0
        assert Xj < 10e-6  # Reasonable range for junction depth

    def test_sheet_resistance(self, cantilever):
        """Test sheet_resistance method."""
        Rs = cantilever.sheet_resistance()
        assert isinstance(Rs, (float, np.floating))
        assert Rs > 0

    def test_Nz(self, cantilever):
        """Test Nz method."""
        Nz = cantilever.Nz()
        assert isinstance(Nz, (float, np.floating))
        assert Nz > 0

    def test_beta(self, cantilever):
        """Test beta calculation."""
        beta = cantilever.beta()
        assert isinstance(beta, (float, np.floating))

    def test_doping_profile(self, cantilever):
        """Test doping_profile method."""
        x, active, total = cantilever.doping_profile()

        # Check return types
        assert isinstance(x, np.ndarray)
        assert isinstance(active, np.ndarray)
        assert isinstance(total, np.ndarray)

        # Check lengths match
        assert len(x) == len(active)
        assert len(x) == len(total)

        # Check values are reasonable
        assert np.all(x >= 0)
        assert np.all(x <= cantilever.t)
        assert np.all(active > 0)
        assert np.all(total > 0)

    def test_doping_profile_truncated_at_thickness(self, cantilever):
        """Test doping profile is truncated at device thickness."""
        x, active, total = cantilever.doping_profile()

        # All x values should be <= thickness
        assert np.all(x <= cantilever.t)


class TestDopingProcessMetrics:
    """Test public process metrics for implanted cantilevers."""

    @pytest.fixture
    def cantilever(self):
        """Create an implanted cantilever for metrics tests."""
        return CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=3600,
            annealing_temp=1173.15,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

    def test_metrics_are_finite(self, cantilever):
        """Test process metrics return finite values."""
        metrics = cantilever.doping_process_metrics()
        values = [
            metrics.alpha_h,
            metrics.nz,
            metrics.beta1,
            metrics.beta2_um,
            metrics.beta,
            metrics.sheet_resistance,
            metrics.junction_depth_m,
            metrics.peak_concentration_cm3,
            metrics.retained_dose_cm2,
        ]

        assert np.all(np.isfinite(values))

    def test_beta_uses_default_thickness(self, cantilever):
        """Test beta combines beta1 and beta2 using cantilever thickness."""
        metrics = cantilever.doping_process_metrics()
        expected_beta = metrics.beta1 - 2 * metrics.beta2_um / (cantilever.t * 1e6)

        assert metrics.beta == pytest.approx(expected_beta)

    def test_thickness_override_changes_only_combined_beta(self, cantilever):
        """Test thickness override affects beta but not lookup beta values."""
        default_metrics = cantilever.doping_process_metrics()
        override_metrics = cantilever.doping_process_metrics(device_thickness_m=4e-6)

        assert override_metrics.beta != pytest.approx(default_metrics.beta)
        assert override_metrics.beta1 == pytest.approx(default_metrics.beta1)
        assert override_metrics.beta2_um == pytest.approx(default_metrics.beta2_um)


class TestDopingOptimization:
    """Test public process-only doping optimization."""

    @pytest.fixture
    def cantilever(self):
        """Create a DopeDealer phosphorus cantilever for optimization tests."""
        return CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="phosphorus",
            annealing_time=3600,
            annealing_temp=273.15 + 950,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
            lookup_source="dopedealer",
        )

    @staticmethod
    def default_objective(metrics):
        """Evaluate the documented default objective."""
        return np.log(metrics.alpha_h) - np.log(metrics.nz) - 2 * np.log(abs(metrics.beta))

    def test_optimize_doping_runs_and_improves_default_objective(self, cantilever):
        """Test default optimizer returns an improved physical process state."""
        original_state = cantilever.doping_current_state().copy()
        starting_objective = self.default_objective(cantilever.doping_process_metrics())

        result = cantilever.optimize_doping_for_hooge_noise(log_dose=True)
        lb, ub = cantilever.doping_optimization_bounds()

        assert result.objective_value <= starting_objective + 1e-8
        assert np.all(result.state >= lb)
        assert np.all(result.state <= ub)
        assert result.optimized is not cantilever
        np.testing.assert_allclose(cantilever.doping_current_state(), original_state)
        assert len(result.all_results) == 1

    def test_boundary_starts_do_not_raise_out_of_bounds(self, cantilever):
        """Test lower and upper boundary starts stay valid for interpolation."""
        lb, ub = cantilever.doping_optimization_bounds()

        result = cantilever.optimize_doping_for_hooge_noise(
            initial_states=[lb, ub],
            log_dose=True,
        )

        assert np.isfinite(result.objective_value)
        assert len(result.all_results) == 2
        assert set(result.boundary_flags) == {
            "annealing_time_at_min",
            "annealing_time_at_max",
            "annealing_temp_at_min",
            "annealing_temp_at_max",
            "implantation_energy_at_min",
            "implantation_energy_at_max",
            "implantation_dose_at_min",
            "implantation_dose_at_max",
        }


class TestMetricConstrainedDopingOptimization:
    """Test metric_constraints in optimize_doping_for_hooge_noise."""

    @pytest.fixture
    def cantilever(self):
        """Create a DopeDealer phosphorus cantilever for constrained optimization."""
        return CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="phosphorus",
            annealing_time=3600,
            annealing_temp=273.15 + 950,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
            lookup_source="dopedealer",
        )

    @pytest.mark.parametrize(
        "constraint",
        [
            MetricConstraint(metric=DopingMetric.SHEET_RESISTANCE),
            MetricConstraint(metric=DopingMetric.SHEET_RESISTANCE, minimum=float("inf")),
            MetricConstraint(metric=DopingMetric.SHEET_RESISTANCE, maximum=float("nan")),
            MetricConstraint(metric=DopingMetric.SHEET_RESISTANCE, minimum=1000.0, maximum=100.0),
        ],
    )
    def test_validation_rejects_bad_bounds(self, cantilever, constraint):
        """Test malformed MetricConstraints raise ValueError."""
        with pytest.raises(ValueError):
            cantilever.optimize_doping_for_hooge_noise(metric_constraints=[constraint])

    def test_validation_rejects_raw_string_metric(self, cantilever):
        """Test passing a raw string instead of DopingMetric raises ValueError."""
        bad = MetricConstraint(metric="sheet_resistance", minimum=300.0)
        with pytest.raises(ValueError):
            cantilever.optimize_doping_for_hooge_noise(metric_constraints=[bad])

    def test_lbfgsb_with_metric_constraints_rejected(self, cantilever):
        """Test L-BFGS-B with metric_constraints raises ValueError."""
        constraints = [MetricConstraint(metric=DopingMetric.SHEET_RESISTANCE, minimum=300.0, maximum=900.0)]
        with pytest.raises(ValueError, match="L-BFGS-B"):
            cantilever.optimize_doping_for_hooge_noise(
                metric_constraints=constraints,
                method="L-BFGS-B",
            )

    def test_sheet_resistance_band_enforced(self, cantilever):
        """Test SLSQP enforces sheet_resistance band on a feasible case."""
        rs_min = 300.0
        rs_max = 900.0
        constraints = [MetricConstraint(metric=DopingMetric.SHEET_RESISTANCE, minimum=rs_min, maximum=rs_max)]

        result = cantilever.optimize_doping_for_hooge_noise(
            metric_constraints=constraints,
            n_random_starts=3,
            random_seed=0,
            log_dose=True,
        )

        tol = 1e-6 * max(rs_max, rs_min)
        assert rs_min - tol <= result.metrics.sheet_resistance <= rs_max + tol
        assert result.scipy_result.success

    def test_infeasible_metric_constraints_raise_runtimeerror(self, cantilever):
        """Test impossible metric constraints raise RuntimeError."""
        impossible = [MetricConstraint(metric=DopingMetric.SHEET_RESISTANCE, minimum=1e8, maximum=1e9)]
        with pytest.raises(RuntimeError):
            cantilever.optimize_doping_for_hooge_noise(
                metric_constraints=impossible,
                n_random_starts=2,
                random_seed=0,
            )

    def test_original_cantilever_not_mutated(self, cantilever):
        """Test the input cantilever is unchanged after a constrained run."""
        snapshot = cantilever.doping_current_state().copy()
        constraints = [MetricConstraint(metric=DopingMetric.SHEET_RESISTANCE, minimum=300.0, maximum=900.0)]

        result = cantilever.optimize_doping_for_hooge_noise(
            metric_constraints=constraints,
            random_seed=0,
        )

        np.testing.assert_allclose(cantilever.doping_current_state(), snapshot)
        assert result.optimized is not cantilever

    def test_device_thickness_threads_into_beta_constraint(self, cantilever):
        """Test device_thickness_m is used by beta-based metric constraints.

        Pin process variables to a narrow window around the current state so the
        optimizer cannot move significantly. Then beta is determined almost
        entirely by device_thickness_m via beta = beta1 - 2*beta2_um/thickness_um.
        A constraint feasible at one thickness must be infeasible at another.
        """
        state = cantilever.doping_current_state()
        narrow = {
            "min_annealing_time": float(state[0]) * 0.999,
            "max_annealing_time": float(state[0]) * 1.001,
            "min_annealing_temp": float(state[1]) - 0.1,
            "max_annealing_temp": float(state[1]) + 0.1,
            "min_implantation_energy": float(state[2]) * 0.999,
            "max_implantation_energy": float(state[2]) * 1.001,
            "min_implantation_dose": float(state[3]) * 0.999,
            "max_implantation_dose": float(state[3]) * 1.001,
        }

        thin_metrics = cantilever.doping_process_metrics(device_thickness_m=2e-6)
        thick_metrics = cantilever.doping_process_metrics(device_thickness_m=10e-6)
        assert thin_metrics.beta != pytest.approx(thick_metrics.beta)

        midpoint = thick_metrics.beta
        half_band = abs(thick_metrics.beta - thin_metrics.beta) * 0.1
        beta_min = midpoint - half_band
        beta_max = midpoint + half_band

        assert not (beta_min <= thin_metrics.beta <= beta_max)
        assert beta_min <= thick_metrics.beta <= beta_max

        constraints = [MetricConstraint(metric=DopingMetric.BETA, minimum=beta_min, maximum=beta_max)]

        result_thick = cantilever.optimize_doping_for_hooge_noise(
            metric_constraints=constraints,
            parameter_constraints=narrow,
            device_thickness_m=10e-6,
            initial_states=[state],
        )
        assert beta_min - 1e-6 <= result_thick.metrics.beta <= beta_max + 1e-6

        with pytest.raises(RuntimeError):
            cantilever.optimize_doping_for_hooge_noise(
                metric_constraints=constraints,
                parameter_constraints=narrow,
                device_thickness_m=2e-6,
                initial_states=[state],
            )

    def test_multistart_prefers_feasible_over_lower_infeasible(self, cantilever, monkeypatch):
        """Test selection picks feasible result over lower-objective infeasible one.

        Monkeypatches scipy.optimize.minimize so two starts return crafted
        OptimizeResults: a lower-objective infeasible one and a higher-objective
        feasible one. Selection must return the feasible one even though its
        objective is worse.
        """
        from scipy.optimize import OptimizeResult

        import piezod.cantilever_implantation as impl_mod

        state_low_rs = np.array([3600.0, 273.15 + 1050.0, 100.0, 5e16])
        state_high_rs = np.array([3600.0, 273.15 + 900.0, 30.0, 5e13])

        metrics_low_rs = cantilever._copy_with_doping_state(state_low_rs).doping_process_metrics()
        metrics_high_rs = cantilever._copy_with_doping_state(state_high_rs).doping_process_metrics()

        rs_min = (metrics_low_rs.sheet_resistance + metrics_high_rs.sheet_resistance) / 2
        rs_max = metrics_high_rs.sheet_resistance * 2.0

        assert metrics_low_rs.sheet_resistance < rs_min
        assert rs_min <= metrics_high_rs.sheet_resistance <= rs_max

        log_low = state_low_rs.copy()
        log_low[3] = np.log10(log_low[3])
        log_high = state_high_rs.copy()
        log_high[3] = np.log10(log_high[3])

        call_counter = {"n": 0}

        def fake_minimize(fun, x0, **kwargs):
            call_counter["n"] += 1
            if call_counter["n"] == 1:
                return OptimizeResult(x=log_low, fun=-100.0, success=True, status=0, message="ok", nit=1)
            return OptimizeResult(x=log_high, fun=-50.0, success=True, status=0, message="ok", nit=1)

        monkeypatch.setattr(impl_mod, "minimize", fake_minimize)

        constraints = [MetricConstraint(metric=DopingMetric.SHEET_RESISTANCE, minimum=rs_min, maximum=rs_max)]
        result = cantilever.optimize_doping_for_hooge_noise(
            metric_constraints=constraints,
            initial_states=[state_low_rs, state_high_rs],
            log_dose=True,
        )

        assert result.objective_value == pytest.approx(-50.0)
        assert rs_min <= result.metrics.sheet_resistance <= rs_max
        assert call_counter["n"] == 2


class TestInheritance:
    """Test that CantileverImplantation properly inherits from Cantilever."""

    def test_inherits_from_cantilever(self):
        """Test CantileverImplantation is a subclass of Cantilever."""
        from piezod.cantilever import Cantilever

        assert issubclass(CantileverImplantation, Cantilever)

    def test_has_cantilever_attributes(self):
        """Test cantilever has inherited attributes."""
        cantilever = CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=3600,
            annealing_temp=1173.15,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

        # Check inherited attributes
        assert hasattr(cantilever, "l")
        assert hasattr(cantilever, "w")
        assert hasattr(cantilever, "t")
        assert hasattr(cantilever, "doping_type")
        assert hasattr(cantilever, "k_b_eV")

    def test_has_dopantNumber_method(self):
        """Test cantilever has inherited dopantNumber method."""
        cantilever = CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=3600,
            annealing_temp=1173.15,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

        assert hasattr(cantilever, "dopantNumber")
        assert cantilever.dopantNumber() == 1  # boron = 1
