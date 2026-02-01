"""Tests for CantileverImplantation class."""

import numpy as np
import pytest

from piezod.cantilever_implantation import CantileverImplantation


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
            "annealing_temp": 1173,  # 900C
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
        assert cantilever.annealing_temp == 1173
        assert cantilever.annealing_time == 3600
        assert cantilever.lookupTableData is None

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
            annealing_temp=1173,
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
            annealing_temp=1173,
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
            "annealing_temp": 1173,  # 900C
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
        T = 1173  # K
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
        T = 1173  # K
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
        T = 1173  # K
        t = 3600  # s
        D = D0 * np.exp(-Ea / k_b_eV / T)
        expected = np.sqrt(D * t)

        assert cantilever.diffusion_length == pytest.approx(expected, rel=1e-10)

    def test_diffusion_length_temperature_dependence(self, cantilever_params):
        """Test diffusion length increases with temperature."""
        cantilever_params["doping_type"] = "boron"

        cantilever_low = CantileverImplantation(**cantilever_params)
        cantilever_params["annealing_temp"] = 1373  # Higher temp
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

        with pytest.raises(ValueError, match="Unknown doping type"):
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
            annealing_temp=1173,
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
            annealing_temp=1173,
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
            annealing_temp=1173,
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
            annealing_temp=1173,
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
            annealing_temp=1173,
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
        assert lb[1] == pytest.approx(273 + 900)
        assert ub[1] == pytest.approx(273 + 1100)

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
        assert lb[1] == pytest.approx(273 + 900)
        assert ub[1] == pytest.approx(273 + 1100)

    def test_multiple_custom_bounds(self, cantilever):
        """Test multiple custom constraints."""
        constraints = {
            "min_annealing_temp": 273 + 950,
            "max_annealing_temp": 273 + 1050,
            "min_implantation_energy": 30,
            "max_implantation_energy": 70,
        }

        lb, ub = cantilever.doping_optimization_bounds(constraints)

        assert lb[1] == pytest.approx(273 + 950)
        assert ub[1] == pytest.approx(273 + 1050)
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
            annealing_temp=1173,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

    def test_current_state(self, cantilever):
        """Test doping_current_state returns correct values."""
        state = cantilever.doping_current_state()

        assert len(state) == 4
        assert state[0] == pytest.approx(3600)
        assert state[1] == pytest.approx(1173)
        assert state[2] == pytest.approx(50)
        assert state[3] == pytest.approx(1e15)

    def test_current_state_is_array(self, cantilever):
        """Test current state returns numpy array."""
        state = cantilever.doping_current_state()
        assert isinstance(state, np.ndarray)

    def test_from_state(self, cantilever):
        """Test doping_cantilever_from_state updates parameters."""
        # Create state vector with 10 elements (6 base + 4 doping)
        new_state = np.array([0, 0, 0, 0, 0, 0, 7200, 1273, 60, 5e15])

        result = cantilever.doping_cantilever_from_state(new_state)

        assert result is cantilever  # Should return self
        assert cantilever.annealing_time == pytest.approx(7200)
        assert cantilever.annealing_temp == pytest.approx(1273)
        assert cantilever.implantation_energy == pytest.approx(60)
        assert cantilever.implantation_dose == pytest.approx(5e15)

    def test_state_round_trip(self, cantilever):
        """Test get/set state round trip preserves values."""
        original_state = cantilever.doping_current_state()

        # Modify the cantilever
        cantilever.annealing_time = 7200
        cantilever.annealing_temp = 1273

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
            annealing_temp=1173,
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


class TestLookupTableMethods:
    """Test methods that require lookup table."""

    @pytest.fixture
    def cantilever(self):
        """Create a test cantilever without lookup table."""
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
            annealing_temp=1173,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

    def test_junction_depth_without_table(self, cantilever):
        """Test junction_depth raises error without lookup table."""
        with pytest.raises(RuntimeError, match="Lookup table not loaded"):
            _ = cantilever.junction_depth

    def test_sheet_resistance_without_table(self, cantilever):
        """Test sheet_resistance raises error without lookup table."""
        with pytest.raises(RuntimeError, match="Lookup table not loaded"):
            cantilever.sheet_resistance()

    def test_Nz_without_table(self, cantilever):
        """Test Nz raises error without lookup table."""
        with pytest.raises(RuntimeError, match="Lookup table not loaded"):
            cantilever.Nz()

    def test_beta_without_table(self, cantilever):
        """Test beta raises error without lookup table."""
        with pytest.raises(RuntimeError, match="Lookup table not loaded"):
            cantilever.beta()

    def test_doping_profile_without_table(self, cantilever):
        """Test doping_profile raises error without lookup table."""
        with pytest.raises(RuntimeError, match="Lookup table not loaded"):
            cantilever.doping_profile()

    def test_load_lookup_table(self, cantilever):
        """Test loading lookup table data."""
        # Create minimal mock lookup table
        mock_data = {
            "z": np.linspace(0, 5, 500),  # microns
            "ImplantDopants": np.array([1, 2, 3]),
            "ImplantDoses": np.array([2e14, 1e15, 2e16]),
            "ImplantEnergies": np.array([20, 50, 80]),
            "AnnealTemps": np.array([900, 1000, 1100]),
            "AnnealTimes": np.array([15, 60, 120]),
            "AnnealOxidation": np.array([1, 2]),
        }

        cantilever.load_lookup_table(mock_data)
        assert cantilever.lookupTableData is not None
        assert cantilever.lookupTableData == mock_data


class TestLookupTableInterpolation:
    """Test methods that use lookup table interpolation with mock data."""

    @pytest.fixture
    def cantilever_with_table(self):
        """Create cantilever with mock lookup table."""
        cantilever = CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=100e-6,
            w=20e-6,
            t=2e-6,
            l_pr_ratio=0.5,
            v_bridge=5.0,
            doping_type="boron",
            annealing_time=3600,  # 60 min
            annealing_temp=1173,  # 900C
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

        # Create comprehensive mock lookup table
        # Grid dimensions: [dopant, dose, energy, temp, time, oxidation]
        dopants = np.array([1, 2, 3])  # B, P, As
        doses = np.array([2e14, 1e15, 2e16])
        energies = np.array([20, 50, 80])
        temps = np.array([900, 1000, 1100])
        times = np.array([15, 60, 120])
        oxidations = np.array([1, 2])

        # Create meshgrid for 6D interpolation
        shape = (len(dopants), len(doses), len(energies), len(temps), len(times), len(oxidations))

        # Create mock data with reasonable values
        mock_data = {
            "z": np.linspace(0, 5, 500),  # depth in microns
            "ImplantDopants": dopants,
            "ImplantDoses": doses,
            "ImplantEnergies": energies,
            "AnnealTemps": temps,
            "AnnealTimes": times,
            "AnnealOxidation": oxidations,
            "Xj": np.random.uniform(0.1e-6, 1e-6, shape),  # junction depth in m
            "Rs": np.random.uniform(10, 1000, shape),  # sheet resistance
            "Nz": np.random.uniform(1e15, 1e20, shape),  # effective doping
            "Beta1": np.random.uniform(1e-10, 1e-11, shape),  # piezoresistive coeff
            "Beta2": np.random.uniform(1e-15, 1e-16, shape),  # piezoresistive coeff
        }

        # Create 7D array for doping profile (z, dopant, dose, energy, temp, time, oxidation)
        shape_7d = (len(mock_data["z"]),) + shape
        mock_data["n"] = np.random.uniform(1e15, 1e20, shape_7d)

        cantilever.load_lookup_table(mock_data)
        return cantilever

    def test_junction_depth_with_table(self, cantilever_with_table):
        """Test junction_depth property with lookup table."""
        Xj = cantilever_with_table.junction_depth
        assert isinstance(Xj, (float, np.floating))
        assert Xj > 0
        assert Xj < 10e-6  # Reasonable range for junction depth

    def test_sheet_resistance_with_table(self, cantilever_with_table):
        """Test sheet_resistance with lookup table."""
        Rs = cantilever_with_table.sheet_resistance()
        assert isinstance(Rs, (float, np.floating))
        assert Rs > 0

    def test_Nz_with_table(self, cantilever_with_table):
        """Test Nz with lookup table."""
        Nz = cantilever_with_table.Nz()
        assert isinstance(Nz, (float, np.floating))
        assert Nz > 0

    def test_beta_with_table(self, cantilever_with_table):
        """Test beta calculation with lookup table."""
        beta = cantilever_with_table.beta()
        assert isinstance(beta, (float, np.floating))

    def test_doping_profile_with_table(self, cantilever_with_table):
        """Test doping_profile with lookup table."""
        x, active, total = cantilever_with_table.doping_profile()

        # Check return types
        assert isinstance(x, np.ndarray)
        assert isinstance(active, np.ndarray)
        assert isinstance(total, np.ndarray)

        # Check lengths match
        assert len(x) == len(active)
        assert len(x) == len(total)

        # Check values are reasonable
        assert np.all(x >= 0)
        assert np.all(x <= cantilever_with_table.t)
        assert np.all(active > 0)
        assert np.all(total > 0)

    def test_doping_profile_truncated_at_thickness(self, cantilever_with_table):
        """Test doping profile is truncated at device thickness."""
        x, active, total = cantilever_with_table.doping_profile()

        # All x values should be <= thickness
        assert np.all(x <= cantilever_with_table.t)


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
            annealing_temp=1173,
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
            annealing_temp=1173,
            annealing_type="inert",
            implantation_energy=50,
            implantation_dose=1e15,
        )

        assert hasattr(cantilever, "dopantNumber")
        assert cantilever.dopantNumber() == 1  # boron = 1
