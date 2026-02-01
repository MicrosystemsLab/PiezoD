"""Tests for CantileverDiffusion class."""

import numpy as np
import pytest

from piezod import CantileverDiffusion


class TestCantileverDiffusionInitialization:
    """Test CantileverDiffusion initialization."""

    def test_init_phosphorus(self) -> None:
        """Test initialization with phosphorus doping."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        assert cant.diffusion_time == 30 * 60
        assert cant.diffusion_temp == 273 + 900
        assert cant.doping_type == "phosphorus"
        assert cant.l == 100e-6
        assert cant.w == 10e-6
        assert cant.t == 1e-6

    def test_init_arsenic(self) -> None:
        """Test initialization with arsenic doping."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="arsenic",
            diffusion_time=45 * 60,
            diffusion_temp=273 + 950,
        )
        assert cant.doping_type == "arsenic"
        assert cant.diffusion_time == 45 * 60

    def test_init_boron(self) -> None:
        """Test initialization with boron doping."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="boron",
            diffusion_time=60 * 60,
            diffusion_temp=273 + 1000,
        )
        assert cant.doping_type == "boron"


class TestDopingProfile:
    """Test doping_profile method for each dopant type."""

    def test_arsenic_profile(self) -> None:
        """Test arsenic diffusion profile calculation."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="arsenic",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 950,
        )
        x, active_doping, total_doping = cant.doping_profile()

        # Check array shapes
        assert len(x) == cant.numZPoints
        assert len(active_doping) == cant.numZPoints
        assert len(total_doping) == cant.numZPoints

        # Check that position array goes from 0 to cantilever thickness
        assert x[0] == pytest.approx(0.0)
        assert x[-1] == pytest.approx(cant.t)

        # Check that doping decreases with depth
        assert active_doping[0] > active_doping[-1]
        assert total_doping[0] > total_doping[-1]

        # Check that active and total are equal for arsenic
        np.testing.assert_array_almost_equal(active_doping, total_doping)

        # Surface concentration should be high
        assert active_doping[0] > 1e19

    def test_boron_profile(self) -> None:
        """Test boron diffusion profile calculation."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="boron",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 950,
        )
        x, active_doping, total_doping = cant.doping_profile()

        # Check array shapes
        assert len(x) == cant.numZPoints
        assert len(active_doping) == cant.numZPoints
        assert len(total_doping) == cant.numZPoints

        # Check that position array is valid
        assert x[0] == pytest.approx(0.0)
        assert x[-1] == pytest.approx(cant.t)

        # Check monotonic decrease
        assert active_doping[0] > active_doping[-1]

        # Check that active and total are equal for boron
        np.testing.assert_array_almost_equal(active_doping, total_doping)

    def test_phosphorus_profile(self) -> None:
        """Test phosphorus diffusion profile calculation."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        x, active_doping, total_doping = cant.doping_profile()

        # Check array shapes
        assert len(x) == cant.numZPoints
        assert len(active_doping) == cant.numZPoints
        assert len(total_doping) == cant.numZPoints

        # Check that position array is valid
        assert x[0] == pytest.approx(0.0)
        assert x[-1] == pytest.approx(cant.t)

        # Check that doping decreases with depth
        assert total_doping[0] > total_doping[-1]

        # Active should be less than or equal to total
        assert np.all(active_doping <= total_doping)

        # Check background doping floor
        assert np.all(total_doping >= 1e15)

        # Surface concentration should be very high for phosphorus
        assert total_doping[0] > 1e20

    def test_phosphorus_low_temp(self) -> None:
        """Test phosphorus diffusion at low temperature (no time offset)."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 750,  # Below 780C threshold
        )
        x, active_doping, total_doping = cant.doping_profile()

        # Should still produce valid profile
        assert len(x) == cant.numZPoints
        assert np.all(total_doping >= 1e15)

    def test_invalid_doping_type(self) -> None:
        """Test that invalid doping type raises ValueError."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="invalid",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        with pytest.raises(ValueError, match="Unknown doping type"):
            cant.doping_profile()


class TestSheetResistance:
    """Test sheet_resistance calculation."""

    def test_sheet_resistance_phosphorus(self) -> None:
        """Test sheet resistance for phosphorus."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        Rs = cant.sheet_resistance()

        # Should be positive and in reasonable range (1-1000 ohms/sq)
        assert Rs > 0
        assert Rs < 1e4

    def test_sheet_resistance_arsenic(self) -> None:
        """Test sheet resistance for arsenic."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="arsenic",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 950,
        )
        Rs = cant.sheet_resistance()
        assert Rs > 0


class TestNz:
    """Test effective carrier density calculation."""

    def test_Nz_phosphorus(self) -> None:
        """Test Nz for phosphorus."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        Nz = cant.Nz()

        # Should be positive and in reasonable range
        assert Nz > 0
        # Should be in units of 1/m^2
        assert Nz > 1e18  # At least 1e14/cm^2 = 1e18/m^2

    def test_Nz_arsenic(self) -> None:
        """Test Nz for arsenic."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="arsenic",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 950,
        )
        Nz = cant.Nz()
        assert Nz > 0

    def test_Nz_boron(self) -> None:
        """Test Nz for boron."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="boron",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 1000,
        )
        Nz = cant.Nz()
        assert Nz > 0


class TestJunctionDepth:
    """Test junction_depth property."""

    def test_junction_depth_phosphorus(self) -> None:
        """Test junction depth for phosphorus."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        xj = cant.junction_depth

        # Should be positive and less than cantilever thickness
        assert xj > 0
        assert xj <= cant.t

    def test_junction_depth_property(self) -> None:
        """Test that junction_depth is a property."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        # Access as property (no parentheses)
        xj1 = cant.junction_depth
        xj2 = cant.junction_depth
        assert xj1 == xj2


class TestAlpha:
    """Test alpha method."""

    def test_alpha_returns_default(self) -> None:
        """Test that alpha returns default_alpha."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        assert cant.alpha() == cant.default_alpha
        assert cant.alpha() == pytest.approx(1e-5)


class TestOptimizationMethods:
    """Test optimization-related methods."""

    def test_doping_optimization_scaling(self) -> None:
        """Test optimization scaling factors."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        scaling = cant.doping_optimization_scaling()
        assert len(scaling) == 2
        assert scaling[0] == 1e-3
        assert scaling[1] == 1e-3

    def test_doping_current_state(self) -> None:
        """Test getting current doping state."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        state = cant.doping_current_state()
        assert len(state) == 2
        assert state[0] == 30 * 60
        assert state[1] == 273 + 900

    def test_doping_cantilever_from_state(self) -> None:
        """Test setting doping parameters from state vector."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        # State vector: [l, w, t, l_pr_ratio, v_bridge, diffusion_time, diffusion_temp]
        new_state = [100e-6, 10e-6, 1e-6, 0.3, 1.0, 45 * 60, 273 + 950]
        cant.doping_cantilever_from_state(new_state)

        assert cant.diffusion_time == 45 * 60
        assert cant.diffusion_temp == 273 + 950

    def test_doping_optimization_bounds_default(self) -> None:
        """Test optimization bounds with default constraints."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        lb, ub = cant.doping_optimization_bounds()

        assert len(lb) == 2
        assert len(ub) == 2

        # Check default bounds
        assert lb[0] == 5 * 60  # min diffusion time
        assert ub[0] == 90 * 60  # max diffusion time
        assert lb[1] == 273 + 800  # min diffusion temp
        assert ub[1] == 273 + 1000  # max diffusion temp

    def test_doping_optimization_bounds_with_constraints(self) -> None:
        """Test optimization bounds with custom constraints."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        constraints = {
            "min_diffusion_time": 10 * 60,
            "max_diffusion_time": 60 * 60,
            "min_diffusion_temp": 273 + 850,
            "max_diffusion_temp": 273 + 950,
        }
        lb, ub = cant.doping_optimization_bounds(constraints)

        assert lb[0] == 10 * 60
        assert ub[0] == 60 * 60
        assert lb[1] == 273 + 850
        assert ub[1] == 273 + 950

    def test_doping_optimization_bounds_partial_constraints(self) -> None:
        """Test optimization bounds with partial constraints."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        constraints = {"min_diffusion_time": 15 * 60}
        lb, ub = cant.doping_optimization_bounds(constraints)

        # Only min_diffusion_time should change
        assert lb[0] == 15 * 60
        assert ub[0] == 90 * 60  # default
        assert lb[1] == 273 + 800  # default
        assert ub[1] == 273 + 1000  # default

    def test_doping_initial_conditions_random_default(self) -> None:
        """Test random initial conditions with default bounds."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        x0 = cant.doping_initial_conditions_random()

        assert len(x0) == 2

        # Check that values are within bounds
        lb, ub = cant.doping_optimization_bounds()
        assert lb[0] <= x0[0] <= ub[0]
        assert lb[1] <= x0[1] <= ub[1]

    def test_doping_initial_conditions_random_with_constraints(self) -> None:
        """Test random initial conditions with custom constraints."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )
        constraints = {
            "min_diffusion_time": 20 * 60,
            "max_diffusion_time": 40 * 60,
            "min_diffusion_temp": 273 + 880,
            "max_diffusion_temp": 273 + 920,
        }
        x0 = cant.doping_initial_conditions_random(constraints)

        # Check that values are within custom bounds
        assert 20 * 60 <= x0[0] <= 40 * 60
        assert 273 + 880 <= x0[1] <= 273 + 920

    def test_doping_initial_conditions_random_generates_different_values(self) -> None:
        """Test that random initial conditions vary between calls."""
        cant = CantileverDiffusion(
            freq_min=1.0,
            freq_max=1e3,
            l=100e-6,
            w=10e-6,
            t=1e-6,
            l_pr_ratio=0.3,
            v_bridge=1.0,
            doping_type="phosphorus",
            diffusion_time=30 * 60,
            diffusion_temp=273 + 900,
        )

        # Generate multiple random initial conditions
        random_values = [cant.doping_initial_conditions_random() for _ in range(10)]

        # Check that not all values are identical (very unlikely with proper randomness)
        # Convert to tuples for comparison
        unique_values = {tuple(v) for v in random_values}
        assert len(unique_values) > 1
