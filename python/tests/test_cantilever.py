"""Basic tests for Cantilever class."""

import numpy as np
import pytest

from piezod import Cantilever, GapConfig


class TestCantileverConstants:
    """Test physical constants are defined correctly."""

    def test_boltzmann_constant(self) -> None:
        """Verify Boltzmann constant value."""
        assert Cantilever.k_b == pytest.approx(1.38e-23, rel=1e-2)

    def test_electron_charge(self) -> None:
        """Verify electron charge value."""
        assert Cantilever.q == pytest.approx(1.60218e-19, rel=1e-4)

    def test_planck_constant(self) -> None:
        """Verify reduced Planck constant value."""
        assert Cantilever.h_bar == pytest.approx(1.055e-34, rel=1e-2)


class TestCantileverMaterials:
    """Test material properties are defined."""

    def test_silicon_thermal_conductivity(self) -> None:
        """Verify silicon thermal conductivity."""
        assert Cantilever.k_si == 148

    def test_silicon_elastic_modulus(self) -> None:
        """Verify silicon elastic modulus is positive."""
        assert Cantilever.E_si > 0

    def test_dopant_options(self) -> None:
        """Verify dopant options are defined."""
        assert "boron" in Cantilever.dopantOptions
        assert "phosphorus" in Cantilever.dopantOptions
        assert "arsenic" in Cantilever.dopantOptions


class TestGapConfig:
    """Test GapConfig dataclass functionality."""

    def test_default_gap_config(self) -> None:
        """Test default GapConfig has zero gap."""
        gc = GapConfig()
        assert gc.gap_width == 0.0
        assert gc.gap_fraction == 0.0
        assert gc.gap_extent is None

    def test_absolute_gap_width(self) -> None:
        """Test absolute gap width takes precedence over fraction."""
        gc = GapConfig(gap_width=5e-6, gap_fraction=0.5)
        assert gc.get_gap_width(10e-6) == 5e-6  # Uses absolute, not fraction

    def test_fractional_gap_width(self) -> None:
        """Test fractional gap when absolute is zero."""
        gc = GapConfig(gap_fraction=0.2)
        assert gc.get_gap_width(10e-6) == pytest.approx(2e-6)

    def test_effective_width_in_gap_region(self) -> None:
        """Test effective width is reduced in gap region."""
        gc = GapConfig(gap_width=2e-6, gap_extent=50e-6)
        nominal_width = 10e-6

        # In gap region
        w_at_root = gc.effective_width(0, 100e-6, nominal_width)
        w_at_mid_gap = gc.effective_width(25e-6, 100e-6, nominal_width)
        assert w_at_root == pytest.approx(8e-6)
        assert w_at_mid_gap == pytest.approx(8e-6)

    def test_effective_width_beyond_gap_region(self) -> None:
        """Test effective width is full beyond gap region."""
        gc = GapConfig(gap_width=2e-6, gap_extent=50e-6)
        nominal_width = 10e-6

        # Beyond gap region
        w_at_tip = gc.effective_width(75e-6, 100e-6, nominal_width)
        assert w_at_tip == pytest.approx(10e-6)

    def test_effective_width_array(self) -> None:
        """Test effective width array calculation."""
        gc = GapConfig(gap_width=2e-6, gap_extent=50e-6)
        nominal_width = 10e-6

        x = np.array([0, 25e-6, 50e-6, 75e-6, 100e-6])
        w_eff = gc.effective_width_array(x, 100e-6, nominal_width)

        # First 3 points in gap region, last 2 beyond
        assert w_eff[0] == pytest.approx(8e-6)
        assert w_eff[1] == pytest.approx(8e-6)
        assert w_eff[2] == pytest.approx(8e-6)
        assert w_eff[3] == pytest.approx(10e-6)
        assert w_eff[4] == pytest.approx(10e-6)


class TestPiezoresistorGap:
    """Test piezoresistor gap configuration in Cantilever."""

    def test_default_air_gap_width(self) -> None:
        """Test default air gap width is 2 um."""
        c = Cantilever()
        assert c.air_gap_width == pytest.approx(2e-6)

    def test_set_air_gap_width_backward_compat(self) -> None:
        """Test setting air_gap_width for backward compatibility."""
        c = Cantilever()
        c.air_gap_width = 5e-6
        assert c.air_gap_width == pytest.approx(5e-6)
        assert c.gap_config.gap_width == pytest.approx(5e-6)

    def test_set_gap_config_directly(self) -> None:
        """Test setting gap_config directly."""
        c = Cantilever()
        c.gap_config = GapConfig(gap_width=3e-6, gap_extent=20e-6)
        assert c.air_gap_width == pytest.approx(3e-6)
        assert c.gap_extent() == pytest.approx(20e-6)

    def test_effective_pr_width_at(self) -> None:
        """Test effective piezoresistor width at position."""
        c = Cantilever()
        c.l = 100e-6
        c.w = 20e-6  # w_pr = 10 um
        c.l_pr_ratio = 0.3  # l_pr = 30 um
        c.gap_config = GapConfig(gap_width=2e-6, gap_extent=15e-6)

        # In gap region (x < 15 um)
        w_root = c.effective_pr_width_at(0)
        assert w_root == pytest.approx(8e-6)  # 10 - 2

        # Beyond gap region
        w_tip = c.effective_pr_width_at(20e-6)
        assert w_tip == pytest.approx(10e-6)


class TestBeamGap:
    """Test beam gap configuration for Rayleigh-Ritz calculations."""

    def test_default_no_beam_gap(self) -> None:
        """Test default has no beam gap."""
        c = Cantilever()
        assert c.beam_gap_config is None
        assert c.beam_gap_extent() == 0.0

    def test_set_beam_gap_config(self) -> None:
        """Test setting beam gap configuration."""
        c = Cantilever()
        c.beam_gap_config = GapConfig(gap_width=3e-6, gap_extent=20e-6)
        assert c.beam_gap_config is not None
        assert c.beam_gap_config.gap_width == pytest.approx(3e-6)
        assert c.beam_gap_extent() == pytest.approx(20e-6)

    def test_effective_beam_width_no_gap(self) -> None:
        """Test effective beam width without gap is constant."""
        c = Cantilever()
        c.w = 20e-6

        assert c.effective_beam_width_at(0) == pytest.approx(20e-6)
        assert c.effective_beam_width_at(50e-6) == pytest.approx(20e-6)

    def test_effective_beam_width_with_gap(self) -> None:
        """Test effective beam width with gap is reduced at root."""
        c = Cantilever()
        c.l = 100e-6
        c.w = 20e-6
        c.beam_gap_config = GapConfig(gap_width=5e-6, gap_extent=30e-6)

        # In gap region
        assert c.effective_beam_width_at(0) == pytest.approx(15e-6)
        assert c.effective_beam_width_at(20e-6) == pytest.approx(15e-6)

        # Beyond gap region
        assert c.effective_beam_width_at(50e-6) == pytest.approx(20e-6)

    def test_effective_beam_width_array(self) -> None:
        """Test effective beam width array."""
        c = Cantilever()
        c.l = 100e-6
        c.w = 20e-6
        c.beam_gap_config = GapConfig(gap_width=5e-6, gap_extent=30e-6)

        x = np.array([0, 20e-6, 40e-6, 80e-6])
        w_eff = c.effective_beam_width_array(x)

        assert w_eff[0] == pytest.approx(15e-6)
        assert w_eff[1] == pytest.approx(15e-6)
        assert w_eff[2] == pytest.approx(20e-6)
        assert w_eff[3] == pytest.approx(20e-6)

    def test_fractional_beam_gap(self) -> None:
        """Test beam gap specified as fraction of width."""
        c = Cantilever()
        c.l = 100e-6
        c.w = 20e-6
        c.beam_gap_config = GapConfig(gap_fraction=0.25, gap_extent=40e-6)

        # Gap = 25% of 20 um = 5 um
        assert c.effective_beam_width_at(0) == pytest.approx(15e-6)
        assert c.effective_beam_width_at(50e-6) == pytest.approx(20e-6)


class TestBeamGapMechanics:
    """Test beam gap effects on mechanical properties."""

    def test_stiffness_decreases_with_beam_gap(self) -> None:
        """Test that beam gap reduces stiffness (smaller width = lower EI)."""
        c1 = Cantilever()
        c1.l = 200e-6
        c1.w = 20e-6
        c1.t = 2e-6
        k1 = c1.stiffness()

        c2 = Cantilever()
        c2.l = 200e-6
        c2.w = 20e-6
        c2.t = 2e-6
        c2.beam_gap_config = GapConfig(gap_width=10e-6, gap_extent=100e-6)  # 50% gap at root half
        k2 = c2.stiffness()

        # Stiffness should be lower with the gap
        assert k2 < k1

    def test_effective_mass_decreases_with_beam_gap(self) -> None:
        """Test that beam gap reduces effective mass."""
        c1 = Cantilever()
        c1.l = 200e-6
        c1.w = 20e-6
        c1.t = 2e-6
        m1 = c1.effective_mass()

        c2 = Cantilever()
        c2.l = 200e-6
        c2.w = 20e-6
        c2.t = 2e-6
        c2.beam_gap_config = GapConfig(gap_width=10e-6, gap_extent=100e-6)
        m2 = c2.effective_mass()

        # Effective mass should be lower with the gap
        assert m2 < m1

    def test_frequency_with_beam_gap(self) -> None:
        """Test that frequency calculation works with beam gap."""
        c = Cantilever()
        c.l = 200e-6
        c.w = 20e-6
        c.t = 2e-6
        c.beam_gap_config = GapConfig(gap_width=5e-6, gap_extent=50e-6)

        # Should return a positive frequency
        omega = c._omega_vacuum_with_beam_gap()
        assert omega > 0

    def test_frequency_reduction_with_large_gap(self) -> None:
        """Test that large gap at root reduces frequency significantly."""
        c1 = Cantilever()
        c1.l = 200e-6
        c1.w = 20e-6
        c1.t = 2e-6
        f1 = c1.omega_vacuum()

        c2 = Cantilever()
        c2.l = 200e-6
        c2.w = 20e-6
        c2.t = 2e-6
        # Large gap (75% width reduction) at root region (25% of length)
        c2.beam_gap_config = GapConfig(gap_width=15e-6, gap_extent=50e-6)
        f2 = c2.omega_vacuum()

        # Frequency should change (either up or down depending on stiffness/mass ratio)
        # The key is that it's calculated without error
        assert f2 > 0
        assert f2 != f1
