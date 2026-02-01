"""Basic tests for Cantilever class."""

import pytest

from piezod import Cantilever


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
