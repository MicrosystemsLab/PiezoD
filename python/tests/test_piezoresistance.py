"""Tests for crystal-orientation pi coefficients and dR/R helpers."""

from __future__ import annotations

import math

import numpy as np
import pytest

from piezod import (
    Cantilever,
    CantileverImplantation,
    CrystalOrientation,
    default_orientation,
    pi_low_doping,
    rotate_in_plane_stress,
)


class TestPiLowDoping:
    """Smith (1954) / Kanda (1982) reference values for pi_l, pi_t."""

    @pytest.mark.parametrize(
        ("doping_type", "orientation", "pi_l_tpa", "pi_t_tpa"),
        [
            ("phosphorus", CrystalOrientation.WAFER_100_DIR_100, -102.2, 53.4),
            ("phosphorus", CrystalOrientation.WAFER_100_DIR_110, -31.2, -17.6),
            ("arsenic", CrystalOrientation.WAFER_100_DIR_100, -102.2, 53.4),
            ("boron", CrystalOrientation.WAFER_100_DIR_100, 6.6, -1.1),
            ("boron", CrystalOrientation.WAFER_100_DIR_110, 71.8, -66.3),
        ],
    )
    def test_reference_values(
        self,
        doping_type: str,
        orientation: CrystalOrientation,
        pi_l_tpa: float,
        pi_t_tpa: float,
    ) -> None:
        pi_l, pi_t = pi_low_doping(doping_type, orientation)
        assert pi_l == pytest.approx(pi_l_tpa * 1e-11, rel=1e-12)
        assert pi_t == pytest.approx(pi_t_tpa * 1e-11, rel=1e-12)

    def test_unknown_doping_type_raises(self) -> None:
        with pytest.raises(ValueError, match="No piezoresistance coefficients"):
            pi_low_doping("germanium", CrystalOrientation.WAFER_100_DIR_100)

    def test_orientation_must_be_enum(self) -> None:
        with pytest.raises(TypeError, match="must be a CrystalOrientation enum"):
            pi_low_doping("phosphorus", "100_along_100")  # type: ignore[arg-type]


class TestDefaultOrientation:
    """Per-doping default orientations match the most-sensitive direction."""

    def test_phosphorus_defaults_to_100(self) -> None:
        assert default_orientation("phosphorus") == CrystalOrientation.WAFER_100_DIR_100

    def test_arsenic_defaults_to_100(self) -> None:
        assert default_orientation("arsenic") == CrystalOrientation.WAFER_100_DIR_100

    def test_boron_defaults_to_110(self) -> None:
        assert default_orientation("boron") == CrystalOrientation.WAFER_100_DIR_110

    def test_unknown_dopant_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown doping_type"):
            default_orientation("xenon")


class TestCantileverPiAccessors:
    """Cantilever exposes signed pi_l / pi_t and a magnitude alias."""

    def _cantilever(self, doping_type: str) -> Cantilever:
        c = Cantilever()
        c.doping_type = doping_type
        return c

    def test_phosphorus_default(self) -> None:
        c = self._cantilever("phosphorus")
        assert c.pi_longitudinal() == pytest.approx(-102.2e-11, rel=1e-12)
        assert c.pi_transverse() == pytest.approx(53.4e-11, rel=1e-12)

    def test_boron_default(self) -> None:
        c = self._cantilever("boron")
        assert c.pi_longitudinal() == pytest.approx(71.8e-11, rel=1e-12)
        assert c.pi_transverse() == pytest.approx(-66.3e-11, rel=1e-12)

    def test_explicit_orientation_override(self) -> None:
        c = self._cantilever("boron")
        c.crystal_orientation = CrystalOrientation.WAFER_100_DIR_100
        assert c.pi_longitudinal() == pytest.approx(6.6e-11, rel=1e-12)
        assert c.pi_transverse() == pytest.approx(-1.1e-11, rel=1e-12)

    def test_max_piezoresistance_factor_equals_abs_pi_longitudinal(self) -> None:
        for doping_type in ("phosphorus", "arsenic", "boron"):
            c = self._cantilever(doping_type)
            assert c.max_piezoresistance_factor() == pytest.approx(abs(c.pi_longitudinal()))


class TestRotateInPlaneStress:
    """In-plane stress rotation."""

    def test_zero_angle_is_identity(self) -> None:
        sxx, syy, sxy = 100e6, -50e6, 30e6
        sigma_l, sigma_t = rotate_in_plane_stress(sxx, syy, sxy, 0.0)
        assert float(sigma_l) == pytest.approx(sxx)
        assert float(sigma_t) == pytest.approx(syy)

    def test_pi_over_two_swaps_l_and_t(self) -> None:
        sxx, syy = 100e6, -50e6
        sigma_l, sigma_t = rotate_in_plane_stress(sxx, syy, 0.0, math.pi / 2)
        assert float(sigma_l) == pytest.approx(syy)
        assert float(sigma_t) == pytest.approx(sxx)

    def test_pure_shear_at_45_deg_is_principal(self) -> None:
        # Pure shear sxy = tau in xy-frame -> principal stresses tau, -tau at 45 deg.
        sigma_l, sigma_t = rotate_in_plane_stress(0.0, 0.0, 1e6, math.pi / 4)
        assert float(sigma_l) == pytest.approx(1e6)
        assert float(sigma_t) == pytest.approx(-1e6)

    def test_vectorized_inputs(self) -> None:
        sxx = np.array([100e6, 200e6, 300e6])
        syy = np.array([-50e6, 0.0, 50e6])
        sxy = np.zeros(3)
        sigma_l, sigma_t = rotate_in_plane_stress(sxx, syy, sxy, math.pi / 2)
        np.testing.assert_allclose(sigma_l, syy, atol=1e-6)
        np.testing.assert_allclose(sigma_t, sxx, atol=1e-6)


class TestCantileverDrOverR:
    """Cantilever.dr_over_r and dr_over_r_from_tensor wiring."""

    @pytest.fixture
    def lightly_doped_phosphorus(self) -> CantileverImplantation:
        # Beta is profile-dependent; at low dose the implant beta still differs
        # from 1 substantially. Tests below check the formula at a known beta.
        return CantileverImplantation(
            freq_min=1e3,
            freq_max=1e6,
            l=200e-6,
            w=20e-6,
            t=2.5e-6,
            l_pr_ratio=0.1,
            v_bridge=1.8,
            doping_type="phosphorus",
            annealing_time=3600,
            annealing_temp=950 + 273.15,
            annealing_type="inert",
            implantation_energy=30,
            implantation_dose=6e13,
            lookup_source="dopedealer",
        )

    def test_zero_stress_gives_zero(self, lightly_doped_phosphorus: CantileverImplantation) -> None:
        c = lightly_doped_phosphorus
        assert c.dr_over_r(0.0, 0.0) == pytest.approx(0.0)
        assert c.dr_over_r_from_tensor(0.0, 0.0, 0.0, 0.0) == pytest.approx(0.0)

    def test_formula_matches_explicit_calculation(self, lightly_doped_phosphorus: CantileverImplantation) -> None:
        c = lightly_doped_phosphorus
        beta = c.beta()
        pi_l = c.pi_longitudinal()
        pi_t = c.pi_transverse()
        sigma_l = 1e6
        sigma_t = -2e6
        expected = beta * (pi_l * sigma_l + pi_t * sigma_t)
        assert c.dr_over_r(sigma_l, sigma_t) == pytest.approx(expected, rel=1e-12)

    def test_tensor_zero_angle_equals_dr_over_r(self, lightly_doped_phosphorus: CantileverImplantation) -> None:
        c = lightly_doped_phosphorus
        sigma = 5e6
        assert c.dr_over_r_from_tensor(sigma, 0.0, 0.0, 0.0) == pytest.approx(c.dr_over_r(sigma, 0.0))

    def test_tensor_pi_over_two_swaps_inputs(self, lightly_doped_phosphorus: CantileverImplantation) -> None:
        c = lightly_doped_phosphorus
        sigma = 5e6
        rotated = c.dr_over_r_from_tensor(sigma, 0.0, 0.0, math.pi / 2)
        assert rotated == pytest.approx(c.dr_over_r(0.0, sigma))

    def test_vectorized_dr_over_r(self, lightly_doped_phosphorus: CantileverImplantation) -> None:
        c = lightly_doped_phosphorus
        sigma_l = np.array([0.0, 1e6, 2e6])
        sigma_t = np.array([0.0, -1e6, 0.0])
        result = c.dr_over_r(sigma_l, sigma_t)
        assert isinstance(result, np.ndarray)
        beta = c.beta()
        pi_l = c.pi_longitudinal()
        pi_t = c.pi_transverse()
        np.testing.assert_allclose(result, beta * (pi_l * sigma_l + pi_t * sigma_t))

    def test_low_doping_known_value(self) -> None:
        # n-Si <100>, beta forced to 1 by setting profile manually -- verify
        # Smith/Kanda magnitude. dr_over_r(1e6, 0) = pi_l_lowdoping * 1e6 = -1.022e-3.
        c = Cantilever()
        c.doping_type = "phosphorus"

        # Force beta = 1 by monkey-patching for this isolated check.
        c.beta = lambda: 1.0  # type: ignore[method-assign]
        assert c.dr_over_r(1e6, 0.0) == pytest.approx(-1.022e-3, rel=1e-3)
