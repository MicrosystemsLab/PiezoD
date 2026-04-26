"""Tests for PiezoresistorFromProfile.

Round-trip check: feeding `CantileverImplantation`'s simulated profile back
into `PiezoresistorFromProfile` should reproduce beta / beta1 / sheet_resistance
to within numerical-integration tolerance, since both use the same mu(n) and
P(n) lookups.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from piezod import (
    CantileverImplantation,
    CrystalOrientation,
    PiezoresistorFromProfile,
    hooge_alpha_from_anneal,
)


@pytest.fixture
def proto1_implantation() -> CantileverImplantation:
    """Proto1 L_IMP2: phosphorus 6e13 cm^-2, 30 keV, 950 C / 60 min, t_dev = 2.5 um."""
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


def _from_simulated_profile(
    c: CantileverImplantation,
    crystal_orientation: CrystalOrientation | None = None,
) -> PiezoresistorFromProfile:
    depth_m, active_cm3, _ = c.doping_profile()
    return PiezoresistorFromProfile(
        depth_m=depth_m,
        active_cm3=active_cm3,
        device_thickness_m=c.t,
        doping_type=c.doping_type,
        annealing_temp=c.annealing_temp,
        annealing_time=c.annealing_time,
        crystal_orientation=crystal_orientation,
    )


class TestConstructorValidation:
    def _kwargs(self) -> dict:
        return {
            "device_thickness_m": 2.5e-6,
            "doping_type": "phosphorus",
            "annealing_temp": 950 + 273.15,
            "annealing_time": 3600,
        }

    def test_mismatched_lengths_raise(self) -> None:
        with pytest.raises(ValueError, match="same length"):
            PiezoresistorFromProfile(
                depth_m=np.linspace(0, 1e-6, 10),
                active_cm3=np.ones(9) * 1e18,
                **self._kwargs(),
            )

    def test_non_monotonic_depth_raises(self) -> None:
        depth = np.array([0.0, 1e-7, 5e-8, 2e-7])
        active = np.full_like(depth, 1e18)
        with pytest.raises(ValueError, match="strictly increasing"):
            PiezoresistorFromProfile(depth_m=depth, active_cm3=active, **self._kwargs())

    def test_non_positive_thickness_raises(self) -> None:
        kwargs = self._kwargs()
        kwargs["device_thickness_m"] = 0.0
        with pytest.raises(ValueError, match="device_thickness_m must be positive"):
            PiezoresistorFromProfile(
                depth_m=np.linspace(0, 1e-6, 10),
                active_cm3=np.full(10, 1e18),
                **kwargs,
            )

    def test_unknown_doping_type_raises(self) -> None:
        kwargs = self._kwargs()
        kwargs["doping_type"] = "germanium"
        with pytest.raises(ValueError, match="Unknown doping_type"):
            PiezoresistorFromProfile(
                depth_m=np.linspace(0, 1e-6, 10),
                active_cm3=np.full(10, 1e18),
                **kwargs,
            )

    def test_too_few_points_after_clip_raises(self) -> None:
        depth = np.linspace(0, 1e-6, 10)
        kwargs = self._kwargs()
        kwargs["device_thickness_m"] = 1e-9  # clips everything but maybe one sample
        with pytest.raises(ValueError, match="fewer than two samples"):
            PiezoresistorFromProfile(
                depth_m=depth,
                active_cm3=np.full(10, 1e18),
                **kwargs,
            )


class TestRoundTripWithSimulatedProfile:
    """Profile-derived metrics should reproduce the lookup-table values."""

    def test_beta_matches_within_few_percent(self, proto1_implantation: CantileverImplantation) -> None:
        # The lookup tables and the trapezoidal integration use slightly
        # different conventions (the lookup is built from a separate process
        # simulator post-process). Agree to ~5%.
        c = proto1_implantation
        pr = _from_simulated_profile(c)
        lookup_beta = c.doping_process_metrics().beta
        profile_beta = pr.beta()
        assert profile_beta == pytest.approx(lookup_beta, rel=0.05)

    def test_alpha_h_matches_exactly(self, proto1_implantation: CantileverImplantation) -> None:
        c = proto1_implantation
        pr = _from_simulated_profile(c)
        assert pr.alpha() == pytest.approx(c.alpha(), rel=1e-12)

    def test_sheet_resistance_matches_within_tolerance(self, proto1_implantation: CantileverImplantation) -> None:
        c = proto1_implantation
        pr = _from_simulated_profile(c)
        # The lookup Rs comes from TSUPREM/DopeDealer's own integration, so
        # exact agreement isn't expected; ~25% accommodates table-vs-profile
        # differences while still validating the formula and units.
        assert pr.sheet_resistance() == pytest.approx(c.sheet_resistance(), rel=0.25)

    def test_peak_concentration_matches(self, proto1_implantation: CantileverImplantation) -> None:
        c = proto1_implantation
        pr = _from_simulated_profile(c)
        depth_m, active_cm3, _ = c.doping_profile()
        assert pr.doping_process_metrics().peak_concentration_cm3 == pytest.approx(np.max(active_cm3))

    def test_metrics_dataclass_shape_matches_implantation(self, proto1_implantation: CantileverImplantation) -> None:
        c = proto1_implantation
        pr = _from_simulated_profile(c)
        impl_fields = set(c.doping_process_metrics().__dataclass_fields__.keys())
        prof_fields = set(pr.doping_process_metrics().__dataclass_fields__.keys())
        assert impl_fields == prof_fields


class TestPiAndDrOverR:
    def test_pi_low_doping_matches_smith_kanda(self, proto1_implantation: CantileverImplantation) -> None:
        pr = _from_simulated_profile(proto1_implantation)
        assert pr.pi_longitudinal() == pytest.approx(-102.2e-11, rel=1e-12)
        assert pr.pi_transverse() == pytest.approx(53.4e-11, rel=1e-12)

    def test_explicit_orientation_override(self, proto1_implantation: CantileverImplantation) -> None:
        pr = _from_simulated_profile(proto1_implantation, crystal_orientation=CrystalOrientation.WAFER_100_DIR_110)
        assert pr.pi_longitudinal() == pytest.approx(-31.2e-11, rel=1e-12)
        assert pr.pi_transverse() == pytest.approx(-17.6e-11, rel=1e-12)

    def test_dr_over_r_zero_stress(self, proto1_implantation: CantileverImplantation) -> None:
        pr = _from_simulated_profile(proto1_implantation)
        assert pr.dr_over_r(0.0, 0.0) == pytest.approx(0.0)

    def test_dr_over_r_formula(self, proto1_implantation: CantileverImplantation) -> None:
        pr = _from_simulated_profile(proto1_implantation)
        beta = pr.beta()
        pi_l = pr.pi_longitudinal()
        pi_t = pr.pi_transverse()
        sigma_l, sigma_t = 3e6, -1e6
        assert pr.dr_over_r(sigma_l, sigma_t) == pytest.approx(beta * (pi_l * sigma_l + pi_t * sigma_t))

    def test_dr_over_r_from_tensor_rotation(self, proto1_implantation: CantileverImplantation) -> None:
        pr = _from_simulated_profile(proto1_implantation)
        sigma = 5e6
        # theta=0: along x, equals dr_over_r(sxx, 0)
        assert pr.dr_over_r_from_tensor(sigma, 0.0, 0.0, 0.0) == pytest.approx(pr.dr_over_r(sigma, 0.0))
        # theta=pi/2: along y, equals dr_over_r(0, sxx)
        assert pr.dr_over_r_from_tensor(sigma, 0.0, 0.0, math.pi / 2) == pytest.approx(pr.dr_over_r(0.0, sigma))

    def test_vectorized_dr_over_r(self, proto1_implantation: CantileverImplantation) -> None:
        pr = _from_simulated_profile(proto1_implantation)
        sigma_l = np.array([0.0, 1e6, 2e6])
        sigma_t = np.array([0.0, -5e5, 1e6])
        result = pr.dr_over_r(sigma_l, sigma_t)
        assert isinstance(result, np.ndarray)
        beta = pr.beta()
        pi_l = pr.pi_longitudinal()
        pi_t = pr.pi_transverse()
        np.testing.assert_allclose(result, beta * (pi_l * sigma_l + pi_t * sigma_t))


class TestSubstrateBackground:
    """Verify the substrate_background_cm3 kwarg drives net active carriers."""

    def _kwargs(self) -> dict:
        return {
            "device_thickness_m": 2.5e-6,
            "doping_type": "phosphorus",
            "annealing_temp": 950 + 273.15,
            "annealing_time": 3600,
        }

    def test_default_substrate_background(self) -> None:
        depth = np.linspace(0, 1e-6, 50)
        active = np.full_like(depth, 1e18)
        pr = PiezoresistorFromProfile(depth_m=depth, active_cm3=active, **self._kwargs())
        assert pr.substrate_background_cm3 == pytest.approx(1e15)

    def test_negative_background_raises(self) -> None:
        depth = np.linspace(0, 1e-6, 50)
        active = np.full_like(depth, 1e18)
        with pytest.raises(ValueError, match="substrate_background_cm3 must be nonnegative"):
            PiezoresistorFromProfile(
                depth_m=depth,
                active_cm3=active,
                substrate_background_cm3=-1.0,
                **self._kwargs(),
            )

    def test_net_active_subtracts_background(self) -> None:
        depth = np.linspace(0, 1e-6, 50)
        active = np.full_like(depth, 1e18)
        pr = PiezoresistorFromProfile(
            depth_m=depth,
            active_cm3=active,
            substrate_background_cm3=2e16,
            **self._kwargs(),
        )
        np.testing.assert_allclose(pr.net_active_cm3(), np.full_like(depth, 1e18 - 2e16))

    def test_total_cm3_floors_at_substrate(self) -> None:
        """`total_cm3()` floors at substrate background where the dopant decays."""
        depth = np.linspace(0, 1e-6, 200)
        # Linear ramp 1e19 -> 1e14
        active = np.linspace(1e19, 1e14, 200)
        pr = PiezoresistorFromProfile(
            depth_m=depth,
            active_cm3=active,
            substrate_background_cm3=1e16,
            **self._kwargs(),
        )
        total = pr.total_cm3()
        assert total[-1] == pytest.approx(1e16)
        # Where dopant > substrate, total == dopant.
        above = active > 1e16
        np.testing.assert_allclose(total[above], active[above])

    def test_junction_depth_tracks_background(self) -> None:
        # Linear ramp from 1e19 down to 1e14 across 1 um
        depth = np.linspace(0, 1e-6, 501)
        active = np.linspace(1e19, 1e14, 501)
        pr = PiezoresistorFromProfile(
            depth_m=depth,
            active_cm3=active,
            substrate_background_cm3=1e17,
            **self._kwargs(),
        )
        # Junction is where the dopant equals 1e17 -- by linear interp,
        # x = (1e19 - 1e17) / (1e19 - 1e14) * 1um ~= 0.99 um.
        expected = (1e19 - 1e17) / (1e19 - 1e14) * 1e-6
        assert pr.junction_depth == pytest.approx(expected, rel=1e-3)

    def test_higher_background_increases_sheet_resistance(self) -> None:
        depth = np.linspace(0, 1e-6, 50)
        active = np.full_like(depth, 1e18)
        pr_low = PiezoresistorFromProfile(
            depth_m=depth, active_cm3=active, substrate_background_cm3=1e15, **self._kwargs()
        )
        pr_high = PiezoresistorFromProfile(
            depth_m=depth, active_cm3=active, substrate_background_cm3=5e17, **self._kwargs()
        )
        assert pr_high.sheet_resistance() > pr_low.sheet_resistance()


class TestHoogeAlphaHelper:
    def test_helper_matches_implantation_alpha(self, proto1_implantation: CantileverImplantation) -> None:
        c = proto1_implantation
        helper = hooge_alpha_from_anneal(c.doping_type, c.annealing_temp, c.annealing_time)
        assert helper == pytest.approx(c.alpha(), rel=1e-12)

    def test_unknown_doping_type_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown doping_type"):
            hooge_alpha_from_anneal("germanium", 1200.0, 3600.0)
