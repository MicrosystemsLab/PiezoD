"""Tests for CantileverPoly class."""

import numpy as np
import pytest

from piezod.cantilever_poly import CantileverPoly, Material


class TestMaterialEnum:
    """Test Material enum."""

    def test_material_values(self) -> None:
        """Test all material values are correct."""
        assert Material.SI.value == "si"
        assert Material.POLY.value == "poly"
        assert Material.TI.value == "ti"
        assert Material.AL.value == "al"
        assert Material.OXIDE.value == "oxide"

    def test_material_enum_members(self) -> None:
        """Test all expected materials exist."""
        materials = [m.value for m in Material]
        assert "si" in materials
        assert "poly" in materials
        assert "ti" in materials
        assert "al" in materials
        assert "oxide" in materials


class TestCantileverPolyInitialization:
    """Test initialization and basic properties."""

    def test_default_initialization(self) -> None:
        """Test cantilever initializes with default values."""
        cant = CantileverPoly()

        assert cant.freq_min == 1.0
        assert cant.freq_max == 1e3
        assert cant.l == 100e-6
        assert cant.w == 10e-6
        assert cant.t_top == 100e-9
        assert cant.t_mid == 500e-9
        assert cant.matl_top == Material.POLY
        assert cant.matl_mid == Material.SI
        assert cant.matl_bot == Material.POLY
        assert cant.l_pr_ratio == 0.5
        assert cant.v_bridge == 1.0
        assert cant.dopant_concentration == 1e19
        assert cant.number_of_piezoresistors == 2
        assert cant.number_of_piezoresistors_on_cantilever == 2

    def test_custom_initialization(self) -> None:
        """Test cantilever initializes with custom values."""
        cant = CantileverPoly(
            freq_min=10.0,
            freq_max=1e4,
            l=200e-6,
            w=20e-6,
            t_top=200e-9,
            t_mid=1e-6,
            t_bot=200e-9,
            matl_top=Material.TI,
            matl_mid=Material.OXIDE,
            matl_bot=Material.AL,
            l_pr_ratio=0.3,
            v_bridge=2.0,
            dopant_concentration=5e19,
            number_of_piezoresistors=4,
            number_of_piezoresistors_on_cantilever=1,
        )

        assert cant.freq_min == 10.0
        assert cant.freq_max == 1e4
        assert cant.l == 200e-6
        assert cant.w == 20e-6
        assert cant.matl_top == Material.TI
        assert cant.matl_mid == Material.OXIDE
        # When number_of_piezoresistors_on_cantilever != 2, t_bot is not forced
        assert cant.number_of_piezoresistors_on_cantilever == 1
        assert cant.t_bot == 200e-9

    def test_symmetric_layers_enforced(self) -> None:
        """Test t_bot is forced equal to t_top for double PR config."""
        cant = CantileverPoly(
            t_top=150e-9,
            t_bot=300e-9,  # Different from t_top
            number_of_piezoresistors_on_cantilever=2,
        )
        # Should be forced equal
        assert cant.t_bot == cant.t_top
        assert cant.t_bot == 150e-9

    def test_asymmetric_layers_allowed(self) -> None:
        """Test t_bot can differ from t_top for single PR config."""
        cant = CantileverPoly(
            t_top=150e-9,
            t_bot=300e-9,
            number_of_piezoresistors_on_cantilever=1,
        )
        # Should remain different
        assert cant.t_top == 150e-9
        assert cant.t_bot == 300e-9

    def test_is_standalone_class(self) -> None:
        """Test CantileverPoly does NOT inherit from Cantilever base class."""
        from piezod.cantilever import Cantilever

        cant = CantileverPoly()
        assert not isinstance(cant, Cantilever)


class TestProperties:
    """Test computed properties."""

    def test_l_pr_property(self) -> None:
        """Test piezoresistor length property."""
        cant = CantileverPoly(l=100e-6, l_pr_ratio=0.5)
        assert cant.l_pr == pytest.approx(50e-6)

    def test_w_pr_property(self) -> None:
        """Test piezoresistor width property."""
        cant = CantileverPoly(w=10e-6)
        assert cant.w_pr == pytest.approx(5e-6)

    def test_t_total_property(self) -> None:
        """Test total thickness property."""
        cant = CantileverPoly(t_top=100e-9, t_mid=500e-9, t_bot=100e-9)
        assert cant.t_total == pytest.approx(700e-9)


class TestMaterialProperties:
    """Test material property lookups."""

    def test_get_elastic_modulus(self) -> None:
        """Test elastic modulus lookup for all materials."""
        cant = CantileverPoly()

        assert cant._get_elastic_modulus(Material.SI) == 169e9
        assert cant._get_elastic_modulus(Material.POLY) == 150e9
        assert cant._get_elastic_modulus(Material.TI) == 85e9
        assert cant._get_elastic_modulus(Material.AL) == 70e9
        assert cant._get_elastic_modulus(Material.OXIDE) == 70e9

    def test_get_density(self) -> None:
        """Test density lookup for all materials."""
        cant = CantileverPoly()

        assert cant._get_density(Material.SI) == 2330
        assert cant._get_density(Material.POLY) == 2330
        assert cant._get_density(Material.TI) == 4506
        assert cant._get_density(Material.AL) == 2700
        assert cant._get_density(Material.OXIDE) == 2634


class TestBeamMechanics:
    """Test beam mechanics calculations."""

    def test_neutral_axis_returns_float(self) -> None:
        """Test neutral axis returns a float."""
        cant = CantileverPoly()
        na = cant.neutral_axis()
        assert isinstance(na, float)

    def test_neutral_axis_positive(self) -> None:
        """Test neutral axis is positive (above bottom)."""
        cant = CantileverPoly()
        na = cant.neutral_axis()
        assert na > 0

    def test_neutral_axis_within_thickness(self) -> None:
        """Test neutral axis is within total thickness."""
        cant = CantileverPoly()
        na = cant.neutral_axis()
        assert na < cant.t_total

    def test_neutral_axis_symmetric_layers(self) -> None:
        """Test neutral axis is at center for symmetric structure."""
        cant = CantileverPoly(
            t_top=100e-9,
            t_mid=500e-9,
            t_bot=100e-9,
            matl_top=Material.POLY,
            matl_mid=Material.SI,
            matl_bot=Material.POLY,
        )
        na = cant.neutral_axis()
        # For symmetric structure, neutral axis should be near center
        center = cant.t_total / 2
        assert na == pytest.approx(center, rel=0.1)

    def test_normalized_curvature_returns_float(self) -> None:
        """Test normalized curvature returns a float."""
        cant = CantileverPoly()
        nc = cant.normalized_curvature()
        assert isinstance(nc, float)

    def test_normalized_curvature_positive(self) -> None:
        """Test normalized curvature is positive."""
        cant = CantileverPoly()
        nc = cant.normalized_curvature()
        assert nc > 0

    def test_stiffness_returns_float(self) -> None:
        """Test stiffness returns a float."""
        cant = CantileverPoly()
        k = cant.stiffness()
        assert isinstance(k, float)

    def test_stiffness_positive(self) -> None:
        """Test stiffness is positive."""
        cant = CantileverPoly()
        k = cant.stiffness()
        assert k > 0

    def test_stiffness_increases_with_thickness(self) -> None:
        """Test stiffness increases with layer thickness."""
        cant1 = CantileverPoly(t_mid=500e-9)
        cant2 = CantileverPoly(t_mid=1000e-9)

        k1 = cant1.stiffness()
        k2 = cant2.stiffness()

        assert k2 > k1

    def test_stiffness_decreases_with_length(self) -> None:
        """Test stiffness decreases with length (k ~ 1/L^3)."""
        cant1 = CantileverPoly(l=100e-6)
        cant2 = CantileverPoly(l=200e-6)

        k1 = cant1.stiffness()
        k2 = cant2.stiffness()

        # k ~ 1/L^3, so doubling L should reduce k by factor of 8
        assert k2 < k1
        assert k1 / k2 == pytest.approx(8, rel=0.2)


class TestResonantFrequency:
    """Test resonant frequency calculations."""

    def test_resonant_frequency_returns_float(self) -> None:
        """Test resonant frequency returns a float."""
        cant = CantileverPoly()
        f = cant.resonant_frequency()
        assert isinstance(f, float)

    def test_resonant_frequency_positive(self) -> None:
        """Test resonant frequency is positive."""
        cant = CantileverPoly()
        f = cant.resonant_frequency()
        assert f > 0

    def test_resonant_frequency_reasonable_range(self) -> None:
        """Test resonant frequency is in reasonable MEMS range."""
        cant = CantileverPoly()
        f = cant.resonant_frequency()
        # Typical MEMS cantilevers: 1 kHz - 1 MHz
        assert 1e3 < f < 1e6

    def test_resonant_frequency_decreases_with_length(self) -> None:
        """Test resonant frequency decreases with length (f ~ 1/L^2)."""
        cant1 = CantileverPoly(l=100e-6)
        cant2 = CantileverPoly(l=200e-6)

        f1 = cant1.resonant_frequency()
        f2 = cant2.resonant_frequency()

        assert f2 < f1
        # f ~ 1/L^2, so doubling L should reduce f by factor of 4
        assert f1 / f2 == pytest.approx(4, rel=0.3)


class TestMobility:
    """Test carrier mobility calculations."""

    def test_mobility_returns_float(self) -> None:
        """Test mobility returns a float."""
        cant = CantileverPoly()
        mu = cant.mobility(1e19)
        assert isinstance(mu, float)

    def test_mobility_positive(self) -> None:
        """Test mobility is positive."""
        cant = CantileverPoly()
        mu = cant.mobility(1e19)
        assert mu > 0

    def test_mobility_decreases_with_concentration(self) -> None:
        """Test mobility decreases at higher doping concentrations."""
        cant = CantileverPoly()
        mu1 = cant.mobility(1e17)
        mu2 = cant.mobility(1e20)

        assert mu1 > mu2

    def test_mobility_poly_correction(self) -> None:
        """Test 50% poly correction is applied."""
        cant = CantileverPoly()
        mu = cant.mobility(1e17)
        # At low concentration, mobility should be ~707 (half of ~1414)
        assert 300 < mu < 750


class TestSheetResistance:
    """Test sheet resistance calculations."""

    def test_sheet_resistance_returns_float(self) -> None:
        """Test sheet_resistance returns a float."""
        cant = CantileverPoly()
        Rs = cant.sheet_resistance()
        assert isinstance(Rs, (float, np.floating))

    def test_sheet_resistance_positive(self) -> None:
        """Test sheet resistance is positive."""
        cant = CantileverPoly()
        Rs = cant.sheet_resistance()
        assert Rs > 0

    def test_sheet_resistance_poly(self) -> None:
        """Test sheet resistance for poly-Si piezoresistor."""
        cant = CantileverPoly(matl_top=Material.POLY)
        Rs = cant.sheet_resistance()
        # Typical poly-Si sheet resistance: 10-1000 ohm/sq
        assert 1 < Rs < 1e4

    def test_sheet_resistance_ti(self) -> None:
        """Test sheet resistance for Ti piezoresistor."""
        cant = CantileverPoly(matl_top=Material.TI, t_top=100e-9)
        Rs = cant.sheet_resistance()
        # Ti resistivity = 42e-6 ohm-cm, t = 100nm
        expected = 42e-6 / (100e-9 * 100)
        assert Rs == pytest.approx(expected, rel=0.01)

    def test_sheet_resistance_al(self) -> None:
        """Test sheet resistance for Al piezoresistor."""
        cant = CantileverPoly(matl_top=Material.AL, t_top=100e-9)
        Rs = cant.sheet_resistance()
        # Al resistivity = 2.7e-6 ohm-cm, t = 100nm
        expected = 2.7e-6 / (100e-9 * 100)
        assert Rs == pytest.approx(expected, rel=0.01)


class TestResistance:
    """Test resistance calculations."""

    def test_resistance_returns_float(self) -> None:
        """Test resistance returns a float."""
        cant = CantileverPoly()
        R = cant.resistance()
        assert isinstance(R, (float, np.floating))

    def test_resistance_positive(self) -> None:
        """Test resistance is positive."""
        cant = CantileverPoly()
        R = cant.resistance()
        assert R > 0

    def test_resistor_length(self) -> None:
        """Test resistor length calculation."""
        cant = CantileverPoly(l=100e-6, l_pr_ratio=0.5)
        rl = cant.resistor_length()
        # 2 * l_pr = 2 * 50e-6 = 100e-6
        assert rl == pytest.approx(100e-6)

    def test_gamma_value(self) -> None:
        """Test gamma returns expected value."""
        cant = CantileverPoly()
        assert cant.gamma() == 0.9


class TestNumberOfCarriers:
    """Test carrier count calculations."""

    def test_number_of_carriers_returns_float(self) -> None:
        """Test number_of_carriers returns a float."""
        cant = CantileverPoly()
        N = cant.number_of_carriers()
        assert isinstance(N, (float, np.floating))

    def test_number_of_carriers_positive(self) -> None:
        """Test number of carriers is positive."""
        cant = CantileverPoly()
        N = cant.number_of_carriers()
        assert N > 0

    def test_number_of_carriers_poly(self) -> None:
        """Test carrier count for poly-Si."""
        cant = CantileverPoly(matl_top=Material.POLY, dopant_concentration=1e19)
        N = cant.number_of_carriers()
        # Should be positive and scale with dopant concentration
        assert N > 0

    def test_number_of_carriers_metal(self) -> None:
        """Test carrier count for metals uses fixed value."""
        cant = CantileverPoly(matl_top=Material.TI)
        N = cant.number_of_carriers()
        # Uses 2e21 carriers/cm^3, should be positive
        assert N > 0


class TestNoise:
    """Test noise calculations."""

    def test_hooge_PSD_returns_array(self) -> None:
        """Test hooge_PSD returns numpy array."""
        cant = CantileverPoly()
        freq = np.logspace(0, 3, 100)
        psd = cant.hooge_PSD(freq)
        assert isinstance(psd, np.ndarray)
        assert len(psd) == len(freq)

    def test_hooge_PSD_positive(self) -> None:
        """Test hooge PSD is positive."""
        cant = CantileverPoly()
        freq = np.logspace(0, 3, 100)
        psd = cant.hooge_PSD(freq)
        assert np.all(psd > 0)

    def test_hooge_PSD_decreases_with_frequency(self) -> None:
        """Test hooge PSD is 1/f."""
        cant = CantileverPoly()
        freq = np.array([1, 10, 100])
        psd = cant.hooge_PSD(freq)
        # 1/f noise: PSD at 10 Hz should be 10x lower than at 1 Hz
        assert psd[0] / psd[1] == pytest.approx(10, rel=0.01)

    def test_hooge_integrated_returns_float(self) -> None:
        """Test hooge_integrated returns a float."""
        cant = CantileverPoly()
        noise = cant.hooge_integrated()
        assert isinstance(noise, float)

    def test_hooge_integrated_positive(self) -> None:
        """Test integrated hooge noise is positive."""
        cant = CantileverPoly()
        noise = cant.hooge_integrated()
        assert noise > 0

    def test_johnson_PSD_returns_array(self) -> None:
        """Test johnson_PSD returns numpy array."""
        cant = CantileverPoly()
        freq = np.logspace(0, 3, 100)
        psd = cant.johnson_PSD(freq)
        assert isinstance(psd, np.ndarray)
        assert len(psd) == len(freq)

    def test_johnson_PSD_flat(self) -> None:
        """Test johnson PSD is flat (white noise)."""
        cant = CantileverPoly()
        freq = np.logspace(0, 3, 100)
        psd = cant.johnson_PSD(freq)
        # All values should be equal
        assert np.allclose(psd, psd[0])

    def test_johnson_integrated_returns_float(self) -> None:
        """Test johnson_integrated returns a float."""
        cant = CantileverPoly()
        noise = cant.johnson_integrated()
        assert isinstance(noise, float)

    def test_thermo_integrated_returns_float(self) -> None:
        """Test thermo_integrated returns a float."""
        cant = CantileverPoly()
        noise = cant.thermo_integrated()
        assert isinstance(noise, float)

    def test_amplifier_PSD_returns_array(self) -> None:
        """Test amplifier_PSD returns numpy array."""
        cant = CantileverPoly()
        freq = np.logspace(0, 3, 100)
        psd = cant.amplifier_PSD(freq)
        assert isinstance(psd, np.ndarray)
        assert len(psd) == len(freq)

    def test_amplifier_integrated_returns_float(self) -> None:
        """Test amplifier_integrated returns a float."""
        cant = CantileverPoly()
        noise = cant.amplifier_integrated()
        assert isinstance(noise, float)

    def test_actuator_noise_integrated_returns_float(self) -> None:
        """Test actuator_noise_integrated returns a float."""
        cant = CantileverPoly()
        noise = cant.actuator_noise_integrated()
        assert isinstance(noise, float)

    def test_integrated_noise_returns_float(self) -> None:
        """Test integrated_noise returns a float."""
        cant = CantileverPoly()
        noise = cant.integrated_noise()
        assert isinstance(noise, float)

    def test_integrated_noise_positive(self) -> None:
        """Test total integrated noise is positive."""
        cant = CantileverPoly()
        noise = cant.integrated_noise()
        assert noise > 0

    def test_knee_frequency_returns_float(self) -> None:
        """Test knee_frequency returns a float."""
        cant = CantileverPoly()
        fk = cant.knee_frequency()
        assert isinstance(fk, float)

    def test_knee_frequency_positive(self) -> None:
        """Test knee frequency is positive."""
        cant = CantileverPoly()
        fk = cant.knee_frequency()
        assert fk > 0


class TestSensitivity:
    """Test sensitivity calculations."""

    def test_piezo_coefficient_poly(self) -> None:
        """Test piezoresistive coefficient for poly-Si."""
        cant = CantileverPoly(matl_top=Material.POLY, dopant_concentration=1e19)
        pi = cant.piezo_coefficient()
        # Should be positive and in reasonable range
        assert pi > 0
        assert pi < 1e-9  # Max ~10^-10 Pa^-1

    def test_piezo_coefficient_ti(self) -> None:
        """Test piezoresistive coefficient for Ti (geometric effect only)."""
        cant = CantileverPoly(matl_top=Material.TI)
        pi = cant.piezo_coefficient()
        # (1 + 2*0.35) / 85e9 ~ 2e-11
        expected = (1 + 2 * 0.35) / 85e9
        assert pi == pytest.approx(expected, rel=0.01)

    def test_piezo_coefficient_al(self) -> None:
        """Test piezoresistive coefficient for Al (geometric effect only)."""
        cant = CantileverPoly(matl_top=Material.AL)
        pi = cant.piezo_coefficient()
        # (1 + 2*0.35) / 70e9 ~ 2.4e-11
        expected = (1 + 2 * 0.35) / 70e9
        assert pi == pytest.approx(expected, rel=0.01)

    def test_piezo_coefficient_richter_model(self) -> None:
        """Test Richter's model reduces coefficient at high doping."""
        cant1 = CantileverPoly(matl_top=Material.POLY, dopant_concentration=1e17)
        cant2 = CantileverPoly(matl_top=Material.POLY, dopant_concentration=1e20)

        pi1 = cant1.piezo_coefficient()
        pi2 = cant2.piezo_coefficient()

        # Higher doping should reduce piezoresistance factor
        assert pi1 > pi2

    def test_force_sensitivity_returns_float(self) -> None:
        """Test force_sensitivity returns a float."""
        cant = CantileverPoly()
        S = cant.force_sensitivity()
        assert isinstance(S, float)

    def test_force_sensitivity_positive(self) -> None:
        """Test force sensitivity is positive."""
        cant = CantileverPoly()
        S = cant.force_sensitivity()
        assert S > 0

    def test_force_sensitivity_increases_with_voltage(self) -> None:
        """Test force sensitivity increases with bridge voltage."""
        cant1 = CantileverPoly(v_bridge=1.0)
        cant2 = CantileverPoly(v_bridge=2.0)

        S1 = cant1.force_sensitivity()
        S2 = cant2.force_sensitivity()

        assert S2 > S1
        # Should scale linearly with voltage
        assert pytest.approx(2.0, rel=0.1) == S2 / S1


class TestPowerAndResolution:
    """Test power dissipation and resolution calculations."""

    def test_power_dissipation_returns_float(self) -> None:
        """Test power_dissipation returns a float."""
        cant = CantileverPoly()
        P = cant.power_dissipation()
        assert isinstance(P, (float, np.floating))

    def test_power_dissipation_positive(self) -> None:
        """Test power dissipation is positive."""
        cant = CantileverPoly()
        P = cant.power_dissipation()
        assert P > 0

    def test_power_scales_with_voltage_squared(self) -> None:
        """Test power scales with V^2."""
        cant1 = CantileverPoly(v_bridge=1.0)
        cant2 = CantileverPoly(v_bridge=2.0)

        P1 = cant1.power_dissipation()
        P2 = cant2.power_dissipation()

        # P ~ V^2, so doubling V should quadruple P
        assert pytest.approx(4.0, rel=0.1) == P2 / P1

    def test_force_resolution_returns_float(self) -> None:
        """Test force_resolution returns a float."""
        cant = CantileverPoly()
        F = cant.force_resolution()
        assert isinstance(F, float)

    def test_force_resolution_positive(self) -> None:
        """Test force resolution is positive."""
        cant = CantileverPoly()
        F = cant.force_resolution()
        assert F > 0

    def test_displacement_resolution_returns_float(self) -> None:
        """Test displacement_resolution returns a float."""
        cant = CantileverPoly()
        D = cant.displacement_resolution()
        assert isinstance(D, float)

    def test_displacement_resolution_positive(self) -> None:
        """Test displacement resolution is positive."""
        cant = CantileverPoly()
        D = cant.displacement_resolution()
        assert D > 0

    def test_displacement_resolution_relation(self) -> None:
        """Test displacement resolution = force resolution / stiffness."""
        cant = CantileverPoly()
        F = cant.force_resolution()
        k = cant.stiffness()
        D = cant.displacement_resolution()

        assert pytest.approx(F / k, rel=1e-10) == D


class TestHoogeAlpha:
    """Test Hooge alpha parameter."""

    def test_hooge_alpha_value(self) -> None:
        """Test Hooge alpha is 1e-3 for polycrystalline."""
        cant = CantileverPoly()
        assert cant.alpha == 1e-3

    def test_hooge_alpha_different_from_single_crystal(self) -> None:
        """Test Hooge alpha differs from single-crystal default (1e-5)."""
        cant = CantileverPoly()
        assert cant.alpha == 1e-3
        # Polycrystalline alpha is 100x higher than single-crystal
        assert cant.alpha > 1e-5


class TestPrintPerformance:
    """Test print_performance method."""

    def test_print_performance_runs(self, capsys: pytest.CaptureFixture) -> None:
        """Test print_performance executes without error."""
        cant = CantileverPoly()
        cant.print_performance()
        captured = capsys.readouterr()
        assert "Cantilever L/W:" in captured.out
        assert "Force resolution:" in captured.out
        assert "Stiffness:" in captured.out


class TestFallbackMaterials:
    """Test fallback behavior for non-standard material configurations."""

    def test_sheet_resistance_fallback(self) -> None:
        """Test sheet resistance fallback for non-standard top material."""
        cant = CantileverPoly(matl_top=Material.SI)  # SI is not POLY, TI, or AL
        Rs = cant.sheet_resistance()
        # Should use poly-Si calculation as fallback
        assert Rs > 0

    def test_piezo_coefficient_fallback(self) -> None:
        """Test piezo coefficient fallback for non-standard top material."""
        cant = CantileverPoly(matl_top=Material.OXIDE)  # OXIDE is not POLY, TI, or AL
        pi = cant.piezo_coefficient()
        # Should use default poly-Si value (103e-11 * 0.6)
        expected = 103e-11 * 0.6
        assert pi == pytest.approx(expected, rel=0.01)


class TestDifferentMaterialConfigurations:
    """Test various material configurations."""

    def test_all_poly_layers(self) -> None:
        """Test cantilever with all poly-Si layers."""
        cant = CantileverPoly(
            matl_top=Material.POLY,
            matl_mid=Material.POLY,
            matl_bot=Material.POLY,
        )
        assert cant.stiffness() > 0
        assert cant.force_resolution() > 0

    def test_oxide_mid_layer(self) -> None:
        """Test cantilever with oxide middle layer."""
        cant = CantileverPoly(
            matl_top=Material.POLY,
            matl_mid=Material.OXIDE,
            matl_bot=Material.POLY,
        )
        assert cant.stiffness() > 0
        assert cant.force_resolution() > 0

    def test_metal_piezoresistor(self) -> None:
        """Test cantilever with metal piezoresistor."""
        cant_ti = CantileverPoly(matl_top=Material.TI)
        cant_al = CantileverPoly(matl_top=Material.AL)

        # Both should work
        assert cant_ti.force_resolution() > 0
        assert cant_al.force_resolution() > 0

        # Metal piezoresistors have lower sensitivity
        cant_poly = CantileverPoly(matl_top=Material.POLY)
        assert cant_ti.force_sensitivity() < cant_poly.force_sensitivity()


class TestPhysicalConstants:
    """Test physical constants are correct."""

    def test_boltzmann_constant(self) -> None:
        """Test Boltzmann constant value."""
        assert CantileverPoly.k_b == pytest.approx(1.38e-23)

    def test_electron_charge(self) -> None:
        """Test electron charge value."""
        assert CantileverPoly.q == pytest.approx(1.60218e-19)

    def test_temperature(self) -> None:
        """Test default temperature."""
        assert CantileverPoly.T == 300


class TestIntegration:
    """Integration tests for full cantilever calculations."""

    def test_full_calculation_workflow(self) -> None:
        """Test complete calculation workflow."""
        cant = CantileverPoly(
            l=200e-6,
            w=20e-6,
            t_top=100e-9,
            t_mid=1000e-9,
            t_bot=100e-9,
            l_pr_ratio=0.4,
            v_bridge=2.0,
            dopant_concentration=5e19,
        )

        # All calculations should complete without error
        na = cant.neutral_axis()
        k = cant.stiffness()
        f = cant.resonant_frequency()
        R = cant.resistance()
        N = cant.number_of_carriers()
        S = cant.force_sensitivity()
        noise = cant.integrated_noise()
        F_res = cant.force_resolution()
        D_res = cant.displacement_resolution()

        # All values should be physical
        assert na > 0
        assert k > 0
        assert f > 0
        assert R > 0
        assert N > 0
        assert S > 0
        assert noise > 0
        assert F_res > 0
        assert D_res > 0

    def test_import_from_package(self) -> None:
        """Test importing from package works."""
        from piezod import CantileverPoly, Material

        cant = CantileverPoly(matl_top=Material.POLY)
        assert cant.stiffness() > 0
