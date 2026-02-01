"""Tests for CantileverPiezoelectric class."""

import numpy as np
import pytest

from piezod.cantilever_piezoelectric import (
    CantileverPiezoelectric,
    FluidType,
    PiezoMaterial,
)


class TestEnums:
    """Test enum definitions."""

    def test_piezo_material_values(self) -> None:
        """Test PiezoMaterial enum values."""
        assert PiezoMaterial.ALN.value == "AlN"
        assert PiezoMaterial.PZT.value == "PZT"

    def test_fluid_type_values(self) -> None:
        """Test FluidType enum values."""
        assert FluidType.VACUUM.value == "vacuum"
        assert FluidType.AIR.value == "air"
        assert FluidType.WATER.value == "water"


class TestInitialization:
    """Test initialization and basic properties."""

    def test_default_initialization(self) -> None:
        """Test cantilever initializes with default values."""
        cant = CantileverPiezoelectric()

        assert cant.freq_min == 1.0
        assert cant.freq_max == 1e3
        assert cant.l_si == 100e-6
        assert cant.w_si == 10e-6
        assert cant.t_si == 1e-6
        assert cant.t_pe == 500e-9
        assert cant.material == PiezoMaterial.ALN
        assert cant.r_shunt == 1e12
        assert cant.fluid == FluidType.AIR
        assert cant.amplifier_noise is True
        assert cant.thermomechanical_noise is True

    def test_custom_initialization(self) -> None:
        """Test cantilever initializes with custom values."""
        cant = CantileverPiezoelectric(
            freq_min=10.0,
            freq_max=1e4,
            l_si=200e-6,
            w_si=20e-6,
            t_si=2e-6,
            t_pe=1e-6,
            material=PiezoMaterial.PZT,
            r_shunt=1e10,
            fluid=FluidType.WATER,
            amplifier_noise=False,
            thermomechanical_noise=False,
        )

        assert cant.freq_min == 10.0
        assert cant.freq_max == 1e4
        assert cant.l_si == 200e-6
        assert cant.w_si == 20e-6
        assert cant.t_si == 2e-6
        assert cant.t_pe == 1e-6
        assert cant.material == PiezoMaterial.PZT
        assert cant.r_shunt == 1e10
        assert cant.fluid == FluidType.WATER
        assert cant.amplifier_noise is False
        assert cant.thermomechanical_noise is False

    def test_aln_material_properties(self) -> None:
        """Test AlN material properties are set correctly."""
        cant = CantileverPiezoelectric(material=PiezoMaterial.ALN)

        assert cant.permittivity == pytest.approx(10 * 8.854e-12)
        assert cant.resistivity == pytest.approx(1e12)
        assert cant.d31 == pytest.approx(2e-12)
        assert cant.E_pe == pytest.approx(396e9)
        assert cant.rho_pe == pytest.approx(3260)

    def test_pzt_material_properties(self) -> None:
        """Test PZT material properties are set correctly."""
        cant = CantileverPiezoelectric(material=PiezoMaterial.PZT)

        assert cant.permittivity == pytest.approx(900 * 8.854e-12)
        assert cant.resistivity == pytest.approx(1e8)
        assert cant.d31 == pytest.approx(70e-12)
        assert cant.E_pe == pytest.approx(55e9)
        assert cant.rho_pe == pytest.approx(7550)

    def test_l_pe_w_pe_equal_si(self) -> None:
        """Test piezoelectric dimensions match silicon."""
        cant = CantileverPiezoelectric(l_si=150e-6, w_si=15e-6)

        assert cant.l_pe == cant.l_si
        assert cant.w_pe == cant.w_si


class TestBeamMechanics:
    """Test beam mechanics calculations."""

    def test_neutral_axis_positive(self) -> None:
        """Test neutral axis is positive."""
        cant = CantileverPiezoelectric()
        assert cant.neutral_axis > 0

    def test_neutral_axis_is_property(self) -> None:
        """Test neutral_axis is a property."""
        cant = CantileverPiezoelectric()
        assert isinstance(cant.neutral_axis, float)

    def test_normalized_curvature_positive(self) -> None:
        """Test normalized curvature is positive."""
        cant = CantileverPiezoelectric()
        assert cant.normalized_curvature > 0

    def test_stiffness_positive(self) -> None:
        """Test stiffness is positive."""
        cant = CantileverPiezoelectric()
        assert cant.stiffness > 0

    def test_stiffness_increases_with_thickness(self) -> None:
        """Test stiffness increases with silicon thickness."""
        cant1 = CantileverPiezoelectric(t_si=1e-6)
        cant2 = CantileverPiezoelectric(t_si=2e-6)

        assert cant2.stiffness > cant1.stiffness

    def test_stiffness_decreases_with_length(self) -> None:
        """Test stiffness decreases with length (cubic)."""
        cant1 = CantileverPiezoelectric(l_si=100e-6)
        cant2 = CantileverPiezoelectric(l_si=200e-6)

        assert cant1.stiffness > cant2.stiffness

    def test_effective_mass_positive(self) -> None:
        """Test effective mass is positive."""
        cant = CantileverPiezoelectric()
        assert cant.effective_mass > 0

    def test_resonant_frequency_positive(self) -> None:
        """Test resonant frequency is positive."""
        cant = CantileverPiezoelectric()
        assert cant.resonant_frequency > 0

    def test_resonant_frequency_increases_with_stiffness(self) -> None:
        """Test resonant frequency increases with stiffness."""
        cant1 = CantileverPiezoelectric(t_si=1e-6)
        cant2 = CantileverPiezoelectric(t_si=2e-6)  # stiffer

        assert cant2.resonant_frequency > cant1.resonant_frequency


class TestRCCalculations:
    """Test resistance and capacitance calculations."""

    def test_Rpe_positive(self) -> None:
        """Test piezo resistance is positive."""
        cant = CantileverPiezoelectric()
        assert cant.Rpe > 0

    def test_Rpe_higher_for_aln(self) -> None:
        """Test AlN has higher resistance than PZT."""
        cant_aln = CantileverPiezoelectric(material=PiezoMaterial.ALN)
        cant_pzt = CantileverPiezoelectric(material=PiezoMaterial.PZT)

        assert cant_aln.Rpe > cant_pzt.Rpe

    def test_Cpe_positive(self) -> None:
        """Test piezo capacitance is positive."""
        cant = CantileverPiezoelectric()
        assert cant.Cpe > 0

    def test_Cpe_higher_for_pzt(self) -> None:
        """Test PZT has higher capacitance than AlN."""
        cant_aln = CantileverPiezoelectric(material=PiezoMaterial.ALN)
        cant_pzt = CantileverPiezoelectric(material=PiezoMaterial.PZT)

        assert cant_pzt.Cpe > cant_aln.Cpe

    def test_R_half_positive(self) -> None:
        """Test half-circuit resistance is positive."""
        cant = CantileverPiezoelectric()
        assert cant.R_half > 0

    def test_R_half_less_than_Rpe(self) -> None:
        """Test half-circuit resistance is less than Rpe."""
        cant = CantileverPiezoelectric()
        assert cant.R_half < cant.Rpe

    def test_C_half_positive(self) -> None:
        """Test half-circuit capacitance is positive."""
        cant = CantileverPiezoelectric()
        assert cant.C_half > 0

    def test_Z_returns_complex_array(self) -> None:
        """Test Z returns complex array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        Z = cant.Z(freq)

        assert isinstance(Z, np.ndarray)
        assert Z.dtype == np.complex128
        assert len(Z) == len(freq)

    def test_Z_half_returns_complex_array(self) -> None:
        """Test Z_half returns complex array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        Z = cant.Z_half(freq)

        assert isinstance(Z, np.ndarray)
        assert Z.dtype == np.complex128
        assert len(Z) == len(freq)


class TestSensitivity:
    """Test sensitivity calculations."""

    def test_q_f_sensitivity_intrinsic_positive(self) -> None:
        """Test intrinsic charge sensitivity is positive."""
        cant = CantileverPiezoelectric()
        assert cant.q_f_sensitivity_intrinsic > 0

    def test_q_f_sensitivity_returns_array(self) -> None:
        """Test q_f_sensitivity returns array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        sens = cant.q_f_sensitivity(freq)

        assert isinstance(sens, np.ndarray)
        assert len(sens) == len(freq)
        assert np.all(sens > 0)

    def test_q_f_sensitivity_constant_with_freq(self) -> None:
        """Test charge sensitivity is constant with frequency."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0, 1000.0])
        sens = cant.q_f_sensitivity(freq)

        assert np.allclose(sens, sens[0])

    def test_q_x_sensitivity_returns_array(self) -> None:
        """Test q_x_sensitivity returns array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        sens = cant.q_x_sensitivity(freq)

        assert isinstance(sens, np.ndarray)
        assert len(sens) == len(freq)
        assert np.all(sens > 0)

    def test_v_f_sensitivity_returns_array(self) -> None:
        """Test v_f_sensitivity returns array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        sens = cant.v_f_sensitivity(freq)

        assert isinstance(sens, np.ndarray)
        assert len(sens) == len(freq)
        assert np.all(sens > 0)

    def test_v_f_sensitivity_increases_with_freq(self) -> None:
        """Test voltage sensitivity increases with frequency at low freq."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        sens = cant.v_f_sensitivity(freq)

        assert sens[1] > sens[0]
        assert sens[2] > sens[1]

    def test_v_x_sensitivity_returns_array(self) -> None:
        """Test v_x_sensitivity returns array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        sens = cant.v_x_sensitivity(freq)

        assert isinstance(sens, np.ndarray)
        assert len(sens) == len(freq)
        assert np.all(sens > 0)


class TestDamping:
    """Test damping and Q factor calculations."""

    def test_kappa_positive(self) -> None:
        """Test kappa is positive."""
        cant = CantileverPiezoelectric()
        assert cant.kappa > 0

    def test_omega_damped_and_Q_vacuum(self) -> None:
        """Test vacuum returns max Q and vacuum frequency."""
        cant = CantileverPiezoelectric(fluid=FluidType.VACUUM)
        omega_d, Q = cant.omega_damped_and_Q()

        assert Q == cant.MAX_Q
        assert omega_d == pytest.approx(2 * np.pi * cant.resonant_frequency)

    def test_omega_damped_and_Q_air(self) -> None:
        """Test air returns damped values."""
        cant = CantileverPiezoelectric(fluid=FluidType.AIR)
        omega_d, Q = cant.omega_damped_and_Q()

        omega_vac = 2 * np.pi * cant.resonant_frequency
        assert omega_d > 0
        assert omega_d <= omega_vac
        assert Q > cant.MIN_Q
        assert Q <= cant.MAX_Q

    def test_omega_damped_and_Q_water(self) -> None:
        """Test water returns lower Q than air."""
        cant_air = CantileverPiezoelectric(fluid=FluidType.AIR)
        cant_water = CantileverPiezoelectric(fluid=FluidType.WATER)

        _, Q_air = cant_air.omega_damped_and_Q()
        _, Q_water = cant_water.omega_damped_and_Q()

        assert Q_water < Q_air

    def test_lookup_tables_shape(self) -> None:
        """Test lookup tables have consistent shapes."""
        cant = CantileverPiezoelectric()

        assert cant.TAU_LOOKUP_REAL.shape == (17, 12)
        assert cant.TAU_LOOKUP_IMAG.shape == (17, 12)
        assert len(cant.REYNOLDS_LOOKUP) == 17
        assert len(cant.KAPPA_LOOKUP) == 12


class TestNoise:
    """Test noise calculations."""

    def test_Vth_returns_array(self) -> None:
        """Test thermomechanical voltage noise returns array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        Vth = cant.Vth(freq)

        assert isinstance(Vth, np.ndarray)
        assert len(Vth) == len(freq)
        assert np.all(Vth > 0)

    def test_Vamp_returns_array(self) -> None:
        """Test amplifier voltage noise returns array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        Vamp = cant.Vamp(freq)

        assert isinstance(Vamp, np.ndarray)
        assert len(Vamp) == len(freq)
        assert np.all(Vamp > 0)

    def test_Vn_returns_array(self) -> None:
        """Test total voltage noise returns array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        Vn = cant.Vn(freq)

        assert isinstance(Vn, np.ndarray)
        assert len(Vn) == len(freq)
        assert np.all(Vn > 0)

    def test_Vn_no_amplifier_noise(self) -> None:
        """Test Vn with amplifier noise disabled."""
        cant1 = CantileverPiezoelectric(amplifier_noise=True)
        cant2 = CantileverPiezoelectric(amplifier_noise=False)
        freq = np.array([10.0, 100.0])

        Vn1 = cant1.Vn(freq)
        Vn2 = cant2.Vn(freq)

        assert np.all(Vn2 < Vn1)

    def test_Vn_no_thermomechanical_noise(self) -> None:
        """Test Vn with thermomechanical noise disabled."""
        cant1 = CantileverPiezoelectric(thermomechanical_noise=True)
        cant2 = CantileverPiezoelectric(thermomechanical_noise=False)
        freq = np.array([10.0, 100.0])

        Vn1 = cant1.Vn(freq)
        Vn2 = cant2.Vn(freq)

        assert np.all(Vn2 < Vn1)

    def test_Qamp_returns_array(self) -> None:
        """Test amplifier charge noise returns array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        Qamp = cant.Qamp(freq)

        assert isinstance(Qamp, np.ndarray)
        assert len(Qamp) == len(freq)
        assert np.all(Qamp > 0)

    def test_Qn_returns_array(self) -> None:
        """Test total charge noise returns array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        Qn = cant.Qn(freq)

        assert isinstance(Qn, np.ndarray)
        assert len(Qn) == len(freq)
        assert np.all(Qn > 0)


class TestResolution:
    """Test resolution calculations."""

    def test_freq_for_plot(self) -> None:
        """Test frequency array generation."""
        cant = CantileverPiezoelectric(freq_min=1.0, freq_max=1000.0)
        freq = cant.freq_for_plot()

        assert len(freq) == cant.N_PLOT_POINTS
        assert freq[0] == pytest.approx(1.0)
        assert freq[-1] == pytest.approx(1000.0)

    def test_Sfminv_returns_array(self) -> None:
        """Test force noise density in voltage mode returns array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        Sf = cant.Sfminv(freq)

        assert isinstance(Sf, np.ndarray)
        assert len(Sf) == len(freq)
        assert np.all(Sf > 0)

    def test_Sxminv_returns_array(self) -> None:
        """Test displacement noise density in voltage mode returns array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        Sx = cant.Sxminv(freq)

        assert isinstance(Sx, np.ndarray)
        assert len(Sx) == len(freq)
        assert np.all(Sx > 0)

    def test_Sfminq_returns_array(self) -> None:
        """Test force noise density in charge mode returns array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        Sf = cant.Sfminq(freq)

        assert isinstance(Sf, np.ndarray)
        assert len(Sf) == len(freq)
        assert np.all(Sf > 0)

    def test_Sxminq_returns_array(self) -> None:
        """Test displacement noise density in charge mode returns array."""
        cant = CantileverPiezoelectric()
        freq = np.array([1.0, 10.0, 100.0])
        Sx = cant.Sxminq(freq)

        assert isinstance(Sx, np.ndarray)
        assert len(Sx) == len(freq)
        assert np.all(Sx > 0)

    def test_Fminv_positive(self) -> None:
        """Test minimum detectable force in voltage mode is positive."""
        cant = CantileverPiezoelectric()
        Fmin = cant.Fminv()

        assert isinstance(Fmin, float)
        assert Fmin > 0

    def test_Xminv_positive(self) -> None:
        """Test minimum detectable displacement in voltage mode is positive."""
        cant = CantileverPiezoelectric()
        Xmin = cant.Xminv()

        assert isinstance(Xmin, float)
        assert Xmin > 0

    def test_Fminq_positive(self) -> None:
        """Test minimum detectable force in charge mode is positive."""
        cant = CantileverPiezoelectric()
        Fmin = cant.Fminq()

        assert isinstance(Fmin, float)
        assert Fmin > 0

    def test_Xminq_positive(self) -> None:
        """Test minimum detectable displacement in charge mode is positive."""
        cant = CantileverPiezoelectric()
        Xmin = cant.Xminq()

        assert isinstance(Xmin, float)
        assert Xmin > 0

    def test_Fminv_relates_to_Xminv(self) -> None:
        """Test Fminv and Xminv are related through stiffness."""
        cant = CantileverPiezoelectric()
        Fmin = cant.Fminv()
        Xmin = cant.Xminv()

        # Xmin = Fmin / k approximately
        ratio = Fmin / cant.stiffness
        assert ratio == pytest.approx(Xmin, rel=0.1)

    def test_Fminq_relates_to_Xminq(self) -> None:
        """Test Fminq and Xminq are related through stiffness."""
        cant = CantileverPiezoelectric()
        Fmin = cant.Fminq()
        Xmin = cant.Xminq()

        # Xmin = Fmin / k approximately
        ratio = Fmin / cant.stiffness
        assert ratio == pytest.approx(Xmin, rel=0.1)

    def test_thermomechanical_force_limit_positive(self) -> None:
        """Test thermomechanical force limit is positive."""
        cant = CantileverPiezoelectric()
        limit = cant.thermomechanical_force_limit()

        assert isinstance(limit, float)
        assert limit > 0

    def test_Vn_total_positive(self) -> None:
        """Test total integrated voltage noise is positive."""
        cant = CantileverPiezoelectric()
        Vn = cant.Vn_total()

        assert isinstance(Vn, float)
        assert Vn > 0

    def test_Qn_total_positive(self) -> None:
        """Test total integrated charge noise is positive."""
        cant = CantileverPiezoelectric()
        Qn = cant.Qn_total()

        assert isinstance(Qn, float)
        assert Qn > 0


class TestPhysicalConstants:
    """Test physical constants are correctly defined."""

    def test_epsilon_0(self) -> None:
        """Test vacuum permittivity."""
        assert pytest.approx(8.854e-12) == CantileverPiezoelectric.EPSILON_0

    def test_k_b(self) -> None:
        """Test Boltzmann constant."""
        assert pytest.approx(1.38e-23) == CantileverPiezoelectric.K_B

    def test_temperature(self) -> None:
        """Test temperature is 300 K."""
        assert CantileverPiezoelectric.T == 300

    def test_fluid_properties(self) -> None:
        """Test fluid properties."""
        assert pytest.approx(1e3) == CantileverPiezoelectric.RHO_WATER
        assert pytest.approx(0.9e-3) == CantileverPiezoelectric.ETA_WATER
        assert pytest.approx(1.2) == CantileverPiezoelectric.RHO_AIR
        assert pytest.approx(17e-6) == CantileverPiezoelectric.ETA_AIR


class TestMaterialComparison:
    """Test comparisons between AlN and PZT materials."""

    def test_pzt_higher_charge_sensitivity(self) -> None:
        """Test PZT has higher charge sensitivity due to higher d31."""
        cant_aln = CantileverPiezoelectric(material=PiezoMaterial.ALN)
        cant_pzt = CantileverPiezoelectric(material=PiezoMaterial.PZT)

        # PZT d31 = 70e-12, AlN d31 = 2e-12
        # However, sensitivity also depends on modulus and geometry
        # PZT has lower modulus, which affects the result
        freq = np.array([100.0])
        sens_aln = cant_aln.q_f_sensitivity(freq)[0]
        sens_pzt = cant_pzt.q_f_sensitivity(freq)[0]

        # PZT should have higher charge sensitivity
        assert sens_pzt > sens_aln


class TestFluidComparison:
    """Test comparisons between fluid environments."""

    def test_vacuum_highest_Q(self) -> None:
        """Test vacuum has highest Q factor."""
        cant_vac = CantileverPiezoelectric(fluid=FluidType.VACUUM)
        cant_air = CantileverPiezoelectric(fluid=FluidType.AIR)
        cant_water = CantileverPiezoelectric(fluid=FluidType.WATER)

        _, Q_vac = cant_vac.omega_damped_and_Q()
        _, Q_air = cant_air.omega_damped_and_Q()
        _, Q_water = cant_water.omega_damped_and_Q()

        assert Q_vac > Q_air > Q_water

    def test_water_lowest_damped_frequency(self) -> None:
        """Test water has lowest damped frequency."""
        cant_air = CantileverPiezoelectric(fluid=FluidType.AIR)
        cant_water = CantileverPiezoelectric(fluid=FluidType.WATER)

        omega_air, _ = cant_air.omega_damped_and_Q()
        omega_water, _ = cant_water.omega_damped_and_Q()

        assert omega_water < omega_air


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_very_small_cantilever(self) -> None:
        """Test with very small cantilever dimensions."""
        cant = CantileverPiezoelectric(
            l_si=10e-6,
            w_si=1e-6,
            t_si=100e-9,
            t_pe=50e-9,
        )

        assert cant.stiffness > 0
        assert cant.resonant_frequency > 0
        assert cant.Fminv() > 0

    def test_very_large_cantilever(self) -> None:
        """Test with large cantilever dimensions."""
        cant = CantileverPiezoelectric(
            l_si=1e-3,
            w_si=100e-6,
            t_si=10e-6,
            t_pe=2e-6,
        )

        assert cant.stiffness > 0
        assert cant.resonant_frequency > 0
        assert cant.Fminv() > 0

    def test_high_frequency_range(self) -> None:
        """Test with high frequency range."""
        cant = CantileverPiezoelectric(freq_min=1e3, freq_max=1e6)

        assert cant.Fminv() > 0
        assert cant.Fminq() > 0

    def test_low_shunt_resistance(self) -> None:
        """Test with low shunt resistance."""
        cant = CantileverPiezoelectric(r_shunt=1e6)

        assert cant.R_half > 0
        assert cant.R_half < cant.Rpe
        assert cant.Fminv() > 0
