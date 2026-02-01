classdef cantileverEpitaxyTest < matlab.unittest.TestCase
    % Tests for cantileverEpitaxy class
    % Reference: Harley 1999 89nm epitaxial cantilever

    properties
        cantilever
    end

    methods (TestClassSetup)
        function addPath(testCase)
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'PiezoD'));
        end
    end

    methods (TestMethodSetup)
        function createCantilever(testCase)
            % Harley 1999 configuration
            freq_min = 10;
            freq_max = 1000;
            l = 300e-6;
            w = 44e-6;
            t = 89e-9;
            l_pr_ratio = 45/300;
            v_bridge = 5;
            doping_type = 'boron';
            concentration = 4e19;
            t_pr_ratio = 30/89;

            testCase.cantilever = cantileverEpitaxy(freq_min, freq_max, ...
                l, w, t, l_pr_ratio, v_bridge, doping_type, concentration, t_pr_ratio);
            testCase.cantilever.fluid = 'vacuum';
            testCase.cantilever.thermal_modeling = 'approx';
            testCase.cantilever.number_of_piezoresistors = 1;
        end
    end

    methods (Test)
        function testConstructor(testCase)
            % Verify object is created with correct properties
            testCase.verifyEqual(testCase.cantilever.l, 300e-6);
            testCase.verifyEqual(testCase.cantilever.w, 44e-6);
            testCase.verifyEqual(testCase.cantilever.t, 89e-9);
            testCase.verifyEqual(testCase.cantilever.doping_type, 'boron');
            testCase.verifyEqual(testCase.cantilever.dopant_concentration, 4e19);
        end

        function testStiffness(testCase)
            k = testCase.cantilever.stiffness();
            testCase.verifyEqual(k, 4.85384e-05, 'RelTol', 0.01);
            testCase.verifyGreaterThan(k, 0);
        end

        function testResistance(testCase)
            R = testCase.cantilever.resistance();
            testCase.verifyEqual(R, 6147.41, 'RelTol', 0.01);
            testCase.verifyGreaterThan(R, 0);
        end

        function testSheetResistance(testCase)
            Rs = testCase.cantilever.sheet_resistance();
            testCase.verifyEqual(Rs, 946.393, 'RelTol', 0.01);
            testCase.verifyGreaterThan(Rs, 0);
        end

        function testNz(testCase)
            Nz = testCase.cantilever.Nz();
            testCase.verifyEqual(Nz, 1.2e18, 'RelTol', 0.01);
            testCase.verifyGreaterThan(Nz, 0);
        end

        function testJunctionDepth(testCase)
            xj = testCase.cantilever.junction_depth();
            testCase.verifyEqual(xj, 3e-08, 'RelTol', 0.01);
            testCase.verifyGreaterThan(xj, 0);
            testCase.verifyLessThan(xj, testCase.cantilever.t);
        end

        function testForceSensitivity(testCase)
            S = testCase.cantilever.force_sensitivity();
            testCase.verifyEqual(S, 800764, 'RelTol', 0.01);
            testCase.verifyGreaterThan(S, 0);
        end

        function testForceResolution(testCase)
            F_min = testCase.cantilever.force_resolution();
            testCase.verifyEqual(F_min, 9.49113e-13, 'RelTol', 0.01);
            testCase.verifyGreaterThan(F_min, 0);
        end

        function testOmegaVacuumHz(testCase)
            f0 = testCase.cantilever.omega_vacuum_hz();
            testCase.verifyEqual(f0, 1359.56, 'RelTol', 0.01);
            testCase.verifyGreaterThan(f0, 0);
        end

        function testIntegratedNoise(testCase)
            V_n = testCase.cantilever.integrated_noise();
            testCase.verifyEqual(V_n, 7.60016e-07, 'RelTol', 0.01);
            testCase.verifyGreaterThan(V_n, 0);
        end

        function testBeta(testCase)
            beta = testCase.cantilever.beta();
            testCase.verifyEqual(beta, 0.315295, 'RelTol', 0.01);
        end

        function testDopingProfile(testCase)
            [z, active, total] = testCase.cantilever.doping_profile();

            % Check array sizes match
            testCase.verifyEqual(length(z), length(active));
            testCase.verifyEqual(length(z), length(total));

            % Check z range
            testCase.verifyEqual(z(1), 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(z(end), testCase.cantilever.t, 'RelTol', 0.01);

            % Check doping is positive
            testCase.verifyGreaterThan(min(active), 0);
            testCase.verifyGreaterThan(min(total), 0);
        end

        function testConsistencyForceResolution(testCase)
            % force_resolution = integrated_noise / force_sensitivity
            F_min = testCase.cantilever.force_resolution();
            V_n = testCase.cantilever.integrated_noise();
            S = testCase.cantilever.force_sensitivity();
            testCase.verifyEqual(F_min, V_n / S, 'RelTol', 0.01);
        end

        function testConsistencyDisplacementResolution(testCase)
            % displacement_resolution = force_resolution / stiffness
            x_min = testCase.cantilever.displacement_resolution();
            F_min = testCase.cantilever.force_resolution();
            k = testCase.cantilever.stiffness();
            testCase.verifyEqual(x_min, F_min / k, 'RelTol', 0.01);
        end
    end
end
