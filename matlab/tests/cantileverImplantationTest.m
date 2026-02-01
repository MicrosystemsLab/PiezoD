classdef cantileverImplantationTest < matlab.unittest.TestCase
    % Tests for cantileverImplantation class

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
            freq_min = 10;
            freq_max = 1000;
            l = 300e-6;
            w = 44e-6;
            t = 89e-9;
            l_pr_ratio = 45/300;
            v_bridge = 5;
            annealing_time = 20*60;
            annealing_temp = 950 + 273;
            implantation_energy = 30;
            implantation_dose = 5e15;
            annealing_type = 'inert';

            testCase.cantilever = cantileverImplantation(freq_min, freq_max, ...
                l, w, t, l_pr_ratio, v_bridge, 'boron', ...
                annealing_time, annealing_temp, annealing_type, ...
                implantation_energy, implantation_dose);
            testCase.cantilever.fluid = 'water';
            testCase.cantilever.number_of_piezoresistors = 4;
        end
    end

    methods (Test)
        function testConstructor(testCase)
            testCase.verifyEqual(testCase.cantilever.implantation_energy, 30);
            testCase.verifyEqual(testCase.cantilever.implantation_dose, 5e15);
            testCase.verifyEqual(testCase.cantilever.annealing_time, 20*60);
            testCase.verifyEqual(testCase.cantilever.annealing_temp, 950+273);
            testCase.verifyEqual(testCase.cantilever.annealing_type, 'inert');
        end

        function testLookupTableLoaded(testCase)
            testCase.verifyNotEmpty(testCase.cantilever.lookupTableData);
        end

        function testStiffness(testCase)
            k = testCase.cantilever.stiffness();
            testCase.verifyEqual(k, 4.85384e-05, 'RelTol', 0.01);
            testCase.verifyGreaterThan(k, 0);
        end

        function testResistance(testCase)
            R = testCase.cantilever.resistance();
            testCase.verifyEqual(R, 763.412, 'RelTol', 0.01);
            testCase.verifyGreaterThan(R, 0);
        end

        function testSheetResistance(testCase)
            Rs = testCase.cantilever.sheet_resistance();
            testCase.verifyEqual(Rs, 69.0004, 'RelTol', 0.01);
            testCase.verifyGreaterThan(Rs, 0);
        end

        function testNz(testCase)
            Nz = testCase.cantilever.Nz();
            testCase.verifyEqual(Nz, 1.86593e19, 'RelTol', 0.01);
            testCase.verifyGreaterThan(Nz, 0);
        end

        function testJunctionDepth(testCase)
            xj = testCase.cantilever.junction_depth();
            testCase.verifyEqual(xj, 1.04889e-06, 'RelTol', 0.01);
            testCase.verifyGreaterThan(xj, 0);
        end

        function testForceSensitivity(testCase)
            % Note: can be negative depending on piezoresistor position
            S = testCase.cantilever.force_sensitivity();
            testCase.verifyEqual(S, -2.47903e06, 'RelTol', 0.01);
        end

        function testForceResolution(testCase)
            % Note: can be negative if sensitivity is negative
            F_min = testCase.cantilever.force_resolution();
            testCase.verifyEqual(F_min, -2.63163e-12, 'RelTol', 0.01);
        end

        function testAlpha(testCase)
            alpha = testCase.cantilever.alpha();
            testCase.verifyEqual(alpha, 5.89717e-07, 'RelTol', 0.01);
            testCase.verifyGreaterThan(alpha, 0);
        end

        function testBeta(testCase)
            beta = testCase.cantilever.beta();
            testCase.verifyEqual(beta, -2.27709, 'RelTol', 0.01);
        end

        function testDopingProfile(testCase)
            [z, active, total] = testCase.cantilever.doping_profile();

            testCase.verifyEqual(length(z), length(active));
            testCase.verifyEqual(length(z), length(total));
            testCase.verifyGreaterThan(min(active), 0);
        end

        function testDiffusionLength(testCase)
            L_d = testCase.cantilever.diffusion_length();
            testCase.verifyGreaterThan(L_d, 0);
        end

        function testAnnealNumberInert(testCase)
            testCase.verifyEqual(testCase.cantilever.annealNumber(), 1);
        end

        function testAnnealNumberOxide(testCase)
            c = cantileverImplantation(10, 1000, 300e-6, 44e-6, 89e-9, ...
                0.15, 5, 'boron', 20*60, 950+273, 'oxide', 30, 5e15);
            testCase.verifyEqual(c.annealNumber(), 2);
        end
    end
end
