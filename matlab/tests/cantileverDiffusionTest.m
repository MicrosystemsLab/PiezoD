classdef cantileverDiffusionTest < matlab.unittest.TestCase
    % Tests for cantileverDiffusion class

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
            diffusion_time = 20*60; % 20 minutes
            diffusion_temp = 800 + 273; % 800C

            testCase.cantilever = cantileverDiffusion(freq_min, freq_max, ...
                l, w, t, l_pr_ratio, v_bridge, 'phosphorus', ...
                diffusion_time, diffusion_temp);
            testCase.cantilever.fluid = 'vacuum';
            testCase.cantilever.number_of_piezoresistors = 4;
        end
    end

    methods (Test)
        function testConstructor(testCase)
            testCase.verifyEqual(testCase.cantilever.diffusion_time, 20*60);
            testCase.verifyEqual(testCase.cantilever.diffusion_temp, 800+273);
            testCase.verifyEqual(testCase.cantilever.doping_type, 'phosphorus');
        end

        function testStiffness(testCase)
            k = testCase.cantilever.stiffness();
            testCase.verifyEqual(k, 3.73372e-05, 'RelTol', 0.01);
            testCase.verifyGreaterThan(k, 0);
        end

        function testResistance(testCase)
            R = testCase.cantilever.resistance();
            testCase.verifyEqual(R, 1442.05, 'RelTol', 0.01);
            testCase.verifyGreaterThan(R, 0);
        end

        function testSheetResistance(testCase)
            Rs = testCase.cantilever.sheet_resistance();
            testCase.verifyEqual(Rs, 179.593, 'RelTol', 0.01);
            testCase.verifyGreaterThan(Rs, 0);
        end

        function testNz(testCase)
            Nz = testCase.cantilever.Nz();
            testCase.verifyEqual(Nz, 4.75808e18, 'RelTol', 0.01);
            testCase.verifyGreaterThan(Nz, 0);
        end

        function testForceSensitivity(testCase)
            S = testCase.cantilever.force_sensitivity();
            testCase.verifyEqual(S, 455495, 'RelTol', 0.01);
            testCase.verifyGreaterThan(S, 0);
        end

        function testForceResolution(testCase)
            F_min = testCase.cantilever.force_resolution();
            testCase.verifyEqual(F_min, 8.34544e-13, 'RelTol', 0.01);
            testCase.verifyGreaterThan(F_min, 0);
        end

        function testDopingProfilePhosphorus(testCase)
            [z, active, total] = testCase.cantilever.doping_profile();

            % Check array sizes
            testCase.verifyEqual(length(z), length(active));
            testCase.verifyEqual(length(z), length(total));

            % Phosphorus profile should have surface concentration
            testCase.verifyGreaterThan(active(1), 1e19);

            % Check z range
            testCase.verifyEqual(z(1), 0, 'AbsTol', 1e-12);
        end

        function testDopingProfileBoron(testCase)
            c = cantileverDiffusion(10, 1000, 300e-6, 44e-6, 500e-9, ...
                0.15, 5, 'boron', 20*60, 900+273);
            [z, active, total] = c.doping_profile();

            % Boron erfc profile
            testCase.verifyGreaterThan(active(1), 1e19);
            % Should decrease with depth
            testCase.verifyGreaterThan(active(1), active(end));
        end

        function testDopingProfileArsenic(testCase)
            c = cantileverDiffusion(10, 1000, 300e-6, 44e-6, 500e-9, ...
                0.15, 5, 'arsenic', 20*60, 900+273);
            [z, active, total] = c.doping_profile();

            % Arsenic erfc profile
            testCase.verifyGreaterThan(active(1), 1e19);
            testCase.verifyGreaterThan(active(1), active(end));
        end

        function testAlphaUsesDefault(testCase)
            alpha = testCase.cantilever.alpha();
            testCase.verifyEqual(alpha, testCase.cantilever.default_alpha);
        end
    end
end
