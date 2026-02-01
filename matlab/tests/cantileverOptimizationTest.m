classdef cantileverOptimizationTest < matlab.unittest.TestCase
    % Tests for cantilever optimization routines
    % Uses 5% tolerance since optimization is inherently approximate

    methods (TestClassSetup)
        function addPath(testCase)
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'PiezoD'));
        end
    end

    methods (Test)
        function testEpitaxyOptimizationFromCurrent(testCase)
            % Create starting cantilever
            c = cantileverEpitaxy(10, 1000, 300e-6, 44e-6, 1e-6, ...
                0.15, 5, 'boron', 4e19, 0.3);
            c.fluid = 'vacuum';
            c.number_of_piezoresistors = 2;

            initial_resolution = c.force_resolution();

            % Setup constraints
            parameter_constraints = {{'min_t', 'max_t', 'min_w', 'max_v_bridge'}, ...
                {1e-6, 1e-6, 5e-6, 10}};
            nonlinear_constraints = {{'omega_min_hz', 'max_power', 'min_k', 'max_k'}, ...
                {1000, 2e-3, 1e-3, 1e1}};
            goal = cantilever.goalForceResolution;

            % Optimize
            c_opt = c.optimize_performance_from_current(parameter_constraints, ...
                nonlinear_constraints, goal);

            % Verify optimization improved or maintained resolution
            final_resolution = c_opt.force_resolution();
            testCase.verifyLessThanOrEqual(final_resolution, initial_resolution * 1.05);

            % Verify constraints are satisfied
            testCase.verifyGreaterThanOrEqual(c_opt.omega_vacuum_hz(), 1000 * 0.99);
            testCase.verifyLessThanOrEqual(c_opt.power_dissipation(), 2e-3 * 1.01);
            testCase.verifyGreaterThanOrEqual(c_opt.stiffness(), 1e-3 * 0.99);
            testCase.verifyLessThanOrEqual(c_opt.stiffness(), 1e1 * 1.01);
        end

        function testDiffusionOptimizationFromCurrent(testCase)
            c = cantileverDiffusion(10, 1000, 300e-6, 44e-6, 1e-6, ...
                0.15, 5, 'phosphorus', 30*60, 850+273);
            c.fluid = 'vacuum';
            c.number_of_piezoresistors = 2;

            initial_resolution = c.force_resolution();

            parameter_constraints = {{'min_t', 'max_t'}, {1e-6, 1e-6}};
            nonlinear_constraints = {{'omega_min_hz', 'max_power'}, {500, 5e-3}};
            goal = cantilever.goalForceResolution;

            c_opt = c.optimize_performance_from_current(parameter_constraints, ...
                nonlinear_constraints, goal);

            final_resolution = c_opt.force_resolution();
            testCase.verifyLessThanOrEqual(final_resolution, initial_resolution * 1.05);
        end

        function testOptimizationBoundsEpitaxy(testCase)
            c = cantileverEpitaxy(10, 1000, 300e-6, 44e-6, 1e-6, ...
                0.15, 5, 'boron', 4e19, 0.3);

            % Test with no constraints
            [lb, ub] = c.optimization_bounds({});
            testCase.verifyEqual(length(lb), length(ub));
            testCase.verifyTrue(all(lb < ub));

            % Test with custom constraints
            constraints = {{'min_l', 'max_l'}, {100e-6, 500e-6}};
            [lb, ub] = c.optimization_bounds(constraints);
            testCase.verifyEqual(lb(1), 100e-6);
            testCase.verifyEqual(ub(1), 500e-6);
        end

        function testOptimizationBoundsDiffusion(testCase)
            c = cantileverDiffusion(10, 1000, 300e-6, 44e-6, 1e-6, ...
                0.15, 5, 'phosphorus', 30*60, 850+273);

            [lb, ub] = c.doping_optimization_bounds({});
            testCase.verifyEqual(length(lb), 2); % time, temp
            testCase.verifyEqual(length(ub), 2);
            testCase.verifyTrue(all(lb < ub));
        end

        function testOptimizationBoundsImplantation(testCase)
            c = cantileverImplantation(10, 1000, 300e-6, 44e-6, 1e-6, ...
                0.15, 5, 'boron', 30*60, 950+273, 'inert', 30, 5e15);

            [lb, ub] = c.doping_optimization_bounds({});
            testCase.verifyEqual(length(lb), 4); % time, temp, energy, dose
            testCase.verifyEqual(length(ub), 4);
            testCase.verifyTrue(all(lb < ub));
        end

        function testRandomInitialConditions(testCase)
            c = cantileverEpitaxy(10, 1000, 300e-6, 44e-6, 1e-6, ...
                0.15, 5, 'boron', 4e19, 0.3);

            [lb, ub] = c.optimization_bounds({});

            % Generate multiple random conditions and verify bounds
            for i = 1:5
                x0 = c.initial_conditions_random({});
                testCase.verifyTrue(all(x0 >= lb));
                testCase.verifyTrue(all(x0 <= ub));
            end
        end

        function testOptimizationScaling(testCase)
            c = cantileverEpitaxy(10, 1000, 300e-6, 44e-6, 1e-6, ...
                0.15, 5, 'boron', 4e19, 0.3);

            scaling = c.optimization_scaling();
            testCase.verifyGreaterThan(length(scaling), 0);
            testCase.verifyTrue(all(scaling > 0));
        end
    end
end
