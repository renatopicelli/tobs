%% --------------------------------------------------------------------- %%
%                           ** TOBS class **                              %
%-------------------------------------------------------------------------%

classdef TOBS

    %% Properties
    properties

        % TOBS parameters
        epsilons
        flip_limits
        symmetry

        % Design variables {0,1}
        design_variables

        % Objective and sensitivities
        objective
        objective_sensitivities

        % Constraints and sensitivities
        % [ constraint_1, constraint_2, ..., constraint_n ]
        constraints
        constraints_limits
        constraints_sensitivities

        % Optimization history
        % [ objective, constraint_1, constraint_2, ..., constraint_n ]
        history

	PythonObjCoeff
	PythonConstCoeff
	PythonRelaxedLimits
	PythonLowerLimits
	PythonUpperLimits
	PythonnDesignVariables
	PythonOptimizerOptions

    end

    %% Methods
    methods

        %% Constructor
        function tobs = TOBS(constraints_limits_in, epsilons_in, flip_limits_in, number_of_variables)

            disp([' '])
            disp(['         Preparing TOBS.'])

            % Optimization constraint limits
            tobs.constraints_limits = constraints_limits_in;

            % Optimization parameters
            tobs.epsilons = epsilons_in;
            tobs.flip_limits = flip_limits_in;

            % Initial design variables
            tobs.design_variables = ones(number_of_variables,1);

            % Default symmetry condition
            tobs.symmetry = 0;

        end % end Constructor

        %% Function to solve optimization problem with ILP
            function tobs = SolveWithILP(tobs)

            % Add CPLEX library.
%             addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64');
%             addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\examples\src\matlab');
            % addpath(genpath("/opt/ibm/ILOG/CPLEX_Studio1271/cplex/matlab/x86-64_linux"));
            % addpath(genpath("/opt/ibm/ILOG/CPLEX_Studio1271/cplex/examples/src/matlab"));

            % Prepare ILP problem.
            clear options
            options.Symmetry = tobs.symmetry;
            % options.Optimizer = 'cplex';
            % options.Optimizer = 'intlinprog';
            options.Optimizer = 'glpk';
            % options.Optimizer = 'python_cplex';
            COptimize = ILP (tobs.epsilons', tobs.constraints_limits', tobs.constraints', tobs.design_variables, tobs.flip_limits, 'Minimize');
            tobs.design_variables = COptimize.Optimize (tobs.objective_sensitivities, tobs.constraints_sensitivities, options);

	    PythonObjCoeff = COptimize.return_values_cplex();

	    disp(PythonObjCoeff);

	    tobs.PythonObjCoeff = COptimize.PythonObjCoeff;
	    tobs.PythonConstCoeff = COptimize.PythonConstCoeff;
	    tobs.PythonRelaxedLimits = COptimize.PythonRelaxedLimits;
	    tobs.PythonLowerLimits = COptimize.PythonLowerLimits;
	    tobs.PythonUpperLimits = COptimize.PythonUpperLimits;
	    tobs.PythonnDesignVariables = COptimize.PythonnDesignVariables;


        end % end SolveWithILP

    end % end methods
end
