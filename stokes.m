function opt = stokes(nvar, x_L, x_U, cst_num, acst_L, acst_U, obj_fun, obj_dfun, cst_fval, jacobian, iteration, epsilons, rho)

addpath(genpath('FEA'))
addpath(genpath('Meshes'))
addpath(genpath('TopOpt'))
disp(' ')
disp('         *****************************')
disp('         **   Topology Pi - START   **')
disp('         *****************************')
disp(' ')

%% --------------------------------------------------------------------- %%
%                              ** Input **                                %
%-------------------------------------------------------------------------%

% Optimization parameters
%radius = 6;                  % Filter radius in length unit
%rho_min = 0.001^3;           % Minimum density (for void elements)
volume_constraint = 3.0e-1; % Compliance (Nm) constraint
flip_limits = 0.05;               % Flip limits

%% --------------------------------------------------------------------- %%
%                         ** Problem set up **                            %
%-------------------------------------------------------------------------%
% Prepare TOBS
tobs = TOBS(volume_constraint, epsilons, flip_limits, nvar);


%% --------------------------------------------------------------------- %%
%                           ** Optimization **                            %
%-------------------------------------------------------------------------%

disp(' ')
disp('         *** Optimization loop *** ')
disp(' ')

tobs.design_variables = rho
tobs.objective = obj_fun;
tobs.objective_sensitivities = obj_dfun;

% Print optimization status on the screen
disp([' It.: ' sprintf('%3i',iteration) '  Obj.: ' sprintf('%5.4f',full(tobs.objective))...
    '  Vol.: ' ])

% Optimization loop %

% Convergence identifiers
is_converged = 0;
difference = 1;

loop = iteration;

% Constraint (compliance) and sensitivities
tobs.constraints = speye (1)*cst_fval;
tobs.constraints_sensitivities = jacobian;

[tobs, PythonObjCoeff, ...
		    PythonConstCoeff, PythonRelaxedLimits, ...
		    PythonLowerLimits, PythonUpperLimits, PythonnDesignVariables] = SolveWithILP(tobs);

% Storing optimization history
tobs.history(loop+1,1) = tobs.objective;
tobs.history(loop+1,2) = tobs.constraints;

% Finite Element analysis
opt = {tobs.design_variables, PythonObjCoeff, ...
		    PythonConstCoeff, PythonRelaxedLimits, ...
		    PythonLowerLimits, PythonUpperLimits, PythonnDesignVariables};

end
