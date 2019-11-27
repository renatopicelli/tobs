%% --------------------------------------------------------------------- %%
%            ** Topology Optimization using the BESO method **            %
%-------------------------------------------------------------------------%

% Renato Picelli (University of São Paulo, Santos, Brazil)

% MATLAB code for structural topology optimization using the BESO method
% Find minimum structural compliance subject to a volume constraint

% 03.06.2019

%function opt = stokes2(funobj, sensibility, design_variables, number_of_variables, current_volfrac, iteration)
function opt = stokes2(nvar, x_L, x_U, cst_num, acst_L, acst_U, obj_fun, obj_dfun, cst_fval, jacobian, iteration)

disp('ola Bruno')
disp(nvar)

% close all; clear all; clc


addpath(genpath(strcat('..', filesep, '..', filesep, 'FEA')))
addpath(genpath(strcat('..', filesep, '..', filesep, 'Meshes')))
addpath(genpath(strcat('..', filesep, '..', filesep, 'TopOpt')))
disp(' ')
disp('         *****************************')
disp('         **   Topology Pi - START   **')
disp('         *****************************')
disp(' ')

%% --------------------------------------------------------------------- %%
%                              ** Input **                                %
%-------------------------------------------------------------------------%

% Material properties
%E = 1.0;   % Young's modulus
%nu = 0.3;  % Poisson's ratio

% Optimization parameters
%radius = 6;                  % Filter radius in length unit
%rho_min = 0.001^3;           % Minimum density (for void elements)
compliance_constraint = 190.3072; % Compliance (Nm) constraint
compliance_constraint = 1.0e-10; % Compliance (Nm) constraint
epsilons = 0.001;                 % Constraint relaxation parameter
epsilons = 0.4;                 % Constraint relaxation parameter
flip_limits = 0.05;               % Flip limits
flip_limits = 10.0;               % Flip limits

%% --------------------------------------------------------------------- %%
%                         ** Problem set up **                            %
%-------------------------------------------------------------------------%
% Prepare TOBS
%tobs = TOBS(compliance_constraint, epsilons, flip_limits, size(fea.mesh.incidence,1));
tobs = TOBS(compliance_constraint, epsilons, flip_limits, nvar);


%% --------------------------------------------------------------------- %%
%                           ** Optimization **                            %
%-------------------------------------------------------------------------%

disp(' ')
disp('         *** Optimization loop *** ')
disp(' ')

% Iterations counter
loop = 0;


% Objective (compliance) and sensitivites
%beso.objective = fea.F'*fea.U;
%beso.objective_sensitivities = ComputeComplianceSensitivities(fea,beso.design_variables);
tobs.objective = obj_fun;
disp('******')
disp('aqui')
disp('******')
tobs.objective_sensitivities = obj_dfun;
sensitivities_previous = tobs.objective_sensitivities;
%tobs.history(loop+1,1) = tobs.objective;
%tobs.history(loop+1,2) = tobs.volume_fraction;

% Print optimization status on the screen
disp([' It.: ' sprintf('%3i',iteration) '  Obj.: ' sprintf('%5.4f',full(tobs.objective))...
    '  Vol.: ' ])

% Optimization loop %

% Convergence identifiers
is_converged = 0;
difference = 1;

loop = iteration;

% Solve with BESO design update scheme
%beso = BESODesignUpdate(beso);

% Constraint (compliance) and sensitivities
tobs.constraints = speye (1)*cst_fval;
tobs.constraints_sensitivities = jacobian;

tobs = SolveWithILP(tobs);

% Stabilization technique (average of the sensitivities history)
tobs.constraints_sensitivities = (tobs.constraints_sensitivities+sensitivities_previous)/2;
% Storing sensitivities for next iteration
sensitivities_previous = tobs.constraints_sensitivities;

% Storing optimization history
tobs.history(loop+1,1) = tobs.objective;
tobs.history(loop+1,2) = tobs.constraints;

% Finite Element analysis
opt = tobs.design_variables;

disp(' ')
disp('         ***************************')
disp('         **   Topology Pi - END   **')
disp('         ***************************')
disp(' ')

end
