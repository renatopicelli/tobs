%% --------------------------------------------------------------------- %%
%         ** Topology Optimization of Binary Structures (TOBS) **         %
%-------------------------------------------------------------------------%

% Renato Picelli (University of São Paulo, Santos, Brazil)
% Raghavendra Sivapuram (University of California, San Diego, USA)

% MATLAB code for structural topology optimization using the TOBS method
% Find minimum structural volume subject to a compliance constraint
% Hydrostatic pressure loads solved via Laplace's equation

% 04.06.2019

function Hydrostatic_load_arch

close all; clear all; clc

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
E = 1;    % Young's modulus
nu = 0.3; % Poisson's ratio

% Optimization parameters
radius = 0.03;               % Filter radius in length unit
compliance_constraint = 7.0; % Final volume fraction
epsilons = 0.001;            % Constraint relaxation parameter
flip_limits = 0.04;          % Flip limits
tau = 1e-4;                  % Convergence tolerance

%% --------------------------------------------------------------------- %%
%                         ** Problem set up **                            %
%-------------------------------------------------------------------------%

% Set up mesh.
mesh = Mesh_FSI('Arch');

% Prepare FEA
fea = FluidStructureInteraction(mesh);
fea = AddSolidMaterial(fea, E, nu, 1.0);
fea.H = BuildFilterMatrix(fea,radius);

% Prepare TOBS
tobs = TOBS(compliance_constraint, epsilons, flip_limits, length(fea.design_domain));

%% --------------------------------------------------------------------- %%
%                           ** Optimization **                            %
%-------------------------------------------------------------------------%

disp(' ')
disp('         *** Optimization loop *** ')
disp(' ')

% Iteration 0 %
loop = 0;

% Finite Element analysis
% Assemble F
fea = AssemblePointLoads(fea);
% Assemble K
fea = AssembleCoupledStructuralAndHydrostaticK(fea);
% Solve coupled structural-hydrostatic analysis
fea = SolveFEA(fea);

% Objective and sensitivites
tobs.objective = mean(tobs.design_variables);
tobs.objective_sensitivities = ones(length(fea.design_domain),1)/length(fea.design_domain);

% Constraints and sensitivities
[tobs.constraints_sensitivities,tobs.constraints] = ComputeComplianceWithCoupledPressureSensitivities(fea);
% Filtering sensitivities
tobs.constraints_sensitivities = fea.H*tobs.constraints_sensitivities;

% Storing sensitivities for the next iteration averaging
sensitivities_previous = tobs.constraints_sensitivities;

% Storing optimization history
tobs.history(loop+1,1) = tobs.objective;
tobs.history(loop+1,2) = tobs.constraints;

% Print optimization status on the screen
disp([' It.: ' sprintf('%3i',loop) '  Obj.: ' sprintf('%5.4f',full(tobs.objective))...
    '  Comp.: ' sprintf('%3.3f',full(tobs.constraints))])

% Optimization loop %

% Convergence identifier
is_converged = 0;
difference = 1;

while (is_converged == 0)

    % Counter update
    loop = loop + 1;

    % Solve with ILP
    tobs = SolveWithILP(tobs);

    % Fluid flooding process (update element types)
    fea = FluidFlooding(fea,tobs.design_variables);

    % Finite Element analysis
    % Assemble F
    fea = AssemblePointLoads(fea);
    % Assemble K
    fea = AssembleCoupledStructuralAndHydrostaticK(fea);
    % Solve coupled structural-hydrostatic analysis
    fea = SolveFEA(fea);

    % Objective and sensitivites
    tobs.objective = mean(tobs.design_variables);
    tobs.objective_sensitivities = ones(length(fea.design_domain),1)/length(fea.design_domain);

    % Constraints and sensitivities
    [tobs.constraints_sensitivities,tobs.constraints] = ComputeComplianceWithCoupledPressureSensitivities(fea);
    % Filtering sensitivities
    tobs.constraints_sensitivities = fea.H*tobs.constraints_sensitivities;

    % Stabilization technique (average of the sensitivities history)
    tobs.constraints_sensitivities = (tobs.constraints_sensitivities+sensitivities_previous)/2;
    % Storing sensitivities for next iteration
    sensitivities_previous = tobs.constraints_sensitivities;

    % Storing optimization history
    tobs.history(loop+1,1) = tobs.objective;
    tobs.history(loop+1,2) = tobs.constraints;

    % Convergence analysis [Huang and Xie, 2007]
    N = 5;
    if (loop >= 2*N) % Analyzing 2*N consecutive iterations
        error_1 = zeros(N,1);
        error_2 = zeros(N,1);
        for i = 1:N
            error_1(i) = tobs.history(loop-i+1,1)-tobs.history(loop-N-i+1,1);
            error_2(i) = tobs.history(loop-i+1,1);
        end
        % Evaluating convergence
        difference = abs(sum(error_1))/sum(error_2);
        % Verifying error tolerance and if constraint is satisfied
        if ((difference <= tau) && (tobs.history(loop,2) <= 1.001*tobs.constraints_limits))
            is_converged = 1;
        end
    else
        difference = 1;
    end

    % Print results
    disp([' It.: ' sprintf('%3i',loop) '  Obj.: ' sprintf('%5.4f',full(tobs.objective))...
        '  Comp.: ' sprintf('%3.3f',tobs.constraints)...
        '  Conv.: ' sprintf('%4.4f',difference)])

    % Plot densities
    PlotFluidStructureTopology(fea); pause(1e-6);

    % Stop at maximum iterations
    if (loop == 300)
        break
    end

end

% Finite Element analysis of the final solution
fea = AssemblePointLoads(fea);
fea = AssembleCoupledStructuralAndHydrostaticK(fea);
fea = SolveFEA(fea);

% Storing objective and constraint
tobs.history(loop+1,1) = mean(tobs.design_variables);
[~,tobs.history(loop+1,2)] = ComputeComplianceWithCoupledPressureSensitivities(fea);

% Plot optimization history
PlotVolumeHistory(tobs.history(:,1))
PlotComplianceHistory(tobs.history(:,2))

disp(' ')
disp('         ***************************')
disp('         **   Topology Pi - END   **')
disp('         ***************************')
disp(' ')

end
