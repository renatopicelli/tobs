%% --------------------------------------------------------------------- %%
%         ** Topology Optimization of Binary Structures (TOBS) **         %
%-------------------------------------------------------------------------%

% Renato Picelli (University of São Paulo, Santos, Brazil)
% Raghavendra Sivapuram (University of California, San Diego, USA)

% MATLAB code for structural topology optimization using the TOBS method
% Find minimum structural volume subject to a compliance constraint

% 03.06.2019

function MBB_120x40_minvol_meshgen

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
E = 1.0;   % Young's modulus
nu = 0.3;  % Poisson's ratio

% Optimization parameters
radius = 6;                       % Filter radius in length unit
rho_min = 0.001^3;                % Minimum density (for void elements)
compliance_constraint = 190.3072; % Compliance (Nm) constraint
epsilons = 0.001;                 % Constraint relaxation parameter
flip_limits = 0.05;               % Flip limits
tau = 0.001;                      % Convergence tolerance

%% --------------------------------------------------------------------- %%
%                         ** Problem set up **                            %
%-------------------------------------------------------------------------%

% Set up mesh. Mesh() means no input mesh will be read.
mesh = Mesh();

% Set up discretization size.
mesh.Lx = 120.0;
mesh.Ly = 40.0;
mesh.nelx = 120;
mesh.nely = 40;
% Generate regular quadrilateral grid.
mesh = GenerateQuad4RectangleMesh(mesh);

% Add Dirichlet boundary conditions
dirichlet.point_1 = [-0.001, -0.001];
dirichlet.point_2 = [0.001, mesh.Ly*1.001];
dirichlet.dofs = [1];
dirichlet.value = 0;
mesh = AddDirichletBC(mesh,dirichlet);
dirichlet.point_1 = [0.999*mesh.Lx, -0.001];
dirichlet.point_2 = [1.001*mesh.Lx, 0.001];
dirichlet.dofs = [2];
dirichlet.value = 0;
mesh = AddDirichletBC(mesh,dirichlet);

% Add Neumann boundary conditions
neumann.point_1 = [-0.001, -0.001];
neumann.point_2 = [0.001, 0.001];
neumann.dofs = [2];
neumann.value = -1;
mesh = AddNeumannBC(mesh,neumann);

% Prepare FEA
fea = FEA(mesh);
fea = AddSolidMaterial(fea, E, nu, 1.0);
fea = AssemblePointLoads(fea);
fea = BuildFilterMatrix(fea,radius);

% Prepare TOBS
tobs = TOBS(compliance_constraint, epsilons, flip_limits, size(fea.mesh.incidence,1));

%% --------------------------------------------------------------------- %%
%                           ** Optimization **                            %
%-------------------------------------------------------------------------%

disp(' ')
disp('         *** Optimization loop *** ')
disp(' ')

% Iteration 0 %
loop = 0;

% Finite Element analysis
densities = tobs.design_variables;
densities(densities == 0) = rho_min;
fea = AssembleStructuralK(fea,densities);
fea = SolveStaticStructuralFEA(fea);

% Objective (volume) and sensitivites
tobs.objective = mean(tobs.design_variables);
tobs.objective_sensitivities = ones(size(fea.mesh.incidence,1),1)/size(fea.mesh.incidence,1);

% Constraint (compliance) and sensitivities
tobs.constraints = fea.F'*fea.U;
tobs.constraints_sensitivities = ComputeComplianceSensitivities(fea,tobs.design_variables);
% Filtering sensitivities
tobs.constraints_sensitivities = fea.H*tobs.constraints_sensitivities;

% PlotScalarPerElement(fea,tobs.constraints_sensitivities,1/5);

% Storing sensitivities for the next iteration averaging
sensitivities_previous = tobs.constraints_sensitivities;

% Storing optimization history
tobs.history(loop+1,1) = tobs.objective;
tobs.history(loop+1,2) = tobs.constraints;

% % Saving iteration results
% directory = strcat('Results/iter',num2str(loop));
% save(directory);

% Print optimization status on the screen
disp([' It.: ' sprintf('%3i',loop) '  Obj.: ' sprintf('%5.4f',full(tobs.objective))...
    '  Comp.: ' sprintf('%3.3f',full(tobs.constraints))])

% Optimization loop %

% Convergence identifier
is_converged = 0;
difference = 1;

while (is_converged == 0)

    % Iteration counter update
    loop = loop + 1;

    % Solve with ILP
    tobs = SolveWithILP(tobs);

    % Finite Element analysis
    densities = tobs.design_variables;
    densities(densities == 0) = rho_min;
    fea = AssembleStructuralK(fea,densities);
    fea = SolveStaticStructuralFEA(fea);

    % Objective (volume) and sensitivites
    tobs.objective = mean(tobs.design_variables);
    tobs.objective_sensitivities = ones(size(fea.mesh.incidence,1),1)/size(fea.mesh.incidence,1);

    % Constraint (compliance) and sensitivities
    tobs.constraints = fea.F'*fea.U;
    tobs.constraints_sensitivities = ComputeComplianceSensitivities(fea,tobs.design_variables);
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
    end

%     % Saving iteration results
%     directory = strcat('Results/iter',num2str(loop));
%     history = tobs.history;
%     save(directory,'densities','history');

    % Print optimization status on the screen
    disp([' It.: ' sprintf('%3i',loop) '  Obj.: ' sprintf('%5.4f',full(tobs.objective))...
        '  Comp.: ' sprintf('%3.3f',full(tobs.constraints))...
        '  Conv.: ' sprintf('%4.4f',difference)])

    % Plot densities
    PlotScalarPerElement(fea,-tobs.design_variables,1); colormap(gray); pause(1e-6)

    % Stop at maximum iterations
    if (loop == 300)
        break
    end

end

% Plot convergence history
PlotVolumeHistory(full(tobs.history(:,1)));
PlotComplianceHistory(full(tobs.history(:,2)));

disp(' ')
disp('         ***************************')
disp('         **   Topology Pi - END   **')
disp('         ***************************')
disp(' ')

end
