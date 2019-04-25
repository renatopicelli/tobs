%% --------------------------------------------------------------------- %%
%         ** Topology Optimization of Binary Structures (TOBS) **         %
%-------------------------------------------------------------------------%

% Renato Picelli (University of São Paulo, Santos, Brazil)
% Raghavendra Sivapuram (University of California, San Diego, USA)

% MATLAB code for structural topology optimization using the TOBS method
% Find minimum structural compliance subject to a volume constraint

% 14.03.2019

function MBB_120x40

close all; clear all; clc

addpath(genpath('../../FEA'))
addpath(genpath('../../Meshes'))
addpath(genpath('../../TopOpt'))

disp(' ')
disp('         *****************************')
disp('         **   Topology Pi - START   **')
disp('         *****************************')
disp(' ')

%% --------------------------------------------------------------------- %%
%                              ** Input **                                %
%-------------------------------------------------------------------------%

% Material properties
E = 1.0; % Young's modulus
v = 0.3;   % Poisson's ratio

% Optimization parameters
radius = 6;                  % Filter radius in length unit
rho_min = 0.001^3;           % Minimum density (for void elements)
final_volume_fraction = 0.5; % Final volume fraction
epsilons = 0.01;             % Constraint relaxation parameter
flip_limits = 0.05;          % Flip limits
tau = 0.001;                 % Convergence tolerance

%% --------------------------------------------------------------------- %%
%                         ** Problem set up **                            %
%-------------------------------------------------------------------------%

% % Prepare FEA
fea = FEA('MBB_120x40_120x40ele');
fea.E = E;
fea.v = v;
fea.F = AssemblePointLoads(fea);
fea.H = BuildFilterMatrix(fea,radius);

% Prepare TOBS
tobs = TOBS();
% Assign constraints relaxation parameters
tobs.epsilons = epsilons;
% Assign flip limits
tobs.flip_limits = flip_limits;
% Assign initial design variables distribution
tobs.design_variables = ones(size(fea.mesh.incidence,1),1);
% % Define constraints
tobs.constraints_limits = final_volume_fraction;

%% --------------------------------------------------------------------- %%
%                           ** Optimization **                            %
%-------------------------------------------------------------------------%

disp(' ')
disp('         *** Optimization loop *** ')
disp(' ')

% Iterations counter
loop = 0;

% Finite Element analysis
densities = tobs.design_variables; densities(densities == 0) = rho_min;
fea.K = AssembleStructuralK(fea,densities);
fea.U = SolveStaticStructuralFEA(fea);

% Objective (compliance) and sensitivites
tobs.objective = fea.F'*fea.U;
[tobs.objective_sensitivities] = ComputeComplianceSensitivities(fea,tobs.design_variables);
% Filtering sensitivities
tobs.objective_sensitivities = fea.H*tobs.objective_sensitivities;

% PlotScalarPerElement(fea,tobs.objective_sensitivities,1/5);

% Storing sensitivities for the next iteration averaging
sensitivities_previous = tobs.objective_sensitivities;

% Constraint (volume) and sensitivities
tobs.constraints = mean(tobs.design_variables);
tobs.constraints_sensitivities = ones(size(fea.mesh.incidence,1),1)/size(fea.mesh.incidence,1);

% Storing optimization history
tobs.history(loop+1,1) = tobs.objective;
tobs.history(loop+1,2) = tobs.constraints;

% % Saving iteration results
% directory = strcat('Results/iter',num2str(loop));
% save(directory);

% Print optimization status on the screen
disp([' It.: ' sprintf('%3i',loop) '  Obj.: ' sprintf('%5.4f',full(tobs.objective))...
    '  Vol.: ' sprintf('%3.3f',full(tobs.constraints))])

% Optimization loop %

% Convergence identifiers
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
    fea.K = AssembleStructuralK(fea,densities);
    fea.U = SolveStaticStructuralFEA(fea);
    
    % Objective and sensitivites
    tobs.objective = fea.F'*fea.U;
    [tobs.objective_sensitivities] = ComputeComplianceSensitivities(fea,tobs.design_variables);
    % Filtering sensitivities
    tobs.objective_sensitivities = fea.H*tobs.objective_sensitivities;
    
    % Evaluating average of alpha history
    tobs.objective_sensitivities = (tobs.objective_sensitivities+sensitivities_previous)/2;
    % Storing alphaKminus1 for next iteration
    sensitivities_previous = tobs.objective_sensitivities;

    % Constraints and sensitivities
    tobs.constraints = mean(tobs.design_variables);
    tobs.constraints_sensitivities = ones(size(fea.mesh.incidence,1),1)/size(fea.mesh.incidence,1);
    
    % Storing optimization history
    tobs.history(loop,1) = tobs.objective;
    tobs.history(loop,2) = tobs.constraints;
    
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
    
    % Print results
    disp([' It.: ' sprintf('%3i',loop) '  Obj.: ' sprintf('%5.4f',full(tobs.objective))...
        '  Vol.: ' sprintf('%3.3f',tobs.constraints)...
        '  Conv.: ' sprintf('%4.4f',difference)])

    % Plot densities
    PlotScalarPerElement(fea,-tobs.design_variables,1); colormap(gray); pause(1e-6)
    
    % Stop at maximum iterations
    if (loop == 300)
        break
    end

end

% Plot convergence history
PlotComplianceHistory(full(tobs.history(:,1)));
PlotVolumeHistory(full(tobs.history(:,2)));

disp(' ')
disp('         ***************************')
disp('         **   Topology Pi - END   **')
disp('         ***************************')
disp(' ')

end
