%% --------------------------------------------------------------------- %%
%            ** Topology Optimization using the BESO method **            %
%-------------------------------------------------------------------------%

% Renato Picelli (University of São Paulo, Santos, Brazil)

% MATLAB code for structural topology optimization using the BESO method
% Find minimum structural compliance subject to a volume constraint

% 03.06.2019

function MBB_120x40_BESO

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
radius = 6;                  % Filter radius in length unit
rho_min = 0.001^3;           % Minimum density (for void elements)
final_volume_fraction = 0.5; % Final volume fraction
ER = 0.01;                   % Evolutionary ratio
ARmax = 0.05;                % Maximum admission ratio
tau = 0.001;                 % Convergence tolerance

%% --------------------------------------------------------------------- %%
%                         ** Problem set up **                            %
%-------------------------------------------------------------------------%

% Set up mesh.
mesh = Mesh('MBB_120x40_120x40ele');

% Prepare FEA
fea = FEA(mesh);
fea = AddSolidMaterial(fea, E, nu, 1.0);
fea = AssemblePointLoads(fea);
fea = BuildFilterMatrix(fea,radius);

% Prepare BESO
beso = BESO(final_volume_fraction, ER, ARmax, size(fea.mesh.incidence,1));

%% --------------------------------------------------------------------- %%
%                           ** Optimization **                            %
%-------------------------------------------------------------------------%

disp(' ')
disp('         *** Optimization loop *** ')
disp(' ')

% Iterations counter
loop = 0;

% Finite Element analysis
densities = beso.design_variables;
densities(densities == 0) = rho_min;
fea = AssembleStructuralK(fea,densities);
fea = SolveStaticStructuralFEA(fea);

% Objective (compliance) and sensitivites
beso.objective = fea.F'*fea.U;
beso.objective_sensitivities = ComputeComplianceSensitivities(fea,beso.design_variables);
% Filtering sensitivities
beso.objective_sensitivities = fea.H*beso.objective_sensitivities;

% PlotScalarPerElement(fea,beso.objective_sensitivities,1/5);

% Storing sensitivities for the next iteration averaging
sensitivities_previous = beso.objective_sensitivities;

% Current volume fraction
beso.volume_fraction = mean(beso.design_variables);

% Storing optimization history
beso.history(loop+1,1) = beso.objective;
beso.history(loop+1,2) = beso.volume_fraction;

% % Saving iteration results
% directory = strcat('Results/iter',num2str(loop));
% save(directory);

% Print optimization status on the screen
disp([' It.: ' sprintf('%3i',loop) '  Obj.: ' sprintf('%5.4f',full(beso.objective))...
    '  Vol.: ' sprintf('%3.3f',beso.volume_fraction)])

% Optimization loop %

% Convergence identifiers
is_converged = 0;
difference = 1;

while (is_converged == 0)

    % Iteration counter update
    loop = loop + 1;

    % Solve with BESO design update scheme
    beso = BESODesignUpdate(beso);

    % Finite Element analysis
    densities = beso.design_variables;
    densities(densities == 0) = rho_min;
    fea = AssembleStructuralK(fea,densities);
    fea = SolveStaticStructuralFEA(fea);

    % Objective and sensitivites
    beso.objective = fea.F'*fea.U;
    beso.objective_sensitivities = ComputeComplianceSensitivities(fea,beso.design_variables);
    % Filtering sensitivities
    beso.objective_sensitivities = fea.H*beso.objective_sensitivities;

    % Evaluating average of alpha history
    beso.objective_sensitivities = (beso.objective_sensitivities+sensitivities_previous)/2;
    % Storing alphaKminus1 for next iteration
    sensitivities_previous = beso.objective_sensitivities;

    % Current volume fraction
    beso.volume_fraction = mean(beso.design_variables);

    % Storing optimization history
    beso.history(loop+1,1) = beso.objective;
    beso.history(loop+1,2) = beso.volume_fraction;

    % Convergence analysis [Huang and Xie, 2007]
    N = 5;
    if (loop >= 2*N) % Analyzing 2*N consecutive iterations
        error_1 = zeros(N,1);
        error_2 = zeros(N,1);
        for i = 1:N
            error_1(i) = beso.history(loop-i+1,1)-beso.history(loop-N-i+1,1);
            error_2(i) = beso.history(loop-i+1,1);
        end
        % Evaluating convergence
        difference = abs(sum(error_1))/sum(error_2);
        % Verifying error tolerance and if constraint is satisfied
        if ((difference <= tau) && (beso.volume_fraction <= 1.001*beso.final_volume_fraction))
            is_converged = 1;
        end
    end

%     % Saving iteration results
%     directory = strcat('Results/iter',num2str(loop));
%     history = tobs.history;
%     save(directory,'densities','history');

    % Print results
    disp([' It.: ' sprintf('%3i',loop) '  Obj.: ' sprintf('%5.4f',full(beso.objective))...
        '  Vol.: ' sprintf('%3.3f',beso.volume_fraction)...
        '  Conv.: ' sprintf('%4.4f',difference)])

    % Plot densities
    PlotScalarPerElement(fea,-beso.design_variables,1); colormap(gray); pause(1e-6)

    % Stop at maximum iterations
    if (loop == 300)
        break
    end

end

% Plot convergence history
PlotComplianceHistory(full(beso.history(:,1)));
PlotVolumeHistory(full(beso.history(:,2)));

disp(' ')
disp('         ***************************')
disp('         **   Topology Pi - END   **')
disp('         ***************************')
disp(' ')

end
