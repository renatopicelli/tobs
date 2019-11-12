%% --------------------------------------------------------------------- %%
%                       ** Structural analysis **                         %
%-------------------------------------------------------------------------%

% Renato Picelli (University of São Paulo, Santos, Brazil)

% 06.06.2019

function structural_analysis

close all; clear all; clc

addpath(genpath(strcat('..', filesep, '..', filesep, 'FEA')))

disp(' ')
disp('         ***************')
disp('         **   START   **')
disp('         ***************')
disp(' ')

%% --------------------------------------------------------------------- %%
%                              ** Input **                                %
%-------------------------------------------------------------------------%

% Material properties
E = 1.0;   % Young's modulus
nu = 0.3;  % Poisson's ratio

%% --------------------------------------------------------------------- %%
%                         ** Problem set up **                            %
%-------------------------------------------------------------------------%

% Set up mesh. Mesh() means no input mesh will be read.
mesh = Mesh();

% Set up discretization size.
mesh.Lx = 20.0;
mesh.Ly = 60.0;
mesh.nelx = 20;
mesh.nely = 60;
% Generate regular quadrilateral grid.
mesh = GenerateQuad4RectangleMesh(mesh);

% Add Dirichlet boundary conditions
dirichlet.point_1 = [-0.001, -0.001];
dirichlet.point_2 = [mesh.Lx*1.001, 0.001];
dirichlet.dofs = [1, 2];
dirichlet.value = 0.0;
mesh = AddDirichletBC(mesh,dirichlet);

% Add Neumann boundary conditions
neumann.point_1 = [9.999, 59.999];
neumann.point_2 = [10.001, 60.0001];
neumann.dofs = [1];
neumann.value = 1.0;
mesh = AddNeumannBC(mesh,neumann);

% Prepare FEA
fea = FEA(mesh);
fea = AddSolidMaterial(fea, E, nu, 1.0);

%% --------------------------------------------------------------------- %%
%                     ** Finite Element Analysis **                       %
%-------------------------------------------------------------------------%

% Assemble force vector
fea = AssemblePointLoads(fea);

% Assemble stiffness matrix
densities = ones(length(fea.mesh.incidence),1);
fea = AssembleStructuralK(fea,densities);

% Solve equilibrium equation
fea = SolveStaticStructuralFEA(fea);

% Plot structural displacements
PlotStructuralDisplacements(fea,1e-1)

disp(' ')
disp('         *************')
disp('         **   END   **')
disp('         *************')
disp(' ')

end
