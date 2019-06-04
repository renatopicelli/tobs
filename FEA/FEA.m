%% --------------------------------------------------------------------- %%
%                            ** FEA class **                              %
%-------------------------------------------------------------------------%

classdef FEA
    
    %% Properties
    properties
        
        % Finite element mesh
        mesh
        
        % Material properties
        E         % Young's modulus
        nu        % Poisson's ratio
        rho       % Solid material density
        th = 1.0; % thickness
        
        % Amount of added materials
        material_counter;
        % List of element materials
        materials_list;
        
        % ID matrix (relation node-DOF)
        % (lines = DOF's, columns = nodes)
        ID
        number_of_equations
        
        % List of elements inside design domain
        design_domain
        
        % Equation matrices
        F % forces
        K % stiffness
        M % mass
        U % displacements (state variable)
        
        % Sensitivity numerical filter matrix
        H
        
    end
    
    %% Methods
    methods
        
        %% Constructor
        function fea = FEA(mesh)
            
            % Generate finite element mesh
            fea.mesh = mesh;
            
            disp([' '])
            disp(['         Preparing FEA.'])
            
            % Creating ID matrix (relation node-DOF)% Initializing ID matrix: dofs = (u_x, u_y)
            fea.ID = ones(2,size(fea.mesh.coordinates,1));
            fea.ID(1,:) = 1+2*[0:(size(fea.mesh.coordinates,1)-1)];
            fea.ID(2,:) = 2+2*[0:(size(fea.mesh.coordinates,1)-1)];

            % Number of equations
            fea.number_of_equations = fea.ID(2,end);
            
            disp(['         Number of equations: ', sprintf('%10i',fea.number_of_equations)])
            
            % Default design domain (all solid elements)
            fea.design_domain = find(fea.mesh.incidence(:,1) == 1);
            
            % Initializing material_counter and materials_list
            fea.material_counter = 0;
            fea.materials_list = ones(size(fea.mesh.incidence,1),1);
            
        end % end Constructor
        
        %% Add solid material
        function fea = AddSolidMaterial(fea, E_in, nu_in, rho_in)
            
            % Material counter
            fea.material_counter = fea.material_counter+1;
            
            % Assign material
            fea.E(fea.material_counter) = E_in;
            fea.nu(fea.material_counter) = nu_in;
            fea.rho(fea.material_counter) = rho_in;
            
        end % end AddSolidMaterial
        
        %% Build filter matrix
        function fea = BuildFilterMatrix(fea,radius)
            
            disp(['         Building filter matrix.'])
            
            % Estimating the number of elements in filter matrix: elements within filter
            % area multiplied by the number of elements.
            nfilterel = length(fea.design_domain)*4*round(radius/(fea.mesh.coordinates(fea.mesh.incidence(1,3),1)-fea.mesh.coordinates(fea.mesh.incidence(1,2),1)))^2;

            % Allocating vectors
            % Element weights
            R = zeros(nfilterel,1);
            % Index vectors
            ifilt = zeros(nfilterel,1);
            jfilt = zeros(nfilterel,1);

            % Selecting coordinates of design nodes
            X_den = fea.mesh.centroids(fea.design_domain,1);
            Y_den = fea.mesh.centroids(fea.design_domain,2);

            % Auxiliary counter
            k = 0;

            % Looping at the elements inside the design domain
            for i = 1:length(fea.design_domain)

                % Element
                el = fea.design_domain(i);

                % Coordinates of the centroid of the element
                x = fea.mesh.centroids(el,1);
                y = fea.mesh.centroids(el,2);

                % Weight function (linear with radial distance)
                r = max(0,1-(((X_den-x).^2+(Y_den-y).^2).^(1/2)/radius));

                % Auxiliary counter update
                k = k+nnz(r);

                % Index vectors building
                ifilt(k-nnz(r)+1:k) = i*ones(nnz(r),1);
                jfilt(k-nnz(r)+1:k) = find(r);

                % Weights update
                R(k-nnz(r)+1:k) = r(r>0)/sum(r);

            end

            % H filter matrix
            fea.H = sparse(ifilt(1:k),jfilt(1:k),R(1:k),length(fea.design_domain),length(fea.design_domain));

        end % end BuildFilterMatrix
        
        %% Assemble point loads
        function fea = AssemblePointLoads(fea)
            
            % Load vector assembly
            fea.F = sparse(zeros(2,1));
            fea.F(fea.number_of_equations) = 0;
            % Applied forces
            for i = 1:size(fea.mesh.neumann_boundary,1)
                if (fea.mesh.neumann_boundary(i,1) ~= 0)
                    fea.F(fea.ID(fea.mesh.neumann_boundary(i,2),fea.mesh.neumann_boundary(i,1)),1) = fea.mesh.neumann_boundary(i,3);
                end
            end
            
        end % end AssemblePointLoads
        
        %% Compute structural element stiffness matrix
        function Ke = ComputeStructuralKe(fea,stress_case,material_index)

            % Initializing matrix and getting FEA info
            Ke = zeros(8);
            el = 1; % Structure grid, element 1
            coord = fea.mesh.coordinates;
            inci = fea.mesh.incidence;
            E = fea.E(material_index);
            nu = fea.nu(material_index);

            % Selecting stress case
            switch stress_case
                case 1 % Plane stress elasticity matrix
                    D = (E/(1-nu*nu))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
                case 2
                    disp('Plane-strain not implemented yet!')
            end
            
            % Gaussian points (2 x 2)
            W = (1/3)*[sqrt(3) sqrt(3); sqrt(3) -sqrt(3); -sqrt(3) sqrt(3); -sqrt(3) -sqrt(3)];

            % Integral
            for j = 1:4
                N = W(j,1);  % Eta
                Q = W(j,2);  % Xi

                % Shape functions derivatives with respect to Xi(Q)
                N1r = (1/4)*(N-1);
                N2r = (1/4)*(1-N);
                N3r = (1/4)*(1+N);
                N4r = -(1/4)*(1+N);

                % Shape functions derivatives with respect to Eta(N)
                N1s = (1/4)*(Q-1);
                N2s = -(1/4)*(1+Q);
                N3s = (1/4)*(1+Q);
                N4s = (1/4)*(1-Q);

                % Coeficientes da Matriz Jacobiana
                J11 = N1r*coord(inci(el,2),1)+N2r*coord(inci(el,3),1)+N3r*coord(inci(el,4),1)+N4r*coord(inci(el,5),1);
                J12 = N1r*coord(inci(el,2),2)+N2r*coord(inci(el,3),2)+N3r*coord(inci(el,4),2)+N4r*coord(inci(el,5),2);
                J21 = N1s*coord(inci(el,2),1)+N2s*coord(inci(el,3),1)+N3s*coord(inci(el,4),1)+N4s*coord(inci(el,5),1);
                J22 = N1s*coord(inci(el,2),2)+N2s*coord(inci(el,3),2)+N3s*coord(inci(el,4),2)+N4s*coord(inci(el,5),2);

                % Jacobian matrix
                J = [J11 J12; J21 J22];

                % Auxiliary matrix [G]
                jacob = det(J);
                % solid
                G = (1/jacob)*[J22 0 -J12 0; 0 -J21 0 J11; -J21 J22 J11 -J12];

                % Auxiliary matrix [P]
                % Solid
                P = [N1r 0 N2r 0 N3r 0 N4r 0; 0 N1r 0 N2r 0 N3r 0 N4r; N1s 0 N2s 0 N3s 0 N4s 0; 0 N1s 0 N2s 0 N3s 0 N4s];

                % Matrix [B] with the derivatives of shape functions
                B = G*P;

                % Stiffness matrix integral
                Ke = Ke + fea.th*B'*D*B*jacob;

            end
            
        end % end ComputeStructuralKe
        
        %% Assemble stiffness matrix
        function fea = AssembleStructuralK(fea,densities)
            
            clear K
            
            % Compute element stiffness matrix
            Ke = ComputeStructuralKe(fea,1,1);
            
            % Number of elements
            nel = size(fea.mesh.incidence,1);
            
            % Assemble global stiffness matrix
            % Solid elements DOF's
            loc = [fea.ID(1,fea.mesh.incidence([1:nel]',2));
                   fea.ID(2,fea.mesh.incidence([1:nel]',2));
                   fea.ID(1,fea.mesh.incidence([1:nel]',3));
                   fea.ID(2,fea.mesh.incidence([1:nel]',3));
                   fea.ID(1,fea.mesh.incidence([1:nel]',4));
                   fea.ID(2,fea.mesh.incidence([1:nel]',4));
                   fea.ID(1,fea.mesh.incidence([1:nel]',5));
                   fea.ID(2,fea.mesh.incidence([1:nel]',5))];

            % Index vectors
            I = reshape(repmat(loc,8,1),nel*64,1);
            J = kron(loc(:),ones(8,1));
            Kg = repmat(Ke(:),nel,1);

            % Mulyiplying physical densities
            for i = 1:nel
                Kg([1:64]+64*(i-1),1) = Kg([1:64]+64*(i-1),1)*densities(i);
            end

            % Assembly
            fea.K = sparse(I,J,Kg);
            
        end % end AssembleStructuralK
        
        %% Solve static case of FEA (Ku = F)
        function fea = SolveStaticStructuralFEA(fea)
            
            % Identifying active DOF's
            dofn = [1:fea.number_of_equations]';
            % Blocking DOF's with dirichlet boundary conditions = 0
            for i = 1:size(fea.mesh.dirichlet_boundary,1)
                if (fea.mesh.dirichlet_boundary(i,1) ~= 0) && (fea.mesh.dirichlet_boundary(i,3) == 0)
                    dofn(fea.ID(fea.mesh.dirichlet_boundary(i,2),fea.mesh.dirichlet_boundary(i,1))) = 0;
                end
            end
            % Active DOF's
            dofa = nonzeros(dofn);
            
            % Solving system responses
            fea.U = sparse(zeros(10,1));
            fea.U(fea.number_of_equations) = 0;

            % Matlab linear solver
            fea.U(dofa) = fea.K(dofa,dofa)\fea.F(dofa);
            
        end % end SolveStaticFEA
        
        %% Compute binary compliance sensitivities
        function [sensitivities,objective] = ComputeComplianceSensitivities(fea,design_variables)

            % Initializing sensitivity vector
            sensitivities = zeros(length(fea.design_domain),1);
            
            % Objective function
            objective = 0;
            
            % Compute element stiffness matrix
            Ke = ComputeStructuralKe(fea,1,1);

            % For elements in the design domain
            for i = 1:length(fea.design_domain)

                if (design_variables(fea.design_domain(i)) == 1)
                    % Element from design domain
                    eleD = fea.design_domain(i);

                    % Element nodes
                    no1 = fea.mesh.incidence(eleD,2);
                    no2 = fea.mesh.incidence(eleD,3);
                    no3 = fea.mesh.incidence(eleD,4);
                    no4 = fea.mesh.incidence(eleD,5);

                    % Vector with element DOF's
                    loc = [fea.ID(1,no1),fea.ID(2,no1),fea.ID(1,no2),fea.ID(2,no2),...
                           fea.ID(1,no3),fea.ID(2,no3),fea.ID(1,no4),fea.ID(2,no4)];

                    % Element's displacement vector
                    Un = full(fea.U(loc));

                    % Sensitivity number with loading sensitivities
                    sensitivities(eleD) = -Un'*Ke*Un;

                    % Summing up objective
                    objective = objective + Un'*Ke*Un;
                end

            end
            
        end % end ComputeComplianceSensitivities
        
    end % end methods
end