%% --------------------------------------------------------------------- %%
%          ** FEA class with fluid-structure interaction **               %
%-------------------------------------------------------------------------%

classdef FluidStructureInteraction
    
    %% Properties
    properties
        
        % Finite element mesh
        mesh
        
        % Material properties
        E         % Young's modulus
        nu        % Poisson's ratio
        rho_s     % Solid material density
        th = 1.0; % thickness
        
        % Amount of added solid materials
        material_counter;
        % List of element solid materials
        materials_list;
        
        % ID matrix (relation node-DOF)
        % (lines = DOF's, columns = nodes)
        ID
        number_of_equations
        
        % Element type (0 = void, 1 = solid, 2 = fluid)
        element_types
        
        % List of element inside design domain
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
        function fea = FluidStructureInteraction(mesh)
            
            % Generate finite element mesh
            fea.mesh = mesh;
            
            disp([' '])
            disp(['         Preparing FEA.'])
            
            % Creating ID matrix (relation node-DOF)
            % Initializing ID matrix: dofs = (u_x, u_y, u_z, p_f, v_x, v_y, v_z)
            fea.ID = ones(2,size(fea.mesh.coordinates,1));
            fea.ID(1,:) = 1+2*[0:(size(fea.mesh.coordinates,1)-1)];
            fea.ID(2,:) = 2+2*[0:(size(fea.mesh.coordinates,1)-1)];
            fea.ID(3,:) = zeros(1,size(fea.mesh.coordinates,1));
            fea.ID(4,:) = [(fea.ID(2,end)+1):(fea.ID(2,end)+size(fea.mesh.coordinates,1))];

            % Number of equations
            fea.number_of_equations = fea.ID(4,end);
            
            % Initial set of element types (1 = solid, 2 = fluid)
            fea.element_types = fea.mesh.incidence(:,1);
            
            % Default design domain (all solid elements)
            fea.design_domain = find(fea.mesh.incidence(:,1) == 1);
            
            % Initializing material_counter and materials_list of solid elements
            fea.material_counter = 0;
            fea.materials_list = ones(size(fea.mesh.incidence,1),1);
            
            disp(['         Number of equations: ', sprintf('%10i',fea.number_of_equations)])
            
        end % end Constructor
        
        %% Add solid material
        function fea = AddSolidMaterial(fea, E_in, nu_in, rho_in)
            
            % Material counter
            fea.material_counter = fea.material_counter+1;
            
            % Assign material
            fea.E(fea.material_counter) = E_in;
            fea.nu(fea.material_counter) = nu_in;
            fea.rho_s(fea.material_counter) = rho_in;
            
        end % end AddSolidMaterial
        
        %% Build filter matrix
        function H = BuildFilterMatrix(fea,radius)
            
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
            H = sparse(ifilt(1:k),jfilt(1:k),R(1:k),length(fea.design_domain),length(fea.design_domain));

        end % end BuildFilterMatrix
        
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
        
        %% Compute hydrostatic fluid element stiffness matrix
        function Kef = ComputeHydrostaticKef(fea)

            % Initializing matrix and getting FEA info
            Kef = zeros(4);
            el = 1; % Structure grid, element 1
            coord = fea.mesh.coordinates;
            inci = fea.mesh.incidence;
            
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
                Gf = (1/jacob)*[J22 -J12; -J21 J11];

                % Auxiliary matrix [P]
                Pf = [N1r N2r N3r N4r; N1s N2s N3s N4s;];

                % Matrix [B] with the derivatives of shape functions
                Bf = Gf*Pf;

                % Stiffness matrix integral
                Kef = Kef + (Bf')*Bf*jacob;

            end
            
        end % end ComputeHydrostaticKef
        
        %% Assemble stiffness matrix
        function fea = AssembleCoupledStructuralAndHydrostaticK(fea)
            
            clear fea.K
            
            % Local incidence matrix
            inci = fea.mesh.incidence;
            
            % Compute element stiffness matrices
            Ke = ComputeStructuralKe(fea,1,1);
            Kef = ComputeHydrostaticKef(fea);
            
            % Number of solid and fluid elements
            nelsol = nnz(fea.element_types == 1);
            nelflu = nnz(fea.element_types == 2);

            % Fluid elements DOF's
            loc = [fea.ID(4,inci(find(fea.element_types == 2),2));
                   fea.ID(4,inci(find(fea.element_types == 2),3));
                   fea.ID(4,inci(find(fea.element_types == 2),4));
                   fea.ID(4,inci(find(fea.element_types == 2),5))];

            % Index vectors
            If = reshape(repmat(loc,4,1),nelflu*16,1);
            Jf = kron(loc(:),ones(4,1));
            Kf = repmat(Kef(:),nelflu,1);

            % Solid elements DOF's
            loc = [fea.ID(1,inci(find(fea.element_types == 1),2));
                   fea.ID(2,inci(find(fea.element_types == 1),2));
                   fea.ID(1,inci(find(fea.element_types == 1),3));
                   fea.ID(2,inci(find(fea.element_types == 1),3));
                   fea.ID(1,inci(find(fea.element_types == 1),4));
                   fea.ID(2,inci(find(fea.element_types == 1),4));
                   fea.ID(1,inci(find(fea.element_types == 1),5));
                   fea.ID(2,inci(find(fea.element_types == 1),5))];

            % Index vectors
            Is = reshape(repmat(loc,8,1),nelsol*64,1);
            Js = kron(loc(:),ones(8,1));
            Ks = repmat(Ke(:),nelsol,1);

            % Global index vectors
            I = [If; Is];
            J = [Jf; Js];
            Kg = [Kf; Ks];

            % Assembly
            K = sparse(I,J,Kg);
            
            % Forcing K to store total number of equations
            if (length(K) < fea.number_of_equations)
                K(fea.number_of_equations,fea.number_of_equations) = 0;
            end
            
            % Couple K matrix
            fea.K = CoupleStructuralAndHydrostaticK(fea,K);
            
        end % end AssembleStructuralAndHydrostaticCoupledK
        
        %% Couple structural and hydrostatic global stiffness matrices
        function K = CoupleStructuralAndHydrostaticK(fea,K)

            % Identifying fluid elements with solid element neighbours
            elfluid = find(fea.element_types == 2);
            for i = 1:length(elfluid)
                % Fluid element to be analyzed
                fluid = elfluid(i);
                % Initializing vector with element types
                ele_type = fea.element_types(nonzeros(fea.mesh.element_neighbours(fluid,2:5)));
                % Identifying and blocking DOF's
                if isempty(find(ele_type == 1))
                    % There is not any solid neighbour element
                    elfluid(i) = 0;
                end
            end
            % fluid elements with solid element neighbours
            elfluid = nonzeros(elfluid);

            % Interface identification, coupling matrix evaluation and coupling of
            % global matrix
            for i = 1:length(elfluid)

                % Fluid element
                fluid = elfluid(i);

                % Finding solid elements neighbours of fluid element i 
                elinterf = fea.mesh.element_neighbours(fluid,find(fea.element_types(nonzeros(fea.mesh.element_neighbours(fluid,2:5))) == 1)+1);

                % For each interface element
                for j = 1:length(elinterf)

                    % Auxiliary variable
                    aux_int = 0;

                    % Initializing interface
                    nodes_interface = 0;

                    % For the 4 nodes of connectivity of the fluid element i
                    for k = 2:5
                        % Searching interface node
                        aux_no = fea.mesh.incidence(fluid,k);
                        % If aux_no is in elinterf connectivity
                        if (length(find(fea.mesh.incidence(elinterf(j),2:5) == aux_no)) >= 1)  
                            % Interface nodes
                            aux_int = aux_int+1;
                            nodes_interface(aux_int) = aux_no;
                        end
                    end

                    % Verifying if the interface is between the last and the first node
                    % of the fluid connectivity. Forcing the interface to be
                    % counter-clockwise with respect to the fluid connectivity
                    if (fea.mesh.incidence(fluid,[2 5]) == nodes_interface)
                        nodes_interface = [ nodes_interface(2) nodes_interface(1) ];
                    end

                    % If there is interface
                    if (nodes_interface ~= 0)
                        % Coupling matrix evaluation
                        Le = ComputeInterfaceCouplingMatrix(nodes_interface,fea);
                        % Interface element DOF's
                        loc = [fea.ID(1,nodes_interface(1)) fea.ID(2,nodes_interface(1)) fea.ID(1,nodes_interface(2))...
                               fea.ID(2,nodes_interface(2)) fea.ID(4,nodes_interface(1)) fea.ID(4,nodes_interface(2))];
                        % Coupling matrix assembly
                        K(loc(1:4),loc) = K(loc(1:4),loc)-Le;
                    end
                end
            end
        end % end CoupleStructuralAndHydrostaticK
        
        %% Compute interface coupling matrix
        function Le = ComputeInterfaceCouplingMatrix(nodes_interface,fea)

            % Interface nodes
            node1 = nodes_interface(1);
            node2 = nodes_interface(2);

            % Nodal coordinates
            x1 = fea.mesh.coordinates(node1,1);
            y1 = fea.mesh.coordinates(node1,2);
            x2 = fea.mesh.coordinates(node2,1);
            y2 = fea.mesh.coordinates(node2,2);

            % Interface element length
            Lface = sqrt(((x2-x1)^2)+((y2-y1)^2));

            % Switched sinus and cossinus of the normal vector
            sinA = (y2-y1)/Lface;    % Cossinus
            cosA = -(x2-x1)/Lface;   % Sinus

            % Coupling matrix
            Le = (Lface/6)*[ 0 0 0 0  2*sinA    sinA;
                             0 0 0 0  2*cosA    cosA;
                             0 0 0 0    sinA  2*sinA;
                             0 0 0 0    cosA  2*cosA];

        end % end ComputeInterfaceCouplingMatrix
        
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

        %% Solve static case of FEA (Ku = F)
        function fea = SolveFEA(fea)
            
            % Impose non-zero DOF's
            for i = 1:size(fea.mesh.dirichlet_boundary,1)
                if (fea.mesh.dirichlet_boundary(i,1) ~= 0) % If there is an imposed DOF
                    if (fea.mesh.dirichlet_boundary(i,3) ~= 0) % If the imposed DOF is non-zero
                        % Penalty
                        Pform = 1e10;
                        % For imposed pressures
                        if (fea.mesh.dirichlet_boundary(i,2) == 4)
                            % Imposed DOF
                            dofpen = fea.ID(fea.mesh.dirichlet_boundary(i,2),fea.mesh.dirichlet_boundary(i,1));
                        else
                            disp('ERROR: not ready for prescribed u.')
                            break;
                        end
                        % Imposing DOF by changing F and K
                        fea.F(dofpen) = fea.F(dofpen) + Pform*fea.K(dofpen,dofpen)*fea.mesh.dirichlet_boundary(i,3);
                        fea.K(dofpen,dofpen) = Pform*fea.K(dofpen,dofpen);
                    end
                end
            end
            
            % Identifying active DOF's
            dofn = [1:fea.number_of_equations]';
            for i = 1:size(fea.mesh.coordinates,1)
                % Types of elements conected to the node i
                ele_type = fea.element_types(nonzeros(fea.mesh.nodal_connectivity(i,2:5)));
                % Identifying and blocking DOF's
                if (isempty(find(ele_type == 1))) && (isempty(find(ele_type == 2)))
                    % Only void elements connected to the node i
                    dofn(fea.ID(4,i)) = 0;
                    dofn(fea.ID(1,i)) = 0;
                    dofn(fea.ID(2,i)) = 0;
                end
                if (isempty(find(ele_type == 1)))
                    % Only fluid elements connected to the node i
                    dofn(fea.ID(1,i)) = 0;
                    dofn(fea.ID(2,i)) = 0;
                elseif (isempty(find(ele_type == 2)))
                    % Only solid elements connected to the node i
                    dofn(fea.ID(4,i)) = 0;
                end
            end
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
        
        %% Compute binary compliance sensitivities with design-dependent coupled hydrostatic pressure loads
        function [sensitivities,objective] = ComputeComplianceWithCoupledPressureSensitivities(fea)

            % Initializing sensitivity vector
            sensitivities = zeros(length(fea.design_domain),1);
            
            % Objective function
            objective = 0;
            
            % Size of the elements
            ele_size = (fea.mesh.coordinates(fea.mesh.incidence(1,3))-fea.mesh.coordinates(fea.mesh.incidence(1,2)));

            % Lc matrix
            Lc = (ele_size/6)*[ 2  0  0  1;
                                2  1  0  0;
                                0 -2 -1  0;
                                1  2  0  0;
                                0 -1 -2  0;
                                0  0 -2 -1;
                                1  0  0  2;
                                0  0 -1 -2];
            
            % Compute element stiffness matrix
            Ke = ComputeStructuralKe(fea,1,1);

            % For elements in the design domain
            for i = 1:length(fea.design_domain)

                if (fea.element_types(fea.design_domain(i)) == 1)
                    
                    % Element from design domain
                    eleD = fea.design_domain(i);

                    % Element nodes
                    node1 = fea.mesh.incidence(eleD,2);
                    node2 = fea.mesh.incidence(eleD,3);
                    node3 = fea.mesh.incidence(eleD,4);
                    node4 = fea.mesh.incidence(eleD,5);

                    % Vector with displacement DOF's
                    loc = [fea.ID(1,node1),fea.ID(2,node1),fea.ID(1,node2),fea.ID(2,node2),...
                           fea.ID(1,node3),fea.ID(2,node3),fea.ID(1,node4),fea.ID(2,node4)];

                    % Element's displacement vector
                    Un = full(fea.U(loc));
                    
                    % Vector with pressure DOF's
                    loc = [fea.ID(4,node1),fea.ID(4,node2),fea.ID(4,node3),fea.ID(4,node4)];

                    % Hydrostatic pressure
                    P = full(max(fea.U(loc)));
                    P = [P; P; P; P];

                    % Sensitivity number with loading sensitivities
                    sensitivities(i) = -0.5*Un'*Ke*Un+Un'*Lc*P;

                    % Summing up objective
                    objective = objective + 0.5*Un'*Ke*Un;
                end

            end
            
        end % end ComputeComplianceWithCoupledPressureSensitivities
        
        %% Fluid flooding process
        function fea = FluidFlooding(fea,design_variables)
            
            % Return design domain to 0/1 configuration
            fea.element_types(fea.design_domain) = design_variables;

            % Propagate fluid region
            % Setting auxiliary counter
            aux = 1;
            % While fluid region still advances
            while (aux ~= 0)
                % Resetting auxiliary counter
                aux = 0;
                % Advancing fluid region
                for i = 1:length(fea.element_types)
                    % For fluid elements
                    if (fea.element_types(i) == 2)
                        % Identifying type of neighbour elements
                        for j = 2:5
                            % If neighbour element exists
                            if (fea.mesh.element_neighbours(i,j) ~= 0)
                                % If neighbour elements is void
                                if (fea.element_types(fea.mesh.element_neighbours(i,j),1) == 0)
                                    % Neighbour element becomes fluid
                                    fea.element_types(fea.mesh.element_neighbours(i,j),1) = 2;
                                    % Auxiliary counter
                                    aux = aux+1;
                                end
                            end
                        end
                    end
                end
            end
        end % end FluidFlooding
        
        %% Identify list of solid elements with fluid neighbours
        function solids_at_interface = ListSolidsAtInterface(fea)

            % List of solid elements
            solids = find(fea.element_types == 1);
            
            % loop through solid elements
            for i = 1:length(solids)
                % Element types of the neighbouring elements
                ele_type = fea.element_types(nonzeros(fea.mesh.element_neighbours(solids(i),2:5)));
                % If there is no fluid element
                if isempty(find(ele_type == 2))
                    % Deleting element from the list
                    solids(i) = 0;
                end
            end
            
            % Solid element with fluid neighbours
            solids_at_interface = nonzeros(solids);
                
        end % end ListSolidsAtInterface
        
    end % end methods
end