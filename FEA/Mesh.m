%% --------------------------------------------------------------------- %%
%                            ** Mesh class **                             %
%-------------------------------------------------------------------------%

classdef Mesh
    
    %% Properties
    properties
        
        % Nodal coordinates (lines = nodes, columns = coordinates)
        % [  coordinate_X1, coordinate_Y1, coordinate_Z1,
        %    coordinate_X2, coordinate_Y2, coordinate_Z2,
        %    ...
        %    coordinate_Xn, coordinate_Yn, coordinate_Zn  ];
        coordinates
        
        % Element incidence - connectivity (lines = elements, columns = element type and connected to node)
        % [  element_1_type, node_1, node_2, node_3, node_4,
        %    element_2_type, node_1, node_2, node_3, node_4,
        %    ...
        %    element_e_type, node_1, node_2, node_3, node_4  ];
        incidence
        
        % Dirichlet boundary conditions
        % [  node, DOF type, value  ];
        dirichlet_boundary
        
        % Neumann boundary conditions
        % [  node, DOF type, value  ];
        neumann_boundary
        
        % Centroids coordinates (lines = elements, columns = coordinates)
        % [  centroid_X1, centroid_Y1
        %    centroid_X2, centroid_Y2
        %    ...
        %    centroid_Xn, centroid_Yn  ]; NOT READY FOR 3D YET!
        centroids
        
        % Nodal incidence - connectivity (lines = nodes, columns = element connected)
        % [  node_1_number, element_1, element_2, element_3, element_4,
        %    node_2_number, element_1, element_2, element_3, element_4,
        %    ...
        %    node_n_number, element_1, element_2, element_3, element_4  ];
        nodal_connectivity
        
        % Element neighbours (lines = elements, columns = side neighbours)
        % [  element_1_number, neighbour_1, neighbour_2, neighbour_3, neighbour_4,
        %    element_2_number, neighbour_1, neighbour_2, neighbour_3, neighbour_4,
        %    ...
        %    element_n_number, neighbour_1, neighbour_2, neighbour_3, neighbour_4  ];
        element_neighbours
        
        % Matrix that relates element number and its position in a
        % rectangular matrix (m x n). Useful for plotting densities
        plot_matrix
        
    end
    
    %% Methods
    methods
        
        %% Constructor
        function mesh = Mesh(example)
            
            % Read mesh 'example.dat'
            [mesh.coordinates, mesh.incidence,...
                mesh.dirichlet_boundary, mesh.neumann_boundary] = ReadMesh(example);
            
            disp(['         Number of nodes:    ', sprintf('%10i',size(mesh.coordinates,1))])
            disp(['         Number of elements: ', sprintf('%10i',size(mesh.incidence,1))])
            
            % Compute elements centroids
            mesh.centroids = ComputeCentroids(mesh);
            
            % Build nodal connectivity
            mesh.nodal_connectivity = BuildNodalConnectivity(mesh);
            
            % Find and store element neighbours
            mesh.element_neighbours = FindElementNeighbours(mesh);
            
            % Building plot_matrix
            mesh.plot_matrix = BuildPlotMatrix(mesh);
            
        end % end Constructor
        
        %% Compute elements centroids
        function centroids = ComputeCentroids(mesh)
            
            disp([' '])
            disp(['         Computing centroids.'])
            
            % Local name for matrices.
            coordinates = mesh.coordinates;
            incidence = mesh.incidence;

            % Initializing centroids matrix
            centroids = zeros(length(incidence),2);
            
            % Calculating centroids
            for i = 1:length(incidence)

                % Maximum coordinate X of element i
                xmax = max([coordinates(incidence(i,2),1) coordinates(incidence(i,3),1) coordinates(incidence(i,4),1) coordinates(incidence(i,5),1)]);
                % Minimum coordinate X of element i
                xmin = min([coordinates(incidence(i,2),1) coordinates(incidence(i,3),1) coordinates(incidence(i,4),1) coordinates(incidence(i,5),1)]);
                % Maximum coordinate Y of element i
                ymax = max([coordinates(incidence(i,2),2) coordinates(incidence(i,3),2) coordinates(incidence(i,4),2) coordinates(incidence(i,5),2)]);
                % Minimum coordinate Y of element i
                ymin = min([coordinates(incidence(i,2),2) coordinates(incidence(i,3),2) coordinates(incidence(i,4),2) coordinates(incidence(i,5),2)]);

                % X component of element i centroid
                dx = (xmax+xmin)/2;
                % Y component of element i centroid
                dy = (ymax+ymin)/2;

                % Centroids matrix
                centroids(i,1) = dx;
                centroids(i,2) = dy;

            end
            
        end % end ComputeCentroids
        
        %% Build matrix with nodal connectivity
        function nodal_connectivity = BuildNodalConnectivity(mesh)
            
            disp(['         Building nodal connectivity.'])

            % Initializing matrix
            nodal_connectivity = zeros(size(mesh.coordinates,1),5);

            % Build matrix
            for j = 1:size(mesh.coordinates,1)

                % Elements connected to the node j
                connectel = [find(mesh.incidence(:,2) == j); find(mesh.incidence(:,3) == j);
                             find(mesh.incidence(:,4) == j); find(mesh.incidence(:,5) == j)];

                % Sorting conectel
                connectel = sort(connectel);

                % Allocating conectel
                nodal_connectivity(j,1) = j;                     
                nodal_connectivity(j,(2:length(connectel)+1)) = connectel;
            end

        end %end BuildNodalConnectivity
        
        %% Build matrix with element neighbours
        function element_neighbours = FindElementNeighbours(mesh)
            
            disp(['         Finding element neighbours.'])

            % Initializing matrix
            element_neighbours = zeros(size(mesh.incidence,1),5);
            element_neighbours(:,1) = [1:size(mesh.incidence,1)]';

            % Find and store element neighbours
            for i = 1:size(mesh.incidence,1)

                % Counter
                aux_element_neighbours = 0;

                for j = 2:5  % (j = element's node)

                    node = mesh.incidence(i,j);

                    for k = 2:5  % (k = node's connectivity)

                        %Se a conectividade de do nó j é diferente de i e de 0 e o elemento i existe
                        if ((mesh.nodal_connectivity(node,k) ~= i) && (mesh.nodal_connectivity(node,k) ~= 0));

                            % Shared nodes counter
                            nodeshared = 0;

                            for l = 2:5  %(l = ith element's incidence)
                                for m = 2:5  %(m = neighbour element's incidence)
                                    % If there is a shared node
                                    if (mesh.incidence(i,l) == mesh.incidence(mesh.nodal_connectivity(node,k),m));
                                        nodeshared = nodeshared+1;
                                    end
                                end
                            end

                            % If elements share 2 nodes, they are neighbours
                            if (nodeshared == 2);
                                
                                % Checking if neighbour was already stored
                                if (length(find(element_neighbours(i,2:5) == mesh.nodal_connectivity(node,k))) ~= 1)
                                    aux_element_neighbours = aux_element_neighbours+1;
                                    % Store neighbour
                                    element_neighbours(i,aux_element_neighbours+1) = mesh.nodal_connectivity(node,k);
                                end

                            end

                        end
                    end
                end

            end 
            
        end % end FindElementNeighbours
        
        %% Build matrix to relate the element number and it is position in a
        % rectangular matrix (m x n)
        function plot_matrix = BuildPlotMatrix(mesh)
            
            disp(['         Building plot matrix.'])
            
            % Internal auxiliary matrix
            centerplot = [mesh.centroids];

            % Maximum of centroids
            max_x = max(max(centerplot(:,1)));
            max_y = max(max(centerplot(:,2)));

            % Elements with maximum of centroid
            elemax_x = find(centerplot(:,1) == max_x);
            elemax_y = find(centerplot(:,2) == max_y);

            % Finding height
            height =  length(find(centerplot(:,1) == centerplot(elemax_y(1),1)));
            
            % Building matrix
            for i = height:-1:1

                % Finding a line of elements in x
                line = find(centerplot(:,2) == min(centerplot(:,2)));

                % Centroid values in y
                yc(i,1) = min(centerplot(:,2));

                % Sorting line of elements according to y coord
                [ordenado,ordem] = sort(centerplot(line,1));

                % Inserting line of elements in plot_matrix
                plot_matrix(i,1:length(line)) = line(ordem);

                % Avoiding same line to be selected two times
                centerplot(line,2) = 2*max([max_x; max_y]);

            end
            
        end % end BuildPlotMatrix
    
    end % end methods
end