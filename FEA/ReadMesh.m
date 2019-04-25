% MATLAB code to read the mesh (generated on ANSYS) of FEM
function [coordinates,incidence,dirichlet_boundary,neumann_boundary] = ReadMesh(example)

    % Reading mesh of a .dat file generated from ANSYS
    fid = fopen(strcat(example,'.dat'));
    
    disp(['         Reading ',example,'.dat mesh ...'])

    %-------------------------------------------------------------------------%
    % Number of nodes and elements

    % Getting line from .dat
    tline = fgetl(fid);
    % Numer of nodes
    nnos = str2double(tline(5:end));

    % Getting line from .dat
    tline = fgetl(fid);
    % Numer of elements
    nel = str2double(tline(5:end));

    %-------------------------------------------------------------------------%
    % Coordinates matrix

    disp(['           ... nodal coordinates,'])

    % Initializing coord matrix
    coord = zeros(nnos,4);

    % Getting line from .dat until coord matrix
    while (strcmp(tline(1:3), 'noh') == 0)
        % Getting line from .dat
        tline = fgetl(fid);
    end

    % Finding commas' positions
    aux = 0;
    for i = 1:length(tline)
        if strcmp(tline(i), ',')
            aux = aux+1;
            pos(aux) = i;
        end
    end

    % Declaring first line of coord
    coord(1,1) = str2double(tline(5:(pos(2)-1)));
    coord(1,2) = str2double(tline((pos(2)+1):(pos(3)-1)));
    coord(1,3) = str2double(tline((pos(3)+1):(pos(4)-1)));
    coord(1,4) = str2double(tline((pos(4)+1):end));

    % Reading the rest of coord matrix
    for i = 2:nnos
        % Getting line from .dat
        tline = fgetl(fid);
        % Declaring lines of coord
        coord(i,1) = str2double(tline(5:(pos(2)-1)));
        coord(i,2) = str2double(tline((pos(2)+1):(pos(3)-1)));
        coord(i,3) = str2double(tline((pos(3)+1):(pos(4)-1)));
        coord(i,4) = str2double(tline((pos(4)+1):end));
    end

    % Adjusting coordinate matrix
    coord = coord(:,2:4);

    %-------------------------------------------------------------------------%
    % Incidence matrix

    disp(['           ... element incidence,'])

    % Initializing inci matrix
    inci = zeros(nel,8);

    % Getting line from .dat
    tline = fgetl(fid);

    % Finding commas' positions
    aux = 0;
    clear pos
    for i = 1:length(tline)
        if strcmp(tline(i), ',')
            aux = aux+1;
            pos(aux) = i;
        end
    end

    % Declaring first line of inci
    inci(1,1) = str2double(tline(5:(pos(2)-1)));
    inci(1,2) = str2double(tline((pos(2)+1):(pos(3)-1)));
    inci(1,3) = str2double(tline((pos(3)+1):(pos(4)-1)));
    inci(1,4) = str2double(tline((pos(4)+1):(pos(5)-1)));
    inci(1,5) = str2double(tline((pos(5)+1):(pos(6)-1)));
    inci(1,6) = str2double(tline((pos(6)+1):(pos(7)-1)));
    inci(1,7) = str2double(tline((pos(7)+1):(pos(8)-1)));
    inci(1,8) = str2double(tline((pos(8)+1):(pos(9)-1)));

    % Reading the rest of inci matrix
    for i = 2:nel
        % Getting line from .dat
        tline = fgetl(fid);
        % Declaring lines of inci
        inci(i,1) = str2double(tline(5:(pos(2)-1)));
        inci(i,2) = str2double(tline((pos(2)+1):(pos(3)-1)));
        inci(i,3) = str2double(tline((pos(3)+1):(pos(4)-1)));
        inci(i,4) = str2double(tline((pos(4)+1):(pos(5)-1)));
        inci(i,5) = str2double(tline((pos(5)+1):(pos(6)-1)));
        inci(i,6) = str2double(tline((pos(6)+1):(pos(7)-1)));
        inci(i,7) = str2double(tline((pos(7)+1):(pos(8)-1)));
        inci(i,8) = str2double(tline((pos(8)+1):(pos(9)-1)));
    end

    % Adjusting incidence matrix
    inci = [ ones(size(inci,1),1) inci(:,5:8) ];

    %-------------------------------------------------------------------------%
    % Boundary conditions matrix

    disp(['           ... boundary conditions.'])

    % Getting line from .dat
    tline = fgetl(fid);
    % Getting line from .dat
    tline = fgetl(fid);

    % Finding commas' positions
    aux = 0;
    clear pos
    for i = 1:length(tline)
        if strcmp(tline(i), ',')
            aux = aux+1;
            pos(aux) = i;
        end
    end

    % Declaring first line of cont
    cont(1,1) = str2double(tline(5:(pos(2)-1)));
    if strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'UX')
        cont(1,2) = 1;
    elseif strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'UY')
        cont(1,2) = 2;
    elseif strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'UZ')
        cont(1,2) = 3;
    elseif strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'PR')
        cont(1,2) = 4;
    end
    cont(1,3) = str2double(tline((pos(3)+1):end));

    % Reading the rest of cont matrix
    aux = 1;
    % Getting line from .dat
    tline = fgetl(fid);
    while (strcmp(tline(1:3), 'end') == 0)
        aux = aux+1;
        % Declaring lines of cont
        cont(aux,1) = str2double(tline(5:(pos(2)-1)));
        if strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'UX')
            cont(aux,2) = 1;
        elseif strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'UY')
            cont(aux,2) = 2;
        elseif strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'UZ')
            cont(aux,2) = 3;
        elseif strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'PR')
            cont(aux,2) = 4;
        end
        cont(aux,3) = str2double(tline((pos(3)+1):end));
        % Getting line from .dat
        tline = fgetl(fid);
    end

    % Size of boundary conditions matrix
    ncont = size(cont,1);

    %-------------------------------------------------------------------------%
    % Load vector

    % Getting line from .dat
    tline = fgetl(fid);

    % Finding commas' positions
    aux = 0;
    clear pos
    for i = 1:length(tline)
        if strcmp(tline(i), ',')
            aux = aux+1;
            pos(aux) = i;
        end
    end

    % Declaring first line of forca
    loads(1,1) = str2double(tline(5:(pos(2)-1)));
    if strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'FX')
        loads(1,2) = 1;
    elseif strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'FY')
        loads(1,2) = 2;
    elseif strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'FZ')
        loads(1,2) = 3;
    elseif strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'DP')
        loads(1,2) = 4;
    end
    loads(1,3) = str2double(tline((pos(3)+1):end));

    % Reading the rest of cont matrix
    aux = 1;
    % Getting line from .dat
    tline = fgetl(fid);
    while (strcmp(tline(1:3), 'end') == 0)
        aux = aux+1;
        % Declaring lines of cont
        loads(aux,1) = str2double(tline(5:(pos(2)-1)));
        if strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'FX')
            loads(aux,2) = 1;
        elseif strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'FY')
            loads(aux,2) = 2;
        elseif strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'FZ')
            loads(aux,2) = 3;
        elseif strcmp(tline([(pos(3)-2) (pos(3)-1)]), 'DP')
            loads(aux,2) = 4;
        end
        loads(aux,3) = str2double(tline((pos(3)+1):end));
        % Getting line from .dat
        tline = fgetl(fid);
    end

    % Length of load vector
    nloads = size(loads,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output matrices for OOP code
    coordinates = coord;
    incidence = inci;
    dirichlet_boundary = cont;
    neumann_boundary = loads;

end