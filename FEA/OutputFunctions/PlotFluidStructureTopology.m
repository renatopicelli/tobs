% Function to plot one scalar value per element

function PlotFluidStructureTopology(fea)

    % Defining matplot
    matplot = fea.element_types;
    matrixplot2 = fea.mesh.plot_matrix;
    % Verifying if there are zeros in matrizplot
    if (length(find(matrixplot2 == 0)) >= 1)
        % Defining void for inexistent elements in the plot
        matplot(size(fea.mesh.incidence,1)+1) = 0;
        matrixplot2(find(matrixplot2 == 0)) = size(fea.mesh.incidence,1)+1;
    end

    % Plotting topology with imagesc
    imagesc(matplot(matrixplot2)); axis off; axis equal;
    % Colormap definition
    map = [   1    1    1;
              0    0    0;
           0.58 0.85 0.85;];
    % Colormap when there are only solid and fluid elements
    if (nnz(find(matplot == 0)) == 0)
        map = [  0    0    0;
              0.58 0.85 0.85;];
    end
    colormap(map);

end