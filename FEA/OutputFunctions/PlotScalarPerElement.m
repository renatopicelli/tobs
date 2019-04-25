% Function to plot one scalar value per element

function PlotScalarPerElement(fea,values,exp)

    % Defining matplot
    matplot = real(values.^(exp));
    matrixplot2 = fea.mesh.plot_matrix;
    % Verifying if there are zeros in matrizplot
    if (length(find(matrixplot2 == 0)) >= 1)
        % Defining void for inexistent elements in the plot
        matplot(size(fea.mesh.incidence,1)+1) = 0;
        matrixplot2(find(matrixplot2 == 0)) = size(fea.mesh.incidence,1)+1;
    end

    % Plotting
    colormap(jet); imagesc(matplot(matrixplot2));
    axis equal; axis off;

end