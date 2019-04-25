% Function to plot structural displacements

function PlotStructuralDisplacements(fea,scale)

    disp([' '])
    disp(['         Plotting structural displacements.'])

    hold on

    % Initializing deformed coordinates and displacements
    Dcoord = fea.mesh.coordinates;
    Ux = zeros(size(fea.mesh.coordinates,1),1);
    Uy = zeros(size(fea.mesh.coordinates,1),1);

    for i = 1:size(fea.mesh.coordinates,1)

        % Identifying nodal displacements
        % X axis
        Dcoord(i,1) = Dcoord(i,1)+scale*fea.U(fea.ID(1,i));
        Ux(i) = fea.U(fea.ID(1,i));
        % Y axis
        Dcoord(i,2) = Dcoord(i,2)+scale*fea.U(fea.ID(2,i));
        Uy(i) = fea.U(fea.ID(2,i));

    end

    % Displacement vector magnitude
    Gu = sqrt( Ux.^2 + Uy.^2 );

    for i = 1:size(fea.mesh.incidence,1)

        inci = fea.mesh.incidence;

        % Plotting displacements
        x = [ Dcoord(inci(i,2),1) Dcoord(inci(i,3),1) Dcoord(inci(i,4),1) Dcoord(inci(i,5),1) ];
        y = [ Dcoord(inci(i,2),2) Dcoord(inci(i,3),2) Dcoord(inci(i,4),2) Dcoord(inci(i,5),2) ];
        c = [ Gu(inci(i,2)) Gu(inci(i,3)) Gu(inci(i,4)) Gu(inci(i,5)) ];
        fill(x,y,c,'LineStyle','none')

    end

    axis off; axis equal;

end
            