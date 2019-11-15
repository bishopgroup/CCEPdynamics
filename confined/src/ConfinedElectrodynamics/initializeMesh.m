function [Mesh] = initializeMesh(Parameter)


%% mesh spacing
Parameter.nMesh = 2.^5;
Parameter.domainLength(3) = 2*Parameter.domainLength(3);% changed to pseudo domain
Mesh.nMesh = [2 2 2]'.*Parameter.nMesh;
Mesh.spacing = Parameter.domainLength' ./ Mesh.nMesh;


%% real-space mesh points
Mesh.x1 = Mesh.spacing(1) * [0 : Mesh.nMesh(1) - 1]';
Mesh.x2 = Mesh.spacing(2) * [0 : Mesh.nMesh(2) - 1]';
Mesh.x3 = Mesh.spacing(3) * [0 : Mesh.nMesh(3) - 1]';
[Mesh.X1, Mesh.X2, Mesh.X3] =  ndgrid(Mesh.x1, Mesh.x2, Mesh.x3);


%% fourier-space wave vectors
% [0, 1, ..., N/2, -N/2, ..., -1]
Mesh.k1 = [0:Mesh.nMesh(1)/2, -Mesh.nMesh(1)/2+1:-1]' ...
    / Parameter.domainLength(1);
Mesh.k2 = [0:Mesh.nMesh(2)/2, -Mesh.nMesh(2)/2+1:-1]' ...
    / Parameter.domainLength(2);
Mesh.k3 = [0:Mesh.nMesh(3)/2, -Mesh.nMesh(3)/2+1:-1]' ...
    / Parameter.domainLength(3);
[Mesh.K1, Mesh.K2, Mesh.K3] = ndgrid(Mesh.k1, Mesh.k2, Mesh.k3);
Mesh.Ksquared = Mesh.K1.^2 + Mesh.K2.^2 + Mesh.K3.^2;


%% linear operator use in solving the Poisson equation
Mesh.operator = exp(-pi * Parameter.splittingAlpha * Mesh.Ksquared)...
    ./ (pi * Mesh.Ksquared);
