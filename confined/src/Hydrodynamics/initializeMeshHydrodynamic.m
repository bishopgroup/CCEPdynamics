function Mesh = initializeMeshHydrodynamic(Parameter, Mesh)
%--------------------------------------------------------------------------
% initializeMesh - initializes particle mesh.
%
% Mesh = initializeMesh(Parameter, Mesh)
%--------------------------------------------------------------------------

%% mesh spacing
Mesh.nMesh = 64*ones(3,1); % round to power of two
Mesh.spacing =  Parameter.domainLength ./ Mesh.nMesh ;


%% real-space mesh points
Mesh.x1 = Mesh.spacing(1) * [0 : Mesh.nMesh(1) - 1]';
Mesh.x2 = Mesh.spacing(2) * [0 : Mesh.nMesh(2) - 1]';
Mesh.x3 = Mesh.spacing(3) * [0 : Mesh.nMesh(3) - 1]';
[Mesh.X2, Mesh.X1, Mesh.X3] = meshgrid(Mesh.x2, Mesh.x1, Mesh.x3);


%% fourier-space wave vectors
% [0, 1, ..., N/2, -N/2, ..., -1]
Mesh.k1 = [[0:Mesh.nMesh(1)/2], [-Mesh.nMesh(1)/2+1:-1]]' / Parameter.domainLength(1);  
Mesh.k2 = [[0:Mesh.nMesh(2)/2], [-Mesh.nMesh(2)/2+1:-1]]' / Parameter.domainLength(2);
Mesh.k3 = [[0:Mesh.nMesh(3)/2], [-Mesh.nMesh(3)/2+1:-1]]' / Parameter.domainLength(3);
[Mesh.K2, Mesh.K1, Mesh.K3] = meshgrid(Mesh.k2, Mesh.k1, Mesh.k3);
Mesh.Ksquared = Mesh.K1.^2 + Mesh.K2.^2 + Mesh.K3.^2;
Mesh.K = sqrt(Mesh.Ksquared);
Mesh.K1Unit = Mesh.K1 ./ Mesh.K;
Mesh.K2Unit = Mesh.K2 ./ Mesh.K;
Mesh.K3Unit = Mesh.K3 ./ Mesh.K;

%% linear operator use in solving the Stokes equation
Mesh.hydrodynamicOperator = 1.5 * (1 + pi * Parameter.splittingAlpha * Mesh.Ksquared) ...
    .* exp(-pi * Parameter.splittingAlpha * Mesh.Ksquared) ...
    ./ (pi * Mesh.Ksquared);

%find(Mesh.Ksquared == 0)



