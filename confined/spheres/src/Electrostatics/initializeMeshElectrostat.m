function [Mesh]= initializeMeshElectrostat(Parameter)
global I e3 e3e3 SSH_X11A_ln SSH_X12A_ln SSH_Y11A SSH_Y12A SSH_Y11B ...
    SSH_Y12B SSH_X11C SSH_X12C SSH_Y11C SSH_Y12C ...
    KroneckerDelta Signature SPH_XA_ln SPH_YA SPH_YB SPH_XC SPH_YC 

%% Tensors
I = eye(3); % Identity Tensor
e3 = zeros(3,1); e3(3) = 1;
e3e3 = zeros(3); e3e3(3,3) = 1;


%% mesh spacing
Mesh.nMesh = [round(Parameter.domainLength(1)/Parameter.domainLength(3)),
              round(Parameter.domainLength(2)/Parameter.domainLength(3)),
              2]' .* Parameter.nMesh;
Mesh.spacing = [Parameter.domainLength(1) / Mesh.nMesh(1),
                Parameter.domainLength(2) / Mesh.nMesh(2),
              2*Parameter.domainLength(3) / Mesh.nMesh(3)];
Mesh.spacing

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


%% Sphere-Plane Resistance
data = load('SpherePlaneResistance.txt');
SPH_XA_ln = griddedInterpolant(log(data(:,1)),log(data(:,2)));
SPH_YA = griddedInterpolant(log(data(:,1)),data(:,3));
SPH_YB = griddedInterpolant(log(data(:,1)),data(:,4));
SPH_XC = griddedInterpolant(log(data(:,1)),data(:,5));
SPH_YC = griddedInterpolant(log(data(:,1)),data(:,6));

%% Two Sphere Resistance
data = load('SphereSphereResistance.txt');
SSH_X11A_ln = griddedInterpolant(log(data(:,1)-2),log(data(:,2)));
SSH_X12A_ln = griddedInterpolant(log(data(:,1)-2),log(-data(:,3)));
SSH_Y11A = griddedInterpolant(log(data(:,1)-2),data(:,4));
SSH_Y12A = griddedInterpolant(log(data(:,1)-2),data(:,5));
SSH_Y11B = griddedInterpolant(log(data(:,1)-2),data(:,6));
SSH_Y12B = griddedInterpolant(log(data(:,1)-2),data(:,7));
SSH_X11C = griddedInterpolant(log(data(:,1)-2),data(:,8));
SSH_X12C = griddedInterpolant(log(data(:,1)-2),data(:,9));
SSH_Y11C = griddedInterpolant(log(data(:,1)-2),data(:,10));
SSH_Y12C = griddedInterpolant(log(data(:,1)-2),data(:,11));

KroneckerDelta=eye(3);

Signature=zeros(3,3,3);
Signature(1,2,3)=1;
Signature(2,3,1)=1;
Signature(3,1,2)=1;
Signature(2,1,3)=-1;
Signature(3,2,1)=-1;
Signature(1,3,2)=-1;
