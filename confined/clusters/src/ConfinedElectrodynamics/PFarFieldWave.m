function [potential, potentialGradient] = ...
    PFarFieldWave(Parameter, Mesh, Particle)
%--------------------------------------------------------------------------
% pFarFieldWave - computes the 'global' potential and potential gradient.
%
% [potential, potentialGradient] = pFarFieldWave(Parameter, Particle)
% computes the global potential (due to Gaussian- filtered charges) and the
% potential gradient (due to Gaussian-filtered dipoles) at particle 
% centers.
%
% potential: (nParticle x 1) array of self-potentials.
%
% potentialGradient: (3 x nParticle) array of self-potential gradients.
%--------------------------------------------------------------------------

Particle.position = Particle.position';

% Convert to Pseudo domain
Particle.charge = [Particle.charge;-Particle.charge];
Particle.field = [Particle.field Particle.field];
Particle.potential = [Particle.potential;-Particle.potential];
Particle.farFieldCharge = [Particle.farFieldCharge; -Particle.farFieldCharge];
Particle.farFieldDipole = [Particle.farFieldDipole Particle.farFieldDipole];

[Particle, Parameter] = mapTo2L3(Particle, Parameter);


%% Distribute charge and dipoles onto mesh points
rho = distributeCharge(Parameter, Particle, Mesh);

%% Solve for potential and potential gradient at mesh points
% Fourier transform
RHO = ifftn(rho);

% Fourier coefficients for potential
PHI = Mesh.operator .* RHO;
idx = Mesh.Ksquared == 0;
PHI(idx) = 0;

% Potential at mesh points
phi = real(fftn(PHI)) / (Mesh.spacing(1)*Mesh.spacing(2)*Mesh.spacing(3)); % nMesh x nMesh x nMesh
% Particle.position
% Parameter.domainLength
% Particle.farFieldDipole
% Particle.farFieldCharge
%% Interpolate potential and potential gradient
[potential, potentialGradient] = interpolatePotential(Parameter, Particle, Mesh, phi);
% Parameter.splittingAlpha

potential = potential(1:Particle.nParticle/2);
potentialGradient = potentialGradient(:,1:Particle.nParticle/2);

%Plot Results
% tmp = zeros(Mesh.nMesh(1),Mesh.nMesh(3));
% for i = 1:Mesh.nMesh(2)
%     tmp(:,:) = phi(:,i,:);
%     imagesc(Mesh.spacing(1),Mesh.spacing(3),tmp'); axis xy; axis equal; colorbar;
%     pause;
% end
% disp('XXXXXXXXXXXXXXXXXXXXXX')
% pause;

%--------------------------------------------------------------------------
function chargeDensity = distributeCharge(Parameter, Particle, Mesh)
% Distributes far field charge and dipole onto the mesh points

% Charge Density
chargeDensity = zeros(Mesh.nMesh(1), Mesh.nMesh(2), Mesh.nMesh(3));

% Find indices of nearest grid points
gridIndex = zeros(Particle.nParticle, 3);
gridIndex(:,1) = mod(round( Particle.position(:,1) / Mesh.spacing(1) ),...
    Mesh.nMesh(1)) + 1;
gridIndex(:,2) = mod(round( Particle.position(:,2) / Mesh.spacing(2) ),...
    Mesh.nMesh(2)) + 1;
gridIndex(:,3) = mod(round( Particle.position(:,3) / Mesh.spacing(3) ),...
    Mesh.nMesh(3)) + 1;

% Find positions of nearest grid points
gridPosition = zeros(Particle.nParticle, 3);
gridPosition(:,1) = Mesh.x1( gridIndex(:,1) );
gridPosition(:,2) = Mesh.x2( gridIndex(:,2) );
gridPosition(:,3) = Mesh.x3( gridIndex(:,3) );

% Fractional displacement between particles and proximal grid point
fractionalDisplacement = Particle.position - gridPosition;
fractionalDisplacement(:,1) = fractionalDisplacement(:,1) ... 
    - round(fractionalDisplacement(:,1) / Parameter.domainLength(1)) ...
    * Parameter.domainLength(1);
fractionalDisplacement(:,1) = fractionalDisplacement(:,1) / Mesh.spacing(1);
fractionalDisplacement(:,2) = fractionalDisplacement(:,2) ... 
    - round(fractionalDisplacement(:,2) / Parameter.domainLength(2)) ...
    * Parameter.domainLength(2);
fractionalDisplacement(:,2) = fractionalDisplacement(:,2) / Mesh.spacing(2);
fractionalDisplacement(:,3) = fractionalDisplacement(:,3) ... 
    - round(fractionalDisplacement(:,3) / Parameter.domainLength(3)) ...
    * Parameter.domainLength(3);
fractionalDisplacement(:,3) = fractionalDisplacement(:,3) / Mesh.spacing(3);


% preallocate
wp = zeros(3,3); % NGRID x NDIM
wq = zeros(3,3); % NGRID x NDIM

for iParticle = 1:Particle.nParticle % loop over particles
    % Weighting Function (charge)
    wq(1,:) = 0.5 * (fractionalDisplacement(iParticle,:).^2 ...
        - fractionalDisplacement(iParticle,:));
    wq(2,:) = (1 - fractionalDisplacement(iParticle,:).^2);
    wq(3,:) = 0.5 * (fractionalDisplacement(iParticle,:).^2 ...
        + fractionalDisplacement(iParticle,:));
    
    % Weighting Function (dipole)
    wp(1,:) = 0.5 * (-1 + 2 * fractionalDisplacement(iParticle,:));
    wp(2,:) = -2 * fractionalDisplacement(iParticle,:);
    wp(3,:) = 0.5 * (1 + 2 * fractionalDisplacement(iParticle,:));  
    wp(1,:) = wp(1,:) / Mesh.spacing(1);
    wp(2,:) = wp(2,:) / Mesh.spacing(2);
    wp(3,:) = wp(3,:) / Mesh.spacing(3);
    
    % Distribute Charge
    n = [-1:1];
    i1 = mod( n + gridIndex(iParticle,1) - 1, Mesh.nMesh(1)) + 1;
    i2 = mod( n + gridIndex(iParticle,2) - 1, Mesh.nMesh(2)) + 1;
    i3 = mod( n + gridIndex(iParticle,3) - 1, Mesh.nMesh(3)) + 1;
 
    for i = 1:3
        for j = 1:3
            for k = 1:3
                chargeDensity(i1(i),i2(j),i3(k)) = chargeDensity(i1(i),i2(j),i3(k)) + ...
                    Particle.farFieldCharge(iParticle) * (wq(i,1) * wq(j,2) * wq(k,3)) + ...
                    Particle.farFieldDipole(1,iParticle) * (wp(i,1) * wq(j,2) * wq(k,3))  + ...
                    Particle.farFieldDipole(2,iParticle) * (wq(i,1) * wp(j,2) * wq(k,3))  + ...
                    Particle.farFieldDipole(3,iParticle) * (wq(i,1) * wq(j,2) * wp(k,3)) ;
            end
        end
    end
end


%--------------------------------------------------------------------------
function [potential, potentialGradient] = ...
    interpolatePotential(Parameter, Particle, Mesh, phi)
% Interpolate the potential and potential gradient at particle centers

potential = zeros(Particle.nParticle,1);
potentialGradient = zeros(3,Particle.nParticle);

% Find indices of nearest grid points
gridIndex = zeros(Particle.nParticle, 3);
gridIndex(:,1) = mod(round( Particle.position(:,1) / Mesh.spacing(1) ),...
    Mesh.nMesh(1)) + 1;
gridIndex(:,2) = mod(round( Particle.position(:,2) / Mesh.spacing(2) ),...
    Mesh.nMesh(2)) + 1;
gridIndex(:,3) = mod(round( Particle.position(:,3) / Mesh.spacing(3) ),...
    Mesh.nMesh(3)) + 1;

% Find positions of nearest grid points
gridPosition = zeros(Particle.nParticle, 3);
gridPosition(:,1) = Mesh.x1( gridIndex(:,1) );
gridPosition(:,2) = Mesh.x2( gridIndex(:,2) );
gridPosition(:,3) = Mesh.x3( gridIndex(:,3) );

% Fractional displacement between particles and proximal grid point
fractionalDisplacement = Particle.position - gridPosition;
fractionalDisplacement(:,1) = fractionalDisplacement(:,1) ... 
    - round(fractionalDisplacement(:,1) / Parameter.domainLength(1)) ...
    * Parameter.domainLength(1);
fractionalDisplacement(:,1) = fractionalDisplacement(:,1) / Mesh.spacing(1);
fractionalDisplacement(:,2) = fractionalDisplacement(:,2) ... 
    - round(fractionalDisplacement(:,2) / Parameter.domainLength(2)) ...
    * Parameter.domainLength(2);
fractionalDisplacement(:,2) = fractionalDisplacement(:,2) / Mesh.spacing(2);
fractionalDisplacement(:,3) = fractionalDisplacement(:,3) ... 
    - round(fractionalDisplacement(:,3) / Parameter.domainLength(3)) ...
    * Parameter.domainLength(3);
fractionalDisplacement(:,3) = fractionalDisplacement(:,3) / Mesh.spacing(3);


% preallocate
w = zeros(3,3); % NGRID x NDIM
dw = zeros(3,3); % NGRID x NDIM

for iParticle = 1:Particle.nParticle % loop over particles
    % Weighting Function
    w(1,:) = 0.5*(fractionalDisplacement(iParticle,:).^2 ...
        - fractionalDisplacement(iParticle,:));
    w(2,:) = (1 - fractionalDisplacement(iParticle,:).^2);
    w(3,:) = 0.5*(fractionalDisplacement(iParticle,:).^2 ...
        + fractionalDisplacement(iParticle,:));    
    
    % Weighting Function
    dw(1,:) = 0.5*(-1 + 2*fractionalDisplacement(iParticle,:));
    dw(2,:) = -2*fractionalDisplacement(iParticle,:);
    dw(3,:) = 0.5*(1 + 2*fractionalDisplacement(iParticle,:));
    dw(1,:) = dw(1,:) / Mesh.spacing(1);
    dw(2,:) = dw(2,:) / Mesh.spacing(2);
    dw(3,:) = dw(3,:) / Mesh.spacing(3);
    
    % Interpolate
    n = [-1:1];
    i1 = mod( n + gridIndex(iParticle,1) - 1, Mesh.nMesh(1)) + 1;
    i2 = mod( n + gridIndex(iParticle,2) - 1, Mesh.nMesh(2)) + 1;
    i3 = mod( n + gridIndex(iParticle,3) - 1, Mesh.nMesh(3)) + 1;
    
    for i = 1:3
        for j = 1:3
            for k = 1:3
                potential(iParticle) = potential(iParticle) + ...
                    phi(i1(i),i2(j),i3(k)) * (w(i,1)*w(j,2)*w(k,3));
                potentialGradient(1,iParticle) = potentialGradient(1,iParticle) + ...
                    phi(i1(i),i2(j),i3(k)) * (dw(i,1)*w(j,2)*w(k,3));
                potentialGradient(2,iParticle) = potentialGradient(2,iParticle) + ...
                    phi(i1(i),i2(j),i3(k)) * (w(i,1)*dw(j,2)*w(k,3));
                potentialGradient(3,iParticle) = potentialGradient(3,iParticle) + ...
                    phi(i1(i),i2(j),i3(k)) * (w(i,1)*w(j,2)*dw(k,3));
            end
        end
    end
end
