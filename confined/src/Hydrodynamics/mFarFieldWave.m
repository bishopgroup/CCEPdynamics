function [velocity, velocityGradient, velocityLaplacian, ...
    velocityLaplacianGradient] = mFarFieldWave(Parameter, Particle, Mesh)
%--------------------------------------------------------------------------
% mFarFieldWave - computes the 'global' velocity and velocity gradient.
%
% [velocity, velocityGradient] = mFarFieldWave(Parameter, Particle)
% computes the global velocity (due to Gaussian-filtered forces) and the
% velocity gradient (due to Gaussian-filtered dipoles) at particle 
% centers.
%
% velocity: (3 x nParticle) array of self-velocity.
%
% velocityGradient: (9 x nParticle) array of self-velocity gradients.
%--------------------------------------------------------------------------

%% Distribute force and dipoles onto mesh points
[gridForceX,gridForceY,gridForceZ] = distributeForce(Parameter, Particle, Mesh);

%% Solve for potential at mesh points
% Fourier transform
gridForceWaveX = ifftn(gridForceX);
gridForceWaveY = ifftn(gridForceY);
gridForceWaveZ = ifftn(gridForceZ);

% Fourier coefficients for velocity
gridVelocityWaveX = Mesh.hydrodynamicOperator .* ((1 - Mesh.K1Unit .* Mesh.K1Unit) .* gridForceWaveX ...
    - Mesh.K1Unit .* Mesh.K2Unit .* gridForceWaveY ...
    - Mesh.K1Unit .* Mesh.K3Unit .* gridForceWaveZ);
gridVelocityWaveY = Mesh.hydrodynamicOperator .* (-Mesh.K2Unit .* Mesh.K1Unit .* gridForceWaveX ...
    + (1 - Mesh.K2Unit .* Mesh.K2Unit) .* gridForceWaveY ...
    - Mesh.K2Unit .* Mesh.K3Unit .* gridForceWaveZ);
gridVelocityWaveZ = Mesh.hydrodynamicOperator .* (-Mesh.K3Unit .* Mesh.K1Unit .* gridForceWaveX ...
    - Mesh.K3Unit .* Mesh.K2Unit .* gridForceWaveY ...
    + (1 - Mesh.K3Unit .* Mesh.K3Unit) .* gridForceWaveZ);

% exclude k = 0 contribution
gridVelocityWaveX(1,1,1) = 0;
gridVelocityWaveY(1,1,1) = 0;
gridVelocityWaveZ(1,1,1) = 0;

% Inverse Fourier transform
gridVelocityX = real(fftn(gridVelocityWaveX)) / (Mesh.spacing(1)*Mesh.spacing(2)*Mesh.spacing(3));
gridVelocityY = real(fftn(gridVelocityWaveY)) / (Mesh.spacing(1)*Mesh.spacing(2)*Mesh.spacing(3));
gridVelocityZ = real(fftn(gridVelocityWaveZ)) / (Mesh.spacing(1)*Mesh.spacing(2)*Mesh.spacing(3));


%% Interpolate potential and potential gradient
[velocity, velocityGradient, velocityLaplacian, ...
    velocityLaplacianGradient] = interpolateVelocity(Parameter, Particle, ...
    Mesh, gridVelocityX, gridVelocityY, gridVelocityZ);

%% Plot Results
% tmp = zeros(Mesh.nMesh,Mesh.nMesh);
% for i = 1:Mesh.nMesh
%     tmp(:,:) = gridVelocityZ(:,i,:);
%     imagesc(Mesh.x1,Mesh.x3,tmp'); axis xy; axis equal; colorbar;
%     pause(0.1);
% end


function [gridForceX,gridForceY,gridForceZ] = ...
    distributeForce(Parameter, Particle, Mesh)
% Distributes far field force and dipole onto the mesh points

nGrid = 4; % up to third-order moments

%% Preallocate array for force at grid points
gridForceX = zeros(Mesh.nMesh(1), Mesh.nMesh(2), Mesh.nMesh(3));
gridForceY = zeros(Mesh.nMesh(1), Mesh.nMesh(2), Mesh.nMesh(3));
gridForceZ = zeros(Mesh.nMesh(1), Mesh.nMesh(2), Mesh.nMesh(3));


%% Grid points
% Find indices of nearest grid points (to the left for nGrid = 4)
gridIndex = zeros(Particle.nParticle, 3);
gridIndex(:,1) = mod(floor( Particle.position(:,1) / Mesh.spacing(1) ),...
    Mesh.nMesh(1)) + 1;
gridIndex(:,2) = mod(floor( Particle.position(:,2) / Mesh.spacing(2) ),...
    Mesh.nMesh(2)) + 1;
gridIndex(:,3) = mod(floor( Particle.position(:,3) / Mesh.spacing(3) ),...
    Mesh.nMesh(3)) + 1;

% Find positions of nearest grid points
gridPosition = zeros(Particle.nParticle, 3);
gridPosition(:,1) = Mesh.x1( gridIndex(:,1) );
gridPosition(:,2) = Mesh.x2( gridIndex(:,2) );
gridPosition(:,3) = Mesh.x3( gridIndex(:,3) );

% Fractional displacement between particles and proximal grid point
delta = Particle.position - gridPosition;

delta(:,1) = delta(:,1) - round(delta(:,1) / Parameter.domainLength(1)) ...
    * Parameter.domainLength(1);
delta(:,1) = delta(:,1) / Mesh.spacing(1);
delta(:,2) = delta(:,2) - round(delta(:,2) / Parameter.domainLength(2)) ...
    * Parameter.domainLength(2);
delta(:,2) = delta(:,2) / Mesh.spacing(2);
delta(:,3) = delta(:,3) - round(delta(:,3) / Parameter.domainLength(3)) ...
    * Parameter.domainLength(3);
delta(:,3) = delta(:,3) / Mesh.spacing(3);

% Indices of grid points around each particle
n = [-1:2]; % for nGrid = 4
iX = mod( bsxfun(@plus, n , gridIndex(:,1)) - 1, Mesh.nMesh(1)) + 1; % nParticle x nGrid
iY = mod( bsxfun(@plus, n , gridIndex(:,2)) - 1, Mesh.nMesh(2)) + 1;
iZ = mod( bsxfun(@plus, n , gridIndex(:,3)) - 1, Mesh.nMesh(3)) + 1;


%% Weighting functions
% Zeroth moment
w0 = zeros(Particle.nParticle,3,nGrid); % nParticle x nDim x nGrid
w0(:,:,1) = -(1/6) * (delta - 2) .* (delta - 1) .* delta;
w0(:,:,2) = 0.5 * (delta - 2) .* (delta - 1) .* (1 + delta);
w0(:,:,3) = 0.5 * delta .* (2 + delta - delta.^2);
w0(:,:,4) = (1/6) * delta .* (-1 + delta.^2);

% First moment
w1 = zeros(Particle.nParticle,3,nGrid); % nParticle x nDim x nGrid
w1(:,:,1) = -(1/3) + delta - 0.5*delta.^2;
w1(:,:,2) = 0.5 * (-1 + delta .* (-4 + 3*delta));
w1(:,:,3) = 1 + delta - 1.5*delta.^2;
w1(:,:,4) = -(1/6) + 0.5*delta.^2;
w1 = w1 / Mesh.spacing(3);

% Second moment
w2 = zeros(Particle.nParticle,3,nGrid); % nParticle x nDim x nGrid
w2(:,:,1) = 0.5 * (1 - delta);
w2(:,:,2) = -1 + 1.5*delta;
w2(:,:,3) = 0.5*(1 - 3*delta);
w2(:,:,4) = 0.5*delta;
w2 = w2 / Mesh.spacing(3)^2;

% Third moment
w3 = zeros(Particle.nParticle,3,nGrid); % nParticle x nDim x nGrid
w3(:,:,1) = -(1/6);
w3(:,:,2) = 0.5;
w3(:,:,3) = -0.5;
w3(:,:,4) = (1/6);
w3 = w3 / Mesh.spacing(3)^3;


%% Distribute forces and dipoles
for i = 1:nGrid
    for j = 1:nGrid
        for k = 1:nGrid
            % linear index for gridForce arrays
            idx(:,1) = iX(:,i) + Mesh.nMesh(1)*(iY(:,j) - 1) ...
                + Mesh.nMesh(1)^2 .*(iZ(:,k) - 1); % nParticle x 1
            idx(:,2) = iX(:,i) + Mesh.nMesh(2)*(iY(:,j) - 1) ...
                + Mesh.nMesh(2)^2 .*(iZ(:,k) - 1); % nParticle x 1
            idx(:,3) = iX(:,i) + Mesh.nMesh(3)*(iY(:,j) - 1) ...
                + Mesh.nMesh(3)^2 .*(iZ(:,k) - 1); % nParticle x 
           
            % distribute forces
            gridForceX(idx(:,1)) = gridForceX(idx(:,1)) ...
                + Particle.farFieldForce(:,1) .* (w0(:,1,i) .* w0(:,2,j) .* w0(:,3,k));
            gridForceY(idx(:,2)) = gridForceY(idx(:,2)) ...
                + Particle.farFieldForce(:,2) .* (w0(:,1,i) .* w0(:,2,j) .* w0(:,3,k));
            gridForceZ(idx(:,3)) = gridForceZ(idx(:,3)) ...
                + Particle.farFieldForce(:,3) .* (w0(:,1,i) .* w0(:,2,j) .* w0(:,3,k));
            
            % distribute dipoles
            gridForceX(idx(:,1)) = gridForceX(idx(:,1)) ...
                + Particle.farFieldDipoleHydro(:,1,1) .* (w1(:,1,i) .* w0(:,2,j) .* w0(:,3,k)) ...
                + Particle.farFieldDipoleHydro(:,2,1) .* (w0(:,1,i) .* w1(:,2,j) .* w0(:,3,k)) ...
                + Particle.farFieldDipoleHydro(:,3,1) .* (w0(:,1,i) .* w0(:,2,j) .* w1(:,3,k));
            gridForceY(idx(:,2)) = gridForceY(idx(:,2)) ...
                + Particle.farFieldDipoleHydro(:,1,2) .* (w1(:,1,i) .* w0(:,2,j) .* w0(:,3,k)) ...
                + Particle.farFieldDipoleHydro(:,2,2) .* (w0(:,1,i) .* w1(:,2,j) .* w0(:,3,k)) ...
                + Particle.farFieldDipoleHydro(:,3,2) .* (w0(:,1,i) .* w0(:,2,j) .* w1(:,3,k));
            gridForceZ(idx(:,3)) = gridForceZ(idx(:,3)) ...
                + Particle.farFieldDipoleHydro(:,1,3) .* (w1(:,1,i) .* w0(:,2,j) .* w0(:,3,k)) ...
                + Particle.farFieldDipoleHydro(:,2,3) .* (w0(:,1,i) .* w1(:,2,j) .* w0(:,3,k)) ...
                + Particle.farFieldDipoleHydro(:,3,3) .* (w0(:,1,i) .* w0(:,2,j) .* w1(:,3,k));
            
            % distribute quadrupoles
            gridForceX(idx(:,1)) = gridForceX(idx(:,1)) ...
                + (1/3) * Particle.farFieldForce(:,1) .* (w2(:,1,i) .* w0(:,2,j) .* w0(:,3,k) ...
                + w0(:,1,i) .* w2(:,2,j) .* w0(:,3,k) + w0(:,1,i) .* w0(:,2,j) .* w2(:,3,k));
            gridForceY(idx(:,2)) = gridForceY(idx(:,2)) ...
                + (1/3) * Particle.farFieldForce(:,2) .* (w2(:,1,i) .* w0(:,2,j) .* w0(:,3,k) ...
                + w0(:,1,i) .* w2(:,2,j) .* w0(:,3,k) + w0(:,1,i) .* w0(:,2,j) .* w2(:,3,k));
            gridForceZ(idx(:,3)) = gridForceZ(idx(:,3)) ...
                + (1/3) * Particle.farFieldForce(:,3) .* (w2(:,1,i) .* w0(:,2,j) .* w0(:,3,k) ...
                + w0(:,1,i) .* w2(:,2,j) .* w0(:,3,k) + w0(:,1,i) .* w0(:,2,j) .* w2(:,3,k));
            
            % distribute octuples
            gridForceX(idx(:,1)) = gridForceX(idx(:,1)) ...
                + (3/5) * Particle.farFieldDipoleHydro(:,1,1) .* (w3(:,1,i) .* w0(:,2,j) .* w0(:,3,k) ...
                + w1(:,1,i) .* w2(:,2,j) .* w0(:,3,k) + w1(:,1,i) .* w0(:,2,j) .* w2(:,3,k)) ...
                + (3/5) * Particle.farFieldDipoleHydro(:,2,1) .* (w2(:,1,i) .* w1(:,2,j) .* w0(:,3,k) ...
                + w0(:,1,i) .* w3(:,2,j) .* w0(:,3,k) + w0(:,1,i) .* w1(:,2,j) .* w2(:,3,k)) ...
                + (3/5) * Particle.farFieldDipoleHydro(:,3,1) .* (w2(:,1,i) .* w0(:,2,j) .* w1(:,3,k) ...
                + w0(:,1,i) .* w2(:,2,j) .* w1(:,3,k) + w0(:,1,i) .* w0(:,2,j) .* w3(:,3,k));
            gridForceY(idx(:,2)) = gridForceY(idx(:,2)) ...
                + (3/5) * Particle.farFieldDipoleHydro(:,1,2) .* (w3(:,1,i) .* w0(:,2,j) .* w0(:,3,k) ...
                + w1(:,1,i) .* w2(:,2,j) .* w0(:,3,k) + w1(:,1,i) .* w0(:,2,j) .* w2(:,3,k)) ...
                + (3/5) * Particle.farFieldDipoleHydro(:,2,2) .* (w2(:,1,i) .* w1(:,2,j) .* w0(:,3,k) ...
                + w0(:,1,i) .* w3(:,2,j) .* w0(:,3,k) + w0(:,1,i) .* w1(:,2,j) .* w2(:,3,k)) ...
                + (3/5) * Particle.farFieldDipoleHydro(:,3,2) .* (w2(:,1,i) .* w0(:,2,j) .* w1(:,3,k) ...
                + w0(:,1,i) .* w2(:,2,j) .* w1(:,3,k) + w0(:,1,i) .* w0(:,2,j) .* w3(:,3,k));
            gridForceZ(idx(:,3)) = gridForceZ(idx(:,3)) ...
                + (3/5) * Particle.farFieldDipoleHydro(:,1,3) .* (w3(:,1,i) .* w0(:,2,j) .* w0(:,3,k) ...
                + w1(:,1,i) .* w2(:,2,j) .* w0(:,3,k) + w1(:,1,i) .* w0(:,2,j) .* w2(:,3,k)) ...
                + (3/5) * Particle.farFieldDipoleHydro(:,2,3) .* (w2(:,1,i) .* w1(:,2,j) .* w0(:,3,k) ...
                + w0(:,1,i) .* w3(:,2,j) .* w0(:,3,k) + w0(:,1,i) .* w1(:,2,j) .* w2(:,3,k)) ...
                + (3/5) * Particle.farFieldDipoleHydro(:,3,3) .* (w2(:,1,i) .* w0(:,2,j) .* w1(:,3,k) ...
                + w0(:,1,i) .* w2(:,2,j) .* w1(:,3,k) + w0(:,1,i) .* w0(:,2,j) .* w3(:,3,k));
        end
    end
end

%% Test distribution
% force = zeros(Particle.nParticle,3);
% dipole = zeros(Particle.nParticle,3,3);
% 
% for nP = 1:Particle.nParticle
%     for i = 1:nGrid
%         for j = 1:nGrid
%             for k = 1:nGrid
%             force(nP,1) = force(nP,1) + gridForceX(iX(nP,i),iY(nP,j),iZ(nP,k));
%             force(nP,2) = force(nP,2) + gridForceY(iX(nP,i),iY(nP,j),iZ(nP,k));
%             force(nP,3) = force(nP,3) + gridForceZ(iX(nP,i),iY(nP,j),iZ(nP,k));
%             
%             
%             dipole(nP,1,1) = dipole(nP,1,1) + ...
%                 Mesh.spacing*(n(i) - delta(nP,1)) * ...
%                 gridForceX(iX(nP,i),iY(nP,j),iZ(nP,k));
%             dipole(nP,1,2) = dipole(nP,1,2) + ...
%                 Mesh.spacing*(n(i) - delta(nP,1)) * ...
%                 gridForceY(iX(nP,i),iY(nP,j),iZ(nP,k));
%             dipole(nP,1,3) = dipole(nP,1,3) + ...
%                 Mesh.spacing*(n(i) - delta(nP,1)) * ...
%                 gridForceZ(iX(nP,i),iY(nP,j),iZ(nP,k));
%             
%             dipole(nP,2,1) = dipole(nP,2,1) + ...
%                 Mesh.spacing*(n(j) - delta(nP,2)) * ...
%                 gridForceX(iX(nP,i),iY(nP,j),iZ(nP,k));
%             dipole(nP,2,2) = dipole(nP,2,2) + ...
%                 Mesh.spacing*(n(j) - delta(nP,2)) * ...
%                 gridForceY(iX(nP,i),iY(nP,j),iZ(nP,k));
%             dipole(nP,2,3) = dipole(nP,2,3) + ...
%                 Mesh.spacing*(n(j) - delta(nP,2)) * ...
%                 gridForceZ(iX(nP,i),iY(nP,j),iZ(nP,k));
%             
%             dipole(nP,3,1) = dipole(nP,3,1) + ...
%                 Mesh.spacing*(n(k) - delta(nP,3)) * ...
%                 gridForceX(iX(nP,i),iY(nP,j),iZ(nP,k));
%             dipole(nP,3,2) = dipole(nP,3,2) + ...
%                 Mesh.spacing*(n(k) - delta(nP,3)) * ...
%                 gridForceY(iX(nP,i),iY(nP,j),iZ(nP,k));
%             dipole(nP,3,3) = dipole(nP,3,3) + ...
%                 Mesh.spacing*(n(k) - delta(nP,3)) * ...
%                 gridForceZ(iX(nP,i),iY(nP,j),iZ(nP,k));
%             end
%         end
%     end
% end
% errorForce = max(max(abs(Particle.farFieldForce - force)))    
% errorDipole = max(max(max(abs(Particle.farFieldDipoleHydro - dipole))))


%--------------------------------------------------------------------------
function [velocity, velocityGradient, velocityLaplacian, ...
    velocityLaplacianGradient] = interpolateVelocity(Parameter, Particle,...
    Mesh, gridVelocityX, gridVelocityY, gridVelocityZ)
% Interpolate the velocity and velocity gradient at particle centers

nGrid = 4; % third-order polynomial interpolation

%% Preallocate velocity arrays
velocity = zeros(Particle.nParticle,3);
velocityLaplacian = zeros(Particle.nParticle,3);
velocityGradient = zeros(Particle.nParticle,3,3);
velocityLaplacianGradient = zeros(Particle.nParticle,3,3);

%% Grid points
% Find indices of nearest grid points (to the left for nGrid = 4)
gridIndex = zeros(Particle.nParticle, 3);
gridIndex(:,1) = mod(floor( Particle.position(:,1) / Mesh.spacing(1) ),...
    Mesh.nMesh(1)) + 1;
gridIndex(:,2) = mod(floor( Particle.position(:,2) / Mesh.spacing(2) ),...
    Mesh.nMesh(2)) + 1;
gridIndex(:,3) = mod(floor( Particle.position(:,3) / Mesh.spacing(3) ),...
    Mesh.nMesh(3)) + 1;

% Find positions of nearest grid points
gridPosition = zeros(Particle.nParticle, 3);
gridPosition(:,1) = Mesh.x1( gridIndex(:,1) );
gridPosition(:,2) = Mesh.x2( gridIndex(:,2) );
gridPosition(:,3) = Mesh.x3( gridIndex(:,3) );

% Fractional displacement between particles and proximal grid point
delta = Particle.position - gridPosition;
delta(:,1) = delta(:,1) - round(delta(:,1) / Parameter.domainLength(1)) ...
    * Parameter.domainLength(1);
delta(:,1) = delta(:,1) / Mesh.spacing(1);
delta(:,2) = delta(:,2) - round(delta(:,2) / Parameter.domainLength(2)) ...
    * Parameter.domainLength(2);
delta(:,2) = delta(:,2) / Mesh.spacing(2);
delta(:,3) = delta(:,3) - round(delta(:,3) / Parameter.domainLength(3)) ...
    * Parameter.domainLength(3);
delta(:,3) = delta(:,3) / Mesh.spacing(3);

% Indices of grid points around each particle
n = [-1:2]; % for nGrid = 4
iX = mod( bsxfun(@plus, n , gridIndex(:,1)) - 1, Mesh.nMesh(1)) + 1; % nParticle x nGrid
iY = mod( bsxfun(@plus, n , gridIndex(:,2)) - 1, Mesh.nMesh(2)) + 1;
iZ = mod( bsxfun(@plus, n , gridIndex(:,3)) - 1, Mesh.nMesh(3)) + 1;


%% Weighting functions
% Zeroth moment
w0 = zeros(Particle.nParticle,3,nGrid); % nParticle x nDim x nGrid
w0(:,:,1) = -(1/6) * (delta - 2) .* (delta - 1) .* delta;
w0(:,:,2) = 0.5 * (delta - 2) .* (delta - 1) .* (1 + delta);
w0(:,:,3) = 0.5 * delta .* (2 + delta - delta.^2);
w0(:,:,4) = (1/6) * delta .* (-1 + delta.^2);

% First moment
w1 = zeros(Particle.nParticle,3,nGrid); % nParticle x nDim x nGrid
w1(:,:,1) = -(1/3) + delta - 0.5*delta.^2;
w1(:,:,2) = 0.5 * (-1 + delta .* (-4 + 3*delta));
w1(:,:,3) = 1 + delta - 1.5*delta.^2;
w1(:,:,4) = -(1/6) + 0.5*delta.^2;
w1 = w1 / Mesh.spacing(3);

% Second moment
w2 = zeros(Particle.nParticle,3,nGrid); % nParticle x nDim x nGrid
w2(:,:,1) = 1 - delta;
w2(:,:,2) = -2 + 3*delta;
w2(:,:,3) = 1 - 3*delta;
w2(:,:,4) = delta;
w2 = w2 / Mesh.spacing(3)^2;

% Third moment
w3 = zeros(Particle.nParticle,3,nGrid); % nParticle x nDim x nGrid
w3(:,:,1) = -1;
w3(:,:,2) = 3;
w3(:,:,3) = -3;
w3(:,:,4) = 1;
w3 = w3 / Mesh.spacing(3)^3;


%% Interpolate velocity and velocity gradients
for i = 1:nGrid
    for j = 1:nGrid
        for k = 1:nGrid
            % linear index for gridForce arrays
            idx(:,1) = iX(:,i) + Mesh.nMesh(1)*(iY(:,j) - 1) ...
                + Mesh.nMesh(1)^2 .*(iZ(:,k) - 1); % nParticle x 1
            idx(:,2) = iX(:,i) + Mesh.nMesh(2)*(iY(:,j) - 1) ...
                + Mesh.nMesh(2)^2 .*(iZ(:,k) - 1); % nParticle x 1
            idx(:,3) = iX(:,i) + Mesh.nMesh(3)*(iY(:,j) - 1) ...
                + Mesh.nMesh(3)^2 .*(iZ(:,k) - 1); % nParticle x 
           
            % interpolate velocity
            velocity(:,1) = velocity(:,1) + ...
                gridVelocityX(idx(:,1)) .* (w0(:,1,i) .* w0(:,2,j) .* w0(:,3,k)); % nParticle x 1
            velocity(:,2) = velocity(:,2) + ...
                gridVelocityY(idx(:,2)) .* (w0(:,1,i) .* w0(:,2,j) .* w0(:,3,k)); 
            velocity(:,3) = velocity(:,3) + ...
                gridVelocityZ(idx(:,3)) .* (w0(:,1,i) .* w0(:,2,j) .* w0(:,3,k));
           
            % interpolate velocity gradient
            velocityGradient(:,1,1) = velocityGradient(:,1,1) + ...
                gridVelocityX(idx(:,1)) .* (w1(:,1,i) .* w0(:,2,j) .* w0(:,3,k));
            velocityGradient(:,1,2) = velocityGradient(:,1,2) + ...
                gridVelocityY(idx(:,2)) .* (w1(:,1,i) .* w0(:,2,j) .* w0(:,3,k));
            velocityGradient(:,1,3) = velocityGradient(:,1,3) + ...
                gridVelocityZ(idx(:,3)) .* (w1(:,1,i) .* w0(:,2,j) .* w0(:,3,k));
            
            velocityGradient(:,2,1) = velocityGradient(:,2,1) + ...
                gridVelocityX(idx(:,1)) .* (w0(:,1,i) .* w1(:,2,j) .* w0(:,3,k));
            velocityGradient(:,2,2) = velocityGradient(:,2,2) + ...
                gridVelocityY(idx(:,2)) .* (w0(:,1,i) .* w1(:,2,j) .* w0(:,3,k));
            velocityGradient(:,2,3) = velocityGradient(:,2,3) + ...
                gridVelocityZ(idx(:,3)) .* (w0(:,1,i) .* w1(:,2,j) .* w0(:,3,k));
            
            velocityGradient(:,3,1) = velocityGradient(:,3,1) + ...
                gridVelocityX(idx(:,1)) .* (w0(:,1,i) .* w0(:,2,j) .* w1(:,3,k));
            velocityGradient(:,3,2) = velocityGradient(:,3,2) + ...
                gridVelocityY(idx(:,2)) .* (w0(:,1,i) .* w0(:,2,j) .* w1(:,3,k));
            velocityGradient(:,3,3) = velocityGradient(:,3,3) + ...
                gridVelocityZ(idx(:,3)) .* (w0(:,1,i) .* w0(:,2,j) .* w1(:,3,k));
            
            % interpolate velocity laplacian
            velocityLaplacian(:,1) = velocityLaplacian(:,1) + ...
                gridVelocityX(idx(:,1)) .* (w2(:,1,i) .* w0(:,2,j) .* w0(:,3,k) ...
                + w0(:,1,i) .* w2(:,2,j) .* w0(:,3,k) ...
                + w0(:,1,i) .* w0(:,2,j) .* w2(:,3,k)); 
            velocityLaplacian(:,2) = velocityLaplacian(:,2) + ...
                gridVelocityY(idx(:,2)) .* (w2(:,1,i) .* w0(:,2,j) .* w0(:,3,k) ...
                + w0(:,1,i) .* w2(:,2,j) .* w0(:,3,k) ...
                + w0(:,1,i) .* w0(:,2,j) .* w2(:,3,k)); 
            velocityLaplacian(:,3) = velocityLaplacian(:,3) + ...
                gridVelocityZ(idx(:,3)) .* (w2(:,1,i) .* w0(:,2,j) .* w0(:,3,k) ...
                + w0(:,1,i) .* w2(:,2,j) .* w0(:,3,k) ...
                + w0(:,1,i) .* w0(:,2,j) .* w2(:,3,k)); 
            
            % interpolate velocity laplacian gradient Dk * Dk * Di * uj
            velocityLaplacianGradient(:,1,1) = velocityLaplacianGradient(:,1,1) + ...
                gridVelocityX(idx(:,1)) .* (w3(:,1,i) .* w0(:,2,j) .* w0(:,3,k) + ...
                w1(:,1,i) .* w2(:,2,j) .* w0(:,3,k) + ...
                w1(:,1,i) .* w0(:,2,j) .* w2(:,3,k));
            velocityLaplacianGradient(:,1,2) = velocityLaplacianGradient(:,1,2) + ...
                gridVelocityY(idx(:,2)) .* (w3(:,1,i) .* w0(:,2,j) .* w0(:,3,k) + ...
                w1(:,1,i) .* w2(:,2,j) .* w0(:,3,k) + ...
                w1(:,1,i) .* w0(:,2,j) .* w2(:,3,k));
            velocityLaplacianGradient(:,1,3) = velocityLaplacianGradient(:,1,3) + ...
                gridVelocityZ(idx(:,3)) .* (w3(:,1,i) .* w0(:,2,j) .* w0(:,3,k) + ...
                w1(:,1,i) .* w2(:,2,j) .* w0(:,3,k) + ...
                w1(:,1,i) .* w0(:,2,j) .* w2(:,3,k));
            
            velocityLaplacianGradient(:,2,1) = velocityLaplacianGradient(:,2,1) + ...
                gridVelocityX(idx(:,1)) .* (w2(:,1,i) .* w1(:,2,j) .* w0(:,3,k) + ...
                w0(:,1,i) .* w3(:,2,j) .* w0(:,3,k) + ...
                w0(:,1,i) .* w1(:,2,j) .* w2(:,3,k));
            velocityLaplacianGradient(:,2,2) = velocityLaplacianGradient(:,2,2) + ...
                gridVelocityY(idx(:,2)) .* (w2(:,1,i) .* w1(:,2,j) .* w0(:,3,k) + ...
                w0(:,1,i) .* w3(:,2,j) .* w0(:,3,k) + ...
                w0(:,1,i) .* w1(:,2,j) .* w2(:,3,k));
            velocityLaplacianGradient(:,2,3) = velocityLaplacianGradient(:,2,3) + ...
                gridVelocityZ(idx(:,3)) .* (w2(:,1,i) .* w1(:,2,j) .* w0(:,3,k) + ...
                w0(:,1,i) .* w3(:,2,j) .* w0(:,3,k) + ...
                w0(:,1,i) .* w1(:,2,j) .* w2(:,3,k));
            
            velocityLaplacianGradient(:,3,1) = velocityLaplacianGradient(:,3,1) + ...
                gridVelocityX(idx(:,1)) .* (w2(:,1,i) .* w0(:,2,j) .* w1(:,3,k) + ...
                w0(:,1,i) .* w2(:,2,j) .* w1(:,3,k) + ...
                w0(:,1,i) .* w0(:,2,j) .* w3(:,3,k));
            velocityLaplacianGradient(:,3,2) = velocityLaplacianGradient(:,3,2) + ...
                gridVelocityY(idx(:,2)) .* (w2(:,1,i) .* w0(:,2,j) .* w1(:,3,k) + ...
                w0(:,1,i) .* w2(:,2,j) .* w1(:,3,k) + ...
                w0(:,1,i) .* w0(:,2,j) .* w3(:,3,k));
            velocityLaplacianGradient(:,3,3) = velocityLaplacianGradient(:,3,3) + ...
                gridVelocityZ(idx(:,3)) .* (w2(:,1,i) .* w0(:,2,j) .* w1(:,3,k) + ...
                w0(:,1,i) .* w2(:,2,j) .* w1(:,3,k) + ...
                w0(:,1,i) .* w0(:,2,j) .* w3(:,3,k));
        end
    end
end
