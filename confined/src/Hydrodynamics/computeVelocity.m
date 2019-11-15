function Particle = computeVelocity(Parameter, Particle, Mesh)
%%-------------------------------------------------------------------------
% computeVelocity - computes the far field contribution to the velocity
% at the particle centers using the Ewald method.
%
% Particle = computeVelocity(Parameter, Particle, Mesh)
% %--------------------------------------------------------------------------
% 
% Particle.farFieldForce = zeros(Particle.nParticle,3);
% Particle.farFieldDipoleHydro = ones(Particle.nParticle,3,3);
% [globalVelocity, globalVelocityGradient, globalVelocityLaplacian] = ...
%     mFarFieldWave(Parameter, Particle, Mesh);
% 
% 
% [selfVelocity, selfVelocityGradient,selfVelocityLaplacian] = ...
%     mFarFieldSelf(Parameter, Particle);  % wave-space 'self' contribution
% globalVelocityGradient
% selfVelocityGradient
% Particle.velocity = Particle.farFieldForce ;
% pause
% return
% %for nP = 1:Particle.nParticle
% %    Particle.farFieldDipoleHydro(nP,:,:) = [1 2 3; 4 5 6; 7 8 9];
% %end
% 
% [globalVelocity, globalVelocityGradient, globalVelocityLaplacian] = ...
%     mFarFieldWave(Parameter, Particle, Mesh);
% [selfVelocity, selfVelocityLaplacian] = ...
%     mFarFieldSelf(Parameter, Particle);
% [localVelocity] = mFarFieldReal(Parameter, Particle);
% 
% Particle.velocity = Particle.force + ...
%     globalVelocity - selfVelocity + localVelocity;% + ...
%  %   (globalVelocityLaplacian - selfVelocityLaplacian)/6;
% globalVelocityLaplacian
% selfVelocityLaplacian
% %Particle.velocityLaplacian = globalVelocityLaplacian - selfVelocityLaplacian;
% 

Np = Particle.nParticle; % abbreviated name for number of particles

if Parameter.isNearFieldHydrodynamics
    nearFieldResistance = rNearFieldHydro(Parameter,Particle);
else
    nearFieldResistance = spalloc(6*Particle.nParticle, ...
        6*Particle.nParticle,1);
end

%% Assemble Right-hand-side (Force, Torque, Strain Rate)
rightHandSide = zeros(11*Np, 1);

% Applied force
rightHandSide(1:3*Np) = reshape(Particle.force, 3*Np, 1); 

% Applied torque
rightHandSide(3*Np+1:6*Np) = reshape(Particle.torque, 3*Np, 1); 

% Strain rate (always zero for solid particles)
rightHandSide(6*Np+1:end) = 0; 


%% Solve for Far Field Force, Torque, and Stresslet

% initial guess
x0 = [reshape(Particle.force, 3*Np, 1); ...
      reshape(Particle.torque, 3*Np, 1);...
      zeros(5*Np, 1)];
  
  [lhs,velocity] = leftHandSide(x0, Parameter, Particle, Mesh, nearFieldResistance);

% Solver parameters
tolerance = 1e-1;
restart = 3*Np;
maxint = 20;

% Solve
[x,flag,relres,iter,resvec] = ...
     gmres(@(x)leftHandSide(x,Parameter,Particle,Mesh,nearFieldResistance), ...
     rightHandSide, restart, tolerance, maxint, [], [], x0);

 % Check result
if flag ~= 0
    semilogy(1:length(resvec),resvec,'-o');
    xlabel('Iteration number');
    ylabel('Relative residual');
    error(sprintf('Houston, we have a problem: flag = %d\n',flag));
end


%% Save results
%Particle.farFieldForce = reshape(x(1:3*Particle.nParticle), ...
%    Particle.nParticle, 3);
%Particle.farFieldDipoleHydro = reshape(x(3*Particle.nParticle+1:end),...
%    Particle.nParticle, 3,3);


[lhs,velocity] = leftHandSide(x, Parameter, Particle, Mesh, nearFieldResistance);
Particle.velocity = velocity;
%Particle.dipole = Particle.farFieldDipoleHydro;
%tmp = reshape(Particle.farFieldDipoleHydro(1,:,:),3,3)

%S = 0.5*(tmp + tmp')
%W = 0.5*(tmp - tmp')

Particle.farFieldTorque = reshape( x(3*Np+1:6*Np), Np, 3);
Particle.farFieldStresslet = reshape( x(6*Np+1:end), Np, 5); % xx, xy, xz, yy, yz
%Particle.farFieldStresslet

%--------------------------------------------------------------------------
function [lhs,velocity] = leftHandSide(x, Parameter, Particle, ...
    Mesh, nearFieldResistance)

Np = Particle.nParticle; % abbreviated name for number of particles

%% Set-up force moments

% Far field force, torque, and stresslet
Particle.farFieldForce = reshape( x(1:3*Np), Np, 3);
Particle.farFieldTorque = reshape( x(3*Np+1:6*Np), Np, 3);
Particle.farFieldStresslet = reshape( x(6*Np+1:end), Np, 5); % xx, xy, xz, yy, yz

% Builing the force dipole Dij
%  Li = eijk * Djk = eijk * Ajk with Ajk = 0.5 * (Djk - Dkj) 
%   
%  L1   D23 - D32              0  L3 -L2 
%  L2 = D31 - D13  and 2*A = -L3   0  L1
%  L3   D12 - D21             L2 -L1   0 
%
%  Sij = 0.5*( Dij + Dji) - Dkk/3
%
%  Dij = Sij + Aij
%
Particle.farFieldDipoleHydro = zeros(Np, 3, 3);

Particle.farFieldDipoleHydro(:,1,1) = Particle.farFieldStresslet(:,1);
Particle.farFieldDipoleHydro(:,1,2) = Particle.farFieldStresslet(:,2) ...
    + 0.5 * Particle.farFieldTorque(:,3);
Particle.farFieldDipoleHydro(:,1,3) = Particle.farFieldStresslet(:,3) ...
    - 0.5 * Particle.farFieldTorque(:,2);

Particle.farFieldDipoleHydro(:,2,1) = Particle.farFieldStresslet(:,2) ...
    - 0.5 * Particle.farFieldTorque(:,3);
Particle.farFieldDipoleHydro(:,2,2) = Particle.farFieldStresslet(:,4);
Particle.farFieldDipoleHydro(:,2,3) = Particle.farFieldStresslet(:,5) ...
    + 0.5 * Particle.farFieldTorque(:,1);

Particle.farFieldDipoleHydro(:,3,1) = Particle.farFieldStresslet(:,3) ...
    + 0.5 * Particle.farFieldTorque(:,2);
Particle.farFieldDipoleHydro(:,3,2) = Particle.farFieldStresslet(:,5) ...
    - 0.5 * Particle.farFieldTorque(:,1);
Particle.farFieldDipoleHydro(:,3,3) = -Particle.farFieldDipoleHydro(:,1,1) ...
    - Particle.farFieldDipoleHydro(:,2,2); % stresslet is traceless


%% Compute the disturbance velocity and velocity gradients 

% Ewald components
[globalVelocity, globalVelocityGradient,globalVelocityLaplacian, globalVelocityLaplacianGradient] = ...
    mFarFieldWave(Parameter, Particle, Mesh); % wave-space contribution
[selfVelocity, selfVelocityGradient,selfVelocityLaplacian, selfVelocityLaplacianGradient] = ...
    mFarFieldSelf(Parameter, Particle);       % wave-space 'self' contribution
[localVelocity, localVelocityGradient, localVelocityLaplacian] = ...
    mFarFieldReal(Parameter, Particle);       % real-space contribution

%reshape(localVelocity(1,:),3,1)
% reshape(globalVelocityGradient(1,:,:),3,3)
% globalVelocityGradient(:,1,1)
% globalVelocityGradient(:,2,2)
% globalVelocityGradient(:,3,3)
% (globalVelocityGradient(:,1,1) + globalVelocityGradient(:,2,2) + globalVelocityGradient(:,3,3) ) ./globalVelocityGradient(:,3,3)

%  Disturbance velocity, velocity gradient, and velocity laplacian
disturbanceVelocity = globalVelocity - selfVelocity + localVelocity;
disturbanceVelocityGradient = globalVelocityGradient ...
    - selfVelocityGradient + localVelocityGradient;
disturbanceVelocityLaplacian = globalVelocityLaplacian ...
    - selfVelocityLaplacian + localVelocityLaplacian;
disturbanceVelocityLaplacianGradient = globalVelocityLaplacianGradient ...
    - selfVelocityLaplacianGradient;


%% Compute vorticity and rate of strain
%
% velocity gradient Gij can be decomposed into a symmetric part Eij (the 
% rate of strain tensor) and an antisymmetric part Wij (the vorticity 
% tensor).
%  
%   Eij = 0.5 * (Gij + Gji)
%   Wij = 0.5 * (Gij - Gji)
%   Gij = Eij + Wij
%
% The vorticity tensor Wij is related to the vorticity vector wi as
%   
%   wi = eijk * Wjk  and    Wij = 0.5 * eijk * wk 
%  
%   w1   W23 - W32              0  w3 -w2 
%   w2 = W31 - W13  and 2*W = -w3   0  w1
%   w3   W12 - W21             w2 -w1   0 
%

% Vorticity at particle centers
% disturbanceVorticity = zeros(Np, 3);
% disturbanceVorticity(:,1) = disturbanceVelocityGradient(:,2,3) ...
%     - disturbanceVelocityGradient(:,3,2);
% disturbanceVorticity(:,2) = disturbanceVelocityGradient(:,3,1) ...
%     - disturbanceVelocityGradient(:,1,3);
% disturbanceVorticity(:,3) = disturbanceVelocityGradient(:,1,2) ...
%     - disturbanceVelocityGradient(:,2,1);

% Strain rate at particle centers
disturbanceStrainRate = zeros(Np, 5); % xx, xy, xz, yy, yz
disturbanceStrainRate(:,1) = disturbanceVelocityGradient(:,1,1);
disturbanceStrainRate(:,2) = 0.5 * (disturbanceVelocityGradient(:,1,2)...
    + disturbanceVelocityGradient(:,2,1));
disturbanceStrainRate(:,3) = 0.5 * (disturbanceVelocityGradient(:,1,3)...
    + disturbanceVelocityGradient(:,3,1));
disturbanceStrainRate(:,4) = disturbanceVelocityGradient(:,2,2);
disturbanceStrainRate(:,5) = 0.5 * (disturbanceVelocityGradient(:,2,3)...
    + disturbanceVelocityGradient(:,3,2));

% disturbanceVelocityLaplacianGradient(:,1,1) + disturbanceVelocityLaplacianGradient(:,2,2)...
%     + disturbanceVelocityLaplacianGradient(:,3,3)

disturbanceStrainRate(:,1) = disturbanceStrainRate(:,1) ...
    + 0.1 * disturbanceVelocityLaplacianGradient(:,1,1);
disturbanceStrainRate(:,2) = disturbanceStrainRate(:,2) ...
    + 0.1 * 0.5 * (disturbanceVelocityLaplacianGradient(:,1,2)...
    + disturbanceVelocityLaplacianGradient(:,2,1));
disturbanceStrainRate(:,3) = disturbanceStrainRate(:,3) ...
    + 0.1 * 0.5 * (disturbanceVelocityLaplacianGradient(:,1,3)...
    + disturbanceVelocityLaplacianGradient(:,3,1));
disturbanceStrainRate(:,4) = disturbanceStrainRate(:,4) ...
    + 0.1 * disturbanceVelocityLaplacianGradient(:,2,2);
disturbanceStrainRate(:,5) = disturbanceStrainRate(:,5) ...
    + 0.1 * 0.5 * (disturbanceVelocityLaplacianGradient(:,2,3)...
    + disturbanceVelocityLaplacianGradient(:,3,2));


%% Compute left-hand-side (Faxen Laws)
lhs = zeros(11*Particle.nParticle, 1);
VW = [disturbanceVelocity; zeros(5*Particle.nParticle, 3)];

% % Compute "A matrix"
% bodyFixedPosition = bsxfun(@minus,Particle.position,Particle.clusterTranslation);
% z = zeros(3);
% A = [repmat(eye(3),1,Particle.nParticle),repmat(zeros(3),1,Particle.nParticle);repmat(zeros(3),1,Particle.nParticle),repmat(eye(3),1,Particle.nParticle)];
%      
% for iPart= 1:Particle.nParticle
%     A(4:end,3*(iPart-1)+1:3*iPart) = -(Signature(:,:,1)*bodyFixedPosition(1,iPart) + Signature(:,:,2)*bodyFixedPosition(2,iPart) ...
%         + Signature(:,:,3)*bodyFixedPosition(3,iPart));
% end

nearFieldFT = nearFieldResistance * VW;

% Force law (near field here)
lhs(1:3*Np) = Particle.farFieldForce + nearFieldFT(1:Np,:);

% Torque law (near field here)
lhs(3*Np+1:6*Np) = Particle.farFieldTorque;
%+ reshape(nearFieldFT(Particle.nParticle+1:end), 3, Particle.nParticle);

% Stresslet law (no near field here)
lhs(6*Np+1:end) = Particle.farFieldStresslet + (10/9) * disturbanceStrainRate;

%% Compute particle velocity
velocity = Particle.farFieldForce + disturbanceVelocity ...
    + disturbanceVelocityLaplacian/6;
