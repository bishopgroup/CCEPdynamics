function Particle = initializeForce(Particle)
%--------------------------------------------------------------------------
% initializeForce - initializes particle force, torque, stresslet.
%
% Particle = initializeForce(Particle)
%--------------------------------------------------------------------------

%% Force - unit force in the negative z-direction
%Particle.force = zeros(Particle.nParticle, 3);
%Particle.force(:,3) = -1;

%% Torque and Stresslet
Particle.torque = zeros(Particle.nParticle, 3);
Particle.dipole = zeros(Particle.nParticle, 9);

