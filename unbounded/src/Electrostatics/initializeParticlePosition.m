function [Particle] = initializeParticlePosition(Parameter, Particle )
%%-------------------------------------------------------------------------
%  - initialize particle positions.
%
% [Parameter, Particle] = initializeParticleCubic(Parameter) prepares the 
% data structure Particle, which contains the following fields:
%
% nParticle: total number of particles in the simulation.
%
% position: (nParticle x 3) array of x,y,z coordinates of each particle.
%
% Particles are initialized on a simple cubic lattice with a desired 
% volume fraction. The domain contains nUnitCell unit cells in each of the 
% 3 directions.
%
% This function also specifies the splitting parameter alpha according to
% the criterion that 4*sqrt(alpha/pi) = L/2 where L is the domain length.
% This choice ensures that real space interactions beyond distances greater
% than half the domain length are neglegible.
%--------------------------------------------------------------------------

%% Simulation domain
Particle.nParticle = Parameter.nParticle;
Particle.position = zeros(3, Parameter.nParticle);
for np = 1:Particle.nParticle
    nOverlap = 1;
    while nOverlap > 0
        Particle.position(1:2,np) = rand(1,2) .* Parameter.domainLength(1:2);
        Particle.position(3,np) = 1 + rand(1) * (Parameter.domainLength(3) - 2);
        nOverlap = checkOverlap( Parameter, Particle );
    end
end

Particle.charge = 4*rand(Particle.nParticle, 1) - 2; % nParticle x 1
Particle.field = [zeros(Particle.nParticle, 2), ones(Particle.nParticle,1)]';  % 3 x nParticle
Particle.potential =  Particle.charge;  % initial guess
Particle.farfieldCharge = Particle.charge;  % initial guess 
Particle.farfieldDipole = Particle.field;  % initial guess 
Particle.force = bsxfun(@times,Particle.field,Particle.charge'); % 3 x nParticle

