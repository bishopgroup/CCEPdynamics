function [Parameter] = initializeParticlePosition(Parameter)
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

% update Parameter structure
Parameter.domainLength = Parameter.L3*[1 1 1]';% L1 L2 gets changed in the initializeMesh
%Parameter.splittingAlpha = (1/4) * Parameter.domainLength(3).^2; 
Parameter.splittingAlpha = 10; 

% update Particle structure
%Particle.nParticle = Parameter.nUnitCell^3 * nParticle;
%Particle.position = zeros(Parameter.nUnitCell, 3);


