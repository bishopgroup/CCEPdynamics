function [U,GradU] = pFarFieldSelf(Parameter,Particle)
%{
pFarFieldSelf - computes the self-contribution to the 'global' potential.
%}

U = (1 - 2 / Parameter.splittingAlpha^0.5) * Particle.farfieldCharge; % potential
GradU = (1 - 4 * pi / (3 * Parameter.splittingAlpha^1.5)) * Particle.farfieldDipole; % potential gradient
