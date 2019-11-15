function [velocity, velocityGradient, velocityLaplacian, ...
    velocityLaplacianGradient] = mFarFieldSelf(Parameter, Particle)
%--------------------------------------------------------------------------
% pFarFieldSelf - computes the self-contribution to the 'global' potential.
%
% [potential, potentialGradient] = pFarFieldSelf(Parameter, Particle)
% computes the self-contribution to the global potential (due to Gaussian-
% filtered charges) and the potential gradient (due to Gaussian-filtered 
% dipoles) at particle centers.
%
% potential: (nParticle x 1) array of self-potentials.
%
% potentialGradient: (nParticle x 3) array of self-potential gradients.
%--------------------------------------------------------------------------

%% self-velocity due to forces, dipoles, quadrupoles, octupoles
velocity = (3 / sqrt(Parameter.splittingAlpha)) * Particle.farFieldForce ...
    - (5 * pi / Parameter.splittingAlpha^1.5) * Particle.farFieldForce / 3;

%% self-velocity gradient due to forces, dipoles, quadrupoles, octupoles
velocityGradient = zeros(Particle.nParticle, 3, 3);
velocityLaplacianGradient = zeros(Particle.nParticle, 3, 3);
I = eye(3);
for m = 1:3
    for i = 1:3
        for n = 1:3
            for j = 1:3
                
                Aminj = (2 * I(m,i) * I(i,j) * I(n,j) - ...
                    I(m,i) * I(n,j) * (1-I(i,j)) + ...
                    4 * (1-I(m,i)) * (1-I(n,j)) * I(i,j) * I(m,n) -...
                    (1-I(m,i)) * (1-I(n,j)) * I(i,n) * I(j,m));
                
                velocityGradient(:,m,i) = velocityGradient(:,m,i) ...
                    + (pi/Parameter.splittingAlpha^1.5) * Aminj * ...
                    Particle.farFieldDipoleHydro(:,n,j) ...
                    - (7*pi^2 / Parameter.splittingAlpha^2.5) * Aminj * ...
                    Particle.farFieldDipoleHydro(:,n,j) * (3/5); 
                
                velocityLaplacianGradient(:,m,i) = velocityLaplacianGradient(:,m,i) ...
                    - (42*pi^2 / (5*Parameter.splittingAlpha^2.5)) * Aminj * ...
                    Particle.farFieldDipoleHydro(:,n,j) ...
                    + (18*pi / Parameter.splittingAlpha^3.5) * Aminj * ...
                    Particle.farFieldDipoleHydro(:,n,j) * (3/5);       
            end
        end
    end
end


%% self-velocity laplacian due to forces, dipoles, quadrupoles, octupoles
velocityLaplacian = (-10 * pi / Parameter.splittingAlpha^1.5) * Particle.farFieldForce ...
    + (42 * pi^2 / Parameter.splittingAlpha^2.5) * Particle.farFieldForce / 3;
