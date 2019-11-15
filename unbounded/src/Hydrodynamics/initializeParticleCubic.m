function [Parameter, Particle] = initializeParticleCubic(Parameter)
%%-------------------------------------------------------------------------
% initializeParticleCubic - initialize particle positions.
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

%% Unit cell

% number of particles in the unit cell
nParticle = 1;

% unit cell length
unitCellLength = (4*pi / (3*Parameter.volumeFraction))^(1/3);

% particle positions in unit cell
zDirectionShift = unitCellLength/2;
unitCellPosition = zeros(1, 3);
unitCellPosition(1,:) = [0,0,zDirectionShift]; 


%% Simulation domain

% update Parameter structure
Parameter.domainLength = Parameter.nUnitCell * unitCellLength* ones(3,1); 
Parameter.nParticle = Parameter.nUnitCell^3 * nParticle;
Parameter.splittingAlpha = (1/4) * Parameter.domainLength(3).^2; 

% update Particle structure
Particle.nParticle = Parameter.nUnitCell^3 * nParticle;
Particle.position = zeros(Parameter.nUnitCell, 3);

% iterate through all unit cells and add particles to position
particleCount = 1;

for iUnitCell = 1:Parameter.nUnitCell
    % x-position of unit cell
    xUnitCell = (iUnitCell-1)*unitCellLength;

    for jUnitCell = 1:Parameter.nUnitCell
        % y-position of unit cell
        yUnitCell = (jUnitCell-1)*unitCellLength;
        
        for kUnitCell = 1:Parameter.nUnitCell
            % z-position of unit cell
            zUnitCell = (kUnitCell-1)*unitCellLength;
            
            % particle positions
            Particle.position(particleCount:particleCount,:) ...
                = bsxfun(@plus, unitCellPosition, ...
                [xUnitCell,yUnitCell,zUnitCell] );
                        
            % update counter
            particleCount = particleCount + 1;
        end
    end
end

alpha = Parameter.splittingAlpha;
x = linspace(4*pi/sqrt(alpha),100,1000);
G1 = 1.5 * gammainc(x,0.5,'upper') * gamma(0.5) ./ x.^0.5 / sqrt(alpha);
G2 = -1.5 * gammainc(x,1.5,'upper') * gamma(1.5) ./ x.^0.5 / sqrt(alpha);
Parameter.G1interp = griddedInterpolant(x,G1);
Parameter.G2interp = griddedInterpolant(x,G2);

