function [velocity, velocityGradient, velocityLaplacian, ...
    velocityGradientLaplacian] = mFarFieldReal(Parameter, Particle)
%--------------------------------------------------------------------------
% pFarFieldReal - computes "local" contribution to velocity and velocity 
%                 gradient due to particle forces and dipoles.
%
% [localPotential,localPotentialGradient] = ...
%                                     pFarFieldReal(Parameter, Particle)
%--------------------------------------------------------------------------

Np = Particle.nParticle;
I = eye(3);

%% Initialize velocity arrays
velocity = zeros(Np,3);
velocityLaplacian = zeros(Np,3);
velocityGradient = zeros(Np,3,3);
velocityGradientLaplacian = zeros(Np,3,3);


%% interparticle distance matrix
delta = zeros(Np,Np,3);

% compute components of the separation vector
delta(:,:,1) = bsxfun(@minus, Particle.position(:,1), ...
    Particle.position(:,1)');
delta(:,:,2) = bsxfun(@minus, Particle.position(:,2), ...
    Particle.position(:,2)');
delta(:,:,3) = bsxfun(@minus, Particle.position(:,3), ...
    Particle.position(:,3)');

% account for periodic boundaries
delta(:,:,1) = delta(:,:,1) - round(delta(:,:,1) / Parameter.domainLength(1)) ...
    * Parameter.domainLength(1);
delta(:,:,2) = delta(:,:,2) - round(delta(:,:,2) / Parameter.domainLength(2)) ...
    * Parameter.domainLength(2);
delta(:,:,3) = delta(:,:,3) - round(delta(:,:,3) / Parameter.domainLength(3)) ...
    * Parameter.domainLength(3);

% Compute interparticle distance matrix (Brute Force)
r = sqrt(delta(:,:,1).^2 + delta(:,:,2).^2 + delta(:,:,3).^2);
rho = r / sqrt(Parameter.splittingAlpha / pi);
rhoSquared = rho.^2;

diagonal = logical(eye(size(r)));
r(diagonal) = Inf;

% Compute unit displacement vectors
e = zeros(Np,Np,3);
e(:,:,1) = delta(:,:,1) ./ r;
e(:,:,2) = delta(:,:,2) ./ r;
e(:,:,3) = delta(:,:,3) ./ r;


%% Compute incomplete gamma functions
phi0 = incompleteGamma(-0.5,rhoSquared);
phi1 = incompleteGamma(0.5,rhoSquared);
phi2 = incompleteGamma(1.5,rhoSquared);
phi3 = incompleteGamma(2.5,rhoSquared);
phi4 = incompleteGamma(3.5,rhoSquared);

phi0(diagonal) = 0;
phi1(diagonal) = 0;
phi2(diagonal) = 0;
phi3(diagonal) = 0;
phi4(diagonal) = 0;



%% Contributions due to forces
% velocity and velocity laplacian
for i = 1:3
    for j = 1:3
        velocity(:,i) = velocity(:,i) ...
           + (1.5 / sqrt(Parameter.splittingAlpha)) ...
           * (-rho.^2 .* phi1 .* (I(i,j) - e(:,:,i) .* e(:,:,j)) ...
           + phi0 * I(i,j)) * Particle.force(:,j);

       velocityLaplacian(:,i) = velocityLaplacian(:,i) ...
           + (1.5 * pi / Parameter.splittingAlpha^1.5) ...
           * (-4 * rho.^4 .* phi3 .* (I(i,j) - e(:,:,i) .* e(:,:,j)) ...
           + 2 * rho.^2 .* phi2 .* (9 * I(i,j) - 7 * e(:,:,i) .* e(:,:,j)) ...
           - 10 * phi1 * I(i,j)) * Particle.force(:,j);
    end
end

% velocity gradient and gradient laplacian
for i = 1:3
    for j = 1:3
        for k = 1:3
            velocityGradient(:,k,i) = velocityGradient(:,k,i) + ...
                (1.5 * sqrt(pi) / Parameter.splittingAlpha) * ...
                (-phi1 .* (4*I(i,j)*e(:,:,k) - I(j,k)*e(:,:,i) - I(i,k)*e(:,:,j))...
                + 2 * rho.^3 .* phi2 .* (I(i,j)*e(:,:,k) - e(:,:,i) .* e(:,:,j) .* e(:,:,k))) * Particle.force(:,j);

%             velocityGradientLaplacian(:,i,k) = velocityGradientLaplacian(:,i,k) + ...
%                 (1.5 * pi^1.5 / Parameter.splittingAlpha^2) * ...
%                 (8 * rho.^5 .* phi4 .* (
%             
%             4*I(i,j)*e(:,:,k) - I(j,k)*e(:,:,i) - I(i,k)*e(:,:,j))...
%                 + 2 * rho.^3 .* phi2 .* (I(i,j)*e(:,:,k) - e(:,:,i) .* e(:,:,j) .* e(:,:,k))) * Particle.force(:,j);
        end
    end
end

%velocityLaplacian


%% Compute velocity due to forces, dipoles, quadrupoles, octupoles

% due to forces, Gij*fj 
% %   Gij = G1*dij - G2*(dij - ei*ej)
% for i = 1:3
%     for j = 1:3
%         localVelocity(:,i) = localVelocity(:,i) ...
%             + (Greens.G1 * I(i,j)) * Particle.force(:,j) ...
%             - (Greens.G2 .* (I(i,j) - e(:,:,i) .* e(:,:,j))) * Particle.force(:,1);
%     end
% end
% 
% % due to dipoles, DnGij * Dnj 
% %   DnGij = (dG1 - dG2)*en*dij + dG2*en*ei*ej + *ei*ej) )
% for i = 1:3
%     for n = 1:3
%         for j = 1:3
%             localVelocity(:,i) = localVelocity(:,i) ...
%                 + ((Greens.dG1 - Greens.dG2) .* e(:,:,n) * I(i,j)) * Particle.force(:,j) ...
%                 + (Greens.G2 .* (I(i,j) - e(:,:,i) .* e(:,:,j))) * Particle.force(:,1);
%         end
%     end
% end





%--------------------------------------------------------------------------
function phi = incompleteGamma(v,y)

phi = gamma(1+v) * gammainc(y,1+v,'upper') .* y.^(-1-v);


