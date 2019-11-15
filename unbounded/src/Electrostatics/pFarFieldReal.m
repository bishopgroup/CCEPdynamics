function [Ul,GradUl] = pFarFieldReal(Parameter,Particle)
%{
Computes "local" contribution to potential and potential gradient due to
particle charges and dipole moments.
%}
global I

% Reflect charges and dipoles
idxReal = 1:Particle.nParticle;
idxImage = Particle.nParticle+1:2*Particle.nParticle;
xi = [Particle.position, Particle.position] ; % (3 x 2*nParticle)
xi(3, idxImage) = -xi(3, idxReal); % Reflect
qff = [Particle.farfieldCharge; Particle.farfieldCharge]; % (2*nParticle x 1)
qff(idxImage) = -qff(idxReal); % Reflect
pff = [Particle.farfieldDipole, Particle.farfieldDipole]; % (3 x 2*nParticle)
pff(1:2, idxImage) = -pff(1:2, idxReal);  % Reflect

% Potential and Field
Ul = zeros(Particle.nParticle, 1);
GradUl = zeros(3, Particle.nParticle);

% Compute interparticle distance matrix (i with j)
xj = xi;
dx1 = bsxfun(@minus, xi(1,:)', xj(1,:));  % xi1 - xj1  (2*nParticle x 2*nParticle)
dx2 = bsxfun(@minus, xi(2,:)', xj(2,:));  % xi2 - xj2  (2*nParticle x 2*nParticle)
dx3 = bsxfun(@minus, xi(3,:)', xj(3,:));  % xi3 - xj3  (2*nParticle x 2*nParticle)
dx1 = dx1 - round(dx1 / Parameter.domainLength(1)) * Parameter.domainLength(1); % periodic bc
dx2 = dx2 - round(dx2 / Parameter.domainLength(2)) * Parameter.domainLength(2); % periodic bc
dx3 = dx3 - round(dx3 / (2 * Parameter.domainLength(3))) * 2 * Parameter.domainLength(3); % periodic bc
rij = sqrt(dx1.^2 + dx2.^2 + dx3.^2); % (2*nParticle x 2*nParticle)

% Compute Green's functions
[G, dG, dGor, d2G] = Green(rij, Parameter.splittingAlpha);  % (2*nParticle x 2*nParticle)

% Compute potential at i due to charges j (and their reflections)
Ul = G(idxReal, :) * qff; % (nParticle x 2*nParticle)*(2*nParticle x 1) --> (nParticle x 1) 

for np = 1:Particle.nParticle
    % Normalized Displacement Vectors (3 x 2*nParticle)
    eij = [dx1(np,:) ./ rij(np,:); 
           dx2(np,:) ./ rij(np,:); 
           dx3(np,:) ./ rij(np,:)];
    eij(:,np) = 0;
    
    % Gradient at np due to charges (and their images)
    GradUl(:,np) = bsxfun(@times, dG(np,:), eij) * qff; % (3 x 2*nParticle) * (2*nParticle x 1)

    % Potential at np due to dipoles (and their images)
    Ul(np) = Ul(np) - sum(pff .* eij, 1) * dG(:,np);

    % Field at np due to dipoles (and their images)
    for j = 1:2*Particle.nParticle
        GradUl(:,np) = GradUl(:,np) - (d2G(np,j) * eij(:,j) * eij(:,j)' + ...
            dGor(np,j) * (I - eij(:,j) * eij(:,j)')) * pff(:,j);
    end
end

%--------------------------------------------------------------------------
function [G,dG,dGor,d2G] = Green(r, splittingAlpha)
erfc_of_r = erfc( sqrt( pi / splittingAlpha ) * r);
exp_of_r =  exp(-pi * r.^2 / splittingAlpha);

G = erfc_of_r ./ r;
dG = -( ( 2 * r / sqrt(splittingAlpha)) .* exp_of_r + erfc_of_r ) ./ r.^2;
dGor = dG ./ r;
d2G = 2 * ( (2 * r .* (splittingAlpha + pi * r.^2) / splittingAlpha^1.5) .* exp_of_r + erfc_of_r ) ./ r.^3;

G(logical(eye(size(G)))) = 0;
dG(logical(eye(size(dG)))) = 0;
dGor(logical(eye(size(dG)))) = 0;
d2G(logical(eye(size(d2G)))) = 0;