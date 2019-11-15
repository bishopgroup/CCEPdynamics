function [Ul,GradUl] = PFarFieldReal(Parameter,Particle)
%{
Computes "local" contribution to potential and potential gradient due to
particle charges and dipole moments.
%}
global I e3e3

% Potential and Field
Ul = zeros(Particle.nParticle,1);
GradUl = zeros(3,Particle.nParticle);

% Compute interparticle distance matrix (i with j)
xi = Particle.position; % (3 x nParticle)
xj = Particle.position; % (3 x nParticle)
dx1 = bsxfun(@minus,xi(1,:)',xj(1,:));  % xi1 - xj1  (nParticle x nParticle)
dx2 = bsxfun(@minus,xi(2,:)',xj(2,:));  % xi2 - xj2  (nParticle x nParticle)
dx3 = bsxfun(@minus,xi(3,:)',xj(3,:));  % xi3 - xj3  (nParticle x nParticle)
dx1 = dx1 - round(dx1 / Parameter.L1) * Parameter.L1; % periodic bc
dx2 = dx2 - round(dx2 / Parameter.L2) * Parameter.L2; % periodic bc
rij = sqrt(dx1.^2 + dx2.^2 + dx3.^2); % (nParticle x nParticle)

% Compute interparticle distance matrix (i with lower image of j)
dx3L = bsxfun(@minus,xi(3,:)',-xj(3,:)); % (nParticle x nParticle)
RijL = sqrt(dx1.^2 + dx2.^2 + dx3L.^2);  % (nParticle x nParticle)

% Compute interparticle distance matrix (i with upper image of j)
dx3U = bsxfun(@minus,xi(3,:)',2*Parameter.L3 - xj(3,:)); % (nParticle x nParticle)
RijU = sqrt(dx1.^2 + dx2.^2 + dx3U.^2);              % (nParticle x nParticle)

% Compute Green's functions and reflections off the walls
[G,dG,dGor,d2G] = Green(rij,Parameter.alpha);      % (nParticle x nParticle)
[GL,dGL,dGLor,d2GL] = Green(RijL,Parameter.alpha); % (nParticle x nParticle)
[GU,dGU,dGUor,d2GU] = Green(RijU,Parameter.alpha); % (nParticle x nParticle)

% Compute potential at i due to charges j (and their reflections)
Ul = (G-GL-GU)*Particle.farFieldCharge; % (nParticle x nParticle)*(nParticle x 1) --> (nParticle x 1) 

% Reflected Dipoles
pUL = (2*e3e3-I) * Particle.farFieldDipole; % (3 x 3) * (3 x nParticle) = (3 x nParticle)

for np = 1:Particle.nParticle
    % Normalized Displacement Vectors (3 x nParticle)
    eij = [dx1(np,:)./rij(np,:); dx2(np,:)./rij(np,:); dx3(np,:)./rij(np,:)];
    eij(:,np) = 0;
    eijL = [dx1(np,:)./RijL(np,:); dx2(np,:)./RijL(np,:); dx3L(np,:)./RijL(np,:)];
    eijL(:,np) = 0; 
    eijU = [dx1(np,:)./RijU(np,:); dx2(np,:)./RijU(np,:); dx3U(np,:)./RijU(np,:)];
    eijU(:,np) = 0;
    
    % Gradient at np due to charges (and their images)
    GradUl(:,np) = (bsxfun(@times, dG(np,:), eij) + ...
        bsxfun(@times, dGL(np,:), eijL) + ...
        bsxfun(@times, dGU(np,:), eijU)) * Particle.farFieldCharge; % (3 x nParticle) * (nParticle x 1)

    % Potential at np due to dipoles (and their images)
    Ul(np) = Ul(np) - sum(Particle.farFieldDipole.*eij,1) * dG(:,np) - ...
       sum(pUL.*eijL,1) *  dGL(:,np) - sum(pUL.*eijU,1) * dGU(:,np);

    % Field at np due to dipoles (and their images)
    for j = 1:Particle.nParticle
        GradUl(:,np) = GradUl(:,np) ...
             - (d2G(np,j) * eij(:,j) * eij(:,j)' + dGor(np,j) * (I - eij(:,j) * eij(:,j)')) * Particle.farFieldDipole(:,j) ...
             - (d2GL(np,j) * eijL(:,j) * eijL(:,j)' + dGLor(np,j) * (I - eijL(:,j) * eijL(:,j)')) * pUL(:,j) ...
             - (d2GU(np,j) * eijU(:,j) * eijU(:,j)' + dGUor(np,j) * (I - eijU(:,j) * eijU(:,j)')) * pUL(:,j);
    end
end

%--------------------------------------------------------------------------
function [G,dG,dGor,d2G] = Green(r,alpha)
erfc_of_r = erfc(sqrt(pi/alpha)*r);
exp_of_r =  exp(-pi*r.^2/alpha);

G = erfc_of_r ./ r;
dG = -( (2*r/sqrt(alpha)) .* exp_of_r + erfc_of_r ) ./ r.^2;
dGor = dG ./ r;
d2G = 2*( (2*r.*(alpha+pi*r.^2)/alpha^1.5) .* exp_of_r + erfc_of_r ) ./ r.^3;

G(logical(eye(size(G)))) = 0;
dG(logical(eye(size(dG)))) = 0;
dGor(logical(eye(size(dG)))) = 0;
d2G(logical(eye(size(d2G)))) = 0;