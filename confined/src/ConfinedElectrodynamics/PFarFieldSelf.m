function [U,GradU] = FarFieldSelf(Parameter,Particle)
%{

%}
global I e3 e3e3 puq_f1 pup_f2 pEq_f2 pEp_f3 pEp_g3

U = zeros(Particle.nParticle,1); % potential
GradU = zeros(3,Particle.nParticle); % potential gradient

Xi = Particle.position(3,:) / Parameter.L3;

for np = 1:Particle.nParticle % loop over particles
    U(np) = (1 + puq_f1(Xi(np)) / Parameter.L3) * Particle.farFieldCharge(np) + ...
        (pup_f2(Xi(np)) / Parameter.L3^2) * Particle.farFieldDipole(3,np);
    GradU(:,np) = (pEq_f2(Xi(np)) / Parameter.L3^2) * e3 * Particle.farFieldCharge(np) + ...
        ((1 + pEp_f3(Xi(np)) / Parameter.L3^3) * e3e3 + ...
         (1 + pEp_g3(Xi(np)) / Parameter.L3^3) * (I-e3e3)) * Particle.farFieldDipole(:,np);
end
    