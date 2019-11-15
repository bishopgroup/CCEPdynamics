function Particle = calculateVelocity(Parameter,Particle)
global Signature

% Compute the resistance matrix
Rnf = RNearField(Parameter,Particle);
Rff = inv(farFieldMobility(Parameter,Particle));
ResistanceMatrix = Rnf;
%ResistanceMatrix = Rff;
%ResistanceMatrix
%pause;

% Compute "A matrix"
bodyFixedPosition = bsxfun(@minus,Particle.position,Particle.clusterTranslation);
z = zeros(3);
A = [repmat(eye(3),1,Particle.nParticle),repmat(zeros(3),1,Particle.nParticle);repmat(zeros(3),1,Particle.nParticle),repmat(eye(3),1,Particle.nParticle)];
     
for iPart= 1:Particle.nParticle
    A(4:end,3*(iPart-1)+1:3*iPart) = -(Signature(:,:,1)*bodyFixedPosition(1,iPart) + Signature(:,:,2)*bodyFixedPosition(2,iPart) ...
        + Signature(:,:,3)*bodyFixedPosition(3,iPart));
end 
%A
%pause;

ARB = A*ResistanceMatrix*A';% from stokesian dynamics for swimming paper
%ARB
%pause;

shiftR = (sum(Particle.position, 2)/Particle.nParticle) - (sum(Particle.position(:,find(Particle.conductivityIndex==1)), 2)/Particle.ncParticle);
Particle.torque = Particle.torque - cross(shiftR, Particle.force);
VW = ARB \ [Particle.force; Particle.torque];

Particle.velocity = VW(1:3);
Particle.angularVelocity = VW(4:6);
