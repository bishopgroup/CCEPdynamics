function Particle = Velocity(Parameter,Particle)

Particle.velocity = Particle.force;

Rnf = rNearField(Parameter,Particle);

FT = [reshape(Particle.force,3*Particle.nParticle,1) ; zeros(3*Particle.nParticle,1)];
VW = Rnf \ FT;

Particle.velocity = reshape(VW(1:3*Particle.nParticle),3,Particle.nParticle);
