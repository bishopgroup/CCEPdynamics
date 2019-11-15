function [Particle, Parameter] = mapTo2L3(Particle, Parameter)
%% Maps the confined configuration to a periodic domain in all direction with electrode separation 2*L3

    %% Mapping variables for emulating confinement using 2L3 all direction periodic lattice 
    ParticlePositions = [Particle.position;Particle.position];
    ParticlePositions(1:Particle.nParticle,3) = Parameter.L3 + ParticlePositions(1:Particle.nParticle,3);
    ParticlePositions(1+Particle.nParticle:end,3) = Parameter.L3 - ParticlePositions(1+Particle.nParticle:end,3);
    Particle.position = ParticlePositions;
    
    Particle.nParticle = 2*Particle.nParticle;
    Parameter.domainLength(3) = 2*Parameter.domainLength(3);
    
    
end