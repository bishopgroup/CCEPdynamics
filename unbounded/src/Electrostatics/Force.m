function Particle = Force(Parameter,Mesh,Particle,NearField)
%{
Compute the electrostatic free energy of a collection of charged spheres.
%}

delta = 1e-2; % step size for finite differences

Particle = Energy(Parameter,Mesh,Particle,NearField);

for n = 1:Particle.nParticle
    part_dx1 = Particle; 
    part_dx1.position(1,n) = part_dx1.position(1,n) + delta;
    part_dx1 = Energy(Parameter,Mesh,part_dx1,NearField);

    part_dx2 = Particle; 
    part_dx2.position(2,n) = part_dx2.position(2,n) + delta;
    part_dx2 = Energy(Parameter,Mesh,part_dx2,NearField);

    part_dx3 = Particle; 
    part_dx3.position(3,n) = part_dx3.position(3,n) + delta;
    part_dx3 = Energy(Parameter,Mesh,part_dx3,NearField);

    Particle.force(:,n) = [(Particle.Ues - part_dx1.Ues) / delta, ...
                           (Particle.Ues - part_dx2.Ues) / delta, ...
                           (Particle.Ues - part_dx3.Ues) / delta];
end

