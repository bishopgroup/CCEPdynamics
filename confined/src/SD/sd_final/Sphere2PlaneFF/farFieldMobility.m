function farFieldMobilityTensor = farFieldMobility(Parameter,Particle)
%{
Given Particle forces F and torques T
%}
KroneckerDelta=eye(3);

farFieldMobilityTensor = zeros(6*Particle.nParticle);

for i = 1:Particle.nParticle
        
    %% Loop through Particle pairs
    for j = i+1:Particle.nParticle
        % InterParticleicle separation (r = xj-xi)
        dx1 = Particle.position(j,1)-Particle.position(i,1);
        dx2 = Particle.position(j,2)-Particle.position(i,2);
        dx3 = Particle.position(j,3)-Particle.position(i,3);
        dx1 = dx1 - round(dx1 / Parameter.domainLength(1)) * Parameter.domainLength(1); % periodic bc
        dx2 = dx2 - round(dx2 / Parameter.domainLength(2)) * Parameter.domainLength(2); % periodic bc
        e = [dx1,dx2,dx3]';
        r = norm(e);
        e = e/r;
        
        Mfar = pairGrandMobilitytensor(r,e);

        iA = [3*(i-1) + [1:3],3*(j-1) + [1:3]];
        farFieldMobilityTensor(iA,iA) = farFieldMobilityTensor(iA,iA) + Mfar(1:6,1:6);

        iC = [3*Particle.nParticle + 3*(i-1) + [1:3], 3*Particle.nParticle + 3*(j-1) + [1:3]];
        farFieldMobilityTensor(iC,iC) = farFieldMobilityTensor(iC,iC) + Mfar(7:12,7:12);
        
        farFieldMobilityTensor(iA,iC) = farFieldMobilityTensor(iA,iC) + Mfar(1:6,7:12);
        farFieldMobilityTensor(iC,iA) = farFieldMobilityTensor(iC,iA) + Mfar(7:12,1:6);
        
    end

    iA = [3*(i-1) + [1:3]];
    iC = [3*Particle.nParticle + 3*(i-1) + [1:3]];
    farFieldMobilityTensor(iA,iA) = farFieldMobilityTensor(iA,iA) + KroneckerDelta;
    farFieldMobilityTensor(iC,iC) = farFieldMobilityTensor(iC,iC) + (4/3)*KroneckerDelta;

end

