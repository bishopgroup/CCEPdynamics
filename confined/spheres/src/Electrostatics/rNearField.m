function Rnf = rNearField(Parameter,Particle)
%{
Given particle forces F and torques T
%}
global KroneckerDelta

Rnf = zeros(6*Particle.nParticle);

for i = 1:Particle.nParticle
        
    %% Loop through particle pairs
    for j = i+1:Particle.nParticle
        % Interparticle separation (r = xj-xi)
        dx1 = Particle.position(1,j)-Particle.position(1,i);
        dx2 = Particle.position(2,j)-Particle.position(2,i);
        dx3 = Particle.position(3,j)-Particle.position(3,i);
        dx1 = dx1 - round(dx1 / Parameter.domainLength(1)) * Parameter.domainLength(1); % periodic bc
        dx2 = dx2 - round(dx2 / Parameter.domainLength(2)) * Parameter.domainLength(2); % periodic bc
        e = [dx1,dx2,dx3]';
        r = norm(e);
        lnr = log(r);
        e = e/r;
        ee = e'*e;
        
        R = PairResistance(r,e);

        iA = [3*(i-1) + [1:3],3*(j-1) + [1:3]];
        Rnf(iA,iA) = Rnf(iA,iA) + R(1:6,1:6);

        iC = [3*Particle.nParticle + 3*(i-1) + [1:3], 3*Particle.nParticle + 3*(j-1) + [1:3]];
        Rnf(iC,iC) = Rnf(iC,iC) + R(7:12,7:12);
        
        Rnf(iA,iC) = Rnf(iA,iC) + R(1:6,7:12);
        Rnf(iC,iA) = Rnf(iC,iA) + R(7:12,1:6);
        
    end

    iA = [3*(i-1) + [1:3]];
    iC = [3*Particle.nParticle + 3*(i-1) + [1:3]];
    Rnf(iA,iA) = Rnf(iA,iA) + KroneckerDelta;
    Rnf(iC,iC) = Rnf(iC,iC) + (4/3)*KroneckerDelta;
    
    %% Wall Interactions with Lower Wall 
    e = [0,0,1]';
    xi = Particle.position(3,i) - 1;
    R = PairResistanceWall(xi,e);

    Rnf(iA,iA) = Rnf(iA,iA) + R(1:3,1:3);
    Rnf(iC,iC) = Rnf(iC,iC) + R(4:6,4:6);
    Rnf(iA,iC) = Rnf(iA,iC) + R(1:3,4:6);
    Rnf(iC,iA) = Rnf(iC,iA) + R(4:6,1:3);
    
    
    %% Wall Interactions with Upper Wall
    e = [0,0,-1]';
    xi = Parameter.domainLength(3) - 1 - Particle.position(3,i);
    R = PairResistanceWall(xi,e); 
   
    Rnf(iA,iA) = Rnf(iA,iA) + R(1:3,1:3);
    Rnf(iC,iC) = Rnf(iC,iC) + R(4:6,4:6);
    Rnf(iA,iC) = Rnf(iA,iC) + R(1:3,4:6);
    Rnf(iC,iA) = Rnf(iC,iA) + R(4:6,1:3);
end