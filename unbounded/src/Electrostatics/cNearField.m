function Cnf = cNearField(Parameter,Particle,NearField)
%{
Given Particleicle charges q and field Eo, this function solves for the
Particleicle potentials U and dipole moments p.
%}
global I

Cnf = zeros(4*Particle.nParticle);

for i = 1:Particle.nParticle
    %% Loop through Particleicle pairs
    for j = i+1:Particle.nParticle
        % InterParticleicle separation (r = xj-xi)
        dx1 = Particle.position(1,j)-Particle.position(1,i);
        dx2 = Particle.position(2,j)-Particle.position(2,i);
        dx3 = Particle.position(3,j)-Particle.position(3,i);
        dx1 = dx1 - round(dx1 / Parameter.domainLength(1)) * Parameter.domainLength(1); % periodic bc
        dx2 = dx2 - round(dx2 / Parameter.domainLength(2)) * Parameter.domainLength(2); % periodic bc
        e = [dx1,dx2,dx3];
        r = norm(e);
        lnxi = real(log(r-2));
        e = e/r;
        ee = e'*e;
        
        % A matrix
        A = zeros(2);
        A(1,1) = NearField.XA11(lnxi);
        A(1,2) = NearField.XA12(lnxi);
        A(2,1) = A(1,2);
        A(2,2) = A(1,1);
        
        % Bs matrix
        Bs = zeros(2,6);
        Bs(1,1:3) = NearField.XB11(lnxi)*e;
        Bs(1,4:6) = NearField.XB12(lnxi)*e;
        Bs(2,1:3) = -Bs(1,4:6);
        Bs(2,4:6) = -Bs(1,1:3);
        
        % B matrix
        B = transpose(Bs);
        
        % C matrix
        C = zeros(6,6);
        C(1:3,1:3) = NearField.XC11(lnxi)*ee + NearField.YC11(lnxi)*(I-ee);
        C(1:3,4:6) = NearField.XC12(lnxi)*ee + NearField.YC12(lnxi)*(I-ee);
        C(4:6,1:3) = C(1:3,4:6);
        C(4:6,4:6) = C(1:3,1:3);
        
        iA = [i,j];
        Cnf(iA,iA) = A;

        iC = [Particle.nParticle + 3*(i-1) + [1:3], Particle.nParticle + 3*(j-1) + [1:3]];
        Cnf(iC,iC) = C;

        Cnf(iA,iC) = Bs;
        Cnf(iC,iA) = B;
    end
    
    %% Wall Interactions with Upper / Lower Wall
    dx3 = -2 * Particle.position(3,i);
    dx3 = dx3 - round(dx3 / (2 * Parameter.domainLength(3))) * 2 * Parameter.domainLength(3); % periodic bc
    e = [0,0,dx3];
    r = abs(dx3);
    lnxi = real(log(r-2));
    e = e/r;
    ee = e'*e;
    
    % A matrix
    A = NearField.XA11(lnxi) - NearField.XA12(lnxi);
    
    % Bs matrix
    Bs = zeros(1,3);
    Bs(1,1:3) = (NearField.XB11(lnxi) + NearField.XB12(lnxi))*e;
    
    % B matrix
    B = transpose(Bs);
    
    % C matrix
    C = (NearField.XC11(lnxi) + NearField.XC12(lnxi))*ee;
    
    iA = i;
    Cnf(iA,iA) = Cnf(iA,iA) + A;
    
    iC = [Particle.nParticle + 3*(i-1) + [1:3]];
    Cnf(iC,iC) = Cnf(iC,iC) + C;
    
    Cnf(iA,iC) = Cnf(iA,iC) + Bs;
    Cnf(iC,iA) = Cnf(iC,iA) + B;
end
