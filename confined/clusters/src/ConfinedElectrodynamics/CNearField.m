function Cnf = CNearField(Parameter,Particle)
%{
Given particle charges q and field Eo, this function solves for the
particle potentials U and dipole moments p.
%}
global  XA11 XA12 XB11 XB12 XC11 XC12 YC11 YC12 I

Cnf = zeros(4*Particle.ncParticle);
Particle.position = Particle.position(:,find(Particle.conductivityIndex==1));

for i = 1:Particle.ncParticle
    %% Loop through particle pairs
    for j = i+1:Particle.ncParticle
        % Interparticle separation (r = xj-xi)
        dx1 = Particle.position(1,j)-Particle.position(1,i);
        dx2 = Particle.position(2,j)-Particle.position(2,i);
        dx3 = Particle.position(3,j)-Particle.position(3,i);
        dx1 = dx1 - round(dx1 / Parameter.L1) * Parameter.L1; % periodic bc
        dx2 = dx2 - round(dx2 / Parameter.L2) * Parameter.L2; % periodic bc
        e = [dx1,dx2,dx3];
        r = norm(e);
        lnxi = log(r-2);
        e = e/r;
        ee = e'*e;
        
        % A matrix
        A = zeros(2);
        A(1,1) = XA11(lnxi);
        A(1,2) = XA12(lnxi);
        A(2,1) = A(1,2);
        A(2,2) = A(1,1);
        
        % Bs matrix
        Bs = zeros(2,6);
        Bs(1,1:3) = XB11(lnxi)*e;
        Bs(1,4:6) = XB12(lnxi)*e;
        Bs(2,1:3) = -Bs(1,4:6);
        Bs(2,4:6) = -Bs(1,1:3);
        
        % B matrix
        B = transpose(Bs);
        
        % C matrix
        C = zeros(6,6);
        C(1:3,1:3) = XC11(lnxi)*ee + YC11(lnxi)*(I-ee);
        C(1:3,4:6) = XC12(lnxi)*ee + YC12(lnxi)*(I-ee);
        C(4:6,1:3) = C(1:3,4:6);
        C(4:6,4:6) = C(1:3,1:3);
        
        iA = [i,j];
        Cnf(iA,iA) = A;

        iC = [Particle.ncParticle + 3*(i-1) + [1:3], Particle.ncParticle + 3*(j-1) + [1:3]];
        Cnf(iC,iC) = C;

        Cnf(iA,iC) = Bs;
        Cnf(iC,iA) = B;
    end
    
    %% Wall Interactions with Lower Wall (image approach)
    dx3 = -2*Particle.position(3,i);
    e = [0,0,dx3];
    r = abs(dx3);
    lnxi = real(log(r-2));
    e = e/r;
    ee = e'*e;
    
    % A matrix
    A = XA11(lnxi) - XA12(lnxi);
    
    % Bs matrix
    Bs = zeros(1,3);
    Bs(1,1:3) = (XB11(lnxi) + XB12(lnxi))*e;
    
    % B matrix
    B = transpose(Bs);
    
    % C matrix
    C = (XC11(lnxi) + XC12(lnxi))*ee;
    
    iA = i;
    Cnf(iA,iA) = Cnf(iA,iA) + A;
    
    iC = [Particle.ncParticle + 3*(i-1) + [1:3]];
    Cnf(iC,iC) = Cnf(iC,iC) + C;
    
    Cnf(iA,iC) = Cnf(iA,iC) + Bs;
    Cnf(iC,iA) = Cnf(iC,iA) + B;

    
    %% Wall Interactions with Upper Wall (image approach)
    dx3 = 2*(Parameter.L3-Particle.position(3,i));
    e = [0,0,dx3];
    r = abs(dx3);
    lnxi = real(log(r-2));
    e = e/r;
    ee = e'*e;

    % A matrix
    A = XA11(lnxi) - XA12(lnxi);
    
    % Bs matrix
    Bs = zeros(1,3);
    Bs(1,1:3) = (XB11(lnxi) + XB12(lnxi))*e;
    
    % B matrix
    B = transpose(Bs);
    
    % C matrix
    C = (XC11(lnxi) + XC12(lnxi))*ee;
    
    iA = i;
    Cnf(iA,iA) = Cnf(iA,iA) + A;
    
    iC = [Particle.ncParticle + 3*(i-1) + [1:3]];
    Cnf(iC,iC) = Cnf(iC,iC) + C;
    
    Cnf(iA,iC) = Cnf(iA,iC) + Bs;
    Cnf(iC,iA) = Cnf(iC,iA) + B;
end
