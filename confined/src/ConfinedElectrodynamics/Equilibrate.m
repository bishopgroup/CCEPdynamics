function Particle = Equilibrate(Parameter,Mesh,Particle,adjmat)
%{
%}

% Compute Near Field Capacitance Matrix
Cnf = CNearField(Parameter,Particle);

storeOldNumber = Particle.nParticle;
storeOldPosition = Particle.position;
Particle.nParticle = Particle.ncParticle;
Particle.position = Particle.position(:,find(Particle.conductivityIndex==1));

qMat = zeros(Particle.nParticle);
UMat = zeros(Particle.nParticle);
Urhs = zeros(Particle.nParticle,1);
wall = diag(adjmat);
wall = wall(find(Particle.conductivityIndex==1));
Particle.adjmatCond(1:Particle.ncParticle+1:end) = wall;

[S, C] = graphconncomp(Particle.adjmatCond);
for cc = 1:S % step through connected components
    ccidx = find(C == cc);
    if any(wall(ccidx) == 2) % in contact with lower wall
        UMat(ccidx,ccidx) = eye(length(ccidx));
        Urhs(ccidx) = Particle.position(3,ccidx);
    elseif any(wall(ccidx) == 3) % in contact with upper wall
        UMat(ccidx,ccidx) = eye(length(ccidx));
        Urhs(ccidx) = Particle.position(3,ccidx) - Parameter.L3;
    else % in contact with no walls
        qMat(ccidx(1),ccidx) = 1; % sum of the charges is constant
        for i = 2:length(ccidx)
            UMat(ccidx(i),ccidx(i-1)) = 1;
            UMat(ccidx(i),ccidx(i)) = -1;
            Urhs(ccidx(i)) = Particle.position(3,ccidx(i-1)) - Particle.position(3,ccidx(i));
        end
    end
end

% Assembly Right-hand-side
b = zeros(4*Particle.nParticle,1);
b(1:Particle.nParticle) = qMat*Particle.charge + Urhs;
b(Particle.nParticle+1:end) = reshape(Particle.field,3*Particle.nParticle,1);

% Solve for Far Field Charge and Dipole Moments, farFieldCharge and farFieldDipole
x0 = [Particle.farFieldCharge; reshape(Particle.farFieldDipole,3*Particle.nParticle,1)];  % initial guess
tol = 1e-3;  % tolerance
restart = 2*Particle.nParticle;
maxint = 40;
[x,flag,relres,iter,resvec] = gmres(@(x)afun(x,Parameter,Mesh,Particle,Cnf,qMat,UMat),...
    b,restart,tol,maxint,[],[],x0);

Particle.farFieldCharge = x(1:Particle.nParticle);
Particle.farFieldDipole = reshape(x(Particle.nParticle+1:end),3,Particle.nParticle);

if flag ~= 0
    fprintf('Houston, we have a problem: flag = %d\n',flag);
    semilogy(1:length(resvec),resvec,'-o');
    xlabel('Iteration number');
    ylabel('Relative residual');
    return
end

% Complete Solution for Potential
[Ul,GradUl] = PFarFieldReal(Parameter,Particle);      % Real Space Part
[Ug,GradUg] = PFarFieldWave(Parameter,Mesh,Particle); % Wave Space Part
[Us,GradUs] = PFarFieldSelf(Parameter,Particle);      % Self-Induction Correction
Particle.potential = Ul + Ug + Us;

% Compute "Total" Dipole Moment
UE = [Particle.potential; reshape(Particle.field,3*Particle.nParticle,1)];
qp = x + Cnf*UE;
Particle.charge = qp(1:Particle.nParticle);
charge= Particle.charge;
Particle.p = reshape(qp(Particle.nParticle+1:end),3,Particle.nParticle);

% Compute Electrostatic Energy
Particle.Ues = 0.5*(sum(Particle.charge .* Particle.potential) - sum(dot(Particle.p,Particle.field,1))) ...
    - sum(Particle.position(3,:) .* Particle.charge');
Particle.position = storeOldPosition;
Particle.nParticle = storeOldNumber;

%--------------------------------------------------------------------------
function Ab = afun(x,Parameter,Mesh,Particle,Cnf,qMat,UMat)

Particle.farFieldCharge = x(1:Particle.nParticle);
Particle.farFieldDipole = reshape(x(Particle.nParticle+1:end),3,Particle.nParticle);

[Ul,GradUl] = PFarFieldReal(Parameter,Particle);      % Real Space Part
[Ug,GradUg] = PFarFieldWave(Parameter,Mesh,Particle); % Wave Space Part
[Us,GradUs] = PFarFieldSelf(Parameter,Particle);      % Self-Induction Correction

% Combine
U = Ul + Ug + Us;
GradU = GradUl + GradUg + GradUs;

% Compute "Total" Charge and Dipole Moments
UE = [U; reshape(GradU,3*Particle.nParticle,1)];
charge = Particle.farFieldCharge + Cnf(1:Particle.nParticle,:)*UE;

% Residuals 
Ab = zeros(4*Particle.nParticle,1);
Ab(1:Particle.nParticle) = qMat*charge + UMat*U;
Ab(Particle.nParticle+1:end) = reshape(GradU,3*Particle.nParticle,1);
