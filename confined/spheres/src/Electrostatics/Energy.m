function [Particle] = Energy(Parameter, Mesh, Particle, NearField)
%{
Compute the electrostatic free energy of a collection of charged spheres.
%}

% Compute Near Field Capacitance Matrix
Cnf = cNearField(Parameter, Particle, NearField);

% Assembly Right-hand-side
b = zeros(4 * Particle.nParticle, 1);
b(1:Particle.nParticle) = Particle.charge;
b(Particle.nParticle+1:end) = reshape(Particle.field, 3 * Particle.nParticle, 1);

% Solve for Far Field Charge and Dipole Moments, farfieldCharge and farfieldDipole
x0 = [Particle.farfieldCharge; reshape(Particle.farfieldDipole, 3 * Particle.nParticle, 1)];  % initial guess
tol = 1e-4;  % tolerance

[x,flag,relres,iter,resvec] = gmres(@(x)afun(x, Parameter, Mesh, Particle, Cnf), b,...
    [],tol,[],[],[],x0);

Particle.farfieldCharge = x(1:Particle.nParticle);
Particle.farfieldDipole = reshape(x(Particle.nParticle+1:end),3,Particle.nParticle);

if flag ~= 0
    semilogy(1:length(resvec),resvec,'-o');
    xlabel('Iteration number');
    ylabel('Relative residual');
    error(sprintf('Houston, we have a problem: flag = %d\n',flag));
end

% Complete Solution for Potential
[Ul,GradUl] = pFarFieldReal(Parameter, Particle);      % Real Space Part
[Ug,GradUg] = pFarFieldWave(Parameter, Mesh, Particle); % Wave Space Part
[Us,GradUs] = pFarFieldSelf(Parameter, Particle);      % Self-Induction Correction

Particle.potential = Ul + Ug + Us;

% Compute "Total" Dipole Moment
UE = [Particle.potential; reshape(Particle.field,3*Particle.nParticle,1)];
qp = x + Cnf*UE;
Particle.dipole = reshape(qp(Particle.nParticle+1:end),3,Particle.nParticle);

% Compute Electrostatic Energy
Particle.Ues = 0.5*(sum(Particle.charge .* Particle.potential) - sum(dot(Particle.dipole,Particle.field,1))) ...
     - sum(Particle.position(3,:) .* Particle.charge');


%--------------------------------------------------------------------------
function Ab = afun(x,Parameter,Mesh,Particle,Cnf)

Particle.farfieldCharge = x(1:Particle.nParticle);
Particle.farfieldDipole = reshape(x(Particle.nParticle+1:end),3,Particle.nParticle);

[Ul, GradUl] = pFarFieldReal(Parameter,Particle);      % Real Space Part
[Ug, GradUg] = pFarFieldWave(Parameter,Mesh,Particle); % Wave Space Part
[Us, GradUs] = pFarFieldSelf(Parameter,Particle);      % Self-Induction Correction

% Combine
potential = Ul + Ug + Us;
GradU = GradUl + GradUg + GradUs;

% Compute "Total" Charge and Dipole Moments
UE = [potential; reshape(GradU,3*Particle.nParticle,1)];
charge = Particle.farfieldCharge + Cnf(1:Particle.nParticle,:)*UE;

% Residuals 
Ab = zeros(4*Particle.nParticle,1);
Ab(1:Particle.nParticle) = charge;
Ab(Particle.nParticle+1:end) = reshape(GradU,3*Particle.nParticle,1);
