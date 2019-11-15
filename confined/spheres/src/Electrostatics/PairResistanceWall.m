function R = PairResistanceWall(xi,e)
global SPH_XA_ln SPH_YA SPH_YB SPH_XC SPH_YC KroneckerDelta Signature

% Compute pairwise resistance matrix
ee = e*e';
lnxi = log(xi);

% Compute Pairwise Resistance Matrix
A = exp(SPH_XA_ln(lnxi))*ee + SPH_YA(lnxi)*(KroneckerDelta - ee) - KroneckerDelta;
C = SPH_XC(lnxi)*ee + SPH_YC(lnxi)*(KroneckerDelta - ee) - (4/3)*KroneckerDelta;

DOT = Signature(:,:,1)*e(1) + Signature(:,:,2)*e(2) + Signature(:,:,3)*e(3);
B = SPH_YB(lnxi)*DOT;

% Compute Pairwise Mobility Matrix
R = [A,  B';
     B,  C];