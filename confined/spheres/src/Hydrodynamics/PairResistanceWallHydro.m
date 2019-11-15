function R = PairResistanceWall(xi,e)


%% Sphere-Plane Resistance
data = load('SpherePlaneResistance.txt');
SPH_XA_ln = griddedInterpolant(log(data(:,1)),log(data(:,2)));
SPH_YA = griddedInterpolant(log(data(:,1)),data(:,3));
SPH_YB = griddedInterpolant(log(data(:,1)),data(:,4));
SPH_XC = griddedInterpolant(log(data(:,1)),data(:,5));
SPH_YC = griddedInterpolant(log(data(:,1)),data(:,6));

KroneckerDelta=eye(3);

Signature=zeros(3,3,3);
Signature(1,2,3)=1;
Signature(2,3,1)=1;
Signature(3,1,2)=1;
Signature(2,1,3)=-1;
Signature(3,2,1)=-1;
Signature(1,3,2)=-1;


% Compute pairwise resistance matrix
ee = e*e';
lnxi = real(log(xi));

% Compute Pairwise Resistance Matrix
A = exp(SPH_XA_ln(lnxi))*ee + SPH_YA(lnxi)*(KroneckerDelta - ee) - KroneckerDelta;
C = SPH_XC(lnxi)*ee + SPH_YC(lnxi)*(KroneckerDelta - ee) - (4/3)*KroneckerDelta;

DOT = Signature(:,:,1)*e(1) + Signature(:,:,2)*e(2) + Signature(:,:,3)*e(3);
B = SPH_YB(lnxi)*DOT;

% Compute Pairwise Mobility Matrix
R = [A,  B';
     B,  C];
