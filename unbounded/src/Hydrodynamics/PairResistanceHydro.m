function R = PairResistance(s,e)


%% Two Sphere Resistance
data = load('SphereSphereResistance.txt');
SSH_X11A_ln = griddedInterpolant(log(data(:,1)-2),log(data(:,2)));
SSH_X12A_ln = griddedInterpolant(log(data(:,1)-2),log(-data(:,3)));
SSH_Y11A = griddedInterpolant(log(data(:,1)-2),data(:,4));
SSH_Y12A = griddedInterpolant(log(data(:,1)-2),data(:,5));
SSH_Y11B = griddedInterpolant(log(data(:,1)-2),data(:,6));
SSH_Y12B = griddedInterpolant(log(data(:,1)-2),data(:,7));
SSH_X11C = griddedInterpolant(log(data(:,1)-2),data(:,8));
SSH_X12C = griddedInterpolant(log(data(:,1)-2),data(:,9));
SSH_Y11C = griddedInterpolant(log(data(:,1)-2),data(:,10));
SSH_Y12C = griddedInterpolant(log(data(:,1)-2),data(:,11));

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
lnxi = real(log(s-2));

% Compute Pairwise Resistance Matrix
A11 = exp(SSH_X11A_ln(lnxi))*ee + SSH_Y11A(lnxi)*(KroneckerDelta - ee) - KroneckerDelta;
A12 = -exp(SSH_X12A_ln(lnxi))*ee + SSH_Y12A(lnxi)*(KroneckerDelta - ee);
C11 = SSH_X11C(lnxi)*ee + SSH_Y11C(lnxi)*(KroneckerDelta - ee) - 4/3*KroneckerDelta;
C12 = SSH_X12C(lnxi)*ee + SSH_Y12C(lnxi)*(KroneckerDelta - ee);

DOT = Signature(:,:,1)*e(1) + Signature(:,:,2)*e(2) + Signature(:,:,3)*e(3);
B11 = SSH_Y11B(lnxi)*DOT;
B12 = SSH_Y12B(lnxi)*DOT;

% Compute Pairwise Mobility Matrix
R = [A11,    A12,  B11',  -B12';
     A12',   A11, B12', -B11';
     B11, B12,  C11,  C12;
     -B12, -B11,  C12', C11];
