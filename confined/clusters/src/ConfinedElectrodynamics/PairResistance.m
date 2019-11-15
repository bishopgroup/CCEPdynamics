function R = PairResistance(s,e)
global SSH_X11A_ln SSH_X12A_ln SSH_Y11A SSH_Y12A SSH_Y11B SSH_Y12B ...
    SSH_X11C SSH_X12C SSH_Y11C SSH_Y12C KroneckerDelta Signature

% Compute pairwise resistance matrix
ee = e*e';
lnxi = log(s-2);

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
