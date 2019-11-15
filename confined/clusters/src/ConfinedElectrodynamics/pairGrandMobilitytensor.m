function Mfar = pairGrandMobilitytensor(s,e)


%% Two Sphere scalar mobility functions
X11A = 1;
X22A = 1;
X12A = 1.5*s^-1 - s^-3;
X21A = 1.5*s^-1 - s^-3;

Y11A = 1;
Y22A = 1;
Y12A = (3/4)*s^-1 + 0.5*s^-3;
Y21A = (3/4)*s^-1 + 0.5*s^-3;

Y11B = 0;
Y22B = 0;
Y12B = -(3/4)*s^-2;
Y21B = (3/4)*s^-2;

X11C = (3/4);
X22C = (3/4);
X12C = (3/4)*s^-3;
X21C = (3/4)*s^-3;

Y11C = (3/4);
Y22C = (3/4);
Y12C = -(3/8)*s^-3;
Y21C = -(3/8)*s^-3;

KroneckerDelta=eye(3);

Signature=zeros(3,3,3);
Signature(1,2,3)=1;
Signature(2,3,1)=1;
Signature(3,1,2)=1;
Signature(2,1,3)=-1;
Signature(3,2,1)=-1;
Signature(1,3,2)=-1;

ee = e*e';

% Compute Pairwise mobolity Matrix
A11 = X11A*ee + Y11A*(KroneckerDelta - ee) - KroneckerDelta;
A12 = X12A*ee + Y12A*(KroneckerDelta - ee);
A22 = X22A*ee + Y22A*(KroneckerDelta - ee) - KroneckerDelta;
A21 = X21A*ee + Y21A*(KroneckerDelta - ee);
C11 = X11C*ee + Y11C*(KroneckerDelta - ee) - 0.75*KroneckerDelta;
C12 = X12C*ee + Y12C*(KroneckerDelta - ee);
C22 = X22C*ee + Y22C*(KroneckerDelta - ee) - 0.75*KroneckerDelta;
C21 = X21C*ee + Y21C*(KroneckerDelta - ee);

DOT = Signature(:,:,1)*e(1) + Signature(:,:,2)*e(2) + Signature(:,:,3)*e(3);
B11 = Y11B*DOT;
B12 = Y12B*DOT;
B22 = Y22B*DOT;
B21 = Y21B*DOT;

% Compute Pairwise Mobility Matrix
Mfar = [A11,    A12,  B11',  B12';
        A21,    A22,  B21',  B22';
        B11,    B12,  C11,   C12;
        B21,    B22,  C21,   C22];
