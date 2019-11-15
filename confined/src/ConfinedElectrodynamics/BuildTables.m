function BuildTables()
global I e3 e3e3 puq_f1 pup_f2 pEq_f2 pEp_f3 pEp_g3 ...
    XA11 XA12 XB11 XB12 XC11 XC12 YC11 YC12 ...
    SSH_X11A_ln SSH_X12A_ln SSH_Y11A SSH_Y12A SSH_Y11B SSH_Y12B SSH_X11C SSH_X12C SSH_Y11C SSH_Y12C ...
    KroneckerDelta Signature ...
    SPH_XA_ln SPH_YA SPH_YB SPH_XC SPH_YC

%% Tensors
I = eye(3); % Identity Tensor
e3 = zeros(3,1); e3(3) = 1;
e3e3 = zeros(3); e3e3(3,3) = 1;


%% Self-Induction Coefficients
puq_f = load('puq_f.txt');
puq_f1 = griddedInterpolant(puq_f(:,1),puq_f(:,2));

pup_f = load('pup_f.txt');
pup_f2 = griddedInterpolant(pup_f(:,1),pup_f(:,2));

pEq_f = load('pEq_f.txt');
pEq_f2 = griddedInterpolant(pEq_f(:,1),pEq_f(:,2));

pEp_f = load('pEp_f.txt');
pEp_f3 = griddedInterpolant(pEp_f(:,1),pEp_f(:,2));

pEp_g = load('pEp_g.txt');
pEp_g3 = griddedInterpolant(pEp_g(:,1),pEp_g(:,2));


%% Two Sphere Capacitance 
data = load('TwoSphereCapacitance.txt');
XA11 = griddedInterpolant(log(data(:,1)-2),data(:,2));
XA12 = griddedInterpolant(log(data(:,1)-2),data(:,3));
XB11 = griddedInterpolant(log(data(:,1)-2),data(:,4));
XB12 = griddedInterpolant(log(data(:,1)-2),data(:,5));
XC11 = griddedInterpolant(log(data(:,1)-2),data(:,6));
XC12 = griddedInterpolant(log(data(:,1)-2),data(:,7));
YC11 = griddedInterpolant(log(data(:,1)-2),data(:,8));
YC12 = griddedInterpolant(log(data(:,1)-2),data(:,9));

%% Sphere-Plane Resistance
data = load('SpherePlaneResistance.txt');
SPH_XA_ln = griddedInterpolant(log(data(:,1)),log(data(:,2)));
SPH_YA = griddedInterpolant(log(data(:,1)),data(:,3));
SPH_YB = griddedInterpolant(log(data(:,1)),data(:,4));
SPH_XC = griddedInterpolant(log(data(:,1)),data(:,5));
SPH_YC = griddedInterpolant(log(data(:,1)),data(:,6));

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
