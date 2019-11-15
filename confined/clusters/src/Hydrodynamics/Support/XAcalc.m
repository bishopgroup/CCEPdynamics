function XAcalc(s)
%{
Solves for the resistance function XA for two equally sized spheres
(D. J. Jeffrey and Y. Onishi, J. Fluid Mech., 1984, 139, 261–290.).

The input paramter s = r / a is the distance between spheres scaled by 
their common radius a.
%}

%% Compute the coefficients
N = 20; % even number
k = [0:N-1]';
fk = XAcoeff(N);

%% Compute XA11
keven = k(1:2:end);
XA11 = bsxfun(@power,2*s,-keven') * fk(1:2:end)


xi = s - 2;
XA11nf = 0.25./xi - 0.225*log(xi) - 0.0267857*xi.*log(xi) + 0.9954;
loglog(xi,XA11,'o',xi,XA11nf)

