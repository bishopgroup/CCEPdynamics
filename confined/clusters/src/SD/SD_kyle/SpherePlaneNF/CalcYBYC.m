function [YB,YC] = CalcYBYC(xi)
%{
Solves for the resistance function YC for a sphere rotating parallel
to a plane surface [1 W.R. Dean and M.E. O’Neill, Mathematika 10, 13
(1963)].

The input paramter is xi = (h-a)/a is the dimensionless separation between
the surface of the sphere and the plane.
%}

a = acosh(1+xi);  % bispherical parameter
N = 500; % Number of terms in the sums

%% Solve for the An Coefficients (equation 27)
nA = [1:N]';
kA = (nA+0.5).*coth((nA+0.5)*a) - coth(a);

% Build sparse matrix
Amat = spalloc(N,N,3*N);

% left boundary
n = 1;
k0 = (0.5).*coth(0.5*a) - coth(a);
Amat(n,n) = ((2*n-1)*k0 - (2*n-3)*kA(n))*(-n / (2*n+1)) ...
    - ((2*n+5)*kA(n) - (2*n+3)*kA(n+1))*((n+1) / (2*n+1));
Amat(n,n+1) = -((2*n+5)*kA(n) - (2*n+3)*kA(n+1))*(-(n+2) / (2*n+3));
% interior points
for n = 2:N-1
    Amat(n,n-1) = ((2*n-1)*kA(n-1) - (2*n-3)*kA(n))*((n-1) / (2*n-1));
    Amat(n,n) = ((2*n-1)*kA(n-1) - (2*n-3)*kA(n))*(-n / (2*n+1)) ...
        - ((2*n+5)*kA(n) - (2*n+3)*kA(n+1))*((n+1) / (2*n+1));
    Amat(n,n+1) = -((2*n+5)*kA(n) - (2*n+3)*kA(n+1))*(-(n+2) / (2*n+3));
end
% right boundary
n = N;
kN = (N+1+0.5).*coth((N+1+0.5)*a) - coth(a);
Amat(n,n-1) = ((2*n-1)*kA(n-1) - (2*n-3)*kA(n))*((n-1) / (2*n-1));
Amat(n,n) = ((2*n-1)*kA(n-1) - (2*n-3)*kA(n))*(-n / (2*n+1)) ...
    - ((2*n+5)*kA(n) - (2*n+3)*kN)*((n+1) / (2*n+1));

% Right hande side
rhs = -(sqrt(2)*exp(-(nA+0.5)*a) ./ ((2*nA+1).*sinh(a))) ...
    .* ((2*nA+1).^2 .* (exp(a)./(2*nA-1) + exp(-a)./(2*nA+3)) .* csch((nA+0.5)*a) ...
    - (2*nA-1).*csch((nA-0.5)*a) - (2*nA+3).*csch((nA+1.5)*a));

% Solve linear system
An = Amat \ rhs;

%% Bn Coefficients (equation 18)
nB = [1:N-1]';
kB = (nB+0.5).*coth((nB+0.5)*a) - coth(a);
Bn = -(2*nB+1).*An(1:end-1) + (nB+2).*An(2:end);
Bn(2:end) = Bn(2:end) + (nB(2:end)-1).*An(1:end-2);

%% Cn Coefficients (equation 19)
nC = [1:N-1]';
kC = (nC+0.5).*coth((nC+0.5)*a) - coth(a);
LC = -(1/sqrt(2)) * (exp(-(nC-0.5)*a)./(2*nC-1) - exp(-(nC+1.5)*a)./(2*nC+3));
Cn = 4*LC*csch(a).*csch((nC+0.5)*a) - 2*kC.*(-An(1:end-1) + ((nC+2)./(2*nC+3)).*An(2:end) );
Cn(2:end) = Cn(2:end) - 2*kC(2:end).*((nC(2:end)-1)./(2*nC(2:end)-1)).*An(1:end-2);

%% Dn Coefficients (equation 21)
nD = [0:N-1]';
kD = (nD+0.5).*coth((nD+0.5)*a) - coth(a);
Dn = 0.5*(nD+1).*(nD+2).*An;
Dn(3:end) = Dn(3:end) - 0.5*(nD(3:end)-1).*nD(3:end).*An(1:end-2);

%% En Coefficients (equation 22)
nE = [0:N-1]';
kE = (nE+0.5).*coth((nE+0.5)*a) - coth(a);
LE = -(1/sqrt(2)) * (exp(-(nE-0.5)*a)./(2*nE-1) - exp(-(nE+1.5)*a)./(2*nE+3));
En = (sqrt(2)*(2*nE+1).*exp(-(nE+0.5)*a) - LE*csch(a)).*csch((nE+0.5)*a) ...
    - kE.*((nE+1).*(nE+2)./ (2*nE+3)).*An;
En(3:end) = En(3:end) + kE(3:end).*((nE(3:end)-1).*nE(3:end)./(2*nE(3:end)-1)).*An(1:end-2);

YB = -sqrt(2)/6*sinh(a)^2 * (sum(En) + sum(nC.*(nC+1).*Cn));

YC = (4/3)*(1/3)*(1 - sinh(a)^3/(4*sqrt(2)) * (sum((2+exp(-(2*nA+1)*a)).*nA.*(nA+1)*2.*An) ...
    + sum((2+exp(-(2*nC+1)*a)).*nC.*(nC+1).*Cn*coth(a)) ...
    - sum((2+exp(-(2*nE+1)*a)).*(2*nE+1-coth(a)).*En) ...
    + sum((2-exp(-(2*nB+1)*a)).*nB.*(nB+1).*Bn*coth(a)) ...
    - sum((2-exp(-(2*nD+1)*a)).*(2*nD+1-coth(a)).*Dn)));
