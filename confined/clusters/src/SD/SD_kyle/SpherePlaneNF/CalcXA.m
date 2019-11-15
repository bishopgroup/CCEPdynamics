function [XA] = CalcXA(xi)
%{
Solves for the resistance functions XA a sphere moving normal
to a plane surface (Brenner, Chem. Eng. Sci. 1961, 16, 242).

The input paramter is xi = (h-a)/a is the dimensionless separation between
the surface of the sphere and the plane.
%}

a = acosh(1+xi);  % bispherical parameter
error = 1e-14;
N = ceil(1-log(error)/(2*a)); % Number of terms in the sum

n = 1:N;
terms = (n.*(n+1)) ./ ((2*n-1).*(2*n+3)) .* ( (2*sinh((2*n+1)*a) ...
    + (2*n+1)*sinh(2*a)) ./ (4*sinh((n+0.5)*a).^2 - (2*n+1).^2*sinh(a)^2) - 1);
XA = (4/3)*sinh(a) * sum(terms); % Equation 2.19
    