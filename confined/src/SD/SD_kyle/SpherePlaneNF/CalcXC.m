function [XA] = CalcXC(xi)
%{
Solves for the resistance functions XC a sphere rotating normal
to a plane surface.

G.B. Jeffrey, Proc. London Math. Soc. 14, 327 (1915)

The input paramter is xi = (h-a)/a is the dimensionless separation between
the surface of the sphere and the plane.
%}

a = acosh(1+xi);  % bispherical parameter
error = 1e-14;
N = ceil(-2*log(error)/(3*a)); % Number of terms in the sum

n = 0:N;
XA = (4/3)*sinh(a)^3*sum( csch((n+1)*a).^3 );


    