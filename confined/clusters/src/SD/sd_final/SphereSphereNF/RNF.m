function R = RNF(h)
%{
Computes the 2-body resistance matrix for a single sphere of unit radius
positioned at a distance h above a solid planar wall.
%}

global SPH_XA_ln SPH_YA SPH_YB SPH_XC SPH_YC eps3jk d3d3 I

xi = h-1;
lnxi = log(xi); % surface separation

if xi < 1e-3
   XA = 1/xi - 1/5*lnxi + 0.97128;
   YA = -8/15*lnxi + 0.9588;
   XC = 4/3*(zeta(3)-0.5*xi*lnxi);
   YC = -8/15*lnxi + 0.5089;
   YB = -2/15*lnxi - 0.2526;
else
   XA = exp(SPH_XA_ln(lnxi));
   YA = SPH_YA(lnxi);
   XC = SPH_XC(lnxi);
   YC = SPH_YC(lnxi);
   YB = SPH_YB(lnxi);
end

%% Compute Pairwise Resistance Matrix
Rfu = XA*d3d3 + YA*(I - d3d3);
Rlu = XC*d3d3 + YC*(I - d3d3);
Rlo = -YB*eps3jk;

R = [Rfu,  Rlo';
     Rlo,  Rlu];
