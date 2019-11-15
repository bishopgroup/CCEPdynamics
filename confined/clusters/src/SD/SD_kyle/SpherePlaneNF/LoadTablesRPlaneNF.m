function LoadTablesRPlaneNF(varargin)
%{
Load tabulated functions for sphere plane resistance.
%}


global SPH_XA_ln SPH_YA SPH_YB SPH_XC SPH_YC eps3jk d3d3 I

   
I = eye(3);
eps3jk = [0 1 0; -1 0 0; 0 0 0];
d3d3 = [0 0 0; 0 0 0; 0 0 1];

%% Sphere-Plane Resistance
data = load('SpherePlaneResistance.txt');
SPH_XA_ln = griddedInterpolant(log(data(:,1)),log(data(:,2)));
SPH_YA = griddedInterpolant(log(data(:,1)),data(:,3));
SPH_YB = griddedInterpolant(log(data(:,1)),data(:,4));
SPH_XC = griddedInterpolant(log(data(:,1)),data(:,5));
SPH_YC = griddedInterpolant(log(data(:,1)),data(:,6));

