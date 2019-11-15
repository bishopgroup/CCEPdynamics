function NearField = initializeNearField(Parameter)
%--------------------------------------------------------------------------
% initializeNearField - loads two-sphere capacitance functions.
%
%  NearField = initializeNearField(Parameter)
%--------------------------------------------------------------------------

%% Two Sphere Capacitance 
data = importdata('TwoSphereCapacitance.txt');
NearField.XA11 = griddedInterpolant(log(data(:,1)-2),data(:,2));
NearField.XA12 = griddedInterpolant(log(data(:,1)-2),data(:,3));
NearField.XB11 = griddedInterpolant(log(data(:,1)-2),data(:,4));
NearField.XB12 = griddedInterpolant(log(data(:,1)-2),data(:,5));
NearField.XC11 = griddedInterpolant(log(data(:,1)-2),data(:,6));
NearField.XC12 = griddedInterpolant(log(data(:,1)-2),data(:,7));
NearField.YC11 = griddedInterpolant(log(data(:,1)-2),data(:,8));
NearField.YC12 = griddedInterpolant(log(data(:,1)-2),data(:,9));
