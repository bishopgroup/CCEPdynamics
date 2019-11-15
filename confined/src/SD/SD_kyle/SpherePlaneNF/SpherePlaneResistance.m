%% Domain
Min = -3;
Max = 3;
Npoints = 1000;
xi = logspace(Min,Max,Npoints)';

% Generate Coefficients
XA = zeros(Npoints,1);
YA = XA; YB = XA;
YB2 = XA; XC = XA;
YC = XA;
for i = 1:length(xi)
    [XA(i)] = CalcXA(xi(i));
    [YA(i),YB(i)] = CalcYAYB(xi(i));
    [YB2(i),YC(i)] = CalcYBYC(xi(i));
    [XC(i)] = CalcXC(xi(i));
end

% Combine Data
table = [xi, XA, YA, YB, XC YC];
save('SpherePlaneResistance.txt','table','-ascii','-tabs','-double');
