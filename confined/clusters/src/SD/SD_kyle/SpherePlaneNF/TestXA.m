function [XA] = TestXA()
%{
Test resistance function XA against lubrication result from
Bossis, Meunier, Sherwood, Phys. Fluids A, 1991, 3, 1853.
%}

Min = -3;
Max = 3;
Npoints = 100;

%% Calc XA
xi = logspace(Min,Max,Npoints)';
XA = zeros(size(xi));
for i = 1:length(xi)
    XA(i) = CalcXA(xi(i));
end

%% Near Field
XAnf = 1./xi - (1/5)*log(xi) + 0.97128; % lubrication
loglog(xi,XA,xi,XAnf,'k:');
set(gca,'FontName','Times New Roman');
xlabel('$\xi$','FontSize',16,'Interpreter','latex');
ylabel('$X^{A}$','FontSize',16,'Interpreter','latex');
