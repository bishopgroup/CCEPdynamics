function [XA] = TestXC()
%{
Test resistance function XC against lubrication result.
%}

%% Reproduce results from G.B. Jeffrey, Proc. London Math. Soc. 14, 327 (1915). 
a = [0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0]';
h = cosh(a);
xi = h-1;
XC = zeros(size(xi));
fprintf('\nTable from Jeffrey\n');
for i = 1:length(xi)
    [XC(i)] = CalcXC(xi(i));
    fprintf('%.1f\t%.4f\t%.4f\n',a(i),1/h(i),0.75*XC(i));
end
fprintf('\n');

%% Plot with lubrication
Min = -3;
Max = 3;
Npoints = 100;

%% Calc XC
xi = logspace(Min,Max,Npoints)';
XC = zeros(size(xi));
for i = 1:length(xi)
    XC(i) = CalcXC(xi(i));
end

%% Near Field
xinf = logspace(Min,-1,Npoints);
XCnf = (4/3)*(zeta(3) + 0.5*xinf.*log(xinf)); % lubrication
semilogx(xi,XC,xinf,XCnf,'k:');
set(gca,'FontName','Times New Roman');
xlabel('$\xi$','FontSize',16,'Interpreter','latex');
ylabel('$X^{C}$','FontSize',16,'Interpreter','latex');
