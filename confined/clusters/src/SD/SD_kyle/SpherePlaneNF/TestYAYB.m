function [YA,YB] = TestYAYB()
%{
Test resistance functions YA and YB against lubrication result from
G. Bossis, A. Meunier, and J.D. Sherwood, Phys. Fluids A Fluid Dyn. 3, 1853 (1991).
%}

%% Reproduce results from M.E. O’Neill, Mathematika 11, 67 (1964)
a = [3,2,1,0.5,0.3,0.1,0.08,0.06,0.04,0.03,0.02]';
h = cosh(a);
xi = h-1;
YA = zeros(size(xi));
YB = zeros(size(xi));
fprintf('\nTabel from ONeill, Mathematika, 1964, 11, 67\n');
for i = 1:length(xi)
    [YA(i),YB(i)] = CalcYAYB(xi(i));
    fprintf('%.2f\t%.4f\t%.4f\t%.5f\n',a(i),h(i),YA(i),0.75*YB(i));
end
fprintf('\n');

%% Reproduce results from A.J. Goldman, R.G. Cox, and H. Brenner, Chem. Eng. Sci. 22, 637 (1967)
a = [3,2,1.5,1,0.5,0.3,0.1,0.08]';
h = cosh(a);
xi = h-1;
YA = zeros(size(xi));
YB = zeros(size(xi));
fprintf('\nTabel from Brenner, Chem Eng Sci, 1967, 22, 637\n');
for i = 1:length(xi)
    [YA(i),YB(i)] = CalcYAYB(xi(i));
    fprintf('%.2f\t%.6f\t%.4f\t%.4e\n',a(i),h(i),YA(i),0.75*YB(i));
end
fprintf('\n');

%% Plot with lubrication results
Min = -3;
Max = 3;
Npoints = 100;

xi = logspace(Min,Max,Npoints)';
XA = zeros(size(xi));
for i = 1:length(xi)
    [YA(i),YB(i)] = CalcYAYB(xi(i));
end

% Near Field
xinf = logspace(Min,0,Npoints)';
YAnf = -(8/15)*log(xinf) + 0.9588; % lubrication
YBnf = -(2/15)*log(xinf) - 0.2526; % lubrication

semilogx(xi,YA,xinf,YAnf,'k:');
set(gca,'FontName','Times New Roman');
xlabel('$\xi$','FontSize',16,'Interpreter','latex');
ylabel('$Y^{A}$','FontSize',16,'Interpreter','latex');

figure
semilogx(xi,YB,xinf,YBnf,'k:');
set(gca,'FontName','Times New Roman');
xlabel('$\xi$','FontSize',16,'Interpreter','latex');
ylabel('$Y^{B}$','FontSize',16,'Interpreter','latex');

