function [YA,YB] = TestYBYC()
%{
Test resistance functions YA and YB against lubrication result from
Bossis, Meunier, Sherwood, Phys. Fluids A, 1991, 3, 1853.
%}

%% Reproduce results from A.J. Goldman, R.G. Cox, and H. Brenner, Chem. Eng. Sci. 22, 637 (1967).
a = [3,2,1.5,1,0.5,0.3,0.1,0.08]';
h = cosh(a);
xi = h-1;
YA = zeros(size(xi));
YB = zeros(size(xi));
fprintf('\nTabel from Brenner, Chem Eng Sci, 1967, 22, 637\n');
for i = 1:length(xi)
    [YB(i),YC(i)] = CalcYBYC(xi(i));
    fprintf('%.2f\t%.6f\t%.4e\t%.4f\n',a(i),h(i),YB(i),0.75*YC(i));
end
fprintf('\n');

%% Plot with lubrication results
Min = -3;
Max = 3;
Npoints = 100;

xi = logspace(Min,Max,Npoints)';
YB = zeros(size(xi));
YC = zeros(size(xi));
for i = 1:length(xi)
    [YB(i),YC(i)] = CalcYBYC(xi(i));
end

% Near Field
xinf = logspace(Min,0,Npoints)';
YBnf = -(2/15)*log(xinf) - 0.2526; % lubrication
YCnf = -(8/15)*log(xinf) + 0.5089; % lubrication

semilogx(xi,YB,xinf,YBnf,'k:');
set(gca,'FontName','Times New Roman');
xlabel('$\xi$','FontSize',16,'Interpreter','latex');
ylabel('$Y^{B}$','FontSize',16,'Interpreter','latex');

figure
semilogx(xi,YC,xinf,YCnf,'k:');
set(gca,'FontName','Times New Roman');
xlabel('$\xi$','FontSize',16,'Interpreter','latex');
ylabel('$Y^{C}$','FontSize',16,'Interpreter','latex');

