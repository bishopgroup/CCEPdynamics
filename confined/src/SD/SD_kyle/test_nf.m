% Check consistency of diagonal elements of $R_{FU}$ matrix with their
% analytical approximations, for different channel heights.
%
% References
% G Bossis, A Meunier, JD Sherwood, Phys. Fluids A 3(8), 1853 (1991)

H = [3 4 5 7 10 15 20 30 50 75 100];

% The diagonal cooeficient of the resistance matrix for a single particle
% near a plane can be expressed can be rigorously expressed as a function
% of "slit's" width between particle and the plane, $\varepsilon a = h -
% a$, where $a$ is particle's radius.
eps = linspace(1.0e-03, 1.0e-01, 100);
fx = -8/15 * log(eps) + 0.9588;
fz = 1 ./ eps - 1/5 * log(eps) + 0.97128;

% Calculate diagonal elements of $R_{FU}$ matrix using SD approach for
% different channel heights.
m = length(H);
n = length(eps);

R11 = zeros(n, m);
R22 = zeros(n, m);
R33 = zeros(n, m);
for j = 1:m;
    for i = 1:n
        h = 1 + eps(i);
        xi = h / H(j);

        R = RSD(xi, H(j));
        RFU = R(1:3,1:3);

        R11(i,j) = RFU(1,1);
        R22(i,j) = RFU(2,2);
        R33(i,j) = RFU(3,3);
    end
end

% Compare results obtained from SD approach with exact values.
plot(eps, R11, eps, fx, 'k-')
xlabel('$\varepsilon$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$R_{11}$', 'FontSize', 16, 'Interpreter', 'latex');
pause;

plot(eps, R22, eps, fx, 'k-')
xlabel('$\varepsilon$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$R_{22}$', 'FontSize', 16, 'Interpreter', 'latex');
pause;

plot(eps, R33, eps, fz, 'k-')
xlabel('$\varepsilon$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$R_{33}$', 'FontSize', 16, 'Interpreter', 'latex');
pause;

% Save the data to external files.
tab = [eps' R11; eps' R22; eps' R33];
save('Rii-H.dat', 'tab', '-ascii');

tab = [eps' fx'; eps' fx'; eps' fz'];
save('Rii-exact.dat', 'tab', '-ascii');