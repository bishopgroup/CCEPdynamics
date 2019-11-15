%Plots drift velocity as a function of channel height.
%
% Validates SD code by reproducing dependence of drift velocity on channel
% height (see Fig. 12 from the referenced paper).
%
% References
% JW Swan and JF Brady, Phys Fluid 22, 103301 (2010).

addpath('Sphere2PlaneFF', 'SpherePlaneFF', 'SpherePlaneNF', 'SphereSphereNF')

LoadTablesM2PlaneFF()
LoadTablesMPlaneFF()
LoadTablesRPlaneNF()

% Channel heights investigated in the paper.
H = [3 4 5 7 10 15 20 30 50 75 100];

% Step to use during calculating numerical derivatives for a given $H$.
% Their number must be equal number of investigated channel heights.
dh = [0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02 0.02, 0.05, 0.05, 0.05];

% Estimate maximal number of points.
nmax = H(end) / dh(1);

X = zeros(nmax, length(H));
Y = zeros(nmax, length(H));
for i = 1:length(H)
    [h, v] = vdrift(H(i), dh(i));
    
    % Pad solution vectors with zeros.
    l = nmax - length(h);
    X(:,i) = wextend('1', 'zpd', h - 1, l, 'r');
    Y(:,i) = wextend('1', 'zpd', v, l, 'r');
end

loglog(X, Y)
ylim([1.0e-4, 1]);
xlabel('$h - a$',...
       'FontSize', 16, 'Interpreter', 'latex');
ylabel('$\nabla \cdot \mathbf{R}_{FU}^{-1} \cdot \delta_{3}$',...
       'FontSize', 16, 'Interpreter', 'latex');
leg = legend('$H/a = 3$', '$H/a = 4$', '$H/a = 5$', '$H/a = 7$', ...
             '$H/a = 10$', '$H/a = 15$', '$H/a = 25$', '$H/a = 30$', ...
             '$H/a = 50$', '$H/a = 75$', '$H/a = 100$');
set(leg, 'Location', 'SouthWest', 'Interpreter', 'latex');