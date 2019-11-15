function draw(r, q)
% DRAW  Draws snapshot of the system configuration
%   Detailed explanation goes here

global H L

if q > 0
    style = '-r';
elseif q == 0
    style = '-k';
else 
    style = '-b';
end
    
part(r([1 3]), 1, style);
hold on;
    
label = sprintf('%.2f', q);
text(r(1), r(3), label, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
hold on;

axis equal;
axis ([ 0 L 0 H ]);
hold off;
end

function part(r, radius, style)
% PARTICLE  Draws representation of the particle
%   Fucntion draws a cricle of radius 'rad' centerd in 'r' using
%   supplied user supplied style 'style'.

nsample = 100;
phi = linspace(0, 2 * pi, nsample);
[x y] = pol2cart(phi, radius);
plot(x(:) + r(1), y(:) + r(2), style);
end