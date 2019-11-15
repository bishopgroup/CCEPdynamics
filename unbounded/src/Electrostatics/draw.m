function draw(r, q, Parameter)
% DRAW  Draws snapshot of the system configuration
N = size(r,2);

nsample = 100;
phi = linspace(0, 2 * pi, nsample);
[xs ys] = pol2cart(phi, 1);

for in = 1:N
    if q(in) > 0
        style = [255 102 102]./256;
    elseif q(in) == 0
        style = '-k';     
    else 
        style = [102 102 255]./256;
    end
    
    part(r(:,in),style,xs,ys);
    hold on;
    
end

axis equal;
axis ([ 0 Parameter.domainLength(1) 0 Parameter.domainLength(3) ]);
set(gca,'xtick',[])
set(gca,'ytick',[])
hold off;

%--------------------------------------------------------------------------
function part(r, style,xs,ys)
% PARTICLE  Draws representation of the particle
%   Fucntion draws a cricle of radius 'rad' centerd in 'r' using
%   supplied user supplied style 'style'.

plot(xs + r(1), ys + r(3), 'Color', style);
