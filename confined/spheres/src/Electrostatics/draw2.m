function draw2(r, q, param)
% DRAW  Draws snapshot of the system configuration
N = size(r,2);
totalCharge = sum(q);

for in = 1:N % color based on total charge for representation
    if totalCharge > 0
        style = [1 0 0];
    elseif q(in) == 0
        style = [0 0 0];     
    else 
        style = [0 0 1];
    end
    
    [x,y,z] = sphere(30);
    s1= surf(x+r(1,in),y+r(2,in),z+r(3,in));
    hold on
    camlight
    lighting phong
    set(s1,'facecolor',style);
    s1.EdgeColor = 'none';
end

axis equal;
axis ([ 0 param.L1 0 param.L2 0 param.L3]);
%axis off;

hold on;
threshold = 0; 
x1 = [ 0 0 param.L1 param.L1];
y1 = [ 0 param.L2 param.L2 0]; 
z1 = ones(1,numel(x1))* threshold;
v = patch(x1,y1,z1, 'r');
set(v,'facealpha',0.8);
set(v,'edgealpha',0.8);
set(gcf,'renderer','opengl') ;

hold on;
threshold = 0; 
x1 = [ 0 0 param.L1 param.L1];
y1 = [ 0 param.L2 param.L2 0]; 
z1 = ones(1,numel(x1))* param.L3;
v = patch(x1,y1,z1, 'b');
set(v,'facealpha',0.8);
set(v,'edgealpha',0.8);
set(gcf,'renderer','opengl') ;

%hold off

