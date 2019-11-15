function draw(r, q, param)
% DRAW  Draws snapshot of the system configuration
N = size(r,2);
totalCharge = sum(q);

for in = 1:N % color based on total charge for representation
    if q(in) > 0
        style = [255 102 102]./256;
    elseif q(in) == 0
        style = [0 0 0];     
    else 
        style = [102 102 255]./256;
    end
    
    [x,y,z] = sphere(30);
    s1= surf(x+r(1,in),y+r(2,in),z+r(3,in));
    hold on
    
    %camlight
    %lighting phong
    set(s1,'facecolor',style);
    s1.EdgeColor = 'none';
end

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
axis equal;
axis ([ 0 param.domainLength(1) 0 param.domainLength(2) 0 param.domainLength(3)]);

%axis off;

hold on;
threshold = 0; 
x1 = [ 0 0 param.domainLength(1) param.domainLength(1)];
y1 = [0  param.domainLength(2) param.domainLength(2) 0]; 
z1 = ones(1,numel(x1))* threshold;
v = patch(x1,y1,z1, 'r');
set(v,'facealpha',0.5);
set(v,'edgealpha',0.5);
set(gcf,'renderer','opengl') ;

hold on;
threshold = 0; 
x1 = [ 0 0 param.domainLength(1) param.domainLength(1)];
y1 = [0  param.domainLength(2) param.domainLength(2) 0]; 
z1 = ones(1,numel(x1))* param.domainLength(3);
v = patch(x1,y1,z1, 'b');
set(v,'facealpha',0.5);
set(v,'edgealpha',0.5);
set(gcf,'renderer','opengl') ;

view(180,0)
hold off

