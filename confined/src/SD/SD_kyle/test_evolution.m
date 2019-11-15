addpath('Sphere2PlaneFF', 'SpherePlaneFF', 'SpherePlaneNF', 'SphereSphereNF')
LoadTablesM2PlaneFF()
LoadTablesMPlaneFF()
LoadTablesRPlaneNF()

[t r q] = evolution(10);

for i = 1:length(t)
    draw(r(i,:), q(i))
    pause(0.1)
end