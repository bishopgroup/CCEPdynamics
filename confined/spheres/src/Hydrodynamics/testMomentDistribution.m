% Charge assignment
clear all



% Force
fxDensity = zeros(Mesh.nMesh, Mesh.nMesh, Mesh.nMesh);
fxDensity = zeros(Mesh.nMesh, Mesh.nMesh, Mesh.nMesh);

% Find indices of nearest grid points
gridIndex = zeros(Particle.nParticle, 3);
gridIndex(:,1) = mod(round( Particle.position(:,1) / Mesh.spacing ),...
    Mesh.nMesh) + 1;
gridIndex(:,2) = mod(round( Particle.position(:,2) / Mesh.spacing ),...
    Mesh.nMesh) + 1;
gridIndex(:,3) = mod(round( Particle.position(:,3) / Mesh.spacing ),...
    Mesh.nMesh) + 1;

% Find positions of nearest grid points
gridPosition = zeros(Particle.nParticle, 3);
gridPosition(:,1) = Mesh.x1( gridIndex(:,1) );
gridPosition(:,2) = Mesh.x2( gridIndex(:,2) );
gridPosition(:,3) = Mesh.x3( gridIndex(:,3) );

% Fractional displacement between particles and proximal grid point
fractionalDisplacement = Particle.position - gridPosition;
fractionalDisplacement = fractionalDisplacement ... 
    - round(fractionalDisplacement / Parameter.domainLength) ...
    * Parameter.domainLength;
fractionalDisplacement = fractionalDisplacement / Mesh.spacing;

% preallocate
wp = zeros(3,3); % NGRID x NDIM
wq = zeros(3,3); % NGRID x NDIM

for iParticle = 1:Particle.nParticle % loop over particles
    % Weighting Function (charge)
    wq(1,:) = 0.5*(fractionalDisplacement(iParticle,:).^2 ...
        - fractionalDisplacement(iParticle,:));
    wq(2,:) = (1 - fractionalDisplacement(iParticle,:).^2);
    wq(3,:) = 0.5*(fractionalDisplacement(iParticle,:).^2 ...
        + fractionalDisplacement(iParticle,:));
    
    % Weighting Function (dipole)
    wp(1,:) = 0.5*(-1 + 2*fractionalDisplacement(iParticle,:));
    wp(2,:) = -2*fractionalDisplacement(iParticle,:);
    wp(3,:) = 0.5*(1 + 2*fractionalDisplacement(iParticle,:));  
    wp = wp / Mesh.spacing;
    
    % Distribute Charge
    n = [-1:1];
    i1 = mod( n + gridIndex(iParticle,1) - 1, Mesh.nMesh) + 1;
    i2 = mod( n + gridIndex(iParticle,2) - 1, Mesh.nMesh) + 1;
    i3 = mod( n + gridIndex(iParticle,3) - 1, Mesh.nMesh) + 1;
 
    for i = 1:3
        for j = 1:3
            for k = 1:3
                chargeDensity(i1(i),i2(j),i3(k)) = chargeDensity(i1(i),i2(j),i3(k)) + ...
                    Particle.farFieldCharge(iParticle) * (wq(i,1)*wq(j,2)*wq(k,3)) + ...
                    Particle.farFieldDipole(1,iParticle) * (wp(i,1)*wq(j,2)*wq(k,3))  + ...
                    Particle.farFieldDipole(2,iParticle) * (wq(i,1)*wp(j,2)*wq(k,3))  + ...
                    Particle.farFieldDipole(3,iParticle) * (wq(i,1)*wq(j,2)*wp(k,3)) ;
            end
        end
    end
end

return


%% 2D scalar example (unit grid spacing)
Q = 1;            % zeroth moment
Qi = [0.1, 0.5];  % first moment
Qij = [0.2, 0.6;
       0.6, 0.4]; % second moment
nGrid = 4; % number of grid points
x = 0.1;
y = 0.2;

% weighting
w0x = [-(-2 + x)*(-1 + x)*x/6, ...
    0.5*(-2 + x)*(-1 + x)*(1 + x), ...
    0.5*x*(2 + x - x^2), ...
    x*(-1 + x^2)/6];
w0y = [-(-2 + y)*(-1 + y)*y/6, ...
    0.5*(-2 + y)*(-1 + y)*(1 + y), ...
    0.5*y*(2 + y - y^2), ...
    y*(-1 + y^2)/6];
w1x = [-(1/3) + x - 0.5*x^2,...
    0.5*(-1 + x*(-4 + 3*x)),...
    1 + x - 1.5*x^2, ...
    (-1 + 3*x^2) / 6];
w1y = [-(1/3) + y - 0.5*y^2,...
    0.5*(-1 + y*(-4 + 3*y)),...
    1 + y - 1.5*y^2, ...
    (-1 + 3*y^2) / 6];
w2x = [(1 - x)/2, ...
    -1 + (3*x)/2, ...
    0.5*(1 - 3*x), ...
    x/2];
w2y = [(1 - y)/2, ...
    -1 + (3*y)/2, ...
    0.5*(1 - 3*y), ...
    y/2];

% charge assignment
rho = zeros(4,4);


for i = 1:4
    for j = 1:4
        rho(i,j) = rho(i,j) + w0x(i)*w0y(j)*Q;
        
        rho(i,j) = rho(i,j) + w1x(i)*w0y(j)*Qi(1);
        rho(i,j) = rho(i,j) + w0x(i)*w1y(j)*Qi(2);

        rho(i,j) = rho(i,j) + w2x(i)*w0y(j)*Qij(1,1);
        rho(i,j) = rho(i,j) + w1x(i)*w1y(j)*Qij(2,1);
        rho(i,j) = rho(i,j) + w1x(i)*w1y(j)*Qij(1,2);
        rho(i,j) = rho(i,j) + w0x(i)*w2y(j)*Qij(2,2);
    end
end

xg = [-1:2];
yg = [-1:2];
[Yg,Xg] = meshgrid(yg,xg);

sum(rho(:))
sum(sum((Xg-x).*rho))
sum(sum((Yg-y).*rho))
sum(sum(0.5*(Xg-x).*(Yg-y).*rho))
sum(sum(0.5*(Yg-y).*(Xg-x).*rho))
sum(sum((Xg-x).^2.*rho))
sum(sum((Yg-y).^2.*rho))


