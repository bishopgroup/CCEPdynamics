%--------------------------------------------------------------------------
% linkedSphereMain - Simulates the trajectory of linked rigid spheres 
%
%--------------------------------------------------------------------------

clear all
fprintf('\n Running simulation for given linked sphere parameters...\n');


%% Set Path
rootPath = pwd();
addpath( fullfile( rootPath, 'src', 'ConfinedElectrodynamics'));
addpath( fullfile( rootPath, 'src', 'Hydrodynamics'));

%% Input Parameters
Parameter.L3 = 20; % distance between electrodes (scaled by a)
Parameter.splittingAlpha = (1/4)*Parameter.L3^2; 

Parameter.L1 = 2*Parameter.L3;
Parameter.L2 = 2*Parameter.L3;
Parameter.delta = 0.1; % minimum surface separation

%% Initialize Clusters
%{
Particle.nParticle = 6;
sphereSeparation = 2 + Parameter.delta/2;
Particle.bodyFixedPosition = [-sphereSeparation/2, 0, 0; 
                               sphereSeparation/2, 0, 0;
                               0, 0, -sphereSeparation*cos(pi/6);
                               0, 0, -sphereSeparation-sphereSeparation*cos(pi/6);
                               0, 0, -2*sphereSeparation-sphereSeparation*cos(pi/6)
                               0, 0, -3*sphereSeparation-sphereSeparation*cos(pi/6)]'; % particle position, 3 x nParticle        
Particle.clusterTranslation = [7,7,11]'; % cluster position
axisAngleIC = [0,1,0,0]; % axis-angle representation
Particle.adjmat = sparse([0 1 1 0 0 0;
                          1 0 1 0 0 0;
                          1 1 0 1 0 0;
                          0 0 1 0 1 0;
                          0 0 0 1 0 1;
                          0 0 0 0 1 0]);
Particle.adjmatCond = sparse([0 1;
                              1 0]);% determines connectivity world cluster

Particle.ncParticle = 2;% lower 3 conductive
Particle.conductivityIndex = [1 1 0 0 0 0];
%}

%{
Particle.nParticle = 4;
sphereSeparation = 2 + Parameter.delta/2;
Particle.bodyFixedPosition = [-sphereSeparation/2, 0, 0; 
                               sphereSeparation/2, 0, 0;
                               0, 0, -sphereSeparation*cos(pi/6);
                               0, 0, -sphereSeparation-sphereSeparation*cos(pi/6)]'; % particle position, 3 x nParticle        
Particle.clusterTranslation = [7,7,11]'; % cluster position
axisAngleIC = [0,1,0,0]; % axis-angle representation
Particle.adjmat = sparse([0 1 1 0;
                          1 0 1 0;
                          1 1 0 1;
                          0 0 1 0]);
Particle.adjmatCond = sparse([0 1;
                              1 0]);% determines connectivity world cluster

Particle.ncParticle = 2;% upper 2 conductive
Particle.conductivityIndex = [1 1 0 0];
%}

%
Particle.nParticle = 4;
sphereSeparation = 2 + Parameter.delta/2;
Particle.bodyFixedPosition = [-sphereSeparation/2, 0, 0; 
                               sphereSeparation/2, 0, 0;
                               0, 0, -sphereSeparation*cos(pi/6);
                               0, 0, -sphereSeparation-sphereSeparation*cos(pi/6)]'; % particle position, 3 x nParticle        
Particle.clusterTranslation = [20,20,10]'; % cluster position
axisAngleIC = [0,1,0,0]; % axis-angle representation
Particle.adjmat = sparse([0 1 1 0;
                          1 0 1 0;
                          1 1 0 1;
                          0 0 1 0]);
Particle.adjmatCond = sparse([0 1;
                              1 0]);% determines connectivity world cluster

Particle.ncParticle = 2;% upper 2 conductive
Particle.conductivityIndex = [1 1 0 0];
%}


%% Derived parameters
Particle.orientation = axisAngleToQuaternion(axisAngleIC);% converted into quaternion for calculations
Particle = CalcPos(Particle); % computes positions
Particle.charge = 0.3*ones(Particle.ncParticle,1);                              % nParticle x 1
Particle.field = [zeros(Particle.ncParticle,2),ones(Particle.ncParticle,1)]';  % 3 x nParticle
Particle.potential =  Particle.charge;    % initial guess
Particle.farFieldCharge = Particle.charge;   % initial guess 
Particle.farFieldDipole = Particle.field;  % initial guess 
Particle.force = zeros(3,1); % 3 x 1 only works wth single cluster of particle
Particle.F = bsxfun(@times,Particle.field,Particle.charge'); % 3 x Np for every particle
Particle.velocity = Particle.force;
ParticleOld = Particle;

%% Initialize Matrices Mesh and Tables
Parameter.domainLength = [Parameter.L1 Parameter.L2 Parameter.L3];
Parameter.alpha = Parameter.splittingAlpha;
[Mesh] = initializeMesh(Parameter);
BuildTables();

%% Integrate
% integration parameters
deltaT = 0.1; % time step
maximumTime = 200; % max time
numberSteps = round(maximumTime/deltaT);  % number of steps

% Save Initial Conditions
timeIteration = 1;
positionPool = reshape(Particle.position,1,3*Particle.nParticle);
quaternionPool = reshape(Particle.orientation,1,4);
chargePool = Particle.charge';
timePool = 0;
eulerRestart = true;
tic

while timeIteration < numberSteps
    % Prepare for Integration Step
    Particle = calculateForce(Parameter,Mesh,Particle); % Compute the forces
    Particle = calculateVelocity(Parameter,Particle); % Compute the velocities
    %Particle = calculateVelocityHydro(Parameter,Particle); % Compute the velocities
    ParticleNext = Particle;
        
    % Euler Step
    if eulerRestart % Euler Step
        ParticleNext.clusterTranslation = Particle.clusterTranslation + deltaT*Particle.velocity;
        ParticleNext.orientation = Particle.orientation + deltaT*angularVelocityToQuaternionRate(Particle.orientation,Particle.angularVelocity)';
    else % 2nd Order Adams-Bashforth Step 
        ParticleNext.clusterTranslation = Particle.clusterTranslation + deltaT*(1.5*Particle.velocity - 0.5*ParticleOld.velocity);
        ParticleNext.orientation = Particle.orientation + deltaT*(1.5*angularVelocityToQuaternionRate(Particle.orientation,Particle.angularVelocity)' - 0.5*angularVelocityToQuaternionRate(ParticleOld.orientation,ParticleOld.angularVelocity)');
    end
    ParticleNext = CalcPos(ParticleNext);
    
    % Check for Overlaps
    [fractionalMove,flag,adjmat] = Overlap(Parameter,ParticleNext,Particle);
    if flag % there are overlaps
        % Take a smaller step
        if eulerRestart % Euler Step
            ParticleNext.clusterTranslation = Particle.clusterTranslation +...
                fractionalMove*deltaT*Particle.velocity;
            ParticleNext.orientation = Particle.orientation +...
                fractionalMove*deltaT*angularVelocityToQuaternionRate(Particle.orientation,Particle.angularVelocity)';
        else % 2nd Order Adams-Bashforth Step
            ParticleNext.clusterTranslation = Particle.clusterTranslation +...
                fractionalMove*deltaT*(1.5*Particle.velocity - 0.5*ParticleOld.velocity);
            ParticleNext.orientation = Particle.orientation +...
                fractionalMove*deltaT*(1.5*angularVelocityToQuaternionRate(Particle.orientation,Particle.angularVelocity)'...
                - 0.5*angularVelocityToQuaternionRate(ParticleOld.orientation,ParticleOld.angularVelocity)');
        end
        ParticleNext = CalcPos(ParticleNext);
        
        % Equilibrate Charges
        ParticleNext = Equilibrate(Parameter,Mesh,ParticleNext,adjmat);

        % Save Intermediate Results
        timeIteration = timeIteration + 1;
        positionPool = [positionPool;reshape(ParticleNext.position,1,3*Particle.nParticle)];
        quaternionPool = [quaternionPool;reshape(ParticleNext.orientation,1,4)];
        chargePool = [chargePool; Particle.charge'];
        timePool = [timePool; timePool(end) + fractionalMove*deltaT];    
        eulerRestart = true;
    else
        % Save Intermediate Results
        timeIteration = timeIteration + 1;
        positionPool = [positionPool;reshape(ParticleNext.position,1,3*Particle.nParticle)];
        quaternionPool = [quaternionPool;reshape(ParticleNext.orientation,1,4)];
        chargePool = [chargePool; Particle.charge'];
        timePool = [timePool; timePool(end) + deltaT];
        eulerRestart = false;
    end
        
    ParticleOld = Particle;
    Particle = ParticleNext;
    
    fprintf('%d of %d iterations.\n',timeIteration,numberSteps);
end

fprintf('\ntotal runtime for simulation was %fsec\n', toc);


%% Save Data
save('results/linkedSphereVariables.mat');


%% PLOT THE RESULTS
%--------------------------------------------------------------------------

% Plot 1: plots the 3D trajectory of particles
plot3(positionPool(:,1),positionPool(:,2),positionPool(:,3),'-o'); hold on;
axis([0 Parameter.L1 0 Parameter.L2 0 Parameter.L3]);
grid on;
plot3(positionPool(:,4),positionPool(:,5),positionPool(:,6),'r-o'); hold off

% Plot 2: plots the trajectory of the center of all the particles
figure;
plot(timePool,positionPool(:,3),'-o',timePool,positionPool(:,6),'r-o');

% Ribbon visulazation of trajectory
sol.x= timePool';
sol.y= [mean(positionPool(:,1:3:end),2)'; mean(positionPool(:,2:3:end),2)';...
    mean(positionPool(:,3:3:end),2)'; quaternionPool'];

% Plot 3: plots the 3D/2D representation of spheres in the box
pseudoChargePool = zeros(size(positionPool,1),Particle.nParticle);
pseudoChargePool(:,find(Particle.conductivityIndex==1)) = chargePool;
figure;
M(size(positionPool,1)) = struct('cdata',[],'colormap',[]);
for i = 1:size(positionPool,1)
    draw(reshape(positionPool(i,:),3,Particle.nParticle), ...
        sum(pseudoChargePool(i,:))*Particle.conductivityIndex, Parameter);
   
    for j = 1:Particle.nParticle
    %plots the 3D trajectory of particles
    plot3(positionPool(1:i,1+3*(j-1)),positionPool(1:i,2+3*(j-1)),...
        positionPool(1:i,3+3*(j-1)),'.','MarkerSize',5); hold on; 
    end
    drawnow;
    M(i) = getframe(gcf);
    hold off;
end
% Create AVI file.
VideoWriter(M, 'linkedSpheres3.avi');


%% Ribbon 3d representaion with the cluster movement
% trajectoryVisualize(sol);
% draw02(reshape(positionPool(1,:),3,Particle.nParticle), sum(pseudoChargePool(1,:))*Particle.conductivityIndex, Parameter);
% draw02(reshape(positionPool(size(positionPool,1),:),3,Particle.nParticle), sum(pseudoChargePool(size(positionPool,1),:))*Particle.conductivityIndex, Parameter);
% camlight
% lighting phong
% axis equal


%% Clean path
rmpath( fullfile( rootPath, 'src', 'ConfinedElectrodynamics'));
rmpath( fullfile( rootPath, 'src', 'Hydrodynamics'));
