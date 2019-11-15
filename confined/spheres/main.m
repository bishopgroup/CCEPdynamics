%--------------------------------------------------------------------------
% syncCCEP - Computes the dynamics of synchronizing Particles in CCEP.
%
% Input arguments required are :
%
%--------------------------------------------------------------------------

% clear all
fprintf('\nThis might take few minutes...\n\n');


%% Initialize variables and Set Path
rootPath = pwd();
addpath( fullfile( rootPath, 'src', 'Hydrodynamics'));
addpath( fullfile( rootPath, 'src', 'Electrostatics'));
addpath( fullfile( rootPath, 'src', 'export_fig'));
addpath( fullfile( rootPath, 'src', 'pmkmp'));


%% Initialize random number generation
rng(103);


%% Initialize Parameters
Parameter.delta = 0.1; % minimum surface separation
Parameter.includeHydrodynamic = false; % includes the hydrodynamic farfield
Parameter.domainLength = [16, 16, 8];
Parameter.splittingAlpha = pi * (Parameter.domainLength(3) / erfcinv(0.001))^2;
Parameter.nMesh = 2.^3;  % number of mesh points on 0 to L3
Parameter.nParticle = 2;


%% Initialize Mesh, particle position and near-Field
[Particle] = initializeParticlePosition(Parameter);
[Mesh] = initializeMeshElectrostat(Parameter);
NearField = initializeNearField(Parameter);


%% Draw initial config
% fprintf('\nCheck the figure window to verify your input positions and press any key to resume\n\n');
% draw(Particle.position, Particle.charge, Parameter);pause;


%% Integrate over time
% integration Parameters
deltaT = 0.1; % time step(use smaller step size if you have issues at electrode charging as two collisions with same particles can cause problem)
totalTime = 10; % max time(470 at 0.2 deltaT for single chain formation)
nIteration = totalTime / deltaT;  % number of steps

% save Initial Conditions
iIteration = 1;
conductivity = 0;
positionPool = reshape(Particle.position,1,3*Particle.nParticle);
chargePool = Particle.charge';
timePool = 0;
eulerRestart = true; 
iConfigurations = 1;

while iIteration < nIteration
    % Prepare for Integration Step
    Particle = Force(Parameter,Mesh,Particle,NearField);

    if Parameter.includeHydrodynamic
        Particle = calculateVelocity(Parameter,Particle);
    else
        Particle = Velocity(Parameter,Particle);
    end
    
    % Compute conductivity and confidence interval
    totalDipole = sum(Particle.dipole(3,:));
    conductivity = [ conductivity 1 + (4*pi*totalDipole /(Parameter.domainLength(1)*Parameter.domainLength(2)*Parameter.domainLength(3)))];
    
    ParticleNEXT = Particle;

    if eulerRestart % Euler Step
        ParticleNEXT.position = Particle.position + deltaT*Particle.velocity;
    else % 2nd Order Adams-Bashforth Step 
        ParticleNEXT.position = Particle.position + deltaT*(1.5*Particle.velocity - 0.5*ParticleOLD.velocity);
    end
    
    % Check for Overlaps
    [fractionalStep,flag,adjmat] = Overlap(Parameter,ParticleNEXT,Particle);
    if flag % there are overlaps

        % Take a smaller step
        if eulerRestart % Euler Step
            ParticleNEXT.position = Particle.position + fractionalStep*deltaT*Particle.velocity;
        else % 2nd Order Adams-Bashforth Step
            ParticleNEXT.position = Particle.position + fractionalStep*deltaT*(1.5*Particle.velocity - 0.5*ParticleOLD.velocity);
        end
        
        % Equilibrate Charges
        ParticleNEXT = Equilibrate(Parameter,Mesh,ParticleNEXT,NearField,adjmat);

        % Save Intermediate Results
        iIteration = iIteration + 1;
        positionPool = [positionPool;reshape(ParticleNEXT.position,1,3*Particle.nParticle)];
        chargePool = [chargePool; Particle.charge'];
        timePool = [timePool; timePool(end) + fractionalStep*deltaT];
        eulerRestart = true;
    else
        % Save Intermediate Results
        iIteration = iIteration + 1;
        positionPool = [positionPool;reshape(ParticleNEXT.position,1,3*Particle.nParticle)];
        chargePool = [chargePool; Particle.charge'];
        timePool = [timePool; timePool(end) + deltaT];
        eulerRestart = false;
    end
        
    ParticleOLD = Particle;
    Particle = ParticleNEXT;
    
    if mod(iIteration, 10) == 0
        fprintf('%d of %d iterations.\n',iIteration,nIteration);
    end
end

%% Save workspace variables
save('results/allVariables.mat')

%% Postprocess the data
%--------------------------------------------------------------------------
%% Colors
cmap = pmkmp(5,'CubicYF');
set(groot,'defaultAxesColorOrder',cmap);


%% Fonts
set(0,'defaultAxesFontName', 'TimesNewRoman')
    

%% Plot
for i = 1:Particle.nParticle
    plot3(positionPool(:,3*(i-1)+1),positionPool(:,3*(i-1)+2),positionPool(:,3*(i-1)+3),'-o','color',rand(1,3)); hold on;
    axis([0 Parameter.domainLength(1) 0 Parameter.domainLength(2) 0 Parameter.domainLength(3)]);
    grid on
end
hold off
% export
export_fig('results/particleTrajectories.pdf','-nocrop');


figure;
for i = 1:Particle.nParticle
    plot(timePool,positionPool(:,3*(i-1)+3),'-o','color',rand(1,3)); hold on;
end
hold off
title('Time series plot','Interpreter','LaTex','FontSize', 16);
xlabel('time(s)','Interpreter','LaTex','FontSize', 16);
ylabel('$x_3 / a$','Interpreter','LaTex','FontSize', 16);
% export
export_fig('results/timeSeriesParticleHeight.pdf','-nocrop');


% Plot 3: plots the 3D/2D representation of spheres in the box
figure;
M(size(positionPool,1)) = struct('cdata',[],'colormap',[]);
for i = 1:size(positionPool,1)
    draw(reshape(positionPool(i,:),3,Particle.nParticle), chargePool(i,:), Parameter);
    M(i) = getframe(gcf);
    hold off;
end
% Create AVI file.
movie2avi(M, 'movie.avi', 'compression', 'None');


%% Clean path
rmpath( fullfile( rootPath, 'src', 'Electrostatics'));
rmpath( fullfile( rootPath, 'src', 'Hydrodynamics'));
rmpath( fullfile( rootPath, 'src', 'export_fig'));
rmpath( fullfile( rootPath, 'src', 'pmkmp'));
