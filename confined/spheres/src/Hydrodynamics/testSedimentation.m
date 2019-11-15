%--------------------------------------------------------------------------
% testSedimentation - Computes the sedimentation velocity for a simple 
%   cubic lattice.
%
% Input is provided in testSedimentationInput.txt with the following fields:
%
%  nUnitCell: number of unit cells (in each direction) to include in the 
%      simulation domain. NOTE: if the 'local' potential is negligible
%      (when particle spacing >> sqrt(alpha/pi)), then 1 unit cell may
%      suffice; otherwise, 2 or more are required.
%
%  meshSpacingFraction: Sets the desired spacing between mesh points as a
%      prescribed fraction of sqrt(alpha/pi).  The meshSpacing is then
%      adusted to by domainLength / nMesh where nMesh is a power of two.
%--------------------------------------------------------------------------

clear all

%% Read input parameters
Parameter = parseSedimentationInput('testSedimentationInput.txt');

%% Brady Table 1 - Particles on cubic lattice
tic
fprintf('\nSimple cubic lattice:\n');
fprintf('  No. of particles =\t%d\n',Parameter.nUnitCell^3 * 1);

% Get Bonnecaze data
%tempData = importdata('BradyTable1.txt','\t',0);
tempData = load('BradyTable1.txt');
volumeFraction = tempData(:,1);
table1Data = tempData(:,2:end);

% Generate new data
myTable1Data = zeros(size(volumeFraction));


for iVolumeFraction = 1:length(volumeFraction)
    
    % Initialize particles
    Parameter.volumeFraction = volumeFraction(iVolumeFraction);
    [Parameter, Particle] = initializeParticleCubic(Parameter);
    
    % Initialize particle force / torque / stresslets
    %% Force - unit force in the negative z-direction
    Particle.force = zeros(Particle.nParticle, 3);
    Particle.force(:,3) = -1;
    Particle = initializeForce(Particle);
    
    % Initialize particle mesh
    Mesh = initializeMeshHydrodynamic(Parameter);
    
    % Initialize nearfield
    Parameter.isNearFieldHydrodynamics = false;
    
    % Compute velocity
    Particle = computeVelocity(Parameter, Particle, Mesh);
    
    % Add to table
    myTable1Data(iVolumeFraction) = -mean(Particle.velocity(:,3));
    
%     if iVolumeFraction == 4
%         break
%     end
end


% Plot the results
h(1) = plot(volumeFraction,table1Data(:,2),'k.-'); hold on;
h(2) = plot(volumeFraction,myTable1Data,'ro');

maxVolumeFraction = (4*pi/3) / 2^3;
plot(maxVolumeFraction*[1,1],[0,10],'k--');
hold off;

%ylim([-1,1])

title('Simple Cubic Lattice','Interpreter','LaTex','FontSize', 16);
xlabel('Volume Fraction','Interpreter','LaTex','FontSize', 16);
ylabel('$U / U_0$','Interpreter','LaTex','FontSize', 16);

set(gca,'FontName','Times New Roman',...
    'FontSize',12,...
    'XScale','log');

legend(h,{'Hashimoto', ...
    'My SD'},...
    'Interpreter','LaTex',...
    'Location','NorthWest', ...
    'FontSize', 14);
legend('boxoff');

fprintf('  Run time is %fs\n',toc);
pause(0.1);

