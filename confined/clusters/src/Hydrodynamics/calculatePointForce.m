function Particle = calculateVelocityHydro(Parameter,Particle)   
%{

%}

% transpose particle fields
Particle.position = Particle.position';
tempForce = Particle.force;
tempTorque = Particle.torque;
Parameter.domainLength = [Parameter.L1 Parameter.L2 Parameter.L3]';

% Initialize 
Parameter.splittingAlpha = (1/4) * Parameter.domainLength(3).^2;
alpha = Parameter.splittingAlpha;
x = linspace(4*pi/sqrt(alpha),100,1000);
G1 = 1.5 * gammainc(x,0.5,'upper') * gamma(0.5) ./ x.^0.5 / sqrt(alpha);
G2 = -1.5 * gammainc(x,1.5,'upper') * gamma(1.5) ./ x.^0.5 / sqrt(alpha);
Parameter.G1interp = griddedInterpolant(x,G1);
Parameter.G2interp = griddedInterpolant(x,G2);
    
% Initialize particle force / torque / stresslets
%% Force - unit force in the negative z-direction
Particle.force = zeros(Particle.nParticle, 3);

%% Torque and Stresslet
Particle.torque = zeros(Particle.nParticle, 3);
Particle.dipoleHydro = zeros(Particle.nParticle, 9);% change dipole and farfieldDipole

tempForce
tempTorque
conversionFunction(Particle)
pointForce = gmres(conversionFunction(Particle),[tempForce;tempTorque]);
conversionFunction(Particle)*pointForce
%pointForce(isnan(pointForce))= 0;
pointForce
%pause;
Particle.force = reshape(pointForce,3,Particle.nParticle)';

% Initialize particle mesh
Mesh = initializeMeshHydrodynamic(Parameter);
    
% include nearfield 
Parameter.isNearFieldHydrodynamics = false;
    
% Compute velocity
Particle = computeVelocity(Parameter, Particle, Mesh);

% convert back the point velocities into cluster velocity
rFromCOM = Particle.position'- repmat(sum(Particle.position,1)/Particle.nParticle,Particle.nParticle,1)';
rFromCOM
rSquared = repmat(sqrt(sum(rFromCOM.*rFromCOM)),3,1).^2;
rSquared
Particle.angularVelocity = sum(cross(rFromCOM, Particle.velocity')./rSquared ,2)/Particle.nParticle;
Particle.angularVelocity(isnan(Particle.angularVelocity))=0;
Particle.angularVelocity
Particle.velocity = (sum(Particle.velocity,1)/Particle.nParticle)';
Particle.velocity

% transpose back the particle fields
Particle.position = Particle.position';
Particle.force = tempForce;
Particle.torque = tempTorque;

end

function pointforceConversionMatrix = conversionFunction(Particle)

    concatMatrix = zeros(3);% considering the first particle to be origin about which system rotates
    centerOfMass = sum(Particle.position,1)/Particle.nParticle;
     for iParticle = 1:Particle.nParticle
         
         y= Particle.position(iParticle,:)- centerOfMass;
         % cross product matrix
         newMatrix = [0 -y(3) y(2);
                      y(3) 0 -y(1);
                     -y(2) y(1) 0];
         
         concatMatrix = horzcat(concatMatrix, newMatrix);
     end
     
     concatMatrix = concatMatrix(:,4:end);
     pointforceConversionMatrix = [repmat(eye(3),1,Particle.nParticle);concatMatrix];
end