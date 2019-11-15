function Particle = calculateVelocity(Parameter,Particle)   
%{

%}

% transpose particle fields
Particle.position = Particle.position';
Particle.force  = Particle.force';

% Initialize 
Parameter.splittingAlpha = (1/4) * Parameter.domainLength(3).^2;
alpha = Parameter.splittingAlpha;
x = linspace(4*pi/sqrt(alpha),100,1000);
G1 = 1.5 * gammainc(x,0.5,'upper') * gamma(0.5) ./ x.^0.5 / sqrt(alpha);
G2 = -1.5 * gammainc(x,1.5,'upper') * gamma(1.5) ./ x.^0.5 / sqrt(alpha);
Parameter.G1interp = griddedInterpolant(x,G1);
Parameter.G2interp = griddedInterpolant(x,G2);
    
% Initialize particle force / torque / stresslets
Particle = initializeForce(Particle);

%disp(Particle.force);pause;

% Initialize particle mesh
Mesh = initializeMeshHydrodynamic(Parameter);
    
% include nearfield ?
Parameter.isNearFieldHydrodynamics = false;
    
% Compute velocity
Particle = computeVelocity(Parameter, Particle, Mesh);

%disp(Particle.velocity);pause;

% transpose back the particle fields
Particle.position = Particle.position';
Particle.force  = Particle.force';
Particle.velocity  = Particle.velocity';