function nOverlap = checkOverlap( Parameter, Particle )
% brute-force evaluation of particle overlap

idx = Particle.position(3,:) ~= 0; % exclude these parameters
Particle.position = Particle.position(:,idx)';
Particle.position

% compute components of the separation vector
deltaX = abs( bsxfun(@minus, Particle.position(:,1), ...
    Particle.position(:,1)') );
deltaY = abs( bsxfun(@minus, Particle.position(:,2), ...
    Particle.position(:,2)') );
deltaZ = abs( bsxfun(@minus, Particle.position(:,3), ...
    Particle.position(:,3)') );

% correct for periodic boundaries
deltaX = deltaX - ...
    round(deltaX / Parameter.domainLength(1)) * Parameter.domainLength(1); 
deltaY = deltaY - ...
    round(deltaY / Parameter.domainLength(2)) * Parameter.domainLength(2); 


% check for overlaps
distanceSquared = deltaX.^2 + deltaY.^2 + deltaZ.^2;
nOverlap = sum(distanceSquared(:) < 4) - Particle.nParticle ;
