function Particle = calculateForce(Parameter,mesh,Particle)
%{
Compute the electrostatic free energy of a collection of charged spheres.
%}

Particle = Energy(Parameter,mesh,Particle);

%% Force
delta = 1e-2; % step size for finite differences

part_dx1 = Particle;
part_dx1.clusterTranslation(1) = part_dx1.clusterTranslation(1) + delta;
part_dx1 = CalcPos(part_dx1);
part_dx1 = Energy(Parameter,mesh,part_dx1);

part_dx2 = Particle;
part_dx2.clusterTranslation(2) = part_dx2.clusterTranslation(2) + delta;
part_dx2 = CalcPos(part_dx2);
part_dx2 = Energy(Parameter,mesh,part_dx2);

part_dx3 = Particle;
part_dx3.clusterTranslation(3) = part_dx3.clusterTranslation(3) + delta;
part_dx3 = CalcPos(part_dx3);
part_dx3 = Energy(Parameter,mesh,part_dx3);

Particle.force = [(Particle.Ues - part_dx1.Ues) / delta;
                  (Particle.Ues - part_dx2.Ues) / delta;
                  (Particle.Ues - part_dx3.Ues) / delta];
              
% %% point force for all particles
% for iParticle = 1:Particle.nParticle
%     
%     part_dx1 = Particle; 
%     part_dx1.position(1,iParticle) = part_dx1.position(1,iParticle) + delta;
%     part_dx1 = Energy(Parameter,mesh,part_dx1);
%     
%     part_dx2 = Particle; 
%     part_dx2.position(2,iParticle) = part_dx2.position(2,iParticle) + delta;
%     part_dx2 = Energy(Parameter,mesh,part_dx2);
%     
%     part_dx3 = Particle; 
%     part_dx3.position(3,iParticle) = part_dx3.position(3,iParticle) + delta;
%     part_dx3 = Energy(Parameter,mesh,part_dx3);
% 
%     Particle.F(:,iParticle) = [(Particle.Ues - part_dx1.Ues) / delta;
%                    (Particle.Ues - part_dx2.Ues) / delta;
%                    (Particle.Ues - part_dx3.Ues) / delta];
%     
% end
% Particle.force
% Particle.F
% pause;

      
%% Torque
delta = 1e-2; % step size for finite differences

Rx = [1 0 0;
      0 cos(delta) sin(delta);
     0 -sin(delta)  cos(delta)];
 
Ry = [cos(delta) 0 -sin(delta);
      0 1 0;
     sin(delta) 0  cos(delta)];
 
Rz = [cos(delta) sin(delta) 0;
      -sin(delta) cos(delta) 0;
      0 0 1];

part_da1 = Particle;
part_da1.position = Rx'*(part_da1.position-repmat(part_da1.clusterTranslation,1,Particle.nParticle)) + repmat(part_da1.clusterTranslation,1,Particle.nParticle);
part_da1 = Energy(Parameter,mesh,part_da1);

part_da2 = Particle;
part_da2.position = Ry'*(part_da2.position-repmat(part_da2.clusterTranslation,1,Particle.nParticle)) + repmat(part_da2.clusterTranslation,1,Particle.nParticle);
part_da2 = Energy(Parameter,mesh,part_da2);

part_da3 = Particle;
part_da3.position = Rz'*(part_da3.position-repmat(part_da3.clusterTranslation,1,Particle.nParticle)) + repmat(part_da3.clusterTranslation,1,Particle.nParticle);
part_da3 = Energy(Parameter,mesh,part_da3);

%
% figure
% draw2D(Particle.position, Particle.charge, Parameter);
% figure
% draw2D(part_da3.position, Particle.charge, Parameter);
% figure
% draw2D(part_da4.position, Particle.charge, Parameter);
% pause;
%}

Particle.torque = [(Particle.Ues - part_da1.Ues) / delta;
                   (Particle.Ues - part_da2.Ues) / delta;
                   (Particle.Ues - part_da3.Ues) / delta];
      
               
               
% %% point force for all particles
% for iParticle = 1:Particle.nParticle
%     
%     delta_1 = part_da1.position(:,iParticle) + Particle.position(:,iParticle);
%     delta_2 = part_da2.position(:,iParticle) + Particle.position(:,iParticle);
%     delta_3 = part_da3.position(:,iParticle) + Particle.position(:,iParticle);
%     delta_rot = delta_1 + delta_2 + delta_3;
% 
%     Particle.Frot(:,iParticle) = [(Particle.Ues - part_dx1.Ues) / delta_rot(1);
%                    (Particle.Ues - part_dx2.Ues) / delta_rot(2);
%                    (Particle.Ues - part_dx3.Ues) / delta_rot(3)];
%     
% end



