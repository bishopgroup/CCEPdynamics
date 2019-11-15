function Particle = CalcPos(Particle)

rotationMatrix = quaternionToRotationMatrix(Particle.orientation);

Particle.position = bsxfun(@plus,Particle.clusterTranslation,rotationMatrix'*Particle.bodyFixedPosition);

 %plot(Particle.position(1,:),Particle.position(2,:),'o');
 %axis equal;
 %pause;
