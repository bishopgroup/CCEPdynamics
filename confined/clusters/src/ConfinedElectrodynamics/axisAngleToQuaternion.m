function quaternion = axisAngleToQuaternion(axisAngle)
%

cosHalfAngle = cos(axisAngle(:,4) / 2);
sinHalfAngle = sin(axisAngle(:,4) / 2);

quaternion = [cosHalfAngle, ...
              (axisAngle(:,1)./norm(axisAngle(:,1:3))).*sinHalfAngle, ...
              (axisAngle(:,2)./norm(axisAngle(:,1:3))).*sinHalfAngle, ...
              (axisAngle(:,3)./norm(axisAngle(:,1:3))).*sinHalfAngle];