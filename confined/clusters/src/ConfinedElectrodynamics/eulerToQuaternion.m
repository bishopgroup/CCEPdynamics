function quaternion = eulerToQuaternion(eulerAngle)
% eulerToQuaternion   converts from (3,1,3) Euler angles to quaternion.
%
% quaternion = eulerToQuaternion(eulerAngle) takes the Nx3 array eulerAngle
% of Euler angles and computes the Nx4 array quaternion. See eqn (441) of
% J. Diebel, Matrix 58, 1-35 (2006).
%
% Example
%   eulerAngle  = [2*pi, pi, 2*pi] .* rand(1,3);
%   quaternion = eulerToQuaternion(eulerAngle);
%
% Kyle Bishop, June 26, 2015
%

cosPhi = cos(eulerAngle(:,1) / 2);
sinPhi = sin(eulerAngle(:,1) / 2);
cosTheta = cos(eulerAngle(:,2) / 2);
sinTheta = sin(eulerAngle(:,2) / 2);
cosPsi = cos(eulerAngle(:,3) / 2);
sinPsi = sin(eulerAngle(:,3) / 2);

quaternion = [cosPhi .* cosTheta .* cosPsi ...
            - sinPhi .* cosTheta .* sinPsi, ...
              cosPhi .* sinTheta .* cosPsi ...
            + sinPhi .* sinTheta .* sinPsi, ...
              cosPhi .* sinTheta .* sinPsi ...
            - sinPhi .* sinTheta .* cosPsi, ...
              cosPhi .* cosTheta .* sinPsi ...
            + sinPhi .* cosTheta .* cosPsi];