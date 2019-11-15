function R = eulerToRotationMatrix(eulerAngle)
% eulerToRotationMatrix converts from (3,1,3) Euler angles to rotationMatrix.
%
% quaternion = eulerToQuaternion(eulerAngle) takes the Nx3 array eulerAngle
% of Euler angles and computes the Nx4 array quaternion. See eqn (441) of
% J. Diebel, Matrix 58, 1-35 (2006).
%
% Example
%   eulerAngle  = [2*pi, pi, 2*pi] .* rand(1,3);
%   rotationMatrix = eulerToRotationMatrix(eulerAngle);
%
% Kyle Bishop, June 26, 2015
%

cosPhi = cos(eulerAngle(:,1));
sinPhi = sin(eulerAngle(:,1));
cosTheta = cos(eulerAngle(:,2));
sinTheta = sin(eulerAngle(:,2));
cosPsi = cos(eulerAngle(:,3));
sinPsi = sin(eulerAngle(:,3));

if min(size(eulerAngle)) == 1 % should be either 3x1 or 1x3 
    R = zeros(3,3);
    R(1,1) = cosPhi .* cosPsi - sinPhi .* cosTheta .* sinPsi;
    R(2,1) = -sinPhi .* cosPsi - cosPhi .* cosTheta .* sinPsi;
    R(3,1) = sinTheta .* sinPsi;
    R(1,2) = cosPhi .* sinPsi + sinPhi .* cosTheta .* cosPsi;
    R(2,2) = -sinPhi .* sinPsi + cosPhi .* cosTheta .* cosPsi;
    R(3,2) = -sinTheta .* cosPsi;
    R(1,3) = sinPhi .* sinTheta;
    R(2,3) = cosPhi .* sinTheta;
    R(3,3) = cosTheta;
else % should be Nx3
    R = zeros(size(eulerAngle,1), 9);
    R(:,1) = cosPhi .* cosPsi - sinPhi .* cosTheta .* sinPsi;
    R(:,2) = -sinPhi .* cosPsi - cosPhi .* cosTheta .* sinPsi;
    R(:,3) = sinTheta .* sinPsi;
    R(:,4) = cosPhi .* sinPsi + sinPhi .* cosTheta .* cosPsi;
    R(:,5) = -sinPhi .* sinPsi + cosPhi .* cosTheta .* cosPsi;
    R(:,6) = -sinTheta .* cosPsi;
    R(:,7) = sinPhi .* sinTheta;
    R(:,8) = cosPhi .* sinTheta;
    R(:,9) = cosTheta;
end
