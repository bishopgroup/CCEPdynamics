function eulerAngle = quaternionToEuler(quaternion)
% quaternionToEuler   converts from quaternion to (3,1,3) Euler angles.
%
% eulerAngle = quaternionToEuler(quaternion) takes the Nx4 array quaternion
% and computes the Nx3 array eulerAngle. See eqn (441) of J. Diebel, Matrix
% 58, 1-35 (2006).
%
% Example
%   quaternion  = rand(1,4);
%   eulerAngle = quaternionToEuler(quaternion);
%
% Kyle Bishop, June 26, 2015
%

eulerAngle = zeros(size(quaternion,1), 3);

eulerAngle(:,1) = atan2(2 * quaternion(:,2) .* quaternion(:,4) ...
    - 2 * quaternion(:,1) .* quaternion(:,3), ...
      2 * quaternion(:,3) .* quaternion(:,4) ...
    + 2 * quaternion(:,1) .* quaternion(:,2));
eulerAngle(:,2) = acos( quaternion(:,4).^2 - quaternion(:,3).^2 ...
    - quaternion(:,2).^2 + quaternion(:,1).^2);
eulerAngle(:,3) = atan2(2 * quaternion(:,2) .* quaternion(:,4) ...
    + 2 * quaternion(:,1) .* quaternion(:,3), ...
    - 2 * quaternion(:,3) .* quaternion(:,4) ...
    + 2 * quaternion(:,1) .* quaternion(:,2));

eulerAngle = real(eulerAngle);

