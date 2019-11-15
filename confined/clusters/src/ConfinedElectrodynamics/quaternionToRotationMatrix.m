function R = quaternionToRotationMatrix(q)
% quaternionToRotationMatrix   converts from quaternion to rotation matrix.
%
% R = quaternionToRotationMatrix(q) takes the Nx4 array quaternion
% and computes the Nx9 rotation matrix.  For the special case of N = 1, the 
% rotation matrix is returned as a 3x3 matrix.  For N > 1, the nth rotation
% matrix can be converted to standard 3x3 form as reshape(R(n,:),3,3).
%
% The rotation matrix relates a vector z in the global coordinates to
% the same vector z' in the body-fixed coordinates, z'_i = R_ij * z_j. 
% See eqn (125) of J. Diebel, Matrix 58, 1-35 (2006) for details. 
%
% Example
%   q = rand(3,4);
%   q = bsxfun(@rdivide, q, sqrt(sum(q.^2 , 2))); % normalize
%   R = quaternionToRotationMatrix(q);
%
% Kyle Bishop, June 26, 2015
%

q = q/sqrt(sumsqr(q));% unit quaternion

if min(size(q)) == 1 % should be either 4x1 or 1x4 
    R = zeros(3,3);
    R(1,1) = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
    R(2,1) = 2 * q(2) * q(3) - 2 * q(1) * q(4);
    R(3,1) = 2 * q(2) * q(4) + 2 * q(1) * q(3);
    R(1,2) = 2 * q(2) * q(3) + 2 * q(1) * q(4);
    R(2,2) = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
    R(3,2) = 2 * q(3) * q(4) - 2 * q(1) * q(2);
    R(1,3) = 2 * q(2) * q(4) - 2 * q(1) * q(3);
    R(2,3) = 2 * q(3) * q(4) + 2 * q(1) * q(2);
    R(3,3) = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;
else % should be Nx4
    R = zeros(size(q,1), 9);
    R(:,1) = q(:,1).^2 + q(:,2).^2 - q(:,3).^2 - q(:,4).^2;     % 11
    R(:,2) = 2 * q(:,2) .* q(:,3) - 2 * q(:,1) .* q(:,4);       % 21
    R(:,3) = 2 * q(:,2) .* q(:,4) + 2 * q(:,1) .* q(:,3);       % 31
    R(:,4) = 2 * q(:,2) .* q(:,3) + 2 * q(:,1) .* q(:,4);       % 12
    R(:,5) = q(:,1).^2 - q(:,2).^2 + q(:,3).^2 - q(:,4).^2;     % 22
    R(:,6) = 2 * q(:,3) .* q(:,4) - 2 * q(:,1) .* q(:,2);       % 32
    R(:,7) = 2 * q(:,2) .* q(:,4) - 2 * q(:,1) .* q(:,3);       % 13
    R(:,8) = 2 * q(:,3) .* q(:,4) + 2 * q(:,1) .* q(:,2);       % 23
    R(:,9) = q(:,1).^2 - q(:,2).^2 - q(:,3).^2 + q(:,4).^2;     % 33
end