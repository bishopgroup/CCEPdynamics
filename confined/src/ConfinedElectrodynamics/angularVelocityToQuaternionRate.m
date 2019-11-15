function quaternionRate = angularVelocityToQuaternionRate(quaternion,W)
%% equations of motion

% Preallocate
quaternionRate = zeros(4,1);% quaternion rate R3->R4
quaternion = quaternion/sqrt(sumsqr(quaternion));% unit quaternion

% Diebel, equations 156 & 109
Qmatrix = [-quaternion(2), -quaternion(3), -quaternion(4); ...
            quaternion(1),  quaternion(4), -quaternion(3); ...
           -quaternion(4),  quaternion(1),  quaternion(2); ...
            quaternion(3), -quaternion(2),  quaternion(1)];
         
% Compute quaternion rates (angular velocity)
quaternionRate = 0.5 * Qmatrix * W;  
