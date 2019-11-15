function quaternion3 = quaternionMultiplication(quaternion1,quaternion2)
%% equations of motion

% Diebel, equations 111
conjQmatrix =  [quaternion2(1), -quaternion2(2), -quaternion2(3), -quaternion2(4); ...
            quaternion2(2), quaternion2(1),  -quaternion2(4), quaternion2(3); ...
            quaternion2(3), quaternion2(4),  quaternion2(1),  -quaternion2(2); ...
            quaternion2(4), -quaternion2(3), quaternion2(2),  quaternion2(1)];
         
% Compute quaternion rates (angular velocity)
quaternion3 = (conjQmatrix * quaternion1')';  