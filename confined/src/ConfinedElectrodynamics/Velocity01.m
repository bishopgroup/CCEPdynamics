function part = Velocity(param,part)

% Compute rotation matrix
Rq = eulerToRotationMatrix(part.ac);

% Compute the resistance matrix
Rnf = RNearField(param,part);

% Rotate resistance tensors into the world frame
RNF = tensorTransform3(Rq, Rnf);

VW = RNF \ [part.F; part.T];

part.V = VW(1:3);
part.W = VW(5);
