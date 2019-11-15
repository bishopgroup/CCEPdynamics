function R = RSD(Xi, H)
%RSD Computes resistance tensor coupling velocities with forces/torques
%
%   R = RSD(Xi, H) computes the Stokesian Dynamics approximation of the
%   resistance tensor (6x6 matrix) coupling translational and orientational
%   velocities with forces and toreques , i.e.
%   \[
%     \mathbf{R} = \left[
%        \begin{array}{cc}
%           \mathbf{R}_{FU} & \mathbf{R}_{F\Omega} \\
%           \mathbf{R}_{LU} & \mathbf{R}_{L\Omega}
%        \end{array} 
%   \]
%   for a single sphere between two parallel walls spearated by distance $H$.
%   Position of the sphere is determined by its fractional height
%   $\Xi \in (0; 1)$ above the lower plane.
%
%   References
%   JW Swan and JF Brady, Phys Fluids, 22, 103301 (2010)

global Flip

addpath('Sphere2PlaneFF', 'SpherePlaneFF', 'SpherePlaneNF')

% Calculate the actual height of the particle above the lower wall.
h = Xi * H;

% Calculate far field resitatance matrix, $R^{\infty}$ which is an inverse of
% grand mobility matrix $M^{\infty}$.
Mfar = M2PlaneFF(Xi, H);
Rfar = inv(Mfar);

% Calculate two-body interactions of the particle with both lower and upper
% wall. For each wall, those interactions consists of two components: exact
% pair-wise lubrication interactions, $R^{2B}$, and a two-body, far-field
% correction $R^{2B, \infty}$ which prevents from counting far-field
% interaction twice as they are already included in grand mobility matrix.
R2B_lower = RPlaneNF(h);
M2Bfar_lower = MPlaneFF(h);
R2Bfar_lower = inv(M2Bfar_lower);

R2B_upper = Flip(1:6,1:6)' * RPlaneNF(H-h) * Flip(1:6,1:6);
M2Bfar_upper = Flip' * MPlaneFF(H-h) * Flip;
R2Bfar_upper = inv(M2Bfar_upper);

% Put all toghether to determine $\mathbf{R}_{SD}$.
R = Rfar(1:6,1:6) + ...
    R2B_lower - R2Bfar_lower(1:6,1:6) + ...
    R2B_upper - R2Bfar_upper(1:6,1:6);
end