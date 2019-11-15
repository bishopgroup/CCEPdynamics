function Tw = tensorTransform3(R,Tb)
%% Transform Triadics: Tw_ijk = Tb_pqr * R_pi * R_qj * R_rk
% 
% Tb is a second order tensor in the body-fixed coordinates 
% R is a rotation matrix that maps from the world coordinates to the
% body-fixed coordinates as zb = R * zw
%
% Consider the following equation among vectors & tensors in the body-fixed
% coordinates: yb = Tb * xb (where Tb is a 3rd order tensor, xb is a 2nd 
% order tensor, and yb is a vector).  To convert to world coordinates, we 
% write that R * yw = Tb * R * xw * R' or, equivalently,
%   
%   yw_i = R'_ip * Tb_pqr * (R_qj * xw_jk * R'_kr)
%   yw_i = (Tb_pqr * R_pi * R_qj * R_rk) * xw_jk
%

Tijr = zeros(3,3,3);
Tijr(:,:,1) = R' * Tb(:,:,1) * R; % T_ijr, r = 1
Tijr(:,:,2) = R' * Tb(:,:,2) * R; % T_ijr, r = 2
Tijr(:,:,3) = R' * Tb(:,:,3) * R; % T_ijr, r = 3
Tw = reshape( reshape(Tijr,9,3) * R, 3, 3, 3);
