function [M] = M2PlaneFF(Xi,H)
%{
Computes the far-field approximation to the mobility matrix following the
approach described in Swan and Brady, Phys. Fluid, 22, 103301 (2010).

Lengths are scaled by particle radius, a.
Forces, Torques, etc are scaled by 6*pi*nu*a^n with n = 1,2,etc.

Input: 
  Xi - h / H, height of the sphere above the lower wall scaled by H
  H  - distance between parallel walls (scaled by a)

Output:
  M - 15x15 mobility matrix
%}

global muf_f1 muf_f3 muf_f5 Iuf_f muf_g1 muf_g3 muf_g5 Iuf_g ...
       mul_f2 mul_f4 Iul_f...
       mus_f2 mus_f4 mus_f6 Ius_f mus_g2 mus_g4 mus_g6 Ius_g ...
       mos_f3 mos_f5 Ios_f ...
       mol_f3 Iol_f mol_g3 Iol_g ...
       mes_f3 mes_f5 mes_f7 Ies_f mes_g3 mes_g5 mes_g7 Ies_g mes_h3 mes_h5 mes_h7 Ies_h ...
       ABmat

   
%% H powers
H2 = H^2; H3 = H*H2; H4 = H*H3; H5 = H*H4; H6 = H*H5; H7 = H*H6;
   
%% MUF (Equation 35, 3x3 matrix)
MUF = Iuf_f * (1 - muf_f1(Xi)/H + muf_f3(Xi)/H3 - muf_f5(Xi)/H5) ...
    + Iuf_g * (1 - muf_g1(Xi)/H + muf_g3(Xi)/H3 - muf_g5(Xi)/H5);

%% MUL (Equation 36, 3x3 matrix)
MUL = -Iul_f * (mul_f2(Xi)/H2 + mul_f4(Xi)/H4); % note sign

%% MOL (Equation 38, 3x3 matrix)
MOL = Iol_f * (0.75 - mol_f3(Xi)/H^3) + Iol_g * (0.75 - mol_g3(Xi)/H^3); 

%% MUS (Equation 37, 3x3x3 matrix reshaped into 3x9)
MUS = reshape( Ius_f * (-mus_f2(Xi)/H2 + mus_f4(Xi)/H4 - mus_f6(Xi)/H6) ...
             + Ius_g * (mus_g2(Xi)/H2 - mus_g4(Xi)/H4 + mus_g6(Xi)/H6), 3,9);

%% MOS (Equation 39, 3x3x3 matrix reshaped into 3 x 9
MOS = reshape( Ios_f * (-mos_f3(Xi)/H3 + mos_f5(Xi)/H5), 3,9); % note sign

%% MES (equation 40, 3x3x3x3 matrix reshaped into 9 x 9)
MES = reshape( Ies_f * (-9/30 - mes_f3(Xi)/H3 + mes_f5(Xi)/H5 - mes_f7(Xi)/H7)...
    + Ies_g * (9/20 - mes_g3(Xi)/H3 + mes_g5(Xi)/H5 - mes_g7(Xi)/H7)...
    + Ies_h * (9/20 - mes_h3(Xi)/H3 + mes_h5(Xi)/H5 - mes_h7(Xi)/H7), 9,9); 

%% Putting it all together
M = [  MUF,  MUL, MUS ; ...
      MUL',  MOL, MOS ; ...
      MUS', MOS', MES ];
  
M = ABmat * M * ABmat';
M = M(1:11,1:11);


     
