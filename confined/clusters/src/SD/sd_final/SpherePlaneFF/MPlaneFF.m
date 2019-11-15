function [M] = MPlaneFF(h)
%{
Computes the far-field approximation to the mobility matrix for a single
sphere near a plane wall following the approach described in Swan and
Brady, Phys. Fluid, 19, 113306 (2007).

Lengths are scaled by particle radius, a.
Forces, Torques, etc are scaled by 6*pi*nu*a^n with n = 1,2,etc.

Input: 
  h - height of the sphere above the plane wall

Output:
  M - 11x11 symmetric, mobility matrix

  M = [ MUF, MUL, MUS;
        MOF, MOL, MOS;
        MEF, MEL, MES]
%}

global Iuf1 Iuf2 MUFinf ...
       Iof1...
       Iol1 Iol2 MOLinf ...
       Ief1 Ief2 Ief3 ...
       Iel1 ...
       Ies1 Ies2 Ies3 Ies4 Ies5 MESinf...
       ABmat
   
%% h powers
h2 = h^2; h3 = h*h2; h4 = h*h3; h5 = h*h4; h6 = h*h5; h7 = h*h6;

%% MUF (Equation B1, 3x3)
MUF = Iuf1 * (-1/16)*(9/h - 2/h3 + 1/h5) ...
    + Iuf2 * (-1/8)*(9/h - 4/h3 + 1/h5);
MUF = MUF + MUFinf;

%% MOF = MUL' (Equation B2, 3x3)
MOF = Iof1 * (3/32)/h4; 

%% MOL (Equation B3, 3x3)
MOL = Iol1 * (-15/64)/h3 + Iol2 * (-3/32)/h3;
MOL = MOL + MOLinf;

%% MEF (Equation B4, 3x3x3)
MEF = reshape( ...
      Ief1 * (-3/160)*(15/h2 - 12/h4 + 5/h6) ...
    + Ief2 * (3/32)*(3/h2 - 3/h4 + 1/h6) ...
    + Ief3 * (-3/16)*(3/h2 - 3/h4 + 1/h6) , 9,3); 

%% MEL (Equation B5, 3x3x3)
MEL = reshape( Iel1 * (9/320)*(5/h3 - 4/h5) , 9,3); % note sign

%% MES (Equation B5, 3x3x3x3)
MES = reshape( ...
      Ies1 * (-3/640)*(10/h3 - 24/h5 + 9/h7) ...
    + Ies2 * (-9/640)*(10/h3 - 8/h5 + 3/h7) ...
    + Ies3 * (3/160)*(20/h3 - 24/h5 + 9/h7) ...
    + Ies4 * (-3/160)*(20/h3 - 24/h5 + 9/h7) ... % MODIFIED
    + Ies5 * (-3/80)*(20/h3 - 24/h5 + 9/h7), 9,9);
% MES = reshape( ...
%       Ies1 * (-3/640)*(10/h3 - 24/h5 + 9/h7) ...
%     + Ies2 * (-9/640)*(10/h3 - 8/h5 + 3/h7) ...
%     + Ies3 * (3/160)*(20/h3 - 24/h5 + 9/h7) ...
%     + Ies4 * (-9/320)*(15/h3 - 16/h5 + 6/h7) ... % ORIGINAL
%     + Ies5 * (-3/80)*(20/h3 - 24/h5 + 9/h7), 9,9);
MES = MES + MESinf;

% Putting it all together
M = [  MUF, MOF', MEF' ; ...
       MOF,  MOL, MEL' ; ...
       MEF,  MEL, MES ];
   
M = ABmat * M * ABmat';
M = M(1:11,1:11);
