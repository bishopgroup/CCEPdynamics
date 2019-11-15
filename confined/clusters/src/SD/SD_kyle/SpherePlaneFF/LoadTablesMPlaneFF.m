function LoadTablesMPlaneFF(varargin)
%{

Load tabulated functions from
  Swan and Brady, Phys. Fluid, 19, 103306 (207)

%}

global Iuf1 Iuf2 MUFinf ...
       Iof1...
       Iol1 Iol2 MOLinf...
       Ief1 Ief2 Ief3 ...
       Iel1 ...
       Ies1 Ies2 Ies3 Ies4 Ies5 MESinf ...
       ABmat Flip

%% AB-matrix (accounts for symmetry of stresslet and strain rate)
ABmat = [1/sqrt(6)        0         0         0 -sqrt(2/3)        0         0         0  1/sqrt(6);
                0  1/sqrt(2)        0  1/sqrt(2)        0         0         0         0         0;
                0         0  1/sqrt(2)        0         0         0  1/sqrt(2)        0         0;
                0         0         0         0         0  1/sqrt(2)        0  1/sqrt(2)        0;
         1/sqrt(2)        0         0         0         0         0         0         0 -1/sqrt(2);
                0  1/sqrt(2)        0 -1/sqrt(2)        0         0         0         0         0;
                0         0  1/sqrt(2)        0         0         0 -1/sqrt(2)        0         0;
                0         0         0         0         0  1/sqrt(2)        0 -1/sqrt(2)        0;
         1/sqrt(3)        0         0         0  1/sqrt(3)        0         0         0  1/sqrt(3)]; % 9x9
ABmat = [eye(6), zeros(6,9); zeros(9,6), ABmat];

%% Flip matrix 
Flip = ones(11,1);
Flip([3,4,5,9,10]) = -1; % flip z-direction
Flip = spdiags(Flip,0,11,11);

%% Unit Vectors & Tensors
d3 = [0 0 1]';
I = eye(3);
eps3jk = [0 1 0; -1 0 0; 0 0 0];

%% UF Coupling (Equation B1)
Iuf1 = [1 0 0; 0 1 0; 0 0 0]; % I - d3d3
Iuf2 = [0 0 0; 0 0 0; 0 0 1]; % d3d3
MUFinf = eye(3);

%% OF Coupling (Equation B2)
Iof1 = eps3jk;

%% OL Coupling (Equation B3)
Iol1 = Iuf1;
Iol2 = Iuf2;
MOLinf = 0.75*eye(3);

%%  EF Coupling (Equation B4)
Ief1 = zeros(3,3,3);
Ief2 = zeros(3,3,3);
Ief3 = zeros(3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            Ief1(i,j,k) = (I(i,k) - d3(i)*d3(k))*d3(j) + ...
                (I(j,k) - d3(j)*d3(k))*d3(i);
            Ief2(i,j,k) = (I(i,j) - d3(i)*d3(j))*d3(k);
            Ief3(i,j,k) = d3(i)*d3(j)*d3(k);
        end
    end
end


%% EL Coupling (Equation B5)
Iel1 = zeros(3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            Iel1(i,j,k) = d3(j)*eps3jk(i,k) + d3(i)*eps3jk(j,k);
        end
    end
end

%% ES Coupling (Equation B6)
Ies1 = zeros(3,3,3,3);
Ies2 = zeros(3,3,3,3);
Ies3 = zeros(3,3,3,3);
Ies4 = zeros(3,3,3,3);
Ies5 = zeros(3,3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                Ies1(i,j,k,l) = (I(i,j) - d3(i)*d3(j))*(I(k,l) - d3(k)*d3(l));
                Ies2(i,j,k,l) = (I(i,k) - d3(i)*d3(k))*(I(j,l) - d3(j)*d3(l)) ...
                    + (I(i,l) - d3(i)*d3(l))*(I(j,k) - d3(j)*d3(k));
                Ies3(i,j,k,l) = (I(i,j) - d3(i)*d3(j))*d3(k)*d3(l) ...
                    + (I(k,l) - d3(k)*d3(l))*d3(i)*d3(j);
                Ies4(i,j,k,l) = (I(i,k) - d3(i)*d3(k))*d3(j)*d3(l) ...
                    + (I(i,l) - d3(i)*d3(l))*d3(j)*d3(k) ...
                    + (I(j,k) - d3(j)*d3(k))*d3(i)*d3(l) ...
                    + (I(j,l) - d3(j)*d3(l))*d3(i)*d3(k);
                Ies5(i,j,k,l) = d3(i)*d3(j)*d3(k)*d3(l);
            end
        end
    end
end

% MESinf (Equation A2 of Durlofsky, Brady, Bossis, J. Fluid Mech., 1987,
% 180, 21)
MESinf = zeros(3,3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                MESinf(i,j,k,l) = 1.5*(9/10)*(d3(i)*d3(j) - I(i,j)/3)*(d3(k)*d3(l) - I(k,l)/3) ...
                    + 0.5*(9/10)*(d3(i)*I(j,l)*d3(k) + d3(j)*I(i,l)*d3(k) + d3(i)*I(j,k)*d3(l) ...
                    + d3(j)*I(i,k)*d3(l) - 4*d3(i)*d3(j)*d3(k)*d3(l)) ...
                    +0.5*(9/10)*(I(i,k)*I(j,l) + I(j,k)*I(i,l) - I(i,j)*I(k,l) + ...
                    d3(i)*d3(j)*I(k,l) + I(i,j)*d3(k)*d3(l) + d3(i)*d3(j)*d3(k)*d3(l)...
                    - d3(i)*I(j,l)*d3(k) - d3(j)*I(i,l)*d3(k) - d3(i)*I(j,k)*d3(l) ...
                    - d3(j)*I(i,k)*d3(l));
            end
        end
    end
end
MESinf = reshape(MESinf,9,9);