function TensorProducts()
%{
This function test the equivalence of different methods for computing
tensor products in MATLAB.  The simplest approach to computing tensor
products uses multidimensional arrays and nested for-loops to compute the
products using Einstein notation.  It will be convenient, however, to
translate these tensor products into matrix-vector products.
%}

%% Number of Iterations for timing
N = 1e5;

%% (3) * (2) = (1)
S = rand(3,3); % second order tensor
MUS = rand(3,3,3); % third order tensor

% Using index notation
tic
for n = 1:N
    P1 = zeros(3,1);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                P1(i) = P1(i) + MUS(i,j,k)*S(j,k);
            end
        end
    end
end
fprintf('\n(3)*(2) = (1), index  : %f\n',toc);

% Using vector operations
tic
for n = 1:N
    P2 = reshape(MUS,3,9) * reshape(S,9,1);
end
fprintf('(3)*(2) = (1), vector : %f\n',toc);
fprintf('(3)*(2) = (1), error  : %.2e\n',norm(P1 - P2));

%% (4) * (2) = (2)
S = rand(3,3); % second order tensor
MES = rand(3,3,3,3); % fourth order tensor

% Using index notation
tic
for n = 1:N
    P1 = zeros(3,3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    P1(i,j) = P1(i,j) + MES(i,j,k,l)*S(k,l);
                end
            end
        end
    end
end
fprintf('\n(4)*(2) = (2), index  : %f\n',toc);

tic
for n = 1:N
    P2 = reshape(reshape(MES,9,9) * reshape(S,9,1),3,3);
end
fprintf('(4)*(2) = (2), vector : %f\n',toc);
fprintf('(4)*(2) = (2), error  : %.2e\n',norm(P1 - P2));
