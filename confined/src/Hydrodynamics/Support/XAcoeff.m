function fk = XAcoeff()
%{
Solves for the resistance function XA for two equally sized spheres
(D. J. Jeffrey and Y. Onishi, J. Fluid Mech., 1984, 139, 261–290.).

N = number of terms in the series solution
%}
N = 10;
%% Comput the P1 matrix (equation 3.6-3.9)
PVsize = [N,N,N];
P = zeros(N,N,N);
V = zeros(N,N,N);
M = N^3;
PVmat = spalloc(2*M,2*M,2*N^4);
% PVmat = zeros(2*M);

linear = reshape(1:M,N,N,N);

c1 = zeros(N,N);
c2 = zeros(N,N);
c3 = zeros(N,N);
c4 = zeros(N,N);
for n = 0:N-1
    for s = 0:N-1
        binomial = nchoosek(n+s,n);
        c1(n+1,s+1) = binomial*(-(2*n/((n+1)*(2*n+3))));
        c2(n+1,s+1) = binomial*(n*(2*n+1)*(2*n*s-n-s+2)/(2*(n+1)*(2*s-1)*(n+s)));
        c3(n+1,s+1) = binomial*(-n*(2*n-1)/(2*(n+1)));
        c4(n+1,s+1) = binomial*(-n*(4*n^2-1)/(2*(n+1)*(2*s+1)));
    end
end

for n = 0:N-1
    for p = 0:N-1
        for q = 0:N-1
            i = linear(n+1, p+1, q+1);
            
            % equation 3.8
            PVmat(i+M,i+M) = -1;
            PVmat(i+M,i) = 1;
            for s = 1:q
                if s >= 0 && q-s >= 0 && p-n-1 >= 0
                    j = linear(s+1, q-s+1, p-n);
                    PVmat(i+M,j) = c1(n+1,s+1);
                end
            end
            
            % equation 3.9
            PVmat(i,i) = -1;
            for s = 1:q
                if s >= 0 && q-s >= 0 && p-n+1 >= 0 && p-n+1 < N
                    j = linear(s+1, q-s+1, p-n+2);
                    PVmat(i,j) = c2(n+1,s+1);
                end
                if s >= 0 && q-s >= 0 && p-n-1 >= 0
                    j = linear(s+1, q-s+1, p-n);
                    PVmat(i,j) = c3(n+1,s+1);
                end
                if s >= 0 && q-s-2 >= 0 && p-n+1 >= 0 && p-n+1 < N
                    j = linear(s+1, q-s-1, p-n+2);
                    PVmat(i,j+M) = c4(n+1,s+1);
                end
            end
        end
    end
end
% sum(PVmat(:)~=0)

% equation 3.6
rhs = zeros(2*M,1);
i = linear(2, 1, 1);
rhs(i) = -1;

% solve for coefficients
PV = PVmat \ rhs;
P = reshape(PV(1:M),N,N,N);
P1 = reshape(P(2,:,:),N,N);

%% fk coefficients (lambda = 1)
fk = zeros(N,1);
for k = 0:N-1
    tmp = 0;
    for q = 0:k
        tmp = tmp + P1(k-q+1,q+1);
    end
    fk(k+1) = 2^k * tmp;
end
