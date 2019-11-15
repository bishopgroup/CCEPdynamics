function [P,V] =Recurence()
global markP markV

N = 5;
P = zeros(N,N,N);
V = zeros(N,N,N);
P(2,1,1) = 1;
V(2,1,1) = 1;

markP = zeros(N,N,N);
markV = zeros(N,N,N);
markP(:,1,1) = 1;
markV(:,1,1) = 1;

for n = 0:N-1
    for p = n-1:N-1
        for q = 0:n
            [n+1,p+1,q+1]
            for s = 1:q
                [s,q-s,p-n+1]
                P(n+1,p+1,q+1) = P(n+1,p+1,q+1) + nchoosek(n+s,n) * ( ...
                    (n*(2*n+1)*(2*n*s-n-s+2)/(2*(n+1)*(2*s-1)*(n+s)))*GetVal(P,s,q-s,p-n+1) - ...
                    (n*(2*n-1)/(2*(n+1)))*GetVal(P,s,q-s,p-n-1) - ...
                    (n*(4*n^2-1)/(2*(n+1)*(2*s+1)))*GetVal(V,s,q-s-2,p-n+1));
            end
            markP(n+1,p+1,q+1) = 1;
           
            V(n+1,p+1,q+1) = P(n+1,p+1,q+1);
%             markV(n+1,p+1,q+1) = 1
            for s = 1:q
                V(n+1,p+1,q+1) = V(n+1,p+1,q+1) ...
                    - (2*n/((n+1)*(2*n+3))) * nchoosek(n+s,n) * GetVal(P,s,q-s,p-n-1);  
            end
        end
    end
end

function val = GetVal(M,i,j,k)
global markP 
i = i + 1;
j = j + 1;
k = k + 1;

[i,j,k]

[imax,jmax,kmax] = size(M);
if i<1 || i>imax || j<1 || j>jmax || k<1 || k>kmax
    val = 0;
else
    val = M(i,j,k);
    if markP(i,j,k)==0
        disp('here')
        pause;
    end
end
