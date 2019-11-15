for i=1:NDIM
UNEW(i)=0;
end

for i=1:Ng
for j=1:Ng
MOBMATA(i,j)=0;
MOBMATB(i,j)=0;
MOBMATC(i,j)=0;
end
end

for i=1:NDIM
for j=1:NDIM
RL(i,j)=0;
end
end

for in=1:NSPHR-1
for jn=in+1:NSPHR

DX=X(jn,t-1)-X(in,t-1);
DY=Y(jn,t-1)-Y(in,t-1);
DZ=Z(jn,t-1)-Z(in,t-1);

DR=sqrt(DX^2+DY^2+DZ^2);

EU(1)=DX/DR;
EU(2)=DY/DR;
EU(3)=DZ/DR;

TEMPM
FTMOB
SSI

for k=1:3
for l=1:3
                
N1=l+3*(jn-1);
N2=k+3*(jn-1);
N3=k+3*(in-1);
N4=l+3*(in-1);

MOBMATB(N3,N1)=B1N(k,l);
MOBMATB(N2,N1)=-B11(k,l)+MOBMATB(N2,N1);
MOBMATB(N3,N4)=B11(k,l)+MOBMATB(N3,N4);
MOBMATB(N1,N3)=MOBMATB(N3,N1);
MOBMATA(N3,N1)=A1N(k,l);
MOBMATC(N3,N1)=C1N(k,l);
MOBMATA(N2,N1)=A11(k,l)+MOBMATA(N2,N1);
MOBMATC(N2,N1)=C11(k,l)+MOBMATC(N2,N1);
MOBMATA(N3,N4)=A11(k,l)+MOBMATA(N3,N4);
MOBMATC(N3,N4)=C11(k,l)+MOBMATC(N3,N4);

RL(N3,N4)=TEMP(k,l)+RL(N3,N4);
RL(N3,N4+Ng)=TEMP(k,l+6)+RL(N3,N4+Ng);
RL(N3,N1)=TEMP(k,l+3);
RL(N3,N1+Ng)=TEMP(k,l+9);
RL(N1,N3+Ng)=RL(N3,N1+Ng);
RL(N3+Ng,N4+Ng)=TEMP(k+6,l+6)+RL(N3+Ng,N4+Ng);
RL(N3+Ng,N1+Ng)=TEMP(k+6,l+9);
RL(N2,N1)=TEMP(k,l)+RL(N2,N1);
RL(N2,N1+Ng)=-TEMP(k,l+6)+RL(N2,N1+Ng);
RL(N2+Ng,N1+Ng)=TEMP(k+6,l+6)+RL(N2+Ng,N1+Ng);
end
end
end
end

for i=1:Ng
MOBMATA(i,i)=MOBMATA(i,i)+1;
MOBMATC(i,i)=MOBMATC(i,i)+0.75;
end

for i=1:Ng
for j=1:i 
MOBMATA(i,j)=MOBMATA(j,i);
MOBMATC(i,j)=MOBMATC(j,i);
end
end

for i=1:NDIM
for j=1:i
RL(i,j)=RL(j,i);
end
end

for i=1:Ng
for j=1:Ng
A(i,j)=MOBMATA(i,j);
A(i,j+Ng)=MOBMATB(i,j);
A(i+Ng,j+Ng)=MOBMATC(i,j);
end
end

for i=1:NDIM
for j=1:NDIM
A(j,i)=A(i,j);
end
end

Am=inv(A);

for i=1:NDIM
for j=1:NDIM
R(i,j)=Am(i,j)+RL(i,j);
end
end

Amm=inv(R);

for i=1:NDIM
for j=1:NDIM
SUM=Amm(i,j)*FU(j);
UNEW(i)=UNEW(i)+SUM;
end
end
