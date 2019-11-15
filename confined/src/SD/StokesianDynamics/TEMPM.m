X11A=1.;
X1NA=1.5/DR-1/DR^3;
Y11A=1.;
Y1NA=0.75/DR+0.5/DR^3;
Y11B=0.;
Y1NB=-3/(4*DR^2);
X11C=0.75;
X1NC=0.75/DR^3;
Y11C=0.75;
Y1NC=-0.375/DR^3;

for i=1:3
for j=1:3

A11(i,j)=X11A*EU(i)*EU(j)+Y11A*(KroneckerDelta(i,j)-EU(i)*EU(j));

A1N(i,j)=X1NA*EU(i)*EU(j)+Y1NA*(KroneckerDelta(i,j)-EU(i)*EU(j));

C11(i,j)=X11C*EU(i)*EU(j)+Y11C*(KroneckerDelta(i,j)-EU(i)*EU(j));

C1N(i,j)=X1NC*EU(i)*EU(j)+Y1NC*(KroneckerDelta(i,j)-EU(i)*EU(j));

end
end

for i=1:3
for j=1:3
DOT=0.;
for k=1:3
    SUM1=Signature(j,i,k)*EU(k);DOT=DOT+SUM1;
end
B11(i,j)=Y11B*DOT;B1N(i,j)=-Y1NB*DOT;
end
end


for i=1:3
A11(i,i)=A11(i,i)-1;
C11(i,i)=C11(i,i)-0.75;
end
