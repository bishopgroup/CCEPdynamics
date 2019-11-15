clear XM;

for i=1:12
for j=1:12
XM(i,j)=0;
end
end

for i=1:3
for j=1:3
XM(i,j)=A11(i,j);
XM(i,j+3)=A1N(i,j);
XM(i+3,j+3)=A11(i,j);
XM(i,j+6)=B11(i,j);
XM(i,j+9)=B1N(i,j);
XM(i+3,j+6)=-B1N(i,j);
XM(i+3,j+9)=-B11(i,j);
XM(i+6,j+6)=C11(i,j);
XM(i+6,j+9)=C1N(i,j);
XM(i+9,j+9)=C11(i,j);
end
end

for i=1:12
for j=1:12
XM(j,i)=XM(i,j);
end
end

for i=1:6
XM(i,i)=XM(i,i)+1;
XM(i+6,i+6)=XM(i+6,i+6)+0.75;
end

XM=inv(XM);
