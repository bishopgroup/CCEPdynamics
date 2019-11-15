if (DR>4.)

X11AR= 1.+2.25*1/DR^2+5.8125*1/DR^4+18.703*1/DR^6+77.426*1/DR^8;

Y11AR= 1.+0.5625*1/DR^2+1.8164*1/DR^4+3.5999*1/DR^6+8.6978*1/DR^8;

Y11BR= -(0.75*1/DR^3+1.4219*1/DR^5+5.7373*1/DR^7+20.708*1/DR^9);

Y11CR= 1.33333+1/DR^4+5.8958*1/DR^6+15.566*1/DR^8;

X11CR= 1.33333 + 1.33333*1/DR^6 + 4.*1/DR^8;

X12AR= -(1.5*1/DR+2.375*1/DR^3+12.094*1/DR^5+ ...
41.648*1/DR^7+148.66*1/DR^9);

Y12AR= -(0.75*1/DR+0.92188*1/DR^3+ ...
2.2061*1/DR^5+5.4714*1/DR^7+16.981*1/DR^9);

Y12BR= 1/DR^2+0.5625*1/DR^4+2.1914*1/DR^6+10.498*1/DR^8;

X12CR= -1.333333*1/DR^3-1.33333*1/DR^9;

Y12CR= 0.666667*1/DR^3+0.75*1/DR^5+1.9219*1/DR^7+1.9352*1/DR^9;

end


if (DR<2.02)
        XI = DR - 2.;
        XI1 = 1./XI;
        DLX = Log10(XI1);         
        XDLX = XI*DLX;
        X11AR= 0.25*XI1+0.225*DLX+0.99538+0.026786*XDLX;
        Y11AR= 0.166667*DLX+0.9983;
        Y11BR= -0.166667*DLX+0.1594-0.083333*XDLX;
        X11CR= 1.40240 + 0.166667*XDLX;
        Y11CR= 0.266667*DLX+0.9372+0.25067*XDLX;
        X12AR= -X11AR + 0.64518;
        Y12AR= -Y11AR + 0.7246;
        Y12BR= -Y11BR + 0.1583;
        Y12CR= 0.0666667*DLX-0.03653 + 0.5*0.16533*XDLX;
        X12CR= -0.200343+0.166667*XDLX;
        Y11AR = Y11AR + 0.1108*XI;
        Y12AR = Y12AR - 0.02365*XI;
        Y11BR = Y11BR - 0.1929*XI;
        Y12BR = Y12BR + 0.09418*XI;
end

if (4>DR & DR>2.1)
	    IB = 7+floor(20.*(DR-2.));
        IA=IB +1;
        C1=(DR-RSS(IB))/(RSS(IA)-RSS(IB));
        X11AR= (X11AS(IA)-X11AS(IB))*C1+X11AS(IB);
        Y11AR= (Y11AS(IA)-Y11AS(IB))*C1+Y11AS(IB);
        Y11BR= (Y11BS(IA)-Y11BS(IB))*C1+Y11BS(IB);
        X11CR= (X11CS(IA)-X11CS(IB))*C1+X11CS(IB);
        Y11CR= (Y11CS(IA)-Y11CS(IB))*C1+Y11CS(IB);
        X12AR= (X12AS(IA)-X12AS(IB))*C1+X12AS(IB);
        Y12AR= (Y12AS(IA)-Y12AS(IB))*C1+Y12AS(IB);
        Y12BR= (Y12BS(IA)-Y12BS(IB))*C1+Y12BS(IB);
        X12CR= (X12CS(IA)-X12CS(IB))*C1+X12CS(IB);
        Y12CR= (Y12CS(IA)-Y12CS(IB))*C1+Y12CS(IB);
end
        
if (2.02<DR & DR<2.1)
        IB = -1 + floor(100.*(DR - 2.));
	    IA=IB +1;
        C1=(DR-RSS(IB))/(RSS(IA)-RSS(IB));
        X11AR= (X11AS(IA)-X11AS(IB))*C1+X11AS(IB);
        Y11AR= (Y11AS(IA)-Y11AS(IB))*C1+Y11AS(IB);
        Y11BR= (Y11BS(IA)-Y11BS(IB))*C1+Y11BS(IB);
        X11CR= (X11CS(IA)-X11CS(IB))*C1+X11CS(IB);
        Y11CR= (Y11CS(IA)-Y11CS(IB))*C1+Y11CS(IB);
        X12AR= (X12AS(IA)-X12AS(IB))*C1+X12AS(IB);
        Y12AR= (Y12AS(IA)-Y12AS(IB))*C1+Y12AS(IB);
        Y12BR= (Y12BS(IA)-Y12BS(IB))*C1+Y12BS(IB);
        X12CR= (X12CS(IA)-X12CS(IB))*C1+X12CS(IB);
        Y12CR= (Y12CS(IA)-Y12CS(IB))*C1+Y12CS(IB);
end

for i=1:3
for j=1:3
DOT=0;
for k=1:3
SUM1 = Signature(j,i,k)*EU(k);DOT=DOT+SUM1;
end
TEMP(i,j+6) = Y11BR*DOT;
TEMP(i,j+9) =-Y12BR*DOT;
TEMP(i+3,j+6) = -TEMP(i,j+9);
TEMP(i+3,j+9) = -TEMP(i,j+6);
TEMP(i,j) = X11AR*EU(i)*EU(j)+Y11AR*(KroneckerDelta(i,j)-EU(i)*EU(j));
TEMP(i,j+3) = X12AR*EU(i)*EU(j)+Y12AR*(KroneckerDelta(i,j)-EU(i)*EU(j));
TEMP(i+3,j+3) = TEMP(i,j);
TEMP(i+6,j+6) = X11CR*EU(i)*EU(j)+Y11CR*(KroneckerDelta(i,j)-EU(i)*EU(j));
TEMP(i+6,j+9) = X12CR*EU(i)*EU(j)+Y12CR*(KroneckerDelta(i,j)-EU(i)*EU(j));
TEMP(i+9,j+9) = TEMP(i+6,j+6);
end
end

for i=1:12
for j=1:12
TEMP(j,i) = TEMP(i,j);
end
end

for i=1:12
for j=1:12
TEMP(i,j) = TEMP(i,j)-XM(i,j);
end
end
