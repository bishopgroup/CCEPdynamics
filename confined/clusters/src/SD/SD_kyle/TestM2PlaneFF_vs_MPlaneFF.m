function TestM2PlaneFF_vs_MPlaneFF()
%{
Checks for consistency between M2PlaneFF and MPlaneFF for large wall 
separations H = 1000. 
%}

addpath('./SpherePlaneFF');
addpath('./Sphere2PlaneFF');

%% Load Tables
LoadTablesM2PlaneFF()
LoadTablesMPlaneFF()


%% Prepare Coordinate
H = 1000; % large separation between walls
Xi = (0.001:0.001:0.499)';
h = H*Xi;
N = length(h);


%% Compute Mobility Matrix
M1 = zeros(N,11,11);
M2 = zeros(N,11,11);
D = zeros(N,1);
for i = 1:N
    M1(i,:,:) = MPlaneFF(h(i));
    M2(i,:,:) = M2PlaneFF(Xi(i),H);
    
    D(i) = norm(reshape(M1(i,:,:) - M2(i,:,:) , 11,11));
end


%% Plot Difference
loglog(h,D);
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$|M_{far}^{2B}-M_{far}|$','FontSize',16,'Interpreter','latex');
pause;


% --MUF---------------------------------------------------------------------
%% M(1,1) & M(2,2)
loglog(h,1-M1(:,1,1),h,1-M2(:,1,1),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$1-M_{1,1}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;

%% M(3,3)
loglog(h,1-M1(:,3,3),h,1-M2(:,3,3),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$1-M_{3,3}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause

%--MOL---------------------------------------------------------------------
%% M(4,4) & M(5,5)
loglog(h,0.75-M1(:,4,4),h,0.75-M2(:,4,4),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$3/4 - M_{4,4}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;

%% M(6,6)
loglog(h,0.75-M1(:,6,6),h,0.75-M2(:,6,6),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$3/4 - M_{6,6}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;

%--MUL---------------------------------------------------------------------
%% M(1,5) = -M(2,4)
semilogx(h,M1(:,1,5),h,M2(:,1,5),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$M_{1,5}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','SouthEast');
pause;


%--MES---------------------------------------------------------------------
%% M(7,7)
loglog(h,9/10-M1(:,7,7),h,9/10-M2(:,7,7),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$9/10 - M_{7,7}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;

%% M(8,8)
loglog(h,9/10-M1(:,8,8),h,9/10-M2(:,8,8),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$9/10 - M_{7,7}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;
    
%% M(9,9)
loglog(h,9/10-M1(:,9,9),h,9/10-M2(:,9,9),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$9/10 - M_{9,9}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;

%% M(10,10)
loglog(h,9/10-M1(:,10,10),h,9/10-M2(:,10,10),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$9/10 - M_{10,10}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;

%% M(11,11)
loglog(h,9/10-M1(:,11,11),h,9/10-M2(:,11,11),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$9/10 - M_{11,11}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;

%% M(7,11)
loglog(h,M1(:,7,11),h,M2(:,7,11),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$M_{7,11}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;


%--MUS---------------------------------------------------------------------
%% M(1,9)
loglog(h,-M1(:,1,9),h,-M2(:,1,9),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$-M_{1,9}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;

%% M(2,10)
loglog(h,M1(:,2,10),h,M2(:,2,10),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$M_{2,10}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;

%% M(3,11)
loglog(h,M1(:,3,11),h,M2(:,3,11),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$M_{3,11}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;

%% M(3,7)
semilogx(h,M1(:,3,7),h,M2(:,3,7),'r');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$M_{3,7}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;

%% M(5,9)
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$M_{5,9}$','FontSize',16,'Interpreter','latex');
legend({'$M_{far}^{2B}$','$M_{far}$'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
pause;