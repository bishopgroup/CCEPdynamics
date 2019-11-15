function TestRPlaneNF_vs_MPlaneFF()
%{
Checks for consistency between RPlaneNF and MPlaneFF.
%}

addpath('./SpherePlaneFF');
addpath('./SpherePlaneNF');

%% Load Tables
LoadTablesRPlaneNF()
LoadTablesMPlaneFF()


%% Prepare Coordinate
Min = -3;
Max = 3;
N = 1000;
xi = logspace(Min,Max,N)';
h = 1+xi;

 
%% Compute Resistance Matrix
RFF = zeros(N,6,6);
RFFS = zeros(N,6,6);
R = zeros(N,6,6);
D = zeros(N,1);
for i = 1:N
    % Compute far-field mobility invert
    MFF = MPlaneFF(h(i));
    tmp = inv(MFF);
    RFFS(i,:,:) = tmp(1:6,1:6); % include stresslet

    tmp = inv(MFF(1:6,1:6));
    RFF(i,:,:) = tmp; % without stresslet contribution
        
    % Compute exact sphere-plane resistance matrix
    R(i,:,:) = RPlaneNF(h(i));
    
    % Error
    D(i) = norm(reshape(R(i,:,:) - RFFS(i,:,:) , 6,6));
end

%% Plot Difference
loglog(xi,abs(D),xi,1./xi.^6,'k:');
set(gca,'FontName','Times New Roman');
xlabel('$h$','FontSize',16,'Interpreter','latex');
ylabel('$|R_{\infty }^{2B} - R^{2B}|$','FontSize',16,'Interpreter','latex');
ylim([1e-10,10]);
pause;

%--RFU---------------------------------------------------------------------
%% R(1,1) = R(2,2)
loglog(xi,abs(RFF(:,1,1)-R(:,1,1)),'b',xi,abs(RFFS(:,1,1)-R(:,1,1)),'r',xi,1./xi.^6,'k:');
set(gca,'FontName','Times New Roman');
xlabel('$\xi$','FontSize',16,'Interpreter','latex');
ylabel('$|R_{\infty }^{2B} - R^{2B}|$','FontSize',16,'Interpreter','latex');
title('$R_{FU}^{11}$','FontSize',16,'Interpreter','latex')
legend({'FL','FLS'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
ylim([1e-10,10]);
pause;


%% R(3,3)
loglog(xi,abs(RFF(:,3,3)-R(:,3,3)),'b',xi,abs(RFFS(:,3,3)-R(:,3,3)),'r',xi,1./xi.^6,'k:');
set(gca,'FontName','Times New Roman');
xlabel('$\xi$','FontSize',16,'Interpreter','latex');
ylabel('$|R_{\infty }^{2B} - R^{2B}|$','FontSize',16,'Interpreter','latex');
title('$R_{FU}^{33}$','FontSize',16,'Interpreter','latex')
legend({'FL','FLS'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
ylim([1e-10,10]);
pause;

%--RLU---------------------------------------------------------------------
%% R(1,5) = -R(2,4)
loglog(xi,abs(RFF(:,1,5)-R(:,1,5)),'b',xi,abs(RFFS(:,1,5)-R(:,1,5)),'r',xi,0.1./xi.^6,'k:');
set(gca,'FontName','Times New Roman');
xlabel('$\xi$','FontSize',16,'Interpreter','latex');
ylabel('$|R_{\infty }^{2B} - R^{2B}|$','FontSize',16,'Interpreter','latex');
title('$R_{LU}^{12}$','FontSize',16,'Interpreter','latex')
legend({'FL','FLS'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
ylim([1e-10,10]);
pause;

 
%--RLO---------------------------------------------------------------------
%% R(4,4) = R(5,5)
loglog(xi,abs(RFF(:,4,4)-R(:,4,4)),'b',xi,abs(RFFS(:,4,4)-R(:,4,4)),'r',xi,0.1./xi.^7,'k:');
set(gca,'FontName','Times New Roman');
xlabel('$\xi$','FontSize',16,'Interpreter','latex');
ylabel('$|R_{\infty }^{2B} - R^{2B}|$','FontSize',16,'Interpreter','latex');
title('$R_{L \Omega}^{11}$','FontSize',16,'Interpreter','latex')
legend({'FL','FLS'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
ylim([1e-10,1]);
pause;

%% R(6,6)
loglog(xi,abs(RFF(:,6,6)-R(:,6,6)),'b',xi,abs(RFFS(:,6,6)-R(:,6,6)),'r',xi,0.1./xi.^7,'k:');
set(gca,'FontName','Times New Roman');
xlabel('$\xi$','FontSize',16,'Interpreter','latex');
ylabel('$|R_{\infty }^{2B} - R^{2B}|$','FontSize',16,'Interpreter','latex');
title('$R_{L \Omega}^{33}$','FontSize',16,'Interpreter','latex')
legend({'FL','FLS'},'FontSize',16,'Interpreter','latex',...
    'location','NorthEast');
ylim([1e-10,10]);
pause;

