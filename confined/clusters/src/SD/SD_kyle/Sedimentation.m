function Sedimentation(H)
%{
Computes the sedimentation velocity of a sphere in a channel.
%}

LoadTablesM2PlaneFF()
LoadTablesMPlaneFF()
LoadTablesRPlaneNF()

%% Positions
Np = 50;
Xi = 1/H + logspace(-5,log10(0.5-1/H),Np)';

%% Compute Far Field Mobility
Uparallel = zeros(Np,1);
Unormal = zeros(Np,1);
Wparallel = zeros(Np,1);
for i = 1:Np
    Rsd = RSD(Xi(i),H);
    Msd = inv(Rsd);

    Uparallel(i) = Msd(1,1); 
    Wparallel(i) = Msd(2,4);
    Unormal(i) = Msd(3,3); 
end

%% Figure 9a
loglog(Xi-1/H,Uparallel);
xlim([1e-5,0.5]);
ylim([0.1,1]);
set(gca,'FontName','Times New Roman');
xlabel('$\Xi-a/H$','FontSize',16,'Interpreter','latex');
ylabel('$U_{\parallel}$','FontSize',16,'Interpreter','latex');

%% Figure 9b
figure;
xi = H*Xi-1;
semilogx(Xi-1/H,Wparallel);
xlim([1e-5,0.5]);
ylim([-0.01,0.05]);
set(gca,'FontName','Times New Roman');
xlabel('$\Xi-a/H$','FontSize',16,'Interpreter','latex');
ylabel('$\Omega_{\parallel}$','FontSize',16,'Interpreter','latex');

%% Figure 10
figure;
loglog(Xi-1/H,Unormal);
xlim([1e-3,0.5]);
ylim([1e-3,1]);
set(gca,'FontName','Times New Roman');
xlabel('$\Xi-a/H$','FontSize',16,'Interpreter','latex');
ylabel('$U_{\perp}$','FontSize',16,'Interpreter','latex');

