function LoadTablesM2PlaneFF(varargin)
%{

Load tabulated functions from
  Swan and Brady, Phys. Fluid, 22, 103301 (2010)

%}

global muf_f1 muf_f3 muf_f5 Iuf_f muf_g1 muf_g3 muf_g5 Iuf_g ...
       mul_f2 mul_f4 Iul_f...
       mus_f2 mus_f4 mus_f6 Ius_f mus_g2 mus_g4 mus_g6 Ius_g ...
       mos_f3 mos_f5 Ios_f ...
       mol_f3 Iol_f mol_g3 Iol_g ...
       mes_f3 mes_f5 mes_f7 Ies_f mes_g3 mes_g5 mes_g7 Ies_g mes_h3 mes_h5 mes_h7 Ies_h

%% Parse Input Arguments
if nargin == 1
    PlotFigures = varargin{1};
else
    PlotFigures = false;
end


%% Unit Vectors & Tensors
d3 = [0 0 1]';
I = eye(3);
eps3jk = [0 1 0; -1 0 0; 0 0 0];

%% UF Coupling
muf_f = load('muf_f.txt');
muf_f1 = griddedInterpolant(muf_f(:,1),-muf_f(:,2));
muf_f3 = griddedInterpolant(muf_f(:,1),muf_f(:,3));
muf_f5 = griddedInterpolant(muf_f(:,1),-muf_f(:,4));
Iuf_f = [1 0 0; 0 1 0; 0 0 0]; % I - d3d3

% Reproduce Figure 2a
if PlotFigures
    Xi = 0.001:0.001:0.999;
    semilogy(Xi,muf_f1(Xi),Xi,muf_f3(Xi),Xi,muf_f5(Xi));
    ylim([1,1e4]);
    set(gca,'FontName','Times New Roman');
    xlabel('$\Xi$','FontSize',16,'Interpreter','latex');
    legend({'$f_{1}^{(UF)}(\Xi)$','$f_{3}^{(UF)}(\Xi)$',...
        '$f_{5}^{(UF)}(\Xi)$'},'FontSize',16,'Interpreter','latex',...
        'location','North');
    pause;
end

muf_g = load('muf_g.txt');
muf_g1 = griddedInterpolant(muf_g(:,1),-muf_g(:,2));
muf_g3 = griddedInterpolant(muf_g(:,1),muf_g(:,3));
muf_g5 = griddedInterpolant(muf_g(:,1),-muf_g(:,4));
Iuf_g = [0 0 0; 0 0 0; 0 0 1]; % d3d3

% Reproduce Figure 2b
if PlotFigures
    Xi = 0.001:0.001:0.999;
    semilogy(Xi,muf_g1(Xi),Xi,muf_g3(Xi),Xi,muf_g5(Xi));
    ylim([1,1e4]);
    set(gca,'FontName','Times New Roman');
    xlabel('$\Xi$','FontSize',16,'Interpreter','latex');
    legend({'$g_{1}^{(UF)}(\Xi)$','$g_{3}^{(UF)}(\Xi)$',...
        '$g_{5}^{(UF)}(\Xi)$'},'FontSize',16,'Interpreter','latex',...
        'location','North');
    pause;
end


%% UL Coupling (OR Equivalently OF)
mul_f = load('mof_f.txt');
mul_f2 = griddedInterpolant(mul_f(:,1),-mul_f(:,2));
mul_f4 = griddedInterpolant(mul_f(:,1),-mul_f(:,3));
Iul_f = eps3jk;

% Reproduce Figure 3a
if PlotFigures
    Xi = 0.001:0.001:0.499;
    semilogy(Xi,mul_f4(Xi));
    ylim([1e-2,1e4]);
    set(gca,'FontName','Times New Roman');
    xlabel('$\Xi$','FontSize',16,'Interpreter','latex');
    ylabel('$f_{4}^{(UL)}(\Xi)$','FontSize',16,'Interpreter','latex');
    pause;
end

% Reproduce Figure 3b
if PlotFigures
    Xi = 0.001:0.001:0.999;
    plot(Xi,mul_f2(Xi));
    ylim([-0.5 0.5]);
    set(gca,'FontName','Times New Roman');
    xlabel('$\Xi$','FontSize',16,'Interpreter','latex');
    ylabel('$f_{2}^{(UL)}(\Xi)$','FontSize',16,'Interpreter','latex');
    pause;
end


%% US Coupling 
mus_f = load('mus_f.txt');
mus_f2 = griddedInterpolant(mus_f(:,1),mus_f(:,2));
mus_f4 = griddedInterpolant(mus_f(:,1),-mus_f(:,3));
mus_f6 = griddedInterpolant(mus_f(:,1),mus_f(:,4));

Ius_f = zeros(3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            Ius_f(i,j,k) = (I(i,j) - d3(i)*d3(j))*d3(k) ...
                + (I(i,k) - d3(i)*d3(k))*d3(j);
        end
    end
end

% Reproduce Figure 4a
if PlotFigures
    Xi = 0.001:0.001:0.499;
    semilogy(Xi,mus_f2(Xi),Xi,mus_f4(Xi),Xi,mus_f6(Xi));
    ylim([1e-1 1e4]);
    set(gca,'FontName','Times New Roman');
    xlabel('$\Xi$','FontSize',16,'Interpreter','latex');
    legend({'$f_{2}^{(US)}(\Xi)$','$f_{4}^{(US)}(\Xi)$',...
        '$f_{6}^{(US)}(\Xi)$'},'FontSize',16,'Interpreter','latex',...
        'location','NorthEast');
    pause;
end

mus_g = load('mus_g.txt');
mus_g2 = griddedInterpolant(mus_f(:,1),-mus_g(:,2));
mus_g4 = griddedInterpolant(mus_f(:,1),mus_g(:,3));
mus_g6 = griddedInterpolant(mus_f(:,1),-mus_g(:,4));

Ius_g = zeros(3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            Ius_g(i,j,k) = d3(i)*(I(j,k) - d3(j)*d3(k)) - 2*d3(i)*d3(j)*d3(k);
        end
    end
end

% Reproduce Figure 4b
if PlotFigures
    Xi = 0.001:0.001:0.499;
    semilogy(Xi,mus_g2(Xi),Xi,mus_g4(Xi),Xi,mus_g6(Xi));
    ylim([1e-1 1e4]);
    set(gca,'FontName','Times New Roman');
    xlabel('$\Xi$','FontSize',16,'Interpreter','latex');
    legend({'$g_{2}^{(US)}(\Xi)$','$g_{4}^{(US)}(\Xi)$',...
        '$g_{6}^{(US)}(\Xi)$'},'FontSize',16,'Interpreter','latex',...
        'location','NorthEast');
    pause;
end


%% OS Coupling 
mos_f = load('mos_f.txt');
mos_f3 = griddedInterpolant(mos_f(:,1),-mos_f(:,2));
mos_f5 = griddedInterpolant(mos_f(:,1),mos_f(:,3));

Ios_f = zeros(3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            Ios_f(i,j,k) = eps3jk(i,j)*d3(k) + eps3jk(i,k)*d3(j);
        end
    end
end

% Reproduce Figure 5
if PlotFigures
    Xi = 0.001:0.001:0.999;
    semilogy(Xi,mos_f3(Xi),Xi,mos_f5(Xi));
    ylim([1 1e4]);
    set(gca,'FontName','Times New Roman');
    xlabel('$\Xi$','FontSize',16,'Interpreter','latex');
    legend({'$f_{3}^{(\Omega S)}(\Xi)$','$f_{5}^{(\Omega S)}(\Xi)$'},...
        'FontSize',16,'Interpreter','latex','location','North');
    pause;
end


%% OL Coupling 
mol_fg = load('mol_fg.txt');
mol_f3 = griddedInterpolant(mol_fg(:,1),mol_fg(:,2));
mol_g3 = griddedInterpolant(mol_fg(:,1),mol_fg(:,3));
Iol_f = Iuf_f;
Iol_g = Iuf_g;

% Reproduce Figure 6
if PlotFigures
    Xi = 0.001:0.001:0.999;
    semilogy(Xi,mol_f3(Xi),Xi,mol_g3(Xi));
    ylim([1 1e4]);
    set(gca,'FontName','Times New Roman');
    xlabel('$\Xi$','FontSize',16,'Interpreter','latex');
    legend({'$f_{3}^{(\Omega L)}(\Xi)$','$g_{3}^{(\Omega L)}(\Xi)$'},...
        'FontSize',16,'Interpreter','latex','location','North');
    pause;
end

%% ES Coupling 
mes_f = load('mes_f.txt');
mes_f3 = griddedInterpolant(mes_f(:,1),-mes_f(:,2));
mes_f5 = griddedInterpolant(mes_f(:,1),mes_f(:,3));
mes_f7 = griddedInterpolant(mes_f(:,1),-mes_f(:,4));

Ies_f = zeros(3,3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                Ies_f(i,j,k,l) = (I(i,j)-d3(i)*d3(j)) * (I(k,l)-d3(k)*d3(l)) ...
                    -2 * (I(i,j)-d3(i)*d3(j)) * d3(k) * d3(l) ...
                    -2 * d3(i) * d3(j) * (I(k,l)-d3(k)*d3(l)) ...
                    +4 * d3(i) * d3(j) * d3(k) * d3(l);
            end
        end
    end
end

% Reproduce Figure 7a
if PlotFigures
    Xi = 0.001:0.001:0.999;
    semilogy(Xi,mes_f3(Xi),Xi,mes_f5(Xi),Xi,mes_f7(Xi));
    ylim([1e-1 1e5]);
    set(gca,'FontName','Times New Roman');
    xlabel('$\Xi$','FontSize',16,'Interpreter','latex');
    legend({'$f_{3}^{(ES)}(\Xi)$','$f_{5}^{(ES)}(\Xi)$',...
        '$f_{7}^{(ES)}(\Xi)$'},'FontSize',16,'Interpreter','latex',...
        'location','North');
    pause;
end

mes_g = load('mes_g.txt');
mes_g3 = griddedInterpolant(mes_g(:,1),-mes_g(:,2));
mes_g5 = griddedInterpolant(mes_g(:,1),mes_g(:,3));
mes_g7 = griddedInterpolant(mes_g(:,1),-mes_g(:,4));

Ies_g = zeros(3,3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                Ies_g(i,j,k,l) = (I(i,k)-d3(i)*d3(k)) * (I(j,l)-d3(j)*d3(l)) ...
                    +(I(i,l)-d3(i)*d3(l)) * (I(j,k)-d3(j)*d3(k)) ...
                    -2 * (I(i,j)-d3(i)*d3(j)) * d3(k) * d3(l) ...
                    -2 * d3(i) * d3(j) * (I(k,l)-d3(k)*d3(l)) ...
                    +4 * d3(i) * d3(j) * d3(k) * d3(l);
            end
        end
    end
end

% Reproduce Figure 7b
if PlotFigures
    Xi = 0.001:0.001:0.999;
    semilogy(Xi,mes_g3(Xi),Xi,mes_g5(Xi),Xi,mes_g7(Xi));
    ylim([1 1e5]);
    set(gca,'FontName','Times New Roman');
    xlabel('$\Xi$','FontSize',16,'Interpreter','latex');
    legend({'$g_{3}^{(ES)}(\Xi)$','$g_{5}^{(ES)}(\Xi)$',...
        '$g_{7}^{(ES)}(\Xi)$'},'FontSize',16,'Interpreter','latex',...
        'location','North');
    pause;
end

mes_h = load('mes_h.txt');
mes_h3 = griddedInterpolant(mes_h(:,1),mes_h(:,2));
mes_h5 = griddedInterpolant(mes_h(:,1),-mes_h(:,3));
mes_h7 = griddedInterpolant(mes_h(:,1),mes_h(:,4));

Ies_h = zeros(3,3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                Ies_h(i,j,k,l) = (I(i,k)-d3(i)*d3(k)) * d3(j) * d3(l) ...
                    + (I(i,l)-d3(i)*d3(l)) * d3(j) * d3(k) ...
                    + (I(j,k)-d3(j)*d3(k)) * d3(i) * d3(l) ...
                    + (I(j,l)-d3(j)*d3(l)) * d3(i) * d3(k);
            end
        end
    end
end

% Reproduce Figure 7c
if PlotFigures
    Xi = 0.001:0.001:0.999;
    semilogy(Xi,mes_h3(Xi),Xi,mes_h5(Xi),Xi,mes_h7(Xi));
    ylim([1 1e5]);
    set(gca,'FontName','Times New Roman');
    xlabel('$\Xi$','FontSize',16,'Interpreter','latex');
    legend({'$h_{3}^{(ES)}(\Xi)$','$h_{5}^{(ES)}(\Xi)$',...
        '$h_{7}^{(ES)}(\Xi)$'},'FontSize',16,'Interpreter','latex',...
        'location','North');
    pause;
end