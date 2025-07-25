clc; clear; close all

FS_inside = 8;   % FontSize inside figure
FS_default = 9;  % FontSize default
FS_lab = 10;      % FontSize of (A) (B) (C)

set(0,'defaultaxesfontsize',FS_default,'defaultlinelinewidth',1.5);

rho_Vec = 0.01:0.001:0.15;    % rho
length_r = length(rho_Vec);


hf = figure;

%%
subplot(2,3,1)
plot_Fig5_bifur_delta()
text(-0.035,2.04e9,'\bf{a}','FontSize',FS_lab);

%%
subplot(2,3,2)
File_folder_name = 'output/3D/epl_0p05/'; 
var_Vec = linspace(6,12,length_r)/10;        % delta
C_matrix = load([File_folder_name, 'Bifur_rho_delta_C_1.dat']);
C_matrix2 = load([File_folder_name, 'Bifur_rho_delta_C_2.dat']);

p = pcolor(rho_Vec,var_Vec,(C_matrix+C_matrix2)/2);
p.FaceColor = 'interp';
p.EdgeColor = 'interp';
xlabel('Exhaustion ($\rho$)','Interpreter','latex');
ylabel('PD-1 expression ($r_p$)','Interpreter','latex');
xlim([0.01 0.12])
text(-0.015,1.202,'\bf{b}','FontSize',FS_lab);
text(0.02,0.7,'Tumor-free','FontSize',FS_inside,'Color','w');
text(0.038,1,'Bistable','FontSize',FS_inside,'Color','w');
text(0.08,1.1,'Tumor','FontSize',FS_inside);

sf = subplot(2,3,3);
s = surf(rho_Vec, var_Vec,C_matrix);
s.EdgeColor = 'none';
hold on
s = surf(rho_Vec, var_Vec,C_matrix2);
xlabel('$\rho$','Interpreter','latex');
ylabel('$r_p$','Interpreter','latex');
zlabel('Steady state ($C$)','Interpreter','latex');
s.EdgeColor = 'none';
xlim([0.01,0.15]);  ylim([6 12]/10)
view(37, 14);  % 默认设置方位角-37.5°，仰角30° 
text(0.04,0,3.2e9,'\bf{c}','FontSize',FS_lab);

%%
subplot(2,3,4)
plot_Fig5_bifur_delta_2()
text(-0.035,1e9,'\bf{d}','FontSize',FS_lab);

%%
subplot(2,3,5)
File_folder_name = 'output/3D/epl_0p01/'; 
var_Vec = linspace(6,12,length_r)/10;        % delta
C_matrix = load([File_folder_name, 'Bifur_rho_delta_C_1.dat']);
C_matrix2 = load([File_folder_name, 'Bifur_rho_delta_C_2.dat']);

p = pcolor(rho_Vec,var_Vec,(C_matrix+C_matrix2)/2);
p.FaceColor = 'interp';
p.EdgeColor = 'interp';
xlabel('Exhaustion ($\rho$)','Interpreter','latex');
ylabel('PD-1 expression ($r_p$)','Interpreter','latex');
xlim([0.01 0.12]);
text(-0.015,12/10,'\bf{e}','FontSize',FS_lab);
hold on 
text(0.02,7/10,'Tumor-free','FontSize',FS_inside,'Color','w');
text(0.08,11/10,'Tumor','FontSize',FS_inside);

%%
subplot(2,3,6)
s = surf(rho_Vec, var_Vec,C_matrix);
s.EdgeColor = 'none';
hold on
s = surf(rho_Vec, var_Vec,C_matrix2);
xlabel('$\rho$','Interpreter','latex');
ylabel('$r_p$','Interpreter','latex');
zlabel('Steady state ($C$)','Interpreter','latex');
s.EdgeColor = 'none';
xlim([0.01,0.15]); ylim([6 12]/10); zlim([0 1.2e9]);
view(30, 14);  % 默认设置方位角-37.5°，仰角30°
text(0.025,0,1.7e9,'\bf{f}','FontSize',FS_lab);

%% export figure
set(hf,'Units','centimeters');
screenposition = get(gcf,'Position');
set(hf,'PaperPosition',[0 0 1.2*screenposition(3) 0.9*screenposition(4)],...
    'PaperSize',[1.2*screenposition(3) 0.9*screenposition(4)]);


% print(hf,'-r700','-dpdf','Figure/Fig5.pdf');
