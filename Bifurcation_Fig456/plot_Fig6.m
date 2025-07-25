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
File_folder_name = 'output/3D/epl_0p05/'; 
var_Vec = linspace(0.65,0.8,length_r);       % delta
C_matrix = load([File_folder_name, 'Bifur_rho_lambdaT_C_1.dat']);
C_matrix2 = load([File_folder_name, 'Bifur_rho_lambdaT_C_2.dat']);

p = pcolor(rho_Vec,var_Vec,(C_matrix+C_matrix2)/2);
p.FaceColor = 'interp';
p.EdgeColor = 'interp';

xlabel('Exhaustion ($\rho$)','Interpreter','latex');
ylabel('$\lambda_T$','Interpreter','latex');
xlim([0.01 0.15])
text(-0.027,0.805,'\bf{a}','FontSize',FS_lab);
hold on 
text(0.075,0.665,'Tumor-free','FontSize',FS_inside,'Color','w');
text(0.075,0.695,'Bistable','FontSize',FS_inside,'Color','w');
text(0.04,0.75,'Tumor','FontSize',FS_inside);

%%
subplot(2,3,2)
s = surf(rho_Vec, var_Vec,C_matrix);
s.EdgeColor = 'none';
hold on
s = surf(rho_Vec, var_Vec,C_matrix2);
xlabel('$\rho$','Interpreter','latex');
ylabel('$\lambda_T$','Interpreter','latex');
zlabel('Steady state ($C$)','Interpreter','latex');
s.EdgeColor = 'none';
xlim([0.01,0.15]); %ylim([6 12]); zlim([0 1.2e9]);
view(-48, 10);  % 默认设置方位角-37.5°，仰角30°
text(-0.085,0.8,2.91e9,'\bf{b}','FontSize',FS_lab);

%%
subplot(2,3,3)
File_folder_name = 'output/3D/epl_0p05/'; 
var_Vec = linspace(0.5,0.8,length_r);       % delta
C_matrix = load([File_folder_name, 'Bifur_rho_lambdaC_C_1.dat']);
C_matrix2 = load([File_folder_name, 'Bifur_rho_lambdaC_C_2.dat']);

p = pcolor(rho_Vec,var_Vec,(C_matrix+C_matrix2)/2);
p.FaceColor = 'interp';
p.EdgeColor = 'interp';

ylabel('$\lambda_C$','Interpreter','latex');
xlabel('Exhaustion ($\rho$)','Interpreter','latex');
% xlim([0.01 0.12])
text(-0.025,0.804,'\bf{c}','FontSize',FS_lab);
hold on 
text(0.036,0.54,'Tumor-free','FontSize',FS_inside,'Color','w');
text(0.04,0.7,'Bistable','FontSize',FS_inside,'Color','w');
text(0.1,0.73,'Tumor','FontSize',FS_inside);


%%
subplot(2,3,4)
File_folder_name = 'output/3D/epl_0p01/'; 
var_Vec = linspace(0.65,0.8,length_r);       % delta
C_matrix = load([File_folder_name, 'Bifur_rho_lambdaT_C_1.dat']);
C_matrix2 = load([File_folder_name, 'Bifur_rho_lambdaT_C_2.dat']);

p = pcolor(rho_Vec,var_Vec,(C_matrix+C_matrix2)/2);
p.FaceColor = 'interp';
p.EdgeColor = 'interp';

ylabel('$\lambda_T$','Interpreter','latex');
xlabel('Exhaustion ($\rho$)','Interpreter','latex');
% xlim([0.01 0.12])
text(-0.03,0.805,'\bf{d}','FontSize',FS_lab);
hold on 
text(0.06,0.665,'Tumor-free','FontSize',FS_inside,'Color','w');
text(0.04,0.76,'Tumor','FontSize',FS_inside);

%%
subplot(2,3,5)
s = surf(rho_Vec, var_Vec,C_matrix);
s.EdgeColor = 'none';
hold on
s = surf(rho_Vec, var_Vec,C_matrix2);
xlabel('$\rho$','Interpreter','latex');
ylabel('$\lambda_T$','Interpreter','latex');
zlabel('Steady state $(C)$','Interpreter','latex');
s.EdgeColor = 'none';
% xlim([0.01,0.15]); %ylim([6 12]); zlim([0 1.2e9]);
view(-50, 15); %view(30, 14);  % 默认设置方位角-37.5°，仰角30°
text(-0.09,0.8,2.15e9,'\bf{e}','FontSize',FS_lab);

%%
subplot(2,3,6)
File_folder_name = 'output/3D/epl_0p01/'; 
var_Vec = linspace(0.5,0.8,length_r);       % delta
C_matrix = load([File_folder_name, 'Bifur_rho_lambdaC_C_1.dat']);
C_matrix2 = load([File_folder_name, 'Bifur_rho_lambdaC_C_2.dat']);

p = pcolor(rho_Vec,var_Vec,(C_matrix+C_matrix2)/2);
p.FaceColor = 'interp';
p.EdgeColor = 'interp';

ylabel('$\lambda_C$','Interpreter','latex');
xlabel('Exhaustion ($\rho$)','Interpreter','latex');
% xlim([0.01 0.12])
text(-0.025,0.804,'\bf{f}','FontSize',FS_lab);
hold on 
text(0.04,0.58,'Tumor-free','FontSize',FS_inside,'Color','w');
text(0.078,0.75,'Tumor','FontSize',FS_inside);


%% export figure
set(hf,'Units','centimeters');
screenposition = get(gcf,'Position');
set(hf,'PaperPosition',[0 0 1.2*screenposition(3) 0.85*screenposition(4)],...
    'PaperSize',[1.2*screenposition(3) 0.85*screenposition(4)]);

%print(hf,'-r700','-dpdf','Figure/Fig6.pdf');
