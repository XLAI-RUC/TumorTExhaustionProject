function plot_Fig4_SS_epsilon_rho()

FS_inside = 9;   % FontSize inside figure
FS_default = 10;  % FontSize default
FS_lab = 14;      % FontSize of (A) (B) (C)

set(0,'defaultaxesfontsize',FS_default,'defaultlinelinewidth',1.5);

epl_Vec = 0.01:0.001:0.15;    % epsilon
rho_Vec = 0.01:0.001:0.15;    % rho

File_folder_name = 'output/3D/'; 
C_matrix = load([File_folder_name, 'Bifur_epl_rho_C.dat']);
T1_matrix = load([File_folder_name, 'Bifur_epl_rho_T1.dat']);
C_matrix2 = load([File_folder_name, 'Bifur_epl_rho_C_2.dat']);
T1_matrix2 = load([File_folder_name, 'Bifur_epl_rho_T1_2.dat']);


p = pcolor(epl_Vec,rho_Vec,(C_matrix+C_matrix2)/2);
p.FaceColor = 'interp';
p.EdgeColor = 'interp';

ylabel('PD-L1 expression ($\varepsilon$)','Interpreter','latex');
xlabel('Exhaustion ($\rho$)','Interpreter','latex');
zlabel('T_1')
xlim([0.01 0.12])

text(0.017,0.03,'Tumor-free','FontSize',FS_inside,'Color','w');
text(0.027,0.12,'Bistable','FontSize',FS_inside);
text(0.083,0.1,'Tumor','FontSize',FS_inside);


end