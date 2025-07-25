function plot_Fig5_bifur_delta_2()
clc; 

FS_inside = 8;   % FontSize inside figure
FS_default = 10;  % FontSize default

set(0,'defaultaxesfontsize',FS_default,'defaultlinelinewidth',1.3);
    

par =  setparameter();
CK = 3e9;  C_max = 1e9;  rho_max = 0.15;

sigma_h = par.delta*par.sigma;
alpha0_h = par.alpha0*par.T0;
aa = par.lambdaC - par.beta1*alpha0_h/par.dE
rho_crit2 = par.lambdaT/(sigma_h*(aa/par.beta2)^2 + 1) - par.d1;
rho_crit1 = 0.0365;

%% plot filled area
X1 = [0, rho_crit2,rho_crit2,0];
Y1 = [0, 0, C_max, C_max];
colorA = [229,255,204]/255; 
plot(X1,Y1,'LineWidth',1,'color',colorA)
m = fill(X1,Y1,'g');
m.FaceColor = colorA;
m.LineWidth = 0.7;
m.EdgeColor = colorA;
hold on 

X1 = [rho_crit2,rho_max,rho_max,rho_crit2];
Y1 = [0, 0,C_max, C_max];
colorB = [1 0.92 0.8];
plot(X1,Y1,'LineWidth',1,'color',colorB )
m = fill(X1,Y1,'g');
m.FaceColor = colorB ;
m.LineWidth = 0.7;
m.EdgeColor = colorB ;

%% plot positive steady state
colorA = [38 122 182]/255;
rho_Vec = 0.001:0.0005:0.15; 

%%
File_folder_name = 'output/2D/epl_0p01/'; 
RootC1 = load([File_folder_name, 'Bifur_epl_0p01_delta_10_sol1.dat']);  
RootC2 = load([File_folder_name, 'Bifur_epl_0p01_delta_10_sol2.dat']);
rho_crit1 = 0.0365;
RootC1(rho_Vec<rho_crit1) = nan;
RootC2(rho_Vec<rho_crit1) = nan;
plot(rho_Vec, RootC1(:,1)*CK,'b'); 
hold on 
text(0.085,0.5*C_max,'$r_p = 1$','FontSize',FS_inside,'Interpreter','latex');


%%
RootC1 = load([File_folder_name, 'Bifur_epl_0p01_delta_8_sol1.dat']);  
RootC2 = load([File_folder_name, 'Bifur_epl_0p01_delta_8_sol2.dat']);
rho_crit1 = 0.0365;
RootC1(rho_Vec<rho_crit1) = nan;
RootC2(rho_Vec<rho_crit1) = nan;
plot(rho_Vec, RootC1(:,1)*CK,'Color',colorA); 
hold on 
text(0.11,0.36*C_max,'$r_p = 0.8$','FontSize',FS_inside,'Interpreter','latex');


%%
RootC1 = load([File_folder_name, 'Bifur_epl_0p01_delta_12_sol1.dat']);  
RootC2 = load([File_folder_name, 'Bifur_epl_0p01_delta_12_sol2.dat']);
plot(rho_Vec, RootC1(:,1)*CK,'Color',colorA); 
hold on 
text(0.07,0.65*C_max,'$r_p = 1.2$','FontSize',FS_inside,'Interpreter','latex');

% File_folder_name = 'output/2D/epl_0p05/'; 
% RootC1 = load([File_folder_name, 'Bifur_epl_0p05_delta_6_sol1.dat']);  
% RootC2 = load([File_folder_name, 'Bifur_epl_0p05_delta_6_sol2.dat']);

% rho_crit1 = 0.122;
% RootC1(rho_Vec<rho_crit1) = nan;
% RootC2(rho_Vec<rho_crit1) = nan;
% plot(rho_Vec, RootC1(:,1)*CK,'--','Color',colorA); 
% hold on 
% plot(rho_Vec, RootC2(:,1)*CK,'Color',colorA); 

% %%
% File_folder_name = 'output/2D/'; 
% RootC1 = load([File_folder_name, 'Bifur_epl_0p05_delta_8_sol1.dat']);  
% RootC2 = load([File_folder_name, 'Bifur_epl_0p05_delta_8_sol2.dat']);
% rho_crit1 = 0.0755;
% RootC1(rho_Vec<rho_crit1) = nan;
% RootC2(rho_Vec<rho_crit1) = nan;
% plot(rho_Vec, RootC1(:,1)*CK,'--','Color',colorA); 
% hold on 
% plot(rho_Vec, RootC2(:,1)*CK,'Color',colorA); 

% %%
% File_folder_name = 'output/2D/'; 
% RootC1 = load([File_folder_name, 'Bifur_epl_0p05_delta_12_sol1.dat']);  
% RootC2 = load([File_folder_name, 'Bifur_epl_0p05_delta_12_sol2.dat']);
% rho_crit1 = 0.0035;
% RootC1(rho_Vec<rho_crit1) = nan;
% RootC2(rho_Vec<rho_crit1) = nan;
% plot(rho_Vec, RootC1(:,1)*CK,'b--'); 
% hold on 
% plot(rho_Vec, RootC2(:,1)*CK,'b'); 


ylim([0 C_max]); xlim([0 rho_max]);

xlabel('Exhaustion ($\rho$)','Interpreter','latex');
ylabel('Steady state ($C$)','Interpreter','latex');


end