function plot_Fig4_bifur_epl()
clc; clear; 

FS_inside = 8;   % FontSize inside figure
FS_default = 10;  % FontSize default
FS_lab = 14;      % FontSize of (A) (B) (C)
LineW = 1.3;

set(0,'defaultaxesfontsize',FS_default,'defaultlinelinewidth',1);
  

par =  setparameter();
CK = 3e9;  C_max = 2e9;  rho_max = 0.1;

sigma_h = par.delta*par.sigma;
alpha0_h = par.alpha0*par.T0;
aa = par.lambdaC - par.beta1*alpha0_h/par.dE
rho_crit2 = par.lambdaT/(sigma_h*(aa/par.beta2)^2 + 1) - par.d1;
rho_crit1 = 0.0365;

%% plot filled area
X1 = [0, rho_crit1,rho_crit1,0];
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

File_folder_name = 'output/2D/'; 
RootC1 = load([File_folder_name, 'Bifur_epl_0p05_sol1.dat']);
RootC2 = load([File_folder_name, 'Bifur_epl_0p05_sol2.dat']);
rho_Vec = 0.01:0.0005:0.15; 
rho_crit1 = 0.0365;

RootC1(rho_Vec<rho_crit1) = nan;
RootC2(rho_Vec<rho_crit1) = nan;

plot(rho_Vec, RootC1(:,1)*CK,'b--','LineWidth',LineW); 
hold on 
plot(rho_Vec, RootC2(:,1)*CK,'b','LineWidth',LineW); 


%%
File_folder_name = 'output/2D/'; 
RootC1 = load([File_folder_name, 'Bifur_epl_0p06_sol1.dat']);
RootC2 = load([File_folder_name, 'Bifur_epl_0p06_sol2.dat']);
rho_crit1 = 0.018;

RootC1(rho_Vec<rho_crit1) = nan;
RootC2(rho_Vec<rho_crit1) = nan;

colorA = [38 122 182]/255;
plot(rho_Vec, RootC1(:,1)*CK,'--','LineWidth',LineW,'Color',colorA); 
hold on 
plot(rho_Vec, RootC2(:,1)*CK,'LineWidth',LineW,'Color',colorA); 
ylim([0 C_max]); xlim([0 rho_max]);

%%
File_folder_name = 'output/2D/'; 
RootC1 = load([File_folder_name, 'Bifur_epl_0p04_sol1.dat']);
RootC2 = load([File_folder_name, 'Bifur_epl_0p04_sol2.dat']);
rho_crit1 = 0.053;

RootC1(rho_Vec<rho_crit1) = nan;
RootC2(rho_Vec<rho_crit1) = nan;

plot(rho_Vec, RootC1(:,1)*CK,'--','LineWidth',LineW,'Color',colorA); 
hold on 
plot(rho_Vec, RootC2(:,1)*CK,'LineWidth',LineW,'Color',colorA); 
ylim([0 C_max]); xlim([0 rho_max]);

%%
File_folder_name = 'output/2D/'; 
RootC1 = load([File_folder_name, 'Bifur_epl_0p03_sol1.dat']);
RootC2 = load([File_folder_name, 'Bifur_epl_0p03_sol2.dat']);
rho_crit1 = 0.053;

RootC1(rho_Vec<rho_crit1) = nan;
RootC2(rho_Vec<rho_crit1) = nan;

plot(rho_Vec, RootC1(:,1)*CK,'LineWidth',LineW,'Color',colorA); 

%% plot steady state E_02
plot([0 rho_crit2], [0 0],'r','LineWidth',2); 
plot([rho_crit2 rho_Vec(end)], [0 0],'r--','LineWidth',2); 


ylim([0 C_max]); xlim([0 rho_max]);

xlabel('Exhaustion ($\rho$)','Interpreter','latex');
ylabel('Steady state ($C$)','Interpreter','latex');

text(0.036,0.6*C_max,'$\varepsilon = 0.05$','Interpreter','latex','FontSize',FS_inside);
text(0.02,0.68*C_max,'$\varepsilon = 0.06$','Interpreter','latex','FontSize',FS_inside);
text(0.053,0.5*C_max,'$\varepsilon = 0.04$','Interpreter','latex','FontSize',FS_inside);
text(0.07,0.4*C_max,'$\varepsilon = 0.03$','Interpreter','latex','FontSize',FS_inside);

text(0.01*0.1,0.92*C_max,'Tumor-free','FontSize',FS_inside);
text(0.37*0.1,0.92*C_max,'Bistable','FontSize',FS_inside);
text(0.65*0.1,0.92*C_max,'Tumor','FontSize',FS_inside);

end