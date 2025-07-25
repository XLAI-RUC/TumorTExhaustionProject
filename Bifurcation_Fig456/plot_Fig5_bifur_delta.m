function plot_Fig5_bifur_delta()
clc; clear; 

FS_inside = 8;   % FontSize inside figure
FS_default = 10;  % FontSize default

set(0,'defaultaxesfontsize',FS_default,'defaultlinelinewidth',1.3);
    

par =  setparameter();
CK = 3e9;  C_max = 2e9;  rho_max = 0.15;

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

% plot the arrows
ax1 = gca; 
axPos = ax1.Position;  % [left, bottom, width, height]
colorC = [18,106,178]/255;
headW = 3; 

%% plot positive steady state
colorA = [38 122 182]/255;

File_folder_name = 'output/2D/epl_0p05/'; 
RootC1 = load([File_folder_name, 'Bifur_epl_0p05_delta_6_sol1.dat']);  
RootC2 = load([File_folder_name, 'Bifur_epl_0p05_delta_6_sol2.dat']);
rho_Vec = 0.001:0.0005:0.15; 
rho_crit1 = 0.122;
RootC1(rho_Vec<rho_crit1) = nan;
RootC2(rho_Vec<rho_crit1) = nan;
plot(rho_Vec, RootC1(:,1)*CK,'--','Color',colorA); 
hold on 
plot(rho_Vec, RootC2(:,1)*CK,'Color',colorA); 

par.delta = 6; 
sigma_h = par.delta*par.sigma;
alpha0_h = par.alpha0*par.T0;
aa = par.lambdaC - par.beta1*alpha0_h/par.dE;
rho_crit2 = par.lambdaT/(sigma_h*(aa/par.beta2)^2 + 1) - par.d1;
text(0.11,0.5*C_max,'$r_p = 0.6$','FontSize',FS_inside,'Interpreter','latex');
[x1,x2,y1,y2] = normalize_f(ax1,axPos,rho_crit2,rho_crit2,0.1e9,1.25e9);
annotation('arrow',[x1, x2],[y1,y2],'Color','r','HeadWidth',headW,'HeadLength',headW,'LineWidth',0.3)
[x1,x2,y1,y2] = normalize_f(ax1,axPos,rho_crit1-0.004,rho_crit1-0.004,0.42e9,0.05e9);
annotation('arrow',[x1, x2],[y1,y2],'Color','r','HeadWidth',headW,'HeadLength',headW,'LineWidth',0.3)

%%
RootC1 = load([File_folder_name, 'Bifur_epl_0p05_delta_8_sol1.dat']);  
RootC2 = load([File_folder_name, 'Bifur_epl_0p05_delta_8_sol2.dat']);
rho_crit1 = 0.0755;
RootC1(rho_Vec<rho_crit1) = nan;
RootC2(rho_Vec<rho_crit1) = nan;
plot(rho_Vec, RootC1(:,1)*CK,'--','Color',colorA); 
hold on 
plot(rho_Vec, RootC2(:,1)*CK,'Color',colorA); 

par.delta = 8; 
sigma_h = par.delta*par.sigma;
alpha0_h = par.alpha0*par.T0;
aa = par.lambdaC - par.beta1*alpha0_h/par.dE;
rho_crit2 = par.lambdaT/(sigma_h*(aa/par.beta2)^2 + 1) - par.d1;
text(0.06,0.4*C_max,'$r_p = 0.8$','FontSize',FS_inside,'Interpreter','latex');
[x1,x2,y1,y2] = normalize_f(ax1,axPos,rho_crit2,rho_crit2,0.06e9,1.33e9);
annotation('arrow',[x1, x2],[y1,y2],'Color','r','HeadWidth',headW,'HeadLength',headW,'LineWidth',0.3)
[x1,x2,y1,y2] = normalize_f(ax1,axPos,rho_crit1-0.003,rho_crit1-0.003,0.49e9,0.05e9);
annotation('arrow',[x1, x2],[y1,y2],'Color','r','HeadWidth',headW,'HeadLength',headW,'LineWidth',0.3)


%%
RootC1 = load([File_folder_name, 'Bifur_epl_0p05_delta_10_sol1.dat']);  
RootC2 = load([File_folder_name, 'Bifur_epl_0p05_delta_10_sol2.dat']);
rho_crit1 = 0.0365;
RootC1(rho_Vec<rho_crit1) = nan;
RootC2(rho_Vec<rho_crit1) = nan;
plot(rho_Vec, RootC1(:,1)*CK,'b--'); 
hold on 
plot(rho_Vec, RootC2(:,1)*CK,'b');

par.delta = 10; 
sigma_h = par.delta*par.sigma;
alpha0_h = par.alpha0*par.T0;
aa = par.lambdaC - par.beta1*alpha0_h/par.dE;
rho_crit2 = par.lambdaT/(sigma_h*(aa/par.beta2)^2 + 1) - par.d1;
text(0.02,0.4*C_max,'$r_p = 1$','FontSize',FS_inside,'Interpreter','latex');
[x1,x2,y1,y2] = normalize_f(ax1,axPos,rho_crit2,rho_crit2,0.1e9,2.6e9);
annotation('arrow',[x1, x2],[y1,y2],'Color','r','HeadWidth',headW,'HeadLength',headW,'LineWidth',0.2)
[x1,x2,y1,y2] = normalize_f(ax1,axPos,rho_crit1-0.003,rho_crit1-0.003,0.9e9,0.08e9);
annotation('arrow',[x1, x2],[y1,y2],'Color','r','HeadWidth',headW,'HeadLength',headW,'LineWidth',0.2)


%%
RootC1 = load([File_folder_name, 'Bifur_epl_0p05_delta_12_sol1.dat']);  
RootC2 = load([File_folder_name, 'Bifur_epl_0p05_delta_12_sol2.dat']);
rho_crit1 = 0.0035;
RootC1(rho_Vec<rho_crit1) = nan;
RootC2(rho_Vec<rho_crit1) = nan;
plot(rho_Vec, RootC1(:,1)*CK,'b--'); 
hold on 
plot(rho_Vec, RootC2(:,1)*CK,'b'); 

par.delta = 12; 
sigma_h = par.delta*par.sigma;
alpha0_h = par.alpha0*par.T0;
aa = par.lambdaC - par.beta1*alpha0_h/par.dE;
rho_crit2 = par.lambdaT/(sigma_h*(aa/par.beta2)^2 + 1) - par.d1;
text(0.001,0.6*C_max,'$r_p = 1.2$','FontSize',FS_inside,'Interpreter','latex');
[x1,x2,y1,y2] = normalize_f(ax1,axPos,rho_crit2,rho_crit2,0.1e9,2.6e9);
annotation('arrow',[x1, x2],[y1,y2],'Color','r','HeadWidth',headW,'HeadLength',headW,'LineWidth',0.2)
[x1,x2,y1,y2] = normalize_f(ax1,axPos,rho_crit1-0.002,rho_crit1-0.002,1e9,0.08e9);
annotation('arrow',[x1, x2],[y1,y2],'Color','r','HeadWidth',headW,'HeadLength',headW,'LineWidth',0.2)


ylim([0 C_max]); xlim([0 rho_max]);

xlabel('Exhaustion ($\rho$)','Interpreter','latex');
ylabel('Steady state ($C$)','Interpreter','latex');

% text(0.0415,0.5*C_max,'$\varepsilon = 0.05$','Interpreter','latex','FontSize',FS_inside);
% text(0.022,0.54*C_max,'$\varepsilon = 0.06$','Interpreter','latex','FontSize',FS_inside);
% text(0.061,0.42*C_max,'$\varepsilon = 0.04$','Interpreter','latex','FontSize',FS_inside);
% text(0.078,0.35*C_max,'$\varepsilon = 0.03$','Interpreter','latex','FontSize',FS_inside);
% 
% text(0.1*0.1,0.95*C_max,'Tumor-free','FontSize',FS_inside);
% text(0.42*0.1,0.95*C_max,'Bistable','FontSize',FS_inside);
% text(0.61*0.1,0.95*C_max,'Tumor','FontSize',FS_inside);

end