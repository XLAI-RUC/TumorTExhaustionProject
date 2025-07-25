function plot_Fig4_bifur_T1()
clc; clear; 

FS_inside = 8;   % FontSize inside figure
FS_default = 10;  % FontSize default
FS_lab = 14;      % FontSize of (A) (B) (C)
LineW = 1.5;

set(0,'defaultaxesfontsize',FS_default,'defaultlinelinewidth',1);
    
par =  setparameter();
CK = 3e9;  C_max = 6e7;  rho_max = 0.1;

File_folder_name = 'output/2D/'; 
RootC1 = load([File_folder_name, 'Bifur_epl_0p05_sol1.dat']);
RootC2 = load([File_folder_name, 'Bifur_epl_0p05_sol2.dat']);
rho_Vec = 0.01:0.0005:0.15; 
rho_crit1 = 0.0365;

sigma_h = par.delta*par.sigma;
alpha0_h = par.alpha0*par.T0;
aa = par.lambdaC - par.beta1*alpha0_h/par.dE
rho_crit2 = par.lambdaT/(sigma_h*(aa/par.beta2)^2 + 1) - par.d1;


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

%% plot steady state E_02
rho_Vec0 = 0:0.0005:0.1; 
T1_eq = sqrt((par.lambdaT - rho_Vec0 - par.d1)./(sigma_h*(rho_Vec0+ par.d1)));
rho_Vec1 = rho_Vec0(rho_Vec0<rho_crit2);
T1_eq1 = T1_eq(rho_Vec0<rho_crit2);
rho_Vec2 = rho_Vec0(rho_Vec0>rho_crit2);
T1_eq2 = T1_eq(rho_Vec0>rho_crit2);

hold on 
plot(rho_Vec1, T1_eq1,'r','LineWidth',LineW); 
plot(rho_Vec2, T1_eq2,'r--','LineWidth',LineW); 

%% plot steady state E_12 and E_22
RootT1 = RootC1(:,2); RootT2 = RootC2(:,2); 
RootT1(rho_Vec<rho_crit1) = nan;
RootT2(rho_Vec<rho_crit1) = nan;

plot(rho_Vec, RootT1*CK,'b--','LineWidth',LineW); 
hold on 
plot(rho_Vec, RootT2*CK,'b','LineWidth',LineW); 
ylim([0 C_max]); xlim([0 rho_max]);

% plot the arrows
ax1 = gca; 
axPos = ax1.Position;  % [left, bottom, width, height]
colorC = [18,106,178]/255;
headW = 5; 
[x1,x2,y1,y2] = normalize_f(ax1,axPos,0.033,0.033,3.6e7,4.6e7);
annotation('arrow',[x1, x2],[y1,y2],'Color',colorC,'HeadWidth',headW,'HeadLength',headW)
[x1,x2,y1,y2] = normalize_f(ax1,axPos,0.063,0.063,4.1e7,2.8e7);
annotation('arrow',[x1, x2],[y1,y2],'Color',colorC,'HeadWidth',headW,'HeadLength',headW)
[x1,x2,y1,y2] = normalize_f(ax1,axPos,0.056,0.04,2e7,2.6e7);
annotation('arrow',[x1, x2],[y1,y2],'Color',colorC,'HeadWidth',headW,'HeadLength',headW)
[x1,x2,y1,y2] = normalize_f(ax1,axPos,0.04,0.056,5.3e7,5e7);
annotation('arrow',[x1, x2],[y1,y2],'Color',colorC,'HeadWidth',headW,'HeadLength',headW)


xlabel('Exhaustion ($\rho$)','Interpreter','latex');
ylabel('Steady state ($T_1$)','Interpreter','latex');

text(0.02*0.1,0.1*C_max,'Tumor-free','FontSize',FS_inside);
text(0.37*0.1,0.1*C_max,'Bistable','FontSize',FS_inside);
text(0.65*0.1,0.1*C_max,'Tumor','FontSize',FS_inside);


end