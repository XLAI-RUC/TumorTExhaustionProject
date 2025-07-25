clear; clc;

y0 =[1e6,5000,1.375*10^7,10^7,0,0];
t1_span = (0:0.5:60);
options = odeset('RelTol', 1e-6, 'AbsTol', [1e-6,1e-6,1e-6,1e-6,1e-6,1e-6]); 

par = setparameter();

rho_orig = 0.1498;
beta2_orig = 1.5688e-8;
rho_betterT1 = 0.0972;        % 'Better Effector'
beta2_betterT1 = 1.5981e-8;   % 'Better Effector'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Anti-PD-L1
par.gamma = 6.25e-6;
par.zeta = 0;%
[t, R] = ode45(@(t,y) ODE_treatment_PD1IL2v_AntiPDL1(t,y,par,rho_orig,beta2_orig), t1_span,y0, options);%解微分方程
plot(t1_span,R(:,1),'red-')
hold on
xlabel('Time (days)')
ylabel('Tumor cell volume (mm^3)')

dlmwrite('Data/treatment_AntiPDL1_25ug.dat',[t1_span',R(:,1)],' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PD1IL2v
par.zeta = 6.25e-6;%
par.gamma = 0;

t1_span = (0:0.5:60);
[t, R] = ode45(@(t,y) ODE_treatment_PD1IL2v_AntiPDL1(t,y,par,rho_betterT1,beta2_betterT1), t1_span,y0, options);%解微分方程
plot(t1_span,R(:,1),'red-')
xlabel('Time (days)'); ylabel('Tumor cell volume (mm^3)');

dlmwrite('Data/treatment_PD1IL2v_25ug.dat',[t1_span',R(:,1)],' ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PD1IL2v + AntiPDL1
par.zeta = 6.25e-6;    
par.gamma = 6.25e-6;
[t, R] = ode45(@(t,y) ODE_treatment_PD1IL2v_AntiPDL1(t,y,par,rho_betterT1,beta2_betterT1), t1_span,y0, options); 
plot(t1_span,R(:,1),'red-')
xlabel('Time (days)'); ylabel('Tumor cell volume (mm^3)')

dlmwrite('Data/treatment_PD1IL2v_AntiPDL1.dat',[t1_span',R(:,1)],' ');





