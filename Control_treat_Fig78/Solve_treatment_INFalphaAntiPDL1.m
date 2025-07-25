clear; clc;


par = setparameter();
par.gamma = 6.25e-6; % paramaters for INFalphaAntiPDL1 kinetics
par.Ka = 9.5e-07;
par.dA = 0.872;
par.alphaA = 3.29e-7;

options = odeset('RelTol', 1e-6, 'AbsTol', [1e-6,1e-6,1e-6,1e-6,1e-6]);
t1_span=(0:0.5:60);
y0 =[1e6,5000,1.375*1e7,1e7,0];  %Initial C¡¢S¡¢E¡¢T cells 

[t, R] = ode45(@(t,y) ODE_treatment_INFalphaAntiPDL1(t,y,par), t1_span,y0, options);
plot(t1_span,R(:,1),'red-')
xlabel('Time (days)')
ylabel('Tumor cell volume (mm^3)')
% save data
dlmwrite('Data/treatment_INFalphaAntiPDL1_25ug.dat',[t1_span',R(:,1)],' ');
