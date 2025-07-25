clear; clc;

data = xlsread('Data/data_exp_INFalphaAntiPDL1.xlsx');
para = xlsread('Data/data_BestFitPara_INFalphaAntiPDL1.csv');
par = setparameter();
par.Ka = para(1); par.dA = para(2);  par.alphaA = para(3); T1_int = para(4);
par.gamma = 6.25e-6;

options = odeset('RelTol', 1e-6, 'AbsTol', [1e-6,1e-6,1e-6,1e-6,1e-6]);% 定义ODE求解器的选项
t1_span = (0:0.5:40);
y0 =[1e6,5000,T1_int,1e7,0]; %Initial C、S、E、T cells 
[t, R] = ode45(@(t,y) ODE_fit_INFalphaAntiPDL1(t,y,par), t1_span,y0, options);%解微分方程
x2 = data(:,1); y2 = data(:,2);

plot(t1_span,R(:,1),'red-',x2,y2,'black o')
xlabel('Time (days)')
ylabel('Tumor cell volume (mm^3)')
legend('fitted curve','data point'); 

% Save data
dlmwrite('Data/fit_INFalphaAntiPDL1_sim.dat',[t1_span',R(:,1)],' ');
dlmwrite('Data/fit_INFalphaAntiPDL1_exp.dat',[x2,y2],' ');
