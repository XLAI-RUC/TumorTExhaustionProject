clear; clc;

data = xlsread('Data/data_exp_PD1_IL2v.xlsx');

par = setparameter();
para = xlsread('Data/data_BestFitPara_PD1IL2V.csv');
par.Kb = para(1);
par.beta2 = para(2);
par.zeta = para(3);
par.rho = para(4);
par.dB = para(5);

options = odeset('RelTol', 1e-6, 'AbsTol', [1e-6,1e-6,1e-6,1e-6,1e-6]);% 定义ODE求解器的选项
t1_span = (0:0.5:30);
y0 = [7.8e7,5000,4*10^7,10^8,0]; 

[t, R] = ode45(@(t,y) ODE_fit_PD1IL2v(t,y,par), t1_span,y0, options);
x2 = data(:,1); y2 = data(:,2);
plot(t1_span,R(:,1),'red-',x2,y2,'black o-')
xlabel('Time (days)'); ylabel('Tumor cell volume (mm^3)')
legend('fitted curve','data point'); 

% Save data
dlmwrite('Data/fit_PD1IL2v_sim.dat',[t1_span',R(:,1)],' ');
dlmwrite('Data/fit_PD1IL2v_exp.dat',[x2,y2],' ');