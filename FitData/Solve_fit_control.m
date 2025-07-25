clear; clc;
data = xlsread('Data/data_experiment_Liang.xlsx','ctrl');
para = xlsread('Data/data_BestFitPara_control.csv');
par = setparameter();
par.lambdaC = para(1);  par.alpha1 = para(2);  par.lambdaT = para(3);  par.rho = para(4);

y0 = [1e6,5000,1.375*10^7,10^7];  % Initial values 
options = odeset('RelTol', 1e-3, 'AbsTol', [1e-3,1e-3,1e-3,1e-3]);% 定义ODE求解器的选项
t1_span = (0:0.5:60);
[t, R] = ode45(@(t,y) ODE_control(t,y,par), t1_span,y0, options);%解微分方程
C = R(:,1);  E = R(:,2);  T_1 = R(:,3);  T_2 = R(:,4);
x2 = data(:,1);  y2 = data(:,2);
plot(t1_span,C,'red-',x2,y2,'black o')
xlabel('Time (days)')
ylabel('Tumor cell volume (mm^3)')
legend('fitted curve','data point'); 
save('data1','t1_span','T_1','C','x2','y2');

dlmwrite('Data/Control_simulation.dat',[t1_span',R(:,1),R(:,3)],' ');
dlmwrite('Data/Control_experiment.dat',[x2,y2],' ');