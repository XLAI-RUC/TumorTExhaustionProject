clear;
clc;

par = setparameter();
par.rho = 0.18;

y0 =[1e6,50,1.3*10^7,10^7,0];  %T、N、Mcells细胞的初值
options = odeset('RelTol', 1e-3, 'AbsTol', [1e-3,1e-3,1e-3,1e-3,1e-3]);% 定义ODE求解器的选项
t_span =(0:0.5:60);

gamma0 = 3.2e-6;
[t, R1] = ode23s(@(t,y) ODE_treatment2(t,y,par,gamma0), t_span,y0, options);%解微分方程
plot(t_span,R1(:,1))
xlabel('Time (days)')
ylabel('Tumor cell number ')
hold on 

gamma0 = 5.2e-6;
[t, R2] = ode23s(@(t,y) ODE_treatment2(t,y,par,gamma0), t_span,y0, options);%解微分方程
plot(t_span,R2(:,1))
xlabel('Time (days)')
ylabel('Tumor cell number ')

gamma0 = 10e-6;
[t, R3] = ode23s(@(t,y) ODE_treatment2(t,y,par,gamma0), t_span,y0, options);%解微分方程
plot(t_span,R3(:,1))
xlabel('Time (days)'); ylabel('Tumor cell number ');

dlmwrite('Data/Tumor_gamma_rho_0p18.dat',[t_span',R1(:,1),R2(:,1),R3(:,1)],' ');

