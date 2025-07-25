function dydt=ODE_LHS(t,Y,LHSmatrix,x,runs)
%% PARAMETERS %%
Parameter_settings_LHS;
alpha_1 = LHSmatrix(x,1); 
lambda_t = LHSmatrix(x,2); 
rho = LHSmatrix(x,3);
d_1 = LHSmatrix(x,4);
lambda_c=LHSmatrix(x,5);
sigma=LHSmatrix(x,6);
epl=LHSmatrix(x,7);
beta_2=LHSmatrix(x,8);

detla = 10;
d_2 = 0.72;
T_0 = 5.59e6;
d_e = 0.18;
beta_1 = 1.961e-8;
alpha_0 = 2.15e-4;
d_a = 0.47;
K_a = 2e-5;
k = 3e9; 
gamma=5e-6;


C = Y(1);
E = Y(2);
T_1 = Y(3);
T_2 = Y(4);
A = Y(5);
dydt(1) = lambda_c*C.*(1-C/k) - beta_1*E.*C - beta_2*T_1.*C;  %C
dydt(2) = alpha_0*T_0 - d_e*E; %E
dydt(3) = (alpha_1*T_0*C + lambda_t*T_1)/(1+(detla*sigma*T_1.*(T_1+epl*C))*(1-A/(K_a+A)))-d_1*T_1 - rho*T_1;%T_1
dydt(4) = rho*T_1 - d_2*T_2;%T_2
dydt(5) = gamma - d_a*A;%A
dydt= dydt(:);