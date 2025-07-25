%% The results should be compared to the PRCC results section in
%% Supplementary Material D and Table D.1 for different N (specified by
%% "runs" in the script below
clear;
clc;

%% Sample size N
runs = 20000;

%% LHS MATRIX  %%
Parameter_settings_LHS_INFalphaAntiPDL1_25ug;

par = setparameter();

par.gamma = 6.25e-6; % paramaters for INFalphaAntiPDL1 kinetics
par.Ka = 9.5e-07;
par.dA = 0.872;
par.alphaA = 3.29e-7;

a = 0.5;
b = 2;

alpha_1_LHS  = LHS_Call(a*par.alpha1, par.alpha1, 5*par.alpha1,  0 ,runs,'unif');   % baseline = 1.4637e-6
lambda_t_LHS = LHS_Call(a*par.lambdaT, par.lambdaT, b*par.lambdaT, 0 ,runs,'unif');  % baseline = 0.7369
rho_LHS      = LHS_Call(a*par.rho, par.rho,   b*par.rho,   0 ,runs,'unif');      % baseline = 0.1498
d_1_LHS      = LHS_Call(a*par.d1, par.d1,        b*par.d1,     0 ,runs,'unif');  % baseline = 0.41
lambda_c_LHS = LHS_Call(a*par.lambdaC, par.lambdaC, b*par.lambdaC,  0 ,runs,'unif'); % baseline = 0.7163
sigma_LHS    = LHS_Call(a*par.sigma, par.sigma, b*par.sigma,  0 ,runs,'unif');     % baseline = 2.725e-17;
epl_LHS      = LHS_Call(a*par.epl, par.epl, 15*par.epl,  0 ,runs,'unif');          % baseline = 1.4637e-6;
beta_2_LHS   = LHS_Call(a*par.beta2, par.beta2, b*par.beta2,  0 ,runs,'unif');    % baseline = 1.5688e-8;

%% LHS MATRIX and PARAMETER LABELS
LHSmatrix = [alpha_1_LHS lambda_t_LHS rho_LHS d_1_LHS lambda_c_LHS sigma_LHS epl_LHS beta_2_LHS] ;

Time_index = [time_points1,time_points2,time_points3,time_points4];

for x=1:runs %Run solution x times choosing different values
    f = @ODE_LHS;
    options = odeset('RelTol', 1e-6, 'AbsTol', [1e-6,1e-6,1e-6,1e-6,1e-6]);% 定义ODE求解器的选项
    [t, y] = ode45(@(t,y) ODE_LHS_INFalphaAntiPDL1(t,y,LHSmatrix(x,:),par), tspan,y0, options);%解微分方程
    B = [t y]; % [time y]
  
    %% Save only the outputs at the time points of interest [time_points]:
    C_1(:,x) = B(Time_index+1,2);
    E_1(:,x) = B(Time_index+1,3);
    T1_1(:,x) = B(Time_index+1,4);
    T2_1(:,x) = B(Time_index+1,5);
    A_1(:,x) = B(Time_index+1,6);
end
%% Save the workspace
% CALCULATE PRCC
alpha = 0.05;
for j = 1:length(Time_index)
    [prccC(j,:), pvalueC(j,:), sign_labelC(j,:)] = PRCC(LHSmatrix,C_1(j,:),1:length(Time_index(j)),PRCC_var,alpha);
    [prccT1(j,:), pvalueT1(j,:), sign_labelT1(j,:)] = PRCC(LHSmatrix,T1_1(j,:),1:length(Time_index(j)),PRCC_var,alpha);
end

save("Model_LHS_INFalphaAntiPDL1_25ug_10_20_30_40.mat");

subplot(1,2,1)
bar(prccC(1,:)); 
k = 8;
set(gca,'XTickLabel',PRCC_var,'XTick',(1:k))
title('PRCCs');

subplot(1,2,2)
bar(pvalueC(1,:));
set(gca,'XTickLabel',PRCC_var,'XTick',(1:k))
title('p values');



