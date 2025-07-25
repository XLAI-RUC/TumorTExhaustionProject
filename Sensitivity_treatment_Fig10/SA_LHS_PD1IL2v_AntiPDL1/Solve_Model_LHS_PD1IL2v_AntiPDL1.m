%% The results should be compared to the PRCC results section in
%% Supplementary Material D and Table D.1 for different N (specified by
%% "runs" in the script below
clear;
clc;

%% Sample size N
runs = 20000;
Parameter_settings_LHS_PD1IL2v_AntiPDl1;

%% LHS MATRIX  %%
a = 0.5;
b = 2;
par =  setparameter();
rho_orig = 0.1498;
beta2_orig = 1.5688e-8;
rho_betterT1 = 0.0972;        % 'Better Effector'
beta2_betterT1 = 1.5981e-8;   % 'Better Effector'

%% PD-L1 Monotherapy
par.gamma = 6.25e-6;
par.zeta = 0;

alpha_1_LHS  = LHS_Call(a*par.alpha1, par.alpha1, 10*par.alpha1,  0 ,runs,'unif');   % baseline = 1.4637e-6
lambda_t_LHS = LHS_Call(a*par.lambdaT, par.lambdaT, b*par.lambdaT, 0 ,runs,'unif');  % baseline = 0.7369
rho_LHS      = LHS_Call(a*par.rho, par.rho,   b*par.rho,   0 ,runs,'unif');      % baseline = 0.1498
d_1_LHS      = LHS_Call(a*par.d1, par.d1,        b*par.d1,     0 ,runs,'unif');  % baseline = 0.41
lambda_c_LHS = LHS_Call(a*par.lambdaC, par.lambdaC, b*par.lambdaC,  0 ,runs,'unif'); % baseline = 0.7163
sigma_LHS    = LHS_Call(a*par.sigma, par.sigma, b*par.sigma,  0 ,runs,'unif');     % baseline = 2.725e-17;
epl_LHS      = LHS_Call(a*par.epl, par.epl, 15*par.epl,  0 ,runs,'unif');          % baseline = 1.4637e-6;
beta_2_LHS   = LHS_Call(a*par.beta2, par.beta2, b*par.beta2,  0 ,runs,'unif');    % baseline = 1.5688e-8;
% LHS MATRIX and PARAMETER LABELS
LHSmatrix_1 = [alpha_1_LHS lambda_t_LHS rho_LHS d_1_LHS lambda_c_LHS sigma_LHS epl_LHS beta_2_LHS];

%% PD1IL2v Monotherapy
par.gamma = 0;
par.zeta = 6.25e-6;%

alpha_1_LHS  = LHS_Call(a*par.alpha1, par.alpha1, 5*par.alpha1,  0 ,runs,'unif');   % baseline = 1.4637e-6
lambda_t_LHS = LHS_Call(a*par.lambdaT, par.lambdaT, b*par.lambdaT, 0 ,runs,'unif');  % baseline = 0.7369
rho_LHS      = LHS_Call(a*rho_orig, rho_orig,   b*rho_orig,   0 ,runs,'unif');       % baseline = 0.1498
d_1_LHS      = LHS_Call(a*par.d1, par.d1,        b*par.d1,     0 ,runs,'unif');      % baseline = 0.41
lambda_c_LHS = LHS_Call(a*par.lambdaC, par.lambdaC, b*par.lambdaC,  0 ,runs,'unif'); % baseline = 0.7163
sigma_LHS    = LHS_Call(a*par.sigma, par.sigma, b*par.sigma,  0 ,runs,'unif');       % baseline = 2.725e-17;
epl_LHS      = LHS_Call(a*par.epl, par.epl, 15*par.epl,  0 ,runs,'unif');            % baseline = 1.4637e-6;
beta_2_LHS   = LHS_Call(a*beta2_orig, beta2_orig, b*beta2_orig,  0 ,runs,'unif');    % baseline = 1.5688e-8;
% LHS MATRIX and PARAMETER LABELS
LHSmatrix_2 = [alpha_1_LHS lambda_t_LHS rho_LHS d_1_LHS lambda_c_LHS sigma_LHS epl_LHS beta_2_LHS] ;

%% PD1IL2v + AntiPDL1
par.zeta = 6.25e-6;
par.gamma = 6.25e-6;

alpha_1_LHS  = LHS_Call(a*par.alpha1, par.alpha1, 5*par.alpha1,  0 ,runs,'unif');        % baseline = 1.4637e-6
lambda_t_LHS = LHS_Call(a*par.lambdaT, par.lambdaT, b*par.lambdaT, 0 ,runs,'unif');      % baseline = 0.7369
rho_LHS      = LHS_Call(a*rho_betterT1, rho_betterT1,   b*rho_betterT1,   0 ,runs,'unif');       % baseline = 0.1498
d_1_LHS      = LHS_Call(a*par.d1, par.d1,        b*par.d1,     0 ,runs,'unif');          % baseline = 0.41
lambda_c_LHS = LHS_Call(a*par.lambdaC, par.lambdaC, b*par.lambdaC,  0 ,runs,'unif');     % baseline = 0.7163
sigma_LHS    = LHS_Call(a*par.sigma, par.sigma, b*par.sigma,  0 ,runs,'unif');           % baseline = 2.725e-17;
epl_LHS      = LHS_Call(a*par.epl, par.epl, 15*par.epl,  0 ,runs,'unif');                % baseline = 1.4637e-6;
beta_2_LHS   = LHS_Call(a*beta2_betterT1, beta2_betterT1, b*beta2_betterT1,  0 ,runs,'unif');    % baseline = 1.5688e-8;
% LHS MATRIX and PARAMETER LABELS
LHSmatrix_3 = [alpha_1_LHS lambda_t_LHS rho_LHS d_1_LHS lambda_c_LHS sigma_LHS epl_LHS beta_2_LHS] ;

options = odeset('RelTol', 1e-6, 'AbsTol', [1e-6,1e-6,1e-6,1e-6,1e-6,1e-6]);
Time_index = [time_points1,time_points2,time_points3,time_points4];


%% Solve ODE
for x = 1:runs % Run solution x times choosing different values
    
    %% Anti-PDL1 Monotherapy
    gamma0 = 6.25e-6; zeta0 = 0;
    [t, y] = ode45(@(t,y) ODE_LHS_PD1IL2v_AntiPDL1(t,y,LHSmatrix_1(x,:),par,gamma0,zeta0), tspan,y0, options);%解微分方程
    B = [t y];   % [time y]
    % Save only the outputs at the time points of interest [time_points]:
    C_1(:,x) = B(Time_index+1,2);
    E_1(:,x) = B(Time_index+1,3);
    T1_1(:,x) = B(Time_index+1,4);
    T2_1(:,x) = B(Time_index+1,5);
    A_1(:,x) = B(Time_index+1,6);
    B_1(:,x) = B(Time_index+1,7);
    
    %% PD1IL2v Monotherapy
    gamma0 = 0; zeta0 = 6.25e-6;
    [t, y] = ode45(@(t,y) ODE_LHS_PD1IL2v_AntiPDL1(t,y,LHSmatrix_2(x,:),par,gamma0,zeta0), tspan,y0, options);%解微分方程
    B = [t y];   % [time y]
    % Save only the outputs at the time points of interest [time_points]:
    C_2(:,x) = B(Time_index+1,2);
    E_2(:,x) = B(Time_index+1,3);
    T1_2(:,x) = B(Time_index+1,4);
    T2_2(:,x) = B(Time_index+1,5);
    A_2(:,x) = B(Time_index+1,6);
    B_2(:,x) = B(Time_index+1,7);
    
    %% PD1IL2v + Anti-PD-L1
    gamma0 = 6.25e-6; zeta0 = 6.25e-6;
    [t, y] = ode45(@(t,y) ODE_LHS_PD1IL2v_AntiPDL1(t,y,LHSmatrix_3(x,:),par,gamma0,zeta0), tspan,y0, options);%解微分方程
    B = [t y];   % [time y]
    % Save only the outputs at the time points of interest [time_points]:
    C_3(:,x) = B(Time_index+1,2);
    E_3(:,x) = B(Time_index+1,3);
    T1_3(:,x) = B(Time_index+1,4);
    T2_3(:,x) = B(Time_index+1,5);
    A_3(:,x) = B(Time_index+1,6);
    B_3(:,x) = B(Time_index+1,7);     
end

%% Save the workspace
% CALCULATE PRCC

alpha = 0.05;
for j = 1:length(Time_index)
    [prccC_1(j,:), pvalueC_1(j,:), sign_labelC_1(j,:)] = PRCC(LHSmatrix_1,C_1(j,:),1:length(Time_index(j)),PRCC_var,alpha);
    [prccT1_1(j,:), pvalueT1_1(j,:), sign_labelT1_1(j,:)] = PRCC(LHSmatrix_1,T1_1(j,:),1:length(Time_index(j)),PRCC_var,alpha);
    
    [prccC_2(j,:), pvalueC_2(j,:), sign_labelC_2(j,:)] = PRCC(LHSmatrix_1,C_2(j,:),1:length(Time_index(j)),PRCC_var,alpha);
    [prccT1_2(j,:), pvalueT1_2(j,:), sign_labelT1_2(j,:)] = PRCC(LHSmatrix_1,T1_2(j,:),1:length(Time_index(j)),PRCC_var,alpha);
    
    [prccC_3(j,:), pvalueC_3(j,:), sign_labelC_3(j,:)] = PRCC(LHSmatrix_1,C_3(j,:),1:length(Time_index(j)),PRCC_var,alpha);
    [prccT1_3(j,:), pvalueT1_3(j,:), sign_labelT1_3(j,:)] = PRCC(LHSmatrix_1,T1_3(j,:),1:length(Time_index(j)),PRCC_var,alpha); 
end


%% Save data
save("Model_LHS_PD1IL2v_AntiPDL1_10_20_30_40.mat");

subplot(1,2,1)
bar(prccC_1(1,:)); 
k = 8;
set(gca,'XTickLabel',PRCC_var,'XTick',(1:k))
title('PRCCs');

subplot(1,2,2)
bar(pvalueC_1(1,:));
set(gca,'XTickLabel',PRCC_var,'XTick',(1:k))
title('p values');



