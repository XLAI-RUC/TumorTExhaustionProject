%% The results should be compared to the PRCC results section in
%% Supplementary Material D and Table D.1 for different N (specified by
%% "runs" in the script below
clear; clc;

%% Sample size N
runs = 10000;

%% LHS MATRIX  %%
Parameter_settings_LHS;
a = 0.5;
b = 2;
alpha_1_LHS= LHS_Call( a*7.7181e-11, alpha_1, b*7.7181e-11,  0 ,runs,'unif'); % baseline = 1.4637e-6
lambda_t_LHS=LHS_Call( a*0.7369, lambda_t,   b*0.7369,     0 ,runs,'unif'); % baseline = 0.4134
rho_LHS=     LHS_Call( a*0.1498, rho,        b*0.1498,   0 ,runs,'unif'); % baseline = 0.0342
d_1_LHS=     LHS_Call( a*0.41, d_1,        b*0.41,     0 ,runs,'unif'); % baseline = 0.3917
lambda_c_LHS=LHS_Call( a*0.7163, lambda_c, b*0.7163,  0 ,runs,'unif'); % baseline = 1.4637e-6
sigma_LHS=   LHS_Call( a*2.725e-17, sigma, b*2.725e-17,  0 ,runs,'unif'); % baseline = 1.4637e-6;
epl_LHS=     LHS_Call( a*0.01, epl, b*0.05,  0 ,runs,'unif'); % baseline = 1.4637e-6;
beta_2_LHS=  LHS_Call( a*1.5688e-8, beta_2, b*1.5688e-8,  0 ,runs,'unif'); % baseline = 1.4637e-6 1.5688e-8;

%% LHS MATRIX and PARAMETER LABELS
LHSmatrix = [alpha_1_LHS lambda_t_LHS rho_LHS d_1_LHS lambda_c_LHS sigma_LHS epl_LHS beta_2_LHS] ;

for x=1:runs %Run solution x times choosing different values
    f=@ODE_LHS;
    x;
    LHSmatrix(x,:);
    options = odeset('RelTol', 1e-4, 'AbsTol', [1e-4,1e-4,1e-4,1e-4,1e-4]);% 定义ODE求解器的选项
    [t,y]=ode15s(@(t,y)ODE_LHS(t,y,LHSmatrix,x,runs),tspan,y0,options); 
    B=[t y]; % [time y]
    %% Save the outputs at ALL time points [tspan]
    %T_lhs(:,x)=Anew(:,1);
    %CD4_lhs(:,x)=Anew(:,2);
    %T1_lhs(:,x)=Anew(:,3);
    %T2_lhs(:,x)=Anew(:,4);
    %V_lhs(:,x)=Anew(:,5);
    
    %% Save only the outputs at the time points of interest [time_points]:
    %% MORE EFFICIENT
    C_lhs1(:,x) = B(time_points1 + 1,2);
    E_lhs1(:,x) = B(time_points1 + 1,3);
    T_1_lhs1(:,x) = B(time_points1 + 1,4);
    T_2_lhs1(:,x) = B(time_points1 + 1,5);
    A_lhs1(:,x) = B(time_points1 + 1,6);
    
    C_lhs2(:,x) = B(time_points2 + 1,2);
    E_lhs2(:,x) = B(time_points2 + 1,3);
    T_1_lhs2(:,x) = B(time_points2 + 1,4);
    T_2_lhs2(:,x) = B(time_points2 + 1,5);
    A_lhs2(:,x) = B(time_points2 + 1,6);
    
    C_lhs3(:,x) = B(time_points3 + 1,2);
    E_lhs3(:,x) = B(time_points3 + 1,3);
    T_1_lhs3(:,x) = B(time_points3 + 1,4);
    T_2_lhs3(:,x) = B(time_points3 + 1,5);
    A_lhs3(:,x) = B(time_points3 + 1,6);
end
%% Save the workspace
% CALCULATE PRCC
alpha = 0.05;
[prcc11, pvalue11, sign_label11] = PRCC(LHSmatrix,C_lhs1,1:length(time_points1),PRCC_var,alpha);
[prcc12, pvalue12, sign_label12] = PRCC(LHSmatrix,T_1_lhs1,1:length(time_points1),PRCC_var,alpha);

[prcc21, pvalue21, sign_label21] = PRCC(LHSmatrix,C_lhs2,1:length(time_points2),PRCC_var,alpha);
[prcc22, pvalue22, sign_label22] = PRCC(LHSmatrix,T_1_lhs2,1:length(time_points2),PRCC_var,alpha);

[prcc31, pvalue31, sign_label31] = PRCC(LHSmatrix,C_lhs3,1:length(time_points3),PRCC_var,alpha);
[prcc32, pvalue32, sign_label32] = PRCC(LHSmatrix,T_1_lhs3,1:length(time_points3),PRCC_var,alpha);

subplot(1,2,1)
bar(prcc11(1,:)); 
k = 8;
set(gca,'XTickLabel',PRCC_var,'XTick',(1:k))
title('PRCCs');
subplot(1,2,2)
bar(pvalue11(1,:));
set(gca,'XTickLabel',PRCC_var,'XTick',(1:k))
title('p values');

save Model_LHS_30_60_120.mat;

