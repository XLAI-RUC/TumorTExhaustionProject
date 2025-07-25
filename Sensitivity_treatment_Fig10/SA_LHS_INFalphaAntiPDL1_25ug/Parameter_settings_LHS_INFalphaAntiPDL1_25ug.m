
% Parameter Labels 
PRCC_var={'\alpha_1', '\lambda_T', '\rho','d_1','\lambda_C','\sigma','\epsilon','\beta_2' };%¦Ñ¦Á_1 

%% TIME SPAN OF THE SIMULATION
c = 40;
t_end = c; % length of the simulations
tspan = (0:t_end);   % time points where the output is calculated
time_points1 = 10; % time points of interest for the US analysis
time_points2 = 20; % time points of interest for the US analysis
time_points3 = 30; % time points of interest for the US analysis
time_points4 = 40;
% INITIAL CONDITION FOR THE ODE MODEL
y0 =[1e6,5000,1.375*1e7,1e7,0];%C¡¢S¡¢E¡¢T cellsÏ¸°ûµÄ³õÖµ

% Variables Labels
y_var_label={'C','E','T_1','T_2','A'};