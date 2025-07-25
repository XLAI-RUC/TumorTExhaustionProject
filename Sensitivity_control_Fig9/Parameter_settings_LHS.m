% PARAMETER BASELINE VALUES
alpha_1 = 7.7181e-11; 
lambda_t = 0.7369;
rho = 0.1498;
d_1 = 0.41;
lambda_c = 0.7163;
sigma = 2.725e-17;
epl = 0.01;
beta_2 = 1.5688e-8;

% Parameter Labels 
PRCC_var={'\alpha_1', '\lambda_T', '\rho','d_1','\lambda_C','\sigma','\epsilon','\beta_2' };%¦Ñ¦Á_1 

%% TIME SPAN OF THE SIMULATION
c = 120;
t_end = c; % length of the simulations
tspan = (0:t_end);   % time points where the output is calculated
time_points1 = 30; % time points of interest for the US analysis
time_points2 = 60; % time points of interest for the US analysis
time_points3 = 120; % time points of interest for the US analysis

% INITIAL CONDITION FOR THE ODE MODEL
y0 =[1e6,500,1.375e7,1e7,0];

% Variables Labels
y_var_label={'C','E','T_1','T_2','A'};