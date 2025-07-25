clear; clc;


%%%%%%%%%%%% Fig. 7b Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For control case, consider the effect of exhaustion \rho （Fig.7b）
rho_values = [0.12, 0.15, 0.18];  
num_rho = length(rho_values);  
options = odeset('RelTol', 1e-6, 'AbsTol', [1e-6,1e-6,1e-6,1e-6,1e-6]);  
 
t1 = (0:1:60);  % time 
sol = zeros(length(t1),3);    
for j = 1:num_rho
    par = setparameter();
    par.rho = rho_values(j);
    gamma0 = 0; %3.125e-6;   % No treatment
    y0 = [1e6, 5000, 1.3*10^7, 10^7, 0];
    [t, R] = ode45(@(t,y) ODE_treatment_AntiPDL1(t,y,par,gamma0), t1, y0, options);
    sol(:,j) = R(:,1);  % Tumor cell
end
 
plot(t1, sol(:, 1)); hold on;
plot(t1, sol(:, 2));
plot(t1, sol(:, 3));
xlabel('Time (day)')  
ylabel('Tumor cell (C)')  
% Save data
dlmwrite('Data/Tumor_t_rho.dat',[t1',sol],' ');


%%%%%%%%%%% Fig. 7d,e Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Under different dosage of anti-PD-L1, consider the effect of exhaustion \rho （Fig.7d,e）
gamma0 = (0:1:150)* 1e-6/4;   % anti-PD-L1 dosage 
n = length(gamma0);  
rho_values = [0.12, 0.15, 0.18];  
num_rho = length(rho_values);  
C = zeros(num_rho, n);   % Tumor cells
T1 = zeros(num_rho, n);  % T1 cells

options = odeset('RelTol', 1e-6, 'AbsTol', [1e-6,1e-6,1e-6,1e-6,1e-6]);  
sol = zeros(); 
for j = 1:num_rho  
    par = setparameter(); 
    par.rho = rho_values(j);

    for i = 1:n  
        y0 = [1e6, 5000, 1.3*10^7, 10^7, 0];  
        t1 = (0:1:400);  
        [t, R] = ode45(@(t,y) ODE_treatment_AntiPDL1(t,y,par,gamma0(i)), t1, y0, options);  
        C(j, i) = R(length(t1),1); 
        T1(j, i) = R(length(t1),3); 
    end  
 
    if j == 1  
        plot(gamma0, C(j, :), 'r-') 
        hold on 
    else  
        plot(gamma0, C(j, :)) 
    end  
end  
 
xlabel('Dose of anti-PD-L1 (ug)')  
ylabel('Tumor cell (C)')  
% Save data
dlmwrite('Data/Cell_gamma_rho.dat',[gamma0',C',T1'],' ');

%%%%%%%%%%%%%%  Fig.7f Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par = setparameter();
dx = 0.01; 
gamma0 =(0:dx:2)*1e-5;  %anti-PD-L1的剂量
n = length(gamma0);
rho = linspace(0,0.2,n); 
C = zeros(n,n); T1 = zeros(n,n);
y0 = [1e6,5000,1.3*10^7,10^7,0];
options = odeset('RelTol', 1e-6, 'AbsTol', [1e-6,1e-6,1e-6,1e-6,1e-6]);% 定义ODE求解器的选项
t1 = (0:10:400);

% distcomp.feature( 'LocalUseMpiexec', false);
% poolobj = gcp('nocreate');    % If no pool, do not create new one.
% if isempty(poolobj)
%     poolobj = parpool('local',n);
% end

for i = 1:n
    for j = 1:n
        par = setparameter();
        rho = linspace(0,0.2,n); 
        gamma0 =(0:dx:2)*1e-5;   % anti-PD-L1 dosage
        par.rho = rho(j);
        [t, R] = ode45(@(t,y) ODE_treatment_AntiPDL1(t,y,par,gamma0(i)), t1,y0, options);
        C(i,j) = R(end,1);   % Tumor cells
        T1(i,j) = R(end,3);  % T1 cells
    end
end
% Save data
save('Data_rho_gamma','C','T1','gamma0','rho');


