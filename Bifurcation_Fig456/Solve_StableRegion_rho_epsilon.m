clc; clear; close all;

par = setparameter();

epl_Vec = 0.01:0.001:0.15;    % epsilon
rho_Vec = 0.01:0.001:0.15;    % rho

length_e = length(epl_Vec);
length_r = length(rho_Vec);

C_matrix = zeros(length_e,length_r); 
T1_matrix = zeros(length_e,length_r);  

t_span = [0 1000];
options = odeset('RelTol', 1e-6, 'AbsTol', [1e-5,1e-5,1e-5,1e-5]); 

y0 = [1e9,5e3,1e7,1e5];
% y0 = [1e3,5e3,1e7,1e5];

for ii = 1:length_e
    par.rho = rho_Vec(ii);
    for jj = 1:length_r
        par.epl = epl_Vec(jj);
        
        [t, y] = ode45(@(t,y) ODE_control(t,y,par), t_span,y0, options);  %解微分方程
        
        C_matrix(ii,jj) = y(end,1);
        T1_matrix(ii,jj) = y(end,3);
    end
end

subplot(1,2,1)
surf(rho_Vec, epl_Vec, C_matrix)
xlabel('\epsilon');
ylabel('\rho');
zlabel('C')

subplot(1,2,2)
surf(rho_Vec, epl_Vec,T1_matrix)
xlabel('\epsilon');
ylabel('\rho');
zlabel('T_1')


%% Save the data
File_folder_name = 'output/3D/'; 
dlmwrite([File_folder_name, 'Bifur_epl_rho_C.dat'],C_matrix,' ');
dlmwrite([File_folder_name, 'Bifur_epl_rho_T1.dat'],T1_matrix,' ');




