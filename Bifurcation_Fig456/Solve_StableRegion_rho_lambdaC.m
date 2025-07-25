clc; clear; close all;

par = setparameter();

par.epl = 0.01;

rho_Vec = 0.01:0.001:0.15;    % rho
length_r = length(rho_Vec);
var_Vec = linspace(0.5,0.8,length_r);    % lambda_C

C_matrix = zeros(length_r,length_r); 
T1_matrix = zeros(length_r,length_r);  

t_span = [0 1000];
options = odeset('RelTol', 1e-6, 'AbsTol', [1e-5,1e-5,1e-5,1e-5]); 

%y0 = [2e9,5e3,1e7,1e5];
y0 = [1e2,5e3,1e7,1e2];

for ii = 1:length_r
    par.rho = rho_Vec(ii);
    for jj = 1:length_r
        par.lambdaC = var_Vec(jj);
        
        [t, y] = ode45(@(t,y) ODE_control(t,y,par), t_span,y0, options);  %解微分方程
        
        C_matrix(ii,jj) = y(end,1);
        T1_matrix(ii,jj) = y(end,3);
    end
end

surf(rho_Vec,var_Vec, C_matrix)
xlabel('\epsilon');
ylabel('\rho');
zlabel('C')


figure
surf(rho_Vec, var_Vec,T1_matrix)
xlabel('\epsilon');
ylabel('\rho');
zlabel('T_1')

File_folder_name = 'output/3D/epl_0p01/'; 
dlmwrite([File_folder_name, 'Bifur_rho_lambdaC_C_2.dat'],C_matrix,' ');
dlmwrite([File_folder_name, 'Bifur_rho_lambdaC_T1_2.dat'],T1_matrix,' ');




