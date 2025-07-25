function Solve_Bifurcation_nondim_rho_delta

epl_0 = 0.05; % epl_0 = 0.05;    
% rho_0 = 0.05; 
delta_0 = 10; 

rho_Vec = 0.001:0.0005:0.15; 

length_r = length(rho_Vec);

RootC1 = zeros(length_r,2);
RootC2 = zeros(length_r,2);
RootC3 = zeros(length_r,2);

for ind = 1:length_r
    
rho_0 = rho_Vec(ind);  

fun = @(x)root2d(x,epl_0,rho_0,delta_0);

x1 = [0.05,0.05];
RootC1(ind,1:2) = fsolve(fun,x1);
x2 = [0.2,0.2];
RootC2(ind,1:2) = fsolve(fun,x2);
x3 = [1,1];
RootC3(ind,1:2) = fsolve(fun,x3);

end

figure

RootC1(RootC1<0) = nan;
RootC2(RootC2<0) = nan;
RootC3(RootC3<0) = nan;

C = 1;% 3e9;

plot(rho_Vec, RootC1(:,1)*C, 'b.','LineWidth',2); 
hold on 
plot(rho_Vec, RootC2(:,1)*C,'r.', 'LineWidth',2); 
plot(rho_Vec, RootC3(:,1)*C, 'g.','LineWidth',2); 
% plot([rho_Vec(1) rho_Vec(end)],[Ceq Ceq],'--') 

xlabel('\rho');
ylabel('C');


% File_folder_name = 'output/2D/epl_0p01/'; 
% dlmwrite([File_folder_name, 'Bifur_epl_0p01_delta_12_sol1.dat'],RootC1,' ');
% dlmwrite([File_folder_name, 'Bifur_epl_0p01_delta_12_sol2.dat'],RootC2,' ');
% dlmwrite([File_folder_name, 'Bifur_epl_0p01_delta_12_sol3.dat'],RootC3,' ');

% RootC1 = load([File_folder_name, 'Bifur_epl_0p05_delta_6_sol1.dat'])

end

function Fun = root2d(x,epl_0,rho_0,delta_0)

par = setparameter();

par.epl = epl_0; 
par.rho = rho_0; 
par.delta = delta_0; 

beta1_t = par.K*par.beta1;
beta2_t = par.K*par.beta2;
alpha0_t = par.alpha0*par.T0./par.K;

sigma_h = par.delta*par.sigma*par.K^2;
alpha1_h = par.alpha1*par.T0;

lambdaC_h = par.lambdaC./beta2_t;
Eeq = alpha0_t/par.dE;
Ceq = (par.lambdaC - beta1_t*Eeq)./par.lambdaC;

rho_h = par.rho + par.d1;

Fun(1) = beta2_t*x(2) - par.lambdaC.*(Ceq - x(1));
Fun(2) = (alpha1_h*x(1) + par.lambdaT*x(2))./(1 + sigma_h*x(2).*(x(2) + par.epl*x(1))) - rho_h*x(2);
end