clear; clc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bistable
par = setparameter();
par.rho = 0.05;   par.epl = 0.05;

% Set initial values 
e1 = [5,36,7,5,5,5]*1e7;  
e2 = (1e8:2e8:12e8);
h = ones(1,length(e1)+length(e2))*5000;
s1 = [33,61,22,25,10,14]*1e6; 
s2 = (6.25e7:-0.5e6:6e7);
y0 = [[e1,e2];h;[s1,s2];[s1,s2]]';

t1_span = (0:0.5:4000);
[num,V,S,p,a1,a2,a3,a4,g1,g2,g3,g4,f1,f2,f22,Cps,Cmi,phi_Cps,phi_Cmi,C,T_1,Eeq] = condition(par,t1_span,y0,[s1,s2]);
% V -- Existence of V-th positive equilibrium£¬S -- Existence of S-th boundary equilibrium 

roots = roots(p);  
real_roots = roots(imag(roots) == 0);
C_equ = real_roots(real_roots>0);
T_1_equ = (par.lambdaC*(1 - C_equ/par.K) - par.beta1*Eeq)/par.beta2;

h1 = subplot(231);
PhasePlane(C,T_1,h1)
save('output/data_exist');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
par = setparameter();
par.rho = 0.4;

% Set initial values 
g = [5e7, 1e8,1e9,3e9,1e8,3e9];
h = ones(1,length(g))*5000;
s=[9.6e6,1e7,1e7,1e7,6e6,6e6];
y0 =[g;h;s;s]';

t1_span=(0:0.5:4000);
[num,V,S,p,a1,a2,a3,a4,g1,g2,g3,g4,f1,f2,f22,Cps,Cmi,phi_Cps,phi_Cmi,C,T_1,Eeq]=condition(par,t1_span,y0,s);

% roots = roots(p);
% real_roots = roots(imag(roots) == 0);
% C_equ=real_roots(real_roots>0);
% T_1_equ=(par.lambdaC*(1 - C_equ/par.K) - par.beta1*Eeq)/par.beta2;

h2 = subplot(232);
PhasePlane(C,T_1, h2)
save('output/data_exist1');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

par = setparameter();
e = [1e7,1e8,1e9,1e9,2e9];
h = ones(1,length(e))*5000;
s = [1.1e7,1.15e7,1e7,5e7,3e7];
y0 =[e;h;s;s]';

t1_span = (0:0.5:4000);
[num,V,S,p,a1,a2,a3,a4,g1,g2,g3,g4,f1,f2,f22,Cps,Cmi,phi_Cps,phi_Cmi,C,T_1,Eeq] = condition(par,t1_span,y0,s);

% roots = roots(p);
% real_roots = roots(imag(roots) == 0);
% C_equ = real_roots(real_roots>0);
% T_1_equ = (par.lambdaC*(1 - C_equ/par.K) - par.beta1*Eeq)/par.beta2;

h3 = subplot(233);
PhasePlane(C,T_1,h3)
save('output/data_exist2');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
par = setparameter();
par.epl = 0.05; par.rho = 0.06;

e = [1e7,1e8,1e9,1e9,2e9];
h = ones(1,length(e))*5000;
s = [1.1e7,1e7,1.15e7,5e7,3e7];
y0 = [e;h;s;s]';

t1_span=(0:0.5:4000);
[num,V,S,p,a1,a2,a3,a4,g1,g2,g3,g4,f1,f2,f22,Cps,Cmi,phi_Cps,phi_Cmi,C,T_1,Eeq]=condition(par,t1_span,y0,s);

% roots = roots(p);
% real_roots = roots(imag(roots) == 0);
% C_equ=real_roots(real_roots>0);
% T_1_equ=(par.lambdaC*(1 - C_equ/par.K) - par.beta1*Eeq)/par.beta2;

hf4 = subplot(234);
PhasePlane(C,T_1,hf4)
save('output/data_exist3');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
par = setparameter();
par.rho = 0.02; 

e = [1e8,1e9,1e9,2e9];
h = zeros(1,length(e));
for i = 1:length(e)
h(i)=5000;
end
s = [1e7,1.1e7,5e7,3e7];
y0 = [e;h;s;s]';

t1_span=(0:0.5:4000);
[num,V,S,p,a1,a2,a3,a4,g1,g2,g3,g4,f1,f2,f22,Cps,Cmi,phi_Cps,phi_Cmi,C,T_1,Eeq]=condition(par,t1_span,y0,s);

% roots = roots(p);
% real_roots = roots(imag(roots) == 0);
% C_equ=real_roots(real_roots>0);
% T_1_equ=(par.lambdaC*(1 - C_equ/par.K) - par.beta1*Eeq)/par.beta2;

hf5 = subplot(235);
PhasePlane(C,T_1,hf5)
save('output/data_exist5');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
par = setparameter();
par.lambdaC = 0.01;  par.beta1 = 1e-5;  par.lambdaT = 0.01;

e = [3e6:2e6:12e6];
h = ones(1,length(e))*5000;
s = [9,8.5,7,6,3]*1e+6; %[1e6:2e6:10e6];
y0 =[e;h;s;s]';

t1_span=(0:0.5:4000);
[num,V,S,p,a1,a2,a3,a4,g1,g2,g3,g4,f1,f2,f22,Cps,Cmi,phi_Cps,phi_Cmi,C,T_1,Eeq] = condition(par,t1_span,y0,s);
% roots = roots(p);
% real_roots = roots(imag(roots) == 0);
% C_equ=real_roots(real_roots>0);
% T_1_equ=(par.lambdaC*(1 - C_equ/par.K) - par.beta1*Eeq)/par.beta2;

hf6 = subplot(236);
PhasePlane(C,T_1,hf6)
save('output/data_exist6');

