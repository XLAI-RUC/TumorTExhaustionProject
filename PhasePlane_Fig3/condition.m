function [num,V,S,p,a1,a2,a3,a4,g1,g2,g3,g4,f1,f2,f22,Cps,Cmi,phi_Cps,phi_Cmi,C,T_1,Eeq] = condition(par,t1_span,y0,s)
clc;

Eeq = par.alpha0*par.T0./par.dE;
rho_h = par.rho + par.d1;
sigma_h = par.delta*par.sigma;
lambdaC_h = par.lambdaC./(par.K*par.beta2);
alpha1_h = par.alpha1*par.T0;

f1 = par.lambdaC - par.beta1 * Eeq;
f2 = par.lambdaT - rho_h;
T10_ceq = sqrt((par.lambdaT-rho_h)./(rho_h*sigma_h));
f22 = par.lambdaC-par.beta1*Eeq-par.beta2*T10_ceq;
% f21 = par.lambdaC - par.beta1 * EquilE + par.beta2 * sqrt((par.lambdaT - rho_h)./(rho_h*sigma_h))  % par.lambdaT - rho > 0

Ceq = par.K*(par.lambdaC - par.beta1*Eeq)./par.lambdaC;
ss = sigma_h*rho_h*lambdaC_h.^2;

g1 = par.epl - 2*lambdaC_h  ;
g2 = par.lambdaT - rho_h - ss*Ceq.^2;
g3 =  ss*Ceq.^2*g1-alpha1_h;
g4 = ss *Ceq.^2* par.epl.^2 - 3*(lambdaC_h - par.epl).*(alpha1_h-lambdaC_h*(rho_h - par.lambdaT)) ;

a1 = ss*(lambdaC_h - par.epl);
a2 = 2*ss*Ceq*(par.epl - 1.5*lambdaC_h);
a3 = alpha1_h + ss * Ceq^2 *(3*lambdaC_h - par.epl) + lambdaC_h*(rho_h - par.lambdaT);
a4 = - par.lambdaC *Ceq*(ss*Ceq^2 + rho_h - par.lambdaT);

p = [a1, a2, a3,a4]; 

Cps = (-a2 + sqrt(a2.^2 - 3*a1*a2))./(3*a1);
Cmi = (-a2 + sqrt(a2.^2 - 3*a1*a2))./(3*a1);
phi_Cps = a1*Cps.^3 + a2*Cps.^2 + a3*Cps + a4;
phi_Cmi = a1*Cmi.^3 + a2*Cmi.^2 + a3*Cmi + a4;

options = odeset('RelTol', 1e-3, 'AbsTol', [1e-3,1e-3,1e-3,1e-3]);   % 定义ODE求解器的选项
C = zeros(length(t1_span),length(s));
T_1 = zeros(length(t1_span),length(s));

for i = 1:length(s)
    [~, R] = ode23s(@(t,y) ODE_control(t,y,par), t1_span,y0(i,:), options);  %解微分方程
    C(:,i) = R(:,1);
    T_1(:,i) = R(:,3);
end

%条件判断
num=0;
V=0;
if f2<=0 %满足第一种正平衡点存在条件
    V=1;
    num=num+1;
end
if f2>0 && g2<0   %满足第二种正平衡点存在条件
    V=2;
    num=num+1;
end
if f2>0 && g1>0 && abs(g2)<1e-6 &&g3>0   %满足第三种正平衡点存在条件
    V=3;
    num=num+1;
end
if f2>0 && g1>0 &&g2>0 &&g3>0 &&g4>0    %满足第四种正平衡点存在条件
    V=4;
    num=num+1;
end

S=0;
if f1<0 &&f2<0 %满足第一种正平衡点稳定条件
    S=1;
elseif f2>0 && f22<0 %满足第二种边界平衡点稳定条件
    S=2;
end

end