function dydt = ODE_treatment_PD1IL2v_AntiPDL1(t,Y,par,rho,beta2)
C = Y(1);
E = Y(2);
T_1 = Y(3);
T_2 = Y(4);
A = Y(5);
B = Y(6);

par.rho = rho;
par.beta2 = beta2;
     
dydt(1) = par.lambdaC*C.*(1-C/par.K) - par.beta1*E.*C - par.beta2*T_1.*C;  %C
dydt(2) = par.alpha0.*par.T0 - par.dE.*E; %T
dydt(3) = (par.alpha1.*par.T0.*C + par.lambdaT*T_1)./(1 + par.delta*par.sigma*T_1.*(T_1+par.epl*C)*(1-B/(par.Kb+B))*(1-A/(par.Ka+A)))-par.rho*T_1 - par.d1*T_1; %T_1
dydt(4) = par.rho*T_1 - par.d2*T_2; %T2
dydt(5) = par.gamma - par.dA*A;  %A
dydt(6) = par.zeta - par.dB*B;   %B
dydt= dydt(:);

end