function dydt = ODE_treatment_INFalphaAntiPDL1(t,Y,par)
C = Y(1);
E = Y(2);
T_1 = Y(3);
T_2 = Y(4);
A = Y(5);

dydt(1) = par.lambdaC*C.*(1-C/par.K) - par.beta1*E.*C - par.beta2*T_1.*C;  %C
dydt(2) = par.alpha0.*par.T0 - par.dE.*E; %T
dydt(3) = (par.alpha1*(1+A/par.alphaA).*par.T0.*C + par.lambdaT*T_1)./(1 + par.delta*par.sigma*T_1.*(T_1+par.epl*C)*(1-A/(par.Ka+A)))-par.rho*T_1 - par.d1*T_1;%T_1
dydt(4) = par.rho*T_1 - par.d2*T_2;%H
dydt(5) = par.gamma - par.dA*A;%A(1+detla*sigma*T_1.*(T_1+epl*C))
dydt= dydt(:);
end