function dydt = ODE_control(~,Y,par)
C = Y(1); E = Y(2); T_1 = Y(3); T_2 = Y(4);
dydt(1) = par.lambdaC*C.*(1-C/par.K) - par.beta1*E.*C - par.beta2*T_1.*C;  %C
dydt(2) = par.alpha0.*par.T0 - par.dE.*E; %E
dydt(3) = (par.alpha1.*par.T0.*C + par.lambdaT*T_1)./(1 + par.delta*par.sigma*T_1.*(T_1+par.epl*C))-par.rho*T_1 - par.d1*T_1;%T_1
dydt(4) = par.rho*T_1 - par.d2*T_2;%T_2
dydt= dydt(:);
end