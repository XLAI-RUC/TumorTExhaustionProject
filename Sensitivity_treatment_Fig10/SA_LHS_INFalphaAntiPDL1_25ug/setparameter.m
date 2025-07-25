function par = setparameter()

%%% cancer cell (C)
par.lambdaC = 0.7163;
par.K = 3e9;
par.beta1 = 1.961e-8;
par.beta2 = 1.5688e-8;

%%% CTL (E)
par.alpha0 = 2.15e-4;
par.T0 = 5.59e6;
par.dE = 0.18;

%%% PD-1low and PD-1high T cell (T1, T2)
par.alpha1 = 7.7181e-11;
par.lambdaT = 0.7369;
par.delta = 10;
par.sigma = 2.725e-17;
par.epl = 0.01;
par.rho = 0.1498;
par.d1 = 0.41;
par.d2 = 0.72;

%%% anti_PD-L1
par.dA = 0.47;
par.Ka = 2e-5; 
par.gamma = 3.125e-6;

%%%
par.dB = 0.1932;
par.Kb = 2.9384e-05;
par.zeta = 3.571e-6;
par.alphaA = 2.9e-3;

end 