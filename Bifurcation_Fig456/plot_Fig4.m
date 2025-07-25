
close all;
FS_lab = 10;

hf = figure;

subplot(221);
plot_Fig4_bifur_C()
text(-0.018,2.04e9,'\bf{a}','FontSize',FS_lab);

subplot(222)
plot_Fig4_bifur_T1()
text(-0.013,6.1e7,'\bf{b}','FontSize',FS_lab);

subplot(223)
plot_Fig4_bifur_epl()
text(-0.018,2.1e9,'\bf{c}','FontSize',FS_lab);

subplot(224)
plot_Fig4_SS_epsilon_rho()
text(-0.015,0.155,'\bf{d}','FontSize',FS_lab);

% print(hf,'-r700','-dpdf','Figure/Fig4.pdf');