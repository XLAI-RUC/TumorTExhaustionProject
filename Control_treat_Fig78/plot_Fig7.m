clc; close all; clear all;
    
    FS_inside = 8;   % FontSize inside figure
    FS_default = 10;  % FontSize default
    FS_lab = 12;      % FontSize of (A) (B) (C)
    colorM = [0 0.4470 0.7410; 0.9 0 0; 0.4660 0.6740 0.1880; 0.7 0.5 0; 0 0 1; 0.7350 0.1780 0.1840;1 0.5 0.15];
    set(0,'defaultaxesfontsize',FS_default,'defaultlinelinewidth',1.3);
    
    hf1 = figure;
    %% Tumor growth curve
    subplot(2,3,1);
    sol = load('Data/Control_simulation.dat');  hold on;
    data_exp = load('Data/Control_experiment.dat'); 
    time_sim = sol(:,1);  C_sim = sol(:,2);  
    time_exp = data_exp(:,1);  C_exp = data_exp(:,2);
    hControl = plot(time_sim, C_sim, 'Color', colorM(1,:));  
    hs1 = scatter(time_exp, C_exp, 'MarkerFaceColor', colorM(1,:),'MarkerEdgeColor', colorM(1,:),'Sizedata', 10); 
    xlabel('Time (days)','fontname','Times','FontSize', FS_default);  
    ylabel('Tumor cells (C) ','fontname','Times','FontSize', FS_default);  
    text(-12.2,10.5e+8,'\bf{a}','FontSize',FS_lab);
    leg = legend([hs1,hControl], {'Experiment','Simulation'},'FontSize',6);  legend boxoff;
    leg.ItemTokenSize = [8,7];
    lgd.FontSize = 6;  
    pos = get(leg, 'Position');
    set(leg, 'Position', [pos(1)-0.04 pos(2)+0.04 0 0]);
    box on;
    
    %% Effect of dosage \gamma on tumor growth
    subplot(2,3,3);
    sol = load('Data/Tumor_gamma_rho_0p15.dat'); 
    t1_span = sol(:,1);   C1 = sol(:,2);  C2 = sol(:,3);  C3 = sol(:,4);
    hp1 = plot(t1_span(20:end), C1(20:end), '-^','MarkerIndices',1:20:length(C1),'Color', colorM(3,:),'MarkerSize',3);  hold on 
    hp2 =  plot(t1_span(20:end), C2(20:end),  '-d','MarkerIndices',1:20:length(C2),'Color', colorM(4,:),'MarkerSize',3);  
    hp3 =  plot(t1_span(20:end), C3(20:end),'-h','MarkerIndices',1:20:length(C3),'Color', colorM(2,:),'MarkerSize',3);  
    
    sol = load('Data/Control_simulation.dat'); 
    t1_span = sol(:,1); C = sol(:,2);
    hp4 =  plot(t1_span, C, '-o','MarkerIndices',1:20:length(C), 'Color', colorM(1,:),'MarkerSize',3);  
    leg = legend([hp4,hp1,hp2,hp3], {'Control','\gamma = 3.2\times 10^{-6}','\gamma = 5.2\times 10^{-6}','\gamma = 10\times 10^{-6}'}, 'Location', 'northwest');  
    legend boxoff;
    leg.ItemTokenSize = [10,10];  % shorter line
    pos = get(leg, 'Position');
    set(leg, 'Position', [pos(1) pos(2)+0.05 0.1 0.1]);
    set(leg,'FontSize',7)
    xlabel('Time (days)','fontname','Times','FontSize', FS_default);  
    ylabel('Tumor cells (C) ','fontname','Times','FontSize', FS_default);  
    text(-12.2,10.5e+8,'\bf{c}','FontSize',FS_lab);
  
    
    %% Effect of exhaustion rate \rho on tumor growth
    subplot(2,3,2);
    sol = load('Data/Tumor_t_rho.dat');
    t1 = sol(:,1);  sol = sol(:,(2:4));
    plot(t1, sol(:,1), 'Color', colorM(7,:));  
    hold on;  
    plot(t1, sol(:,2),'-o','MarkerIndices',1:10:size(sol,1), 'Color', colorM(1,:),'MarkerSize',3);  
    plot(t1, sol(:,3), 'Color', colorM(3,:));  
    xlim([0,60]); ylim([0 13e+8]);
    text(40,1.31e+9,'\rho = 0.18','Interpreter','tex','FontSize',FS_inside);
    text(40,1.02e+9,'\rho = 0.15','Interpreter','tex','FontSize',FS_inside);
    text(40,0.65e+9,'\rho = 0.12','Interpreter','tex','FontSize',FS_inside);
    xlabel('Time (days)','fontname','Times','FontSize', FS_default);  
    ylabel('Tumor cells (C) ','fontname','Times','FontSize', FS_default);  
    text(-12,13.5e+8,'\bf{b}','FontSize',FS_lab); 
    
    %%
    subplot(2,3,4);
    sol = load('Data/Cell_gamma_rho.dat');
    gamma0 = sol(:,1);  C = sol(:,(2:4));   T1 = sol(:,(5:7));
    plot(gamma0, C(:,1), 'Color', colorM(7,:),'MarkerSize',3);  
    hold on
    plot(gamma0, C(:,2), '-o','MarkerIndices',1:10:size(C,1),'Color',colorM(1,:),'MarkerSize',3);  
    plot(gamma0, C(:,3), 'Color', colorM(3,:),'MarkerSize',3);  
    xlabel('\gamma','fontname','Times','FontSize', FS_default)
    ylabel('Tumor cells (C)','fontname','Times','FontSize', FS_default)
    text(0.25e-5,0.4e+9,'\rho = 0.12','Interpreter','tex','FontSize',FS_inside,'Rotation', -52);
    text(0.45e-5,0.53e+9,'\rho = 0.15','Interpreter','tex','FontSize',FS_inside,'Rotation', -50);
    text(0.6e-5,0.75e+9,'\rho = 0.18','Interpreter','tex','FontSize',FS_inside,'Rotation', -50);
    xlim([0,1.25e-5]);  ylim([0 13e+8]);
    text(-0.2e-5,13.6e+8,'\bf{d}','FontSize',FS_lab);

    %%
    subplot(2,3,5);
    plot(gamma0, T1(:,1), 'Color', colorM(7,:),'MarkerSize',3);  
    hold on
    plot(gamma0, T1(:,2), '-o','MarkerIndices',1:30:size(C,1),'Color',colorM(1,:),'MarkerSize',3);  
    plot(gamma0, T1(:,3), 'Color', colorM(3,:),'MarkerSize',3);  
    xlabel('\gamma','fontname','Times','FontSize', FS_default)
    ylabel('PD-1^{lo}CD8^+ T_{ex} (T_1)','fontname','Times','FontSize', FS_default)
    text(1e-5,5.8e+7,'\rho = 0.12','Interpreter','tex','FontSize',FS_inside,'Rotation', 37);
    text(1.25e-5,5.5e+7,'\rho = 0.15','Interpreter','tex','FontSize',FS_inside,'Rotation', 33);
    text(1.5e-5,5.1e+7,'\rho = 0.18','Interpreter','tex','FontSize',FS_inside,'Rotation', 30);
    xlim([0,3.7e-5]);  ylim([2.6e+7 8.45e+7]);
    text(-1e-5,8.8e+7,'\bf{e}','FontSize',FS_lab);
    
    %% The tumorelimination vs. immune escape
    subplot(2,3,6);
    load('Data/Data_rho_gamma.mat');  
    [X,Y] = meshgrid(gamma0,rho);
    mymap = [0.75 0.75 0.75;0.45 0.45 0.45];
    pcolor(X,Y,C');
    colormap(mymap); 
    shading interp;  xlim([0, 2e-5]);  caxis([0 1e+5]); 
    xlabel('\gamma','fontname','Times','FontSize', FS_default);  
    ylabel('\rho ','fontname','Times','FontSize', FS_default); 
    text(0.5e-5,0.06,'Elimination','FontSize',FS_inside);
    text(0.15e-5,0.17,'Escape','FontSize',FS_inside);
    text(-0.4e-5,0.205,'\bf{f}','FontSize',FS_lab);
    xlim([0,1.6*1e-5]);
    
    
    %% Export figure
    set(hf1,'Units','centimeters');
    screenposition = get(gcf,'Position');
    set(hf1,'PaperPosition',[0 0 1.2*screenposition(3) 0.9*screenposition(4)],...
            'PaperSize',[1.2*screenposition(3) 0.9*screenposition(4)]);
        
   print(hf1,'-r700','-dpdf','Fig/Fig7.pdf');   