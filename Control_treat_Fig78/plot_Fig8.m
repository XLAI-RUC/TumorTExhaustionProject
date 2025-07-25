clc; close all; clear all;
    
    FS_inside = 9;   % FontSize inside figure
    FS_default = 9;  % FontSize default
    FS_lab = 11;      % FontSize of (A) (B) (C)
    FS = 11;
    
    set(0,'defaultaxesfontsize',FS_default,'defaultlinelinewidth',1.5);
    colorM = [0.9 0 0; 0.4660 0.6740 0.1880;0 0.4470 0.7410; 0.7 0.5 0; 0 0 1; 0.7350 0.1780 0.1840];
    hf1 = figure;
    
    subplot(2,2,1);
    R = load('Data/fit_INFalphaAntiPDL1_sim.dat'); 
    point = load('Data/fit_INFalphaAntiPDL1_exp.dat'); 
    time = R(:,1); C = R(:,2);  x2 = point(:,1);  y2 = point(:,2); 
    hs1 = scatter(x2, y2, 'MarkerFaceColor', colorM(3,:),'MarkerEdgeColor', colorM(3,:),'Sizedata', 8); 
    hold on
    hs2 = plot(time, C,'Color', colorM(1,:)); box on;
    text(3,3.2e+8,'IFN\alpha-anti-PD-L1','Interpreter','tex','FontSize',FS_inside);
    xlabel('Time (days)','Interpreter','tex','FontSize',12, 'fontname', 'Times');  
    ylabel('Tumor cells (C)', 'FontSize',12, 'fontname', 'Times'); 
    text(-6,3.65e+8,'\bf{a}','FontSize',FS_lab);
    ylim([0 3.5e+8]);
    leg = legend([hs1,hs2], {'Experiment','Simulation'});  legend boxoff;
    leg.ItemTokenSize = [14,14];
    pos = get(leg, 'Position');
    set(leg, 'Position', [pos(1)*0.8 pos(2)*0.95 0 0]);
    
    subplot(2,2,2);
    R = load('Data/fit_PD1IL2v_sim.dat'); 
    point = load('Data/fit_PD1IL2v_exp.dat'); 
    time = R(:,1); C = R(:,2); x2 = point(:,1); y2 = point(:,2);
    hs1 = plot(x2, y2, 'o-', 'Color', colorM(3,:), 'LineWidth', 1.5, 'Marker', 'o', 'MarkerEdgeColor', colorM(3,:), 'MarkerFaceColor', colorM(3,:), 'MarkerSize', 3);
    hold on
    hs2 = plot(time, C,'Color', colorM(1,:)); 
    text(1.5,1.66e+8,'PD1-IL2v','Interpreter','tex','FontSize',FS_inside);
    xlabel('Time (days)','Interpreter','tex','FontSize',12, 'fontname', 'Times');  
    ylabel('Tumor cells (C)', 'FontSize',12, 'fontname', 'Times'); 
    text(-4.7,1.87e+8,'\bf{b}','FontSize',FS_lab);  ylim([0 1.8e+8]);
    leg = legend([hs1,hs2], {'Experiment','Simulation'});  legend boxoff;
    leg.ItemTokenSize = [14,14];
    pos = get(leg, 'Position');
    set(leg, 'Position', [pos(1)*1.1 pos(2)*0.8 0 0]);
        
    subplot(2,2,3);
    R = load('Data/treatment_PD1IL2v_AntiPDL1.dat'); 
    time = R(:,1); C = R(:,2);
    hs1=plot(time, C,'-d','MarkerIndices',1:20:length(C),'Color', colorM(1,:),'MarkerSize',3); 
    hold on
    
    R = load('Data/treatment_INFalphaAntiPDL1_25ug.dat'); 
    time = R(:,1); C = R(:,2);
    hs2 = plot(time, C,'-h','MarkerIndices',1:20:length(C),'Color', colorM(3,:),'MarkerSize',3); 
    
    R = load('Data/treatment_AntiPDL1_25ug.dat'); 
    time = R(:,1); C = R(:,2);
    hs3 = plot(time, C,'-o','MarkerIndices',1:20:length(C),'Color', colorM(2,:),'MarkerSize',3); 
    
    R=load('Data/treatment_PD1IL2v_25ug.dat'); 
    time=R(:,1);C = R(:,2);
    hs4 = plot(time, C,'-^','MarkerIndices',1:20:length(C),'Color', colorM(4,:),'MarkerSize',3); 
    R = load('Data/Control_simulation.dat'); 
    time = R(:,1); C = R(:,2);
    hs5 = plot(time, C,'-*','MarkerIndices',1:5:length(C),'Color', colorM(5,:),'MarkerSize',3); 
    
    leg = legend([hs5,hs3,hs2,hs4,hs1], {'control','Anti-PD-L1','IFN{\alpha}-anti-PD-L1','PD1-IL2v','Anti-PD-L1 + PD1-IL2v'});  
    leg.ItemTokenSize = [18,18];
    pos = get(leg, 'Position');
    set(leg, 'Position', [pos(1) pos(2)*0.8 0.4 0]);
    legend boxoff;
    xlabel('Time (days)','Interpreter','tex','FontSize',12, 'fontname', 'Times');  
    ylabel('Tumor cells (C)', 'FontSize',12, 'fontname', 'Times'); 
    text(-4.8,7.6e+7,'\bf{c}','FontSize',FS_lab);
    ylim([0 7.5e+7]); xlim([0 40])
        
    subplot(2,2,4);
    load('Data/data_gamma_zeta'); 
    [zeta_grid, gamma_grid] = ndgrid(zeta, gamma);
    pcolor(zeta_grid* 1e-6,gamma_grid* 1e-6, C);  
    shading flat;  
    mymap = [0.75 0.75 0.75;0.45 0.45 0.45];
    colormap(mymap); 
    set(gca, 'CLim',[0,1]);  
    ylabel('\gamma (Anti-PD-L1)','fontname','Times','FontSize', 12);  
    xlabel('\gamma_B (PD1-IL2v)','fontname','Times','FontSize', 12);  
    xlims = xlim;  
    ylims = ylim;    
    text((xlims(1) + (xlims(2) - xlims(1)) * 0.4), ylims(2) - (ylims(2) - ylims(1)) * 0.4, '\bf{Elimination}', 'FontSize', FS_inside);  
    text(xlims(1) + (xlims(2) - xlims(1))* 0.1, ylims(2) - (ylims(2) - ylims(1)) * 0.8, '\bf{Escape}', 'FontSize', FS_inside); 
    text(-0.07e-5,7e-6,'\bf{d}','FontSize',FS_lab);

    print(hf1,'-r700','-dpdf','Fig/Fig8.pdf');   