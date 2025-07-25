clc; close all; clear all;
    
    FS_inside = 7;   % FontSize inside figure
    FS_default = 10;  % FontSize default
    FS_lab = 11;      % FontSize of (A) (B) (C)

    set(0,'defaultaxesfontsize',FS_default,'defaultlinelinewidth',0.6);
    
    hf1 = figure;
    %%  %% OS curves 
    load Model_LHS_30_60_120.mat
    p1 = prcc11(1,:); p1_t = prcc12(1,:);
    p2 = prcc21(1,:); p2_t = prcc22(1,:);
    p3 = prcc31(1,:); p3_t = prcc32(1,:);
    
    subplot(2,1,1)
    data0 = [p1',p2',p3'];
    [data ind] = sort(data0,'descend');
    PRCC_var1 = PRCC_var(ind(:,1));
    categories = {'Group A', 'Group B', 'Group C'};
    b = bar(data, 'grouped'); % 'grouped'为并列显示
    set(gca,'XTickLabel',PRCC_var1,'XTick',(1:k));
    ylabel('PRCC');
    text(-0.2,1,'\bf{a}','FontSize',FS_lab);  
    text(4,0.75,'Tumor cells'); 
       
    subplot(2,1,2)
    data0 = [p1_t',p2_t',p3_t'];
    [data ind] = sort(data0);
    PRCC_var1 = PRCC_var(ind(:,1));
    categories = {'Group A', 'Group B', 'Group C'};
    b = bar(data, 'grouped'); % 'grouped'为并列显示
    set(gca,'XTickLabel',PRCC_var1,'XTick',(1:k));
    ylabel('PRCC'); 
    text(-0.2,1,'\bf{b}','FontSize',FS_lab);
    text(4,0.7,'T_1 cells'); 
    
    set(hf1,'Units','centimeters');
    screenposition = get(gcf,'Position');
    set(hf1,'PaperPosition',[0 0 1*screenposition(3) 0.8*screenposition(4)],...
    'PaperSize',[1*screenposition(3) 0.8*screenposition(4)]);
    print(hf1,'-r700','-dpdf','Fig9.pdf');  