clc; close all; clear all;
    
    FS_inside = 7;   % FontSize inside figure
    FS_default = 10;  % FontSize default
    FS_lab = 11;      % FontSize of (A) (B) (C)

    set(0,'defaultaxesfontsize',FS_default,'defaultlinelinewidth',0.6);
     k = 8;
     

    %%  %% OS curves 
    load("Model_LHS_PD1IL2v_AntiPDl1_10_20_30_40.mat")
     
    
    %% Anti-PD-L1 Monotherapy 
    p1 = prccC_1(1,:); p1_t = prccT1_1(1,:);
    p2 = prccC_1(2,:); p2_t = prccT1_1(2,:);
    p3 = prccC_1(3,:); p3_t = prccT1_1(3,:);
    p4 = prccC_1(4,:); p4_t = prccT1_1(4,:);
    
    figure
    subplot(2,1,1)
    data0 = [p1',p2',p3',p4'];
    [data ind] = sort(data0,'descend');
    PRCC_var1 = PRCC_var(ind(:,1));
    categories = {'Group A', 'Group B', 'Group C'};
    b = bar(data, 'grouped'); % 'grouped'为并列显示
    set(gca,'XTickLabel',PRCC_var1,'XTick',(1:k));
    ylabel('PRCC');
    text(-0.2,1,'\bf{a}','FontSize',FS_lab);  
    text(4,0.75,'Tumor cells'); 
       
    subplot(2,1,2)
    data0 = [p1_t',p2_t',p3_t',p4_t'];
    [data ind] = sort(data0);
    PRCC_var1 = PRCC_var(ind(:,1));
    categories = {'Group A', 'Group B', 'Group C'};
    b = bar(data, 'grouped'); % 'grouped'为并列显示
    set(gca,'XTickLabel',PRCC_var1,'XTick',(1:k));
    ylabel('PRCC'); 
    text(-0.2,1,'\bf{b}','FontSize',FS_lab);
    text(4,0.7,'T_1 cells'); 
    
    %% PD1IL2v Monotherapy
    p1 = prccC_2(1,:); p1_t = prccT1_2(1,:);
    p2 = prccC_2(2,:); p2_t = prccT1_2(2,:);
    p3 = prccC_2(3,:); p3_t = prccT1_2(3,:);
    p4 = prccC_2(4,:); p4_t = prccT1_2(4,:);
    
    figure
    subplot(2,1,1)
    data0 = [p1',p2',p3',p4'];
    [data ind] = sort(data0,'descend');
    PRCC_var1 = PRCC_var(ind(:,1));
    categories = {'Group A', 'Group B', 'Group C'};
    b = bar(data, 'grouped'); % 'grouped'为并列显示
    set(gca,'XTickLabel',PRCC_var1,'XTick',(1:k));
    ylabel('PRCC');
    text(-0.2,1,'\bf{a}','FontSize',FS_lab);  
    text(4,0.75,'Tumor cells'); 
       
    subplot(2,1,2)
    data0 = [p1_t',p2_t',p3_t',p4_t'];
    [data ind] = sort(data0);
    PRCC_var1 = PRCC_var(ind(:,1));
    categories = {'Group A', 'Group B', 'Group C'};
    b = bar(data, 'grouped'); % 'grouped'为并列显示
    set(gca,'XTickLabel',PRCC_var1,'XTick',(1:k));
    ylabel('PRCC'); 
    text(-0.2,1,'\bf{b}','FontSize',FS_lab);
    text(4,0.7,'T_1 cells'); 
    
   %% PD1IL2v + Anti-PD-L1
    p1 = prccC_3(1,:); p1_t = prccT1_3(1,:);
    p2 = prccC_3(2,:); p2_t = prccT1_3(2,:);
    p3 = prccC_3(3,:); p3_t = prccT1_3(3,:);
    p4 = prccC_3(4,:); p4_t = prccT1_3(4,:);
    
    figure
    subplot(2,1,1)
    data0 = [p1',p2',p3',p4'];
    [data ind] = sort(data0,'descend');
    PRCC_var1 = PRCC_var(ind(:,1));
    categories = {'Group A', 'Group B', 'Group C'};
    b = bar(data, 'grouped'); % 'grouped'为并列显示
    set(gca,'XTickLabel',PRCC_var1,'XTick',(1:k));
    ylabel('PRCC');
    text(-0.2,1,'\bf{a}','FontSize',FS_lab);  
    text(4,0.75,'Tumor cells'); 
       
    subplot(2,1,2)
    data0 = [p1_t',p2_t',p3_t',p4_t'];
    [data ind] = sort(data0);
    PRCC_var1 = PRCC_var(ind(:,1));
    categories = {'Group A', 'Group B', 'Group C'};
    b = bar(data, 'grouped'); % 'grouped'为并列显示
    set(gca,'XTickLabel',PRCC_var1,'XTick',(1:k));
    ylabel('PRCC'); 
    text(-0.2,1,'\bf{b}','FontSize',FS_lab);
    text(4,0.7,'T_1 cells'); 
    

%     print(hf1,'-r700','-dpdf','Fig_LHS_PD1IL2v_AntiPDl1.pdf');  