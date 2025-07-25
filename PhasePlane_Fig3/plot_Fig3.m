    clc; close all; clear;
    
    FS_inside = 9;   % FontSize inside figure
    FS_default = 7;  % FontSize default
    FS_lab = 11;      % FontSize of (A) (B) (C)
    FS=9;

    set(0,'defaultaxesfontsize',FS_default,'defaultlinelinewidth',0.8);
    colorM = [0.9 0 0; 0.4660 0.6740 0.1880;0 0.4470 0.7410; 0.7 0.5 0; 0 0 1; 0.7350 0.1780 0.1840];
    arrowh_size_axis = [190 90]+20;
    arrowh_location_axis = 100;
    
    hf1 = figure;
    
    subplot(2,3,1);  
    load('output/data_exist2.mat');  
    PhasePlane(C, T_1, gca); legend off;
    x1 = zeros(5,4);  z1 = [12,10,12,3,1,5]; z2 = z1+5;
    for k = 1: size(C,2)
        rowC = C(:, k);   
        rowT_1 = T_1(:, k);
        x1(k,1:2) = rowC([z1(k),z2(k)]);
        x1(k,3:4) = rowT_1([z1(k),z2(k)]);
    end
    for k = 1:size(C,2)
       arrowh(x1(k,1:2),x1(k,3:4),colorM(3,:),arrowh_size_axis,arrowh_location_axis);
    end
    text(-4.2e+8,5.1e+7,'\bf{a}','FontSize',FS_lab);
    text(0.9e+9,3.3e+7,'{\it S}_{11}','Interpreter','tex','FontSize',FS);
    
    subplot(2,3,2);  
    load('output/data_exist1.mat');  
    PhasePlane(C, T_1, gca); legend off;
    x1 = zeros(5,4);  z1 = [12,10,5,5,10,1]; z2 = z1+5;
    for k = 1: size(C,2)
        rowC = C(:, k);   
        rowT_1 = T_1(:, k);
        x1(k,1:2) = rowC([z1(k),z2(k)]);
        x1(k,3:4) = rowT_1([z1(k),z2(k)]);
    end
    for k = 1:size(C,2)
       arrowh(x1(k,1:2),x1(k,3:4),colorM(3,:),arrowh_size_axis,arrowh_location_axis);
    end
    text(-7.8e+8,10.2e+6,'\bf{b}','FontSize',FS_lab);
    text(2e+9,7.9e+6,'{\it S}_{11}','Interpreter','tex','FontSize',FS);
    
    subplot(2,3,3);  
    load('output/data_exist3.mat');  
    PhasePlane(C, T_1, gca); legend off;
    x1 = zeros(5,4);  z1 = [10,5,3,1,1,5]; z2 = z1+2;
    for k = 1: size(C,2)
        rowC = C(:, k);   
        rowT_1 = T_1(:, k);
        x1(k,1:2) = rowC([z1(k),z2(k)]);
        x1(k,3:4) = rowT_1([z1(k),z2(k)]);
    end
    for k = 1:size(C,2)
       arrowh(x1(k,1:2),x1(k,3:4),colorM(3,:),arrowh_size_axis,arrowh_location_axis);
    end
    text(-4.2e+8,5.1e+7,'\bf{c}','FontSize',FS_lab);
    text(1.8e+9,1.9e+7,'{\it S}_{11}','Interpreter','tex','FontSize',FS);
    
    subplot(2,3,4); 
    load('output/data_exist.mat');  
    PhasePlane(C, T_1, gca);
    legend off;
    x1 = zeros(5,4);  z1 = [13,10,10,13,12,12,1*ones(1,8)]; z2 = z1+5;
    for k = 1: size(C,2)
        rowC = C(:, k);   
        rowT_1 = T_1(:, k);
        x1(k,1:2) = rowC([z1(k),z2(k)]);
        x1(k,3:4) = rowT_1([z1(k),z2(k)]);
    end
    for k = 1:size(C,2)
       arrowh(x1(k,1:2),x1(k,3:4),colorM(3,:),arrowh_size_axis,arrowh_location_axis);
    end
    plot(2.05e+8,4.3e+7,'rs','MarkerFaceColor','r','MarkerSize',4)
    text(-2.5e+8,6.7e+7,'\bf{d}','FontSize',FS_lab);
    text(1.4e+8,3.85e+7,'{\it S}_{11}','Interpreter','tex','FontSize',8,'Color','r');
    text(11.5e+8,2.2e+7,'{\it S}_{12}','Interpreter','tex','FontSize',FS,'Color','r');
    text(-1e+8,5.15e+7,'{\it S}_{02}','Interpreter','tex','FontSize',FS,'Color','r');
     ylim([0.9e+7 6.5e+7]);

    subplot(2,3,5);  
    load('output/data_exist5.mat');  
    PhasePlane(C, T_1, gca);  legend off;
    x1 = zeros(5,4);  z1 = [5,5,1,1,1,5]; z2 = z1+2;
    for k = 1: size(C,2)
        rowC = C(:, k);   
        rowT_1 = T_1(:, k);
        x1(k,1:2) = rowC([z1(k),z2(k)]);
        x1(k,3:4) = rowT_1([z1(k),z2(k)]);
    end
    for k = 1:size(C,2)
       arrowh(x1(k,1:2),x1(k,3:4),colorM(3,:),arrowh_size_axis,arrowh_location_axis);
    end
    text(-4e+8,6.2e+7,'\bf{e}','FontSize',FS_lab);
    text(0.3e+8,5.4e+7,'{\it S}_{02}','Interpreter','tex','FontSize',FS);

    subplot(2,3,6);  
    load('output/data_exist6.mat');  
    PhasePlane(C, T_1, gca);  legend off;
    x1 = zeros(5,4);  z1 = [1,2,1,1,1,1]; z2 = z1+2;
    for k = 1: size(C,2)
        rowC = C(:, k);   
        rowT_1 = T_1(:, k);
        x1(k,1:2) = rowC([z1(k),z2(k)]);
        x1(k,3:4) = rowT_1([z1(k),z2(k)]);
    end
    for k = 1:size(C,2)
       arrowh(x1(k,1:2),x1(k,3:4),colorM(3,:),arrowh_size_axis,arrowh_location_axis);
    end
    ylim([0 10e+6])
    text(-2.5e+6,10.1e+6,'\bf{f}','FontSize',FS_lab);
    text(0e+6,0.7e+6,'{\it S}_{01}','Interpreter','tex','FontSize',FS);
    
    
    set(hf1,'Units','centimeters');
    screenposition = get(gcf,'Position');
     set(hf1,'PaperPosition',[0 0 1.2*screenposition(3) 0.95*screenposition(4)],...
            'PaperSize',[1.2*screenposition(3) 0.95*screenposition(4)]);
    print(hf1,'-r700','-dpdf','Fig/Fig3.pdf');   