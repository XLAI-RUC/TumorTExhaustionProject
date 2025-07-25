function PhasePlane(C,T_1,ax)
clc;

colorM = [0.9 0 0; 0.4660 0.6740 0.1880;0 0.4470 0.7410; 0.7 0.5 0; 0 0 1; 0.7350 0.1780 0.1840];

hold on
m = size(C); % ����C����������ɵľ���  
for i = 1:m(2) 
   rowC = C(:, i);   
   rowT_1 = T_1(:, i); 
   plot(rowC, rowT_1, 'Color',colorM(3,:)); 
   h1= scatter(rowC(1), rowT_1(1),5, 'filled', 'k'); 
   h2 =scatter(rowC(m(1)), rowT_1(m(1)),18, 'filled', 'r'); 
end  
    % ... ������������xlabel, ylabel, legend��Ҳ��Ҫʹ��ax��  
    xlabel(ax, 'Tumor cells (C)','fontname','Times');  
    ylabel(ax, 'PD-1^{lo}CD8^+ T cells (T_1)','fontname','Times');  
    legendHandles = [h1, h2, plot([],[],'Parent',ax)]; % ָ��ParentΪax  
    legendStrings = {'initial values', 'equilibrium', 'Line'};  
    legend(ax, legendHandles, legendStrings, 'Location', 'best','FontSize', 4);  

end