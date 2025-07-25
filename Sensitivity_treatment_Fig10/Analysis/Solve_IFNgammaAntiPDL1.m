%% Modified from MBS 2025 (Ma,Lai)
%% Step 1: plot the distribution of all parameters, that is, Fig1
%% Step 2: plot the heat plot of some parameters,   that is, Fig2
%% ========================================================================
Num = 20000;
FontS = 7; % Font Size 
LineW = 1;  % Line Width

set(0,'defaultaxesfontsize',FontS ,'defaultlinelinewidth',LineW );
%% read simulation data 
load("../SA_LHS_INFalphaAntiPDL1_25ug/Data/Model_LHS_INFalphaAntiPDL1_25ug_10_20_30_40.mat");
  Y = LHSmatrix;  % LHSmatrix is the LHS parameter set

%% ------------------------------------------------------------------------
Control  = zeros(Num,1);
UnControl = zeros(Num,1);

Control  = (C_1(4,:)<y0(1)); 
sum(Control)
UnControl = (C_1(4,:)>5*y0(1)); % C_1(4,:) is tumor cell count at day 60
sum(UnControl)
h1 = figure;
xlabels = { '\alpha_1','\lambda_T','\rho','d_1','\lambda_C','\sigma','\epsilon','\beta_2'};  

for i = 1:8
    subplot(3,3,i);  
    
    param = Y(:, i);
    death_param = param(UnControl == 1);
    [f_death, x_death] = ksdensity(death_param);  %  probability density estimate, f_death
    fill(x_death, f_death, [0.8500, 0.3250, 0.0980], ...
        'FaceAlpha', 0.3, 'EdgeColor', [0.8500, 0.3250, 0.0980], 'LineWidth', 2); hold on;

    live_param = param(Control == 1);
    [f_live, x_live] = ksdensity(live_param);
    fill(x_live, f_live, [0.3010, 0.7450, 0.9330], ...
        'FaceAlpha', 0.3, 'EdgeColor', [0.3010, 0.7450, 0.9330], 'LineWidth', LineW);

    [~, pval] = kstest2(death_param, live_param);  % A test decision for the null hypothesis that the data in vectors x1 and x2 are from the same continuous distribution, using the two-sample Kolmogorov-Smirnov test

    xlabel(['$', xlabels{i}, '$'], 'Interpreter', 'latex', 'FontSize', FontS);
    
    text(0.05, 0.85, ['p = ', num2str(pval, '%.3g')], ...
        'Units', 'normalized', 'FontSize', FontS, 'Color', 'k');

    if pval < 0.05
        text(0.05, 0.70, '*', 'Units', 'normalized', ...
            'FontSize', FontS, 'Color', 'r', 'FontWeight', 'bold');
    end
end

% save the figure as pdf
set(gcf,'Units','centimeters');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 1.0*screenposition(3) 1.2*screenposition(4)],...
           'PaperSize',[1.0*screenposition(3) 1.2*screenposition(4)]);

 print(h1,'-dpdf','Fig_INFgammaAntiPDL1_dist.pdf');

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %% ---------------------------------------------------------------------- 
selected_idx = [1:6,8];
param_labels = { '\alpha_1','\lambda_T','\rho','d_1','\lambda_C','\sigma','\beta_2'};
%param_labels = { '\alpha_1','\lambda_T','\rho','d_1','\lambda_C','\sigma','\epsilon','\beta_2'};  

% Just consider selected paramters
Y_sub = Y(:, selected_idx);

h1  = figure;
set(0, 'defaultaxesfontsize', FontS, 'defaultlinelinewidth', 1.0);

n = length(selected_idx);  

% %hw = tight_subplot(n,n,[.1 .03],[.05 .05],[.05 .01]);
% hw = tight_subplot(n,n,[.03 .03],[.03 .03],[.03 .03]);
% %%% (Nh, Nw, gap[gap_h gap_w], marg_h[lower upper], marg_w[left right] )

for i = 1:n
    for j = 1:n
        subplot(n, n, (i-1)*n + j);
        if i == j
            % === KDE  ===
            param = Y_sub(:, i);

            % Tumor-control
            dp = param(UnControl == 1);
            [f_d, x_d] = ksdensity(dp);
            fill(x_d, f_d, [0.8500, 0.3250, 0.0980], ...
                'FaceAlpha', 0.3, 'EdgeColor', [0.8500, 0.3250, 0.0980], 'LineWidth', LineW); hold on;

            % Immune-escape
            lp = param(Control == 1);
            [f_l, x_l] = ksdensity(lp);
            fill(x_l, f_l, [0.3010, 0.7450, 0.9330], ...
                'FaceAlpha', 0.3, 'EdgeColor', [0.3010, 0.7450, 0.9330], 'LineWidth', LineW);

            axis tight; axis square;

            % KS
            [~, pval] = kstest2(dp, lp);
            
%             % Show p-value
%             text(0.05, 0.85, ['p = ', num2str(pval, '%.3g')], ...
%                 'Units', 'normalized', 'FontSize', 14, 'Color', 'k');
%         
%             % If p < 0.05, mark "*".
%             if pval < 0.05
%                 text(0.05, 0.70, '*', 'Units', 'normalized', ...
%                     'FontSize', 24, 'Color', 'r', 'FontWeight', 'bold');
%             end
        if i>1 && i<n
           set(gca, 'XTick', [], 'YTick', []);
           xlabel(['$', param_labels{j}, '$'], 'Interpreter', 'latex');
        end  
        if i==1
           set(gca,'XTick', []);         
           ylabel(['$', param_labels{i}, '$'], 'Interpreter', 'latex');
           xlabel(['$', param_labels{j}, '$'], 'Interpreter', 'latex');
           
        end  
        if i==n && j==7 
              set(gca,'YTick', []); 
              xt = get(gca, 'XTick');
              set(gca, 'XTickLabel', xt/1e-8); 
              xlabel(['$', param_labels{j},'(\times 10^{-8})', '$'], 'Interpreter', 'latex');
        end
 
        elseif j < i  %%% Lower triangular 
            % === Tumor-control
            x = Y_sub(Control == 1, j);
            y = Y_sub(Control == 1, i);

            xi = linspace(min(x), max(x), 100);
            yi = linspace(min(y), max(y), 100);
            [X, Yg] = meshgrid(xi, yi);
            [f, ~] = ksdensity([x, y], [X(:), Yg(:)]);
            F = reshape(f, size(X));

            imagesc(xi, yi, F);
            axis xy square;
            xlim([min(xi), max(xi)]);
            ylim([min(yi), max(yi)]);
            
            % Set ticks and x-y labels
            if j==1 && i<n
               set(gca,'XTick', []); 
               if i==2 || i==3 || i==5
                   yt = get(gca, 'YTick');
                   yticks([min(yt),(min(yt)+max(yt))/2,max(yt)]);
               end  
               
               ylabel(['$', param_labels{i}, '$'], 'Interpreter', 'latex');
            end  
            if i==n && j>1 && j ~=6 && j~=7
                set(gca,'YTick', []); 
                xlabel(['$', param_labels{j}, '$'], 'Interpreter', 'latex');
            end 
            if j>1 && i<n
                set(gca,'XTick', [],'YTick', []); 
            end
            if i==n && j==1
               %set(gca,'YTick', []); 
               xt = get(gca, 'XTick');
               set(gca, 'XTickLabel', xt/1e-10); 
               xlabel(['$', param_labels{j},'(\times 10^{-10})', '$'], 'Interpreter', 'latex');
               ylabel(['$', param_labels{i}, '$'], 'Interpreter', 'latex');
            end
            if i==n && j==6
               set(gca,'YTick', []); 
               xt = get(gca, 'XTick');
               set(gca, 'XTickLabel', xt/1e-17); 
               xlabel(['$', param_labels{j},'(\times 10^{-17})', '$'], 'Interpreter', 'latex');
            end 
                           
        elseif j > i   %%% Upper triangular
            % === Immune-escape
            x = Y_sub(UnControl == 1, j);
            y = Y_sub(UnControl == 1, i);

            xi = linspace(min(x), max(x), 100);
            yi = linspace(min(y), max(y), 100);
            [X, Yg] = meshgrid(xi, yi);
            [f, ~] = ksdensity([x, y], [X(:), Yg(:)]);
            F = reshape(f, size(X));

            imagesc(xi, yi, F);
            axis xy square;
            xlim([min(xi), max(xi)]);
            ylim([min(yi), max(yi)]);
            
            set(gca,'XTick', [],'YTick', []); 
        end
    end
end
%% ------------------------------------------------------------------------

legend_center_y = 0.03;
box_w = 0.015;
box_h = 0.015;

y_box = legend_center_y - box_h/2;
y_text = y_box;

% --- Tumor control  ---
annotation('rectangle', [0.18, y_box, box_w, box_h], ...
    'FaceColor', [0.3010, 0.7450, 0.9330], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
annotation('textbox', [0.18 + box_w + 0.005, y_text, 0.1, box_h], ...
    'String', 'Tumor control', 'FontSize', FontS, ...
    'EdgeColor', 'none', 'VerticalAlignment', 'middle');

% --- Immune escape ---
annotation('rectangle', [0.34, y_box, box_w, box_h], ...
    'FaceColor', [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
annotation('textbox', [0.34 + box_w + 0.005, y_text, 0.15, box_h], ...
    'String', 'Immune escape', 'FontSize', FontS, ...
    'EdgeColor', 'none', 'VerticalAlignment', 'middle');

% --- colorbar ---
cmap = parula(256);
cbar_img = reshape(cmap, [1, size(cmap,1), 3]);

cbar_x = 0.665;
cbar_y = 0.02;
cbar_w = 0.18;
cbar_h = 0.025;
cbar_center_y = cbar_y + cbar_h / 2;

axes('Position', [cbar_x, cbar_y, cbar_w, cbar_h]);
image(cbar_img);
axis off;

% Low/High  colorbar 
y_textbox = cbar_center_y - box_h / 2;

annotation('textbox', [0.63, y_box, 0.03, box_h], ...
    'String', 'Low', 'EdgeColor', 'none', ...
    'FontSize', FontS, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

annotation('textbox', [0.85, y_box, 0.05, box_h], ...
    'String', 'High', 'EdgeColor', 'none', ...
    'FontSize', FontS, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% % %% ------------------------------------------------------------------------
% save the figure as pdf
set(gcf,'Units','centimeters');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 1*screenposition(3) 1.45*screenposition(4)],...
           'PaperSize',[1*screenposition(3) 1.45*screenposition(4)]);

print(h1,'-dpdf','Fig_INFgammaAntiPDL1.pdf');