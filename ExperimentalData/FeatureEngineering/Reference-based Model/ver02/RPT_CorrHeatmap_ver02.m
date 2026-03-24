%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Correlation Heatmap (ver02)
% - Pearson Correlation between Features (X) and Labels (y)
% - 3 Heatmaps: Charge-only, Discharge-only, Combined
% - Style: p-value stars, in-cell numbers, academic colormap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ========================================================================
% Section 1: Configuration & Data Loading
% ========================================================================
currentScriptPath = mfilename('fullpath');
[ver02Dir, ~, ~] = fileparts(currentScriptPath);
saveDir = fullfile(ver02Dir, 'FeatureLabelVisual');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

dataDir = fullfile(ver02Dir, 'RPT_FeatureLabelExtractor');
load(fullfile(dataDir, 'Feature_Matrix_ver02.mat'), 'FeatureTable_ver02');
load(fullfile(dataDir, 'Label_Matrix_ver02.mat'), 'LabelTable_ver02');

%% ========================================================================
% Section 2: Build Named Columns
% ========================================================================
% --- Charge Features ---
X_chg_names = {'dQ_chg_S1','dQ_chg_S2','dQ_chg_S3','dQ_chg_S4','dQ_chg_S5', 'Peak_H_chg', 'Peak_A_chg', 'Peak_Pos_chg'};
X_chg = table2array(FeatureTable_ver02(:, X_chg_names));

% --- Discharge Features ---
X_dch_names = {'dQ_dch_S1','dQ_dch_S2','dQ_dch_S3','dQ_dch_S4','dQ_dch_S5', 'Peak_H_dch', 'Peak_A_dch', 'Peak_Pos_dch', 'Energy_dch'};
X_dch = table2array(FeatureTable_ver02(:, X_dch_names));

% --- Common Features ---
X_common_names = {'C_eff_chg', 'C_eff_dch'};
X_common = table2array(FeatureTable_ver02(:, X_common_names));

% --- Target Labels ---
y_names = {'SOH', 'LLI', 'LAM'};
y_targets = table2array(LabelTable_ver02(:, y_names));

%% ========================================================================
% Section 3: Academic Colormap (Blue → White → Red)
% ========================================================================
n = 128;
blue = [0.15 0.25 0.55];
white = [1 1 1];
red = [0.7 0.15 0.15];
c1 = [linspace(blue(1),white(1),n)', linspace(blue(2),white(2),n)', linspace(blue(3),white(3),n)'];
c2 = [linspace(white(1),red(1),n)', linspace(white(2),red(2),n)', linspace(white(3),red(3),n)'];
cmap = [c1; c2];

%% ========================================================================
% Section 4: Plot Heatmaps
% ========================================================================

% ---- [1] Charge Heatmap (Square: Features + Labels) ----
all_chg = [X_chg, X_common, y_targets];
all_chg_names = [X_chg_names, X_common_names, y_names];
n_feat_chg = length(X_chg_names) + length(X_common_names);

[R_chg, P_chg] = corr(all_chg, 'Rows', 'pairwise');
n_chg = length(all_chg_names);

fig1 = figure('Position', [50, 50, 1200, 1200], 'Color', 'w', 'Visible', 'on');
imagesc(R_chg); colormap(cmap); caxis([-1 1]);
cb = colorbar('Location','eastoutside');
cb.Label.String = 'Pearson r'; cb.Label.FontSize = 11; cb.Label.FontWeight = 'bold';
axis square;
set(gca, 'XTick',1:n_chg, 'XTickLabel',all_chg_names, 'FontSize',12, 'FontWeight','bold');
set(gca, 'YTick',1:n_chg, 'YTickLabel',all_chg_names, 'FontSize',12, 'FontWeight','bold');
set(gca, 'XAxisLocation','top'); xtickangle(45);
hold on;
line([0.5,n_chg+0.5],[0.5,n_chg+0.5],'Color',[0.3 0.3 0.3],'LineWidth',1.5);
line([n_feat_chg+0.5,n_feat_chg+0.5],[0.5,n_chg+0.5],'Color','k','LineWidth',2);
line([0.5,n_chg+0.5],[n_feat_chg+0.5,n_feat_chg+0.5],'Color','k','LineWidth',2);

for i=1:n_chg
    for j=1:n_chg
        if i==j, continue; end
        val=R_chg(i,j); pval=P_chg(i,j);
        if pval<0.001, star='***'; elseif pval<0.01, star='**'; elseif pval<0.05, star='*'; else, star=''; end
        if abs(val)>0.6, clr='w'; else, clr='k'; end
        text(j,i,sprintf('%.2f%s',val,star),'HorizontalAlignment','center','FontSize',10,'FontWeight','bold','Color',clr);
    end
end
hold off;
title({'Charge Features & Labels Correlation','(*p<0.05, **p<0.01, ***p<0.001)'},'FontSize',16,'FontWeight','bold');
saveas(fig1,fullfile(saveDir,'CorrHeatmap_Charge.fig'));
close(fig1);
fprintf('Charge Heatmap saved.\n');

% ---- [2] Discharge Heatmap (Square: Features + Labels) ----
all_dch = [X_dch, X_common, y_targets];
all_dch_names = [X_dch_names, X_common_names, y_names];
n_feat_dch = length(X_dch_names) + length(X_common_names);

[R_dch, P_dch] = corr(all_dch, 'Rows', 'pairwise');
n_dch = length(all_dch_names);

fig2 = figure('Position', [100, 50, 1200, 1200], 'Color', 'w', 'Visible', 'on');
imagesc(R_dch); colormap(cmap); caxis([-1 1]);
cb = colorbar('Location','eastoutside');
cb.Label.String = 'Pearson r'; cb.Label.FontSize = 11; cb.Label.FontWeight = 'bold';
axis square;
set(gca, 'XTick',1:n_dch, 'XTickLabel',all_dch_names, 'FontSize',12, 'FontWeight','bold');
set(gca, 'YTick',1:n_dch, 'YTickLabel',all_dch_names, 'FontSize',12, 'FontWeight','bold');
set(gca, 'XAxisLocation','top'); xtickangle(45);
hold on;
line([0.5,n_dch+0.5],[0.5,n_dch+0.5],'Color',[0.3 0.3 0.3],'LineWidth',1.5);
line([n_feat_dch+0.5,n_feat_dch+0.5],[0.5,n_dch+0.5],'Color','k','LineWidth',2);
line([0.5,n_dch+0.5],[n_feat_dch+0.5,n_feat_dch+0.5],'Color','k','LineWidth',2);

for i=1:n_dch
    for j=1:n_dch
        if i==j, continue; end
        val=R_dch(i,j); pval=P_dch(i,j);
        if pval<0.001, star='***'; elseif pval<0.01, star='**'; elseif pval<0.05, star='*'; else, star=''; end
        if abs(val)>0.6, clr='w'; else, clr='k'; end
        text(j,i,sprintf('%.2f%s',val,star),'HorizontalAlignment','center','FontSize',10,'FontWeight','bold','Color',clr);
    end
end
hold off;
title({'Discharge Features & Labels Correlation','(*p<0.05, **p<0.01, ***p<0.001)'},'FontSize',16,'FontWeight','bold');
saveas(fig2,fullfile(saveDir,'CorrHeatmap_Discharge.fig'));
close(fig2);
fprintf('Discharge Heatmap saved.\n');

% ---- [3] Combined Heatmap (Square: All Features + All Labels) ----
all_combined = [X_chg, X_dch, X_common, y_targets];
all_combined_names = [X_chg_names, X_dch_names, X_common_names, y_names];
n_feat_all = length(X_chg_names) + length(X_dch_names) + length(X_common_names);

[R_all, P_all] = corr(all_combined, 'Rows', 'pairwise');
n_all = length(all_combined_names);

fig3 = figure('Position', [50, 30, 1800, 1800], 'Color', 'w', 'Visible', 'on');
imagesc(R_all); colormap(cmap); caxis([-1 1]);
cb = colorbar('Location','eastoutside');
cb.Label.String = 'Pearson r'; cb.Label.FontSize = 12; cb.Label.FontWeight = 'bold';
axis square;
set(gca, 'XTick',1:n_all, 'XTickLabel',all_combined_names, 'FontSize',10, 'FontWeight','bold');
set(gca, 'YTick',1:n_all, 'YTickLabel',all_combined_names, 'FontSize',10, 'FontWeight','bold');
set(gca, 'XAxisLocation','top'); xtickangle(45);

hold on;
line([0.5,n_all+0.5],[0.5,n_all+0.5],'Color',[0.3 0.3 0.3],'LineWidth',1.5);
% Feature / Label boundary
line([n_feat_all+0.5,n_feat_all+0.5],[0.5,n_all+0.5],'Color','k','LineWidth',2);
line([0.5,n_all+0.5],[n_feat_all+0.5,n_feat_all+0.5],'Color','k','LineWidth',2);
% Charge / Discharge Feature boundary
n_chg_feat = length(X_chg_names);
line([n_chg_feat+0.5,n_chg_feat+0.5],[0.5,n_all+0.5],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','--');
line([0.5,n_all+0.5],[n_chg_feat+0.5,n_chg_feat+0.5],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','--');

for i=1:n_all
    for j=1:n_all
        if i==j, continue; end
        val=R_all(i,j); pval=P_all(i,j);
        if pval<0.001, star='***'; elseif pval<0.01, star='**'; elseif pval<0.05, star='*'; else, star=''; end
        if abs(val)>0.6, clr='w'; else, clr='k'; end
        text(j,i,sprintf('%.2f%s',val,star),'HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Color',clr);
    end
end
hold off;
title({'All Features & Labels Correlation Matrix','(*p<0.05, **p<0.01, ***p<0.001)'},'FontSize',16,'FontWeight','bold');
saveas(fig3,fullfile(saveDir,'CorrHeatmap_Combined.fig'));
close(fig3);
fprintf('Combined Heatmap saved.\n');

fprintf('\nAll Correlation Heatmaps saved to: %s\n', saveDir);
