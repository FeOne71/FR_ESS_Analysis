
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Comparison Visualization: MLR vs. Random Forest (No PCA)
% Purpose: Compare baseline MLR with optimized No-PCA RF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% 1. Load Data
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis';
path_main  = fullfile(baseDir, 'Modeling_Results', 'Modeling_Results.mat');
path_nopca = fullfile(baseDir, 'Modeling_Results', 'Modeling_Results_NoPCA.mat');
saveDir    = fullfile(baseDir, 'Modeling_Results', 'Paper_Figures_Comparison');

if ~exist(path_main, 'file') || ~exist(path_nopca, 'file')
    error('Modeling result files missing. Run both RPT_Modeling.m and RPT_Modeling_NoPCA.m first.');
end
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

% Load MLR results from main file
main_data = load(path_main, 'metrics_MLR', 'label_names');
% Load No-PCA RF results from NoPCA file
nopca_data = load(path_nopca, 'metrics_RF_Raw');

label_names = main_data.label_names;
metrics_MLR = main_data.metrics_MLR;
metrics_RF  = nopca_data.metrics_RF_Raw;

% Common Plot Settings
ppt_fontsize_title = 16;
ppt_fontsize_axis  = 14;
ppt_fontsize_tick  = 12;
ppt_fontsize_value = 11;
color_mlr = [0.6 0.6 0.6];    % Neutral gray for baseline
color_rf  = [0.2 0.6 0.2];    % Strong green for the "Winner" (No-PCA RF)
n_labels = length(label_names);

%% 2. Calculate Statistics (Means and SD)
% Index 3 is R2
mlr_r2_means = squeeze(mean(metrics_MLR(:, :, 3), 1));
mlr_r2_stds  = squeeze(std(metrics_MLR(:, :, 3), 0, 1));
rf_r2_means  = squeeze(mean(metrics_RF(:, :, 3), 1));
rf_r2_stds   = squeeze(std(metrics_RF(:, :, 3), 0, 1));

%% 3. Generate Comparison Bar Chart (R2)
fprintf('\nGenerating MLR vs No-PCA RF Comparison Figure...\n');
fig_comp = figure('Position', [100, 100, 900, 600], 'Name', 'Final Comparison: MLR vs No-PCA RF');

data_to_plot = [mlr_r2_means; rf_r2_means]';
b = bar(data_to_plot, 'grouped', 'EdgeColor', 'none');
b(1).FaceColor = color_mlr;
b(2).FaceColor = color_rf;

hold on;

% Add Error Bars (Horizontal placement adjustment for grouped bars)
group_width = min(0.8, n_labels/(n_labels+1.5));
for i = 1:n_labels
    % MLR Error Bar
    x_mlr = i - group_width/4;
    errorbar(x_mlr, mlr_r2_means(i), mlr_r2_stds(i), 'k.', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % RF Error Bar
    x_rf = i + group_width/4;
    errorbar(x_rf, rf_r2_means(i), rf_r2_stds(i), 'k.', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % Value Labels
    text(x_mlr, mlr_r2_means(i) + mlr_r2_stds(i) + 0.02, sprintf('%.2f', mlr_r2_means(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', ppt_fontsize_value, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);
    text(x_rf, rf_r2_means(i) + rf_r2_stds(i) + 0.02, sprintf('%.2f', rf_r2_means(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', ppt_fontsize_value, 'FontWeight', 'bold', 'Color', color_rf);
    
    % Improvement Text (%)
    improvement = (rf_r2_means(i) - mlr_r2_means(i)) / mlr_r2_means(i) * 100;
    if improvement > 0
        text(i, max(rf_r2_means(i), mlr_r2_means(i)) + 0.12, sprintf('+%.1f%%', improvement), ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
        text(i, max(rf_r2_means(i), mlr_r2_means(i)) + 0.16, '▲', 'HorizontalAlignment', 'center', 'Color', 'r');
    end
end

% Styling
ylabel('R² Score (Cross-Validation Avg)', 'FontSize', ppt_fontsize_axis);
set(gca, 'XTickLabel', label_names, 'FontSize', ppt_fontsize_tick);
ylim([0 1.2]);
grid on;
legend({'MLR (Base)', 'Random Forest (No PCA)'}, 'Location', 'northwest', 'FontSize', 12);
title('Final Model Performance Comparison: MLR vs. No-PCA RF', 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');

% Add Annotation Box
annotation('textbox', [0.65, 0.15, 0.25, 0.1], 'String', ...
    {['Avg SOH Improvement: ' sprintf('%.1f%%', (rf_r2_means(1)/mlr_r2_means(1)-1)*100)], ...
     ['Avg LLI Improvement: ' sprintf('%.1f%%', (rf_r2_means(2)/mlr_r2_means(2)-1)*100)], ...
     ['Avg LAM Improvement: ' sprintf('%.1f%%', (rf_r2_means(3)/mlr_r2_means(3)-1)*100)]}, ...
    'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', color_rf, 'LineWidth', 1.5);

saveas(fig_comp, fullfile(saveDir, 'Comparison_MLR_vs_RFNoPCA_R2.fig'));

fprintf('\nFinal Comparison Visualization Completed.\n');
fprintf('Figure saved in: %s\n', saveDir);
