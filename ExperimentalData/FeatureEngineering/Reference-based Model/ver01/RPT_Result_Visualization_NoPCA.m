
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Result Visualization (No PCA) - ENHANCED
% Follows the professional style of RPT_Result_Visualization.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% 1. Load Data
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis';
resultMatPath = fullfile(baseDir, 'Modeling_Results', 'Modeling_Results_NoPCA.mat');
saveDir = fullfile(baseDir, 'Modeling_Results', 'Paper_Figures_NoPCA');

if ~exist(resultMatPath, 'file')
    error('Result file not found! Run RPT_Modeling_NoPCA.m first.');
end
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

fprintf('Loading No-PCA Modeling Results...\n');
load(resultMatPath);
fprintf('Data Loaded.\n');

% Common Settings (Standardized across scripts)
ppt_fontsize_title = 16;
ppt_fontsize_axis  = 14;
ppt_fontsize_tick  = 12;
ppt_fontsize_value = 11;
color_rf  = [0.4 0.2 0.7]; % Dark Purple for NoPCA
label_units = {'(%)', '(%)', '(%)'};
n_labels = length(label_names);

%% 2. RF (No PCA) R2 Score Bar Chart
fprintf('\nGenerating R2 Performance Bar Chart...\n');
fig_rf_r2 = figure('Position', [100, 100, 600, 500], 'Name', 'RF NoPCA R2 Score');
rf_r2_means = squeeze(mean(metrics_RF_Raw(:, :, 3), 1));

b = bar(rf_r2_means, 0.6, 'FaceColor', color_rf, 'EdgeColor', 'none');
hold on;
for i = 1:length(rf_r2_means)
    text(i, rf_r2_means(i)+0.02, sprintf('%.2f', rf_r2_means(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', ppt_fontsize_value, 'FontWeight', 'bold');
end
yline(0, '-k');
ylim([0 1.05]);
set(gca, 'XTickLabel', label_names, 'FontSize', ppt_fontsize_tick);
ylabel('R² Score', 'FontSize', ppt_fontsize_axis);
title('RF Model Performance: No PCA (R²)', 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
grid on; axis square;
saveas(fig_rf_r2, fullfile(saveDir, 'RF_NoPCA_Performance_R2.fig'));

%% 3. Predicted vs Actual Scatter Plots
fprintf('\nGenerating Predicted vs Actual Scatter Plots...\n');
fig_scat = figure('Position', [100, 100, 1500, 450], 'Name', 'RF NoPCA Pred vs Actual');

all_true_mat = vertcat(kfold_Y_true{:});
all_pred_mat = vertcat(kfold_Y_pred_RF{:});

for i = 1:n_labels
    subplot(1, 3, i);
    scatter(all_true_mat(:, i), all_pred_mat(:, i), 35, color_rf, 'filled', 'MarkerFaceAlpha', 0.5);
    hold on;
    
    % Draw 1:1 line
    ax_min = min([all_true_mat(:, i); all_pred_mat(:, i)]);
    ax_max = max([all_true_mat(:, i); all_pred_mat(:, i)]);
    margin = (ax_max - ax_min) * 0.1;
    plot([ax_min-margin, ax_max+margin], [ax_min-margin, ax_max+margin], '--k', 'LineWidth', 1.5);
    
    xlabel(['Actual ' label_names{i} ' ' label_units{i}], 'FontSize', ppt_fontsize_axis);
    ylabel(['Predicted ' label_names{i} ' ' label_units{i}], 'FontSize', ppt_fontsize_axis);
    title(sprintf('%s (R²=%.2f)', label_names{i}, rf_r2_means(i)), 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
    
    grid on; axis square;
    xlim([ax_min-margin, ax_max+margin]);
    ylim([ax_min-margin, ax_max+margin]);
end
sgtitle('Random Forest (No PCA): Predicted vs Actual', 'FontSize', 18, 'FontWeight', 'bold');
saveas(fig_scat, fullfile(saveDir, 'RF_NoPCA_Predicted_vs_Actual.fig'));

%% 4. Error Distribution Histograms
fprintf('\nGenerating Error Distribution Histograms...\n');
fig_err = figure('Position', [100, 100, 1200, 400], 'Name', 'RF NoPCA Error Distribution');

for i = 1:n_labels
    subplot(1, 3, i);
    errors = all_pred_mat(:, i) - all_true_mat(:, i);
    
    histogram(errors, 20, 'FaceColor', color_rf, 'EdgeColor', 'w', 'FaceAlpha', 0.7);
    hold on;
    xline(0, '--k', 'LineWidth', 1.5);
    
    rmse_fold = sqrt(mean(errors.^2));
    xlabel(['Error ' label_units{i}], 'FontSize', ppt_fontsize_axis);
    ylabel('Count', 'FontSize', ppt_fontsize_axis);
    title(sprintf('%s Error (RMSE=%.2f)', label_names{i}, rmse_fold), ...
        'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
    grid on;
end
sgtitle('Random Forest (No PCA): Error Distribution', 'FontSize', 18, 'FontWeight', 'bold');
saveas(fig_err, fullfile(saveDir, 'RF_NoPCA_Error_Distribution.fig'));

%% 5. Feature Importance
fprintf('\nGenerating Feature Importance Figures...\n');
fig_fi = figure('Position', [100, 100, 1500, 500], 'Name', 'RF NoPCA Feature Importance');

for lbl = 1:n_labels
    subplot(1, 3, lbl);
    
    % Use saved X, Y for consistency
    t = templateTree('MinLeafSize', 1);
    rf_temp = fitrensemble(X, Y(:, lbl), 'Method', 'Bag', 'NumLearningCycles', 500, 'Learners', t);
    imp = predictorImportance(rf_temp);
    
    [imp_sorted, sort_idx] = sort(imp, 'descend');
    barh(imp_sorted, 'FaceColor', color_rf, 'EdgeColor', 'none');
    set(gca, 'YTick', 1:length(feature_names), 'YTickLabel', feature_names(sort_idx), ...
        'YDir', 'reverse', 'FontSize', 9);
    xlabel('Importance Score', 'FontSize', 11);
    title(sprintf('%s Critical Features', label_names{lbl}), 'FontSize', 13, 'FontWeight', 'bold');
    grid on;
end
sgtitle('Random Forest (No PCA): Feature Importance Ranking', 'FontSize', 18, 'FontWeight', 'bold');
saveas(fig_fi, fullfile(saveDir, 'RF_NoPCA_Feature_Importance.fig'));

%% 5. [NEW] Optimization & Convergence Detail (No PCA)
fprintf('\nGenerating Convergence Figures...\n');

if exist('rmse_heatmap', 'var')
    for lbl = 1:n_labels
        fig_conv = figure('Position', [100, 100, 700, 500], 'Name', ['RF Convergence (NoPCA): ' label_names{lbl}]);
        hold on;
        line_styles = {'-o', '-s', '-d', '-^'};
        colors_conv = lines(length(hp_leaf_grid));
        
        for i = 1:length(hp_leaf_grid)
            plot(hp_trees_grid, rmse_heatmap(i, :, lbl), line_styles{mod(i-1,4)+1}, ...
                'Color', colors_conv(i,:), 'LineWidth', 2, 'MarkerSize', 8, ...
                'DisplayName', sprintf('MinLeafSize=%d', hp_leaf_grid(i)));
        end
        
        xlabel('Number of Trees', 'FontSize', ppt_fontsize_axis);
        ylabel(['OOB RMSE (%)'], 'FontSize', ppt_fontsize_axis);
        title(['RF Optimization & Convergence (No PCA): ' label_names{lbl}], 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 10);
        grid on;
        
        % Highlight selected point (500, 1)
        plot(500, rmse_heatmap(1, 4, lbl), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'y', 'HandleVisibility', 'off');
        text(500, rmse_heatmap(1, 4, lbl)*0.95, 'Selected (Optimal)', 'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        
        saveas(fig_conv, fullfile(saveDir, ['RF_NoPCA_Convergence_Plot_' label_names{lbl} '.fig']));
    end
end

%% 6. [NEW] Model Diagnostics (Residual Analysis)
fprintf('\nGenerating Residual Figures...\n');

fig_resid = figure('Position', [100, 100, 1200, 400], 'Name', 'RF NoPCA Residual Analysis');
all_true_mat = vertcat(kfold_Y_true{:});
all_pred_mat = vertcat(kfold_Y_pred_RF{:});

for i = 1:n_labels
    subplot(1, 3, i);
    residuals = all_pred_mat(:, i) - all_true_mat(:, i);
    
    scatter(all_true_mat(:, i), residuals, 30, color_rf, 'filled', 'MarkerFaceAlpha', 0.5);
    hold on;
    yline(0, '--r', 'LineWidth', 2);
    
    xlabel(['Actual ' label_names{i} ' (%)'], 'FontSize', ppt_fontsize_axis);
    ylabel('Residual (Pred - Act)', 'FontSize', ppt_fontsize_axis);
    title(sprintf('%s Residuals (Std=%.2f)', label_names{i}, std(residuals)), ...
        'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
    grid on;
end
sgtitle('Random Forest (No PCA): Residual Analysis', 'FontSize', 18, 'FontWeight', 'bold');
saveas(fig_resid, fullfile(saveDir, 'RF_NoPCA_Residual_Analysis.fig'));


fprintf('\nAll Professional Visualizations (No PCA) Completed.\n');
fprintf('Figures saved in: %s\n', saveDir);
