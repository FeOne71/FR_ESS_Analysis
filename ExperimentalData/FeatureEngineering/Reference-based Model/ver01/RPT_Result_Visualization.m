
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Result Visualization Script
% Generates Professional-Grade Plots for MLR vs Random Forest
% Loads data from 'Modeling_Results.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% 1. Load Data
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis';
resultMatPath = fullfile(baseDir, 'Modeling_Results', 'Modeling_Results.mat');
saveDir = fullfile(baseDir, 'Modeling_Results', 'Paper_Figures');

if ~exist(resultMatPath, 'file')
    error('Result file not found! Run RPT_Modeling.m first.');
end
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

fprintf('Loading Modeling Results...\n');
load(resultMatPath);
fprintf('Data Loaded.\n');

% Common Settings
ppt_fontsize_title = 16;
ppt_fontsize_axis  = 14;
ppt_fontsize_tick  = 12;
ppt_fontsize_value = 11;
color_mlr = [0.2 0.5 0.85]; % Blue
color_rf  = [0.85 0.33 0.1]; % Orange
label_units = {'(%)', '(%)', '(%)'}; % Since we updated units to %

%% 2. MLR Specific Visualization
fprintf('\nGenerating MLR Figures...\n');

% 2.1 MLR R2 Score Bar Chart
fig_mlr_r2 = figure('Position', [100, 100, 600, 500], 'Name', 'MLR R2 Score');
mlr_r2_means = squeeze(mean(metrics_MLR(:, :, 3), 1));

b = bar(mlr_r2_means, 0.6, 'FaceColor', color_mlr, 'EdgeColor', 'none');
hold on;
for i = 1:length(mlr_r2_means)
    text(i, mlr_r2_means(i)+0.02, sprintf('%.2f', mlr_r2_means(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', ppt_fontsize_value, 'FontWeight', 'bold');
end
yline(0, '-k');
ylim([0 1.05]);
set(gca, 'XTickLabel', label_names, 'FontSize', ppt_fontsize_tick);
ylabel('R² Score', 'FontSize', ppt_fontsize_axis);
title('MLR Model Performance (R²)', 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
grid on; axis square;
saveas(fig_mlr_r2, fullfile(saveDir, 'MLR_Performance_R2.fig'));

% 2.2 MLR Predicted vs Actual Scatter
fig_mlr_scat = figure('Position', [100, 100, 1200, 400], 'Name', 'MLR Pred vs Actual');
for i = 1:n_labels
    subplot(1, 3, i);
    all_true = vertcat(kfold_Y_true{:});
    all_pred = vertcat(kfold_Y_pred_MLR{:});
    
    scatter(all_true(:, i), all_pred(:, i), 30, color_mlr, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    
    % Draw y=x line
    ax_min = min([all_true(:, i); all_pred(:, i)]);
    ax_max = max([all_true(:, i); all_pred(:, i)]);
    margin = (ax_max - ax_min) * 0.1;
    plot([ax_min-margin, ax_max+margin], [ax_min-margin, ax_max+margin], '--k', 'LineWidth', 1.5);
    
    xlabel(['Actual ' label_names{i}], 'FontSize', ppt_fontsize_axis);
    ylabel(['Predicted ' label_names{i}], 'FontSize', ppt_fontsize_axis);
    title(sprintf('%s (R²=%.2f)', label_names{i}, mlr_r2_means(i)), ...
        'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
    
    grid on; axis square;
    xlim([ax_min-margin, ax_max+margin]);
    ylim([ax_min-margin, ax_max+margin]);
end
sgtitle('MLR: Predicted vs Actual', 'FontSize', 18, 'FontWeight', 'bold');
saveas(fig_mlr_scat, fullfile(saveDir, 'MLR_Predicted_vs_Actual.fig'));


%% 3. RF Specific Visualization
fprintf('\nGenerating RF Figures...\n');

% 3.1 RF R2 Score Bar Chart
fig_rf_r2 = figure('Position', [100, 100, 600, 500], 'Name', 'RF R2 Score');
rf_r2_means = squeeze(mean(metrics_RF(:, :, 3), 1));

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
title('Random Forest Performance (R²)', 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
grid on; axis square;
saveas(fig_rf_r2, fullfile(saveDir, 'RF_Performance_R2.fig'));

% 3.2 RF Predicted vs Actual Scatter
fig_rf_scat = figure('Position', [100, 100, 1200, 400], 'Name', 'RF Pred vs Actual');
for i = 1:n_labels
    subplot(1, 3, i);
    all_true = vertcat(kfold_Y_true{:});
    all_pred = vertcat(kfold_Y_pred_RF{:});
    
    scatter(all_true(:, i), all_pred(:, i), 30, color_rf, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    
    ax_min = min([all_true(:, i); all_pred(:, i)]);
    ax_max = max([all_true(:, i); all_pred(:, i)]);
    margin = (ax_max - ax_min) * 0.1;
    plot([ax_min-margin, ax_max+margin], [ax_min-margin, ax_max+margin], '--k', 'LineWidth', 1.5);
    
    xlabel(['Actual ' label_names{i}], 'FontSize', ppt_fontsize_axis);
    ylabel(['Predicted ' label_names{i}], 'FontSize', ppt_fontsize_axis);
    title(sprintf('%s (R²=%.2f)', label_names{i}, rf_r2_means(i)), ...
        'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
    
    grid on; axis square;
    xlim([ax_min-margin, ax_max+margin]);
    ylim([ax_min-margin, ax_max+margin]);
end
sgtitle('Random Forest: Predicted vs Actual', 'FontSize', 18, 'FontWeight', 'bold');
saveas(fig_rf_scat, fullfile(saveDir, 'RF_Predicted_vs_Actual.fig'));

% 3.3 RF Error Distribution Histogram
fig_err = figure('Position', [100, 100, 1200, 400], 'Name', 'Prediction Error Distribution');
for i = 1:n_labels
    subplot(1, 3, i);
    all_true = vertcat(kfold_Y_true{:});
    all_pred = vertcat(kfold_Y_pred_RF{:});
    errors = all_pred(:, i) - all_true(:, i);
    
    histogram(errors, 20, 'FaceColor', color_rf, 'EdgeColor', 'w', 'FaceAlpha', 0.7);
    hold on;
    xline(0, '--k', 'LineWidth', 1.5);
    
    xlabel(['Error ' label_units{i}], 'FontSize', ppt_fontsize_axis);
    ylabel('Count', 'FontSize', ppt_fontsize_axis);
    title(sprintf('%s Error (RMSE=%.2f)', label_names{i}, sqrt(mean(errors.^2))), ...
        'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
    grid on;
end
sgtitle('Random Forest Prediction Error Distribution', 'FontSize', 18, 'FontWeight', 'bold');
saveas(fig_err, fullfile(saveDir, 'RF_Error_Distribution.fig'));


%% 4. Comparison Visualization (MLR vs RF)
fprintf('\nGenerating Comparison Figures...\n');

% 4.1 Comparison Bar Chart (R2 Score)
fig_comp = figure('Position', [100, 100, 800, 500], 'Name', 'MLR vs RF Comparison');
comparison_data = [mlr_r2_means; rf_r2_means]';

b = bar(comparison_data, 'grouped');
b(1).FaceColor = color_mlr; b(1).EdgeColor = 'none';
b(2).FaceColor = color_rf; b(2).EdgeColor = 'none';

hold on;
% Add value labels
for i = 1:n_labels
    % MLR Label
    text(b(1).XEndPoints(i), mlr_r2_means(i)+0.03, sprintf('%.2f', mlr_r2_means(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', ppt_fontsize_value, 'Color', color_mlr, 'FontWeight', 'bold');
    % RF Label
    text(b(2).XEndPoints(i), rf_r2_means(i)+0.03, sprintf('%.2f', rf_r2_means(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', ppt_fontsize_value, 'Color', color_rf, 'FontWeight', 'bold');
    
    % Highlight Improvement for LLI
    if strcmp(label_names{i}, 'LLI')
        x_center = (b(1).XEndPoints(i) + b(2).XEndPoints(i))/2;
        y_top = max(mlr_r2_means(i), rf_r2_means(i)) + 0.15;
        improvement = (rf_r2_means(i) - mlr_r2_means(i)) / mlr_r2_means(i) * 100;
        
        annotation('textarrow', [0.5 0.5], [0.5 0.4], ... % Dummy coordinates, actual arrow logic below
            'String', sprintf('+%.0f%% Improved', improvement), ...
            'Color', [0.2 0.7 0.2], 'FontSize', 12, 'FontWeight', 'bold');
        
        text(x_center, y_top-0.05, '★', 'Color', 'r', 'FontSize', 20, 'HorizontalAlignment', 'center');
    end
end
ylabel('R² Score', 'FontSize', ppt_fontsize_axis);
ylim([0 1.2]); % Give space for labels
set(gca, 'XTickLabel', label_names, 'FontSize', ppt_fontsize_tick);
legend({'MLR (Linear)', 'Random Forest (Ensemble)'}, 'Location', 'northwest', 'FontSize', 12);
title('Model Performance Comparison (MLR vs RF)', 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
grid on;

saveas(fig_comp, fullfile(saveDir, 'Comparison_MLR_vs_RF.fig'));


%% 5. [NEW] Optimization & Convergence Detail
fprintf('\nGenerating Convergence Figures...\n');

% 5.1 Convergence Plot (Trees vs. OOB Error)
if exist('rmse_heatmap', 'var')
    for lbl = 1:n_labels
        fig_conv = figure('Position', [100, 100, 700, 500], 'Name', ['RF Convergence: ' label_names{lbl}]);
        hold on;
        line_styles = {'-o', '-s', '-d', '-^'};
        colors_conv = lines(length(hp_leaf_grid));
        
        for i = 1:length(hp_leaf_grid)
            plot(hp_trees_grid, rmse_heatmap(i, :, lbl), line_styles{mod(i-1,4)+1}, ...
                'Color', colors_conv(i,:), 'LineWidth', 2, 'MarkerSize', 8, ...
                'DisplayName', sprintf('MinLeafSize=%d', hp_leaf_grid(i)));
        end
        
        xlabel('Number of Trees', 'FontSize', ppt_fontsize_axis);
        ylabel(['OOB RMSE ' label_units{lbl}], 'FontSize', ppt_fontsize_axis);
        title(['RF Optimization & Convergence: ' label_names{lbl}], 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 10);
        grid on;
        
        % Highlight selected point (500, 1) - Assuming index 1 for Leaf 1 (ri=1) and index 4 for Tree 500 (ci=4)
        % Note: Index might vary if grid changes, but currently [1,3,5,10] and [50,100,200,500]
        plot(500, rmse_heatmap(1, 4, lbl), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'y', 'HandleVisibility', 'off');
        text(500, rmse_heatmap(1, 4, lbl)*0.95, 'Selected (Optimal)', 'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        
        saveas(fig_conv, fullfile(saveDir, ['RF_Convergence_Plot_' label_names{lbl} '.fig']));
    end
end

%% 6. [NEW] Model Interpretability (Feature Importance)
fprintf('\nGenerating Feature Importance Figures...\n');

fig_fi = figure('Position', [100, 100, 1500, 450], 'Name', 'Feature Importance');
for lbl = 1:n_labels
    % Modeling.m에서 predictorImportance를 계산해서 저장하지 않았으므로, 
    % 여기서 가볍게 전체 데이터로 다시 계산 (또는 저장된 값이 있다면 사용)
    subplot(1, 3, lbl);
    
    % Simple Fit to get Importance
    t = templateTree('MinLeafSize', 1);
    rf_temp = fitrensemble(X, Y(:, lbl), 'Method', 'Bag', 'NumLearningCycles', 100, 'Learners', t);
    imp = predictorImportance(rf_temp);
    
    [imp_sorted, sort_idx] = sort(imp, 'descend');
    barh(imp_sorted, 'FaceColor', [0.3 0.6 0.9]);
    set(gca, 'YTick', 1:length(feature_names), 'YTickLabel', feature_names(sort_idx), ...
        'YDir', 'reverse', 'FontSize', 9);
    xlabel('Importance', 'FontSize', 11);
    title(sprintf('%s Feature Importance', label_names{lbl}), 'FontSize', 13, 'FontWeight', 'bold');
    grid on;
end
saveas(fig_fi, fullfile(saveDir, 'RF_Feature_Importance_Detail.fig'));

%% 7. [NEW] Model Diagnostics (Residual analysis)
fprintf('\nGenerating Residual Figures...\n');

fig_resid = figure('Position', [100, 100, 1200, 400], 'Name', 'Residual Analysis');
for i = 1:n_labels
    subplot(1, 3, i);
    all_true = vertcat(kfold_Y_true{:});
    all_pred = vertcat(kfold_Y_pred_RF{:});
    residuals = all_pred(:, i) - all_true(:, i);
    
    scatter(all_true(:, i), residuals, 30, [0.2 0.6 0.4], 'filled', 'MarkerFaceAlpha', 0.5);
    hold on;
    yline(0, '--r', 'LineWidth', 2);
    
    xlabel(['Actual ' label_names{i}], 'FontSize', ppt_fontsize_axis);
    ylabel('Residual (Pred - Act)', 'FontSize', ppt_fontsize_axis);
    title(sprintf('%s Residuals (Std=%.2f)', label_names{i}, std(residuals)), ...
        'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
    grid on;
end
sgtitle('Random Forest: Residual Analysis', 'FontSize', 18, 'FontWeight', 'bold');
saveas(fig_resid, fullfile(saveDir, 'RF_Residual_Analysis.fig'));


fprintf('\nAll Visualization Completed.\n');
fprintf('Figures saved in: %s\n', saveDir);
