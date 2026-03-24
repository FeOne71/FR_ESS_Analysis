% Vis_RF.m: Random Forest 결과 시각화 스크립트
% 예측 성능(Scatter, Error Dist) 및 최종 투입 변수 중요도(Feature Importance) 분석
clear; clc; close all;

set(0, 'DefaultFigureWindowStyle', 'docked');

base_dir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
save_dir = fullfile(base_dir, 'ML_RandomForest');
res_path = fullfile(save_dir, 'Result_RandomForest.mat');
fprintf('Loading Results from: %s\n', res_path);
if ~exist(res_path, 'file')
    error('Result file not found! Please run Train_RF.m first.');
end
load(res_path, 'Results_RF');

label_names = Results_RF.label_names;
n_labels = length(label_names);
feature_names = Results_RF.Merged.feature_names;

% We now only have the 'Merged' feature set and Group K-Fold CV.
mod_name = 'Merged';
fprintf('Generating figures for Integrated Model (%s) ...\n', mod_name);

%% 1. Parity Plot (Actual vs Predicted Scatter)
fig_scatter = figure('Name', 'RF Parity Plot', 'Position', [50, 50, 1600, 400]);
sgtitle('Random Forest - 1:1 Parity Plot', 'FontSize', 14, 'FontWeight', 'bold');

for lbl = 1:n_labels
    y_t = Results_RF.(mod_name).(label_names{lbl}).Y_true;
    y_p = Results_RF.(mod_name).(label_names{lbl}).Y_pred;
    
    subplot(1, n_labels, lbl); hold on; box on; grid on;
    
    unit_str = '%';
    rmse = sqrt(mean((y_t - y_p).^2));
    r2 = 1 - sum((y_t - y_p).^2) / sum((y_t - mean(y_t)).^2);
    
    scatter(y_t, y_p, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
    
    min_v = min([y_t; y_p]); max_v = max([y_t; y_p]);
    range_margin = max(abs((max_v - min_v) * 0.1), 0.1);
    plot([min_v-range_margin, max_v+range_margin], [min_v-range_margin, max_v+range_margin], 'k--', 'LineWidth', 1.5);
    
    axis([min_v-range_margin max_v+range_margin min_v-range_margin max_v+range_margin]);
    title(sprintf('%s (R^2=%.3f, RMSE=%.3f %s)', label_names{lbl}, r2, rmse, unit_str), 'Interpreter', 'none');
    xlabel(sprintf('Actual (%s)', unit_str)); ylabel(sprintf('Predicted (%s)', unit_str));
    axis square;
end

%% 2. Error Distribution (Histogram)
fig_err = figure('Name', 'RF Error Dist', 'Position', [100, 100, 1600, 300]);
sgtitle('Random Forest - Error Distribution', 'FontSize', 14, 'FontWeight', 'bold');

for lbl = 1:n_labels
    y_t = Results_RF.(mod_name).(label_names{lbl}).Y_true;
    y_p = Results_RF.(mod_name).(label_names{lbl}).Y_pred;
    err = y_p - y_t;
    
    subplot(1, n_labels, lbl);
    histogram(err, 20, 'FaceColor', [0.85 0.33 0.10], 'EdgeColor', 'w');
    title(label_names{lbl}, 'Interpreter', 'none');
    xlabel('Error (Pred - Actual) [%]'); ylabel('Frequency');
    grid on;
end

%% 3. Residual Plot (Residuals vs Actuals)
fig_res = figure('Name', 'RF Residual Plot', 'Position', [150, 150, 1600, 300]);
sgtitle('Random Forest - Residuals vs Actual', 'FontSize', 14, 'FontWeight', 'bold');

for lbl = 1:n_labels
    y_t = Results_RF.(mod_name).(label_names{lbl}).Y_true;
    y_p = Results_RF.(mod_name).(label_names{lbl}).Y_pred;
    err = y_p - y_t;
    
    subplot(1, n_labels, lbl); hold on; box on; grid on;
    scatter(y_t, err, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
    yline(0, 'k--', 'LineWidth', 1.5);
    title(label_names{lbl}, 'Interpreter', 'none');
    xlabel('Actual [%]'); ylabel('Residuals [%]');
end

%% 4. Hyperparameter Convergence Profile
fig_hp = figure('Name', 'RF Hyperparameter Tuning', 'Position', [200, 200, 1600, 300]);
sgtitle('Random Forest - OOB Error vs. Num Trees', 'FontSize', 14, 'FontWeight', 'bold');

for lbl = 1:n_labels
    hp_hist = Results_RF.(mod_name).(label_names{lbl}).hp_history;
    
    subplot(1, n_labels, lbl); hold on; box on; grid on;
    if ~isempty(hp_hist)
        leaf_sizes = unique(hp_hist(:, 2));
        colors = lines(length(leaf_sizes));
        for i = 1:length(leaf_sizes)
            idx = hp_hist(:, 2) == leaf_sizes(i);
            plot(hp_hist(idx, 1), hp_hist(idx, 3), '-o', 'LineWidth', 1.5, 'Color', colors(i,:), ...
                'DisplayName', sprintf('MinLeafSize = %d', leaf_sizes(i)));
        end
        legend('Location', 'best');
    end
    title(label_names{lbl}, 'Interpreter', 'none');
    xlabel('Number of Trees'); ylabel('OOB RMSE');
end

%% 5. Feature Importance
fig_imp = figure('Name', 'RF Feature Importance', 'Position', [250, 250, 1800, 450]);
sgtitle('Random Forest - Feature Importance (Final Model)', 'FontSize', 14, 'FontWeight', 'bold');

for lbl = 1:n_labels
    rf_final = Results_RF.(mod_name).(label_names{lbl}).Model;
    imp = oobPermutedPredictorImportance(rf_final);
    
    [sorted_imp, sort_idx] = sort(imp, 'ascend');
    sorted_names = feature_names(sort_idx);
    
    subplot(1, n_labels, lbl);
    barh(sorted_imp, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'none');
    set(gca, 'YTick', 1:length(feature_names), 'YTickLabel', sorted_names, 'TickLabelInterpreter', 'none', 'FontSize', 8);
    title(label_names{lbl}, 'Interpreter', 'none');
    xlabel('Importance Score');
    grid on;
end

%% Save Figures
saveas(fig_scatter, fullfile(save_dir, 'RF_ParityPlot.fig'));
saveas(fig_err,     fullfile(save_dir, 'RF_ErrorDist.fig'));
saveas(fig_res,     fullfile(save_dir, 'RF_ResidualPlot.fig'));
saveas(fig_hp,      fullfile(save_dir, 'RF_Hyperparams.fig'));
saveas(fig_imp,     fullfile(save_dir, 'RF_Importance.fig'));

fprintf('\n=== Visualization Complete ===\n');
fprintf('All figures successfully saved to:\n  %s\n', save_dir);
