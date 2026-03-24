% Vis_SVM.m: Support Vector Machine 결과 시각화 스크립트
% 예측 성능(Scatter, Error Dist) 분석 (SVM은 변수 중요도 추출 방식이 달라 생략)
clear; clc; close all;

set(0, 'DefaultFigureWindowStyle', 'docked');

base_dir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
save_dir = fullfile(base_dir, 'ML_SVM');
res_path = fullfile(save_dir, 'Result_SVM.mat');
fprintf('Loading Results from: %s\n', res_path);
if ~exist(res_path, 'file')
    error('Result file not found! Please run Train_SVM.m first.');
end
load(res_path, 'Results_SVM');

label_names = Results_SVM.label_names;
n_labels = length(label_names);
feature_names = Results_SVM.feature_names;

fprintf('Generating figures for SVM Model (Gaussian) ...\n');

Y_true = Results_SVM.Y_True;
Y_pred = Results_SVM.Y_Pred;

%% 1. Parity Plot (Actual vs Predicted Scatter)
fig_scatter = figure('Name', 'SVM Parity Plot', 'Position', [50, 50, 1600, 400]);
sgtitle('Support Vector Machine - 1:1 Parity Plot', 'FontSize', 14, 'FontWeight', 'bold');

for lbl = 1:n_labels
    subplot(1, n_labels, lbl); hold on; box on; grid on;
    y_t = Y_true(:, lbl);
    y_p = Y_pred(:, lbl);
    
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
fig_err = figure('Name', 'SVM Error Dist', 'Position', [100, 100, 1600, 300]);
sgtitle('Support Vector Machine - Error Distribution', 'FontSize', 14, 'FontWeight', 'bold');

for lbl = 1:n_labels
    err = Y_pred(:, lbl) - Y_true(:, lbl);
    
    subplot(1, n_labels, lbl);
    histogram(err, 20, 'FaceColor', [0.85 0.33 0.10], 'EdgeColor', 'w');
    title(label_names{lbl}, 'Interpreter', 'none');
    xlabel('Error (Pred - Actual) [%]'); ylabel('Frequency');
    grid on;
end

%% 3. Residual Plot (Residuals vs Actuals)
fig_res = figure('Name', 'SVM Residual Plot', 'Position', [150, 150, 1600, 300]);
sgtitle('Support Vector Machine - Residuals vs Actual', 'FontSize', 14, 'FontWeight', 'bold');

for lbl = 1:n_labels
    y_t = Y_true(:, lbl);
    y_p = Y_pred(:, lbl);
    err = y_p - y_t;
    
    subplot(1, n_labels, lbl); hold on; box on; grid on;
    scatter(y_t, err, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
    yline(0, 'k--', 'LineWidth', 1.5);
    title(label_names{lbl}, 'Interpreter', 'none');
    xlabel('Actual [%]'); ylabel('Residuals [%]');
end

%% 4. Hyperparameter Convergence Profile
fig_hp = figure('Name', 'SVM Hyperparameter Tuning', 'Position', [200, 200, 1600, 300]);
sgtitle('Support Vector Machine - Bayesian Optimization Trace', 'FontSize', 14, 'FontWeight', 'bold');

for lbl = 1:n_labels
    hp_hist = Results_SVM.(label_names{lbl}).hp_history;
    
    subplot(1, n_labels, lbl); hold on; box on; grid on;
    if ~isempty(hp_hist)
        plot(1:length(hp_hist), hp_hist, '-o', 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880]);
    end
    title(label_names{lbl}, 'Interpreter', 'none');
    xlabel('Iteration'); ylabel('Min Objective');
end

%% Save Figures
saveas(fig_scatter, fullfile(save_dir, 'SVM_ParityPlot.fig'));
saveas(fig_err,     fullfile(save_dir, 'SVM_ErrorDist.fig'));
saveas(fig_res,     fullfile(save_dir, 'SVM_ResidualPlot.fig'));
saveas(fig_hp,      fullfile(save_dir, 'SVM_Hyperparams.fig'));

fprintf('\n=== Visualization Complete ===\n');
fprintf('All figures successfully saved to:\n  %s\n', save_dir);
