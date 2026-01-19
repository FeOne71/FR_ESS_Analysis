%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 06_DriveCycle_XGBoost_Evaluation.m
% XGBoost 모델 평가 및 시각화
% 
% 목적: 
% - 4-Fold CV 결과 종합 및 분석
% - 예측 vs 실제 산점도
% - 잔차 플롯
% - Feature Importance 시각화
% - 성능 리포트 생성
%
% 입력:
% - XGBoost_Training_Results.mat
% - XGBoost_Fold*_Model.mat
%
% 출력:
% - 시각화 그래프들 (figures 폴더)
% - 평가 리포트
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('=== Drive Cycle XGBoost Evaluation ===\n');

%% Configuration - User Settings
% =========================================================================
inputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
outputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
figuresDir = fullfile(outputDir, 'figures', 'XGBoost_Evaluation');
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

n_folds = 4;
% =========================================================================

%% Load training results (Multi-Target Support)
fprintf('\n=== Loading Training Results ===\n');

% Try to load overall summary first (for multi-target evaluation)
overall_summary_file = fullfile(inputDir, 'XGBoost_AllTargets_Summary.mat');
if exist(overall_summary_file, 'file')
    load(overall_summary_file);
    fprintf('Loaded overall summary (Multi-Target): %s\n', overall_summary_file);
    fprintf('Found %d target variables: %s\n', length(targetVars), strjoin(targetVars, ', '));
    
    % Check if targetEventType exists (for separated models)
    if exist('targetEventType', 'var')
        fprintf('Model type: %s\n', targetEventType);
    end
    
    % Load individual target results for detailed evaluation
    use_multi_target = true;
else
    % Fallback: Try to load single target result (backward compatibility)
    results_file = fullfile(inputDir, 'XGBoost_Training_Results.mat');
    if ~exist(results_file, 'file')
        error('Training results not found. Please run DriveCycle_XGBoost_Training_05.m first');
    end
    
    load(results_file);
    fprintf('Loaded single target training results: %s\n', results_file);
    
    % Check if targetEventType exists (for separated models)
    if exist('targetEventType', 'var')
        fprintf('Model type: %s\n', targetEventType);
    end
    
    % Convert to multi-target format for consistency
    targetVars = {targetVar};
    all_target_results = struct();
    all_target_results.(targetVar) = struct();
    all_target_results.(targetVar).fold_results = fold_results;
    all_target_results.(targetVar).overall_mae = overall_mae;
    all_target_results.(targetVar).overall_rmse = overall_rmse;
    all_target_results.(targetVar).overall_mape = overall_mape;
    all_target_results.(targetVar).overall_r2 = overall_r2;
    all_target_results.(targetVar).overall_bias = overall_bias;
    all_target_results.(targetVar).overall_relative_bias = overall_relative_bias;
    all_target_results.(targetVar).fold_mae_mean = mean(fold_mae);
    all_target_results.(targetVar).fold_mae_std = std(fold_mae);
    all_target_results.(targetVar).fold_rmse_mean = mean(fold_rmse);
    all_target_results.(targetVar).fold_rmse_std = std(fold_rmse);
    all_target_results.(targetVar).fold_r2_mean = mean(fold_r2);
    all_target_results.(targetVar).fold_r2_std = std(fold_r2);
    
    use_multi_target = false;
end

% Note: SOC data loading will be done after loading individual target results
% (moved to after target results are loaded to access all_fold_ids and all_actuals)

%% Process each target variable (Multi-Target Support)
if use_multi_target && exist('targetVars', 'var') && length(targetVars) > 1
    % Multi-target: Process each target separately
    fprintf('\n=== Multi-Target Evaluation ===\n');
    fprintf('Found %d target variables. Processing each separately...\n', length(targetVars));
    
    for targetIdx = 1:length(targetVars)
        targetVar = targetVars{targetIdx};
        fprintf('\n');
        fprintf('########################################################################\n');
        fprintf('### Evaluating Model for Target Variable: %s (%d/%d) ###\n', targetVar, targetIdx, length(targetVars));
        fprintf('########################################################################\n');
        
        % Load individual target results
        target_results_file = fullfile(inputDir, sprintf('XGBoost_%s_Training_Results.mat', targetVar));
        if ~exist(target_results_file, 'file')
            fprintf('WARNING: Results file not found for %s: %s\n', targetVar, target_results_file);
            fprintf('Skipping evaluation for this target...\n');
            continue;
        end
        
        % Load target-specific results
        load(target_results_file);
        fprintf('Loaded results for target: %s\n', targetVar);
        
        % Load SOC data from test sets (for SOC-wise analysis)
        % This must be done after loading target results to access all_fold_ids and all_actuals
        all_soc_data = [];
        if exist('all_fold_ids', 'var') && exist('all_actuals', 'var')
            for fold_num = 1:n_folds
                test_file = fullfile(inputDir, sprintf('DriveCycle_Fold%d_Test.mat', fold_num));
                if exist(test_file, 'file')
                    load(test_file, 'testData');
                    if istable(testData) && ismember('SOC', testData.Properties.VariableNames)
                        fold_soc = testData.SOC(~isnan(all_actuals(all_fold_ids == fold_num)));
                        all_soc_data = [all_soc_data; fold_soc(:)];
                    end
                end
            end
        end
        
        % If SOC data is not available, try to extract from fold results
        if isempty(all_soc_data) && exist('fold_results', 'var') && isfield(fold_results, 'Fold1')
            % Extract SOC from test data stored in fold results
            all_soc_data = [];
            for fold_num = 1:n_folds
                fold_name = sprintf('Fold%d', fold_num);
                if isfield(fold_results, fold_name) && isfield(fold_results.(fold_name), 'testData')
                    testData = fold_results.(fold_name).testData;
                    if istable(testData) && ismember('SOC', testData.Properties.VariableNames)
                        valid_idx = ~isnan(fold_results.(fold_name).y_test);
                        all_soc_data = [all_soc_data; testData.SOC(valid_idx)];
                    end
                end
            end
        end
        
        % Create target-specific figures directory
        target_figuresDir = fullfile(figuresDir, targetVar);
        if ~exist(target_figuresDir, 'dir')
            mkdir(target_figuresDir);
        end
        original_figuresDir = figuresDir;
        figuresDir = target_figuresDir;
        
        % Continue with evaluation code below (rest of script)
        % [The rest of the script will process this target]
    end
else
    % Single target: Use existing code
    if ~use_multi_target
        % Already loaded single target results
        targetVar = targetVars{1};  % Use first (and only) target
    end
    
    % Load SOC data from test sets (for SOC-wise analysis)
    % This must be done after loading target results to access all_fold_ids and all_actuals
    all_soc_data = [];
    if exist('all_fold_ids', 'var') && exist('all_actuals', 'var')
        for fold_num = 1:n_folds
            test_file = fullfile(inputDir, sprintf('DriveCycle_Fold%d_Test.mat', fold_num));
            if exist(test_file, 'file')
                load(test_file, 'testData');
                if istable(testData) && ismember('SOC', testData.Properties.VariableNames)
                    fold_soc = testData.SOC(~isnan(all_actuals(all_fold_ids == fold_num)));
                    all_soc_data = [all_soc_data; fold_soc(:)];
                end
            end
        end
    end
    
    % If SOC data is not available, try to extract from fold results
    if isempty(all_soc_data) && exist('fold_results', 'var') && isfield(fold_results, 'Fold1')
        % Extract SOC from test data stored in fold results
        all_soc_data = [];
        for fold_num = 1:n_folds
            fold_name = sprintf('Fold%d', fold_num);
            if isfield(fold_results, fold_name) && isfield(fold_results.(fold_name), 'testData')
                testData = fold_results.(fold_name).testData;
                if istable(testData) && ismember('SOC', testData.Properties.VariableNames)
                    valid_idx = ~isnan(fold_results.(fold_name).y_test);
                    all_soc_data = [all_soc_data; testData.SOC(valid_idx)];
                end
            end
        end
    end
    
    % Create target-specific figures directory
    target_figuresDir = fullfile(figuresDir, targetVar);
    if ~exist(target_figuresDir, 'dir')
        mkdir(target_figuresDir);
    end
    original_figuresDir = figuresDir;
    figuresDir = target_figuresDir;
end

%% Summary Statistics
fprintf('\n=== Cross-Validation Summary ===\n');
fprintf('Number of folds: %d\n', n_folds);
fprintf('Target variable: %s\n', targetVar);
if exist('targetEventType', 'var')
    fprintf('Event type: %s\n', targetEventType);
    if strcmp(targetEventType, 'Charge')
        fprintf('  -> CHARGE-only model evaluation\n');
    elseif strcmp(targetEventType, 'Discharge')
        fprintf('  -> DISCHARGE-only model evaluation\n');
    else
        fprintf('  -> INTEGRATED model evaluation\n');
    end
end

fprintf('\nPerformance Metrics (Mean ± Std):\n');
fprintf('  MAE:  %.4f ± %.4f\n', mean(fold_mae), std(fold_mae));
fprintf('  RMSE: %.4f ± %.4f\n', mean(fold_rmse), std(fold_rmse));
fprintf('  MAPE: %.2f%% ± %.2f%%\n', mean(fold_mape), std(fold_mape));
fprintf('  R²:   %.4f ± %.4f\n', mean(fold_r2), std(fold_r2));

fprintf('\nOverall Performance (All Test Data):\n');
fprintf('  MAE:  %.4f\n', overall_mae);
fprintf('  RMSE: %.4f\n', overall_rmse);
fprintf('  MAPE: %.2f%%\n', overall_mape);
fprintf('  R²:   %.4f\n', overall_r2);

%% SOC-wise Performance Analysis
fprintf('\n=== Performance by SOC Level ===\n');

% Extract SOC data from fold results
all_soc_data = [];
for fold_num = 1:n_folds
    fold_name = sprintf('Fold%d', fold_num);
    if isfield(fold_results, fold_name)
        if isfield(fold_results.(fold_name), 'soc_test')
            fold_soc = fold_results.(fold_name).soc_test;
            fold_valid = ~isnan(fold_results.(fold_name).y_test);
            all_soc_data = [all_soc_data; fold_soc(fold_valid)];
        elseif isfield(fold_results.(fold_name), 'testData')
            testData = fold_results.(fold_name).testData;
            if istable(testData) && ismember('SOC', testData.Properties.VariableNames)
                fold_valid = ~isnan(fold_results.(fold_name).y_test);
                all_soc_data = [all_soc_data; testData.SOC(fold_valid)];
            end
        end
    end
end

% If still empty, try loading from test files
if isempty(all_soc_data)
    for fold_num = 1:n_folds
        test_file = fullfile(inputDir, sprintf('DriveCycle_Fold%d_Test.mat', fold_num));
        if exist(test_file, 'file')
            load(test_file, 'testData');
            if ismember('SOC', testData.Properties.VariableNames)
                % Match indices with predictions
                fold_start = sum(all_fold_ids < fold_num) + 1;
                fold_end = sum(all_fold_ids <= fold_num);
                if fold_end >= fold_start
                    all_soc_data = [all_soc_data; testData.SOC(fold_start:fold_end)];
                end
            end
        end
    end
end

if ~isempty(all_soc_data) && length(all_soc_data) == length(all_actuals)
    uniqueSOCs = unique(all_soc_data(~isnan(all_soc_data)));
    uniqueSOCs = sort(uniqueSOCs);
    
    soc_performance = struct();
    
    for s = 1:length(uniqueSOCs)
        soc = uniqueSOCs(s);
        idx = (all_soc_data == soc) & ~isnan(all_actuals) & ~isnan(all_predictions);
        
        if sum(idx) > 0
            soc_actuals = all_actuals(idx);
            soc_predictions = all_predictions(idx);
            
            soc_mae = mean(abs(soc_actuals - soc_predictions));
            soc_rmse = sqrt(mean((soc_actuals - soc_predictions).^2));
            soc_mape = mean(abs((soc_actuals - soc_predictions) ./ soc_actuals)) * 100;
            soc_r2 = 1 - sum((soc_actuals - soc_predictions).^2) / sum((soc_actuals - mean(soc_actuals)).^2);
            
            fprintf('SOC %d%%: N=%d, MAE=%.4f, RMSE=%.4f, MAPE=%.2f%%, R²=%.4f\n', ...
                soc, sum(idx), soc_mae, soc_rmse, soc_mape, soc_r2);
            
            soc_performance.(sprintf('SOC%d', soc)) = struct();
            soc_performance.(sprintf('SOC%d', soc)).n_samples = sum(idx);
            soc_performance.(sprintf('SOC%d', soc)).mae = soc_mae;
            soc_performance.(sprintf('SOC%d', soc)).rmse = soc_rmse;
            soc_performance.(sprintf('SOC%d', soc)).mape = soc_mape;
            soc_performance.(sprintf('SOC%d', soc)).r2 = soc_r2;
        end
    end
    
    % Visualize SOC-wise performance
    if length(uniqueSOCs) > 1
        figure('Position', [100, 100, 1400, 600]);
        
        % Subplot 1: RMSE by SOC
        subplot(1, 3, 1);
        soc_rmse_values = [];
        for s = 1:length(uniqueSOCs)
            soc = uniqueSOCs(s);
            if isfield(soc_performance, sprintf('SOC%d', soc))
                soc_rmse_values(s) = soc_performance.(sprintf('SOC%d', soc)).rmse;
            else
                soc_rmse_values(s) = NaN;
            end
        end
        bar(uniqueSOCs, soc_rmse_values);
        xlabel('SOC (%)', 'FontSize', 12);
        ylabel('RMSE (Ah)', 'FontSize', 12);
        title('RMSE by SOC Level', 'FontSize', 14);
        grid on;
        
        % Subplot 2: R² by SOC
        subplot(1, 3, 2);
        soc_r2_values = [];
        for s = 1:length(uniqueSOCs)
            soc = uniqueSOCs(s);
            if isfield(soc_performance, sprintf('SOC%d', soc))
                soc_r2_values(s) = soc_performance.(sprintf('SOC%d', soc)).r2;
            else
                soc_r2_values(s) = NaN;
            end
        end
        bar(uniqueSOCs, soc_r2_values);
        xlabel('SOC (%)', 'FontSize', 12);
        ylabel('R²', 'FontSize', 12);
        title('R² by SOC Level', 'FontSize', 14);
        ylim([min(0, min(soc_r2_values)-0.1), max(1, max(soc_r2_values)+0.1)]);
        grid on;
        
        % Subplot 3: Sample count by SOC
        subplot(1, 3, 3);
        soc_n_values = [];
        for s = 1:length(uniqueSOCs)
            soc = uniqueSOCs(s);
            if isfield(soc_performance, sprintf('SOC%d', soc))
                soc_n_values(s) = soc_performance.(sprintf('SOC%d', soc)).n_samples;
            else
                soc_n_values(s) = 0;
            end
        end
        bar(uniqueSOCs, soc_n_values);
        xlabel('SOC (%)', 'FontSize', 12);
        ylabel('Number of Samples', 'FontSize', 12);
        title('Sample Count by SOC Level', 'FontSize', 14);
        grid on;
        
        sgtitle('Model Performance by SOC Level', 'FontSize', 16);
        
        savefig(fullfile(figuresDir, 'XGBoost_Performance_BySOC.fig'));
        fprintf('Saved: XGBoost_Performance_BySOC.fig\n');
    end
else
    fprintf('WARNING: SOC data not available for SOC-wise analysis\n');
    fprintf('Make sure SOC column exists in test data\n');
    soc_performance = struct();
    all_soc_data = [];
end

%% Generate Standard Evaluation Plots
fprintf('\n=== Generating Standard Evaluation Plots ===\n');

% 전체 예측 결과 모음 (All Folds)
Y_true = all_actuals;
Y_pred = all_predictions;

% 1. Actual vs Estimated (Regression Plot)
fig1 = figure('Name', 'Actual vs Estimated', 'Position', [100, 100, 800, 700], 'Visible', 'on');
h1 = scatter(Y_true, Y_pred, 40, 'b', 'filled');
try
    h1.MarkerFaceAlpha = 0.5;
catch
end
hold on;
% 기준선 (Perfect Fit)
plot([min(Y_true), max(Y_true)], [min(Y_true), max(Y_true)], 'r--', 'LineWidth', 2);
% ±5% 오차선 (선택사항)
plot([min(Y_true), max(Y_true)], [min(Y_true)*0.95, max(Y_true)*0.95], 'k:', 'LineWidth', 1);
plot([min(Y_true), max(Y_true)], [min(Y_true)*1.05, max(Y_true)*1.05], 'k:', 'LineWidth', 1);
xlabel('Actual Capacity (Ah)', 'FontSize', 12);
ylabel('Estimated Capacity (Ah)', 'FontSize', 12);
title(sprintf('SOH Estimation Performance\nRMSE: %.3f Ah, R²: %.3f', overall_rmse, overall_r2), 'FontSize', 14);
grid on;
axis square;
savefig(fullfile(figuresDir, 'Actual_vs_Estimated.fig'));
close(fig1);
fprintf('  Saved: Actual_vs_Estimated.fig\n');

% Time-series Plot: Cycle vs Capacity (Actual vs Estimated)
fprintf('\n=== Generating Time-series Plot ===\n');

% Load Cycle and Channel information from test data
all_cycles = [];
all_channels = [];
for fold_num = 1:n_folds
    test_file = fullfile(inputDir, sprintf('DriveCycle_Fold%d_Test.mat', fold_num));
    if exist(test_file, 'file')
        load(test_file, 'testData');
        if istable(testData) && ismember('Cycle', testData.Properties.VariableNames) && ...
           ismember('Channel', testData.Properties.VariableNames)
            % Match indices with predictions
            fold_start = sum(all_fold_ids < fold_num) + 1;
            fold_end = sum(all_fold_ids <= fold_num);
            if fold_end >= fold_start && fold_end <= height(testData)
                fold_cycles = testData.Cycle(fold_start:fold_end);
                fold_channels = testData.Channel(fold_start:fold_end);
                all_cycles = [all_cycles; fold_cycles(:)];
                all_channels = [all_channels; fold_channels(:)];
            end
        end
    end
end

% If still empty, try to extract from fold_results
if isempty(all_cycles) || length(all_cycles) ~= length(all_actuals)
    all_cycles = [];
    all_channels = [];
    for fold_num = 1:n_folds
        fold_name = sprintf('Fold%d', fold_num);
        if isfield(fold_results, fold_name) && isfield(fold_results.(fold_name), 'testData')
            testData = fold_results.(fold_name).testData;
            if istable(testData) && ismember('Cycle', testData.Properties.VariableNames) && ...
               ismember('Channel', testData.Properties.VariableNames)
                valid_idx = ~isnan(fold_results.(fold_name).y_test);
                all_cycles = [all_cycles; testData.Cycle(valid_idx)];
                all_channels = [all_channels; testData.Channel(valid_idx)];
            end
        end
    end
end

% Create time-series plot if we have cycle and channel data
if ~isempty(all_cycles) && length(all_cycles) == length(all_actuals) && ...
   ~isempty(all_channels) && length(all_channels) == length(all_actuals)
    
    % Group by Channel and Cycle
    unique_channels = unique(all_channels);
    
    % Create figure with subplots for each channel (or combined)
    fig_ts = figure('Name', 'Time-series: Capacity vs Cycle', 'Position', [150, 150, 1600, 900], 'Visible', 'on');
    
    % If too many channels, create a combined plot with different colors
    if length(unique_channels) <= 8
        % Create subplot grid (2x4 for 8 channels max)
        nRows = 2;
        nCols = 4;
        
        for ch_idx = 1:length(unique_channels)
            ch = unique_channels(ch_idx);
            subplot(nRows, nCols, ch_idx);
            
            % Filter data for this channel
            ch_mask = (all_channels == ch);
            ch_cycles = all_cycles(ch_mask);
            ch_actuals = all_actuals(ch_mask);
            ch_predictions = all_predictions(ch_mask);
            
            % Sort by cycle
            [ch_cycles_sorted, sort_idx] = sort(ch_cycles);
            ch_actuals_sorted = ch_actuals(sort_idx);
            ch_predictions_sorted = ch_predictions(sort_idx);
            
            % Plot actual (black line)
            plot(ch_cycles_sorted, ch_actuals_sorted, 'k-', 'LineWidth', 2, 'DisplayName', 'Actual (RPT)');
            hold on;
            
            % Plot estimated (blue dots)
            plot(ch_cycles_sorted, ch_predictions_sorted, 'b.', 'MarkerSize', 12, 'DisplayName', 'Estimated');
            
            xlabel('Cycle', 'FontSize', 11, 'FontWeight', 'bold');
            ylabel('Capacity (Ah)', 'FontSize', 11, 'FontWeight', 'bold');
            title(sprintf('Channel %d', ch), 'FontSize', 12, 'FontWeight', 'bold');
            legend('Location', 'best', 'FontSize', 9);
            grid on;
        end
        
        sgtitle('SOH Estimation Time-series: Capacity vs Cycle (by Channel)', 'FontSize', 16, 'FontWeight', 'bold');
    else
        % Too many channels - create combined plot
        hold on;
        
        % Use different colors for different channels
        colors = lines(length(unique_channels));
        
        for ch_idx = 1:length(unique_channels)
            ch = unique_channels(ch_idx);
            
            % Filter data for this channel
            ch_mask = (all_channels == ch);
            ch_cycles = all_cycles(ch_mask);
            ch_actuals = all_actuals(ch_mask);
            ch_predictions = all_predictions(ch_mask);
            
            % Sort by cycle
            [ch_cycles_sorted, sort_idx] = sort(ch_cycles);
            ch_actuals_sorted = ch_actuals(sort_idx);
            ch_predictions_sorted = ch_predictions(sort_idx);
            
            % Plot actual (line)
            plot(ch_cycles_sorted, ch_actuals_sorted, '-', 'Color', colors(ch_idx, :), ...
                 'LineWidth', 1.5, 'DisplayName', sprintf('Ch%d Actual', ch));
            
            % Plot estimated (dots)
            plot(ch_cycles_sorted, ch_predictions_sorted, '.', 'Color', colors(ch_idx, :), ...
                 'MarkerSize', 10, 'DisplayName', sprintf('Ch%d Est', ch));
        end
        
        xlabel('Cycle', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Capacity (Ah)', 'FontSize', 12, 'FontWeight', 'bold');
        title('SOH Estimation Time-series: Capacity vs Cycle (All Channels)', 'FontSize', 14, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 9, 'NumColumns', 2);
        grid on;
    end
    
    savefig(fullfile(figuresDir, 'Time_Series_Capacity_vs_Cycle.fig'));
    close(fig_ts);
    fprintf('  Saved: Time_Series_Capacity_vs_Cycle.fig\n');
else
    fprintf('  WARNING: Cycle/Channel data not available for time-series plot\n');
    fprintf('  Make sure Cycle and Channel columns exist in test data\n');
end

% 2. Residual Plot (오차 분석)
% 오차가 특정 패턴을 가지면 모델에 문제가 있다는 뜻 (랜덤하게 퍼져야 좋음)
residuals = Y_true - Y_pred;
fig2 = figure('Name', 'Residual Analysis', 'Position', [150, 150, 1000, 500], 'Visible', 'on');

subplot(1, 2, 1);
h2 = scatter(Y_true, residuals, 40, 'k', 'filled');
try
    h2.MarkerFaceAlpha = 0.5;
catch
end
yline(0, 'r--', 'LineWidth', 2);
xlabel('Actual Capacity (Ah)', 'FontSize', 12);
ylabel('Residual Error (Ah)', 'FontSize', 12);
title('Residuals vs Actual', 'FontSize', 14);
grid on;

subplot(1, 2, 2);
histogram(residuals, 30);
xlabel('Residual Error (Ah)', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('Error Distribution', 'FontSize', 14);
grid on;

sgtitle('Model Error Analysis', 'FontSize', 16);
savefig(fullfile(figuresDir, 'Residual_Analysis.fig'));
close(fig2);
fprintf('  Saved: Residual_Analysis.fig\n');

%% Figure 1: Estimation vs Actual (All Folds Combined)
fprintf('\n=== Creating Additional Visualizations ===\n');

figure('Position', [100, 100, 1200, 800], 'Visible', 'on');

% Subplot 1: Estimation vs Actual (All data)
subplot(2, 2, 1);
% Use scatter without Alpha for compatibility (older MATLAB versions)
h1 = scatter(all_actuals, all_predictions, 50, 'filled');
try
    h1.MarkerFaceAlpha = 0.6;  % Try to set alpha if supported
catch
    % Alpha not supported, use default
end
hold on;
plot([min(all_actuals), max(all_actuals)], [min(all_actuals), max(all_actuals)], ...
    'r--', 'LineWidth', 2);
xlabel('Actual Capacity (Ah)', 'FontSize', 12);
ylabel('Estimated Capacity (Ah)', 'FontSize', 12);
title(sprintf('Estimation vs Actual (All Folds)\nR² = %.4f, RMSE = %.4f', ...
    overall_r2, overall_rmse), 'FontSize', 14);
grid on;
axis equal;

% Subplot 2: Residual Plot
subplot(2, 2, 2);
residuals = all_actuals - all_predictions;
h2 = scatter(all_predictions, residuals, 50, 'filled');
try
    h2.MarkerFaceAlpha = 0.6;  % Try to set alpha if supported
catch
    % Alpha not supported, use default
end
hold on;
plot([min(all_predictions), max(all_predictions)], [0, 0], ...
    'r--', 'LineWidth', 2);
xlabel('Estimated Capacity (Ah)', 'FontSize', 12);
ylabel('Residuals (Ah)', 'FontSize', 12);
title('Residual Plot', 'FontSize', 14);
grid on;

% Subplot 3: Fold-wise Performance Comparison
subplot(2, 2, 3);
fold_nums = 1:n_folds;
bar_data = [fold_mae, fold_rmse];
bar(fold_nums, bar_data);
xlabel('Fold', 'FontSize', 12);
ylabel('Error', 'FontSize', 12);
title('Fold-wise Performance', 'FontSize', 14);
legend('MAE', 'RMSE', 'Location', 'best');
grid on;
xticks(1:n_folds);

% Subplot 4: R² by Fold
subplot(2, 2, 4);
bar(fold_nums, fold_r2);
xlabel('Fold', 'FontSize', 12);
ylabel('R²', 'FontSize', 12);
title('R² by Fold', 'FontSize', 14);
ylim([0, 1]);
grid on;
xticks(1:n_folds);

sgtitle(sprintf('XGBoost Model Evaluation - %s', targetVar), 'FontSize', 16);

savefig(fullfile(figuresDir, 'XGBoost_Evaluation_Overall.fig'));
fprintf('Saved: XGBoost_Evaluation_Overall.fig\n');
close(gcf);

%% Figure 2: Estimation vs Actual (Individual Folds)
figure('Position', [100, 100, 1400, 1000]);

for fold_num = 1:n_folds
    subplot(2, 2, fold_num);
    
    if isfield(fold_results, sprintf('Fold%d', fold_num))
        fold_data = fold_results.(sprintf('Fold%d', fold_num));
        y_test = fold_data.y_test;
        y_test_pred = fold_data.y_test_pred;
        
        h = scatter(y_test, y_test_pred, 50, 'filled');
        try
            h.MarkerFaceAlpha = 0.6;  % Try to set alpha if supported
        catch
            % Alpha not supported, use default
        end
        hold on;
        plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], ...
            'r--', 'LineWidth', 2);
        xlabel('Actual Capacity (Ah)', 'FontSize', 11);
        ylabel('Estimated Capacity (Ah)', 'FontSize', 11);
        title(sprintf('Fold %d (R² = %.4f, RMSE = %.4f)', ...
            fold_num, fold_data.test_r2, fold_data.test_rmse), 'FontSize', 12);
        grid on;
        axis equal;
    end
end

sgtitle('Estimation vs Actual - Individual Folds', 'FontSize', 16);

savefig(fullfile(figuresDir, 'XGBoost_Evaluation_ByFold.fig'));
fprintf('Saved: XGBoost_Evaluation_ByFold.fig\n');
close(gcf);

%% Figure 3: Feature Importance (Average across folds)
fprintf('\n=== Feature Importance Analysis ===\n');

% Collect feature importance from all folds
all_importance = [];
all_feature_names = {};

for fold_num = 1:n_folds
    if isfield(fold_results, sprintf('Fold%d', fold_num))
        fold_data = fold_results.(sprintf('Fold%d', fold_num));
        if isfield(fold_data, 'feature_importance') && ~isempty(fold_data.feature_importance)
            all_importance = [all_importance; fold_data.feature_importance(:)'];
            if isempty(all_feature_names)
                all_feature_names = fold_data.available_features;
            end
        end
    end
end

if ~isempty(all_importance)
    % Average importance across folds
    avg_importance = mean(all_importance, 1);
    
    % Sort by importance
    [sorted_importance, sort_idx] = sort(avg_importance, 'descend');
    sorted_features = all_feature_names(sort_idx);
    
    % Top 20 features
    n_top = min(20, length(sorted_features));
    top_features = sorted_features(1:n_top);
    top_importance = sorted_importance(1:n_top);
    
    figure('Position', [100, 100, 1200, 800]);
    barh(flipud(top_importance));
    set(gca, 'YTickLabel', flipud(top_features), 'FontSize', 10);
    xlabel('Feature Importance', 'FontSize', 12);
    ylabel('Features', 'FontSize', 12);
    title(sprintf('Top %d Feature Importance (Average across Folds)', n_top), 'FontSize', 14);
    grid on;
    
    savefig(fullfile(figuresDir, 'XGBoost_FeatureImportance.fig'));    
    fprintf('Saved: XGBoost_FeatureImportance.fig\n');
    close(gcf);
    
    % Display top features
    fprintf('\nTop 10 Most Important Features:\n');
    for i = 1:min(10, length(top_features))
        fprintf('  %d. %s (importance: %.4f)\n', i, top_features{i}, top_importance(i));
    end
end

%% Generate Report
fprintf('\n=== Generating Report ===\n');

report_file = fullfile(outputDir, 'XGBoost_Evaluation_Report.txt');
fid = fopen(report_file, 'w');

fprintf(fid, '========================================\n');
fprintf(fid, 'XGBoost Model Evaluation Report\n');
fprintf(fid, '========================================\n\n');

fprintf(fid, 'Target Variable: %s\n', targetVar);
fprintf(fid, 'Number of Folds: %d\n', n_folds);
fprintf(fid, 'Date: %s\n\n', datestr(now));

fprintf(fid, 'Cross-Validation Performance (Mean ± Std):\n');
fprintf(fid, '  MAE:  %.4f ± %.4f\n', mean(fold_mae), std(fold_mae));
fprintf(fid, '  RMSE: %.4f ± %.4f\n', mean(fold_rmse), std(fold_rmse));
fprintf(fid, '  MAPE: %.2f%% ± %.2f%%\n', mean(fold_mape), std(fold_mape));
fprintf(fid, '  R²:   %.4f ± %.4f\n\n', mean(fold_r2), std(fold_r2));

fprintf(fid, 'Overall Performance (All Test Data Combined):\n');
fprintf(fid, '  MAE:  %.4f\n', overall_mae);
fprintf(fid, '  RMSE: %.4f\n', overall_rmse);
fprintf(fid, '  MAPE: %.2f%%\n', overall_mape);
fprintf(fid, '  R²:   %.4f\n\n', overall_r2);

fprintf(fid, 'Fold-wise Performance:\n');
for fold_num = 1:n_folds
    if isfield(fold_results, sprintf('Fold%d', fold_num))
        fold_data = fold_results.(sprintf('Fold%d', fold_num));
        fprintf(fid, '  Fold %d: MAE=%.4f, RMSE=%.4f, MAPE=%.2f%%, R²=%.4f\n', ...
            fold_num, fold_data.test_mae, fold_data.test_rmse, ...
            fold_data.test_mape, fold_data.test_r2);
    end
end

if ~isempty(all_importance)
    fprintf(fid, '\nTop 10 Most Important Features:\n');
    for i = 1:min(10, length(top_features))
        fprintf(fid, '  %d. %s (importance: %.4f)\n', i, top_features{i}, top_importance(i));
    end
end

% Add SOC-wise performance to report
if exist('soc_performance', 'var') && ~isempty(fieldnames(soc_performance))
    fprintf(fid, '\n=== Performance by SOC Level ===\n');
    uniqueSOCs = [];
    for s = 1:100
        if isfield(soc_performance, sprintf('SOC%d', s))
            uniqueSOCs(end+1) = s;
        end
    end
    uniqueSOCs = sort(uniqueSOCs);
    
    for s = 1:length(uniqueSOCs)
        soc = uniqueSOCs(s);
        if isfield(soc_performance, sprintf('SOC%d', soc))
            perf = soc_performance.(sprintf('SOC%d', soc));
            fprintf(fid, 'SOC %d%%: N=%d, MAE=%.4f, RMSE=%.4f, MAPE=%.2f%%, R²=%.4f\n', ...
                soc, perf.n_samples, perf.mae, perf.rmse, perf.mape, perf.r2);
        end
    end
end

fclose(fid);
fprintf('Saved report: %s\n', report_file);

fprintf('\n=== Evaluation Complete ===\n');
fprintf('All figures saved to: %s\n', figuresDir);

