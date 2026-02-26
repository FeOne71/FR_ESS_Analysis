%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT_Modeling.m
% PCA + MLR & Random Forest Modeling with 5-Fold CV
%
% 목적: 14개 배터리 피처를 PCA로 차원 축소 후,
%       SOH/LLI/LAM을 추정하는 MLR 및 RF 모델을 구축/검증한다.
%
% 검증 방식: 5-Fold Cross-Validation (셀 혼합)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
warning off;

%% ========================================================================
% Section 1: 데이터 로드 & 전처리
% ========================================================================
fprintf('=== Section 1: 데이터 로드 & 전처리 ===\n');

% --- 경로 설정 ---
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis';
dataPath = fullfile(baseDir, 'Dataset', 'Feature_Matrix_Final.mat');
saveDir = fullfile(baseDir, 'Modeling_Results');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

% --- 데이터 로드 ---
fprintf('Loading Feature Matrix...\n');
load(dataPath, 'FeatureTable');
fprintf('  Total samples: %d\n', height(FeatureTable));

% --- 입력/라벨 분리 ---
X = FeatureTable.X_Normalized;   % 14개 피처 (이미 C-rate별 정규화됨)
Y = FeatureTable.Y_Labels;       % [SOH, LLI, LAM]
cellIDs = FeatureTable.CellID;   % 셀 ID (그룹 변수)

% --- 그룹 정보 ---
unique_cells = unique(cellIDs);
n_cells = length(unique_cells);
n_features = size(X, 2);
n_labels = size(Y, 2);
label_names = {'SOH', 'LLI', 'LAM'};

% --- NaN 제거 ---
nan_rows = any(isnan(X), 2) | any(isnan(Y), 2);
if sum(nan_rows) > 0
    fprintf('  Removing %d rows with NaN values.\n', sum(nan_rows));
    X(nan_rows, :) = [];
    Y(nan_rows, :) = [];
    cellIDs(nan_rows) = [];
end

% --- 데이터 요약 출력 ---
fprintf('\n--- 데이터 요약 ---\n');
fprintf('  유효 샘플 수: %d\n', size(X, 1));
fprintf('  피처 수: %d\n', n_features);
fprintf('  라벨 수: %d (%s)\n', n_labels, strjoin(label_names, ', '));
fprintf('  셀 그룹: %d개 (%s)\n', n_cells, strjoin(unique_cells', ', '));

% 셀별 샘플 수
fprintf('\n  셀별 샘플 수:\n');
for i = 1:n_cells
    mask = strcmp(cellIDs, unique_cells{i});
    fprintf('    %s: %d samples\n', unique_cells{i}, sum(mask));
end

fprintf('\n=== Section 1 완료 ===\n');

%% ========================================================================
% Section 2: PCA (Scree Plot, 95% 기준 PC 선택)
% ========================================================================
fprintf('\n=== Section 2: PCA 분석 ===\n');

% --- PCA 수행 (전체 데이터) ---
[coeff, score, latent, ~, explained] = pca(X);
cumulative_var = cumsum(explained);

% --- 95% 누적 설명력 기준 PC 개수 결정 ---
threshold = 95;
n_pc = find(cumulative_var >= threshold, 1, 'first');
fprintf('  누적 설명력 %.0f%% 기준 → PC 개수: %d개 (총 14개 중)\n', threshold, n_pc);
fprintf('  선택된 PC의 누적 설명력: %.2f%%\n', cumulative_var(n_pc));

% --- 개별 PC 설명력 출력 ---
fprintf('\n  PC별 설명력:\n');
for k = 1:n_pc
    fprintf('    PC%d: %.2f%% (누적: %.2f%%)\n', k, explained(k), cumulative_var(k));
end

% --- Scree Plot 생성 ---
fig_scree = figure('Position', [100, 100, 800, 400], 'Name', 'PCA Scree Plot');

subplot(1, 2, 1);
bar(explained, 'FaceColor', [0.3 0.6 0.9]);
hold on;
bar(1:n_pc, explained(1:n_pc), 'FaceColor', [0.2 0.4 0.8]);
xlabel('Principal Component', 'FontSize', 12);
ylabel('Variance Explained (%)', 'FontSize', 12);
title('Individual Variance', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
hold off;

subplot(1, 2, 2);
plot(cumulative_var, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0.2 0.4 0.8]);
hold on;
yline(threshold, '--r', sprintf('%.0f%%', threshold), 'LineWidth', 1.5, 'FontSize', 11);
xline(n_pc, '--k', sprintf('PC%d', n_pc), 'LineWidth', 1.5, 'FontSize', 11);
plot(n_pc, cumulative_var(n_pc), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Number of PCs', 'FontSize', 12);
ylabel('Cumulative Variance (%)', 'FontSize', 12);
title('Cumulative Variance', 'FontSize', 13, 'FontWeight', 'bold');
ylim([0 100]);
grid on;
hold off;

sgtitle('PCA Scree Plot (14 Features)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig_scree, fullfile(saveDir, 'PCA_Scree_Plot.fig'));
fprintf('  Scree Plot 저장 완료\n');

% --- PC score 저장 ---
X_pc = score(:, 1:n_pc);

fprintf('\n=== Section 2 완료 ===\n');

%% ========================================================================
% Section 3: [1단계] 전체 데이터 MLR Sanity Check
% ========================================================================
fprintf('\n=== Section 3: 전체 데이터 MLR Sanity Check ===\n');
fprintf('  (Train = Test, 과적합 허용 → 피처가 라벨을 설명 가능한지 확인)\n\n');

R2_sanity = zeros(1, n_labels);
Y_pred_sanity = zeros(size(Y));

for lbl = 1:n_labels
    y = Y(:, lbl);
    mdl = fitlm(X_pc, y);
    Y_pred_sanity(:, lbl) = predict(mdl, X_pc);
    R2_sanity(lbl) = mdl.Rsquared.Ordinary;
    fprintf('  %s: R² = %.4f (Adjusted R² = %.4f)\n', ...
        label_names{lbl}, mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted);
end

% --- Predicted vs Actual 산점도 ---
fig_sanity = figure('Position', [50, 50, 1500, 450], 'Name', 'MLR Sanity Check');
label_units_plot = {'Capacity (Ah)', 'LLI (%)', 'LAM (%)'};
colors_cell = lines(n_cells);

for lbl = 1:n_labels
    subplot(1, 3, lbl);
    hold on;
    for ci = 1:n_cells
        mask_ci = strcmp(cellIDs, unique_cells{ci});
        scatter(Y(mask_ci, lbl), Y_pred_sanity(mask_ci, lbl), 30, ...
            colors_cell(ci,:), 'filled', 'MarkerFaceAlpha', 0.7, ...
            'DisplayName', unique_cells{ci});
    end
    ax_min = min([Y(:,lbl); Y_pred_sanity(:,lbl)]);
    ax_max = max([Y(:,lbl); Y_pred_sanity(:,lbl)]);
    margin = (ax_max - ax_min) * 0.05;
    plot([ax_min-margin, ax_max+margin], [ax_min-margin, ax_max+margin], ...
        '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel(sprintf('Actual %s', label_units_plot{lbl}), 'FontSize', 11);
    ylabel(sprintf('Predicted %s', label_units_plot{lbl}), 'FontSize', 11);
    title(sprintf('%s (R² = %.4f)', label_names{lbl}, R2_sanity(lbl)), ...
        'FontSize', 13, 'FontWeight', 'bold');
    xlim([ax_min-margin, ax_max+margin]);
    ylim([ax_min-margin, ax_max+margin]);
    grid on; axis square;
    if lbl == 1, legend('Location', 'southeast', 'FontSize', 7); end
    hold off;
end

sgtitle(sprintf('MLR Sanity Check (All Data, %d PCs)', n_pc), ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig_sanity, fullfile(saveDir, 'MLR_Sanity_Check.fig'));

fprintf('\n=== Section 3 완료 ===\n');

%% ========================================================================
% Section 4: 5-Fold CV (MLR + RF + HP 튜닝)
% ========================================================================
fprintf('\n=== Section 4: 5-Fold CV (셀 혼합) ===\n');
fprintf('  8개 셀 데이터를 섞어서 랜덤 5-Fold 분할\n\n');

K = 5;
rng(42);
n_samples = size(X, 1);
cv_indices = crossvalind('Kfold', n_samples, K);

% 결과 저장
metrics_MLR = zeros(K, n_labels, 3);  % (fold, label, [RMSE MAE R2])
metrics_RF  = zeros(K, n_labels, 3);
kfold_Y_true = cell(K, 1);
kfold_Y_pred_MLR = cell(K, 1);
kfold_Y_pred_RF  = cell(K, 1);

% HP 그리드
hp_numTrees = [50, 100, 200, 500];
hp_minLeaf  = [1, 3, 5, 10];
best_hp_log = zeros(K, n_labels, 2);

for fold = 1:K
    fprintf('  [Fold %d/%d] Train: %d, Test: %d samples\n', ...
        fold, K, sum(cv_indices ~= fold), sum(cv_indices == fold));
    
    test_mask  = (cv_indices == fold);
    train_mask = ~test_mask;
    
    X_train_raw = X(train_mask, :);
    X_test_raw  = X(test_mask, :);
    Y_train = Y(train_mask, :);
    Y_test  = Y(test_mask, :);
    
    % Z-score (Train 기준)
    mu_tr  = mean(X_train_raw);
    std_tr = std(X_train_raw);
    std_tr(std_tr == 0) = 1;
    X_train_z = (X_train_raw - mu_tr) ./ std_tr;
    X_test_z  = (X_test_raw  - mu_tr) ./ std_tr;
    
    % PCA (Train 축)
    [coeff_k, score_k, ~, ~, expl_k] = pca(X_train_z);
    cumvar_k = cumsum(expl_k);
    n_pc_k = find(cumvar_k >= threshold, 1, 'first');
    
    X_train_pc = score_k(:, 1:n_pc_k);
    X_test_pc  = X_test_z * coeff_k(:, 1:n_pc_k);
    
    kfold_Y_true{fold} = Y_test;
    kfold_Y_pred_MLR{fold} = zeros(size(Y_test));
    kfold_Y_pred_RF{fold}  = zeros(size(Y_test));
    
    for lbl = 1:n_labels
        y_train = Y_train(:, lbl);
        y_test  = Y_test(:, lbl);
        
        % --- MLR ---
        mdl_mlr = fitlm(X_train_pc, y_train);
        y_pred_mlr = predict(mdl_mlr, X_test_pc);
        
        metrics_MLR(fold, lbl, 1) = sqrt(mean((y_test - y_pred_mlr).^2));
        metrics_MLR(fold, lbl, 2) = mean(abs(y_test - y_pred_mlr));
        ss_res = sum((y_test - y_pred_mlr).^2);
        ss_tot = sum((y_test - mean(y_test)).^2);
        if ss_tot == 0, metrics_MLR(fold, lbl, 3) = NaN;
        else, metrics_MLR(fold, lbl, 3) = 1 - ss_res / ss_tot; end
        
        % --- RF: HP 그리드 서치 (OOB 기반) ---
        best_oob = Inf;
        best_nt = hp_numTrees(1);
        best_ml = hp_minLeaf(1);
        
        for nt = hp_numTrees
            for ml = hp_minLeaf
                t = templateTree('MinLeafSize', ml);
                rf_temp = fitrensemble(X_train_pc, y_train, ...
                    'Method', 'Bag', ...
                    'NumLearningCycles', nt, ...
                    'Learners', t);
                y_oob = oobPredict(rf_temp);
                oob_rmse = sqrt(mean((y_train - y_oob).^2));
                
                if oob_rmse < best_oob
                    best_oob = oob_rmse;
                    best_nt = nt;
                    best_ml = ml;
                end
            end
        end
        
        best_hp_log(fold, lbl, :) = [best_nt, best_ml];
        
        % 최적 HP로 최종 RF 학습 → Test 예측
        t_best = templateTree('MinLeafSize', best_ml);
        rf_final = fitrensemble(X_train_pc, y_train, ...
            'Method', 'Bag', ...
            'NumLearningCycles', best_nt, ...
            'Learners', t_best);
        y_pred_rf = predict(rf_final, X_test_pc);
        
        metrics_RF(fold, lbl, 1) = sqrt(mean((y_test - y_pred_rf).^2));
        metrics_RF(fold, lbl, 2) = mean(abs(y_test - y_pred_rf));
        ss_res_rf = sum((y_test - y_pred_rf).^2);
        if ss_tot == 0, metrics_RF(fold, lbl, 3) = NaN;
        else, metrics_RF(fold, lbl, 3) = 1 - ss_res_rf / ss_tot; end
        
        kfold_Y_pred_MLR{fold}(:, lbl) = y_pred_mlr;
        kfold_Y_pred_RF{fold}(:, lbl)  = y_pred_rf;
    end
    
    fprintf('    PCs: %d | RF HP(SOH): [Trees=%d, MinLeaf=%d]\n', ...
        n_pc_k, best_hp_log(fold, 1, 1), best_hp_log(fold, 1, 2));
    fprintf('    MLR RMSE(SOH)=%.4f | RF RMSE(SOH)=%.4f\n', ...
        metrics_MLR(fold, 1, 1), metrics_RF(fold, 1, 1));
end

fprintf('\n=== Section 4 완료 ===\n');

%% ========================================================================
% Section 5: 성능 요약 테이블 출력
% ========================================================================
fprintf('\n=== Section 5: 5-Fold CV 성능 요약 ===\n');

for lbl = 1:n_labels
    fprintf('\n────────────────────────────────────────────\n');
    fprintf(' [%s] Performance Summary\n', label_names{lbl});
    fprintf('────────────────────────────────────────────\n');
    fprintf('%-8s | %12s %12s %12s | %12s %12s %12s\n', ...
        'Fold', 'MLR_RMSE', 'MLR_MAE', 'MLR_R²', ...
        'RF_RMSE', 'RF_MAE', 'RF_R²');
    fprintf('%s\n', repmat('-', 1, 86));
    
    for fold = 1:K
        fprintf('Fold %d   | %12.4f %12.4f %12.4f | %12.4f %12.4f %12.4f\n', ...
            fold, ...
            metrics_MLR(fold, lbl, 1), metrics_MLR(fold, lbl, 2), metrics_MLR(fold, lbl, 3), ...
            metrics_RF(fold, lbl, 1),  metrics_RF(fold, lbl, 2),  metrics_RF(fold, lbl, 3));
    end
    
    fprintf('%s\n', repmat('-', 1, 86));
    fprintf('%-8s | %12.4f %12.4f %12.4f | %12.4f %12.4f %12.4f\n', ...
        'Mean', ...
        mean(metrics_MLR(:, lbl, 1)), mean(metrics_MLR(:, lbl, 2)), mean(metrics_MLR(:, lbl, 3)), ...
        mean(metrics_RF(:, lbl, 1)),  mean(metrics_RF(:, lbl, 2)),  mean(metrics_RF(:, lbl, 3)));
    fprintf('%-8s | %12.4f %12.4f %12.4f | %12.4f %12.4f %12.4f\n', ...
        'Std', ...
        std(metrics_MLR(:, lbl, 1)), std(metrics_MLR(:, lbl, 2)), std(metrics_MLR(:, lbl, 3)), ...
        std(metrics_RF(:, lbl, 1)),  std(metrics_RF(:, lbl, 2)),  std(metrics_RF(:, lbl, 3)));
end

fprintf('\n=== Section 5 완료 ===\n');

%% ========================================================================
% Section 6: 시각화 (분석 + PPT + 프로세스 증명)
% ========================================================================
fprintf('\n=== Section 6: 시각화 ===\n');

% PPT 폰트 설정
ppt_fontsize_title = 18;
ppt_fontsize_axis  = 14;
ppt_fontsize_tick  = 12;
ppt_fontsize_value = 10;

% ---------- 6.1: Predicted vs Actual (RF) ----------
fprintf('  6.1: Predicted vs Actual (RF)...\n');
fig_scatter = figure('Position', [50, 50, 1500, 450], 'Name', 'K-Fold Predicted vs Actual');

for lbl = 1:n_labels
    subplot(1, 3, lbl);
    hold on;
    
    all_true = vertcat(kfold_Y_true{:});
    all_pred = vertcat(kfold_Y_pred_RF{:});
    
    scatter(all_true(:, lbl), all_pred(:, lbl), 25, ...
        [0.85 0.33 0.1], 'filled', 'MarkerFaceAlpha', 0.6);
    
    ax_min = min([all_true(:, lbl); all_pred(:, lbl)]);
    ax_max = max([all_true(:, lbl); all_pred(:, lbl)]);
    margin = (ax_max - ax_min) * 0.05;
    plot([ax_min-margin, ax_max+margin], [ax_min-margin, ax_max+margin], ...
        '--k', 'LineWidth', 1.5);
    
    ss_res_all = sum((all_true(:,lbl) - all_pred(:,lbl)).^2);
    ss_tot_all = sum((all_true(:,lbl) - mean(all_true(:,lbl))).^2);
    r2_all = 1 - ss_res_all / ss_tot_all;
    
    xlabel('Actual', 'FontSize', ppt_fontsize_axis);
    ylabel('Predicted (RF)', 'FontSize', ppt_fontsize_axis);
    title(sprintf('%s (R²=%.4f)', label_names{lbl}, r2_all), ...
        'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
    xlim([ax_min-margin, ax_max+margin]);
    ylim([ax_min-margin, ax_max+margin]);
    grid on; axis square;
    hold off;
end
sgtitle('5-Fold CV: RF Predicted vs Actual', 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
saveas(fig_scatter, fullfile(saveDir, 'KFold_RF_Predicted_vs_Actual.fig'));

% ---------- 6.2: Feature Importance ----------
fprintf('  6.2: Feature Importance...\n');
feature_names = {'Chg_dQ1','Chg_dQ2','Chg_dQ3','Chg_dQ4','Chg_dQ5', ...
                 'Dch_dQ1','Dch_dQ2','Dch_dQ3','Dch_dQ4','Dch_dQ5', ...
                 'Chg_PkH','Chg_PkA','Dch_PkH','Dch_PkA'};

fig_fi = figure('Position', [50, 50, 1500, 450], 'Name', 'Feature Importance');
for lbl = 1:n_labels
    t_fi = templateTree('MinLeafSize', 1);
    rf_fi = fitrensemble(X, Y(:, lbl), ...
        'Method', 'Bag', 'NumLearningCycles', 200, 'Learners', t_fi);
    imp = predictorImportance(rf_fi);
    
    subplot(1, 3, lbl);
    [imp_sorted, sort_idx] = sort(imp, 'descend');
    barh(imp_sorted, 'FaceColor', [0.3 0.6 0.9]);
    set(gca, 'YTick', 1:n_features, 'YTickLabel', feature_names(sort_idx), ...
        'YDir', 'reverse', 'FontSize', 9);
    xlabel('Importance', 'FontSize', 11);
    title(sprintf('%s', label_names{lbl}), 'FontSize', 13, 'FontWeight', 'bold');
    grid on;
end
sgtitle('RF Feature Importance (All Data)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig_fi, fullfile(saveDir, 'RF_Feature_Importance.fig'));

% ---------- 6.3: Error Distribution Histogram ----------
fprintf('  6.3: Error Distribution...\n');
fig_errhist = figure('Position', [50, 50, 1500, 450], 'Name', 'Error Distribution');
for lbl = 1:n_labels
    all_errors = [];
    for fold = 1:K
        err = kfold_Y_pred_RF{fold}(:, lbl) - kfold_Y_true{fold}(:, lbl);
        all_errors = [all_errors; err]; %#ok<AGROW>
    end
    subplot(1, 3, lbl);
    histogram(all_errors, 25, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'w', 'FaceAlpha', 0.85);
    hold on;
    xline(0, '--r', 'LineWidth', 2);
    xline(mean(all_errors), '-k', sprintf('Mean=%.4f', mean(all_errors)), ...
        'LineWidth', 1.5, 'FontSize', 10, 'LabelOrientation', 'horizontal');
    xlabel('Prediction Error', 'FontSize', 11);
    ylabel('Count', 'FontSize', 11);
    title(sprintf('%s (Std=%.4f)', label_names{lbl}, std(all_errors)), ...
        'FontSize', 13, 'FontWeight', 'bold');
    grid on; hold off;
end
sgtitle('5-Fold CV: RF Prediction Error Distribution', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig_errhist, fullfile(saveDir, 'RF_Error_Distribution.fig'));

% ---------- 6.4: MLR vs RF SOH RMSE Grouped Bar ----------
fprintf('  6.4: MLR vs RF SOH 비교...\n');
cell_labels_fold = arrayfun(@(x) sprintf('Fold %d', x), 1:K, 'UniformOutput', false);
soh_rmse_mlr = metrics_MLR(:, 1, 1);
soh_rmse_rf  = metrics_RF(:, 1, 1);

fig_comp = figure('Position', [50, 50, 900, 500], 'Name', 'MLR vs RF SOH');
comparison_data = [soh_rmse_mlr, soh_rmse_rf];
b = bar(comparison_data, 'grouped');
b(1).FaceColor = [0.2 0.5 0.85]; b(1).EdgeColor = 'none';
b(2).FaceColor = [0.85 0.33 0.1]; b(2).EdgeColor = 'none';
hold on;
yline(mean(soh_rmse_mlr), '--', sprintf('MLR Mean=%.3f', mean(soh_rmse_mlr)), ...
    'Color', [0.2 0.5 0.85], 'LineWidth', 1.5, 'FontSize', ppt_fontsize_value);
yline(mean(soh_rmse_rf), '--', sprintf('RF Mean=%.3f', mean(soh_rmse_rf)), ...
    'Color', [0.85 0.33 0.1], 'LineWidth', 1.5, 'FontSize', ppt_fontsize_value);
set(gca, 'XTickLabel', cell_labels_fold, 'FontSize', ppt_fontsize_tick);
xlabel('Fold', 'FontSize', ppt_fontsize_axis);
ylabel('RMSE (Ah)', 'FontSize', ppt_fontsize_axis);
title('SOH: MLR vs Random Forest', 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
legend({'MLR', 'Random Forest'}, 'FontSize', ppt_fontsize_tick, 'Location', 'northwest');
grid on; hold off;
saveas(fig_comp, fullfile(saveDir, 'MLR_vs_RF_SOH_Comparison.fig'));

% ---------- 6.5: Overall Label Comparison ----------
fprintf('  6.5: Overall Label 비교...\n');
fig_overall = figure('Position', [50, 50, 900, 500], 'Name', 'Overall Label Comparison');
label_units_ppt = {'RMSE (Ah)', 'RMSE (V)', 'RMSE (Ah/V)'};
bar_colors = {[0.2 0.5 0.85], [0.85 0.33 0.1]};

for lbl = 1:n_labels
    subplot(1, 3, lbl);
    mean_rmse = [mean(metrics_MLR(:, lbl, 1)), mean(metrics_RF(:, lbl, 1))];
    b = bar(mean_rmse, 0.6);
    b.FaceColor = 'flat';
    b.CData(1,:) = bar_colors{1}; b.CData(2,:) = bar_colors{2};
    b.EdgeColor = 'none';
    hold on;
    text(1, mean_rmse(1)*1.05, sprintf('%.4f', mean_rmse(1)), ...
        'HorizontalAlignment', 'center', 'FontSize', ppt_fontsize_tick, 'FontWeight', 'bold');
    text(2, mean_rmse(2)*1.05, sprintf('%.4f', mean_rmse(2)), ...
        'HorizontalAlignment', 'center', 'FontSize', ppt_fontsize_tick, 'FontWeight', 'bold');
    set(gca, 'XTickLabel', {'MLR', 'RF'}, 'FontSize', ppt_fontsize_tick);
    ylabel(label_units_ppt{lbl}, 'FontSize', ppt_fontsize_axis);
    title(label_names{lbl}, 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
    grid on; hold off;
end
sgtitle('5-Fold CV Mean RMSE: MLR vs RF', 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
saveas(fig_overall, fullfile(saveDir, 'Overall_Label_Comparison.fig'));

% ---------- 6.6: HP 튜닝 히트맵 ----------
fprintf('  6.6: HP 튜닝 히트맵 (All Labels)...\n');
hp_trees_grid = [50, 100, 200, 500];
hp_leaf_grid  = [1, 3, 5, 10];
rmse_heatmap  = zeros(length(hp_leaf_grid), length(hp_trees_grid), n_labels);

% PCA for HP Tuning (Using all data once)
mu_all = mean(X); std_all = std(X); std_all(std_all==0) = 1;
X_z_all = (X - mu_all) ./ std_all;
[coeff_all, score_all, ~, ~, expl_all] = pca(X_z_all);
cumvar_all = cumsum(expl_all);
n_pc_all = find(cumvar_all >= threshold, 1, 'first');
X_pc_all = score_all(:, 1:n_pc_all);

for lbl = 1:n_labels
    fprintf('    Processing Label: %s...\n', label_names{lbl});
    for ri = 1:length(hp_leaf_grid)
        for ci = 1:length(hp_trees_grid)
            t_hp = templateTree('MinLeafSize', hp_leaf_grid(ri));
            rf_hp = fitrensemble(X_pc_all, Y(:, lbl), ...
                'Method', 'Bag', 'NumLearningCycles', hp_trees_grid(ci), 'Learners', t_hp);
            y_oob = oobPredict(rf_hp);
            rmse_heatmap(ri, ci, lbl) = sqrt(mean((Y(:, lbl) - y_oob).^2));
        end
    end
    
    % 개별 라벨 시각화 (Optional, figure만 생성하되 대표로 하나만 표시하거나 다 저장)
    fig_hp = figure('Position', [50, 50, 700, 550], 'Name', ['HP Tuning Heatmap: ' label_names{lbl}], 'Visible', 'off');
    imagesc(rmse_heatmap(:,:,lbl)); colormap(flipud(hot));
    cb = colorbar; cb.Label.String = sprintf('%s OOB RMSE', label_names{lbl}); cb.Label.FontSize = ppt_fontsize_tick;
    set(gca, 'XTick', 1:length(hp_trees_grid), 'XTickLabel', hp_trees_grid, ...
        'YTick', 1:length(hp_leaf_grid), 'YTickLabel', hp_leaf_grid, 'FontSize', ppt_fontsize_tick);
    xlabel('Number of Trees', 'FontSize', ppt_fontsize_axis);
    ylabel('Min Leaf Size', 'FontSize', ppt_fontsize_axis);
    title(['RF Hyperparameter Tuning: ' label_names{lbl}], 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
    
    [~, min_idx] = min(reshape(rmse_heatmap(:,:,lbl), [], 1));
    [min_r, min_c] = ind2sub([length(hp_leaf_grid), length(hp_trees_grid)], min_idx);
    hold on;
    plot(min_c, min_r, 'wo', 'MarkerSize', 20, 'LineWidth', 3);
    for ri = 1:length(hp_leaf_grid)
        for ci = 1:length(hp_trees_grid)
            if ri == min_r && ci == min_c
                text(ci, ri, sprintf('%.3f', rmse_heatmap(ri,ci,lbl)), ...
                    'Color', 'w', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
            else
                text(ci, ri, sprintf('%.3f', rmse_heatmap(ri,ci,lbl)), ...
                    'Color', 'k', 'HorizontalAlignment', 'center');
            end
        end
    end
    hold off;
    saveas(fig_hp, fullfile(saveDir, ['HP_Tuning_Heatmap_' label_names{lbl} '.fig']));
end

% ---------- 6.7: Residual Plot ----------
fprintf('  6.7: Residual Plot...\n');
fig_resid = figure('Position', [50, 50, 1500, 450], 'Name', 'Residual Plot');
for lbl = 1:n_labels
    subplot(1, 3, lbl);
    all_true_r = vertcat(kfold_Y_true{:});
    all_pred_r = vertcat(kfold_Y_pred_RF{:});
    residuals = all_pred_r(:, lbl) - all_true_r(:, lbl);
    
    scatter(all_true_r(:, lbl), residuals, 25, [0.2 0.5 0.85], 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    yline(0, '--r', 'LineWidth', 2);
    sigma = std(residuals);
    yline(sigma, ':k', sprintf('+1σ=%.4f', sigma), 'FontSize', 9);
    yline(-sigma, ':k', sprintf('-1σ=%.4f', -sigma), 'FontSize', 9);
    xlabel('Actual Value', 'FontSize', ppt_fontsize_axis);
    ylabel('Residual (Pred - Actual)', 'FontSize', ppt_fontsize_axis);
    title(sprintf('%s (Bias=%.4f)', label_names{lbl}, mean(residuals)), ...
        'FontSize', ppt_fontsize_title - 2, 'FontWeight', 'bold');
    grid on; hold off;
end
sgtitle('5-Fold CV: RF Residual Analysis', 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
saveas(fig_resid, fullfile(saveDir, 'RF_Residual_Plot.fig'));

% ---------- 6.8: MLR vs RF R² 개선 비교 ----------
fprintf('  6.8: MLR vs RF R² 비교...\n');
fig_r2comp = figure('Position', [50, 50, 800, 500], 'Name', 'MLR vs RF R2');
r2_mlr_kf = squeeze(mean(metrics_MLR(:, :, 3), 1));
r2_rf_kf  = squeeze(mean(metrics_RF(:, :, 3), 1));
r2_data = [r2_mlr_kf; r2_rf_kf]';

b = bar(r2_data, 'grouped');
b(1).FaceColor = [0.2 0.5 0.85]; b(1).EdgeColor = 'none';
b(2).FaceColor = [0.85 0.33 0.1]; b(2).EdgeColor = 'none';
hold on;
for lbl = 1:n_labels
    x_mlr = b(1).XEndPoints(lbl);
    x_rf  = b(2).XEndPoints(lbl);
    text(x_mlr, r2_mlr_kf(lbl)+0.02, sprintf('%.3f', r2_mlr_kf(lbl)), ...
        'HorizontalAlignment', 'center', 'FontSize', ppt_fontsize_value, 'FontWeight', 'bold');
    text(x_rf, r2_rf_kf(lbl)+0.02, sprintf('%.3f', r2_rf_kf(lbl)), ...
        'HorizontalAlignment', 'center', 'FontSize', ppt_fontsize_value, 'FontWeight', 'bold');
    improvement = r2_rf_kf(lbl) - r2_mlr_kf(lbl);
    x_mid = (x_mlr + x_rf) / 2;
    y_top = max(r2_mlr_kf(lbl), r2_rf_kf(lbl)) + 0.06;
    text(x_mid, y_top, sprintf('+%.1f%%p', improvement*100), ...
        'HorizontalAlignment', 'center', 'FontSize', ppt_fontsize_tick, ...
        'FontWeight', 'bold', 'Color', [0.1 0.6 0.1]);
end
yline(0, '-k', 'LineWidth', 0.5);
set(gca, 'XTickLabel', label_names, 'FontSize', ppt_fontsize_tick);
xlabel('Label', 'FontSize', ppt_fontsize_axis);
ylabel('R²', 'FontSize', ppt_fontsize_axis);
title('5-Fold CV: MLR vs RF (R²)', 'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
legend({'MLR', 'Random Forest'}, 'FontSize', ppt_fontsize_tick, 'Location', 'northwest');
grid on; hold off;
saveas(fig_r2comp, fullfile(saveDir, 'MLR_vs_RF_R2_Comparison.fig'));

% ---------- 6.9: SOH Trace Plot (셀별, K-Fold) ----------
fprintf('  6.9: SOH Trace Plot (셀별)...\n');

% K-Fold 결과에서 셀별 예측 추출
fig_trace = figure('Position', [50, 50, 1600, 800], 'Name', 'SOH Trace per Cell');
all_true_cat = vertcat(kfold_Y_true{:});
all_pred_mlr_cat = vertcat(kfold_Y_pred_MLR{:});
all_pred_rf_cat  = vertcat(kfold_Y_pred_RF{:});

% cv_indices 순서대로 원래 cellIDs 유지
cellIDs_all = cellIDs;  % NaN 제거 후의 cellIDs

for ci = 1:n_cells
    subplot(2, 4, ci);
    cell_mask = strcmp(cellIDs_all, unique_cells{ci});
    y_true_cell = Y(cell_mask, 1);  % SOH
    
    % K-Fold에서 이 셀의 예측값 찾기 (test에 해당하는 샘플)
    % 원본 인덱스 기반으로 예측값 모으기
    pred_rf_cell = zeros(sum(cell_mask), 1);
    pred_mlr_cell = zeros(sum(cell_mask), 1);
    cell_indices = find(cell_mask);
    
    for s = 1:length(cell_indices)
        orig_idx = cell_indices(s);
        fold_of_sample = cv_indices(orig_idx);
        % 이 fold의 test set에서 이 샘플의 위치 찾기
        test_indices_in_fold = find(cv_indices == fold_of_sample);
        pos_in_fold = find(test_indices_in_fold == orig_idx);
        pred_rf_cell(s) = kfold_Y_pred_RF{fold_of_sample}(pos_in_fold, 1);
        pred_mlr_cell(s) = kfold_Y_pred_MLR{fold_of_sample}(pos_in_fold, 1);
    end
    
    x_t = 1:sum(cell_mask);
    plot(x_t, y_true_cell, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
        'MarkerFaceColor', 'k', 'DisplayName', 'Actual');
    hold on;
    plot(x_t, pred_mlr_cell, 'b--s', 'LineWidth', 1, 'MarkerSize', 5, ...
        'MarkerFaceColor', [0.2 0.5 0.85], 'DisplayName', 'MLR');
    plot(x_t, pred_rf_cell, 'r--^', 'LineWidth', 1, 'MarkerSize', 5, ...
        'MarkerFaceColor', [0.85 0.33 0.1], 'DisplayName', 'RF');
    
    % R² 계산
    ss_res_t = sum((y_true_cell - pred_rf_cell).^2);
    ss_tot_t = sum((y_true_cell - mean(y_true_cell)).^2);
    if ss_tot_t == 0, r2_t = NaN; else, r2_t = 1 - ss_res_t / ss_tot_t; end
    
    title(sprintf('%s (R²=%.2f)', unique_cells{ci}, r2_t), ...
        'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Sample', 'FontSize', 10);
    ylabel('SOH (Ah)', 'FontSize', 10);
    grid on;
    if ci == 1, legend('Location', 'best', 'FontSize', 7); end
    hold off;
end

sgtitle('5-Fold CV: SOH Trace per Cell (Actual vs MLR vs RF)', ...
    'FontSize', ppt_fontsize_title, 'FontWeight', 'bold');
saveas(fig_trace, fullfile(saveDir, 'SOH_Trace_Per_Cell.fig'));

fprintf('\n=== Section 6 완료 ===\n');

% --- 최종 파일 목록 ---
fprintf('\n=== 전체 저장 파일 목록 ===\n');
fprintf('  저장 경로: %s\n', saveDir);
fprintf('  [분석]\n');
fprintf('    PCA_Scree_Plot.fig\n');
fprintf('    MLR_Sanity_Check.fig\n');
fprintf('    RF_Feature_Importance.fig\n');
fprintf('    RF_Error_Distribution.fig\n');
fprintf('    RF_Residual_Plot.fig\n');
fprintf('  [PPT 요약]\n');
fprintf('    KFold_RF_Predicted_vs_Actual.fig\n');
fprintf('    MLR_vs_RF_SOH_Comparison.fig\n');
fprintf('    Overall_Label_Comparison.fig\n');
fprintf('    MLR_vs_RF_R2_Comparison.fig\n');
fprintf('  [프로세스 증명]\n');
fprintf('    HP_Tuning_Heatmap.fig\n');
fprintf('    SOH_Trace_Per_Cell.fig\n');
fprintf('\n=== 전체 모델링 파이프라인 종료. ===\n');

%% [ADD] Save Modeling Results for Visualization Script
resultMatPath = fullfile(saveDir, 'Modeling_Results.mat');
fprintf('\nSaving Modeling Results to %s ...\n', resultMatPath);

save(resultMatPath, ...
    'metrics_MLR', 'metrics_RF', ...            % Performance Metrics
    'kfold_Y_true', ...                         % Ground Truth (Fold-wise)
    'kfold_Y_pred_MLR', 'kfold_Y_pred_RF', ...  % Predictions (Fold-wise)
    'label_names', 'feature_names', ...         % Metadata
    'X', 'Y', 'cellIDs', 'cv_indices', ...      % Raw Data used for modeling
    'n_features', 'n_labels', 'n_samples', ...  % Dimensions
    'rmse_heatmap', 'hp_trees_grid', 'hp_leaf_grid'); % [ADD] Grid Search Results

fprintf('Save Complete.\n');
