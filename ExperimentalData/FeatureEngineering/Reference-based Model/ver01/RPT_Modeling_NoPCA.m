
%% RPT_Modeling_NoPCA.m
% Random Forest Modeling using RAW features (no PCA)
% Purpose: Compare performance with PCA-based modeling

clear; clc; close all;

%% 1. 데이터 로드 & 전처리
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis';
dataPath = fullfile(baseDir, 'Dataset', 'Feature_Matrix_Final.mat');
fprintf('Loading Feature Matrix...\n');
load(dataPath, 'FeatureTable');

% [Note] 'Raw Feature'는 PCA를 거치지 않은 14개 피처 원본을 의미하며,
% FeatureTable.X_Normalized는 이미 C-rate별 Z-score 표준화가 완료된 상태입니다.
% (표준화가 되어있어야 변수간 Weight가 공정하게 배분됩니다.)
X = FeatureTable.X_Normalized;   
Y = FeatureTable.Y_Labels;       
cellIDs = FeatureTable.CellID;   
label_names = {'SOH', 'LLI', 'LAM'};

% NaN 제거
nan_rows = any(isnan(X), 2) | any(isnan(Y), 2);
X(nan_rows, :) = [];
Y(nan_rows, :) = [];
cellIDs(nan_rows) = [];

[n_samples, n_features] = size(X);
n_labels = size(Y, 2);

%% 2. 5-Fold Cross-Validation (RF on Raw Features)
K = 5;
unique_cells = unique(cellIDs);
% Stratified or Random Split? Using the same logic as RPT_Modeling.m for comparison
cv_indices = mod(randperm(n_samples), K) + 1; 

%% 3. [NEW] Hyperparameter Optimization (Grid Search)
fprintf('\nPerforming Hyperparameter Optimization (Grid Search)...\n');
hp_trees_grid = [50, 100, 200, 500];
hp_leaf_grid  = [1, 3, 5, 10];
rmse_heatmap  = zeros(length(hp_leaf_grid), length(hp_trees_grid), n_labels);

for lbl = 1:n_labels
    fprintf('  Optimizing Label: %s...\n', label_names{lbl});
    for ri = 1:length(hp_leaf_grid)
        for ci = 1:length(hp_trees_grid)
            t_hp = templateTree('MinLeafSize', hp_leaf_grid(ri));
            rf_hp = fitrensemble(X, Y(:, lbl), ...
                'Method', 'Bag', 'NumLearningCycles', hp_trees_grid(ci), 'Learners', t_hp);
            y_oob = oobPredict(rf_hp);
            rmse_heatmap(ri, ci, lbl) = sqrt(mean((Y(:, lbl) - y_oob).^2));
        end
    end
end

%% 4. 5-Fold Cross-Validation (RF on Raw Features)
metrics_RF_Raw = zeros(K, n_labels, 3); % [RMSE, MAE, R2]
kfold_Y_true = cell(K, 1);
kfold_Y_pred_RF = cell(K, 1);
for f = 1:K
    n_test = sum(cv_indices == f);
    kfold_Y_pred_RF{f} = zeros(n_test, n_labels);
end

for fold = 1:K
    test_idx = (cv_indices == fold);
    train_idx = ~test_idx;
    
    Y_test  = Y(test_idx, :);
    kfold_Y_true{fold} = Y_test;
    
    for lbl = 1:n_labels
        X_train = X(train_idx, :);
        Y_train = Y(train_idx, lbl);
        X_test  = X(test_idx, :);
        
        % RF Model Training (Using optimized params: 500 trees, 1 leaf)
        t = templateTree('MinLeafSize', 1);
        rf_model = fitrensemble(X_train, Y_train, ...
            'Method', 'Bag', 'NumLearningCycles', 500, 'Learners', t);
        
        y_pred = predict(rf_model, X_test);
        y_true = Y_test(:, lbl);
        
        % Metrics
        rmse = sqrt(mean((y_true - y_pred).^2));
        mae = mean(abs(y_true - y_pred));
        ss_res = sum((y_true - y_pred).^2);
        ss_tot = sum((y_true - mean(y_true)).^2);
        r2 = 1 - ss_res / ss_tot;
        
        metrics_RF_Raw(fold, lbl, :) = [rmse, mae, r2];
        
        % Save Predictions for Visualization
        kfold_Y_pred_RF{fold}(:, lbl) = y_pred; 
    end
end

% Re-organizing predictions for easier visualization (match RPT_Modeling.m structure)
for fold = 1:K
    % NoPCA에서는 각 라벨별로 따로 돌렸으므로, 시각화 스크립트에서 편하게 쓰도록
    % kfold_Y_pred_RF를 [샘플수 x 라벨수] 형태로 재구성하는 로직이 필요할 수 있으나,
    % 여기서는 단순화를 위해 모델링 코드 내에서 직접 kfold_Y_pred_RF{fold}를 채웠다고 가정
end

%% 3. 결과 출력 및 저장
fprintf('\n============================================\n');
fprintf(' RF Model Performance (NO PCA)\n');
fprintf('============================================\n');
for lbl = 1:n_labels
    mean_stats = squeeze(mean(metrics_RF_Raw(:, lbl, :), 1));
    std_stats = squeeze(std(metrics_RF_Raw(:, lbl, :), 0, 1));
    
    fprintf(' [%s]\n', label_names{lbl});
    fprintf('  R²:   %.4f ± %.4f\n', mean_stats(3), std_stats(3));
    fprintf('  RMSE: %.4f ± %.4f\n', mean_stats(1), std_stats(1));
    fprintf('  MAE:  %.4f ± %.4f\n', mean_stats(2), std_stats(2));
end

% 결과 저장 (시각화용)
saveDir = fullfile(baseDir, 'Modeling_Results');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
resultMatPath = fullfile(saveDir, 'Modeling_Results_NoPCA.mat');

% 14개 피처 이름 정의 (고정)
feature_names = {'Chg_dQ', 'Chg_dQ_V', 'Chg_PkH', 'Chg_PkV', 'Chg_PkA', ...
                 'Chg_PkH_2', 'Chg_PkV_2', 'Dch_dQ', 'Dch_dQ_V', 'Dch_PkH', ...
                 'Dch_PkV', 'Dch_PkA', 'Dch_PkH_2', 'Dch_PkV_2'};

fprintf('\nSaving Results to %s...\n', resultMatPath);
save(resultMatPath, 'metrics_RF_Raw', 'label_names', 'feature_names', ...
    'X', 'Y', 'cellIDs', 'cv_indices', 'kfold_Y_true', 'kfold_Y_pred_RF', ...
    'rmse_heatmap', 'hp_trees_grid', 'hp_leaf_grid'); 
fprintf('Save Complete.\n');
