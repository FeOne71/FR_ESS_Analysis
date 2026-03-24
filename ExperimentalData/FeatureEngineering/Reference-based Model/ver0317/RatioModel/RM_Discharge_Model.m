% RM_Discharge_Model.m
% 방전 전용 마스킹 증강 모델 (XGBoost, LightGBM)
% 피처: dQ_d_ratio_01~11 (11개) + C_eff_dch (1개) = 12 피처
% NaN 네이티브 처리 (0 임퓨테이션 금지)

clear; clc; close all;

%% 1. Paths & Setup
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
rmDir   = fullfile(verDir, 'RatioModel');
visDir  = fullfile(rmDir, 'Visualization');
if ~exist(visDir,'dir'), mkdir(visDir); end

pyenv('Version', 'C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

%% 2. Load Lab Data
d = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat'));
FM = d.FM;
Q_nom = 64;

dQ_d_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false);
dQ_d_raw = FM{:, dQ_d_cols};   % 215 x 11
C_eff_d  = FM.C_eff_dch;
y_lab    = FM.Static_Capacity / Q_nom * 100;
cellIDs  = FM.CellID;

nRows = size(dQ_d_raw, 1);
nD = size(dQ_d_raw, 2);  % 11

%% 3. Generate contiguous window patterns
dch_patterns = {};
for w = 1:nD
    for s = 1:(nD-w+1)
        dch_patterns{end+1} = s:(s+w-1);
    end
end
fprintf('Contiguous discharge patterns: %d\n', length(dch_patterns));

%% 4. Masking Augmented Data
rng(42);
N_aug = 30;

X_aug = []; y_aug = []; cellID_aug = string([]);

for i = 1:nRows
    dQ = dQ_d_raw(i,:);
    
    % Original (no masking)
    ratio_orig = dQ / sum(dQ, 'omitnan');
    X_aug(end+1,:) = [ratio_orig, C_eff_d(i)];
    y_aug(end+1,1) = y_lab(i);
    cellID_aug(end+1,1) = string(cellIDs(i));
    
    % N_aug masked versions
    for a = 1:N_aug
        dp = dch_patterns{randi(length(dch_patterns))};
        dQ_masked = nan(1, nD);
        dQ_masked(dp) = dQ(dp);
        ratio_m = dQ_masked / sum(dQ_masked, 'omitnan');
        X_aug(end+1,:) = [ratio_m, C_eff_d(i)];
        y_aug(end+1,1) = y_lab(i);
        cellID_aug(end+1,1) = string(cellIDs(i));
    end
end

fprintf('Augmented: %d rows (%.1f%% NaN)\n', size(X_aug,1), sum(isnan(X_aug(:)))/numel(X_aug)*100);

%% 5. Standardize (NaN stays NaN)
mu_dch = mean(X_aug, 1, 'omitnan');
sig_dch = std(X_aug, 0, 1, 'omitnan');
sig_dch(sig_dch==0) = 1;
X_aug_s = (X_aug - mu_dch) ./ sig_dch;

%% 6. LOCO-CV
uniqueCells = unique(cellIDs);
nCells = length(uniqueCells);
modelTypes = {'XGBoost', 'LightGBM'};

fprintf('\n=== DISCHARGE MODEL LOCO-CV ===\n');

for m = 1:2
    mName = modelTypes{m};
    y_true_all = []; y_pred_all = [];
    
    for fold = 1:nCells
        testCell = uniqueCells(fold);
        
        trainMask = ~strcmp(cellID_aug, testCell);
        X_tr = X_aug_s(trainMask,:); y_tr = y_aug(trainMask);
        
        testMask_orig = strcmp(cellIDs, testCell);
        X_te_orig = dQ_d_raw(testMask_orig,:);
        C_te = C_eff_d(testMask_orig);
        y_te = y_lab(testMask_orig);
        
        ratio_te = X_te_orig ./ sum(X_te_orig, 2, 'omitnan');
        X_te = [ratio_te, C_te];
        X_te_s = (X_te - mu_dch) ./ sig_dch;
        
        switch mName
            case 'XGBoost'
                mdl = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
                    'learning_rate',0.05,'min_child_weight',int32(3),'subsample',0.8,...
                    'colsample_bytree',0.8,'random_state',int32(42)));
                mdl.fit(py.numpy.array(X_tr), py.numpy.array(y_tr));
                y_p = double(mdl.predict(py.numpy.array(X_te_s)));
            case 'LightGBM'
                mdl = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
                    'learning_rate',0.05,'min_child_samples',int32(5),'subsample',0.8,...
                    'colsample_bytree',0.8,'random_state',int32(42),'verbose',int32(-1)));
                mdl.fit(py.numpy.array(X_tr), py.numpy.array(y_tr));
                y_p = double(mdl.predict(py.numpy.array(X_te_s)));
        end
        y_true_all = [y_true_all; y_te]; y_pred_all = [y_pred_all; y_p(:)];
    end
    
    err = y_true_all - y_pred_all;
    rmse = sqrt(mean(err.^2));
    mape = mean(abs(err./y_true_all))*100;
    r2 = 1 - sum(err.^2)/sum((y_true_all-mean(y_true_all)).^2);
    fprintf('  %s: RMSE=%.3f%%, MAPE=%.3f%%, R2=%.4f\n', mName, rmse, mape, r2);
    CV_Results.(mName) = struct('RMSE', rmse, 'MAPE', mape, 'R2', r2);
end

%% 7. Train final models
fprintf('\n=== Training final DISCHARGE models ===\n');

FinalModels.XGB = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
    'learning_rate',0.05,'min_child_weight',int32(3),'subsample',0.8,...
    'colsample_bytree',0.8,'random_state',int32(42)));
FinalModels.XGB.fit(py.numpy.array(X_aug_s), py.numpy.array(y_aug));
fprintf('  XGBoost trained.\n');

FinalModels.LGB = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
    'learning_rate',0.05,'min_child_samples',int32(5),'subsample',0.8,...
    'colsample_bytree',0.8,'random_state',int32(42),'verbose',int32(-1)));
FinalModels.LGB.fit(py.numpy.array(X_aug_s), py.numpy.array(y_aug));
fprintf('  LightGBM trained.\n');

%% 8. Save
save(fullfile(rmDir, 'RM_Discharge_Final.mat'), 'FinalModels', 'CV_Results', 'mu_dch', 'sig_dch');
fprintf('Saved: RM_Discharge_Final.mat\n');
fprintf('=== Discharge Model Complete! ===\n');
