% DM_01_Train_Models.m
% Discharge-Only Model: LOCO-CV Training
% Features: dQ_d_01~11 + C_eff_dch = 12 features (dQ_d_12 excluded due to NaN)
% Target: SOH% = Static_Capacity / 64 * 100
% Models: LASSO, RF, GBM, XGBoost, LightGBM

clear; clc; close all;

%% 1. Load Data
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
dmDir   = fullfile(verDir, 'DischargeModel');
visDir  = fullfile(dmDir, 'Visualization');
if ~exist(visDir,'dir'), mkdir(visDir); end

pyenv('Version', 'C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

d = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat'));
FM = d.FM;

%% 2. Define Features & Target
Q_nom = 64;
feat_names = [arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false), {'C_eff_dch'}];
nFeat = length(feat_names);
targetCol = 'Static_Capacity';

X_lab = FM{:, feat_names};
y_lab = FM.(targetCol) / Q_nom * 100;  % SOH%

% Standardize
mu_lab = mean(X_lab, 1);
sigma_lab = std(X_lab, 0, 1);
sigma_lab(sigma_lab == 0) = 1;
X_lab_s = (X_lab - mu_lab) ./ sigma_lab;

%% 3. LOCO-CV
cellIDs = unique(FM.CellID);
nCells = length(cellIDs);
modelNames = {'LASSO','RF','GBM','XGBoost','LightGBM'};
nModels = length(modelNames);

Results = table();
resultIdx = 0;

fprintf('=== Discharge-Only Model LOCO-CV (%d features) ===\n', nFeat);

for m = 1:nModels
    mName = modelNames{m};
    fprintf('  Model: %s ... ', mName);
    
    y_true_all = [];
    y_pred_all = [];
    
    for fold = 1:nCells
        testIdx = FM.CellID == cellIDs(fold);
        trainIdx = ~testIdx;
        
        X_tr = X_lab_s(trainIdx,:); y_tr = y_lab(trainIdx);
        X_te = X_lab_s(testIdx,:);  y_te = y_lab(testIdx);
        
        switch mName
            case 'LASSO'
                [B, FI] = lasso(X_tr, y_tr, 'CV', 4);
                idx = FI.IndexMinMSE;
                y_p = X_te * B(:,idx) + FI.Intercept(idx);
            case 'RF'
                mdl = fitrensemble(X_tr, y_tr, 'Method', 'Bag', ...
                    'NumLearningCycles', 100, ...
                    'Learners', templateTree('MaxNumSplits', 20, 'MinLeafSize', 5, 'Surrogate', 'on'));
                y_p = predict(mdl, X_te);
            case 'GBM'
                mdl = fitrensemble(X_tr, y_tr, 'Method', 'LSBoost', ...
                    'NumLearningCycles', 100, 'LearnRate', 0.1, ...
                    'Learners', templateTree('MaxNumSplits', 10, 'MinLeafSize', 5, 'Surrogate', 'on'));
                y_p = predict(mdl, X_te);
            case 'XGBoost'
                xgb = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(100),...
                    'max_depth',int32(5),'learning_rate',0.1,'min_child_weight',int32(3),...
                    'subsample',0.8,'colsample_bytree',0.8,'random_state',int32(42)));
                xgb.fit(py.numpy.array(X_tr), py.numpy.array(y_tr));
                y_p = double(xgb.predict(py.numpy.array(X_te)));
            case 'LightGBM'
                lgb = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(100),...
                    'max_depth',int32(5),'learning_rate',0.1,'min_child_samples',int32(5),...
                    'subsample',0.8,'colsample_bytree',0.8,'random_state',int32(42),'verbose',int32(-1)));
                lgb.fit(py.numpy.array(X_tr), py.numpy.array(y_tr));
                y_p = double(lgb.predict(py.numpy.array(X_te)));
        end
        
        y_true_all = [y_true_all; y_te];
        y_pred_all = [y_pred_all; y_p(:)];
    end
    
    err = y_true_all - y_pred_all;
    rmse = sqrt(mean(err.^2));
    mape = mean(abs(err ./ y_true_all)) * 100;
    r2 = 1 - sum(err.^2) / sum((y_true_all - mean(y_true_all)).^2);
    
    fprintf('RMSE=%.3f%%, MAPE=%.2f%%, R2=%.4f\n', rmse, mape, r2);
    resultIdx = resultIdx + 1;
    Results(resultIdx,:) = table({mName}, nFeat, rmse, mape, r2, ...
        'VariableNames', {'Model','NumFeatures','RMSE_pct','MAPE','R2'});
end

%% 4. Display & Save
fprintf('\n========== DISCHARGE MODEL LOCO-CV RESULTS ==========\n');
disp(Results);
writetable(Results, fullfile(dmDir, 'DM_LOCO_CV_Results.csv'));
save(fullfile(dmDir, 'DM_LOCO_CV_Results.mat'), 'Results', 'mu_lab', 'sigma_lab', 'feat_names');

%% 5. Bar Chart
fig = figure('Position', [50 100 800 400]);
b = bar(Results.RMSE_pct);
set(gca, 'XTickLabel', Results.Model, 'FontSize', 12);
ylabel('RMSE (SOH %)', 'FontSize', 13, 'FontWeight', 'bold');
title('Discharge-Only Model: LOCO-CV RMSE Comparison', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
for i = 1:height(Results)
    text(i, Results.RMSE_pct(i)+0.05, sprintf('%.2f%%', Results.RMSE_pct(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end
saveas(fig, fullfile(visDir, 'DM_RMSE_Comparison.png'));
fprintf('Saved: DM_RMSE_Comparison.png\n');
fprintf('=== Discharge Model Training Complete! ===\n');
