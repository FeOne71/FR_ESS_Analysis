% ML_01_Train_Models.m
% LOCO-CV (Leave-One-Cell-Out) ML Training Pipeline
% 5 Models × 4 Feature Sets = 20 experiments
% Models: LASSO, RF, GBM, XGBoost, LightGBM
% Feature Sets: Full(23), Field(5), Top2(3), Top2Field(3)

clear; clc; close all;

%% 1. Load Data
baseDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData';
verDir = fullfile(baseDir, 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
d = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat'));
FM = d.FM;

% Setup Python
pyenv('Version', 'C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

%% 2. Define Feature Sets
% A: Full (23 features) - exclude dQ_c_01, dQ_c_02, dQ_d_12
featA = [arrayfun(@(i) sprintf('dQ_c_%02d',i), 3:12, 'UniformOutput', false), ...
         arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false), ...
         {'C_eff_chg', 'C_eff_dch'}];

% B: Field-Available (Seg07~10 charge + C_eff)
featB = {'dQ_c_07','dQ_c_08','dQ_c_09','dQ_c_10','C_eff_chg','C_eff_dch'};

% C: Top-2 (PCC strongest 2 segments + C_eff)
%    Charge Seg07 (r=0.98 at 0.1C), Discharge Seg09 (r=0.93 at 3C, 0.91 across)
featC = {'dQ_c_07','dQ_d_07','C_eff_chg','C_eff_dch'};

% D: Top-2 Field (field-available PCC strongest 2 + C_eff)
%    Charge Seg08 (r=0.95 at 2C), Charge Seg09 (r=0.85 at 1C)
featD = {'dQ_c_08','dQ_c_09','C_eff_chg','C_eff_dch'};

featureSetNames = {'A_Full', 'B_Field', 'C_Top2', 'D_Top2Field'};
featureSets = {featA, featB, featC, featD};
nFeatSets = length(featureSets);

%% 3. Define Models
modelNames = {'LASSO', 'RF', 'GBM', 'XGBoost', 'LightGBM'};
nModels = length(modelNames);

%% 4. Cell IDs for LOCO-CV
cellIDs = unique(FM.CellID);
nCells = length(cellIDs);
targetCol = 'Static_Capacity';

%% 5. Run LOCO-CV for all combinations
Results = table();
resultIdx = 0;

for fs = 1:nFeatSets
    featNames = featureSets{fs};
    nFeat = length(featNames);
    fsName = featureSetNames{fs};
    fprintf('\n=== Feature Set: %s (%d features) ===\n', fsName, nFeat);
    
    for m = 1:nModels
        mName = modelNames{m};
        fprintf('  Model: %s ... ', mName);
        
        y_true_all = [];
        y_pred_all = [];
        
        for fold = 1:nCells
            % Split
            testCell = cellIDs(fold);
            testIdx = FM.CellID == testCell;
            trainIdx = ~testIdx;
            
            X_train = FM{trainIdx, featNames};
            y_train = FM{trainIdx, targetCol};
            X_test = FM{testIdx, featNames};
            y_test = FM{testIdx, targetCol};
            
            % Standardize (train stats only) - applied to ALL models per paper
            mu = mean(X_train, 1);
            sigma = std(X_train, 0, 1);
            sigma(sigma == 0) = 1;
            X_train_s = (X_train - mu) ./ sigma;
            X_test_s = (X_test - mu) ./ sigma;
            
            % Train & Predict
            switch mName
                case 'LASSO'
                    [B, FitInfo] = lasso(X_train_s, y_train, 'CV', 4);
                    idxLambda = FitInfo.IndexMinMSE;
                    coef = B(:, idxLambda);
                    intercept = FitInfo.Intercept(idxLambda);
                    y_pred = X_test_s * coef + intercept;
                    
                case 'RF'
                    mdl = fitrensemble(X_train_s, y_train, ...
                        'Method', 'Bag', ...
                        'NumLearningCycles', 100, ...
                        'Learners', templateTree('MaxNumSplits', 20, 'MinLeafSize', 5));
                    y_pred = predict(mdl, X_test_s);
                    
                case 'GBM'
                    mdl = fitrensemble(X_train_s, y_train, ...
                        'Method', 'LSBoost', ...
                        'NumLearningCycles', 100, ...
                        'LearnRate', 0.1, ...
                        'Learners', templateTree('MaxNumSplits', 10, 'MinLeafSize', 5));
                    y_pred = predict(mdl, X_test_s);
                    
                case 'XGBoost'
                    xgb = py.xgboost.XGBRegressor(...
                        pyargs('n_estimators', int32(100), ...
                               'max_depth', int32(5), ...
                               'learning_rate', 0.1, ...
                               'min_child_weight', int32(3), ...
                               'subsample', 0.8, ...
                               'colsample_bytree', 0.8, ...
                               'random_state', int32(42)));
                    X_tr_py = py.numpy.array(X_train_s);
                    y_tr_py = py.numpy.array(y_train);
                    xgb.fit(X_tr_py, y_tr_py);
                    X_te_py = py.numpy.array(X_test_s);
                    y_pred_py = xgb.predict(X_te_py);
                    y_pred = double(y_pred_py);
                    
                case 'LightGBM'
                    lgb = py.lightgbm.LGBMRegressor(...
                        pyargs('n_estimators', int32(100), ...
                               'max_depth', int32(5), ...
                               'learning_rate', 0.1, ...
                               'min_child_samples', int32(5), ...
                               'subsample', 0.8, ...
                               'colsample_bytree', 0.8, ...
                               'random_state', int32(42), ...
                               'verbose', int32(-1)));
                    X_tr_py = py.numpy.array(X_train_s);
                    y_tr_py = py.numpy.array(y_train);
                    lgb.fit(X_tr_py, y_tr_py);
                    X_te_py = py.numpy.array(X_test_s);
                    y_pred_py = lgb.predict(X_te_py);
                    y_pred = double(y_pred_py);
            end
            
            y_true_all = [y_true_all; y_test];
            y_pred_all = [y_pred_all; y_pred(:)];
        end
        
        % Compute metrics
        err = y_true_all - y_pred_all;
        rmse = sqrt(mean(err.^2));
        mape = mean(abs(err ./ y_true_all)) * 100;
        ss_res = sum(err.^2);
        ss_tot = sum((y_true_all - mean(y_true_all)).^2);
        r2 = 1 - ss_res / ss_tot;
        
        fprintf('RMSE=%.3f Ah, MAPE=%.2f%%, R2=%.4f\n', rmse, mape, r2);
        
        resultIdx = resultIdx + 1;
        Results(resultIdx,:) = table({fsName}, {mName}, nFeat, rmse, mape, r2, ...
            'VariableNames', {'FeatureSet','Model','NumFeatures','RMSE','MAPE','R2'});
    end
end

%% 6. Display and Save Results
fprintf('\n\n========== FINAL RESULTS ==========\n');
disp(Results);

writetable(Results, fullfile(verDir, 'ML_01_Results.csv'));
save(fullfile(verDir, 'ML_01_Results.mat'), 'Results');
fprintf('Saved: ML_01_Results.csv, ML_01_Results.mat\n');

%% 7. Quick Visualization - RMSE comparison bar chart
fig = figure('Position', [50 100 1400 500]);

% Reshape for grouped bar
rmse_matrix = NaN(nModels, nFeatSets);
for fs = 1:nFeatSets
    for m = 1:nModels
        idx = strcmp(Results.FeatureSet, featureSetNames{fs}) & ...
              strcmp(Results.Model, modelNames{m});
        rmse_matrix(m, fs) = Results.RMSE(idx);
    end
end

b = bar(rmse_matrix);
set(gca, 'XTickLabel', modelNames, 'FontSize', 11);
ylabel('RMSE (Ah)', 'FontSize', 12, 'FontWeight', 'bold');
title('SOH Estimation RMSE: Model × Feature Set Comparison (LOCO-CV)', ...
    'FontSize', 14, 'FontWeight', 'bold');
legend(strrep(featureSetNames, '_', ' '), 'Location', 'northwest', 'FontSize', 10);
grid on;

% Annotate values
for fs = 1:nFeatSets
    for m = 1:nModels
        val = rmse_matrix(m, fs);
        text(b(fs).XEndPoints(m), b(fs).YEndPoints(m) + 0.02, ...
            sprintf('%.3f', val), 'HorizontalAlignment', 'center', ...
            'FontSize', 7, 'FontWeight', 'bold');
    end
end

visDir = fullfile(verDir, 'Visualization');
if ~exist(visDir, 'dir'), mkdir(visDir); end
saveas(fig, fullfile(visDir, 'ML_RMSE_Comparison.png'));
fprintf('Saved: ML_RMSE_Comparison.png\n');

fprintf('\n=== ML Training Pipeline Complete! ===\n');
