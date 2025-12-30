%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 05_DriveCycle_XGBoost_Training.m
% XGBoost 모델 학습 (4-Fold Cross Validation)
% 
% 목적: 
% - 각 Fold별로 XGBoost 모델 학습
% - 하이퍼파라미터 튜닝 (Grid Search)
% - Early Stopping으로 과적합 방지
% - 각 Fold별 모델 저장 및 평가
%
% 입력:
% - DriveCycle_Fold*_Train.mat, DriveCycle_Fold*_Test.mat
% - Selected_Features.mat+
%+
% 출력:
% - XGBoost_Fold*_Model.mat (각 Fold별 모델)
% - XGBoost_Training_Results.mat (학습 결과)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

fprintf('=== Drive Cycle XGBoost Training (4-Fold CV) ===\n');

%% Configuration - User Settings
% =========================================================================
inputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
outputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Target variables (Multi-Target Regression: each target trained separately)
% XGBoost supports only single target variable, so we train one model per target
targetVars = {'Capacity_C3', 'Capacity_OCV'};  % List of target variables to train
% Note: Each target variable will be trained with a separate XGBoost model

% Event Type Filtering (충전/방전 분리 모델링)
% [방법 1] 충전 모델 만들 때
targetEventType = 'Discharge'; 

% [방법 2] 방전 모델 만들 때 (현재 방전 이벤트 데이터 없음)
% targetEventType = 'Discharge';  % WARNING: No discharge events in data!

% [방법 3] 통합 모델 (현재는 Charge만 있으므로 Charge와 동일)
% targetEventType = 'All';

% Hyperparameter grid for tuning (Improved for better performance)
hyperparams = struct();
hyperparams.max_depth = [4, 5, 6, 7];  % Increased max depth
hyperparams.learning_rate = [0.05, 0.1, 0.15];  % Multiple learning rates
hyperparams.n_estimators = [200, 300, 500];  % Increased n_estimators for better learning
hyperparams.subsample = [0.8, 1.0];  % Add subsampling to prevent overfitting
hyperparams.colsample_bytree = [0.8, 1.0];  % Add column subsampling
hyperparams.reg_alpha = [0, 0.1];  % L1 regularization
hyperparams.reg_lambda = [1, 1.5];  % L2 regularization

% Early stopping
early_stopping_rounds = 50;
validation_fraction = 0.2;  % 20% of train data for validation

% Number of folds
n_folds = 4;

fprintf('\n=== Model Configuration ===\n');
fprintf('Model Type: XGBoost (via xgboost_train.m using XGBoost DLL)\n');
fprintf('Target variables: %d variables (Multi-Target Regression)\n', length(targetVars));
fprintf('  Target list: %s\n', strjoin(targetVars, ', '));
fprintf('  [INFO] XGBoost supports only single target variable.\n');
fprintf('         Each target will be trained with a separate model.\n');
fprintf('         Total models to train: %d (one per target variable)\n', length(targetVars));
fprintf('Event type filter: %s\n', targetEventType);
if strcmp(targetEventType, 'Charge')
    fprintf('  -> Training CHARGE-only model\n');
elseif strcmp(targetEventType, 'Discharge')
    fprintf('  -> Training DISCHARGE-only model\n');
else
    fprintf('  -> Training INTEGRATED model (Charge + Discharge)\n');
end
fprintf('Cross-Validation: %d-Fold CV\n', n_folds);
fprintf('Validation split: %.1f%% of training data\n', validation_fraction * 100);
fprintf('Early stopping rounds: %d\n', early_stopping_rounds);
fprintf('\nHyperparameter Grid:\n');
fprintf('  max_depth: %s\n', mat2str(hyperparams.max_depth));
fprintf('  learning_rate: %s\n', mat2str(hyperparams.learning_rate));
fprintf('  n_estimators: %s\n', mat2str(hyperparams.n_estimators));
fprintf('  subsample: %s\n', mat2str(hyperparams.subsample));
fprintf('  colsample_bytree: %s\n', mat2str(hyperparams.colsample_bytree));
% =========================================================================

%% Load selected features
fprintf('\n=== Loading Selected Features ===\n');
features_file = fullfile(inputDir, 'Selected_Features.mat');
if ~exist(features_file, 'file')
    error('Selected features not found. Please run DriveCycle_FeatureEngineering_04.m first');
end

load(features_file, 'selected_features');
fprintf('Loaded %d selected features\n', length(selected_features));

%% Initialize results storage for all targets
all_target_results = struct();

%% Process each target variable (Multi-Target Regression)
for targetIdx = 1:length(targetVars)
    targetVar = targetVars{targetIdx};
    
    fprintf('\n');
    fprintf('########################################################################\n');
    fprintf('### Training Model for Target Variable: %s (%d/%d) ###\n', targetVar, targetIdx, length(targetVars));
    fprintf('########################################################################\n');
    
    % Initialize results storage for this target
    fold_results = struct();
    all_predictions = [];
    all_actuals = [];
    all_fold_ids = [];
    
    %% Process each fold
    for fold_num = 1:n_folds
    fprintf('\n========================================\n');
    fprintf('=== Processing Fold %d ===\n', fold_num);
    fprintf('========================================\n');
    
    % Load fold data
    train_file = fullfile(inputDir, sprintf('DriveCycle_Fold%d_Train.mat', fold_num));
    test_file = fullfile(inputDir, sprintf('DriveCycle_Fold%d_Test.mat', fold_num));
    
    if ~exist(train_file, 'file') || ~exist(test_file, 'file')
        fprintf('WARNING: Fold %d data not found, skipping...\n', fold_num);
        continue;
    end
    
    load(train_file, 'trainData');
    load(test_file, 'testData');
    
    fprintf('Original train data: %d rows\n', height(trainData));
    fprintf('Original test data: %d rows\n', height(testData));
    
    % Debug: Check available columns
    fprintf('Available columns in train data: %s\n', strjoin(trainData.Properties.VariableNames, ', '));
    
    % Debug: Check EventType values before filtering
    if ismember('EventType', trainData.Properties.VariableNames)
        unique_event_types = unique(trainData.EventType);
        fprintf('Available EventType values in train data: %s\n', strjoin(string(unique_event_types), ', '));
        
        % Count by EventType
        for i = 1:length(unique_event_types)
            et = unique_event_types(i);
            count = sum(strcmpi(trainData.EventType, et));
            fprintf('  %s: %d rows\n', string(et), count);
        end
    else
        fprintf('WARNING: EventType column not found in train data!\n');
    end
    
    % Debug: Check SOC column
    if ismember('SOC', trainData.Properties.VariableNames)
        unique_soc = unique(trainData.SOC(~isnan(trainData.SOC)));
        fprintf('Available SOC values in train data: %s\n', mat2str(unique_soc));
    else
        fprintf('WARNING: SOC column not found in train data!\n');
        fprintf('Please re-run DriveCycle_DataAggregation_02.m to add SOC column.\n');
    end
    
    % Filter by EventType if specified
    if exist('targetEventType', 'var') && ~strcmp(targetEventType, 'All')
        if ismember('EventType', trainData.Properties.VariableNames)
            % Use case-insensitive matching and handle potential whitespace
            train_mask = strcmpi(strtrim(cellstr(trainData.EventType)), targetEventType);
            test_mask = strcmpi(strtrim(cellstr(testData.EventType)), targetEventType);
            
            trainData = trainData(train_mask, :);
            testData = testData(test_mask, :);
            
            fprintf('Filtered by EventType: %s\n', targetEventType);
            fprintf('Train data: %d rows (after filtering, from %d original)\n', height(trainData), sum(train_mask == false) + height(trainData));
            fprintf('Test data: %d rows (after filtering, from %d original)\n', height(testData), sum(test_mask == false) + height(testData));
            
            if height(trainData) == 0
                fprintf('  ERROR: No training data after filtering!\n');
                fprintf('  Target EventType: "%s"\n', targetEventType);
                if ismember('EventType', trainData.Properties.VariableNames)
                    fprintf('  Available EventTypes: %s\n', strjoin(unique(trainData.EventType), ', '));
                end
                fprintf('  Skipping fold %d\n', fold_num);
                continue;
            end
            if height(testData) == 0
                fprintf('  ERROR: No test data after filtering!\n');
                fprintf('  Target EventType: "%s"\n', targetEventType);
                if ismember('EventType', testData.Properties.VariableNames)
                    fprintf('  Available EventTypes: %s\n', strjoin(unique(testData.EventType), ', '));
                end
                fprintf('  Skipping fold %d\n', fold_num);
                continue;
            end
        else
            fprintf('  WARNING: EventType column not found. Using all data.\n');
        end
    else
        fprintf('Using all event types (integrated model)\n');
    end
    
    fprintf('Final train data: %d rows\n', height(trainData));
    fprintf('Final test data: %d rows\n', height(testData));
    
    % Debug: Data distribution analysis
    fprintf('\n=== Data Distribution Analysis ===\n');
    if ismember(targetVar, trainData.Properties.VariableNames)
        y_train_all = trainData.(targetVar);
        y_test_all = testData.(targetVar);
        fprintf('Target variable: %s\n', targetVar);
        fprintf('Train data - Target statistics:\n');
        fprintf('  Mean: %.4f, Std: %.4f\n', mean(y_train_all(~isnan(y_train_all))), std(y_train_all(~isnan(y_train_all))));
        fprintf('  Min: %.4f, Max: %.4f\n', min(y_train_all(~isnan(y_train_all))), max(y_train_all(~isnan(y_train_all))));
        fprintf('  Valid samples: %d / %d (%.1f%%)\n', sum(~isnan(y_train_all)), length(y_train_all), 100*sum(~isnan(y_train_all))/length(y_train_all));
        fprintf('Test data - Target statistics:\n');
        fprintf('  Mean: %.4f, Std: %.4f\n', mean(y_test_all(~isnan(y_test_all))), std(y_test_all(~isnan(y_test_all))));
        fprintf('  Min: %.4f, Max: %.4f\n', min(y_test_all(~isnan(y_test_all))), max(y_test_all(~isnan(y_test_all))));
        fprintf('  Valid samples: %d / %d (%.1f%%)\n', sum(~isnan(y_test_all)), length(y_test_all), 100*sum(~isnan(y_test_all))/length(y_test_all));
    end
    
    % SOC distribution analysis
    if ismember('SOC', trainData.Properties.VariableNames)
        fprintf('\nSOC distribution in train data:\n');
        soc_train_all = trainData.SOC;
        unique_soc = unique(soc_train_all(~isnan(soc_train_all)));
        for s = 1:length(unique_soc)
            soc_val = unique_soc(s);
            count = sum(soc_train_all == soc_val);
            fprintf('  SOC %d: %d samples (%.1f%%)\n', soc_val, count, 100*count/length(soc_train_all));
        end
    end
    
    % Prepare features and target
    % Check which features are available in the data
    available_features = intersect(selected_features, trainData.Properties.VariableNames);
    fprintf('Available features: %d / %d\n', length(available_features), length(selected_features));
    
    % [핵심 방어] 혹시라도 04번에서 안 걸러진 RPT 변수가 있다면 최후의 방어 (삭제)
    rpt_in_features = available_features(startsWith(available_features, 'RPT_'));
    if ~isempty(rpt_in_features)
        fprintf('  WARNING: Found %d RPT variables in selected features (removing to prevent data leakage)\n', length(rpt_in_features));
        available_features = available_features(~startsWith(available_features, 'RPT_'));
        fprintf('  Removed RPT variables: %s\n', strjoin(rpt_in_features, ', '));
    end
    
    % [핵심] 통합 모델을 위해 'SOC' 변수 강제 주입
    % (이유: SOC 50과 90의 저항 차이를 모델이 알아야 하므로)
    if ismember('SOC', trainData.Properties.VariableNames) && ~ismember('SOC', available_features)
        available_features{end+1} = 'SOC';
        fprintf('  Added "SOC" to feature list for integrated modeling.\n');
    end
    
    fprintf('Final feature count: %d\n', length(available_features));
    
    % Extract X and y
    X_train = trainData{:, available_features};
    y_train = trainData.(targetVar);
    
    X_test = testData{:, available_features};
    y_test = testData.(targetVar);
    
    % Store SOC data for later analysis (if available)
    if ismember('SOC', trainData.Properties.VariableNames)
        soc_train = trainData.SOC;
        soc_test = testData.SOC;
    else
        soc_train = [];
        soc_test = [];
        fprintf('  WARNING: SOC column not found in data!\n');
        fprintf('  SOC information is critical for unified model performance.\n');
    end
    
    % Remove rows with NaN in target
    train_valid = ~isnan(y_train);
    test_valid = ~isnan(y_test);
    
    X_train = X_train(train_valid, :);
    y_train = y_train(train_valid);
    if exist('soc_train', 'var') && ~isempty(soc_train)
        soc_train = soc_train(train_valid);
    end
    
    X_test = X_test(test_valid, :);
    y_test = y_test(test_valid);
    if exist('soc_test', 'var') && ~isempty(soc_test)
        soc_test = soc_test(test_valid);
    end
    if exist('testData', 'var')
        testData = testData(test_valid, :);
    end
    
    fprintf('Valid train samples: %d\n', size(X_train, 1));
    fprintf('Valid test samples: %d\n', size(X_test, 1));
    
    % Debug: Feature statistics
    fprintf('\n=== Feature Statistics ===\n');
    fprintf('Feature count: %d\n', size(X_train, 2));
    fprintf('Feature names: %s\n', strjoin(available_features, ', '));
    fprintf('\nFeature missing value analysis:\n');
    for f = 1:size(X_train, 2)
        missing_train = sum(isnan(X_train(:, f)));
        missing_test = sum(isnan(X_test(:, f)));
        if missing_train > 0 || missing_test > 0
            fprintf('  %s: Train missing=%d (%.1f%%), Test missing=%d (%.1f%%)\n', ...
                available_features{f}, missing_train, 100*missing_train/size(X_train,1), ...
                missing_test, 100*missing_test/size(X_test,1));
        end
    end
    
    % Remove rows with NaN in features (no imputation)
    % Check for NaN in features
    train_nan_mask = any(isnan(X_train), 2);
    test_nan_mask = any(isnan(X_test), 2);
    
    if sum(train_nan_mask) > 0
        fprintf('  WARNING: Removing %d training samples with NaN in features (%.1f%%)\n', ...
            sum(train_nan_mask), 100*sum(train_nan_mask)/size(X_train,1));
        X_train = X_train(~train_nan_mask, :);
        y_train = y_train(~train_nan_mask);
        if exist('soc_train', 'var') && ~isempty(soc_train)
            soc_train = soc_train(~train_nan_mask);
        end
    end
    
    if sum(test_nan_mask) > 0
        fprintf('  WARNING: Removing %d test samples with NaN in features (%.1f%%)\n', ...
            sum(test_nan_mask), 100*sum(test_nan_mask)/size(X_test,1));
        X_test = X_test(~test_nan_mask, :);
        y_test = y_test(~test_nan_mask);
        if exist('soc_test', 'var') && ~isempty(soc_test) && length(soc_test) == length(test_nan_mask)
            soc_test = soc_test(~test_nan_mask);
        end
        if exist('testData', 'var') && height(testData) == length(test_nan_mask)
            testData = testData(~test_nan_mask, :);
        end
    end
    
    fprintf('After NaN removal - Train: %d samples, Test: %d samples\n', ...
        size(X_train, 1), size(X_test, 1));
    
    % Z-score normalization (Standardization)
    % Calculate statistics from training data only
    fprintf('\n=== Feature Normalization (Z-score) ===\n');
    feature_means = mean(X_train, 1);
    feature_stds = std(X_train, 0, 1);  % 0 = sample standard deviation (N-1)
    
    % Handle zero standard deviation (constant features)
    zero_std_mask = feature_stds < eps;
    if sum(zero_std_mask) > 0
        fprintf('  WARNING: %d features have zero standard deviation (constant features)\n', sum(zero_std_mask));
        fprintf('    Features: %s\n', strjoin(available_features(zero_std_mask), ', '));
        % Set std to 1 to avoid division by zero (feature will be centered but not scaled)
        feature_stds(zero_std_mask) = 1;
    end
    
    % Normalize training data
    X_train_normalized = (X_train - feature_means) ./ feature_stds;
    
    % Normalize test data using training statistics (important!)
    X_test_normalized = (X_test - feature_means) ./ feature_stds;
    
    % Replace original data with normalized data
    X_train = X_train_normalized;
    X_test = X_test_normalized;
    
    % Debug: Print normalization statistics
    fprintf('Normalization statistics (from training data):\n');
    for f = 1:length(available_features)
        fprintf('  %s: mean=%.4f, std=%.4f\n', available_features{f}, feature_means(f), feature_stds(f));
    end
    
    % Debug: Verify normalization (training data should have mean≈0, std≈1)
    fprintf('\nVerification - Training data after normalization:\n');
    fprintf('  Mean per feature (should be ~0): %s\n', mat2str(mean(X_train, 1), 4));
    fprintf('  Std per feature (should be ~1): %s\n', mat2str(std(X_train, 0, 1), 4));
    
    % Split train into train and validation for hyperparameter tuning
    n_train = size(X_train, 1);
    n_val = round(n_train * validation_fraction);
    val_idx = randperm(n_train, n_val);
    train_idx = setdiff(1:n_train, val_idx);
    
    X_train_split = X_train(train_idx, :);
    y_train_split = y_train(train_idx);
    X_val_split = X_train(val_idx, :);
    y_val_split = y_train(val_idx);
    
    fprintf('Train split: %d samples (%.1f%%)\n', length(train_idx), 100*length(train_idx)/n_train);
    fprintf('Validation split: %d samples (%.1f%%)\n', length(val_idx), 100*length(val_idx)/n_train);
    
    % Debug: Data split statistics
    fprintf('\nData split target statistics:\n');
    fprintf('  Train split - Mean: %.4f, Std: %.4f\n', mean(y_train_split), std(y_train_split));
    fprintf('  Val split - Mean: %.4f, Std: %.4f\n', mean(y_val_split), std(y_val_split));
    
    %% Hyperparameter Tuning (Grid Search)
    fprintf('\n=== Hyperparameter Tuning (Grid Search) ===\n');
    tuning_start_time = tic;
    
    best_score = Inf;
    best_params = struct();
    best_model = [];
    
    % Generate all parameter combinations
    param_combinations = generateParamCombinations(hyperparams);
    fprintf('Total parameter combinations: %d\n', length(param_combinations));
    
    % Limit search space if too large (random sample)
    max_combinations = 50;  % Limit to 50 combinations for speed
    if length(param_combinations) > max_combinations
        fprintf('Sampling %d random combinations from %d total\n', ...
            max_combinations, length(param_combinations));
        sample_idx = randperm(length(param_combinations), max_combinations);
        param_combinations = param_combinations(sample_idx);
    end
    
    % Debug: Store all parameter combinations and their scores
    all_param_scores = [];
    
    for comb_idx = 1:length(param_combinations)
        params = param_combinations{comb_idx};
        
        fprintf('  Testing combination %d/%d: ', comb_idx, length(param_combinations));
        fprintf('depth=%d, lr=%.3f, n_est=%d, subsample=%.1f, colsample=%.1f\n', ...
            params.max_depth, params.learning_rate, params.n_estimators, ...
            params.subsample, params.colsample_bytree);
        
        % Train model with these parameters
        try
            model = trainXGBoost(X_train_split, y_train_split, X_val_split, y_val_split, ...
                params, early_stopping_rounds);
            
            % Evaluate on validation set
            y_val_pred = predictXGBoost(model, X_val_split);
            val_rmse = sqrt(mean((y_val_split - y_val_pred).^2));
            val_mae = mean(abs(y_val_split - y_val_pred));
            val_mape = mean(abs((y_val_split - y_val_pred) ./ y_val_split)) * 100;
            val_r2 = 1 - sum((y_val_split - y_val_pred).^2) / sum((y_val_split - mean(y_val_split)).^2);
            
            % Store score
            all_param_scores = [all_param_scores; struct('params', params, 'rmse', val_rmse, ...
                'mae', val_mae, 'mape', val_mape, 'r2', val_r2)];
            
            if val_rmse < best_score
                best_score = val_rmse;
                best_params = params;
                best_model = model;
                fprintf('    -> New best! RMSE: %.4f, MAE: %.4f, MAPE: %.2f%%, R²: %.4f\n', ...
                    val_rmse, val_mae, val_mape, val_r2);
            else
                fprintf('    -> RMSE: %.4f, MAE: %.4f, MAPE: %.2f%%, R²: %.4f\n', ...
                    val_rmse, val_mae, val_mape, val_r2);
            end
        catch ME
            fprintf('    -> Error: %s\n', ME.message);
            continue;
        end
    end
    
    % Debug: Print top 5 parameter combinations
    if ~isempty(all_param_scores)
        [~, sort_idx] = sort([all_param_scores.rmse]);
        fprintf('\n=== Top 5 Parameter Combinations (by Validation RMSE) ===\n');
        fprintf('Rank | Depth | LR    | N_Est | RMSE    | MAE     | MAPE    | R²      \n');
        fprintf('-----|-------|-------|-------|---------|---------|---------|---------\n');
        for i = 1:min(5, length(all_param_scores))
            idx = sort_idx(i);
            p = all_param_scores(idx).params;
            fprintf('%4d | %5d | %5.3f | %5d | %8.4f | %7.4f | %7.2f%% | %7.4f\n', ...
                i, p.max_depth, p.learning_rate, p.n_estimators, ...
                all_param_scores(idx).rmse, all_param_scores(idx).mae, ...
                all_param_scores(idx).mape, all_param_scores(idx).r2);
        end
    end
    
    tuning_time = toc(tuning_start_time);
    fprintf('\n=== Hyperparameter Tuning Summary ===\n');
    fprintf('Total combinations tested: %d\n', length(param_combinations));
    fprintf('Tuning time: %.2f seconds (%.2f minutes)\n', tuning_time, tuning_time/60);
    
    if isempty(best_params) || isempty(fieldnames(best_params))
        fprintf('\n[ERROR] No valid model was trained. All parameter combinations failed.\n');
        fprintf('Possible causes:\n');
        fprintf('  1. XGBoost DLL not loaded (compiler required for loadlibrary)\n');
        fprintf('  2. DLL path incorrect\n');
        fprintf('  3. Missing dependencies\n');
        fprintf('\nSolution: Install MinGW-w64 compiler or use LSBoost instead.\n');
        error('Hyperparameter tuning failed. No valid models were trained.');
    end
    
    fprintf('Best parameters found:\n');
    fprintf('  max_depth: %d\n', best_params.max_depth);
    fprintf('  learning_rate: %.4f\n', best_params.learning_rate);
    fprintf('  n_estimators: %d\n', best_params.n_estimators);
    fprintf('  subsample: %.2f\n', best_params.subsample);
    fprintf('  colsample_bytree: %.2f\n', best_params.colsample_bytree);
    fprintf('Best validation RMSE: %.4f\n', best_score);
    
    %% Train final model on full training set
    fprintf('\n=== Training Final Model ===\n');
    training_start_time = tic;
    final_model = trainXGBoost(X_train, y_train, X_val_split, y_val_split, ...
        best_params, early_stopping_rounds);
    training_time = toc(training_start_time);
    fprintf('Model training completed in %.2f seconds (%.2f minutes)\n', training_time, training_time/60);
    
    % Debug: Model information
    fprintf('\n=== Model Information (Fold %d, Target: %s) ===\n', fold_num, targetVar);
    fprintf('Model type: %s\n', final_model.type);
    
    if strcmp(final_model.type, 'xgboost_dll')
        % XGBoost DLL model structure
        fprintf('XGBoost DLL Model:\n');
        if isfield(final_model, 'iters_optimal')
            fprintf('  Optimal iterations: %d\n', final_model.iters_optimal);
        end
        if isfield(final_model, 'params')
            fprintf('  Objective: %s\n', final_model.params.objective);
            fprintf('  Booster: %s\n', final_model.params.booster);
        end
        fprintf('Best hyperparameters:\n');
        fprintf('  max_depth: %d\n', best_params.max_depth);
        fprintf('  learning_rate (eta): %.4f\n', best_params.learning_rate);
        fprintf('  n_estimators: %d\n', best_params.n_estimators);
        fprintf('  subsample: %.2f\n', best_params.subsample);
        fprintf('  colsample_bytree: %.2f\n', best_params.colsample_bytree);
    elseif strcmp(final_model.type, 'fitrensemble')
        % LSBoost model structure
        fprintf('LSBoost Model (fitrensemble):\n');
        fprintf('  Number of trees: %d\n', final_model.model.NumTrained);
        fprintf('  Learning rate: %.4f\n', final_model.model.ModelParameters.LearnRate);
        fprintf('  Method: %s\n', final_model.model.ModelParameters.Method);
        fprintf('Best hyperparameters:\n');
        fprintf('  max_depth: %d\n', best_params.max_depth);
        fprintf('  learning_rate: %.4f\n', best_params.learning_rate);
        fprintf('  n_estimators: %d\n', best_params.n_estimators);
        fprintf('  subsample: %.2f\n', best_params.subsample);
        fprintf('  colsample_bytree: %.2f\n', best_params.colsample_bytree);
    else
        fprintf('Unknown model type\n');
    end
    
    %% Generate Training Visualizations
    fprintf('\n=== Generating Training Visualizations ===\n');
    figuresDir = fullfile(outputDir, 'figures', 'Training');
    if ~exist(figuresDir, 'dir')
        mkdir(figuresDir);
    end
    
    % Feature Importance Plot
    feature_importance = getFeatureImportance(final_model);
    if ~isempty(feature_importance) && length(feature_importance) == length(available_features)
        [sorted_imp, idx] = sort(feature_importance, 'descend');
        sorted_names = available_features(idx);
        
        % Print feature importance table
        fprintf('\n=== Feature Importance (Fold %d) ===\n', fold_num);
        fprintf('Rank | Feature | Importance | Relative Importance\n');
        fprintf('-----|---------|------------|-------------------\n');
        total_imp = sum(sorted_imp);
        for i = 1:length(sorted_names)
            rel_imp = 100 * sorted_imp(i) / total_imp;
            fprintf('%4d | %-7s | %10.4f | %17.2f%%\n', i, sorted_names{i}, sorted_imp(i), rel_imp);
        end
        
        % 상위 15개 중요도 시각화
        n_show = min(15, length(sorted_names));
        fig = figure('Name', sprintf('Feature Importance - Fold %d', fold_num), ...
            'Position', [100, 100, 800, 600], 'Visible', 'on');
        barh(sorted_imp(1:n_show));
        yticks(1:n_show);
        yticklabels(strrep(sorted_names(1:n_show), '_', '\_'));
        xlabel('Importance Score', 'FontSize', 12);
        ylabel('Features', 'FontSize', 12);
        title(sprintf('XGBoost Feature Importance (Fold %d)', fold_num), 'FontSize', 14);
        grid on;
        savefig(fullfile(figuresDir, sprintf('Feature_Imp_Fold%d.fig', fold_num)));
        saveas(fig, fullfile(figuresDir, sprintf('Feature_Imp_Fold%d.png', fold_num)));
        close(fig);
        fprintf('  Saved: Feature_Imp_Fold%d.fig\n', fold_num);
    else
        if strcmp(final_model.type, 'xgboost_dll')
            fprintf('  INFO: XGBoost DLL model - feature importance not directly available\n');
            fprintf('    XGBoost C API does not provide feature importance extraction.\n');
            fprintf('    To get feature importance, use permutation importance or similar methods.\n');
        else
            fprintf('  WARNING: Feature importance extraction failed or dimension mismatch\n');
            fprintf('    Expected: %d features, Got: %d importance values\n', ...
                length(available_features), length(feature_importance));
        end
    end
    
    %% Evaluate on test set
    fprintf('\n=== Evaluating on Test Set ===\n');
    prediction_start_time = tic;
    y_test_pred = predictXGBoost(final_model, X_test);
    prediction_time = toc(prediction_start_time);
    fprintf('Prediction completed in %.4f seconds\n', prediction_time);
    
    % Calculate metrics
    residuals = y_test - y_test_pred;
    mae = mean(abs(residuals));
    rmse = sqrt(mean(residuals.^2));
    mape = mean(abs(residuals ./ y_test)) * 100;
    r2 = 1 - sum(residuals.^2) / sum((y_test - mean(y_test)).^2);
    
    % Additional metrics
    mean_actual = mean(y_test);
    mean_pred = mean(y_test_pred);
    bias = mean_pred - mean_actual;
    relative_bias = 100 * bias / mean_actual;
    
    fprintf('\n=== Test Set Performance (Fold %d, Target: %s) ===\n', fold_num, targetVar);
    fprintf('Sample size: %d\n', length(y_test));
    fprintf('Metrics:\n');
    fprintf('  MAE:  %.4f\n', mae);
    fprintf('  RMSE: %.4f\n', rmse);
    fprintf('  MAPE: %.2f%%\n', mape);
    fprintf('  R²:   %.4f\n', r2);
    fprintf('  Bias: %.4f (%.2f%% relative)\n', bias, relative_bias);
    
    % Debug: Additional performance metrics
    fprintf('\nDetailed Performance Metrics:\n');
    fprintf('  MAE relative to mean: %.2f%%\n', 100 * mae / mean_actual);
    fprintf('  RMSE relative to mean: %.2f%%\n', 100 * rmse / mean_actual);
    fprintf('  Explained Variance: %.2f%%\n', 100 * r2);
    if r2 < 0
        fprintf('  [WARNING] R² < 0: Model performs worse than predicting the mean!\n');
    elseif r2 < 0.3
        fprintf('  [WARNING] R² < 0.3: Weak model performance\n');
    elseif r2 < 0.7
        fprintf('  [INFO] 0.3 ≤ R² < 0.7: Moderate model performance\n');
    else
        fprintf('  [INFO] R² ≥ 0.7: Strong model performance\n');
    end
    
    fprintf('\nTarget statistics:\n');
    fprintf('  Actual - Mean: %.4f, Std: %.4f, Range: [%.4f, %.4f]\n', ...
        mean_actual, std(y_test), min(y_test), max(y_test));
    fprintf('  Predicted - Mean: %.4f, Std: %.4f, Range: [%.4f, %.4f]\n', ...
        mean_pred, std(y_test_pred), min(y_test_pred), max(y_test_pred));
    fprintf('  Variance Ratio (Pred/Actual): %.4f\n', var(y_test_pred) / var(y_test));
    if var(y_test_pred) / var(y_test) < 0.5
        fprintf('    [WARNING] Predicted variance is much lower than actual (underfitting risk)\n');
    elseif var(y_test_pred) / var(y_test) > 1.5
        fprintf('    [WARNING] Predicted variance is much higher than actual (overfitting risk)\n');
    end
    
    fprintf('\nResidual statistics:\n');
    fprintf('  Mean: %.4f (should be close to 0)\n', mean(residuals));
    fprintf('  Std: %.4f\n', std(residuals));
    fprintf('  Min: %.4f, Max: %.4f\n', min(residuals), max(residuals));
    fprintf('  Skewness: %.4f (0 = symmetric, >0 = right-skewed, <0 = left-skewed)\n', skewness(residuals));
    fprintf('  Kurtosis: %.4f (3 = normal, >3 = heavy-tailed, <3 = light-tailed)\n', kurtosis(residuals));
    
    % Debug: Residual distribution analysis
    residual_percentiles = prctile(residuals, [5, 25, 50, 75, 95]);
    fprintf('  Percentiles: 5th=%.4f, 25th=%.4f, 50th=%.4f, 75th=%.4f, 95th=%.4f\n', ...
        residual_percentiles(1), residual_percentiles(2), residual_percentiles(3), ...
        residual_percentiles(4), residual_percentiles(5));
    
    % Debug: Outlier detection in residuals
    iqr_residuals = iqr(residuals);
    q1_residuals = prctile(residuals, 25);
    q3_residuals = prctile(residuals, 75);
    outlier_lower = q1_residuals - 1.5 * iqr_residuals;
    outlier_upper = q3_residuals + 1.5 * iqr_residuals;
    n_outliers = sum(residuals < outlier_lower | residuals > outlier_upper);
    fprintf('  Outliers (IQR method): %d (%.2f%%)\n', n_outliers, 100*n_outliers/length(residuals));
    
    % SOC별 성능 분석
    if exist('soc_test', 'var') && ~isempty(soc_test) && length(soc_test) == length(y_test)
        fprintf('\n=== Performance by SOC (Fold %d, Target: %s) ===\n', fold_num, targetVar);
        unique_soc_test = unique(soc_test(~isnan(soc_test)));
        fprintf('SOC | N    | MAE     | RMSE    | MAPE    | R²      \n');
        fprintf('----|------|---------|---------|---------|---------\n');
        for s = 1:length(unique_soc_test)
            soc_val = unique_soc_test(s);
            soc_mask = (soc_test == soc_val);
            if sum(soc_mask) > 0
                soc_y_test = y_test(soc_mask);
                soc_y_pred = y_test_pred(soc_mask);
                soc_residuals = soc_y_test - soc_y_pred;
                soc_mae = mean(abs(soc_residuals));
                soc_rmse = sqrt(mean(soc_residuals.^2));
                soc_mape = mean(abs(soc_residuals ./ soc_y_test)) * 100;
                soc_r2 = 1 - sum(soc_residuals.^2) / sum((soc_y_test - mean(soc_y_test)).^2);
                fprintf(' %2d | %4d | %7.4f | %7.4f | %7.2f%% | %7.4f\n', ...
                    soc_val, sum(soc_mask), soc_mae, soc_rmse, soc_mape, soc_r2);
            end
        end
    end
    
    % Get feature importance
    feature_importance = getFeatureImportance(final_model);
    
    %% Store results
    fold_results.(sprintf('Fold%d', fold_num)) = struct();
    fold_results.(sprintf('Fold%d', fold_num)).best_params = best_params;
    fold_results.(sprintf('Fold%d', fold_num)).test_mae = mae;
    fold_results.(sprintf('Fold%d', fold_num)).test_rmse = rmse;
    fold_results.(sprintf('Fold%d', fold_num)).test_mape = mape;
    fold_results.(sprintf('Fold%d', fold_num)).test_r2 = r2;
    fold_results.(sprintf('Fold%d', fold_num)).test_bias = bias;
    fold_results.(sprintf('Fold%d', fold_num)).test_relative_bias = relative_bias;
    fold_results.(sprintf('Fold%d', fold_num)).tuning_time = tuning_time;
    fold_results.(sprintf('Fold%d', fold_num)).training_time = training_time;
    fold_results.(sprintf('Fold%d', fold_num)).prediction_time = prediction_time;
    fold_results.(sprintf('Fold%d', fold_num)).y_test = y_test;
    fold_results.(sprintf('Fold%d', fold_num)).y_test_pred = y_test_pred;
    fold_results.(sprintf('Fold%d', fold_num)).feature_importance = feature_importance;
    fold_results.(sprintf('Fold%d', fold_num)).available_features = available_features;
    
    % Store SOC data for evaluation
    % Note: soc_test is already filtered to match y_test size after NaN removal
    if exist('soc_test', 'var') && ~isempty(soc_test) && length(soc_test) == length(y_test)
        fold_results.(sprintf('Fold%d', fold_num)).soc_test = soc_test;
        % testData is also already filtered, but we need to ensure it matches
        if exist('testData', 'var') && height(testData) == length(y_test)
            fold_results.(sprintf('Fold%d', fold_num)).testData = testData;
        end
    end
    
    % Store predictions for overall evaluation
    all_predictions = [all_predictions; y_test_pred];
    all_actuals = [all_actuals; y_test];
    all_fold_ids = [all_fold_ids; repmat(fold_num, length(y_test), 1)];
    
    % Save individual fold model (target-specific filename)
    savePath_model = fullfile(outputDir, sprintf('XGBoost_%s_Fold%d_Model.mat', targetVar, fold_num));
    save(savePath_model, 'final_model', 'best_params', 'available_features', ...
        'mae', 'rmse', 'mape', 'r2', 'y_test', 'y_test_pred', 'feature_importance', ...
        'fold_num', 'targetVar');
    fprintf('Saved model: %s\n', savePath_model);
    end  % End of fold loop
    
    %% Overall CV Results for this target
    fprintf('\n========================================\n');
    fprintf('=== Overall Cross-Validation Results (%s) ===\n', targetVar);
    fprintf('========================================\n');

% Calculate mean and std across folds
fold_mae = [];
fold_rmse = [];
fold_mape = [];
fold_r2 = [];

for fold_num = 1:n_folds
    if isfield(fold_results, sprintf('Fold%d', fold_num))
        fold_mae = [fold_mae; fold_results.(sprintf('Fold%d', fold_num)).test_mae];
        fold_rmse = [fold_rmse; fold_results.(sprintf('Fold%d', fold_num)).test_rmse];
        fold_mape = [fold_mape; fold_results.(sprintf('Fold%d', fold_num)).test_mape];
        fold_r2 = [fold_r2; fold_results.(sprintf('Fold%d', fold_num)).test_r2];
    end
end

    fprintf('\n=== Cross-Validation Performance Summary (%s) ===\n', targetVar);
    fprintf('Number of folds completed: %d / %d\n', length(fold_mae), n_folds);
    fprintf('\nPerformance Metrics (Mean ± Std across folds):\n');
    fprintf('  MAE:  %.4f ± %.4f\n', mean(fold_mae), std(fold_mae));
    fprintf('  RMSE: %.4f ± %.4f\n', mean(fold_rmse), std(fold_rmse));
    fprintf('  MAPE: %.2f%% ± %.2f%%\n', mean(fold_mape), std(fold_mape));
    fprintf('  R²:   %.4f ± %.4f\n', mean(fold_r2), std(fold_r2));
    
    % Debug: Detailed performance analysis
    fprintf('\n=== Detailed Performance Analysis (%s) ===\n', targetVar);
    fprintf('Performance Stability:\n');
    fprintf('  MAE CV (Coefficient of Variation): %.2f%%\n', 100*std(fold_mae)/mean(fold_mae));
    fprintf('  RMSE CV: %.2f%%\n', 100*std(fold_rmse)/mean(fold_rmse));
    fprintf('  R² Range: [%.4f, %.4f]\n', min(fold_r2), max(fold_r2));
    fprintf('  R² Spread: %.4f\n', max(fold_r2) - min(fold_r2));

    % Per-fold performance table
    fprintf('\nPer-fold Performance (%s):\n', targetVar);
    fprintf('Fold | MAE     | RMSE    | MAPE    | R²      | Bias    | Best Depth\n');
    fprintf('-----|---------|---------|---------|---------|---------|------------\n');
    for fold_num = 1:n_folds
        if isfield(fold_results, sprintf('Fold%d', fold_num))
            fold_data = fold_results.(sprintf('Fold%d', fold_num));
            best_depth = fold_data.best_params.max_depth;
            fprintf('  %d  | %7.4f | %7.4f | %7.2f%% | %7.4f | %7.4f | %d\n', ...
                fold_num, fold_data.test_mae, fold_data.test_rmse, ...
                fold_data.test_mape, fold_data.test_r2, fold_data.test_bias, best_depth);
        end
    end

    % Timing statistics
    if isfield(fold_results, 'Fold1')
        total_tuning_time = 0;
        total_training_time = 0;
        total_prediction_time = 0;
        for fold_num = 1:n_folds
            fold_name = sprintf('Fold%d', fold_num);
            if isfield(fold_results, fold_name)
                if isfield(fold_results.(fold_name), 'tuning_time')
                    total_tuning_time = total_tuning_time + fold_results.(fold_name).tuning_time;
                end
                if isfield(fold_results.(fold_name), 'training_time')
                    total_training_time = total_training_time + fold_results.(fold_name).training_time;
                end
                if isfield(fold_results.(fold_name), 'prediction_time')
                    total_prediction_time = total_prediction_time + fold_results.(fold_name).prediction_time;
                end
            end
        end
        fprintf('\nTiming Statistics (%s):\n', targetVar);
        fprintf('  Total tuning time: %.2f seconds (%.2f minutes)\n', total_tuning_time, total_tuning_time/60);
        fprintf('  Total training time: %.2f seconds (%.2f minutes)\n', total_training_time, total_training_time/60);
        fprintf('  Total prediction time: %.2f seconds\n', total_prediction_time);
        fprintf('  Average per fold - Tuning: %.2f sec, Training: %.2f sec, Prediction: %.4f sec\n', ...
            total_tuning_time/length(fold_mae), total_training_time/length(fold_mae), total_prediction_time/length(fold_mae));
    end

    % Overall metrics (all test data combined)
    overall_mae = mean(abs(all_actuals - all_predictions));
    overall_rmse = sqrt(mean((all_actuals - all_predictions).^2));
    overall_mape = mean(abs((all_actuals - all_predictions) ./ all_actuals)) * 100;
    overall_r2 = 1 - sum((all_actuals - all_predictions).^2) / sum((all_actuals - mean(all_actuals)).^2);
    overall_bias = mean(all_predictions - all_actuals);
    overall_relative_bias = 100 * overall_bias / mean(all_actuals);
    
    fprintf('\nOverall Performance - All Test Data Combined (%s):\n', targetVar);
    fprintf('  MAE:  %.4f\n', overall_mae);
    fprintf('  RMSE: %.4f\n', overall_rmse);
    fprintf('  MAPE: %.2f%%\n', overall_mape);
    fprintf('  R²:   %.4f\n', overall_r2);
    fprintf('  Bias: %.4f (%.2f%% relative)\n', overall_bias, overall_relative_bias);
    fprintf('  Sample size: %d\n', length(all_actuals));
    
    % Debug: Residual analysis
    overall_residuals = all_actuals - all_predictions;
    fprintf('\nResidual Analysis (%s):\n', targetVar);
    fprintf('  Residual Mean: %.4f (should be close to 0)\n', mean(overall_residuals));
    fprintf('  Residual Std: %.4f\n', std(overall_residuals));
    fprintf('  Residual Range: [%.4f, %.4f]\n', min(overall_residuals), max(overall_residuals));
    fprintf('  Residual Skewness: %.4f\n', skewness(overall_residuals));
    fprintf('  Residual Kurtosis: %.4f\n', kurtosis(overall_residuals));
    
    % Debug: Prediction vs Actual statistics
    fprintf('\nPrediction vs Actual Statistics (%s):\n', targetVar);
    fprintf('  Actual - Mean: %.4f, Std: %.4f, Range: [%.4f, %.4f]\n', ...
        mean(all_actuals), std(all_actuals), min(all_actuals), max(all_actuals));
    fprintf('  Predicted - Mean: %.4f, Std: %.4f, Range: [%.4f, %.4f]\n', ...
        mean(all_predictions), std(all_predictions), min(all_predictions), max(all_predictions));
    fprintf('  Variance Ratio (Pred/Actual): %.4f\n', var(all_predictions) / var(all_actuals));
    
    %% Save results for this target
    savePath_results = fullfile(outputDir, sprintf('XGBoost_%s_Training_Results.mat', targetVar));
    if exist('targetEventType', 'var')
        save(savePath_results, 'fold_results', 'all_predictions', 'all_actuals', 'all_fold_ids', ...
            'fold_mae', 'fold_rmse', 'fold_mape', 'fold_r2', ...
            'overall_mae', 'overall_rmse', 'overall_mape', 'overall_r2', ...
            'overall_bias', 'overall_relative_bias', ...
            'targetVar', 'targetEventType', 'n_folds', 'hyperparams', ...
            'available_features', 'selected_features');
    else
        save(savePath_results, 'fold_results', 'all_predictions', 'all_actuals', 'all_fold_ids', ...
            'fold_mae', 'fold_rmse', 'fold_mape', 'fold_r2', ...
            'overall_mae', 'overall_rmse', 'overall_mape', 'overall_r2', ...
            'overall_bias', 'overall_relative_bias', ...
            'targetVar', 'n_folds', 'hyperparams', ...
            'available_features', 'selected_features');
    end
    
    fprintf('\nSaved results: %s\n', savePath_results);
    
    % Store results for overall summary
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
    
    fprintf('\n=== Training Complete for Target: %s ===\n', targetVar);
end  % End of target variable loop

%% Overall Summary (All Targets)
fprintf('\n');
fprintf('########################################################################\n');
fprintf('### Overall Summary - All Target Variables ###\n');
fprintf('########################################################################\n');

fprintf('\n=== Performance Comparison Across Targets ===\n');
fprintf('Target Variable | MAE (Mean±Std) | RMSE (Mean±Std) | MAPE (Mean±Std) | R² (Mean±Std)\n');
fprintf('----------------|-----------------|------------------|------------------|---------------\n');
for targetIdx = 1:length(targetVars)
    targetVar = targetVars{targetIdx};
    if isfield(all_target_results, targetVar)
        result = all_target_results.(targetVar);
        fprintf('%-15s | %6.4f±%6.4f | %7.4f±%7.4f | %6.2f%%±%5.2f%% | %7.4f±%7.4f\n', ...
            targetVar, result.fold_mae_mean, result.fold_mae_std, ...
            result.fold_rmse_mean, result.fold_rmse_std, ...
            result.overall_mape, 0, ...  % MAPE std not calculated per fold
            result.fold_r2_mean, result.fold_r2_std);
    end
end

fprintf('\n=== Overall Test Performance (All Targets Combined) ===\n');
for targetIdx = 1:length(targetVars)
    targetVar = targetVars{targetIdx};
    if isfield(all_target_results, targetVar)
        result = all_target_results.(targetVar);
        fprintf('\n%s:\n', targetVar);
        fprintf('  MAE:  %.4f\n', result.overall_mae);
        fprintf('  RMSE: %.4f\n', result.overall_rmse);
        fprintf('  MAPE: %.2f%%\n', result.overall_mape);
        fprintf('  R²:   %.4f\n', result.overall_r2);
        fprintf('  Bias: %.4f (%.2f%% relative)\n', result.overall_bias, result.overall_relative_bias);
    end
end

%% Save overall results
savePath_overall = fullfile(outputDir, 'XGBoost_AllTargets_Summary.mat');
if exist('targetEventType', 'var')
    save(savePath_overall, 'all_target_results', 'targetVars', 'targetEventType', ...
        'n_folds', 'hyperparams', 'selected_features');
else
    save(savePath_overall, 'all_target_results', 'targetVars', ...
        'n_folds', 'hyperparams', 'selected_features');
end
fprintf('\nSaved overall summary: %s\n', savePath_overall);

fprintf('\n=== XGBoost Training Complete (All Targets) ===\n');
if exist('targetEventType', 'var')
    fprintf('Model type: %s\n', targetEventType);
    fprintf('Target variables trained: %d\n', length(targetVars));
    fprintf('  %s\n', strjoin(targetVars, ', '));
    if ~strcmp(targetEventType, 'All')
        fprintf('\n[INFO] This is a %s-only model.\n', targetEventType);
        fprintf('       To compare with other event types, run the script again with different targetEventType.\n');
        fprintf('       Options: ''Charge'', ''Discharge'', or ''All''\n');
    end
end

%% Helper Functions

function combinations = generateParamCombinations(hyperparams)
    % Generate all combinations of hyperparameters
    fields = fieldnames(hyperparams);
    values = cell(length(fields), 1);
    
    for i = 1:length(fields)
        values{i} = hyperparams.(fields{i});
    end
    
    % Create meshgrid-like combinations
    [grids{1:length(fields)}] = ndgrid(values{:});
    
    combinations = {};
    n_comb = numel(grids{1});
    
    for comb = 1:n_comb
        params = struct();
        for i = 1:length(fields)
            idx = cell(1, length(fields));
            [idx{:}] = ind2sub(size(grids{1}), comb);
            params.(fields{i}) = grids{i}(idx{i});
        end
        combinations{comb} = params;
    end
end

function model = trainXGBoost(X_train, y_train, X_val, y_val, params, early_stopping)
    % Train XGBoost model using xgboost_train.m (XGBoost DLL)
    % This function wraps xgboost_train.m for regression
    
    % Add path to xgboost functions if not already added
    xgboost_path = 'C:\Users\Chulwon Jung\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\Functions to run xgboost in Matlab';
    if ~contains(path, xgboost_path)
        addpath(xgboost_path);
    end
    
    % Try to use XGBoost DLL, fallback to LSBoost if it fails
    use_xgboost_dll = true;
    
    try
        % Convert parameters to xgboost format
        xgb_params = struct();
        xgb_params.booster = 'gbtree';
        xgb_params.objective = 'reg:squarederror';  % Regression objective
        xgb_params.max_depth = params.max_depth;
        xgb_params.eta = params.learning_rate;  % learning_rate in XGBoost is called 'eta'
        xgb_params.min_child_weight = 1;
        xgb_params.subsample = params.subsample;
        xgb_params.colsample_bytree = params.colsample_bytree;
        xgb_params.num_parallel_tree = 1;
        
        % Set max number of iterations
        max_num_iters = params.n_estimators;
        
        % For regression, use 'None' eval_metric (no internal CV)
        % We handle CV externally in the main script
        eval_metric = 'None';
        
        % Don't save model to file (keep in memory)
        model_filename = [];
        
        % Train XGBoost model (this will try to load DLL inside xgboost_train.m)
        model = xgboost_train(X_train, y_train, xgb_params, max_num_iters, eval_metric, model_filename);
        
        % Store original params for reference
        model.original_params = params;
        model.type = 'xgboost_dll';
        
    catch ME
        % If XGBoost DLL fails, fallback to LSBoost
        warning('XGBoost DLL failed: %s\nFalling back to LSBoost (MATLAB fitrensemble)', ME.message);
        use_xgboost_dll = false;
        
        % Fallback to LSBoost
        max_splits = min(2^params.max_depth - 1, 1000);
        t = templateTree('MaxNumSplits', max_splits, 'MinLeafSize', 5);
        ensemble_model = fitrensemble(X_train, y_train, ...
            'Method', 'LSBoost', ...
            'NumLearningCycles', params.n_estimators, ...
            'LearnRate', params.learning_rate, ...
            'Learners', t);
        
        model = struct();
        model.type = 'fitrensemble';
        model.model = ensemble_model;
        model.params = params;
    end
end

function y_pred = predictXGBoost(model, X)
    % Predict using XGBoost model (xgboost_test.m)
    if isstruct(model) && isfield(model, 'type')
        if strcmp(model.type, 'xgboost_dll')
            % XGBoost DLL prediction
            loadmodel = 0;  % Use model from memory (not from file)
            y_pred = xgboost_test(X, model, loadmodel);
        elseif strcmp(model.type, 'fitrensemble')
            % fitrensemble prediction (backward compatibility)
            y_pred = predict(model.model, X);
        else
            error('Unknown model type: %s', model.type);
        end
    elseif isa(model, 'RegressionEnsemble')
        % Direct RegressionEnsemble (backward compatibility)
        y_pred = predict(model, X);
    else
        error('Unknown model type');
    end
end

function importance = getFeatureImportance(model)
    % Get feature importance from model (XGBoost DLL or fitrensemble)
    if isstruct(model) && isfield(model, 'type')
        if strcmp(model.type, 'xgboost_dll')
            % XGBoost DLL model - feature importance not directly available
            % XGBoost C API doesn't provide easy feature importance extraction
            % Return empty or use params as proxy (not ideal, but informative)
            % Note: XGBoost DLL doesn't provide feature importance in the same way
            % For now, return empty array - feature importance would need to be
            % calculated separately using permutation importance or similar methods
            importance = [];
            % Alternative: could use model.params to infer, but not accurate
        elseif strcmp(model.type, 'fitrensemble')
            % fitrensemble importance (OOBPermutedPredictorDeltaError)
            try
                importance = model.model.OOBPermutedPredictorDeltaError;
                if isempty(importance) || all(isnan(importance))
                    % Fallback: use predictor importance
                    try
                        importance = predictorImportance(model.model);
                    catch
                        importance = [];
                    end
                end
            catch
                % Try predictorImportance as fallback
                try
                    importance = predictorImportance(model.model);
                catch
                    importance = [];
                end
            end
        else
            importance = [];
        end
    elseif isa(model, 'RegressionEnsemble')
        % Direct RegressionEnsemble (backward compatibility)
        try
            importance = model.OOBPermutedPredictorDeltaError;
            if isempty(importance)
                importance = predictorImportance(model);
            end
        catch
            importance = predictorImportance(model);
        end
    else
        importance = [];
    end
end

