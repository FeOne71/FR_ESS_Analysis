%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 04_DriveCycle_FeatureEngineering.m
% Feature Engineering & Selection
% 
% 목적: 
% - 결측치/이상치 처리
% - Target과의 상관분석
% - 다중공선성 해결 (VIF 계산)
% - Feature Selection (XGBoost Feature Importance)
% - 선택된 feature 리스트 저장
%
% 입력:
% - DriveCycle_Fold*_Train.mat (각 Fold의 Train 데이터)
%
% 출력:
% - Selected_Features.mat (선택된 feature 리스트)
% - Feature_Engineering_Results.mat (상관분석, VIF 결과 등)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

fprintf('=== Drive Cycle Feature Engineering ===\n');

%% Configuration - User Settings
% =========================================================================
inputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
outputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Target variables (Capacity + RPT DCIR values)
targetVars = {'Capacity_C3', 'Capacity_OCV'};
% Note: RPT variables will be added to targetVars after loading data

% Feature selection thresholds (Relaxed for more features)
missing_threshold = 0.95;     % Missing ratio > 95%: remove (relaxed - keep features available in some groups)
correlation_threshold = 0.01;  % |correlation| < 0.01: remove (very relaxed - keep most features)
vif_threshold = 100;          % VIF > 100: remove (very relaxed - only remove extreme multicollinearity)
use_vif_filtering = false;     % Skip VIF filtering (use correlation-based multicollinearity check only)
multicollinearity_threshold = 0.999;  % Correlation > 0.999: remove one (extremely relaxed - almost no removal)
min_features_to_keep = 10;    % Minimum number of features to keep (even if correlation is low)

% Final feature count
target_feature_count = 30;   % Aim for ~30 features
% =========================================================================

%% Load first fold's train data to get feature list
fprintf('\n=== Loading Data ===\n');
fold1_train_file = fullfile(inputDir, 'DriveCycle_Fold1_Train.mat');
if ~exist(fold1_train_file, 'file')
    error('Fold 1 train data not found. Please run DriveCycle_DataSplit_03.m first');
end

load(fold1_train_file, 'trainData');
fprintf('Loaded Fold 1 train data: %d rows, %d columns\n', height(trainData), width(trainData));

%% Identify feature columns (exclude metadata, target, and RPT data)
% [핵심 원칙] RPT 데이터는 Target/Reference이지 Feature가 아님 (Data Leakage 방지)
fprintf('\n=== Feature Identification ===\n');
allVars = trainData.Properties.VariableNames;

% Metadata variables (not used as features)
metadataVars = {'Cycle', 'Channel', 'DC_Profile', 'EventType', 'SOC'};  % SOC will be added later in training

% Target variables (what we want to predict)
targetVars = {'Capacity_C3', 'Capacity_OCV'};

% [핵심] RPT로 시작하는 모든 변수 찾기 (Target/Reference이므로 Feature에서 제외)
rptVars = allVars(startsWith(allVars, 'RPT_'));
fprintf('Found %d RPT variables (will be excluded from features)\n', length(rptVars));

% Add RPT variables to target variables (RPT-DCIR values are also targets)
targetVars = [targetVars, rptVars];
fprintf('Total target variables: %d (Capacity: 2, RPT-DCIR: %d)\n', length(targetVars), length(rptVars));

% Excluded variables: metadata + target + RPT data (RPT is now in targetVars, so don't double exclude)
excludedVars = [metadataVars, targetVars];

% Final feature candidates (only drive cycle resistance values: R_1s, R_3s, etc.)
featureVars = setdiff(allVars, excludedVars);

fprintf('\nTotal features: %d\n', length(featureVars));
fprintf('Metadata variables: %s\n', strjoin(metadataVars, ', '));
fprintf('Target variables: %s\n', strjoin(targetVars, ', '));
fprintf('RPT variables (excluded): %d variables\n', length(rptVars));
if length(rptVars) > 0 && length(rptVars) <= 10
    fprintf('  RPT variables: %s\n', strjoin(rptVars, ', '));
elseif length(rptVars) > 10
    fprintf('  RPT variables (first 10): %s\n', strjoin(rptVars(1:10), ', '));
end
fprintf('\n[INFO] RPT data is used as Target/Reference only, not as input features.\n');
fprintf('       This prevents data leakage and ensures realistic SOH estimation.\n');

%% Combine all folds' train data for feature engineering
fprintf('\n=== Combining All Folds Train Data ===\n');
allTrainData = trainData;

for fold_num = 2:4
    fold_train_file = fullfile(inputDir, sprintf('DriveCycle_Fold%d_Train.mat', fold_num));
    if exist(fold_train_file, 'file')
        load(fold_train_file, 'trainData');
        allTrainData = [allTrainData; trainData];
        fprintf('Added Fold %d: %d rows\n', fold_num, height(trainData));
    end
end

fprintf('Total combined train data: %d rows\n', height(allTrainData));

%% Step 1: Missing Value Analysis (Group-wise, same as outlier removal)
fprintf('\n=== Step 1: Missing Value Analysis (Group-wise) ===\n');

% Use same grouping as outlier removal
missing_group_by = {'Cycle', 'DC_Profile', 'SOC', 'EventType'};  % Grouping variables (DC별 * 사이클별 * SOC별 * 이벤트타입별)
fprintf('Grouping by: %s\n', strjoin(missing_group_by, ', '));

% Check if grouping variables exist
available_group_vars_missing = {};
for g = 1:length(missing_group_by)
    if ismember(missing_group_by{g}, allTrainData.Properties.VariableNames)
        available_group_vars_missing{end+1} = missing_group_by{g};
    else
        fprintf('WARNING: Grouping variable "%s" not found in data. Skipping.\n', missing_group_by{g});
    end
end

if isempty(available_group_vars_missing)
    fprintf('WARNING: No grouping variables available. Using global missing value analysis.\n');
    % Global missing value analysis (fallback)
    X_all = allTrainData{:, featureVars};
    missing_ratio = sum(isnan(X_all), 1) / size(X_all, 1);
    remove_missing = missing_ratio > missing_threshold;
    features_after_missing = featureVars(~remove_missing);
    
    fprintf('Features removed (missing > %.0f%%): %d\n', missing_threshold*100, sum(remove_missing));
    fprintf('Remaining features: %d\n', length(features_after_missing));
    
    if sum(remove_missing) > 0
        fprintf('Removed features: %s\n', strjoin(featureVars(remove_missing), ', '));
    end
else
    % Group-wise missing value analysis
    fprintf('\n--- Group-wise Missing Value Analysis ---\n');
    
    % Get unique groups
    group_data = allTrainData(:, available_group_vars_missing);
    [unique_groups, ~, group_idx] = unique(group_data, 'rows');
    n_groups = height(unique_groups);
    
    fprintf('Number of groups: %d\n', n_groups);
    
    % For each feature, calculate missing ratio per group and overall
    feature_missing_ratios = zeros(length(featureVars), n_groups);
    
    for f_idx = 1:length(featureVars)
        feat_name = featureVars{f_idx};
        feat_data = allTrainData.(feat_name);
        
        % Calculate missing ratio for each group
        for g = 1:n_groups
            group_mask = (group_idx == g);
            group_feat_data = feat_data(group_mask);
            
            if sum(group_mask) > 0
                group_missing_ratio = sum(isnan(group_feat_data)) / sum(group_mask);
                feature_missing_ratios(f_idx, g) = group_missing_ratio;
            else
                feature_missing_ratios(f_idx, g) = 1.0;  % No data = 100% missing
            end
        end
    end
    
    % For each feature, use the maximum missing ratio across all groups
    % (if a feature is missing in most groups, it should be removed)
    max_missing_ratio_per_feature = max(feature_missing_ratios, [], 2);
    
    % Also calculate minimum missing ratio (best case scenario)
    % This helps keep features that are available in at least some groups
    min_missing_ratio_per_feature = min(feature_missing_ratios, [], 2);
    
    % Also calculate overall missing ratio for comparison
    X_all = allTrainData{:, featureVars};
    overall_missing_ratio = sum(isnan(X_all), 1) / size(X_all, 1);
    
    % Remove features where:
    % 1. Maximum missing ratio > threshold (too many groups have missing data)
    % AND
    % 2. Minimum missing ratio > threshold (even best group has too much missing)
    % This keeps features that are available in at least some groups
    remove_missing = (max_missing_ratio_per_feature > missing_threshold) & ...
                     (min_missing_ratio_per_feature > missing_threshold);
    features_after_missing = featureVars(~remove_missing);
    
    fprintf('\nFeatures removed (max missing ratio across groups > %.0f%%): %d\n', missing_threshold*100, sum(remove_missing));
    fprintf('Remaining features: %d\n', length(features_after_missing));
    
    if sum(remove_missing) > 0
        fprintf('Removed features: %s\n', strjoin(featureVars(remove_missing), ', '));
        fprintf('\nDetailed missing ratios (removed features):\n');
        for f_idx = find(remove_missing)'
            fprintf('  %s: Overall=%.1f%%, Max across groups=%.1f%%\n', ...
                featureVars{f_idx}, 100*overall_missing_ratio(f_idx), 100*max_missing_ratio_per_feature(f_idx));
        end
    end
    
    if n_groups > 20
        fprintf('  (Too many groups to display individually)\n');
    end
end

%% Step 2: Outlier Detection (Group-wise IQR Method)
fprintf('\n=== Step 2: Outlier Detection (Group-wise IQR) ===\n');

% Configuration for outlier removal
outlier_method = 'IQR';  % 'IQR' or 'Zscore' or 'ModifiedZscore'
outlier_iqr_multiplier = 1.5;  % IQR multiplier (1.5 = standard, 2.0 = more conservative)
outlier_group_by = {'Cycle', 'DC_Profile', 'SOC', 'EventType'};  % Grouping variables (DC별 * 사이클별 * SOC별 * 이벤트타입별)
outlier_remove_enabled = true;  % Set to false to skip outlier removal

if ~outlier_remove_enabled
    fprintf('Outlier removal is disabled. Skipping.\n');
    features_after_outlier = features_after_missing;
    outlier_removed_count = 0;
else
    fprintf('Method: %s (Multiplier: %.1f)\n', outlier_method, outlier_iqr_multiplier);
    fprintf('Grouping by: %s\n', strjoin(outlier_group_by, ', '));
    
    % Check if grouping variables exist
    available_group_vars = {};
    for g = 1:length(outlier_group_by)
        if ismember(outlier_group_by{g}, allTrainData.Properties.VariableNames)
            available_group_vars{end+1} = outlier_group_by{g};
        else
            fprintf('WARNING: Grouping variable "%s" not found in data. Skipping.\n', outlier_group_by{g});
        end
    end
    
    if isempty(available_group_vars)
        fprintf('WARNING: No grouping variables available. Using global outlier removal.\n');
        available_group_vars = {};
    end
    
    % Get feature columns for outlier detection (only resistance features, not capacity)
    outlier_features = features_after_missing;
    % Remove capacity-related features from outlier detection (they are targets)
    outlier_features = outlier_features(~contains(outlier_features, 'Capacity'));
    
    fprintf('Features to check for outliers: %d\n', length(outlier_features));
    fprintf('  Examples: %s\n', strjoin(outlier_features(1:min(5, length(outlier_features))), ', '));
    
    % Initialize outlier mask (all false = no outliers initially)
    outlier_mask = false(height(allTrainData), 1);
    outlier_stats = struct();
    
    if isempty(available_group_vars)
        % Global outlier removal (no grouping)
        fprintf('\n--- Global Outlier Removal ---\n');
        
        for f_idx = 1:length(outlier_features)
            feat_name = outlier_features{f_idx};
            feat_data = allTrainData.(feat_name);
            valid_data = feat_data(~isnan(feat_data));
            
            if length(valid_data) < 10
                fprintf('  %s: Insufficient data (%d valid points). Skipping.\n', feat_name, length(valid_data));
                continue;
            end
            
            % Calculate IQR
            Q1 = prctile(valid_data, 25);
            Q3 = prctile(valid_data, 75);
            IQR = Q3 - Q1;
            lower_bound = Q1 - outlier_iqr_multiplier * IQR;
            upper_bound = Q3 + outlier_iqr_multiplier * IQR;
            
            % Find outliers
            feat_outliers = (feat_data < lower_bound) | (feat_data > upper_bound);
            feat_outliers(isnan(feat_data)) = false;  % Don't count NaN as outlier
            
            outlier_mask = outlier_mask | feat_outliers;
            
            n_outliers = sum(feat_outliers);
            if n_outliers > 0
                fprintf('  %s: %d outliers (%.1f%%) [Bounds: %.2f - %.2f]\n', ...
                    feat_name, n_outliers, 100*n_outliers/length(valid_data), lower_bound, upper_bound);
            end
        end
    else
        % Group-wise outlier removal
        fprintf('\n--- Group-wise Outlier Removal ---\n');
        
        % Get unique groups
        group_data = allTrainData(:, available_group_vars);
        [unique_groups, ~, group_idx] = unique(group_data, 'rows');
        n_groups = height(unique_groups);
        
        fprintf('Number of groups: %d\n', n_groups);
        
        % Process each group
        for g = 1:n_groups
            group_mask = (group_idx == g);
            group_data_subset = allTrainData(group_mask, :);
            n_group_samples = height(group_data_subset);
            
            if n_group_samples < 5
                % Too few samples in group, skip outlier detection
                continue;
            end
            
            % Create group label for display
            group_label_parts = {};
            for v = 1:length(available_group_vars)
                var_name = available_group_vars{v};
                var_value = unique_groups.(var_name)(g);
                if isnumeric(var_value)
                    group_label_parts{end+1} = sprintf('%s=%d', var_name, var_value);
                else
                    group_label_parts{end+1} = sprintf('%s=%s', var_name, string(var_value));
                end
            end
            group_label = strjoin(group_label_parts, ', ');
            
            % Check each feature in this group
            group_outlier_mask = false(n_group_samples, 1);
            
            for f_idx = 1:length(outlier_features)
                feat_name = outlier_features{f_idx};
                feat_data = group_data_subset.(feat_name);
                valid_data = feat_data(~isnan(feat_data));
                
                if length(valid_data) < 5
                    % Too few valid points in this group for this feature
                    continue;
                end
                
                % Calculate IQR for this group and feature
                Q1 = prctile(valid_data, 25);
                Q3 = prctile(valid_data, 75);
                IQR = Q3 - Q1;
                
                if IQR == 0
                    % All values are the same, skip
                    continue;
                end
                
                lower_bound = Q1 - outlier_iqr_multiplier * IQR;
                upper_bound = Q3 + outlier_iqr_multiplier * IQR;
                
                % Find outliers in this group for this feature
                feat_outliers = (feat_data < lower_bound) | (feat_data > upper_bound);
                feat_outliers(isnan(feat_data)) = false;  % Don't count NaN as outlier
                
                group_outlier_mask = group_outlier_mask | feat_outliers;
            end
            
            % Mark outliers in the main mask
            group_indices = find(group_mask);
            outlier_mask(group_indices(group_outlier_mask)) = true;
            
            n_group_outliers = sum(group_outlier_mask);
            if n_group_outliers > 0 && n_groups <= 20  % Only print if not too many groups
                fprintf('  Group %d (%s): %d outliers (%.1f%%)\n', ...
                    g, group_label, n_group_outliers, 100*n_group_outliers/n_group_samples);
            end
        end
        
        if n_groups > 20
            fprintf('  (Too many groups to display individually)\n');
        end
    end
    
    % Count total outliers
    outlier_removed_count = sum(outlier_mask);
    outlier_removed_ratio = 100 * outlier_removed_count / height(allTrainData);
    
    fprintf('\n--- Outlier Removal Summary ---\n');
    fprintf('Total outliers detected: %d (%.2f%%)\n', outlier_removed_count, outlier_removed_ratio);
    
    % Remove outliers from data
    if outlier_removed_count > 0
        allTrainData = allTrainData(~outlier_mask, :);
        fprintf('Removed %d outlier rows. Remaining rows: %d\n', ...
            outlier_removed_count, height(allTrainData));
        
        % Update feature list (features might have been removed if all values became NaN)
        X_after_outlier = allTrainData{:, features_after_missing};
        missing_ratio_after = sum(isnan(X_after_outlier), 1) / size(X_after_outlier, 1);
        features_after_outlier = features_after_missing(missing_ratio_after < 1.0);  % Keep features with at least some data
        
        if length(features_after_outlier) < length(features_after_missing)
            fprintf('WARNING: %d features removed due to 100%% missing after outlier removal\n', ...
                length(features_after_missing) - length(features_after_outlier));
        end
    else
        fprintf('No outliers detected. No rows removed.\n');
        features_after_outlier = features_after_missing;
    end
    
    % Store outlier statistics
    outlier_stats.method = outlier_method;
    outlier_stats.iqr_multiplier = outlier_iqr_multiplier;
    outlier_stats.grouping_vars = available_group_vars;
    outlier_stats.removed_count = outlier_removed_count;
    outlier_stats.removed_ratio = outlier_removed_ratio;
    outlier_stats.features_checked = outlier_features;
end

%% Step 3: Correlation Analysis with Target
fprintf('\n=== Step 3: Target Correlation Analysis ===\n');

% Handle missing values for correlation calculation
X_clean = allTrainData{:, features_after_outlier};
valid_rows = ~any(isnan(X_clean), 2);  % Rows with no missing values
X_clean = X_clean(valid_rows, :);

% Calculate correlation with each target
target_correlations = struct();
for t = 1:length(targetVars)
    targetVar = targetVars{t};
    y = allTrainData.(targetVar)(valid_rows);
    
    % Remove NaN from target
    valid_target = ~isnan(y);
    X_valid = X_clean(valid_target, :);
    y_valid = y(valid_target);
    
    if sum(valid_target) > 0
        % Calculate correlation
        corr_vals = zeros(size(X_valid, 2), 1);
        for f = 1:size(X_valid, 2)
            x_feat = X_valid(:, f);
            valid_feat = ~isnan(x_feat);
            if sum(valid_feat) > 10  % Need at least 10 valid points
                corr_vals(f) = corr(x_feat(valid_feat), y_valid(valid_feat), 'rows', 'complete');
            else
                corr_vals(f) = 0;
            end
        end
        
        target_correlations.(targetVar) = corr_vals;
        
        % Remove features with low correlation
        low_corr = abs(corr_vals) < correlation_threshold;
        fprintf('\nTarget: %s\n', targetVar);
        fprintf('Features with |correlation| < %.2f: %d\n', correlation_threshold, sum(low_corr));
        
        % Display top correlations
        [sorted_corr, idx] = sort(abs(corr_vals), 'descend');
        fprintf('Top 10 feature correlations:\n');
        for i = 1:min(10, length(sorted_corr))
            fprintf('  %d. %s: %.4f\n', i, features_after_outlier{idx(i)}, corr_vals(idx(i)));
        end
    end
end

% Keep features that have correlation > threshold with at least one target
if length(targetVars) > 1
    max_corr = max([abs(target_correlations.(targetVars{1})), ...
                    abs(target_correlations.(targetVars{2}))], [], 2);
else
    max_corr = abs(target_correlations.(targetVars{1}));
end

keep_corr = max_corr >= correlation_threshold;
features_after_corr = features_after_outlier(keep_corr);

fprintf('\nFeatures removed (low correlation): %d\n', sum(~keep_corr));
fprintf('Remaining features: %d\n', length(features_after_corr));

% Ensure minimum number of features (keep top N by correlation if too few)
if length(features_after_corr) < min_features_to_keep && length(features_after_outlier) >= min_features_to_keep
    fprintf('\nWARNING: Only %d features passed correlation threshold (%.2f)\n', ...
        length(features_after_corr), correlation_threshold);
    fprintf('Keeping top %d features by correlation to ensure minimum feature count\n', min_features_to_keep);
    
    % Sort by max correlation and keep top N
    [sorted_corr, sorted_idx] = sort(max_corr, 'descend');
    top_n_idx = sorted_idx(1:min(min_features_to_keep, length(features_after_outlier)));
    features_after_corr = features_after_outlier(top_n_idx);
    
    fprintf('Selected top %d features:\n', length(features_after_corr));
    for i = 1:length(features_after_corr)
        fprintf('  %d. %s (corr=%.4f)\n', i, features_after_corr{i}, max_corr(top_n_idx(i)));
    end
end

%% Step 3-2: SOC별 상관계수 분석 및 표 출력
fprintf('\n=== Step 3-2: SOC별 상관계수 분석 (Capacity_C3) ===\n');

% SOC별 상관계수 계산
soc_levels = [90, 70, 50];
soc_labels = {'SOC90', 'SOC70', 'SOC50'};
event_types = {'Charge', 'Discharge'};

% 상관계수 저장용 구조체
soc_corr_table = struct();

% Capacity_C3와의 상관계수만 계산
if ismember('Capacity_C3', targetVars) && ismember('SOC', allTrainData.Properties.VariableNames)
    fprintf('\n--- SOC별 상관계수 표 (Capacity_C3) ---\n');
    fprintf('\n');
    
    % 표 헤더 출력
    fprintf('%-12s | %-10s | %10s | %10s | %10s | %10s\n', ...
        'SOC', 'EventType', 'R_1s', 'R_3s', 'R_5s', 'N (Samples)');
    fprintf('%s\n', repmat('-', 1, 70));
    
    for soc_idx = 1:length(soc_levels)
        soc_level = soc_levels(soc_idx);
        soc_label = soc_labels{soc_idx};
        
        for event_idx = 1:length(event_types)
            event_type = event_types{event_idx};
            
            % SOC 필터링
            if isnumeric(allTrainData.SOC)
                soc_mask = allTrainData.SOC == soc_level;
            elseif iscategorical(allTrainData.SOC) || iscellstr(allTrainData.SOC)
                soc_str = allTrainData.SOC;
                if iscategorical(soc_str)
                    soc_str = cellstr(soc_str);
                end
                soc_mask = false(height(allTrainData), 1);
                for i = 1:height(allTrainData)
                    soc_val = soc_str{i};
                    num_match = regexp(soc_val, '(\d+)', 'tokens');
                    if ~isempty(num_match)
                        soc_num = str2double(num_match{1}{1});
                        if soc_num == soc_level
                            soc_mask(i) = true;
                        end
                    end
                end
            else
                try
                    soc_numeric = str2double(allTrainData.SOC);
                    soc_mask = soc_numeric == soc_level;
                catch
                    continue;
                end
            end
            
            % EventType 필터링
            if ismember('EventType', allTrainData.Properties.VariableNames)
                event_mask = strcmp(allTrainData.EventType, event_type);
                combined_mask = soc_mask & event_mask;
            else
                combined_mask = soc_mask;
            end
            
            soc_event_data = allTrainData(combined_mask, :);
            
            if height(soc_event_data) < 10
                continue;
            end
            
            % 상관계수 계산
            X_soc = soc_event_data{:, features_after_corr};
            y_soc = soc_event_data.Capacity_C3;
            
            % 유효한 행 찾기
            valid_mask = ~any(isnan(X_soc), 2) & ~isnan(y_soc);
            X_soc_valid = X_soc(valid_mask, :);
            y_soc_valid = y_soc(valid_mask);
            
            if sum(valid_mask) < 10
                continue;
            end
            
            % 각 특성별 상관계수 계산 (R_1s, R_3s, R_5s)
            corr_r1s = NaN;
            corr_r3s = NaN;
            corr_r5s = NaN;
            
            % features_after_corr에서 R_1s, R_3s, R_5s 찾기
            r1s_idx = find(strcmp(features_after_corr, 'R_1s'), 1);
            r3s_idx = find(strcmp(features_after_corr, 'R_3s'), 1);
            r5s_idx = find(strcmp(features_after_corr, 'R_5s'), 1);
            
            % R_1s 상관계수 계산
            if ~isempty(r1s_idx) && r1s_idx <= size(X_soc_valid, 2)
                r1s_data = X_soc_valid(:, r1s_idx);
                valid_r1s = ~isnan(r1s_data) & ~isnan(y_soc_valid);
                if sum(valid_r1s) >= 10
                    corr_r1s = corr(r1s_data(valid_r1s), y_soc_valid(valid_r1s), 'rows', 'complete');
                end
            end
            
            % R_3s 상관계수 계산
            if ~isempty(r3s_idx) && r3s_idx <= size(X_soc_valid, 2)
                r3s_data = X_soc_valid(:, r3s_idx);
                valid_r3s = ~isnan(r3s_data) & ~isnan(y_soc_valid);
                if sum(valid_r3s) >= 10
                    corr_r3s = corr(r3s_data(valid_r3s), y_soc_valid(valid_r3s), 'rows', 'complete');
                end
            end
            
            % R_5s 상관계수 계산
            if ~isempty(r5s_idx) && r5s_idx <= size(X_soc_valid, 2)
                r5s_data = X_soc_valid(:, r5s_idx);
                valid_r5s = ~isnan(r5s_data) & ~isnan(y_soc_valid);
                if sum(valid_r5s) >= 10
                    corr_r5s = corr(r5s_data(valid_r5s), y_soc_valid(valid_r5s), 'rows', 'complete');
                end
            end
            
            % 표 출력
            fprintf('%-12s | %-10s | %10.4f | %10.4f | %10.4f | %10d\n', ...
                soc_label, event_type, corr_r1s, corr_r3s, corr_r5s, sum(valid_mask));
            
            % 저장
            key = sprintf('%s_%s', soc_label, event_type);
            soc_corr_table.(key).R_1s = corr_r1s;
            soc_corr_table.(key).R_3s = corr_r3s;
            soc_corr_table.(key).R_5s = corr_r5s;
            soc_corr_table.(key).N = sum(valid_mask);
        end
    end
    
    fprintf('\n');
end

%% Generate Correlation Visualizations
fprintf('\n=== Generating Correlation Visualizations ===\n');
figuresDir = fullfile(outputDir, 'figures', 'FeatureEng');
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

% 1. Top Features Correlation Bar Plot
if isfield(target_correlations, 'Capacity_C3')
    [sorted_corr, idx] = sort(abs(target_correlations.Capacity_C3), 'descend');
    sorted_feats = features_after_outlier(idx);
    
    % 상위 10개만 표시
    n_top = min(10, length(sorted_feats));
    top_corr = target_correlations.Capacity_C3(idx(1:n_top)); % 부호 포함
    top_feats = sorted_feats(1:n_top);
    
    fig = figure('Name', 'Top Feature Correlations', 'Position', [100, 100, 800, 600], 'Visible', 'on');
    b = barh(top_corr);
    b.FaceColor = 'flat';
    % 양의 상관관계는 파랑, 음의 상관관계(저항 등)는 빨강
    for k = 1:length(top_corr)
        if top_corr(k) < 0
            b.CData(k,:) = [1 0 0]; % Red for negative correlation
        else
            b.CData(k,:) = [0 0 1]; % Blue for positive correlation
        end
    end
    yticks(1:n_top);
    yticklabels(strrep(top_feats, '_', '\_'));
    xlabel('Pearson Correlation Coefficient (r)', 'FontSize', 12);
    ylabel('Features', 'FontSize', 12);
    title('Top Features Correlated with Capacity (C/3)', 'FontSize', 14);
    grid on;
    savefig(fullfile(figuresDir, 'Top_Correlations.fig'));
    fprintf('  Saved: Top_Correlations.fig\n');
    
    % 2. SOC별 상관분석 히트맵 (SOC 90, 50, 70) - 참고: DriveCycle_CorrelationAnalysis_04.m
    % SOC별로 데이터를 필터링하여 각각 상관분석 수행
    % X: 주행부하에서 추출한 충전/방전 시간별 저항값 (R_1s, R_3s, R_5s, R_10s, R_30s, R_60s)
    % Y: 용량 2개 (Capacity_C3, Capacity_OCV) + RPT DCIR (충전/방전 시간별 저항값)
    % NOTE: 히트맵 시각화에서는 원본 feature 리스트를 사용하여 모든 저항값(R_1s~R_60s)을 포함
    if length(featureVars) > 0 && ismember('Capacity_C3', allTrainData.Properties.VariableNames) && ismember('SOC', allTrainData.Properties.VariableNames)
        % === DEBUG: Feature Variables Analysis ===
        fprintf('\n=== [DEBUG] Feature Variables Analysis for Heatmap ===\n');
        fprintf('Original featureVars: %d개\n', length(featureVars));
        fprintf('  featureVars list: %s\n', strjoin(featureVars, ', '));
        
        % Check missing ratio for each feature
        fprintf('\n[DEBUG] Missing ratio for each feature:\n');
        for f_idx = 1:length(featureVars)
            feat = featureVars{f_idx};
            if ismember(feat, allTrainData.Properties.VariableNames)
                feat_data = allTrainData.(feat);
                missing_ratio = 100 * sum(isnan(feat_data)) / length(feat_data);
                valid_count = sum(~isnan(feat_data));
                fprintf('  %s: %.1f%% missing (%d valid / %d total)\n', ...
                    feat, missing_ratio, valid_count, length(feat_data));
            else
                fprintf('  %s: NOT FOUND in table\n', feat);
            end
        end
        
        % X: 주행부하 저항값들 (원본 feature 리스트 사용 - 모든 저항값 포함)
        x_vars = featureVars;
        
        % Y: 용량 2개 + RPT DCIR 값들
        y_vars = {'Capacity_C3', 'Capacity_OCV'};
        
        % RPT DCIR 변수들 찾기 (RPT_Chg_*, RPT_Dch_*)
        allVars_in_data = allTrainData.Properties.VariableNames;
        rpt_vars = allVars_in_data(startsWith(allVars_in_data, 'RPT_'));
        y_vars = [y_vars, rpt_vars];
        
        % 전체 변수 (X + Y)
        selected_vars = [x_vars, y_vars];
        
        % Check if all variables exist
        available_vars = intersect(selected_vars, allTrainData.Properties.VariableNames);
        
        % X와 Y 구분
        available_x_vars = intersect(x_vars, allTrainData.Properties.VariableNames);
        available_y_vars = intersect(y_vars, allTrainData.Properties.VariableNames);
        
        % === DEBUG: Available Variables ===
        fprintf('\n=== [DEBUG] Available Variables for Heatmap ===\n');
        fprintf('X variables requested: %d개\n', length(x_vars));
        fprintf('  Requested X vars: %s\n', strjoin(x_vars, ', '));
        fprintf('X variables available in table: %d개\n', length(available_x_vars));
        if length(available_x_vars) > 0
            fprintf('  Available X vars: %s\n', strjoin(available_x_vars, ', '));
        else
            fprintf('  WARNING: No X variables available in table!\n');
        end
        
        % Check which X variables are missing
        missing_x_vars = setdiff(x_vars, available_x_vars);
        if ~isempty(missing_x_vars)
            fprintf('  Missing X vars (not in table): %s\n', strjoin(missing_x_vars, ', '));
        end
        
        fprintf('\n=== SOC별 상관분석 히트맵 (X vs Y) ===\n');
        fprintf('X (주행부하 저항값): %d개 (원본 feature 리스트 사용 - 모든 저항값 포함)\n', length(available_x_vars));
        if length(available_x_vars) > 0
            fprintf('  X 변수 예시: %s\n', strjoin(available_x_vars(1:min(5, end)), ', '));
        end
        fprintf('Y (용량 + RPT DCIR): %d개\n', length(available_y_vars));
        if length(available_y_vars) > 0
            fprintf('  Y 변수 예시: %s\n', strjoin(available_y_vars(1:min(5, end)), ', '));
        end
        fprintf('히트맵에 포함될 변수 수: %d (X: %d + Y: %d)\n', length(available_vars), length(available_x_vars), length(available_y_vars));
        
        if length(available_x_vars) > 0 && length(available_y_vars) > 0
            % SOC별 + EventType별로 분석 (SOC 90, 50, 70) × (Charge, Discharge)
            soc_levels = [90, 70, 50];
            soc_labels = {'SOC90', 'SOC70', 'SOC50'};
            event_types = {'Charge', 'Discharge'};
            
            % === Initialize report data structure ===
            report_data = struct();
            report_data.feature_availability = struct();
            report_data.summary_table = [];
            
            fprintf('\n=== SOC별 + EventType별 상관분석 시각화 ===\n');
            
            for soc_idx = 1:length(soc_levels)
                soc_level = soc_levels(soc_idx);
                soc_label = soc_labels{soc_idx};
                
                for event_idx = 1:length(event_types)
                    event_type = event_types{event_idx};
                    
                    fprintf('\n--- %s - %s 분석 ---\n', soc_label, event_type);
                    
                    % SOC별 + EventType별 데이터 필터링
                    % SOC 필터링
                    if isnumeric(allTrainData.SOC)
                        soc_mask = allTrainData.SOC == soc_level;
                    elseif iscategorical(allTrainData.SOC) || iscellstr(allTrainData.SOC)
                        soc_str = allTrainData.SOC;
                        if iscategorical(soc_str)
                            soc_str = cellstr(soc_str);
                        end
                        soc_mask = false(height(allTrainData), 1);
                        for i = 1:height(allTrainData)
                            soc_val = soc_str{i};
                            num_match = regexp(soc_val, '(\d+)', 'tokens');
                            if ~isempty(num_match)
                                soc_num = str2double(num_match{1}{1});
                                if soc_num == soc_level
                                    soc_mask(i) = true;
                                end
                            end
                        end
                    else
                        try
                            soc_numeric = str2double(allTrainData.SOC);
                            soc_mask = soc_numeric == soc_level;
                        catch
                            fprintf('  WARNING: %s - SOC 컬럼 형식을 인식할 수 없습니다.\n', soc_label);
                            continue;
                        end
                    end
                    
                    % EventType 필터링
                    if ismember('EventType', allTrainData.Properties.VariableNames)
                        event_mask = strcmp(allTrainData.EventType, event_type);
                        combined_mask = soc_mask & event_mask;
                    else
                        fprintf('  WARNING: EventType 컬럼이 없습니다. SOC만 필터링합니다.\n');
                        combined_mask = soc_mask;
                    end
                    
                    soc_data = allTrainData(combined_mask, :);
                    
                    fprintf('  %s - %s 데이터: %d rows\n', soc_label, event_type, height(soc_data));
                    
                    % SOC 값 확인 (디버깅)
                    if height(soc_data) == 0
                        unique_soc = unique(allTrainData.SOC);
                        fprintf('  디버깅: 전체 데이터의 SOC 값들: %s\n', mat2str(unique_soc(1:min(10, length(unique_soc)))));
                        if ismember('EventType', allTrainData.Properties.VariableNames)
                            unique_events = unique(allTrainData.EventType);
                            fprintf('  디버깅: 전체 데이터의 EventType 값들: %s\n', strjoin(unique_events, ', '));
                        end
                    end
                
                if height(soc_data) < 3
                    fprintf('  WARNING: %s 데이터가 부족합니다 (n=%d). 건너뜁니다.\n', soc_label, height(soc_data));
                    continue;
                end
                
                    if height(soc_data) < 3
                        fprintf('  WARNING: %s - %s 데이터가 부족합니다 (n=%d). 건너뜁니다.\n', soc_label, event_type, height(soc_data));
                        continue;
                    end
                    
                    % Extract X and Y data separately
                    X_data = soc_data{:, available_x_vars};
                    Y_data = soc_data{:, available_y_vars};
                    
                    % === DEBUG: Detailed NaN analysis for each X variable ===
                    fprintf('\n  [DEBUG] NaN analysis for each X variable (%s - %s):\n', soc_label, event_type);
                    
                    % Store data for report
                    feature_stats = struct();
                    for x_idx = 1:length(available_x_vars)
                        x_var = available_x_vars{x_idx};
                        x_col_data = X_data(:, x_idx);
                        x_nan_count = sum(isnan(x_col_data));
                        x_valid_count = sum(~isnan(x_col_data));
                        x_nan_ratio = 100 * x_nan_count / length(x_col_data);
                        fprintf('    %s: %d valid / %d total (%.1f%% NaN)\n', ...
                            x_var, x_valid_count, length(x_col_data), x_nan_ratio);
                        
                        % Store for report
                        feature_stats.(x_var).valid_count = x_valid_count;
                        feature_stats.(x_var).total_count = length(x_col_data);
                        feature_stats.(x_var).nan_ratio = x_nan_ratio;
                        feature_stats.(x_var).available = (x_valid_count > 0);
                    end
                    
                    % Store in report_data
                    report_key = sprintf('%s_%s', soc_label, event_type);
                    if ~exist('report_data', 'var')
                        report_data = struct();
                        report_data.feature_availability = struct();
                        report_data.correlation_analysis = struct();
                    end
                    if ~isfield(report_data, 'feature_availability')
                        report_data.feature_availability = struct();
                    end
                    if ~isfield(report_data, 'correlation_analysis')
                        report_data.correlation_analysis = struct();
                    end
                    report_data.feature_availability.(report_key) = feature_stats;
                    
                    % Store correlation analysis results
                    if exist('corr_mat', 'var') && exist('r2_mat', 'var') && exist('p_mat', 'var') && ...
                       exist('selected_x_vars', 'var') && exist('selected_y_vars', 'var')
                        % Store correlation results for Capacity_C3 and Capacity_OCV
                        corr_results = struct();
                        corr_results.x_vars = selected_x_vars;
                        corr_results.y_vars = selected_y_vars;
                        corr_results.correlation_matrix = corr_mat;
                        corr_results.r2_matrix = r2_mat;
                        corr_results.pvalue_matrix = p_mat;
                        corr_results.n_samples_matrix = n_samples_mat;
                        report_data.correlation_analysis.(report_key) = corr_results;
                    end
                    
                    % Find valid rows (pairwise deletion: 각 X-Y 쌍별로 유효한 데이터만 사용)
                    % NaN 비율: 모든 값이 NaN인 행의 비율
                    x_all_nan = all(isnan(X_data), 2);
                    y_all_nan = all(isnan(Y_data), 2);
                    fprintf('  X 데이터: %d rows, 모든 X가 NaN인 행: %d (%.1f%%)\n', size(X_data, 1), sum(x_all_nan), 100*sum(x_all_nan)/size(X_data, 1));
                    fprintf('  Y 데이터: %d rows, 모든 Y가 NaN인 행: %d (%.1f%%)\n', size(Y_data, 1), sum(y_all_nan), 100*sum(y_all_nan)/size(Y_data, 1));
                    
                    % 최소한 하나의 X와 하나의 Y가 유효한 행 찾기
                    valid_rows_x = ~x_all_nan;  % 모든 X가 NaN이 아닌 행
                    valid_rows_y = ~y_all_nan;  % 모든 Y가 NaN이 아닌 행
                    valid_rows = valid_rows_x | valid_rows_y;  % 하나라도 유효하면 포함
                    
                    X_data_clean = X_data(valid_rows, :);
                    Y_data_clean = Y_data(valid_rows, :);
                    
                    if size(X_data_clean, 1) < 3
                        fprintf('  WARNING: %s - %s 유효한 데이터가 부족합니다. 건너뜁니다.\n', soc_label, event_type);
                        continue;
                    end
                    
                    % Calculate preliminary correlations to select top 20
                    % X와 Capacity_C3 간의 상관계수로 상위 20개 선택
                    if ismember('Capacity_C3', available_y_vars)
                        cap_idx = find(strcmp(available_y_vars, 'Capacity_C3'));
                        if ~isempty(cap_idx) && cap_idx <= size(Y_data_clean, 2)
                            y_cap = Y_data_clean(:, cap_idx);
                            x_corrs = zeros(size(X_data_clean, 2), 1);
                            
                            for x_idx = 1:size(X_data_clean, 2)
                                x_vec = X_data_clean(:, x_idx);
                                valid_pair = ~isnan(y_cap) & ~isnan(x_vec);
                                if sum(valid_pair) >= 3
                                    [R_pair, ~] = corrcoef(y_cap(valid_pair), x_vec(valid_pair));
                                    x_corrs(x_idx) = abs(R_pair(1, 2));
                                else
                                    x_corrs(x_idx) = 0;
                                end
                            end
                            
                            % 상위 20개 X 선택
                            [~, x_sort_idx] = sort(x_corrs, 'descend');
                            n_x_select = min(20, length(available_x_vars));
                            selected_x_indices = x_sort_idx(1:n_x_select);
                            selected_x_vars = available_x_vars(selected_x_indices);
                            
                            % Y는 용량 2개 + 상위 18개 RPT 변수
                            capacity_vars = {'Capacity_C3', 'Capacity_OCV'};
                            rpt_vars_list = setdiff(available_y_vars, capacity_vars);
                            
                            if length(rpt_vars_list) > 0
                                % RPT 변수들과 Capacity_C3 간의 상관계수 계산
                                rpt_corrs = zeros(length(rpt_vars_list), 1);
                                for rpt_idx = 1:length(rpt_vars_list)
                                    rpt_var = rpt_vars_list{rpt_idx};
                                    rpt_col_idx = find(strcmp(available_y_vars, rpt_var));
                                    if ~isempty(rpt_col_idx) && rpt_col_idx <= size(Y_data_clean, 2)
                                        rpt_vec = Y_data_clean(:, rpt_col_idx);
                                        valid_pair = ~isnan(y_cap) & ~isnan(rpt_vec);
                                        if sum(valid_pair) >= 3
                                            [R_pair, ~] = corrcoef(y_cap(valid_pair), rpt_vec(valid_pair));
                                            rpt_corrs(rpt_idx) = abs(R_pair(1, 2));
                                        end
                                    end
                                end
                                
                                [~, rpt_sort_idx] = sort(rpt_corrs, 'descend');
                                n_rpt_select = min(18, length(rpt_vars_list));
                                selected_rpt_indices = rpt_sort_idx(1:n_rpt_select);
                                selected_rpt_vars = rpt_vars_list(selected_rpt_indices);
                            else
                                selected_rpt_vars = {};
                            end
                            
                            selected_y_vars = [capacity_vars, selected_rpt_vars];
                            selected_y_vars = intersect(selected_y_vars, available_y_vars);
                            
                            fprintf('  선택된 X 변수: %d개 (상위 %d개)\n', length(selected_x_vars), n_x_select);
                            fprintf('    X 변수 목록: %s\n', strjoin(selected_x_vars, ', '));
                            fprintf('  선택된 Y 변수: %d개 (용량 2개 + RPT 상위 %d개)\n', length(selected_y_vars), length(selected_rpt_vars));
                            fprintf('    Y 변수 목록: %s\n', strjoin(selected_y_vars, ', '));
                        else
                            % Capacity_C3가 없으면 모든 변수 사용
                            selected_x_vars = available_x_vars(1:min(20, end));
                            selected_y_vars = available_y_vars(1:min(20, end));
                        end
                    else
                        % Capacity_C3가 없으면 모든 변수 사용 (최대 20개)
                        selected_x_vars = available_x_vars(1:min(20, end));
                        selected_y_vars = available_y_vars(1:min(20, end));
                    end
                    
                    % 선택된 변수로 데이터 추출
                    X_data_selected = soc_data{:, selected_x_vars};
                    Y_data_selected = soc_data{:, selected_y_vars};
                    
                    % 유효한 행 찾기 (최소한 하나의 X와 하나의 Y가 유효해야 함)
                    x_selected_all_nan = all(isnan(X_data_selected), 2);
                    y_selected_all_nan = all(isnan(Y_data_selected), 2);
                    valid_rows = ~x_selected_all_nan & ~y_selected_all_nan;  % X와 Y 모두 유효한 행만
                    
                    X_data_clean = X_data_selected(valid_rows, :);
                    Y_data_clean = Y_data_selected(valid_rows, :);
                    
                    fprintf('  유효한 데이터: %d rows (X와 Y 모두 유효한 행)\n', size(X_data_clean, 1));
                    fprintf('  최종 X 변수: %d개, Y 변수: %d개\n', size(X_data_clean, 2), size(Y_data_clean, 2));
                    
                    % 실제 유효한 데이터 포인트 수 확인
                    if size(X_data_clean, 1) > 0
                        x_valid_count = sum(sum(~isnan(X_data_clean)));
                        y_valid_count = sum(sum(~isnan(Y_data_clean)));
                        fprintf('  X 유효한 값 개수: %d / %d (%.1f%%)\n', x_valid_count, numel(X_data_clean), 100*x_valid_count/numel(X_data_clean));
                        fprintf('  Y 유효한 값 개수: %d / %d (%.1f%%)\n', y_valid_count, numel(Y_data_clean), 100*y_valid_count/numel(Y_data_clean));
                    end
                    
                    if size(X_data_clean, 1) < 3
                        fprintf('  WARNING: %s - %s 유효한 데이터가 부족합니다. 건너뜁니다.\n', soc_label, event_type);
                        continue;
                    end
                    
                    % Calculate X-Y correlation matrix (X vs Y only)
                    % corr_mat: rows = Y variables, columns = X variables
                    corr_mat = zeros(size(Y_data_clean, 2), size(X_data_clean, 2));
                    p_mat = zeros(size(Y_data_clean, 2), size(X_data_clean, 2));
                    r2_mat = zeros(size(Y_data_clean, 2), size(X_data_clean, 2));  % R² (결정계수)
                    n_samples_mat = zeros(size(Y_data_clean, 2), size(X_data_clean, 2));  % 유효 샘플 수
                    
                    for y_idx = 1:size(Y_data_clean, 2)
                        for x_idx = 1:size(X_data_clean, 2)
                            y_vec = Y_data_clean(:, y_idx);
                            x_vec = X_data_clean(:, x_idx);
                            
                            % Remove NaN pairs
                            valid_pair = ~isnan(y_vec) & ~isnan(x_vec);
                            n_valid = sum(valid_pair);
                            n_samples_mat(y_idx, x_idx) = n_valid;
                            
                            if n_valid >= 3
                                [R_pair, P_pair] = corrcoef(y_vec(valid_pair), x_vec(valid_pair));
                                corr_mat(y_idx, x_idx) = R_pair(1, 2);
                                p_mat(y_idx, x_idx) = P_pair(1, 2);
                                r2_mat(y_idx, x_idx) = R_pair(1, 2)^2;  % R² = r²
                            else
                                corr_mat(y_idx, x_idx) = NaN;
                                p_mat(y_idx, x_idx) = NaN;
                                r2_mat(y_idx, x_idx) = NaN;
                            end
                        end
                    end
                    
                    % === DEBUG: X-Y Correlation Analysis Report (for Documentation) ===
                    fprintf('\n  === [DEBUG] X-Y Correlation Analysis Report (%s - %s) ===\n', soc_label, event_type);
                    
                    % 주요 Y 변수 (Capacity_C3, Capacity_OCV)와 모든 X 변수 간의 상관분석 결과 출력
                    main_y_vars = {'Capacity_C3', 'Capacity_OCV'};
                    fprintf('\n  --- Correlation Analysis: X (Drive Cycle Resistance) vs Y (Capacity) ---\n');
                    fprintf('  X Variable | Y Variable    | Correlation (r) | R² (Coefficient) | p-value    | N (Samples)\n');
                    fprintf('  -----------|---------------|-----------------|-------------------|------------|------------\n');
                    
                    for main_y_idx = 1:length(main_y_vars)
                        main_y_var = main_y_vars{main_y_idx};
                        y_idx_in_selected = find(strcmp(selected_y_vars, main_y_var));
                        
                        if ~isempty(y_idx_in_selected) && y_idx_in_selected <= size(corr_mat, 1)
                            for x_idx = 1:size(X_data_clean, 2)
                                x_var = selected_x_vars{x_idx};
                                r_val = corr_mat(y_idx_in_selected, x_idx);
                                r2_val = r2_mat(y_idx_in_selected, x_idx);
                                p_val = p_mat(y_idx_in_selected, x_idx);
                                n_val = n_samples_mat(y_idx_in_selected, x_idx);
                                
                                if ~isnan(r_val)
                                    % p-value 표시 형식
                                    if p_val < 0.001
                                        p_str = sprintf('%.2e', p_val);
                                    elseif p_val < 0.01
                                        p_str = sprintf('%.3e', p_val);
                                    else
                                        p_str = sprintf('%.4f', p_val);
                                    end
                                    
                                    % 유의성 표시
                                    if p_val < 0.001
                                        sig_str = '***';
                                    elseif p_val < 0.01
                                        sig_str = '**';
                                    elseif p_val < 0.05
                                        sig_str = '*';
                                    else
                                        sig_str = '';
                                    end
                                    
                                    fprintf('  %-10s | %-13s | %14.4f | %16.4f | %10s | %11d\n', ...
                                        x_var, main_y_var, r_val, r2_val, [p_str, sig_str], n_val);
                                else
                                    fprintf('  %-10s | %-13s | %14s | %16s | %10s | %11s\n', ...
                                        x_var, main_y_var, 'NaN', 'NaN', 'NaN', '-');
                                end
                            end
                        end
                    end
                    
                    fprintf('\n  Significance: *** p<0.001, ** p<0.01, * p<0.05\n');
                    
                    % 상관계수 요약 (X 변수별로 Capacity_C3와의 상관계수 순위)
                    if ismember('Capacity_C3', selected_y_vars)
                        cap3_idx = find(strcmp(selected_y_vars, 'Capacity_C3'));
                        if ~isempty(cap3_idx) && cap3_idx <= size(corr_mat, 1)
                            fprintf('\n  --- Summary: X Variables Ranked by Correlation with Capacity_C3 ---\n');
                            x_corr_with_cap3 = corr_mat(cap3_idx, :);
                            x_r2_with_cap3 = r2_mat(cap3_idx, :);
                            x_p_with_cap3 = p_mat(cap3_idx, :);
                            
                            % 유효한 상관계수만 필터링
                            valid_corr_mask = ~isnan(x_corr_with_cap3);
                            valid_corrs = x_corr_with_cap3(valid_corr_mask);
                            valid_x_vars = selected_x_vars(valid_corr_mask);
                            valid_r2 = x_r2_with_cap3(valid_corr_mask);
                            valid_p = x_p_with_cap3(valid_corr_mask);
                            
                            % 절댓값 기준으로 정렬
                            [~, sort_idx] = sort(abs(valid_corrs), 'descend');
                            
                            fprintf('  Rank | X Variable | Correlation (r) | R²        | p-value    | Interpretation\n');
                            fprintf('  -----|------------|-----------------|-----------|------------|----------------\n');
                            for rank = 1:length(sort_idx)
                                idx = sort_idx(rank);
                                x_var = valid_x_vars{idx};
                                r_val = valid_corrs(idx);
                                r2_val = valid_r2(idx);
                                p_val = valid_p(idx);
                                
                                % 해석
                                if abs(r_val) >= 0.7
                                    interp_str = 'Strong';
                                elseif abs(r_val) >= 0.5
                                    interp_str = 'Moderate';
                                elseif abs(r_val) >= 0.3
                                    interp_str = 'Weak';
                                else
                                    interp_str = 'Very Weak';
                                end
                                
                                if r_val < 0
                                    interp_str = [interp_str, ' (Negative)'];
                                else
                                    interp_str = [interp_str, ' (Positive)'];
                                end
                                
                                if p_val < 0.05
                                    interp_str = [interp_str, ', Significant'];
                                else
                                    interp_str = [interp_str, ', Not Significant'];
                                end
                                
                                p_str = sprintf('%.4f', p_val);
                                if p_val < 0.001
                                    p_str = sprintf('%.2e', p_val);
                                end
                                
                                fprintf('  %4d | %-10s | %14.4f | %9.4f | %10s | %s\n', ...
                                    rank, x_var, r_val, r2_val, p_str, interp_str);
                            end
                        end
                    end
                    
                    % === DEBUG: Highlight all cells with |r| > 0.6 ===
                    fprintf('\n  --- All High Correlations (|r| > 0.6) in Heatmap ---\n');
                    high_corr_count = 0;
                    for y_idx = 1:size(corr_mat, 1)
                        for x_idx = 1:size(corr_mat, 2)
                            r_val = corr_mat(y_idx, x_idx);
                            if ~isnan(r_val) && abs(r_val) > 0.6
                                high_corr_count = high_corr_count + 1;
                                y_var = selected_y_vars{y_idx};
                                x_var = selected_x_vars{x_idx};
                                r2_val = r2_mat(y_idx, x_idx);
                                p_val = p_mat(y_idx, x_idx);
                                n_val = n_samples_mat(y_idx, x_idx);
                                
                                p_str = sprintf('%.2e', p_val);
                                if p_val < 0.001
                                    sig_str = '***';
                                elseif p_val < 0.01
                                    sig_str = '**';
                                elseif p_val < 0.05
                                    sig_str = '*';
                                else
                                    sig_str = '';
                                end
                                
                                fprintf('  [%d] %s ↔ %s: r=%.4f, R²=%.4f, p=%s%s, n=%d\n', ...
                                    high_corr_count, x_var, y_var, r_val, r2_val, p_str, sig_str, n_val);
                            end
                        end
                    end
                    if high_corr_count == 0
                        fprintf('  No correlations with |r| > 0.6 found.\n');
                    else
                        fprintf('  Total: %d high correlation pairs (|r| > 0.6)\n', high_corr_count);
                    end
                
                    % Create X-Y correlation heatmap (Y vs X)
                    fig2 = figure('Name', sprintf('Correlation Heatmap (X vs Y) - %s %s', soc_label, event_type), ...
                        'Position', [150 + soc_idx*50 + event_idx*30, 150 + soc_idx*50 + event_idx*30, max(1200, size(X_data_clean, 2)*50), max(800, size(Y_data_clean, 2)*50)], 'Visible', 'on');
                
                % Create cell labels with r, p-value
                cellLabels = cell(size(corr_mat));
                for i = 1:size(corr_mat, 1)  % Y variables (rows)
                    for j = 1:size(corr_mat, 2)  % X variables (columns)
                        r = corr_mat(i, j);
                        p = p_mat(i, j);
                        
                        if isnan(r)
                            cellLabels{i, j} = 'NaN';
                        else
                            if p < 0.001
                                pStr = sprintf('p=%.2e', p);
                            elseif p < 0.01
                                pStr = sprintf('p=%.3e', p);
                            else
                                pStr = sprintf('p=%.3f', p);
                            end
                            cellLabels{i, j} = sprintf('r=%.2f\n%s', r, pStr);
                        end
                    end
                end
                
                % Create heatmap using imagesc (참고: DriveCycle_CorrelationAnalysis_04.m)
                ax = axes('Position', [0.15, 0.1, 0.7, 0.75]);
                imagesc(corr_mat);
                
                % Custom red-blue colormap (참고: DriveCycle_CorrelationAnalysis_04.m)
                nColors = 256;
                blueVals = linspace(1, 0, nColors/2)';
                greenVals = linspace(0, 1, nColors/2)';
                redVals = linspace(0, 1, nColors/2)';
                blueToWhite = [redVals, greenVals, blueVals];
                
                redVals2 = ones(nColors/2, 1);
                greenVals2 = linspace(1, 0, nColors/2)';
                blueVals2 = linspace(1, 0, nColors/2)';
                whiteToRed = [redVals2, greenVals2, blueVals2];
                
                redBlueMap = [blueToWhite; whiteToRed];
                colormap(ax, redBlueMap);
                caxis([-1 1]);
                colorbar;
                
                    % Set axis labels
                    xLabels = strrep(selected_x_vars, '_', '\_');
                    yLabels = strrep(selected_y_vars, '_', '\_');
                
                set(ax, 'XTick', 1:length(selected_x_vars), 'XTickLabel', xLabels, ...
                    'YTick', 1:length(selected_y_vars), 'YTickLabel', yLabels, ...
                    'XTickLabelRotation', 45, 'FontSize', 9);
                
                xlabel('X: 주행부하 저항값', 'FontSize', 12, 'FontWeight', 'bold');
                ylabel('Y: 용량 + RPT DCIR', 'FontSize', 12, 'FontWeight', 'bold');
                
                % Add cell labels and highlight high correlations (|r| > 0.6)
                high_corr_threshold = 0.6;
                for i = 1:size(corr_mat, 1)
                    for j = 1:size(corr_mat, 2)
                        if ~isnan(corr_mat(i, j))
                            % Add text label
                            text(j, i, cellLabels{i, j}, ...
                                'HorizontalAlignment', 'center', ...
                                'VerticalAlignment', 'middle', ...
                                'FontSize', 8, ...
                                'Color', 'white', ...
                                'FontWeight', 'bold');
                            
                            % Highlight high correlations (|r| > 0.6) with thick border
                            if abs(corr_mat(i, j)) > high_corr_threshold
                                % Draw thick border around the cell
                                % Rectangle coordinates: [x_left, y_bottom, width, height]
                                % Note: imagesc uses pixel coordinates, so we need to adjust
                                rectangle('Position', [j-0.5, i-0.5, 1, 1], ...
                                    'EdgeColor', 'yellow', ...
                                    'LineWidth', 3, ...
                                    'LineStyle', '-');
                            end
                        end
                    end
                end
                
                    title(sprintf('Correlation Matrix (X vs Y) - %s %s\nX: 주행부하 저항값 (SOC=%d%%에서 추출한 모든 %s 이벤트의 저항값), Y: 용량+RPT DCIR (n=%d)', ...
                        soc_label, event_type, soc_level, event_type, size(X_data_clean, 1)), ...
                        'FontSize', 14, 'FontWeight', 'bold');
                    axis square;
                    
                    savePath_heatmap = fullfile(figuresDir, sprintf('Corr_Heatmap_%s_%s.fig', soc_label, event_type));
                    savefig(savePath_heatmap);
                    fprintf('  Saved: Corr_Heatmap_%s_%s.fig\n', soc_label, event_type);
                    close(fig2);  % Close figure after saving
                end
            end
        else
            fprintf('WARNING: X 또는 Y 변수가 없어서 SOC별 히트맵을 생성할 수 없습니다.\n');
            fprintf('  available_x_vars: %d개, available_y_vars: %d개\n', length(available_x_vars), length(available_y_vars));
        end
    else
        fprintf('WARNING: SOC별 히트맵 생성을 위한 조건이 충족되지 않았습니다.\n');
        fprintf('  featureVars: %d개\n', length(featureVars));
        fprintf('  Capacity_C3 존재: %d\n', ismember('Capacity_C3', allTrainData.Properties.VariableNames));
        fprintf('  SOC 존재: %d\n', ismember('SOC', allTrainData.Properties.VariableNames));
    end
    
end

fprintf('Feature engineering visualizations saved to: %s\n', figuresDir);

%% Step 4: Multicollinearity Analysis (VIF) - Iterative Removal
fprintf('\n=== Step 4: Multicollinearity Analysis (VIF) ===\n');

% Force check use_vif_filtering variable (in case it was modified elsewhere)
if ~exist('use_vif_filtering', 'var') || isempty(use_vif_filtering)
    use_vif_filtering = false;  % Default to false
end

if use_vif_filtering
    % Prepare data for VIF calculation
    X_vif = allTrainData{:, features_after_corr};
    valid_rows_vif = ~any(isnan(X_vif), 2);
    X_vif_clean = X_vif(valid_rows_vif, :);
    
    features_after_vif = features_after_corr;
    vif_values = [];
    
    if size(X_vif_clean, 1) > 0 && size(X_vif_clean, 2) > 0
        % Check if we have enough samples (need at least as many samples as features)
        n_samples = size(X_vif_clean, 1);
        n_features = size(X_vif_clean, 2);
        
        fprintf('Samples: %d, Features: %d\n', n_samples, n_features);
        
        if n_samples < n_features
            fprintf('WARNING: More features than samples. Using correlation-based filtering instead of VIF.\n');
            % Use correlation matrix instead
            corr_matrix = corr(X_vif_clean, 'rows', 'complete');
            
            % Find highly correlated pairs and remove one from each pair
            [row_idx, col_idx] = find(triu(corr_matrix, 1) > multicollinearity_threshold);
            
            features_to_remove = [];
            for p = 1:length(row_idx)
                % Keep the feature with higher index (arbitrary choice)
                features_to_remove = [features_to_remove, col_idx(p)];
            end
            features_to_remove = unique(features_to_remove);
            
            features_after_vif(features_to_remove) = [];
            fprintf('Removed %d features due to high correlation (correlation > %.2f)\n', ...
                length(features_to_remove), multicollinearity_threshold);
            fprintf('Remaining features: %d\n', length(features_after_vif));
        else
            % Iterative VIF removal: remove one feature at a time
            fprintf('Using iterative VIF removal...\n');
            max_iterations = 50;  % Prevent infinite loop
            iteration = 0;
            
            current_features = features_after_corr;
            current_X = X_vif_clean;
            
            while iteration < max_iterations
                iteration = iteration + 1;
                
                % Calculate VIF
                try
                    current_vif = calculateVIF(current_X);
                catch ME
                    fprintf('VIF calculation failed: %s\n', ME.message);
                    break;
                end
                
                % Find features with VIF > threshold
                high_vif_idx = find(current_vif > vif_threshold);
                
                if isempty(high_vif_idx)
                    % No more high VIF features
                    fprintf('No features with VIF > %d. Stopping.\n', vif_threshold);
                    break;
                end
                
                % Remove the feature with highest VIF
                [~, max_vif_idx] = max(current_vif);
                feature_to_remove = current_features(max_vif_idx);
                
                fprintf('Iteration %d: Removing %s (VIF = %.2f)\n', ...
                    iteration, feature_to_remove{1}, current_vif(max_vif_idx));
                
                % Remove from current lists
                current_features(max_vif_idx) = [];
                current_X(:, max_vif_idx) = [];
                
                % Check if we have enough samples left
                if size(current_X, 1) < size(current_X, 2)
                    fprintf('Not enough samples for remaining features. Stopping.\n');
                    break;
                end
                
                % Ensure minimum feature count
                if length(current_features) < min_features_to_keep
                    fprintf('WARNING: VIF removal would reduce features below minimum (%d). Stopping.\n', min_features_to_keep);
                    break;
                end
            end
            
            features_after_vif = current_features;
            vif_values = current_vif;
            
            fprintf('After iterative VIF removal: %d features remaining\n', length(features_after_vif));
        end
    else
        fprintf('WARNING: Not enough valid data for VIF calculation\n');
        features_after_vif = features_after_corr;
        vif_values = NaN(size(features_after_corr));
    end
else
    % Skip VIF filtering completely - use correlation-based multicollinearity check only
    fprintf('VIF filtering is DISABLED (use_vif_filtering = false)\n');
    fprintf('Skipping Step 4. Will use correlation-based multicollinearity check in Step 5 instead.\n');
    features_after_vif = features_after_corr;
    vif_values = [];
    fprintf('Features after Step 4 (VIF skipped): %d\n', length(features_after_vif));
end

%% Step 5: Feature Correlation Matrix (Additional multicollinearity check)
fprintf('\n=== Step 5: Feature Correlation Matrix Analysis ===\n');

X_corr = allTrainData{:, features_after_vif};
valid_rows_corr = ~any(isnan(X_corr), 2);
X_corr_clean = X_corr(valid_rows_corr, :);

if size(X_corr_clean, 1) > 10 && size(X_corr_clean, 2) > 1
    corr_matrix = corr(X_corr_clean, 'rows', 'complete');
    
    % Display feature correlation matrix
    fprintf('\n--- Feature-to-Feature Correlation Matrix ---\n');
    fprintf('Features: %s\n', strjoin(features_after_vif, ', '));
    fprintf('\nCorrelation Matrix (upper triangle):\n');
    fprintf('%12s', '');
    for i = 1:length(features_after_vif)
        fprintf('%12s', features_after_vif{i});
    end
    fprintf('\n');
    for i = 1:size(corr_matrix, 1)
        fprintf('%12s', features_after_vif{i});
        for j = 1:size(corr_matrix, 2)
            if j <= i
                fprintf('%12s', '-');
            else
                fprintf('%12.4f', corr_matrix(i, j));
            end
        end
        fprintf('\n');
    end
    
    % Find highly correlated feature pairs
    [row_idx, col_idx] = find(triu(corr_matrix, 1) > multicollinearity_threshold);
    
    % Also report correlations > 0.8 (moderate-high) for information
    [moderate_row, moderate_col] = find(triu(corr_matrix, 1) > 0.8);
    if ~isempty(moderate_row)
        fprintf('\n--- Moderate-High Correlations (|r| > 0.8) ---\n');
        for p = 1:length(moderate_row)
            fprintf('  %s ↔ %s: r = %.4f\n', ...
                features_after_vif{moderate_row(p)}, ...
                features_after_vif{moderate_col(p)}, ...
                corr_matrix(moderate_row(p), moderate_col(p)));
        end
    end
    
    if ~isempty(row_idx)
        fprintf('Found %d highly correlated feature pairs (correlation > %.2f)\n', ...
            length(row_idx), multicollinearity_threshold);
        
        % Remove one feature from each pair (keep the one with higher target correlation)
        features_to_remove = [];
        
        % Calculate correlation with targets for features_after_vif
        % Find indices of features_after_vif in the original features_after_outlier list
        [~, feat_indices] = ismember(features_after_vif, features_after_outlier);
        valid_indices = feat_indices > 0;
        feature_corrs = NaN(length(features_after_vif), 1);
        feature_corrs(valid_indices) = max_corr(feat_indices(valid_indices));
        
        for p = 1:length(row_idx)
            feat1_idx = row_idx(p);
            feat2_idx = col_idx(p);
            
            % Get correlations for these features
            corr1 = feature_corrs(feat1_idx);
            corr2 = feature_corrs(feat2_idx);
            
            % If correlations are available, keep the one with higher correlation
            if ~isnan(corr1) && ~isnan(corr2)
                if corr1 < corr2
                    features_to_remove = [features_to_remove, feat1_idx];
                else
                    features_to_remove = [features_to_remove, feat2_idx];
                end
            else
                % If correlation not available, remove the second one (arbitrary)
                features_to_remove = [features_to_remove, feat2_idx];
            end
        end
        
        features_to_remove = unique(features_to_remove);
        features_after_multicoll = features_after_vif;
        features_after_multicoll(features_to_remove) = [];
        
        fprintf('Removed %d features due to high correlation\n', length(features_to_remove));
        fprintf('Remaining features: %d\n', length(features_after_multicoll));
    else
        features_after_multicoll = features_after_vif;
        fprintf('No highly correlated pairs found\n');
    end
else
    features_after_multicoll = features_after_vif;
    fprintf('Not enough data for correlation matrix analysis\n');
end

%% Step 6: XGBoost Feature Importance (Final Selection)
fprintf('\n=== Step 6: XGBoost Feature Importance ===\n');
fprintf('This step will be performed during model training\n');
fprintf('For now, keeping %d features after previous steps\n', length(features_after_multicoll));

% Final feature list
selected_features = features_after_multicoll;

% IMPORTANT: Protect SOC variable (if exists) - Required for unified model
% SOC information is critical for distinguishing between aging and SOC effects
if ismember('SOC', allVars) && ~ismember('SOC', selected_features)
    fprintf('\nWARNING: SOC variable was removed during feature selection!\n');
    fprintf('Adding SOC back to selected features (required for unified model)\n');
    selected_features{end+1} = 'SOC';
end

%% Save results
fprintf('\n=== Saving Results ===\n');

% Save selected features
savePath_features = fullfile(outputDir, 'Selected_Features.mat');
save(savePath_features, 'selected_features', 'featureVars', 'targetVars', ...
    'missing_threshold', 'correlation_threshold', 'vif_threshold', ...
    'target_feature_count');

% Save detailed results
savePath_results = fullfile(outputDir, 'Feature_Engineering_Results.mat');
save(savePath_results, 'selected_features', 'features_after_missing', ...
    'features_after_outlier', 'features_after_corr', 'features_after_vif', 'features_after_multicoll', ...
    'missing_ratio', 'target_correlations', 'vif_values', ...
    'remove_missing', 'keep_corr', 'allTrainData', 'outlier_stats');

fprintf('Saved selected features: %s\n', savePath_features);
fprintf('Saved detailed results: %s\n', savePath_results);

fprintf('\n=== Feature Engineering Summary ===\n');
fprintf('Initial features: %d\n', length(featureVars));
fprintf('After missing value removal: %d\n', length(features_after_missing));
if exist('outlier_removed_count', 'var') && outlier_removed_count > 0
    fprintf('After outlier removal: %d features, %d rows removed (%.2f%%)\n', ...
        length(features_after_outlier), outlier_removed_count, ...
        100*outlier_removed_count/(outlier_removed_count + height(allTrainData)));
end
fprintf('After correlation filtering: %d\n', length(features_after_corr));
fprintf('After VIF filtering: %d\n', length(features_after_vif));
fprintf('After multicollinearity check: %d\n', length(features_after_multicoll));
fprintf('Final selected features: %d\n', length(selected_features));

% === Generate Feature Availability Report (Debug Output for Documentation) ===
if exist('report_data', 'var') && isfield(report_data, 'feature_availability') && ...
   ~isempty(fieldnames(report_data.feature_availability))
    fprintf('\n=== Feature Availability Report (for Documentation) ===\n');
    
    soc_levels = [90, 70, 50];
    soc_labels = {'SOC90', 'SOC70', 'SOC50'};
    event_types = {'Charge', 'Discharge'};
    
    % Get available_x_vars from report_data or use featureVars
    if exist('available_x_vars', 'var') && ~isempty(available_x_vars)
        report_features = available_x_vars;
    else
        report_features = featureVars;
    end
    
    % Display summary grouped by feature
    fprintf('\n--- Feature Availability Summary (for Report Documentation) ---\n');
    fprintf('Total features analyzed: %d\n', length(report_features));
    fprintf('Total combinations: %d (3 SOC × 2 EventType)\n', length(soc_levels) * length(event_types));
    
    for f_idx = 1:length(report_features)
        feat = report_features{f_idx};
        fprintf('\n[Feature: %s]\n', feat);
        fprintf('  SOC Level | Event Type  | Status      | Valid Count | Total Count | NaN Ratio\n');
        fprintf('  ----------|-------------|-------------|-------------|-------------|----------\n');
        
        available_count = 0;
        for soc_idx = 1:length(soc_levels)
            soc_label = soc_labels{soc_idx};
            for event_idx = 1:length(event_types)
                event_type = event_types{event_idx};
                report_key = sprintf('%s_%s', soc_label, event_type);
                
                if isfield(report_data.feature_availability, report_key)
                    feature_stats = report_data.feature_availability.(report_key);
                    if isfield(feature_stats, feat)
                        stats = feature_stats.(feat);
                        if stats.available
                            available_count = available_count + 1;
                            status_str = 'Available  ';
                        else
                            status_str = 'NOT Avail.';
                        end
                        fprintf('  %-9s | %-11s | %-11s | %11d | %11d | %8.1f%%\n', ...
                            soc_label, event_type, status_str, ...
                            stats.valid_count, stats.total_count, stats.nan_ratio);
                    else
                        fprintf('  %-9s | %-11s | Not Found   |            - |            - |        -\n', ...
                            soc_label, event_type);
                    end
                else
                    fprintf('  %-9s | %-11s | No Data     |            - |            - |        -\n', ...
                        soc_label, event_type);
                end
            end
        end
        fprintf('  Summary: Available in %d/6 combinations\n', available_count);
    end
    
    fprintf('\n--- Interpretation for Report ---\n');
    fprintf('This table shows which resistance features (R_1s, R_3s, R_5s, R_10s, R_30s, R_60s)\n');
    fprintf('are available for each SOC level (90%%, 70%%, 50%%) and event type (Charge, Discharge).\n');
    fprintf('Features with high NaN ratios (>70%%) are typically excluded from feature selection.\n');
end

fprintf('\n=== Feature Engineering Complete ===\n');

%% Helper function: Calculate VIF
function vif = calculateVIF(X)
    % X: n x p matrix (n samples, p features)
    % Returns: p x 1 vector of VIF values
    % Improved version with better handling of rank deficiency
    
    [n, p] = size(X);
    vif = zeros(p, 1);
    
    % Normalize features to avoid numerical issues
    X_norm = (X - mean(X, 1)) ./ (std(X, 0, 1) + eps);
    
    for i = 1:p
        % Regress feature i on all other features
        y = X_norm(:, i);
        X_other = X_norm;
        X_other(:, i) = [];  % Remove feature i
        
        % Use pinv (pseudo-inverse) to handle rank deficiency
        try
            % Add intercept (already normalized, so mean is 0)
            % Use least squares with regularization for stability
            beta = pinv(X_other) * y;
            y_pred = X_other * beta;
            
            % Calculate R-squared
            ss_res = sum((y - y_pred).^2);
            ss_tot = sum((y - mean(y)).^2);
            
            if ss_tot > eps
                r_squared = 1 - (ss_res / ss_tot);
                r_squared = max(0, min(1, r_squared));  % Clamp to [0, 1]
            else
                r_squared = 0;
            end
            
            % VIF = 1 / (1 - R^2)
            if r_squared < 1 - eps
                vif(i) = 1 / (1 - r_squared);
            else
                vif(i) = Inf;  % Perfect multicollinearity
            end
        catch ME
            % If calculation fails, set to high value
            vif(i) = Inf;
        end
    end
end

