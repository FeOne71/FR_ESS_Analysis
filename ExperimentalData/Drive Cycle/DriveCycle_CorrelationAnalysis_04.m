%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 04_DriveCycle_CorrelationAnalysis.m
% 상관분석 스크립트
% 
% 목적: 
% - DriveCycle_Summary_Table.mat를 불러와서 상관분석 수행
% - 6개 저항성분(R_1s, R_3s, R_5s, R_10s, R_30s, R_60s)과 용량(Capacity_C3) 간의 상관관계 분석
% - 전체 통합 분석 (DC Profile 구분 없이)
% - Pearson 상관계수 계산
%
% 입력:
% - DriveCycle_Summary_Table.mat (02번 스크립트 출력)
%
% 출력:
% - 상관계수 결과 테이블 (콘솔 출력)
% - 상관분석 결과 저장 (선택적)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('=== Drive Cycle Correlation Analysis ===\n');

%% Configuration - Paths
inputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
outputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Configuration - Analysis Settings
% Event Type Selection
targetEventType = 'All';  % 'All': both charge and discharge, 'Charge': only charge, 'Discharge': only discharge

%% Load Summary Table
summaryTablePath = fullfile(inputDir, 'DriveCycle_Summary_Table.mat');
if ~exist(summaryTablePath, 'file')
    fprintf('ERROR: DriveCycle_Summary_Table.mat not found!\n');
    fprintf('Expected path: %s\n', summaryTablePath);
    fprintf('Please run 02_DriveCycle_DataAggregation.m first to generate the summary table.\n');
    return;
end

load(summaryTablePath, 'summaryTable');
fprintf('Loaded summary table: %d rows\n', height(summaryTable));

% Debug: Check table structure
fprintf('\n=== Table Structure Check ===\n');
fprintf('Table columns: %s\n', strjoin(summaryTable.Properties.VariableNames, ', '));
if ismember('EventType', summaryTable.Properties.VariableNames)
    fprintf('Unique Event Types: %s\n', strjoin(unique(summaryTable.EventType), ', '));
else
    fprintf('WARNING: EventType column not found in table!\n');
end
fprintf('Unique DC Profiles: %s\n', strjoin(unique(summaryTable.DC_Profile), ', '));
fprintf('Unique Channels: %s\n', mat2str(unique(summaryTable.Channel)'));
fprintf('Unique Cycles: %s\n', mat2str(unique(summaryTable.Cycle)'));
fprintf('Data range - Capacity_C3: [%.2f, %.2f] Ah\n', min(summaryTable.Capacity_C3), max(summaryTable.Capacity_C3));

% Data composition analysis
fprintf('\n=== Data Composition Analysis ===\n');
fprintf('Total rows in summary table: %d\n', height(summaryTable));
fprintf('Unique Cycle × Channel combinations: %d\n', ...
    height(unique(summaryTable(:, {'Cycle', 'Channel'}))));
fprintf('Unique Cycle × Channel × DC_Profile combinations: %d\n', ...
    height(unique(summaryTable(:, {'Cycle', 'Channel', 'DC_Profile'}))));

% Count by Cycle × Channel (RPT × Cell combinations)
fprintf('\nCycle × Channel combinations (RPT × Cell):\n');
cycleChannelCombos = unique(summaryTable(:, {'Cycle', 'Channel'}));
fprintf('  Total unique combinations: %d\n', height(cycleChannelCombos));
for i = 1:min(10, height(cycleChannelCombos))
    fprintf('  Cycle %d × Channel %d\n', cycleChannelCombos.Cycle(i), cycleChannelCombos.Channel(i));
end
if height(cycleChannelCombos) > 10
    fprintf('  ... (showing first 10)\n');
end

% Count by Cycle
fprintf('\nData points by Cycle (RPT):\n');
uniqueCycles = unique(summaryTable.Cycle);
for i = 1:length(uniqueCycles)
    cycleNum = uniqueCycles(i);
    count = sum(summaryTable.Cycle == cycleNum);
    fprintf('  Cycle %d: %d rows\n', cycleNum, count);
end

% Count by Channel
fprintf('\nData points by Channel (Cell):\n');
uniqueChannels = unique(summaryTable.Channel);
for i = 1:length(uniqueChannels)
    chNum = uniqueChannels(i);
    count = sum(summaryTable.Channel == chNum);
    fprintf('  Channel %d: %d rows\n', chNum, count);
end

% Filter by EventType if specified
if ismember('EventType', summaryTable.Properties.VariableNames)
    if strcmp(targetEventType, 'Charge')
        summaryTable = summaryTable(strcmp(summaryTable.EventType, 'Charge'), :);
        fprintf('\nFiltered to Charge events only: %d rows\n', height(summaryTable));
    elseif strcmp(targetEventType, 'Discharge')
        summaryTable = summaryTable(strcmp(summaryTable.EventType, 'Discharge'), :);
        fprintf('\nFiltered to Discharge events only: %d rows\n', height(summaryTable));
    else
        fprintf('\nUsing all event types: %d rows\n', height(summaryTable));
    end
end

%% Outlier Removal using IQR method (Cycle-wise)
fprintf('\n\n========================================\n');
fprintf('=== Outlier Removal (IQR Method, Cycle-wise) ===\n');
fprintf('========================================\n');
fprintf('Method: IQR (Interquartile Range) with 1.5 multiplier\n');
fprintf('Applied: Separately for each Cycle\n');
fprintf('Variables: All variables (Capacity_C3, R_1s, R_3s, R_5s, R_10s, R_30s, R_60s)\n');
fprintf('Criterion: Values outside [Q1 - 1.5*IQR, Q3 + 1.5*IQR] are considered outliers\n\n');

% Store original data count
originalRowCount = height(summaryTable);
fprintf('Original data rows: %d\n', originalRowCount);

% Define variables for outlier removal
variableNames = {'Capacity_C3', 'R_1s', 'R_3s', 'R_5s', 'R_10s', 'R_30s', 'R_60s'};
variableLabels = {'Capacity_C3', 'Rchg 1s', 'Rchg 3s', 'Rchg 5s', 'Rchg 10s', 'Rchg 30s', 'Rchg 60s'};

% Check which variables exist
availableVarsForOutlier = {};
availableLabelsForOutlier = {};
for idx = 1:length(variableNames)
    if ismember(variableNames{idx}, summaryTable.Properties.VariableNames)
        availableVarsForOutlier{end+1} = variableNames{idx};
        availableLabelsForOutlier{end+1} = variableLabels{idx};
    end
end

fprintf('Variables to check for outliers: %s\n', strjoin(availableLabelsForOutlier, ', '));
fprintf('\n');

% Initialize outlier mask (true = keep, false = remove)
outlierMask = true(height(summaryTable), 1);

% Get unique cycles
uniqueCycles = unique(summaryTable.Cycle);
fprintf('Processing %d cycles for outlier removal...\n', length(uniqueCycles));

% Statistics for debugging
outlierStats = struct();
outlierStats.ByCycle = struct();
outlierStats.ByVariable = struct();
totalOutliersRemoved = 0;

% Process each cycle separately
for cycleIdx = 1:length(uniqueCycles)
    cycleNum = uniqueCycles(cycleIdx);
    cycleMask = summaryTable.Cycle == cycleNum;
    cycleData = summaryTable(cycleMask, :);
    cycleRowIndices = find(cycleMask);
    
    fprintf('\n--- Cycle %d (n=%d) ---\n', cycleNum, height(cycleData));
    
    cycleOutlierCount = 0;
    cycleOutlierDetails = {};
    
    % Check each variable for outliers in this cycle
    for varIdx = 1:length(availableVarsForOutlier)
        varName = availableVarsForOutlier{varIdx};
        varLabel = availableLabelsForOutlier{varIdx};
        
        % Extract valid data for this variable in this cycle
        varData = cycleData.(varName);
        validMask = ~isnan(varData);
        validData = varData(validMask);
        
        if length(validData) < 4  % Need at least 4 points for IQR
            continue;
        end
        
        % Calculate IQR
        Q1 = prctile(validData, 25);
        Q3 = prctile(validData, 75);
        IQR_val = Q3 - Q1;
        lowerBound = Q1 - 1.5 * IQR_val;
        upperBound = Q3 + 1.5 * IQR_val;
        
        % Find outliers
        outlierMask_var = (varData < lowerBound) | (varData > upperBound);
        outlierCount = sum(outlierMask_var);
        
        if outlierCount > 0
            cycleOutlierCount = cycleOutlierCount + outlierCount;
            outlierValues = varData(outlierMask_var);
            
            % Mark outliers in the main table
            for i = 1:length(cycleRowIndices)
                if outlierMask_var(i)
                    outlierMask(cycleRowIndices(i)) = false;
                end
            end
            
            fprintf('  %s: %d outliers (Q1=%.4f, Q3=%.4f, IQR=%.4f, Range=[%.4f, %.4f])\n', ...
                varLabel, outlierCount, Q1, Q3, IQR_val, lowerBound, upperBound);
            fprintf('    Outlier values: %s\n', mat2str(outlierValues', 4));
            
            % Store details
            cycleOutlierDetails{end+1} = sprintf('%s: %d outliers', varLabel, outlierCount);
        end
    end
    
    if cycleOutlierCount > 0
        fprintf('  Total outliers in Cycle %d: %d rows\n', cycleNum, cycleOutlierCount);
        totalOutliersRemoved = totalOutliersRemoved + cycleOutlierCount;
    else
        fprintf('  No outliers detected in Cycle %d\n', cycleNum);
    end
    
    % Store statistics
    outlierStats.ByCycle.(sprintf('Cycle%d', cycleNum)) = cycleOutlierCount;
end

% Remove outliers from summary table
summaryTable_original = summaryTable;  % Keep original for comparison
summaryTable = summaryTable(outlierMask, :);

fprintf('\n========================================\n');
fprintf('Outlier Removal Summary:\n');
fprintf('  Original rows: %d\n', originalRowCount);
fprintf('  Rows after removal: %d\n', height(summaryTable));
fprintf('  Rows removed: %d (%.2f%%)\n', ...
    originalRowCount - height(summaryTable), ...
    100 * (originalRowCount - height(summaryTable)) / originalRowCount);
fprintf('========================================\n');

%% Define variables for correlation matrix
% All variables: Capacity_C3, R_1s, R_3s, R_5s, R_10s, R_30s, R_60s
% Note: Using outlier-removed data
variableNames = {'Capacity_C3', 'R_1s', 'R_3s', 'R_5s', 'R_10s', 'R_30s', 'R_60s'};
variableLabels = {'Capacity_C3', 'Rchg 1s', 'Rchg 3s', 'Rchg 5s', 'Rchg 10s', 'Rchg 30s', 'Rchg 60s'};

%% Correlation Matrix Analysis
fprintf('\n\n========================================\n');
fprintf('=== Correlation Matrix Analysis ===\n');
fprintf('========================================\n');
fprintf('Analysis Type: Overall (전체 데이터셋 통합 분석)\n');
fprintf('Event Type: %s\n', targetEventType);
fprintf('\n');
fprintf('[분석 방식 설명]\n');
fprintf('- 모든 셀(Channel)의 데이터를 통합하여 하나의 큰 데이터셋으로 분석\n');
fprintf('- Cycle, Channel, DC_Profile 구분 없이 모든 행을 하나로 합쳐서 상관분석 수행\n');
fprintf('- 각 행 = 하나의 데이터 포인트 (Cycle × Channel × DC_Profile × EventType 조합)\n');
fprintf('- 목표: 셀 조건 상관없이 전체 데이터셋에서 용량과 저항의 상관관계 파악\n');
fprintf('- 이상치 제거: IQR 방법 (Cycle별 독립 적용, 1.5*IQR 기준)\n');
fprintf('\n');

% Check which variables exist in the table
availableVars = {};
availableLabels = {};
for idx = 1:length(variableNames)
    if ismember(variableNames{idx}, summaryTable.Properties.VariableNames)
        availableVars{end+1} = variableNames{idx};
        availableLabels{end+1} = variableLabels{idx};
    end
end

if length(availableVars) < 2
    fprintf('ERROR: Insufficient variables for correlation matrix analysis!\n');
    return;
end

fprintf('Available variables: %s\n', strjoin(availableLabels, ', '));
fprintf('\n');

% Extract data for correlation matrix
% Use pairwise deletion: each variable pair uses only rows where both variables are non-NaN
dataMatrix = [];
for idx = 1:length(availableVars)
    dataMatrix = [dataMatrix, summaryTable.(availableVars{idx})];
end

% Calculate correlation matrix with pairwise deletion
% This way, R_60s missing doesn't cause R_1s data to be discarded
[R, P] = corrcoef(dataMatrix, 'rows', 'pairwise');  % Pearson correlation

% Calculate Spearman correlation (rank-based, for monotonic relationships)
[R_spearman, P_spearman] = corr(dataMatrix, 'Type', 'Spearman', 'rows', 'pairwise');

% Count valid pairs for each correlation
nMatrix = zeros(size(R));
for i = 1:size(dataMatrix, 2)
    for j = 1:size(dataMatrix, 2)
        validPair = ~isnan(dataMatrix(:, i)) & ~isnan(dataMatrix(:, j));
        nMatrix(i, j) = sum(validPair);
    end
end

fprintf('Correlation matrix calculated using pairwise deletion\n');
fprintf('(각 변수 쌍별로 유효한 데이터만 사용하여 계산)\n');
fprintf('Method: Pearson (linear) and Spearman (monotonic)\n');
fprintf('Sample sizes for each pair:\n');
fprintf('%12s', '');
for idx = 1:length(availableLabels)
    fprintf('%12s', availableLabels{idx});
end
fprintf('\n');
for i = 1:size(nMatrix, 1)
    fprintf('%12s', availableLabels{i});
    for j = 1:size(nMatrix, 2)
        fprintf('%12d', nMatrix(i, j));
    end
    fprintf('\n');
end
fprintf('\n');

% Calculate 95% confidence intervals for Pearson correlation coefficients
% Using Fisher's z-transformation
alpha = 0.05;
z_crit = norminv(1 - alpha/2);
corrMatrix = R;  % Pearson
pValueMatrix = P;  % Pearson
ciLowerMatrix = zeros(size(R));
ciUpperMatrix = zeros(size(R));

for i = 1:size(R, 1)
    for j = 1:size(R, 2)
        if i == j
            ciLowerMatrix(i, j) = 1;
            ciUpperMatrix(i, j) = 1;
        else
            r = R(i, j);
            n_pair = nMatrix(i, j);  % Use pairwise sample size
            if abs(r) < 1 && n_pair > 3
                % Fisher's z-transformation
                z = 0.5 * log((1 + r) / (1 - r));
                se = 1 / sqrt(n_pair - 3);
                z_lower = z - z_crit * se;
                z_upper = z + z_crit * se;
                % Transform back
                ciLowerMatrix(i, j) = (exp(2*z_lower) - 1) / (exp(2*z_lower) + 1);
                ciUpperMatrix(i, j) = (exp(2*z_upper) - 1) / (exp(2*z_upper) + 1);
            else
                ciLowerMatrix(i, j) = r;
                ciUpperMatrix(i, j) = r;
            end
        end
    end
end

% Calculate 95% confidence intervals for Spearman correlation coefficients
ciLowerMatrix_spearman = zeros(size(R_spearman));
ciUpperMatrix_spearman = zeros(size(R_spearman));

for i = 1:size(R_spearman, 1)
    for j = 1:size(R_spearman, 2)
        if i == j
            ciLowerMatrix_spearman(i, j) = 1;
            ciUpperMatrix_spearman(i, j) = 1;
        else
            r = R_spearman(i, j);
            n_pair = nMatrix(i, j);
            if abs(r) < 1 && n_pair > 3
                % Fisher's z-transformation for Spearman
                z = 0.5 * log((1 + r) / (1 - r));
                se = 1 / sqrt(n_pair - 3);
                z_lower = z - z_crit * se;
                z_upper = z + z_crit * se;
                % Transform back
                ciLowerMatrix_spearman(i, j) = (exp(2*z_lower) - 1) / (exp(2*z_lower) + 1);
                ciUpperMatrix_spearman(i, j) = (exp(2*z_upper) - 1) / (exp(2*z_upper) + 1);
            else
                ciLowerMatrix_spearman(i, j) = r;
                ciUpperMatrix_spearman(i, j) = r;
            end
        end
    end
end

% Store results
correlationResults = struct();
correlationResults.VariableNames = availableVars;
correlationResults.VariableLabels = availableLabels;
% Pearson correlation
correlationResults.CorrelationMatrix = corrMatrix;  % Pearson
correlationResults.PValueMatrix = pValueMatrix;  % Pearson
correlationResults.CILowerMatrix = ciLowerMatrix;  % Pearson
correlationResults.CIUpperMatrix = ciUpperMatrix;  % Pearson
% Spearman correlation
correlationResults.CorrelationMatrix_Spearman = R_spearman;
correlationResults.PValueMatrix_Spearman = P_spearman;
correlationResults.CILowerMatrix_Spearman = ciLowerMatrix_spearman;
correlationResults.CIUpperMatrix_Spearman = ciUpperMatrix_spearman;
correlationResults.SampleSizeMatrix = nMatrix;  % Pairwise sample sizes
correlationResults.DataMatrix = dataMatrix;  % Store outlier-removed data for scatter plots
correlationResults.OutlierRemoval = struct();
correlationResults.OutlierRemoval.Method = 'IQR (Cycle-wise)';
correlationResults.OutlierRemoval.Multiplier = 1.5;
correlationResults.OutlierRemoval.OriginalRowCount = originalRowCount;
correlationResults.OutlierRemoval.FinalRowCount = height(summaryTable);
correlationResults.OutlierRemoval.RowsRemoved = originalRowCount - height(summaryTable);
correlationResults.OutlierRemoval.RemovalPercentage = 100 * (originalRowCount - height(summaryTable)) / originalRowCount;

%% Display correlation matrix summary
fprintf('\n========================================\n');
fprintf('=== Correlation Matrix Summary (Pearson) ===\n');
fprintf('========================================\n');
fprintf('Note: Sample sizes vary by variable pair (pairwise deletion)\n');
fprintf('Minimum sample size: %d\n', min(nMatrix(nMatrix > 0)));
fprintf('Maximum sample size: %d\n', max(nMatrix(:)));
fprintf('Method: Pearson (measures linear relationships)\n');
fprintf('\nCorrelation Matrix (Pearson R):\n');
fprintf('%12s', '');
for idx = 1:length(availableLabels)
    fprintf('%12s', availableLabels{idx});
end
fprintf('\n');
for i = 1:size(corrMatrix, 1)
    fprintf('%12s', availableLabels{i});
    for j = 1:size(corrMatrix, 2)
        fprintf('%12.4f', corrMatrix(i, j));
    end
    fprintf('\n');
end

fprintf('\n========================================\n');
fprintf('=== Correlation Matrix Summary (Spearman) ===\n');
fprintf('========================================\n');
fprintf('Method: Spearman (measures monotonic relationships)\n');
fprintf('\nCorrelation Matrix (Spearman ρ):\n');
fprintf('%12s', '');
for idx = 1:length(availableLabels)
    fprintf('%12s', availableLabels{idx});
end
fprintf('\n');
for i = 1:size(R_spearman, 1)
    fprintf('%12s', availableLabels{i});
    for j = 1:size(R_spearman, 2)
        fprintf('%12.4f', R_spearman(i, j));
    end
    fprintf('\n');
end

% Comparison: Pearson vs Spearman
fprintf('\n========================================\n');
fprintf('=== Pearson vs Spearman Comparison ===\n');
fprintf('========================================\n');
fprintf('%-20s %12s %12s %12s\n', 'Variable Pair', 'Pearson R', 'Spearman ρ', 'Difference');
fprintf('%s\n', repmat('-', 1, 60));
for i = 1:size(corrMatrix, 1)
    for j = i+1:size(corrMatrix, 2)
        pearson_r = corrMatrix(i, j);
        spearman_r = R_spearman(i, j);
        diff = abs(spearman_r - pearson_r);
        fprintf('%-20s %12.4f %12.4f %12.4f\n', ...
            sprintf('%s vs %s', availableLabels{i}, availableLabels{j}), ...
            pearson_r, spearman_r, diff);
    end
end
fprintf('\nInterpretation:\n');
fprintf('- If Spearman > Pearson: Non-linear but monotonic relationship\n');
fprintf('- If Spearman ≈ Pearson: Linear relationship\n');
fprintf('- If Spearman < Pearson: Possible outliers affecting Pearson\n');
fprintf('- Large difference: Check for non-linearity or outliers\n');

%% Detailed Statistics for Each Variable
fprintf('\n\n========================================\n');
fprintf('=== Variable Descriptive Statistics ===\n');
fprintf('========================================\n');
fprintf('%-15s %12s %12s %12s %12s %12s\n', 'Variable', 'Mean', 'Std', 'Min', 'Max', 'Median');
fprintf('%s\n', repmat('-', 1, 80));
for idx = 1:length(availableVars)
    varData = dataMatrix(:, idx);
    fprintf('%-15s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
        availableLabels{idx}, ...
        mean(varData), ...
        std(varData), ...
        min(varData), ...
        max(varData), ...
        median(varData));
end

%% Detailed Pairwise Correlation Analysis
fprintf('\n\n========================================\n');
fprintf('=== Detailed Pairwise Correlations ===\n');
fprintf('========================================\n');
fprintf('%-20s %-20s %8s %8s %8s %15s %8s\n', ...
    'Variable 1', 'Variable 2', 'R', 'R²', 'p-value', '95%% CI', 'Sig');
fprintf('%s\n', repmat('-', 1, 100));

% Significance levels
alpha = 0.05;
sigLevels = {'***', '**', '*', 'ns'};  % p<0.001, p<0.01, p<0.05, not significant

for i = 1:size(corrMatrix, 1)
    for j = i+1:size(corrMatrix, 2)  % Only upper triangle
        r = corrMatrix(i, j);
        r_squared = r^2;
        p = pValueMatrix(i, j);
        ci_low = ciLowerMatrix(i, j);
        ci_up = ciUpperMatrix(i, j);
        
        % Determine significance
        if p < 0.001
            sig = sigLevels{1};
        elseif p < 0.01
            sig = sigLevels{2};
        elseif p < 0.05
            sig = sigLevels{3};
        else
            sig = sigLevels{4};
        end
        
        fprintf('%-20s %-20s %8.4f %8.4f %8.4e %7.4f-%.4f %8s\n', ...
            availableLabels{i}, ...
            availableLabels{j}, ...
            r, r_squared, p, ci_low, ci_up, sig);
    end
end

%% Significant Correlations Summary
fprintf('\n\n========================================\n');
fprintf('=== Significant Correlations (p < 0.05) ===\n');
fprintf('========================================\n');
sigCount = 0;
for i = 1:size(corrMatrix, 1)
    for j = i+1:size(corrMatrix, 2)
        p = pValueMatrix(i, j);
        if p < 0.05
            sigCount = sigCount + 1;
            r = corrMatrix(i, j);
            r_squared = r^2;
            ci_low = ciLowerMatrix(i, j);
            ci_up = ciUpperMatrix(i, j);
            
            % Strength interpretation
            if abs(r) >= 0.9
                strength = 'Very Strong';
            elseif abs(r) >= 0.7
                strength = 'Strong';
            elseif abs(r) >= 0.5
                strength = 'Moderate';
            elseif abs(r) >= 0.3
                strength = 'Weak';
            else
                strength = 'Very Weak';
            end
            
            % Direction
            if r > 0
                direction = 'Positive';
            else
                direction = 'Negative';
            end
            
            fprintf('[%d] %s vs %s\n', sigCount, availableLabels{i}, availableLabels{j});
            fprintf('    R = %.4f (R² = %.4f), %s %s correlation\n', r, r_squared, direction, strength);
            fprintf('    p = %.4e, 95%% CI: [%.4f, %.4f]\n', p, ci_low, ci_up);
            fprintf('\n');
        end
    end
end
if sigCount == 0
    fprintf('No significant correlations found (p < 0.05)\n');
end

%% Strong Correlations Summary (|r| > 0.7)
fprintf('\n\n========================================\n');
fprintf('=== Strong Correlations (|r| > 0.7) ===\n');
fprintf('========================================\n');
strongCount = 0;
for i = 1:size(corrMatrix, 1)
    for j = i+1:size(corrMatrix, 2)
        r = corrMatrix(i, j);
        if abs(r) > 0.7
            strongCount = strongCount + 1;
            p = pValueMatrix(i, j);
            r_squared = r^2;
            ci_low = ciLowerMatrix(i, j);
            ci_up = ciUpperMatrix(i, j);
            
            fprintf('[%d] %s vs %s\n', strongCount, availableLabels{i}, availableLabels{j});
            fprintf('    R = %.4f (R² = %.4f)\n', r, r_squared);
            fprintf('    p = %.4e, 95%% CI: [%.4f, %.4f]\n', p, ci_low, ci_up);
            if p < 0.05
                fprintf('    Status: Significant (p < 0.05)\n');
            else
                fprintf('    Status: Not significant (p >= 0.05)\n');
            end
            fprintf('\n');
        end
    end
end
if strongCount == 0
    fprintf('No strong correlations found (|r| > 0.7)\n');
end

%% Capacity_C3 vs Rchg Correlations (Most Important)
fprintf('\n\n========================================\n');
fprintf('=== Capacity_C3 vs Rchg Correlations ===\n');
fprintf('========================================\n');
Capacity_C3Idx = find(strcmp(availableVars, 'Capacity_C3'));
if ~isempty(Capacity_C3Idx)
    fprintf('Pearson Correlation:\n');
    fprintf('%-15s %8s %8s %12s %15s %8s\n', ...
        'Rchg Interval', 'R', 'R²', 'p-value', '95%% CI', 'Sig');
    fprintf('%s\n', repmat('-', 1, 70));
    
    for idx = 1:length(availableVars)
        if idx ~= Capacity_C3Idx
            r = corrMatrix(Capacity_C3Idx, idx);
            r_squared = r^2;
            p = pValueMatrix(Capacity_C3Idx, idx);
            ci_low = ciLowerMatrix(Capacity_C3Idx, idx);
            ci_up = ciUpperMatrix(Capacity_C3Idx, idx);
            
            % Determine significance
            if p < 0.001
                sig = '***';
            elseif p < 0.01
                sig = '**';
            elseif p < 0.05
                sig = '*';
            else
                sig = 'ns';
            end
            
            fprintf('%-15s %8.4f %8.4f %12.4e %7.4f-%.4f %8s\n', ...
                availableLabels{idx}, r, r_squared, p, ci_low, ci_up, sig);
        end
    end
    
    fprintf('\nSpearman Correlation:\n');
    fprintf('%-15s %8s %8s %12s %15s %8s\n', ...
        'Rchg Interval', 'ρ', 'ρ²', 'p-value', '95%% CI', 'Sig');
    fprintf('%s\n', repmat('-', 1, 70));
    
    for idx = 1:length(availableVars)
        if idx ~= Capacity_C3Idx
            r = R_spearman(Capacity_C3Idx, idx);
            r_squared = r^2;
            p = P_spearman(Capacity_C3Idx, idx);
            ci_low = ciLowerMatrix_spearman(Capacity_C3Idx, idx);
            ci_up = ciUpperMatrix_spearman(Capacity_C3Idx, idx);
            
            % Determine significance
            if p < 0.001
                sig = '***';
            elseif p < 0.01
                sig = '**';
            elseif p < 0.05
                sig = '*';
            else
                sig = 'ns';
            end
            
            fprintf('%-15s %8.4f %8.4f %12.4e %7.4f-%.4f %8s\n', ...
                availableLabels{idx}, r, r_squared, p, ci_low, ci_up, sig);
        end
    end
    
    fprintf('\nComparison (Capacity_C3 vs Rchg):\n');
    fprintf('%-15s %12s %12s %12s\n', 'Rchg Interval', 'Pearson R', 'Spearman ρ', 'Difference');
    fprintf('%s\n', repmat('-', 1, 55));
    for idx = 1:length(availableVars)
        if idx ~= Capacity_C3Idx
            pearson_r = corrMatrix(Capacity_C3Idx, idx);
            spearman_r = R_spearman(Capacity_C3Idx, idx);
            diff = abs(spearman_r - pearson_r);
            fprintf('%-15s %12.4f %12.4f %12.4f\n', ...
                availableLabels{idx}, pearson_r, spearman_r, diff);
        end
    end
end

%% Missing Data Information (결측치 정보)
fprintf('\n\n========================================\n');
fprintf('=== Missing Data Information (결측치 정보) ===\n');
fprintf('========================================\n');
fprintf('설명: Pairwise Deletion 방식을 사용합니다.\n');
fprintf('각 변수 쌍별로 유효한 데이터만 사용하므로, R_60s가 없어도 R_1s 분석에는 영향 없음.\n\n');
fprintf('필터링된 테이블 총 행 수: %d\n', height(summaryTable));

% Find Capacity_C3 index for later use
Capacity_C3Idx = find(strcmp(availableVars, 'Capacity_C3'));

% Check for each variable
fprintf('\n변수별 결측치 현황:\n');
fprintf('%-15s %10s %10s %12s\n', 'Variable', 'Missing', 'Complete', 'Missing %%');
fprintf('%s\n', repmat('-', 1, 50));
for idx = 1:length(availableVars)
    varData = summaryTable.(availableVars{idx});
    missing = sum(isnan(varData));
    complete = sum(~isnan(varData));
    missingPct = 100 * missing / height(summaryTable);
    fprintf('%-15s %10d %10d %11.2f%%\n', availableLabels{idx}, missing, complete, missingPct);
end

% Explain why some rows are excluded
fprintf('\n[참고] Pairwise Deletion 방식:\n');
fprintf('- 각 변수 쌍별로 유효한 데이터만 사용하여 상관계수를 계산합니다.\n');
fprintf('- 예: Capacity_C3 vs R_1s는 R_1s가 있는 모든 행을 사용 (R_60s가 없어도 OK)\n');
fprintf('- 예: Capacity_C3 vs R_60s는 R_60s가 있는 행만 사용\n');
fprintf('- 이 방식으로 데이터 손실을 최소화합니다.\n');
if ~isempty(Capacity_C3Idx)
    fprintf('\n변수 쌍별 샘플 크기:\n');
    r1sIdx = find(strcmp(availableVars, 'R_1s'));
    if ~isempty(r1sIdx)
        fprintf('  Capacity_C3 vs R_1s: %d개\n', nMatrix(Capacity_C3Idx, r1sIdx));
    end
    if ismember('R_60s', availableVars)
        r60sIdx = find(strcmp(availableVars, 'R_60s'));
        if ~isempty(r60sIdx)
            fprintf('  Capacity_C3 vs R_60s: %d개\n', nMatrix(Capacity_C3Idx, r60sIdx));
            if ~isempty(r1sIdx)
                fprintf('  → R_60s가 없는 %d개 행도 R_1s 분석에 사용됨\n', ...
                    nMatrix(Capacity_C3Idx, r1sIdx) - nMatrix(Capacity_C3Idx, r60sIdx));
            end
        end
    end
end

%% Summary Statistics
fprintf('\n\n========================================\n');
fprintf('=== Analysis Summary ===\n');
fprintf('========================================\n');
fprintf('Event Type: %s\n', targetEventType);
fprintf('Outlier Removal: IQR method (Cycle-wise, 1.5*IQR)\n');
fprintf('  Original rows: %d\n', originalRowCount);
fprintf('  Final rows: %d\n', height(summaryTable));
fprintf('  Rows removed: %d (%.2f%%)\n', ...
    originalRowCount - height(summaryTable), ...
    correlationResults.OutlierRemoval.RemovalPercentage);
fprintf('Total Variables: %d\n', length(availableVars));
fprintf('Sample Size Range: %d ~ %d (varies by variable pair)\n', ...
    min(nMatrix(nMatrix > 0)), max(nMatrix(:)));
fprintf('Total Correlation Pairs: %d\n', nchoosek(length(availableVars), 2));
fprintf('Significant Correlations (p < 0.05): %d\n', sigCount);
fprintf('Strong Correlations (|r| > 0.7): %d\n', strongCount);
fprintf('\n');

% Save results
savePath = fullfile(outputDir, 'DriveCycle_CorrelationResults.mat');
save(savePath, 'correlationResults', 'targetEventType');
fprintf('\nResults saved to: %s\n', savePath);

%% Visualization
    fprintf('\n=== Creating Visualization ===\n');
    
    % Create figures directory
    figuresDir = fullfile(outputDir, 'figures', 'CorrelationAnalysis');
    if ~exist(figuresDir, 'dir')
        mkdir(figuresDir);
    end
    
    % Event type label for titles
    eventTypeLabel = '';
    if strcmp(targetEventType, 'Charge')
        eventTypeLabel = ' (Charge)';
    elseif strcmp(targetEventType, 'Discharge')
        eventTypeLabel = ' (Discharge)';
    else
        eventTypeLabel = ' (All Events)';
    end
    
    %% Create correlation matrix heatmap (like the reference image)
    fig = figure('Name', sprintf('Correlation Heatmap%s', eventTypeLabel), ...
                 'Position', [100, 100, 1200, 1000], 'Visible', 'on');
    
    % Prepare correlation matrix data (Pearson)
    corrData = correlationResults.CorrelationMatrix;  % Pearson
    pData = correlationResults.PValueMatrix;  % Pearson
    ciLower = correlationResults.CILowerMatrix;  % Pearson
    ciUpper = correlationResults.CIUpperMatrix;  % Pearson
    varLabels = correlationResults.VariableLabels;
    
    % Create cell labels with r, p, CI, and * for |r|>0.9
    cellLabels = cell(size(corrData));
    for i = 1:size(corrData, 1)
        for j = 1:size(corrData, 2)
            r = corrData(i, j);
            p = pData(i, j);
            ci_low = ciLower(i, j);
            ci_up = ciUpper(i, j);
            
            if i == j
                % Diagonal: self-correlation
                cellLabels{i, j} = sprintf('r=%.2f\np=%.2f\nCI[%.2f,%.2f]', r, p, ci_low, ci_up);
            else
                % Off-diagonal: correlation between different variables
                % Format: r=X.XX, p=X.XXe-YY or p=X.XXX, CI[X.XX, Y.YY]
                if p < 0.001
                    pStr = sprintf('p=%.2e', p);
                elseif p < 0.01
                    pStr = sprintf('p=%.3e', p);
                else
                    pStr = sprintf('p=%.3f', p);
                end
                
                % Add asterisk if |r| > 0.9
                asterisk = '';
                if abs(r) > 0.9
                    asterisk = '*';
                end
                
                cellLabels{i, j} = sprintf('%sr=%.2f\n%s\nCI[%.2f,%.2f]', ...
                    asterisk, r, pStr, ci_low, ci_up);
            end
        end
    end
    
    % Create heatmap using imagesc (more control over cell labels)
    ax = axes('Position', [0.1, 0.1, 0.75, 0.75]);
    imagesc(corrData);
    
    % Create red-blue colormap (red for positive, blue for negative)
    % Custom red-blue colormap: blue (negative) -> white (zero) -> red (positive)
    nColors = 256;
    % Create colormap: blue -> cyan -> white -> yellow -> red
    % First half: blue to white (negative to zero)
    blueVals = linspace(1, 0, nColors/2)';  % Blue decreases
    greenVals = linspace(0, 1, nColors/2)';  % Green increases
    redVals = linspace(0, 1, nColors/2)';   % Red increases
    blueToWhite = [redVals, greenVals, blueVals];
    
    % Second half: white to red (zero to positive)
    redVals2 = ones(nColors/2, 1);           % Red stays at 1
    greenVals2 = linspace(1, 0, nColors/2)'; % Green decreases
    blueVals2 = linspace(1, 0, nColors/2)';  % Blue decreases
    whiteToRed = [redVals2, greenVals2, blueVals2];
    
    redBlueMap = [blueToWhite; whiteToRed];
    colormap(ax, redBlueMap);
    caxis([-1 1]);
    colorbar;
    
    % Set axis labels
    set(ax, 'XTick', 1:length(varLabels), 'XTickLabel', varLabels, ...
            'YTick', 1:length(varLabels), 'YTickLabel', varLabels, ...
            'XTickLabelRotation', 45, 'FontSize', 10);
    
    % Add cell labels
    for i = 1:size(corrData, 1)
        for j = 1:size(corrData, 2)
            text(j, i, cellLabels{i, j}, ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', ...
                 'FontSize', 8, ...
                 'Color', 'white', ...
                 'FontWeight', 'bold');
        end
    end
    
    % Title
    title(sprintf('Correlation Heatmap (Pearson, r, p, 95%% CI; if |r|>0.9)%s\n(IQR Outlier Removal: %d rows removed, %.1f%%)', ...
          eventTypeLabel, correlationResults.OutlierRemoval.RowsRemoved, ...
          correlationResults.OutlierRemoval.RemovalPercentage), ...
          'FontSize', 14, 'FontWeight', 'bold');
    
    % Add note about asterisk
    text(0.5, 0.02, '* indicates |r| > 0.9', ...
         'Units', 'normalized', 'HorizontalAlignment', 'center', ...
         'FontSize', 10, 'FontAngle', 'italic');
    
    % Save figure
    savePath_fig = fullfile(figuresDir, sprintf('Correlation_Heatmap_Pearson%s.fig', eventTypeLabel));
    saveas(fig, savePath_fig);
    fprintf('Saved: %s\n', savePath_fig);
    
    %% Create Spearman correlation matrix heatmap
    fprintf('\n=== Creating Spearman Correlation Heatmap ===\n');
    fig_spearman = figure('Name', sprintf('Correlation Heatmap (Spearman)%s', eventTypeLabel), ...
                 'Position', [100, 100, 1200, 1000], 'Visible', 'on');
    
    % Prepare Spearman correlation matrix data
    corrData_spearman = correlationResults.CorrelationMatrix_Spearman;
    pData_spearman = correlationResults.PValueMatrix_Spearman;
    ciLower_spearman = correlationResults.CILowerMatrix_Spearman;
    ciUpper_spearman = correlationResults.CIUpperMatrix_Spearman;
    
    % Create cell labels with r, p, CI, and * for |r|>0.9
    cellLabels_spearman = cell(size(corrData_spearman));
    for i = 1:size(corrData_spearman, 1)
        for j = 1:size(corrData_spearman, 2)
            r = corrData_spearman(i, j);
            p = pData_spearman(i, j);
            ci_low = ciLower_spearman(i, j);
            ci_up = ciUpper_spearman(i, j);
            
            if i == j
                cellLabels_spearman{i, j} = sprintf('ρ=%.2f\np=%.2f\nCI[%.2f,%.2f]', r, p, ci_low, ci_up);
            else
                if p < 0.001
                    pStr = sprintf('p=%.2e', p);
                elseif p < 0.01
                    pStr = sprintf('p=%.3e', p);
                else
                    pStr = sprintf('p=%.3f', p);
                end
                
                asterisk = '';
                if abs(r) > 0.9
                    asterisk = '*';
                end
                
                cellLabels_spearman{i, j} = sprintf('%sρ=%.2f\n%s\nCI[%.2f,%.2f]', ...
                    asterisk, r, pStr, ci_low, ci_up);
            end
        end
    end
    
    % Create heatmap using imagesc
    ax_spearman = axes('Position', [0.1, 0.1, 0.75, 0.75]);
    imagesc(corrData_spearman);
    
    % Custom red-blue colormap
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
    colormap(ax_spearman, redBlueMap);
    caxis([-1 1]);
    colorbar;
    
    % Set axis labels
    set(ax_spearman, 'XTick', 1:length(varLabels), 'XTickLabel', varLabels, ...
            'YTick', 1:length(varLabels), 'YTickLabel', varLabels, ...
            'XTickLabelRotation', 45, 'FontSize', 10);
    
    % Add cell labels
    for i = 1:size(corrData_spearman, 1)
        for j = 1:size(corrData_spearman, 2)
            text(j, i, cellLabels_spearman{i, j}, ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', ...
                 'FontSize', 8, ...
                 'Color', 'white', ...
                 'FontWeight', 'bold');
        end
    end
    
    % Title
    title(sprintf('Correlation Heatmap (Spearman, ρ, p, 95%% CI; if |ρ|>0.9)%s\n(IQR Outlier Removal: %d rows removed, %.1f%%)', ...
          eventTypeLabel, correlationResults.OutlierRemoval.RowsRemoved, ...
          correlationResults.OutlierRemoval.RemovalPercentage), ...
          'FontSize', 14, 'FontWeight', 'bold');
    
    % Add note about asterisk
    text(0.5, 0.02, '* indicates |ρ| > 0.9', ...
         'Units', 'normalized', 'HorizontalAlignment', 'center', ...
         'FontSize', 10, 'FontAngle', 'italic');
    
    % Save figure
    savePath_fig_spearman = fullfile(figuresDir, sprintf('Correlation_Heatmap_Spearman%s.fig', eventTypeLabel));
    saveas(fig_spearman, savePath_fig_spearman);
    fprintf('Saved: %s\n', savePath_fig_spearman);
    
    %% Create Scatter Plots: Capacity_C3 vs Rchg (with regression line)
    fprintf('\n=== Creating Scatter Plots: Capacity_C3 vs Rchg ===\n');
    
    % Find Capacity_C3 index
    Capacity_C3Idx = find(strcmp(availableVars, 'Capacity_C3'));
    if isempty(Capacity_C3Idx)
        fprintf('WARNING: Capacity_C3 variable not found. Skipping scatter plots.\n');
    else
        % Create scatter plots for each Rchg interval
        for idx = 1:length(availableVars)
            if idx == Capacity_C3Idx
                continue;  % Skip Capacity_C3 vs Capacity_C3
            end
            
            rchgVar = availableVars{idx};
            rchgLabel = availableLabels{idx};
            
            % Extract valid pairs (pairwise deletion)
            validMask = ~isnan(dataMatrix(:, Capacity_C3Idx)) & ~isnan(dataMatrix(:, idx));
            capData = dataMatrix(validMask, Capacity_C3Idx);
            rchgData = dataMatrix(validMask, idx);
            
            if length(capData) < 3
                fprintf('  Skipping %s: insufficient data (n=%d)\n', rchgLabel, length(capData));
                continue;
            end
            
            % Calculate correlation for this pair
            [R_pair, P_pair] = corrcoef(capData, rchgData);
            r_val = R_pair(1, 2);
            p_val = P_pair(1, 2);
            r_squared = r_val^2;
            
            % Linear regression
            p_fit = polyfit(capData, rchgData, 1);
            y_fit = polyval(p_fit, capData);
            
            % Create figure
            fig_scatter = figure('Name', sprintf('Capacity_C3 vs %s%s', rchgLabel, eventTypeLabel), ...
                                'Position', [100, 100, 800, 600], 'Visible', 'on');
            
            % Get Cycle and Channel information for coloring
            validRows = find(validMask);
            cycleData = summaryTable.Cycle(validRows);
            channelData = summaryTable.Channel(validRows);
            
            % Scatter plot with Cycle-based coloring
            uniqueCycles = unique(cycleData);
            cycleColors = lines(length(uniqueCycles));
            hold on;
            
            for cycIdx = 1:length(uniqueCycles)
                cycleNum = uniqueCycles(cycIdx);
                cycleMask = cycleData == cycleNum;
                scatter(capData(cycleMask), rchgData(cycleMask), 80, ...
                    cycleColors(cycIdx, :), 'filled', 'MarkerFaceAlpha', 0.7, ...
                    'DisplayName', sprintf('Cycle %d', cycleNum));
            end
            
            % Regression line
            [capData_sorted, sortIdx] = sort(capData);
            y_fit_sorted = polyval(p_fit, capData_sorted);
            plot(capData_sorted, y_fit_sorted, 'k--', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
            
            % Labels and title
            xlabel('Capacity_C3 (Ah)', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel(sprintf('%s (mΩ)', rchgLabel), 'FontSize', 12, 'FontWeight', 'bold');
            title(sprintf('Capacity_C3 vs %s%s\n(IQR Outlier Removal Applied)', rchgLabel, eventTypeLabel), ...
                  'FontSize', 14, 'FontWeight', 'bold');
            
            % Add statistics text
            statsText = sprintf('R = %.4f\nR² = %.4f\np = %.4e\nn = %d', ...
                              r_val, r_squared, p_val, length(capData));
            if p_val < 0.001
                sigText = '***';
            elseif p_val < 0.01
                sigText = '**';
            elseif p_val < 0.05
                sigText = '*';
            else
                sigText = 'ns';
            end
            statsText = sprintf('%s\nSig: %s', statsText, sigText);
            
            % Position text box (top right)
            text(0.98, 0.98, statsText, ...
                 'Units', 'normalized', ...
                 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'top', ...
                 'FontSize', 10, ...
                 'BackgroundColor', 'white', ...
                 'EdgeColor', 'black', ...
                 'LineWidth', 1);
            
            legend('Location', 'best', 'FontSize', 9);
            grid on;
            hold off;
            
            % Save figure
            savePath_scatter = fullfile(figuresDir, sprintf('Scatter_Capacity_C3_vs_%s%s.fig', ...
                strrep(rchgLabel, ' ', '_'), eventTypeLabel));
            saveas(fig_scatter, savePath_scatter);
            fprintf('  Saved: %s (n=%d, R=%.4f, R²=%.4f, p=%.4e)\n', ...
                rchgLabel, length(capData), r_val, r_squared, p_val);
        end
    end
    
fprintf('\n=== Visualization Complete ===\n');

fprintf('\n=== Correlation Analysis Complete ===\n');

