%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 04_FieldData_CorrelationAnalysis.m
% 상관분석 스크립트 (필드 데이터용)
% 
% 목적: 
% - FieldData_Summary_Table.mat를 불러와서 상관분석 수행
% - 동일한 그룹의 이벤트들로만 상관분석 (각 그룹별로 별도 분석)
% - SOH vs Rchg 상관분석 (그룹별)
% - Year vs Rchg 상관분석 (그룹별)
% - Pearson 상관계수 계산
%
% 입력:
% - FieldData_Summary_Table.mat (02번 스크립트 출력)
%
% 출력:
% - 상관계수 결과 테이블 (콘솔 출력)
% - 상관분석 결과 저장 (그룹별)
% - figures/Group_*/Correlation_*.fig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('=== Field Data Correlation Analysis ===\n');

%% Configuration - Paths
% =========================================================================
scriptDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(scriptDir, 'EventsResults');
inputDir = fullfile(resultsDir, 'FieldData_Analysis');
outputDir = fullfile(resultsDir, 'FieldData_Analysis');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Event Type Selection (선택적)
targetEventType = '';  % 빈 문자열: 모든 이벤트 타입

% 필터링: 모든 연도에 데이터가 있는 그룹만 사용
require_all_years = true;  % true: 모든 연도에 데이터가 있는 그룹만 사용
% =========================================================================

%% Load Summary Table (이상치 제거된 테이블 우선 사용)
summaryTablePath_cleaned = fullfile(inputDir, 'FieldData_Summary_Table_Cleaned.mat');
summaryTablePath_original = fullfile(inputDir, 'FieldData_Summary_Table.mat');

if exist(summaryTablePath_cleaned, 'file')
    load(summaryTablePath_cleaned, 'summaryTable');
    fprintf('Loaded cleaned summary table (IQR outliers removed): %d rows\n', height(summaryTable));
elseif exist(summaryTablePath_original, 'file')
    load(summaryTablePath_original, 'summaryTable');
    fprintf('Loaded original summary table: %d rows\n', height(summaryTable));
    fprintf('WARNING: Cleaned table not found. Using original table.\n');
    fprintf('Please run 03_FieldData_Visualization.m first to generate cleaned table.\n');
else
    fprintf('ERROR: FieldData_Summary_Table.mat not found!\n');
    fprintf('Expected path: %s\n', summaryTablePath_original);
    fprintf('Please run 02_FieldData_DataAggregation.m first to generate the summary table.\n');
    return;
end

% Debug: Check table structure
fprintf('\n=== Table Structure Check ===\n');
fprintf('Table columns: %s\n', strjoin(summaryTable.Properties.VariableNames, ', '));
if ismember('EventType', summaryTable.Properties.VariableNames)
    fprintf('Unique Event Types: %s\n', strjoin(unique(summaryTable.EventType), ', '));
else
    fprintf('WARNING: EventType column not found in table!\n');
end
fprintf('Unique Groups: %d\n', length(unique(summaryTable.Group)));
fprintf('Unique Years: %s\n', mat2str(unique(summaryTable.Year)'));
if ismember('SOH', summaryTable.Properties.VariableNames)
    validSOH = ~isnan(summaryTable.SOH);
    if sum(validSOH) > 0
        fprintf('Data range - SOH: [%.2f, %.2f] %%\n', min(summaryTable.SOH(validSOH)), max(summaryTable.SOH(validSOH)));
    else
        fprintf('WARNING: No valid SOH data found!\n');
    end
end

% Data composition analysis
fprintf('\n=== Data Composition Analysis ===\n');
fprintf('Total rows in summary table: %d\n', height(summaryTable));
fprintf('Unique Group × Year combinations: %d\n', ...
    height(unique(summaryTable(:, {'Group', 'Year'}))));

% Count by Year
fprintf('\nData points by Year:\n');
uniqueYears = unique(summaryTable.Year);
for i = 1:length(uniqueYears)
    yearNum = uniqueYears(i);
    count = sum(summaryTable.Year == yearNum);
    fprintf('  Year %d: %d rows\n', yearNum, count);
end

% Count by Group
fprintf('\nData points by Group:\n');
uniqueGroups_all = unique(summaryTable.Group);
for i = 1:min(10, length(uniqueGroups_all))
    groupName = uniqueGroups_all{i};
    count = sum(strcmp(summaryTable.Group, groupName));
    fprintf('  %s: %d rows\n', groupName, count);
end
if length(uniqueGroups_all) > 10
    fprintf('  ... (showing first 10, total: %d groups)\n', length(uniqueGroups_all));
end

% Event Type 필터링
if ~isempty(targetEventType) && ismember('EventType', summaryTable.Properties.VariableNames)
    originalCount = height(summaryTable);
    summaryTable = summaryTable(strcmp(summaryTable.EventType, targetEventType), :);
    fprintf('\nFiltered to %s events only: %d rows (removed %d rows)\n', ...
        targetEventType, height(summaryTable), originalCount - height(summaryTable));
elseif ismember('EventType', summaryTable.Properties.VariableNames)
    fprintf('\nUsing all event types: %d rows\n', height(summaryTable));
end

%% Create figures directory
figuresDir = fullfile(outputDir, 'figures');
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

%% Get unique groups (모든 연도에 데이터가 있는 그룹만)
all_uniqueGroups = unique(summaryTable.Group);
fprintf('\nTotal unique groups: %d\n', length(all_uniqueGroups));

% 모든 연도 확인
all_required_years = sort(unique(summaryTable.Year));

% 모든 연도에 데이터가 있는 그룹만 필터링
uniqueGroups = {};
for g_idx = 1:length(all_uniqueGroups)
    groupName = all_uniqueGroups{g_idx};
    groupMask = strcmp(summaryTable.Group, groupName);
    groupData = summaryTable(groupMask, :);
    
    % 모든 연도에 데이터가 있는지 확인
    if require_all_years
        groupYears = unique(groupData.Year);
        has_all_years = true;
        missing_years = {};
        for y_idx = 1:length(all_required_years)
            req_year = all_required_years(y_idx);
            if ~any(groupYears == req_year)
                has_all_years = false;
                missing_years{end+1} = sprintf('Y%d', req_year);
            end
        end
        
        if ~has_all_years
            fprintf('  Skipping group %s: missing years %s\n', groupName, strjoin(missing_years, ', '));
            continue;
        end
    end
    
    uniqueGroups{end+1} = groupName;
end
uniqueGroups = uniqueGroups';

fprintf('Qualified groups (all years have data): %d\n', length(uniqueGroups));

%% Define variables for correlation
timeIntervals = {'R_chg_1s', 'R_chg_3s', 'R_chg_5s', 'R_chg_10s', 'R_chg_30s', 'R_chg_60s'};
timeIntervalLabels = {'Rchg 1s', 'Rchg 3s', 'Rchg 5s', 'Rchg 10s', 'Rchg 30s', 'Rchg 60s'};

%% Process each group separately
allCorrelationResults = struct();

% Create red-blue colormap (한 번만 생성)
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

for g_idx = 1:length(uniqueGroups)
    groupName = uniqueGroups{g_idx};
    % Convert group name for display: _ to space (keep decimal points)
    groupName_display = strrep(groupName, '_', ' ');
    
    fprintf('\n\n========================================\n');
    fprintf('=== Processing Group %d/%d: %s ===\n', g_idx, length(uniqueGroups), groupName_display);
    fprintf('========================================\n');
    
    % 그룹별 데이터 필터링
    groupMask = strcmp(summaryTable.Group, groupName);
    groupTable = summaryTable(groupMask, :);
    
    fprintf('Events in this group: %d\n', height(groupTable));
    
    %% Debug: Group data check
    fprintf('\n--- Group Data Check ---\n');
    fprintf('Group: %s\n', groupName_display);
    fprintf('Total events in group: %d\n', height(groupTable));
    fprintf('Years in group: %s\n', mat2str(unique(groupTable.Year)'));
    
    % Check data availability for each variable
    fprintf('Data availability:\n');
    fprintf('  SOH: %d valid / %d total\n', sum(~isnan(groupTable.SOH)), height(groupTable));
    for t_idx = 1:length(timeIntervals)
        timeInterval = timeIntervals{t_idx};
        timeLabel = timeIntervalLabels{t_idx};
        validCount = sum(~isnan(groupTable.(timeInterval)));
        fprintf('  %s: %d valid / %d total\n', timeLabel, validCount, height(groupTable));
    end
    
    % 그룹별 figures 폴더 생성 (안전한 파일명으로 변환)
    groupFolderName = strrep(groupName, '.', '_');
    groupFolderName = strrep(groupFolderName, '-', '_N');  % 마이너스 기호를 _N으로 변환
    groupFolderName = strrep(groupFolderName, ' ', '_');
    groupFigDir = fullfile(figuresDir, ['Group_' groupFolderName]);
    if ~exist(groupFigDir, 'dir')
        mkdir(groupFigDir);
    end
    
    %% SOH vs Rchg Correlation Analysis
    fprintf('\n--- SOH vs Rchg Correlations ---\n');
    
    validMaskSOH = ~isnan(groupTable.SOH);
    if sum(validMaskSOH) == 0
        fprintf('No SOH data available for this group. Skipping SOH correlation analysis.\n');
    else
        fprintf('Valid SOH data points: %d\n', sum(validMaskSOH));
        
        % 상관계수 결과 저장
        sohCorrResults = struct();
        sohCorrResults.Group = groupName;
        sohCorrResults.TimeInterval = {};
        sohCorrResults.R = [];
        sohCorrResults.R_squared = [];
        sohCorrResults.P = [];
        sohCorrResults.CI_lower = [];
        sohCorrResults.CI_upper = [];
        sohCorrResults.N = [];
        
        fprintf('\nCorrelation Method: Pearson (linear relationship)\n');
        fprintf('Sample size: %d (pairwise deletion for each interval)\n', sum(validMaskSOH));
        fprintf('\n%-15s %8s %8s %12s %15s %8s\n', ...
            'Rchg Interval', 'R', 'R²', 'p-value', '95%% CI', 'Sig');
        fprintf('%s\n', repmat('-', 1, 70));
        
        for t_idx = 1:length(timeIntervals)
            timeInterval = timeIntervals{t_idx};
            timeLabel = timeIntervalLabels{t_idx};
            
            validMask = validMaskSOH & ~isnan(groupTable.(timeInterval));
            if sum(validMask) < 3
                fprintf('%-15s %8s\n', timeLabel, 'Insufficient data');
                continue;
            end
            
            soh_data = groupTable.SOH(validMask);
            rchg_data = groupTable.(timeInterval)(validMask);
            
            % Correlation
            [R, P] = corrcoef(soh_data, rchg_data);
            r_val = R(1, 2);
            p_val = P(1, 2);
            r_squared = r_val^2;
            
            % 95% Confidence Interval (Fisher's z-transformation)
            n = length(soh_data);
            if abs(r_val) < 1 && n > 3
                z = 0.5 * log((1 + r_val) / (1 - r_val));
                se = 1 / sqrt(n - 3);
                z_crit = norminv(0.975);
                z_lower = z - z_crit * se;
                z_upper = z + z_crit * se;
                ci_lower = (exp(2*z_lower) - 1) / (exp(2*z_lower) + 1);
                ci_upper = (exp(2*z_upper) - 1) / (exp(2*z_upper) + 1);
            else
                ci_lower = r_val;
                ci_upper = r_val;
            end
            
            % Significance
            if p_val < 0.001
                sig = '***';
            elseif p_val < 0.01
                sig = '**';
            elseif p_val < 0.05
                sig = '*';
            else
                sig = 'ns';
            end
            
            fprintf('%-15s %8.4f %8.4f %12.4e %7.4f-%.4f %8s\n', ...
                timeLabel, r_val, r_squared, p_val, ci_lower, ci_upper, sig);
            
            % 결과 저장
            sohCorrResults.TimeInterval{end+1} = timeLabel;
            sohCorrResults.R(end+1) = r_val;
            sohCorrResults.R_squared(end+1) = r_squared;
            sohCorrResults.P(end+1) = p_val;
            sohCorrResults.CI_lower(end+1) = ci_lower;
            sohCorrResults.CI_upper(end+1) = ci_upper;
            sohCorrResults.N(end+1) = n;
        end
        
        % 그룹별 결과 저장 (안전한 필드명으로 변환)
        group_key_safe = strrep(groupName, '.', '_');
        group_key_safe = strrep(group_key_safe, '-', '_N');  % 마이너스 기호를 _N으로 변환
        group_key_safe = strrep(group_key_safe, ' ', '_');
        allCorrelationResults.(group_key_safe).SOH = sohCorrResults;
    end
    
    %% Year vs Rchg Correlation Analysis
    fprintf('\n--- Year vs Rchg Correlations ---\n');
    
    validMaskYear = ~isnan(groupTable.Year);
    if sum(validMaskYear) == 0
        fprintf('No Year data available for this group. Skipping Year correlation analysis.\n');
    else
        fprintf('Valid Year data points: %d\n', sum(validMaskYear));
        fprintf('Year range: [%d, %d]\n', min(groupTable.Year(validMaskYear)), max(groupTable.Year(validMaskYear)));
        
        % 상관계수 결과 저장
        yearCorrResults = struct();
        yearCorrResults.Group = groupName;
        yearCorrResults.TimeInterval = {};
        yearCorrResults.R = [];
        yearCorrResults.R_squared = [];
        yearCorrResults.P = [];
        yearCorrResults.CI_lower = [];
        yearCorrResults.CI_upper = [];
        yearCorrResults.N = [];
        
        fprintf('\n%-15s %8s %8s %12s %15s %8s\n', ...
            'Rchg Interval', 'R', 'R²', 'p-value', '95%% CI', 'Sig');
        fprintf('%s\n', repmat('-', 1, 70));
        
        for t_idx = 1:length(timeIntervals)
            timeInterval = timeIntervals{t_idx};
            timeLabel = timeIntervalLabels{t_idx};
            
            validMask = validMaskYear & ~isnan(groupTable.(timeInterval));
            if sum(validMask) < 3
                fprintf('%-15s %8s\n', timeLabel, 'Insufficient data');
                continue;
            end
            
            year_data = groupTable.Year(validMask);
            rchg_data = groupTable.(timeInterval)(validMask);
            
            % Correlation
            [R, P] = corrcoef(year_data, rchg_data);
            r_val = R(1, 2);
            p_val = P(1, 2);
            r_squared = r_val^2;
            
            % 95% Confidence Interval
            n = length(year_data);
            if abs(r_val) < 1 && n > 3
                z = 0.5 * log((1 + r_val) / (1 - r_val));
                se = 1 / sqrt(n - 3);
                z_crit = norminv(0.975);
                z_lower = z - z_crit * se;
                z_upper = z + z_crit * se;
                ci_lower = (exp(2*z_lower) - 1) / (exp(2*z_lower) + 1);
                ci_upper = (exp(2*z_upper) - 1) / (exp(2*z_upper) + 1);
            else
                ci_lower = r_val;
                ci_upper = r_val;
            end
            
            % Significance
            if p_val < 0.001
                sig = '***';
            elseif p_val < 0.01
                sig = '**';
            elseif p_val < 0.05
                sig = '*';
            else
                sig = 'ns';
            end
            
            fprintf('%-15s %8.4f %8.4f %12.4e %7.4f-%.4f %8s\n', ...
                timeLabel, r_val, r_squared, p_val, ci_lower, ci_upper, sig);
            
            % 결과 저장
            yearCorrResults.TimeInterval{end+1} = timeLabel;
            yearCorrResults.R(end+1) = r_val;
            yearCorrResults.R_squared(end+1) = r_squared;
            yearCorrResults.P(end+1) = p_val;
            yearCorrResults.CI_lower(end+1) = ci_lower;
            yearCorrResults.CI_upper(end+1) = ci_upper;
            yearCorrResults.N(end+1) = n;
        end
        
        % 그룹별 결과 저장 (안전한 필드명으로 변환)
        group_key_safe = strrep(groupName, '.', '_');
        group_key_safe = strrep(group_key_safe, '-', '_N');  % 마이너스 기호를 _N으로 변환
        group_key_safe = strrep(group_key_safe, ' ', '_');
        allCorrelationResults.(group_key_safe).Year = yearCorrResults;
    end
    
    %% Create Correlation Heatmap (SOH + Year vs Rchg - 8x8 matrix)
    if sum(validMaskSOH) > 0
        fprintf('\n--- Creating Correlation Heatmap (8x8 matrix: SOH + Year + 6 Rchg intervals) ---\n');
        
        % Prepare data matrix: SOH + Year + 6 Rchg intervals
        dataMatrix = [];
        varLabels = {'SOH'};
        
        % SOH data
        soh_data = groupTable.SOH(validMaskSOH);
        dataMatrix = [dataMatrix, soh_data];
        fprintf('  Added SOH: %d data points\n', length(soh_data));
        
        % Year data (align with SOH)
        validMaskYear = ~isnan(groupTable.Year);
        if sum(validMaskYear) > 0
            % Align Year with SOH: use same validMaskSOH indices
            year_data = groupTable.Year(validMaskSOH);
            % Check if Year is valid where SOH is valid
            year_valid_in_soh = ~isnan(year_data);
            if sum(year_valid_in_soh) >= 3
                dataMatrix = [dataMatrix, year_data];
                varLabels{end+1} = 'Year';
                fprintf('  Added Year: %d valid data points\n', sum(year_valid_in_soh));
            else
                fprintf('  Skipped Year: insufficient data (%d valid points, need >= 3)\n', sum(year_valid_in_soh));
            end
        else
            fprintf('  Skipped Year: no valid Year data\n');
        end
        
        % Rchg data for each interval
        % Get indices where SOH is valid (for alignment)
        soh_valid_indices = find(validMaskSOH);
        
        for t_idx = 1:length(timeIntervals)
            timeInterval = timeIntervals{t_idx};
            timeLabel = timeIntervalLabels{t_idx};
            validMask = validMaskSOH & ~isnan(groupTable.(timeInterval));
            
            if sum(validMask) >= 3
                % Extract Rchg data where both SOH and Rchg are valid
                rchg_data = groupTable.(timeInterval)(validMask);
                % Align with SOH data: create array same size as soh_data
                aligned_rchg = NaN(size(soh_data));
                % Find positions in soh_data where this Rchg interval is also valid
                valid_in_soh_data = validMask(validMaskSOH);
                aligned_rchg(valid_in_soh_data) = rchg_data;
                dataMatrix = [dataMatrix, aligned_rchg];
                varLabels{end+1} = timeLabel;
                fprintf('  Added %s: %d valid data points\n', timeLabel, sum(validMask));
            else
                fprintf('  Skipped %s: insufficient data (%d valid points, need >= 3)\n', ...
                    timeLabel, sum(validMask));
            end
        end
        
        fprintf('  Total variables in correlation matrix: %d\n', size(dataMatrix, 2));
        
        % Calculate correlation matrix (pairwise deletion)
        nVars = size(dataMatrix, 2);
        fprintf('  Calculating %dx%d correlation matrix using pairwise deletion...\n', nVars, nVars);
        
        % Pearson correlation
        corrMatrix = NaN(nVars, nVars);
        pMatrix = NaN(nVars, nVars);
        ciLowerMatrix = NaN(nVars, nVars);
        ciUpperMatrix = NaN(nVars, nVars);
        
        % Spearman correlation
        corrMatrix_Spearman = NaN(nVars, nVars);
        pMatrix_Spearman = NaN(nVars, nVars);
        ciLowerMatrix_Spearman = NaN(nVars, nVars);
        ciUpperMatrix_Spearman = NaN(nVars, nVars);
        
        nMatrix = zeros(nVars, nVars);  % Sample size for each pair
        
        for i = 1:nVars
            for j = 1:nVars
                % Pairwise deletion
                validPair = ~isnan(dataMatrix(:, i)) & ~isnan(dataMatrix(:, j));
                n_pair = sum(validPair);
                nMatrix(i, j) = n_pair;
                
                if n_pair >= 3
                    x = dataMatrix(validPair, i);
                    y = dataMatrix(validPair, j);
                    
                    % Pearson correlation
                    [R, P] = corrcoef(x, y);
                    r_val = R(1, 2);
                    p_val = P(1, 2);
                    
                    corrMatrix(i, j) = r_val;
                    pMatrix(i, j) = p_val;
                    
                    % 95% Confidence Interval (Pearson)
                    if abs(r_val) < 1 && n_pair > 3
                        z = 0.5 * log((1 + r_val) / (1 - r_val));
                        se = 1 / sqrt(n_pair - 3);
                        z_crit = norminv(0.975);
                        z_lower = z - z_crit * se;
                        z_upper = z + z_crit * se;
                        ci_lower = (exp(2*z_lower) - 1) / (exp(2*z_lower) + 1);
                        ci_upper = (exp(2*z_upper) - 1) / (exp(2*z_upper) + 1);
                    else
                        ci_lower = r_val;
                        ci_upper = r_val;
                    end
                    ciLowerMatrix(i, j) = ci_lower;
                    ciUpperMatrix(i, j) = ci_upper;
                    
                    % Spearman correlation
                    [R_spear, P_spear] = corr(x, y, 'Type', 'Spearman');
                    if isscalar(R_spear)
                        r_spear_val = R_spear;
                        p_spear_val = P_spear;
                    else
                        r_spear_val = R_spear(1, 2);
                        p_spear_val = P_spear(1, 2);
                    end
                    
                    corrMatrix_Spearman(i, j) = r_spear_val;
                    pMatrix_Spearman(i, j) = p_spear_val;
                    
                    % 95% Confidence Interval (Spearman)
                    if abs(r_spear_val) < 1 && n_pair > 3
                        z_spear = 0.5 * log((1 + r_spear_val) / (1 - r_spear_val));
                        se_spear = 1 / sqrt(n_pair - 3);
                        z_crit = norminv(0.975);
                        z_lower_spear = z_spear - z_crit * se_spear;
                        z_upper_spear = z_spear + z_crit * se_spear;
                        ci_lower_spear = (exp(2*z_lower_spear) - 1) / (exp(2*z_lower_spear) + 1);
                        ci_upper_spear = (exp(2*z_upper_spear) - 1) / (exp(2*z_upper_spear) + 1);
                    else
                        ci_lower_spear = r_spear_val;
                        ci_upper_spear = r_spear_val;
                    end
                    ciLowerMatrix_Spearman(i, j) = ci_lower_spear;
                    ciUpperMatrix_Spearman(i, j) = ci_upper_spear;
                end
            end
        end
        
        fprintf('  Correlation matrix calculated. Sample sizes: min=%d, max=%d\n', ...
            min(nMatrix(nMatrix > 0)), max(nMatrix(:)));
        
        % Display sample sizes matrix (레퍼런스 스타일)
        fprintf('\n  Sample sizes for each pair:\n');
        fprintf('%12s', '');
        for idx = 1:length(varLabels)
            fprintf('%12s', varLabels{idx});
        end
        fprintf('\n');
        for i = 1:size(nMatrix, 1)
            fprintf('%12s', varLabels{i});
            for j = 1:size(nMatrix, 2)
                fprintf('%12d', nMatrix(i, j));
            end
            fprintf('\n');
        end
        
        %% Display correlation matrix summary (like reference)
        fprintf('\n========================================\n');
        fprintf('=== Correlation Matrix Summary (Pearson) - %s ===\n', groupName_display);
        fprintf('========================================\n');
        fprintf('Note: Sample sizes vary by variable pair (pairwise deletion)\n');
        fprintf('Minimum sample size: %d\n', min(nMatrix(nMatrix > 0)));
        fprintf('Maximum sample size: %d\n', max(nMatrix(:)));
        fprintf('Method: Pearson (measures linear relationships)\n');
        fprintf('\nCorrelation Matrix (Pearson R):\n');
        fprintf('%12s', '');
        for idx = 1:length(varLabels)
            fprintf('%12s', varLabels{idx});
        end
        fprintf('\n');
        for i = 1:size(corrMatrix, 1)
            fprintf('%12s', varLabels{i});
            for j = 1:size(corrMatrix, 2)
                if isnan(corrMatrix(i, j))
                    fprintf('%12s', 'N/A');
                else
                    fprintf('%12.4f', corrMatrix(i, j));
                end
            end
            fprintf('\n');
        end
        
        %% Spearman Correlation Matrix Summary
        fprintf('\n========================================\n');
        fprintf('=== Correlation Matrix Summary (Spearman) - %s ===\n', groupName_display);
        fprintf('========================================\n');
        fprintf('Method: Spearman (measures monotonic relationships)\n');
        fprintf('\nCorrelation Matrix (Spearman ρ):\n');
        fprintf('%12s', '');
        for idx = 1:length(varLabels)
            fprintf('%12s', varLabels{idx});
        end
        fprintf('\n');
        for i = 1:size(corrMatrix_Spearman, 1)
            fprintf('%12s', varLabels{i});
            for j = 1:size(corrMatrix_Spearman, 2)
                if isnan(corrMatrix_Spearman(i, j))
                    fprintf('%12s', 'N/A');
                else
                    fprintf('%12.4f', corrMatrix_Spearman(i, j));
                end
            end
            fprintf('\n');
        end
        
        %% Pearson vs Spearman Comparison
        fprintf('\n========================================\n');
        fprintf('=== Pearson vs Spearman Comparison - %s ===\n', groupName_display);
        fprintf('========================================\n');
        fprintf('%-20s %12s %12s %12s\n', 'Variable Pair', 'Pearson R', 'Spearman ρ', 'Difference');
        fprintf('%s\n', repmat('-', 1, 60));
        for i = 1:size(corrMatrix, 1)
            for j = i+1:size(corrMatrix, 2)
                if isnan(corrMatrix(i, j)) || isnan(corrMatrix_Spearman(i, j))
                    continue;
                end
                pearson_r = corrMatrix(i, j);
                spearman_r = corrMatrix_Spearman(i, j);
                diff = abs(spearman_r - pearson_r);
                fprintf('%-20s %12.4f %12.4f %12.4f\n', ...
                    sprintf('%s vs %s', varLabels{i}, varLabels{j}), ...
                    pearson_r, spearman_r, diff);
            end
        end
        fprintf('\nInterpretation:\n');
        fprintf('- If Spearman > Pearson: Non-linear but monotonic relationship\n');
        fprintf('- If Spearman ≈ Pearson: Linear relationship\n');
        fprintf('- If Spearman < Pearson: Possible outliers affecting Pearson\n');
        fprintf('- Large difference: Check for non-linearity or outliers\n');
        
        %% Variable Descriptive Statistics
        fprintf('\n========================================\n');
        fprintf('=== Variable Descriptive Statistics - %s ===\n', groupName_display);
        fprintf('========================================\n');
        fprintf('%-15s %12s %12s %12s %12s %12s\n', 'Variable', 'Mean', 'Std', 'Min', 'Max', 'Median');
        fprintf('%s\n', repmat('-', 1, 80));
        for idx = 1:length(varLabels)
            varData = dataMatrix(:, idx);
            validData = varData(~isnan(varData));
            if ~isempty(validData)
                fprintf('%-15s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
                    varLabels{idx}, ...
                    mean(validData), ...
                    std(validData), ...
                    min(validData), ...
                    max(validData), ...
                    median(validData));
            else
                fprintf('%-15s %12s\n', varLabels{idx}, 'No data');
            end
        end
        
        %% Detailed Pairwise Correlations
        fprintf('\n========================================\n');
        fprintf('=== Detailed Pairwise Correlations - %s ===\n', groupName_display);
        fprintf('========================================\n');
        fprintf('%-20s %-20s %8s %8s %12s %15s %8s\n', ...
            'Variable 1', 'Variable 2', 'R', 'R²', 'p-value', '95%% CI', 'Sig');
        fprintf('%s\n', repmat('-', 1, 100));
        
        sigLevels = {'***', '**', '*', 'ns'};  % p<0.001, p<0.01, p<0.05, not significant
        
        for i = 1:size(corrMatrix, 1)
            for j = i+1:size(corrMatrix, 2)  % Only upper triangle
                if isnan(corrMatrix(i, j))
                    continue;
                end
                r = corrMatrix(i, j);
                r_squared = r^2;
                p = pMatrix(i, j);
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
                
                fprintf('%-20s %-20s %8.4f %8.4f %12.4e %7.4f-%.4f %8s\n', ...
                    varLabels{i}, ...
                    varLabels{j}, ...
                    r, r_squared, p, ci_low, ci_up, sig);
            end
        end
        
        %% Significant Correlations Summary
        fprintf('\n========================================\n');
        fprintf('=== Significant Correlations (p < 0.05) - %s ===\n', groupName_display);
        fprintf('========================================\n');
        sigCount = 0;
        for i = 1:size(corrMatrix, 1)
            for j = i+1:size(corrMatrix, 2)
                if isnan(corrMatrix(i, j)) || isnan(pMatrix(i, j))
                    continue;
                end
                p = pMatrix(i, j);
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
                    
                    fprintf('[%d] %s vs %s\n', sigCount, varLabels{i}, varLabels{j});
                    fprintf('    R = %.4f (R² = %.4f), %s %s correlation\n', r, r_squared, direction, strength);
                    fprintf('    p = %.4e, 95%% CI: [%.4f, %.4f], n = %d\n', p, ci_low, ci_up, nMatrix(i, j));
                    fprintf('\n');
                end
            end
        end
        if sigCount == 0
            fprintf('No significant correlations found (p < 0.05)\n');
        end
        
        %% Strong Correlations Summary (|r| > 0.7)
        fprintf('\n========================================\n');
        fprintf('=== Strong Correlations (|r| > 0.7) - %s ===\n', groupName_display);
        fprintf('========================================\n');
        strongCount = 0;
        for i = 1:size(corrMatrix, 1)
            for j = i+1:size(corrMatrix, 2)
                if isnan(corrMatrix(i, j))
                    continue;
                end
                r = corrMatrix(i, j);
                if abs(r) > 0.7
                    strongCount = strongCount + 1;
                    p = pMatrix(i, j);
                    r_squared = r^2;
                    ci_low = ciLowerMatrix(i, j);
                    ci_up = ciUpperMatrix(i, j);
                    
                    fprintf('[%d] %s vs %s\n', strongCount, varLabels{i}, varLabels{j});
                    fprintf('    R = %.4f (R² = %.4f)\n', r, r_squared);
                    fprintf('    p = %.4e, 95%% CI: [%.4f, %.4f], n = %d\n', p, ci_low, ci_up, nMatrix(i, j));
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
        
        %% SOH vs Rchg Correlations (Most Important)
        sohIdx = find(strcmp(varLabels, 'SOH'));
        if ~isempty(sohIdx)
            fprintf('\n========================================\n');
            fprintf('=== SOH vs Rchg Correlations - %s ===\n', groupName_display);
            fprintf('========================================\n');
            fprintf('Pearson Correlation:\n');
            fprintf('%-15s %8s %8s %12s %15s %8s %8s\n', ...
                'Rchg Interval', 'R', 'R²', 'p-value', '95%% CI', 'Sig', 'n');
            fprintf('%s\n', repmat('-', 1, 80));
            
            for idx = 1:length(varLabels)
                if idx ~= sohIdx
                    if isnan(corrMatrix(sohIdx, idx))
                        continue;
                    end
                    r = corrMatrix(sohIdx, idx);
                    r_squared = r^2;
                    p = pMatrix(sohIdx, idx);
                    ci_low = ciLowerMatrix(sohIdx, idx);
                    ci_up = ciUpperMatrix(sohIdx, idx);
                    n = nMatrix(sohIdx, idx);
                    
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
                    
                    fprintf('%-15s %8.4f %8.4f %12.4e %7.4f-%.4f %8s %8d\n', ...
                        varLabels{idx}, r, r_squared, p, ci_low, ci_up, sig, n);
                end
            end
        end
        
        %% Missing Data Information
        fprintf('\n========================================\n');
        fprintf('=== Missing Data Information - %s ===\n', groupName_display);
        fprintf('========================================\n');
        fprintf('설명: Pairwise Deletion 방식을 사용합니다.\n');
        fprintf('각 변수 쌍별로 유효한 데이터만 사용하므로, Rchg_60s가 없어도 Rchg_1s 분석에는 영향 없음.\n\n');
        fprintf('그룹 총 행 수: %d\n', height(groupTable));
        
        fprintf('\n변수별 결측치 현황:\n');
        fprintf('%-15s %10s %10s %12s\n', 'Variable', 'Missing', 'Complete', 'Missing %%');
        fprintf('%s\n', repmat('-', 1, 50));
        for idx = 1:length(varLabels)
            varData = dataMatrix(:, idx);
            missing = sum(isnan(varData));
            complete = sum(~isnan(varData));
            total = length(varData);
            missing_pct = (missing / total) * 100;
            fprintf('%-15s %10d %10d %11.2f%%\n', varLabels{idx}, missing, complete, missing_pct);
        end
        
        fprintf('\n변수 쌍별 샘플 크기:\n');
        fprintf('%-20s %-20s %10s\n', 'Variable 1', 'Variable 2', 'Sample Size');
        fprintf('%s\n', repmat('-', 1, 55));
        for i = 1:size(nMatrix, 1)
            for j = i+1:size(nMatrix, 2)
                fprintf('%-20s %-20s %10d\n', varLabels{i}, varLabels{j}, nMatrix(i, j));
            end
        end
        
        %% Analysis Summary for this group
        fprintf('\n========================================\n');
        fprintf('=== Analysis Summary - %s ===\n', groupName_display);
        fprintf('========================================\n');
        fprintf('  Total rows in group: %d\n', height(groupTable));
        fprintf('  Variables in correlation matrix: %d\n', nVars);
        fprintf('  Total Correlation Pairs: %d\n', nchoosek(nVars, 2));
        fprintf('  Significant Correlations (p < 0.05): %d\n', sigCount);
        fprintf('  Strong Correlations (|r| > 0.7): %d\n', strongCount);
        fprintf('  Minimum sample size: %d\n', min(nMatrix(nMatrix > 0)));
        fprintf('  Maximum sample size: %d\n', max(nMatrix(:)));
        
        % Create cell labels
        cellLabels = cell(size(corrMatrix));
        for i = 1:size(corrMatrix, 1)
            for j = 1:size(corrMatrix, 2)
                if isnan(corrMatrix(i, j))
                    cellLabels{i, j} = 'N/A';
                elseif i == j
                    % Diagonal: self-correlation
                    cellLabels{i, j} = sprintf('r=%.2f\np=%.2f\nCI[%.2f,%.2f]', ...
                        corrMatrix(i, j), pMatrix(i, j), ciLowerMatrix(i, j), ciUpperMatrix(i, j));
                else
                    % Off-diagonal
                    r = corrMatrix(i, j);
                    p = pMatrix(i, j);
                    ci_low = ciLowerMatrix(i, j);
                    ci_up = ciUpperMatrix(i, j);
                    
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
                    
                    cellLabels{i, j} = sprintf('%sr=%.2f\n%s\nCI[%.2f,%.2f]', ...
                        asterisk, r, pStr, ci_low, ci_up);
                end
            end
        end
        
        % Convert group name for display: _ to space (keep decimal points)
        groupName_display = strrep(groupName, '_', ' ');
        
        % Create individual figure for this group's heatmap
        fig_heatmap_group = figure('Name', sprintf('Correlation Heatmap - %s', groupName_display), ...
                                   'Position', [100 + (g_idx-1)*50, 100 + (g_idx-1)*50, 1200, 1000], ...
                                   'Visible', 'on');
        
        % Create axes
        ax = axes('Parent', fig_heatmap_group);
        
        % Create heatmap using imagesc
        imagesc(ax, corrMatrix);
        
        % Apply colormap
        colormap(ax, redBlueMap);
        caxis(ax, [-1 1]);
        
        % Set axis labels
        set(ax, 'XTick', 1:length(varLabels), 'XTickLabel', varLabels, ...
                'YTick', 1:length(varLabels), 'YTickLabel', varLabels, ...
                'XTickLabelRotation', 45, 'FontSize', 10);
        
        % Add cell labels with dynamic text color based on correlation value
        for i = 1:size(corrMatrix, 1)
            for j = 1:size(corrMatrix, 2)
                if ~isnan(corrMatrix(i, j))
                    % Determine text color based on correlation value
                    % Low |r| values (near 0) have white/light background -> use black text
                    % High |r| values (near ±1) have dark background -> use white text
                    abs_r = abs(corrMatrix(i, j));
                    if abs_r < 0.3
                        textColor = 'black';  % Light background, use black text
                    else
                        textColor = 'white';  % Dark background, use white text
                    end
                    
                    text(ax, j, i, cellLabels{i, j}, ...
                         'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'middle', ...
                         'FontSize', 7, ...
                         'Color', textColor, ...
                         'FontWeight', 'bold');
                end
            end
        end
        
        % Add overall title (include group name, remove axes title to avoid overlap)
        sgtitle(fig_heatmap_group, sprintf('Correlation Heatmap (r, p, 95%% CI; if |r|>0.9) - %s', groupName_display), ...
                'FontSize', 16, 'FontWeight', 'bold');
        
        % Add note about asterisk
        annotation(fig_heatmap_group, 'textbox', [0.4, 0.02, 0.2, 0.03], ...
                   'String', '* indicates |r| > 0.9', ...
                   'FontSize', 10, 'FontAngle', 'italic', ...
                   'HorizontalAlignment', 'center', 'EdgeColor', 'none');
        
        % Add colorbar
        cb = colorbar(ax);
        cb.Position = [0.92, 0.1, 0.02, 0.8];
        colormap(ax, redBlueMap);
        caxis(ax, [-1 1]);
        
        % Save individual heatmap figure
        savePath_heatmap_group = fullfile(groupFigDir, sprintf('Correlation_Heatmap_%s.fig', groupFolderName));
        saveas(fig_heatmap_group, savePath_heatmap_group);
        fprintf('\nSaved heatmap: %s\n', savePath_heatmap_group);
        close(fig_heatmap_group);
    end
    
    %% Create Scatter Plots with Regression (SOH vs Rchg)
    if sum(validMaskSOH) > 0
        fprintf('\n--- Creating Scatter Plots (SOH vs Rchg) ---\n');
        
        % Convert group name for display: _ to space (keep decimal points)
        groupName_display = strrep(groupName, '_', ' ');
        
        fig = figure('Name', sprintf('SOH vs Rchg Scatter - %s', groupName_display), ...
                     'Position', [100, 100, 1600, 1000], 'Visible', 'on');
        
        for t_idx = 1:length(timeIntervals)
            timeInterval = timeIntervals{t_idx};
            timeLabel = timeIntervalLabels{t_idx};
            
            % 서브플롯 그리드: 6개 시간 간격 (2x3)
            subplot(2, 3, t_idx);
            validMask = validMaskSOH & ~isnan(groupTable.(timeInterval));
            
            if sum(validMask) >= 3
                soh_data = groupTable.SOH(validMask);
                rchg_data = groupTable.(timeInterval)(validMask);
                
                scatter(soh_data, rchg_data, 50, 'filled', 'MarkerFaceAlpha', 0.6);
                hold on;
                
                % Linear regression
                p_fit = polyfit(soh_data, rchg_data, 1);
                x_fit = sort(unique(soh_data));
                y_fit = polyval(p_fit, x_fit);
                plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
                
                % Correlation
                [R, P] = corrcoef(soh_data, rchg_data);
                r_val = R(1, 2);
                p_val = P(1, 2);
                r_squared = r_val^2;
                
                % Statistics text
                statsText = sprintf('R = %.4f\nR² = %.4f\np = %.4e\nn = %d', ...
                                  r_val, r_squared, p_val, length(soh_data));
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
                
                text(0.05, 0.95, statsText, 'Units', 'normalized', ...
                     'VerticalAlignment', 'top', 'FontSize', 10, ...
                     'BackgroundColor', 'white', 'EdgeColor', 'black');
                
                xlabel('SOH (%)', 'FontSize', 12, 'FontWeight', 'bold');
                ylabel(sprintf('%s (mΩ)', timeLabel), 'FontSize', 12, 'FontWeight', 'bold');
                title(timeLabel, 'FontSize', 14, 'FontWeight', 'bold');
                grid on;
                hold off;
            else
                text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'FontSize', 14);
                title(timeLabel, 'FontSize', 14, 'FontWeight', 'bold');
            end
        end
        
        sgtitle(sprintf('SOH vs Rchg Correlation - %s', groupName_display), ...
                'FontSize', 16, 'FontWeight', 'bold');
        
        % 저장
        savePath = fullfile(groupFigDir, sprintf('Correlation_Scatter_SOH_%s.fig', groupFolderName));
        saveas(fig, savePath);
        fprintf('Saved: %s\n', savePath);
        close(fig);
    end
end

% Note: Each group's heatmap is now saved as a separate figure file
fprintf('\n=== All heatmaps saved as individual figures ===\n');

%% Save All Correlation Results
fprintf('\n\n=== Saving Correlation Results ===\n');
savePath = fullfile(outputDir, 'FieldData_CorrelationResults.mat');
save(savePath, 'allCorrelationResults');
fprintf('Results saved to: %s\n', savePath);

fprintf('\n=== Correlation Analysis Complete ===\n');

