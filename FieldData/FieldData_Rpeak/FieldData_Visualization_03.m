%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 03_FieldData_Visualization.m
% 통합 시각화 (필드 데이터용)
% 
% 목적: 
% - FieldData_Summary_Table.mat를 로드
% - Year vs Rchg 통합 그래프 생성 (레퍼런스 스타일)
% - Group별로 subplot 생성
% - Color/Marker: Year별로 구분
%
% 입력:
% - FieldData_Summary_Table.mat (02번 스크립트 출력)
%
% 출력:
% - figures/Integrated_SOH_Rchg/SOH_vs_R_*.fig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('=== Field Data Visualization ===\n');

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
targetEventType = '';  % 빈 문자열: 모든 이벤트 타입, 'charge' 또는 'discharge': 특정 타입만

% 최대 그룹 수 (상위 N개 그룹만 시각화, 0이면 모든 그룹)
maxGroupsToPlot = 10;  % 상위 10개 그룹만 시각화 (이벤트 수 기준)

% 필터링: 모든 연도에 데이터가 있는 그룹만 사용
require_all_years = true;  % true: 모든 연도에 데이터가 있는 그룹만 사용
% =========================================================================

%% Load Summary Table
summaryTablePath = fullfile(inputDir, 'FieldData_Summary_Table.mat');
if ~exist(summaryTablePath, 'file')
    fprintf('ERROR: FieldData_Summary_Table.mat not found!\n');
    fprintf('Expected path: %s\n', summaryTablePath);
    fprintf('Please run 02_FieldData_DataAggregation.m first to generate the summary table.\n');
    return;
end

load(summaryTablePath, 'summaryTable');
fprintf('Loaded summary table: %d rows\n', height(summaryTable));

% Event Type 필터링
if ~isempty(targetEventType) && ismember('EventType', summaryTable.Properties.VariableNames)
    summaryTable = summaryTable(strcmp(summaryTable.EventType, targetEventType), :);
    fprintf('Filtered to %s events only: %d rows\n', targetEventType, height(summaryTable));
end

% Debug: Check table structure
fprintf('\n=== Table Structure Check ===\n');
fprintf('Table columns: %s\n', strjoin(summaryTable.Properties.VariableNames, ', '));
if ismember('EventType', summaryTable.Properties.VariableNames)
    fprintf('Unique Event Types: %s\n', strjoin(unique(summaryTable.EventType), ', '));
end
fprintf('Unique Groups: %d\n', length(unique(summaryTable.Group)));
fprintf('Unique Years: %s\n', mat2str(unique(summaryTable.Year)'));
if ismember('SOH', summaryTable.Properties.VariableNames)
    validSOH = ~isnan(summaryTable.SOH);
    fprintf('Data range - SOH: [%.2f, %.2f] %%\n', min(summaryTable.SOH(validSOH)), max(summaryTable.SOH(validSOH)));
end

%% Create figures directory
figuresDir = fullfile(outputDir, 'figures', 'Integrated_SOH_Rchg');
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

%% Define time intervals to plot
timeIntervals = {'R_chg_1s', 'R_chg_3s', 'R_chg_5s', 'R_chg_10s', 'R_chg_30s', 'R_chg_60s'};
timeIntervalLabels = {'Rchg 1s', 'Rchg 3s', 'Rchg 5s', 'Rchg 10s', 'Rchg 30s', 'Rchg 60s'};

%% Get unique groups and select top N groups by event count
uniqueGroups = unique(summaryTable.Group);
fprintf('\n=== Group Selection ===\n');
fprintf('Total unique groups: %d\n', length(uniqueGroups));

% 모든 연도 확인
all_required_years = sort(unique(summaryTable.Year));

% 그룹별 이벤트 수 계산 및 모든 연도에 데이터가 있는지 확인
groupCounts = zeros(length(uniqueGroups), 1);
qualifiedGroups = {};
for g_idx = 1:length(uniqueGroups)
    groupName = uniqueGroups{g_idx};
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
            fprintf('  Skipping group %s: missing years %s (has %d years: %s)\n', ...
                groupName, strjoin(missing_years, ', '), length(groupYears), mat2str(sort(groupYears)'));
            continue;
        end
    end
    
    groupCounts(length(qualifiedGroups)+1) = height(groupData);
    qualifiedGroups{end+1} = groupName;
end

% Qualified groups만 사용
uniqueGroups = qualifiedGroups';
groupCounts = groupCounts(1:length(uniqueGroups));
fprintf('Qualified groups (all years have data): %d\n', length(uniqueGroups));

% 상위 N개 그룹 선택
if maxGroupsToPlot > 0 && maxGroupsToPlot < length(uniqueGroups)
    [~, sortIdx] = sort(groupCounts, 'descend');
    selectedGroups = uniqueGroups(sortIdx(1:maxGroupsToPlot));
    fprintf('Selected top %d groups by event count:\n', maxGroupsToPlot);
    for g_idx = 1:length(selectedGroups)
        groupIdx = find(strcmp(uniqueGroups, selectedGroups{g_idx}), 1);
        fprintf('  %s: %d events\n', selectedGroups{g_idx}, groupCounts(groupIdx));
    end
else
    selectedGroups = uniqueGroups;
    fprintf('Using all %d groups\n', length(selectedGroups));
end

%% Get unique Years for marker mapping (레퍼런스: Channel 마커)
uniqueYears = sort(unique(summaryTable.Year));
fprintf('\nUnique Years: %s\n', mat2str(uniqueYears'));

%% IQR 이상치 제거 (그룹별, 연도별로 적용)
fprintf('\n=== IQR Outlier Removal (by Group and Year) ===\n');
fprintf('Original table size: %d rows\n', height(summaryTable));

% 저항값 임계값 설정
resistance_threshold_mOhm = 5.0;  % 5mΩ 이상인 저항값 제거

% 그룹별로 이상치 제거
all_cleaned_tables = {};
total_outliers_removed = 0;

for g_idx = 1:length(selectedGroups)
    groupName = selectedGroups{g_idx};
    groupMask = strcmp(summaryTable.Group, groupName);
    groupTable = summaryTable(groupMask, :);
    
    if height(groupTable) == 0
        continue;
    end
    
    fprintf('\n  Processing group: %s (%d rows)\n', groupName, height(groupTable));
    
    % 그룹 내에서 연도별 IQR 이상치 제거
    outlier_idx = false(height(groupTable), 1);
    unique_years = unique(groupTable.Year);
    
    for t_idx = 1:length(timeIntervals)
        timeInterval = timeIntervals{t_idx};
        if ~ismember(timeInterval, groupTable.Properties.VariableNames)
            continue;
        end
        
        for y_idx = 1:length(unique_years)
            year_val = unique_years(y_idx);
            year_mask = groupTable.Year == year_val;
            year_data = groupTable.(timeInterval)(year_mask);
            year_data_valid = year_data(~isnan(year_data));
            
            % 레퍼런스와 동일: 최소 4개 이상의 데이터가 있을 때만 이상치 제거 (IQR 계산을 위해 최소 4개 필요)
            if length(year_data_valid) >= 4
                Q1 = prctile(year_data_valid, 25);
                Q3 = prctile(year_data_valid, 75);
                IQR = Q3 - Q1;
                lower_bound = Q1 - 1.5 * IQR;
                upper_bound = Q3 + 1.5 * IQR;
                
                year_outliers = year_mask & (groupTable.(timeInterval) < lower_bound | groupTable.(timeInterval) > upper_bound);
                outlier_idx = outlier_idx | year_outliers;
                if sum(year_outliers) > 0
                    fprintf('    %s, Year %d, %s: Removed %d outliers (bounds: [%.4f, %.4f], n=%d)\n', ...
                        groupName, year_val, timeInterval, sum(year_outliers), lower_bound, upper_bound, length(year_data_valid));
                end
            else
                fprintf('    %s, Year %d, %s: Skipped (insufficient data: %d points, need >= 4)\n', ...
                    groupName, year_val, timeInterval, length(year_data_valid));
            end
        end
    end
    
    % 추가 필터: 저항값이 5mΩ 이상인 행 제거
    resistance_outlier_idx = false(height(groupTable), 1);
    for t_idx = 1:length(timeIntervals)
        timeInterval = timeIntervals{t_idx};
        if ismember(timeInterval, groupTable.Properties.VariableNames)
            % 해당 Rchg 값이 5mΩ 이상인 행 표시
            resistance_outlier_idx = resistance_outlier_idx | (groupTable.(timeInterval) >= resistance_threshold_mOhm);
        end
    end
    
    if sum(resistance_outlier_idx) > 0
        fprintf('    %s: Removed %d rows with resistance >= %.1f mΩ\n', ...
            groupName, sum(resistance_outlier_idx), resistance_threshold_mOhm);
    end
    
    % IQR 이상치와 저항값 임계값 이상치 모두 제거
    outlier_idx = outlier_idx | resistance_outlier_idx;
    
    % 그룹별 이상치 제거된 테이블 생성
    groupTable_cleaned = groupTable(~outlier_idx, :);
    outliers_removed = sum(outlier_idx);
    total_outliers_removed = total_outliers_removed + outliers_removed;
    
    fprintf('    After outlier removal: %d rows remaining (removed %d outliers)\n', ...
        height(groupTable_cleaned), outliers_removed);
    
    all_cleaned_tables{end+1} = groupTable_cleaned;
end

% 모든 그룹의 이상치 제거된 테이블 합치기
if ~isempty(all_cleaned_tables)
    summaryTable_cleaned = vertcat(all_cleaned_tables{:});
    fprintf('\nTotal outliers removed: %d\n', total_outliers_removed);
    fprintf('After outlier removal: %d rows remaining (from %d original rows)\n', ...
        height(summaryTable_cleaned), height(summaryTable));
    
    % 이상치 제거된 테이블로 교체
    summaryTable = summaryTable_cleaned;
else
    fprintf('WARNING: No cleaned tables generated!\n');
end

%% Save Cleaned Tables by Group (그룹별 이상치 제거된 테이블 저장)
fprintf('\n=== Saving Cleaned Tables by Group ===\n');
cleanedTablesDir = fullfile(outputDir, 'Cleaned_Tables_ByGroup');
if ~exist(cleanedTablesDir, 'dir')
    mkdir(cleanedTablesDir);
end

for g_idx = 1:length(selectedGroups)
    groupName = selectedGroups{g_idx};
    groupMask = strcmp(summaryTable.Group, groupName);
    groupTable_cleaned = summaryTable(groupMask, :);
    
    if height(groupTable_cleaned) > 0
        % 안전한 파일명으로 변환
        groupFileName = strrep(groupName, '.', '_');
        groupFileName = strrep(groupFileName, '-', '_N');
        groupFileName = strrep(groupFileName, ' ', '_');
        
        savePath = fullfile(cleanedTablesDir, sprintf('CleanedTable_%s.mat', groupFileName));
        save(savePath, 'groupTable_cleaned');
        fprintf('  Saved: %s (%d rows)\n', groupFileName, height(groupTable_cleaned));
    end
end

% 전체 이상치 제거된 테이블 저장
savePath_all = fullfile(outputDir, 'FieldData_Summary_Table_Cleaned.mat');
save(savePath_all, 'summaryTable');
fprintf('\nSaved overall cleaned table: %s\n', savePath_all);

%% Define color map for Groups (레퍼런스: DC Profile 색상)
groupColors = lines(length(selectedGroups));

%% Define marker styles for Years (레퍼런스: Channel 마커)
yearMarkers = {'o', 's', '^', 'd', 'v', '>', '<', 'p'};
if length(uniqueYears) > length(yearMarkers)
    % If more years than markers, cycle through markers
    yearMarkers = repmat(yearMarkers, 1, ceil(length(uniqueYears) / length(yearMarkers)));
end

%% Create figures (one for each time interval)
for fig_idx = 1:length(timeIntervals)
    timeInterval = timeIntervals{fig_idx};
    timeLabel = timeIntervalLabels{fig_idx};
    
    % Check if column exists
    if ~ismember(timeInterval, summaryTable.Properties.VariableNames)
        fprintf('WARNING: Column %s not found in summary table. Skipping.\n', timeInterval);
        continue;
    end
    
    % Debug: Check data availability for this time interval
    validData = ~isnan(summaryTable.(timeInterval)) & ~isnan(summaryTable.Year);
    fprintf('\n--- Creating plot for %s ---\n', timeLabel);
    fprintf('Valid data points: %d / %d\n', sum(validData), height(summaryTable));
    
    if sum(validData) == 0
        fprintf('WARNING: No valid data for %s. Skipping plot.\n', timeLabel);
        continue;
    end
    
    % Calculate subplot grid size (레퍼런스: 2x4 그리드)
    nGroups = length(selectedGroups);
    if nGroups <= 4
        nRows = 2;
        nCols = 2;
    elseif nGroups <= 6
        nRows = 2;
        nCols = 3;
    elseif nGroups <= 8
        nRows = 2;
        nCols = 4;
    elseif nGroups <= 9
        nRows = 3;
        nCols = 3;
    elseif nGroups <= 12
        nRows = 3;
        nCols = 4;
    else
        nRows = 4;
        nCols = 4;
    end
    
    % Create figure with subplots (Group별로 subplot, 레퍼런스: DC Profile별 subplot)
    fig = figure('Name', sprintf('Year vs %s (Group별)', timeLabel), ...
                 'Position', [100 + (fig_idx-1)*50, 100 + (fig_idx-1)*50, 1600, 900], ...
                 'Visible', 'on');
    
    totalPointsPlotted = 0;
    
    % Plot each Group in a separate subplot (레퍼런스: DC Profile별 subplot)
    for g_idx = 1:min(nGroups, nRows*nCols)
        groupName = selectedGroups{g_idx};
        groupColor = groupColors(g_idx, :);
        
        % Convert group name for display: _ to space (keep decimal points)
        groupName_display = strrep(groupName, '_', ' ');
        
        subplot(nRows, nCols, g_idx);
        hold on;
        
        legendEntries = {};
        legendHandles = [];
        
        % Calculate year-wise means for trend line
        yearMeans = [];
        yearMeanYears = [];
        yearCounts = [];  % 개수 저장용
        
        % Plot all years for this group (레퍼런스: Channel별 마커)
        for y_idx = 1:length(uniqueYears)
            year = uniqueYears(y_idx);
            marker = yearMarkers{y_idx};
            
            % Filter data for this Group and Year
            mask = strcmp(summaryTable.Group, groupName) & ...
                   summaryTable.Year == year & ...
                   ~isnan(summaryTable.(timeInterval)) & ...
                   ~isnan(summaryTable.Year);
            
            if sum(mask) > 0
                yearData = summaryTable.Year(mask);
                rchgData = summaryTable.(timeInterval)(mask);
                
                % Sort by Year for line connection (though all should be same year)
                [yearData_sorted, sortIdx] = sort(yearData);
                rchgData_sorted = rchgData(sortIdx);
                
                % Plot scatter points with group-specific color, year-specific marker
                % 레퍼런스: DC 색상 기반, Channel별 마커
                yearColor = groupColor;  % Use group color as base
                % Make year colors slightly different by adjusting brightness (레퍼런스 방식)
                yearColor = yearColor * (0.7 + 0.3 * (y_idx / length(uniqueYears)));
                
                h = scatter(yearData_sorted, rchgData_sorted, 80, yearColor, ...
                           'filled', 'Marker', marker, ...
                           'MarkerEdgeColor', 'black', 'LineWidth', 1.2, ...
                           'DisplayName', sprintf('Y%d', year));
                
                % Connect points with line for same year (trend visualization)
                if length(yearData_sorted) > 1
                    plot(yearData_sorted, rchgData_sorted, '--', 'Color', yearColor, ...
                         'LineWidth', 1.5, 'HandleVisibility', 'off');
                end
                
                % Calculate mean for this year (for trend line)
                yearMean = mean(rchgData_sorted);
                yearCount = length(rchgData_sorted);
                yearMeans(end+1) = yearMean;
                yearMeanYears(end+1) = year;
                yearCounts(end+1) = yearCount;
                
                totalPointsPlotted = totalPointsPlotted + length(yearData);
                
                legendEntries{end+1} = sprintf('Y%d', year);
                legendHandles(end+1) = h;
            end
        end
        
        % Plot trend line connecting year-wise means (점선으로)
        if length(yearMeans) >= 2
            [yearMeanYears_sorted, sortIdx] = sort(yearMeanYears);
            yearMeans_sorted = yearMeans(sortIdx);
            plot(yearMeanYears_sorted, yearMeans_sorted, '--', 'Color', groupColor, ...
                 'LineWidth', 2.5, 'HandleVisibility', 'off');
        end
        
        % Display mean values as text (각 연도별 평균값 및 개수 텍스트 표시)
        for mean_idx = 1:length(yearMeans)
            year_val = yearMeanYears(mean_idx);
            mean_val = yearMeans(mean_idx);
            count_val = yearCounts(mean_idx);
            text(year_val, mean_val, sprintf('%.3f\nn=%d', mean_val, count_val), ...
                 'FontSize', 7, 'Color', groupColor, ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', ...
                 'FontWeight', 'bold', ...
                 'BackgroundColor', 'white', ...
                 'EdgeColor', groupColor, ...
                 'Margin', 1);
        end
        
        % Subplot labels and title (소수점과 띄어쓰기 포함하여 표시)
        xlabel('Year', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel(sprintf('%s (mΩ)', timeLabel), 'FontSize', 10, 'FontWeight', 'bold');
        title(groupName_display, 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        
        % Add legend for years (레퍼런스: Channel legend)
        if ~isempty(legendHandles)
            legend(legendHandles, legendEntries, 'Location', 'best', 'FontSize', 8, 'NumColumns', 2);
        end
    end
    
    fprintf('Total points plotted: %d\n', totalPointsPlotted);
    
    if totalPointsPlotted == 0
        fprintf('WARNING: No data points plotted for %s!\n', timeLabel);
        close(fig);
        continue;
    end
    
    % Set common axis limits for all subplots
    allYearData = summaryTable.Year(validData);
    allRchgData = summaryTable.(timeInterval)(validData);
    xlim_range = [min(allYearData) - 0.5, max(allYearData) + 0.5];  % Year는 정수이므로 0.5씩 여유
    ylim_range = [0, 2.5];  % Y축 범위 고정: 0 ~ 2.5 mΩ
    
    for g_idx = 1:min(nGroups, nRows*nCols)
        subplot(nRows, nCols, g_idx);
        xlim(xlim_range);
        ylim(ylim_range);
        set(gca, 'XTick', uniqueYears);  % 연도만 표시
    end
    
    % Overall title
    sgtitle(sprintf('Year vs %s (Group별)', timeLabel), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Save figure (Group별 subplot)
    savePath = fullfile(figuresDir, sprintf('Year_vs_%s_Group별.fig', timeInterval));
    saveas(fig, savePath);
    fprintf('Saved: %s\n', savePath);
    close(fig);
    
    %% Create integrated plot (all groups and years in one graph)
    fig_integrated = figure('Name', sprintf('Year vs %s (통합)', timeLabel), ...
                            'Position', [100 + (fig_idx-1)*50 + 100, 100 + (fig_idx-1)*50, 1200, 800], ...
                            'Visible', 'on');
    
    hold on;
    
    % Plot data for each Group and Year combination
    legendEntries = {};
    legendHandles = [];
    totalPointsPlotted_integrated = 0;
    
    for g_idx = 1:length(selectedGroups)
        groupName = selectedGroups{g_idx};
        groupColor = groupColors(g_idx, :);
        
        % Convert group name for display: _ to space (keep decimal points)
        groupName_display = strrep(groupName, '_', ' ');
        
        % Calculate year-wise means for trend line
        yearMeans_integrated = [];
        yearMeanYears_integrated = [];
        yearCounts_integrated = [];  % 개수 저장용
        
        for y_idx = 1:length(uniqueYears)
            year = uniqueYears(y_idx);
            marker = yearMarkers{y_idx};
            
            % Filter data for this Group and Year
            mask = strcmp(summaryTable.Group, groupName) & ...
                   summaryTable.Year == year & ...
                   ~isnan(summaryTable.(timeInterval)) & ...
                   ~isnan(summaryTable.Year);
            
            if sum(mask) > 0
                yearData = summaryTable.Year(mask);
                rchgData = summaryTable.(timeInterval)(mask);
                
                % Sort by Year for line connection
                [yearData_sorted, sortIdx] = sort(yearData);
                rchgData_sorted = rchgData(sortIdx);
                
                % Plot scatter points with group color (레퍼런스: DC color)
                h = scatter(yearData_sorted, rchgData_sorted, 80, groupColor, ...
                           'filled', 'Marker', marker, ...
                           'MarkerEdgeColor', 'black', 'LineWidth', 1.2, ...
                           'DisplayName', sprintf('%s-Y%d', groupName_display, year));
                
                % Connect points with line for same Group-Year combination
                if length(yearData_sorted) > 1
                    plot(yearData_sorted, rchgData_sorted, '--', 'Color', groupColor, ...
                         'LineWidth', 1.0, 'HandleVisibility', 'off');
                end
                
                % Calculate mean for this year (for trend line)
                yearMean = mean(rchgData_sorted);
                yearCount = length(rchgData_sorted);
                yearMeans_integrated(end+1) = yearMean;
                yearMeanYears_integrated(end+1) = year;
                yearCounts_integrated(end+1) = yearCount;
                
                totalPointsPlotted_integrated = totalPointsPlotted_integrated + length(yearData);
                
                % Add to legend (only once per Group, 레퍼런스: DC Profile만 legend에)
                if y_idx == 1
                    legendEntries{end+1} = groupName_display;
                    legendHandles(end+1) = h;
                end
            end
        end
        
        % Plot trend line connecting year-wise means for this group (점선으로)
        if length(yearMeans_integrated) >= 2
            [yearMeanYears_sorted, sortIdx] = sort(yearMeanYears_integrated);
            yearMeans_sorted = yearMeans_integrated(sortIdx);
            plot(yearMeanYears_sorted, yearMeans_sorted, '--', 'Color', groupColor, ...
                 'LineWidth', 2.5, 'HandleVisibility', 'off');
        end
        
        % Display mean values as text (각 연도별 평균값 및 개수 텍스트 표시)
        for mean_idx = 1:length(yearMeans_integrated)
            year_val = yearMeanYears_integrated(mean_idx);
            mean_val = yearMeans_integrated(mean_idx);
            count_val = yearCounts_integrated(mean_idx);
            text(year_val, mean_val, sprintf('%.3f\nn=%d', mean_val, count_val), ...
                 'FontSize', 7, 'Color', groupColor, ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', ...
                 'FontWeight', 'bold', ...
                 'BackgroundColor', 'white', ...
                 'EdgeColor', groupColor, ...
                 'Margin', 1);
        end
    end
    
    if totalPointsPlotted_integrated > 0
        % Set axis limits
        allYearData = summaryTable.Year(validData);
        allRchgData = summaryTable.(timeInterval)(validData);
        xlim_range = [min(allYearData) - 0.5, max(allYearData) + 0.5];
        ylim_range = [0, 2.5];  % Y축 범위 고정: 0 ~ 2.5 mΩ
        
        xlim(xlim_range);
        ylim(ylim_range);
        
        % Add legend
        if ~isempty(legendHandles)
            legend(legendHandles, legendEntries, 'Location', 'best', 'FontSize', 10);
        end
        
        % Labels and title
        xlabel('Year', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(sprintf('%s (mΩ)', timeLabel), 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Year vs %s (All Groups & Years)', timeLabel), ...
              'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        
        % Add text annotation for year markers (레퍼런스: Channel marker legend)
        if length(uniqueYears) <= 8
            markerLegend = cell(length(uniqueYears), 1);
            for y_idx = 1:length(uniqueYears)
                markerLegend{y_idx} = sprintf('Y%d: %s', uniqueYears(y_idx), yearMarkers{y_idx});
            end
            text(0.02, 0.98, strjoin(markerLegend, '\n'), ...
                 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                 'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');
        end
        
        % Save integrated figure
        if isempty(targetEventType)
            savePath = fullfile(figuresDir, sprintf('Year_vs_%s_통합.fig', timeInterval));
        else
            savePath = fullfile(figuresDir, sprintf('Year_vs_%s_통합_%s.fig', timeInterval, targetEventType));
        end
        saveas(fig_integrated, savePath);
        fprintf('Saved: %s\n', savePath);
    else
        fprintf('WARNING: No data points for integrated plot of %s!\n', timeLabel);
    end
    close(fig_integrated);
    
    %% Create Box Plots (Year vs Rchg) - 이상치 제거된 데이터 사용
    fprintf('\n--- Creating Box Plots for %s ---\n', timeLabel);
    
    % Group별 박스플롯
    fig_box = figure('Name', sprintf('Year vs %s Box Plot (Group별)', timeLabel), ...
                     'Position', [100 + (fig_idx-1)*50, 100 + (fig_idx-1)*50, 1600, 900], ...
                     'Visible', 'on');
    
    for g_idx = 1:min(nGroups, nRows*nCols)
        groupName = selectedGroups{g_idx};
        groupColor = groupColors(g_idx, :);
        
        % Convert group name for display: _ to space (keep decimal points)
        groupName_display = strrep(groupName, '_', ' ');
        
        subplot(nRows, nCols, g_idx);
        hold on;
        
        % Filter data for this Group (이상치 제거된 데이터 사용)
        groupMask = strcmp(summaryTable.Group, groupName) & ...
                   ~isnan(summaryTable.(timeInterval)) & ...
                   ~isnan(summaryTable.Year);
        
        if sum(groupMask) >= 3
            year_data = summaryTable.Year(groupMask);
            rchg_data = summaryTable.(timeInterval)(groupMask);
            
            % Year별로 데이터 그룹화
            uniqueYears_group = sort(unique(year_data));
            
            % 박스플롯을 위한 데이터 준비: 모든 데이터를 하나의 벡터로, 그룹 변수 생성
            all_rchg = [];
            all_groups = [];
            yearLabels = cell(length(uniqueYears_group), 1);
            means = zeros(length(uniqueYears_group), 1);
            stds = zeros(length(uniqueYears_group), 1);
            medians = zeros(length(uniqueYears_group), 1);
            mins = zeros(length(uniqueYears_group), 1);
            maxs = zeros(length(uniqueYears_group), 1);
            counts = zeros(length(uniqueYears_group), 1);
            
            for y_idx = 1:length(uniqueYears_group)
                year = uniqueYears_group(y_idx);
                yearMask = year_data == year;
                year_rchg = rchg_data(yearMask);
                
                if ~isempty(year_rchg)
                    all_rchg = [all_rchg; year_rchg(:)];
                    all_groups = [all_groups; repmat(y_idx, length(year_rchg), 1)];
                    yearLabels{y_idx} = sprintf('%d', year);
                    means(y_idx) = mean(year_rchg);
                    stds(y_idx) = std(year_rchg);
                    medians(y_idx) = median(year_rchg);
                    mins(y_idx) = min(year_rchg);
                    maxs(y_idx) = max(year_rchg);
                    counts(y_idx) = length(year_rchg);
                end
            end
            
            % 박스플롯 생성
            if ~isempty(all_rchg) && length(unique(all_groups)) > 0
                bp = boxplot(all_rchg, all_groups, 'Labels', yearLabels, 'Colors', 'k', ...
                            'Symbol', 'k+', 'OutlierSize', 4);
                
                % 박스플롯 색상 설정 (그룹 색상 사용)
                h = findobj(bp, 'Tag', 'Box');
                for i = 1:length(h)
                    patch(get(h(i), 'XData'), get(h(i), 'YData'), groupColor, ...
                          'FaceAlpha', 0.5, 'EdgeColor', groupColor, 'LineWidth', 1.5);
                end
                
                % 평균값을 점으로 표시
                for y_idx = 1:length(uniqueYears_group)
                    if counts(y_idx) > 0
                        plot(y_idx, means(y_idx), 'ro', 'MarkerSize', 8, ...
                             'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
                    end
                end
                
                % 각 Year별 통계 정보 텍스트 표시 (박스플롯 위에)
                for y_idx = 1:length(uniqueYears_group)
                    if counts(y_idx) > 0
                        % Y축 범위를 고려하여 텍스트 위치 설정 (박스플롯 위쪽)
                        y_pos = 2.5 - 0.1;  % Y축 최대값 근처
                        
                        % 통계 정보 텍스트 생성
                        statsText_year = sprintf('n=%d\nμ=%.3f\nσ=%.3f\nmed=%.3f', ...
                            counts(y_idx), means(y_idx), stds(y_idx), medians(y_idx));
                        
                        text(y_idx, y_pos, statsText_year, ...
                             'FontSize', 7, 'Color', 'black', ...
                             'HorizontalAlignment', 'center', ...
                             'VerticalAlignment', 'top', ...
                             'FontWeight', 'bold', ...
                             'BackgroundColor', 'white', ...
                             'EdgeColor', 'black', ...
                             'Margin', 2);
                    end
                end
                
                % 통계 정보 텍스트 추가 (Year vs Rchg 상관계수)
                validMaskYear = ~isnan(summaryTable.Year);
                validMaskBoth = groupMask & validMaskYear;
                if sum(validMaskBoth) >= 3
                    year_full = summaryTable.Year(validMaskBoth);
                    rchg_full = summaryTable.(timeInterval)(validMaskBoth);
                    [R, P] = corrcoef(year_full, rchg_full);
                    r_val = R(1, 2);
                    p_val = P(1, 2);
                    
                    statsText = sprintf('R = %.4f\np = %.4e\nn = %d', ...
                                      r_val, p_val, length(rchg_full));
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
                         'VerticalAlignment', 'top', 'FontSize', 9, ...
                         'BackgroundColor', 'white', 'EdgeColor', 'black');
                end
                
                xlabel('Year', 'FontSize', 10, 'FontWeight', 'bold');
                ylabel(sprintf('%s (mΩ)', timeLabel), 'FontSize', 10, 'FontWeight', 'bold');
                title(groupName_display, 'FontSize', 11, 'FontWeight', 'bold');
                grid on;
                
                % Y축 범위 고정
                ylim([0, 2.5]);
            else
                text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'FontSize', 14);
                title(groupName_display, 'FontSize', 11, 'FontWeight', 'bold');
            end
        else
            text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', 'FontSize', 14);
            title(groupName_display, 'FontSize', 11, 'FontWeight', 'bold');
        end
    end
    
    % Overall title
    sgtitle(sprintf('Year vs %s Box Plot (Group별)', timeLabel), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Save figure
    savePath_box = fullfile(figuresDir, sprintf('Year_vs_%s_BoxPlot_Group별.fig', timeInterval));
    saveas(fig_box, savePath_box);
    fprintf('Saved: %s\n', savePath_box);
    close(fig_box);
    
    %% Create integrated box plot (all groups in one graph)
    fig_box_integrated = figure('Name', sprintf('Year vs %s Box Plot (통합)', timeLabel), ...
                                'Position', [100 + (fig_idx-1)*50 + 100, 100 + (fig_idx-1)*50, 1200, 800], ...
                                'Visible', 'on');
    
    hold on;
    
    % 모든 그룹의 데이터를 Year별로 그룹화
    validMask_all = ~isnan(summaryTable.(timeInterval)) & ~isnan(summaryTable.Year);
    if sum(validMask_all) >= 3
        year_data_all = summaryTable.Year(validMask_all);
        rchg_data_all = summaryTable.(timeInterval)(validMask_all);
        
        % Year별로 데이터 그룹화
        uniqueYears_all = sort(unique(year_data_all));
        
        % 박스플롯을 위한 데이터 준비: 모든 데이터를 하나의 벡터로, 그룹 변수 생성
        all_rchg_all = [];
        all_groups_all = [];
        yearLabels_all = cell(length(uniqueYears_all), 1);
        means_all = zeros(length(uniqueYears_all), 1);
        stds_all = zeros(length(uniqueYears_all), 1);
        medians_all = zeros(length(uniqueYears_all), 1);
        mins_all = zeros(length(uniqueYears_all), 1);
        maxs_all = zeros(length(uniqueYears_all), 1);
        counts_all = zeros(length(uniqueYears_all), 1);
        
        for y_idx = 1:length(uniqueYears_all)
            year = uniqueYears_all(y_idx);
            yearMask = year_data_all == year;
            year_rchg = rchg_data_all(yearMask);
            
            if ~isempty(year_rchg)
                all_rchg_all = [all_rchg_all; year_rchg(:)];
                all_groups_all = [all_groups_all; repmat(y_idx, length(year_rchg), 1)];
                yearLabels_all{y_idx} = sprintf('%d', year);
                means_all(y_idx) = mean(year_rchg);
                stds_all(y_idx) = std(year_rchg);
                medians_all(y_idx) = median(year_rchg);
                mins_all(y_idx) = min(year_rchg);
                maxs_all(y_idx) = max(year_rchg);
                counts_all(y_idx) = length(year_rchg);
            end
        end
        
        % 박스플롯 생성
        if ~isempty(all_rchg_all) && length(unique(all_groups_all)) > 0
            bp = boxplot(all_rchg_all, all_groups_all, 'Labels', yearLabels_all, 'Colors', 'k', ...
                        'Symbol', 'k+', 'OutlierSize', 4);
            
            % 각 Year별 평균값 계산 및 표시
            % 평균값을 점으로 표시
            for y_idx = 1:length(uniqueYears_all)
                if counts_all(y_idx) > 0
                    plot(y_idx, means_all(y_idx), 'ro', 'MarkerSize', 10, ...
                         'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
                end
            end
            
            % 각 Year별 통계 정보 텍스트 표시 (박스플롯 위에)
            for y_idx = 1:length(uniqueYears_all)
                if counts_all(y_idx) > 0
                    % Y축 범위를 고려하여 텍스트 위치 설정 (박스플롯 위쪽)
                    y_pos = 2.5 - 0.1;  % Y축 최대값 근처
                    
                    % 통계 정보 텍스트 생성
                    statsText_year = sprintf('n=%d\nμ=%.3f\nσ=%.3f\nmed=%.3f', ...
                        counts_all(y_idx), means_all(y_idx), stds_all(y_idx), medians_all(y_idx));
                    
                    text(y_idx, y_pos, statsText_year, ...
                         'FontSize', 8, 'Color', 'black', ...
                         'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'top', ...
                         'FontWeight', 'bold', ...
                         'BackgroundColor', 'white', ...
                         'EdgeColor', 'black', ...
                         'Margin', 2);
                end
            end
            
            % 통계 정보 텍스트 추가
            [R, P] = corrcoef(year_data_all, rchg_data_all);
            r_val = R(1, 2);
            p_val = P(1, 2);
            
            statsText = sprintf('R = %.4f\np = %.4e\nn = %d', ...
                              r_val, p_val, length(rchg_data_all));
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
            
            xlabel('Year', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel(sprintf('%s (mΩ)', timeLabel), 'FontSize', 12, 'FontWeight', 'bold');
            title(sprintf('Year vs %s Box Plot (All Groups & Years)', timeLabel), ...
                  'FontSize', 14, 'FontWeight', 'bold');
            grid on;
            
            % Y축 범위 고정
            ylim([0, 2.5]);
        else
            text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', 'FontSize', 14);
        end
    else
        text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', 14);
    end
    
    % Save integrated box plot
    if isempty(targetEventType)
        savePath_box_integrated = fullfile(figuresDir, sprintf('Year_vs_%s_BoxPlot_통합.fig', timeInterval));
    else
        savePath_box_integrated = fullfile(figuresDir, sprintf('Year_vs_%s_BoxPlot_통합_%s.fig', timeInterval, targetEventType));
    end
    saveas(fig_box_integrated, savePath_box_integrated);
    fprintf('Saved: %s\n', savePath_box_integrated);
    close(fig_box_integrated);
end

fprintf('\n=== Visualization Complete ===\n');
fprintf('All figures saved to figures/Integrated_SOH_Rchg/\n');
fprintf('Total figures created: %d\n', length(timeIntervals) * 4);  % Group별 산점도 + 통합 산점도 + Group별 박스플롯 + 통합 박스플롯
