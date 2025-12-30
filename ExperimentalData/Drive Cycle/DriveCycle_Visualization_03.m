%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 03_DriveCycle_Visualization.m
% 통합 시각화
% 
% 목적: 
% - DriveCycle_Summary_Table.mat를 로드
% - Capacity vs Resistance 통합 그래프 생성 (6개: R_1s, R_3s, R_5s, R_10s, R_30s, R_60s)
% - 충전 및 방전 저항 모두 시각화
% - Color: DC Profile (DC1~DC8)
% - Marker: Channel (Ch9~Ch16)
%
% 입력:
% - DriveCycle_Summary_Table.mat (02번 스크립트 출력)
%
% 출력:
% - figures/Integrated_Capacity_Rchg/Capacity_vs_R_*.fig (3개)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('=== Drive Cycle DCIR Visualization ===\n');

%% Configuration - Paths
inputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
outputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Configuration - Event Type Selection
% Options: 'Charge', 'Discharge', or 'All'
% Note: 'All' will create separate visualizations for Charge and Discharge
targetEventType = 'All';  % 'All': both charge and discharge (separate plots), 'Charge': only charge, 'Discharge': only discharge

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
% Check for Capacity_C3 and Capacity_OCV
if ismember('Capacity_C3', summaryTable.Properties.VariableNames)
    fprintf('Data range - Capacity_C3: [%.2f, %.2f] Ah\n', min(summaryTable.Capacity_C3(~isnan(summaryTable.Capacity_C3))), max(summaryTable.Capacity_C3(~isnan(summaryTable.Capacity_C3))));
end
if ismember('Capacity_OCV', summaryTable.Properties.VariableNames)
    fprintf('Data range - Capacity_OCV: [%.2f, %.2f] Ah\n', min(summaryTable.Capacity_OCV(~isnan(summaryTable.Capacity_OCV))), max(summaryTable.Capacity_OCV(~isnan(summaryTable.Capacity_OCV))));
end
% Legacy support for 'Capacity' column
if ismember('Capacity', summaryTable.Properties.VariableNames)
    fprintf('Data range - Capacity (legacy): [%.2f, %.2f] Ah\n', min(summaryTable.Capacity(~isnan(summaryTable.Capacity))), max(summaryTable.Capacity(~isnan(summaryTable.Capacity))));
end

% Store original summary table for 'All' case
originalSummaryTable = summaryTable;

% Determine which event types to visualize
if ismember('EventType', summaryTable.Properties.VariableNames)
    if strcmp(targetEventType, 'Charge')
        summaryTable = summaryTable(strcmp(summaryTable.EventType, 'Charge'), :);
        eventTypesToPlot = {'Charge'};
        fprintf('\nFiltered to Charge events only: %d rows\n', height(summaryTable));
    elseif strcmp(targetEventType, 'Discharge')
        summaryTable = summaryTable(strcmp(summaryTable.EventType, 'Discharge'), :);
        eventTypesToPlot = {'Discharge'};
        fprintf('\nFiltered to Discharge events only: %d rows\n', height(summaryTable));
    else
        % 'All': visualize both Charge and Discharge separately
        eventTypesToPlot = {'Charge', 'Discharge'};
        fprintf('\nWill visualize both Charge and Discharge events separately\n');
        fprintf('Total rows: %d\n', height(summaryTable));
    end
else
    eventTypesToPlot = {'All'};  % No EventType column
    fprintf('\nNo EventType column found. Visualizing all data.\n');
end

%% Create figures directory
figuresDir = fullfile(outputDir, 'figures', 'Integrated_Capacity_Rchg');
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

%% Define time intervals to plot (all 6 intervals)
timeIntervals = {'R_1s', 'R_3s', 'R_5s', 'R_10s', 'R_30s', 'R_60s'};
timeIntervalLabels = {'R 1s', 'R 3s', 'R 5s', 'R 10s', 'R 30s', 'R 60s'};

%% Get unique DC Profiles and Channels for color/marker mapping
uniqueDCProfiles = unique(summaryTable.DC_Profile);
uniqueChannels = unique(summaryTable.Channel);

% Ensure DC1~DC8 are always included (even if no data)
allDCProfiles = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};
for d = 1:length(allDCProfiles)
    if ~ismember(allDCProfiles{d}, uniqueDCProfiles)
        uniqueDCProfiles{end+1} = allDCProfiles{d};
    end
end

% Sort DC profiles (DC1, DC2, ...)
dc_nums = zeros(length(uniqueDCProfiles), 1);
for d = 1:length(uniqueDCProfiles)
    dc_str = uniqueDCProfiles{d};
    num_match = regexp(dc_str, 'DC(\d+)', 'tokens');
    if ~isempty(num_match)
        dc_nums(d) = str2double(num_match{1}{1});
    else
        dc_nums(d) = 999;
    end
end
[~, dc_sort_idx] = sort(dc_nums);
sortedDCProfiles = uniqueDCProfiles(dc_sort_idx);

% Sort channels
sortedChannels = sort(uniqueChannels);

%% Define color map for DC Profiles (8 colors)
dcColors = lines(length(sortedDCProfiles));

%% Define marker styles for Channels (8 different markers)
channelMarkers = {'o', 's', '^', 'd', 'v', '>', '<', 'p'};
if length(sortedChannels) > length(channelMarkers)
    % If more channels than markers, cycle through markers
    channelMarkers = repmat(channelMarkers, 1, ceil(length(sortedChannels) / length(channelMarkers)));
end

%% Create figures (one for each time interval and event type)
% Use subplot grid: DC Profile별로 subplot, Channel은 마커로 구분
% Loop through event types first, then time intervals
for eventTypeIdx = 1:length(eventTypesToPlot)
    currentEventType = eventTypesToPlot{eventTypeIdx};
    
    % Filter data for current event type
    if strcmp(currentEventType, 'All')
        currentSummaryTable = originalSummaryTable;
        eventTypeLabel = '';
        eventTypeSuffix = '';
    else
        if ismember('EventType', originalSummaryTable.Properties.VariableNames)
            currentSummaryTable = originalSummaryTable(strcmp(originalSummaryTable.EventType, currentEventType), :);
            eventTypeLabel = sprintf(' (%s)', currentEventType);
            eventTypeSuffix = sprintf('_%s', currentEventType);
        else
            currentSummaryTable = originalSummaryTable;
            eventTypeLabel = '';
            eventTypeSuffix = '';
        end
    end
    
    fprintf('\n=== Visualizing %s Events ===\n', currentEventType);
    
    % Loop through all time intervals
    for fig_idx = 1:length(timeIntervals)
        timeInterval = timeIntervals{fig_idx};
        timeLabel = timeIntervalLabels{fig_idx};
        
        % Check if column exists
        if ~ismember(timeInterval, currentSummaryTable.Properties.VariableNames)
            fprintf('WARNING: Column %s not found in summary table. Skipping.\n', timeInterval);
            continue;
        end
    
    % Determine which capacity column to use
    % Use Capacity_C3 if available, otherwise fall back to Capacity_OCV or Capacity
    capacityCol = 'Capacity_C3';
    if ~ismember('Capacity_C3', summaryTable.Properties.VariableNames)
        if ismember('Capacity_OCV', summaryTable.Properties.VariableNames)
            capacityCol = 'Capacity_OCV';
        elseif ismember('Capacity', summaryTable.Properties.VariableNames)
            capacityCol = 'Capacity';
        else
            fprintf('WARNING: No capacity column found. Skipping plot.\n');
            continue;
        end
    end
    
        % Debug: Check data availability for this time interval
        validData = ~isnan(currentSummaryTable.(timeInterval)) & ~isnan(currentSummaryTable.(capacityCol));
        fprintf('\n--- Creating plot for %s%s ---\n', timeLabel, eventTypeLabel);
        fprintf('Using capacity column: %s\n', capacityCol);
        fprintf('Valid data points: %d / %d\n', sum(validData), height(currentSummaryTable));
        
        if sum(validData) == 0
            fprintf('WARNING: No valid data for %s%s. Skipping plot.\n', timeLabel, eventTypeLabel);
            continue;
        end
    
        % Create figure with subplots (2 rows x 4 columns for 8 DC profiles)
        fig = figure('Name', sprintf('Capacity vs %s (DC별)%s', timeLabel, eventTypeLabel), ...
                     'Position', [100 + (fig_idx-1)*50 + eventTypeIdx*20, 100 + (fig_idx-1)*50, 1600, 900], ...
                     'Visible', 'on');
        
        % Create 2x4 subplot grid
        nRows = 2;
        nCols = 4;
        
        totalPointsPlotted = 0;
        
        % Get unique DC profiles and channels from current summary table
        if height(currentSummaryTable) > 0
            currentDCProfiles = unique(currentSummaryTable.DC_Profile);
            currentChannels = unique(currentSummaryTable.Channel);
        else
            currentDCProfiles = sortedDCProfiles;
            currentChannels = sortedChannels;
        end
        
        % Plot each DC Profile in a separate subplot
        for dc_idx = 1:length(sortedDCProfiles)
            dcProfile = sortedDCProfiles{dc_idx};
            
            subplot(nRows, nCols, dc_idx);
            hold on;
            
            legendEntries = {};
            legendHandles = [];
            
            % Plot all channels for this DC profile
            for ch_idx = 1:length(sortedChannels)
                channel = sortedChannels(ch_idx);
                marker = channelMarkers{ch_idx};
                
                % Filter data for this DC Profile and Channel
                mask = strcmp(currentSummaryTable.DC_Profile, dcProfile) & ...
                       currentSummaryTable.Channel == channel & ...
                       ~isnan(currentSummaryTable.(timeInterval)) & ...
                       ~isnan(currentSummaryTable.(capacityCol));
            
            if sum(mask) > 0
                capData = summaryTable.(capacityCol)(mask);
                rchgData = summaryTable.(timeInterval)(mask);
                
                    % Sort by capacity for line connection
                    [capData_sorted, sortIdx] = sort(capData);
                    rchgData_sorted = rchgData(sortIdx);
                    
                    % Plot scatter points with channel-specific color
                    chColor = dcColors(dc_idx, :);  % Use DC color as base
                    % Make channel colors slightly different by adjusting brightness
                    chColor = chColor * (0.7 + 0.3 * (ch_idx / length(sortedChannels)));
                    
                    h = scatter(capData_sorted, rchgData_sorted, 80, chColor, ...
                               'filled', 'Marker', marker, ...
                               'MarkerEdgeColor', 'black', 'LineWidth', 1.2, ...
                               'DisplayName', sprintf('Ch%d', channel));
                    
                    % Connect points with line for same channel (trend visualization)
                    if length(capData_sorted) > 1
                        plot(capData_sorted, rchgData_sorted, '--', 'Color', chColor, ...
                             'LineWidth', 1.5, 'HandleVisibility', 'off');
                    end
                    
                    totalPointsPlotted = totalPointsPlotted + length(capData);
                    
                    legendEntries{end+1} = sprintf('Ch%d', channel);
                    legendHandles(end+1) = h;
                end
            end
        
        % Subplot labels and title
        xlabel('Capacity (Ah)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel(sprintf('%s (mΩ)', timeLabel), 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('%s%s', dcProfile, eventTypeLabel), 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        
        % Add legend for channels
        if ~isempty(legendHandles)
            legend(legendHandles, legendEntries, 'Location', 'best', 'FontSize', 8, 'NumColumns', 2);
        else
            % No data for this DC profile - display "No Data" message
            text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 14, ...
                'Color', 'red', 'Units', 'normalized');
        end
        
        % Set consistent axis limits across all subplots (will be set after all plots)
    end
    
    fprintf('Total points plotted: %d\n', totalPointsPlotted);
    
    % Always create and save figure, even if no data points
    % (This ensures all DC profiles are shown, even if empty)
    
        % Set common axis limits for all subplots (only if valid data exists)
        if sum(validData) > 0
            allCapData = currentSummaryTable.(capacityCol)(validData);
            allRchgData = currentSummaryTable.(timeInterval)(validData);
        xlim_range = [min(allCapData) - 0.05*range(allCapData), max(allCapData) + 0.05*range(allCapData)];
        ylim_range = [min(allRchgData) - 0.05*range(allRchgData), max(allRchgData) + 0.05*range(allRchgData)];
        
        % Reverse x-axis so initial capacity (higher value) is on the left
        % This makes degradation trend more intuitive (left to right = degradation)
        for dc_idx = 1:length(sortedDCProfiles)
            subplot(nRows, nCols, dc_idx);
            xlim(xlim_range);  % Set normal order [min, max]
            ylim(ylim_range);
            set(gca, 'XDir', 'reverse');  % Reverse x-axis direction (high value on left)
        end
    else
        % No valid data - set default axis limits
        for dc_idx = 1:length(sortedDCProfiles)
            subplot(nRows, nCols, dc_idx);
            xlim([0, 100]);  % Default capacity range
            ylim([0, 100]);  % Default resistance range
        end
    end
    
    % Overall title
    sgtitle(sprintf('Capacity vs %s (DC Profile별)%s', timeLabel, eventTypeLabel), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
        % Save figure (DC별 subplot) - Always save, even if no data
        savePath = fullfile(figuresDir, sprintf('Capacity_vs_%s_DC별%s.fig', timeInterval, eventTypeSuffix));
        saveas(fig, savePath);
        fprintf('Saved: %s (DC별, %s events)\n', savePath, currentEventType);
        close(fig);
        
        %% Create integrated plot (all DC profiles and channels in one graph)
        fig_integrated = figure('Name', sprintf('Capacity vs %s (통합)%s', timeLabel, eventTypeLabel), ...
                                'Position', [100 + (fig_idx-1)*50 + 100 + eventTypeIdx*20, 100 + (fig_idx-1)*50, 1200, 800], ...
                                'Visible', 'on');
        
        hold on;
        
        % Plot data for each DC Profile and Channel combination
        legendEntries = {};
        legendHandles = [];
        totalPointsPlotted_integrated = 0;
        
        for dc_idx = 1:length(sortedDCProfiles)
            dcProfile = sortedDCProfiles{dc_idx};
            dcColor = dcColors(dc_idx, :);
            
            for ch_idx = 1:length(sortedChannels)
                channel = sortedChannels(ch_idx);
                marker = channelMarkers{ch_idx};
                
                % Filter data for this DC Profile and Channel
                mask = strcmp(currentSummaryTable.DC_Profile, dcProfile) & ...
                       currentSummaryTable.Channel == channel & ...
                       ~isnan(currentSummaryTable.(timeInterval)) & ...
                       ~isnan(currentSummaryTable.(capacityCol));
                
                if sum(mask) > 0
                    capData = currentSummaryTable.(capacityCol)(mask);
                    rchgData = currentSummaryTable.(timeInterval)(mask);
                
                % Sort by capacity for line connection
                [capData_sorted, sortIdx] = sort(capData);
                rchgData_sorted = rchgData(sortIdx);
                
                % Plot scatter points with DC color
                h = scatter(capData_sorted, rchgData_sorted, 80, dcColor, ...
                           'filled', 'Marker', marker, ...
                           'MarkerEdgeColor', 'black', 'LineWidth', 1.2, ...
                           'DisplayName', sprintf('%s-Ch%d', dcProfile, channel));
                
                % Connect points with line for same DC-Channel combination
                if length(capData_sorted) > 1
                    plot(capData_sorted, rchgData_sorted, '--', 'Color', dcColor, ...
                         'LineWidth', 1.0, 'HandleVisibility', 'off');
                end
                
                totalPointsPlotted_integrated = totalPointsPlotted_integrated + length(capData);
                
                % Add to legend (only once per DC Profile)
                if ch_idx == 1
                    legendEntries{end+1} = dcProfile;
                    legendHandles(end+1) = h;
                end
            end
        end
    end
    
        if totalPointsPlotted_integrated > 0
            % Set axis limits
            allCapData = currentSummaryTable.(capacityCol)(validData);
            allRchgData = currentSummaryTable.(timeInterval)(validData);
            xlim_range = [min(allCapData) - 0.05*range(allCapData), max(allCapData) + 0.05*range(allCapData)];
            ylim_range = [min(allRchgData) - 0.05*range(allRchgData), max(allRchgData) + 0.05*range(allRchgData)];
            
            xlim(xlim_range);
            ylim(ylim_range);
            set(gca, 'XDir', 'reverse');  % Reverse x-axis (initial capacity on left)
            
            % Add legend
            if ~isempty(legendHandles)
                legend(legendHandles, legendEntries, 'Location', 'best', 'FontSize', 10);
            end
        else
            % No data - set default limits and show message
            xlim([0, 100]);
            ylim([0, 100]);
            text(0.5, 0.5, sprintf('No data for %s%s', timeLabel, eventTypeLabel), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 14, 'Color', 'red', 'Units', 'normalized');
        end
        
        % Labels and title
        xlabel('Capacity (Ah)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(sprintf('%s (mΩ)', timeLabel), 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Capacity vs %s (All Channels & DC Profiles)%s', timeLabel, eventTypeLabel), ...
              'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        
        % Add text annotation for channel markers (only if data exists)
        if totalPointsPlotted_integrated > 0 && length(sortedChannels) <= 8
            markerLegend = cell(length(sortedChannels), 1);
            for ch_idx = 1:length(sortedChannels)
                markerLegend{ch_idx} = sprintf('Ch%d: %s', sortedChannels(ch_idx), channelMarkers{ch_idx});
            end
            text(0.02, 0.98, strjoin(markerLegend, '\n'), ...
                 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                 'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');
        end
        
        % Save integrated figure - Always save, even if no data
        savePath = fullfile(figuresDir, sprintf('Capacity_vs_%s_통합%s.fig', timeInterval, eventTypeSuffix));
        saveas(fig_integrated, savePath);
        if totalPointsPlotted_integrated > 0
            fprintf('Saved: %s (통합, %s events, %d points)\n', savePath, currentEventType, totalPointsPlotted_integrated);
        else
            fprintf('Saved: %s (통합, %s events, NO DATA)\n', savePath, currentEventType);
        end
        close(fig_integrated);
    end  % End of timeIntervals loop
end  % End of eventTypesToPlot loop

fprintf('\n=== Visualization Complete ===\n');
fprintf('All figures saved to figures/Integrated_Capacity_Rchg/\n');
fprintf('Total time intervals: %d\n', length(timeIntervals));
fprintf('Total event types: %d\n', length(eventTypesToPlot));
fprintf('Total figures created: %d (per event type)\n', length(timeIntervals) * length(eventTypesToPlot));

