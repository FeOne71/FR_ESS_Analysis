%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experimental Data Capacity Retention Summary Code - Dynamic Version
% Required .csv files: RPT / Aging 
% Updated date: 2025-01-27 
% Features: Dynamic RPT and Aging data loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% =========================================================================
% Configuration
% =========================================================================
% Folder directory
dataDir   = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\Experimental Data';
rptFolder    = fullfile(dataDir, 'RPT');
agingFolder  = fullfile(dataDir, 'Aging');
saveFolder   = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Capacity_Trend_Figures';     

% Channel (can be single channel or array of channels)
chList = {'09', '10', '11', '12', '13', '14', '15', '16'};  % Modify as needed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define RPT cycles to analyze (can be easily modified)
rptCycles = [0, 200, 400, 600, 800];  % Add more cycles here as needed

%% Define Aging cycles to analyze (can be easily modified)
agingCycles = [0, 200, 400, 600, 800];  % Start and end cycles for each aging period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fig color
color_static = '#CD534C';
color_ocv    = '#0073C2';
color_aging  = '#20854E';

% =========================================================================
% Dynamic Data Loading Functions
% =========================================================================

% Function to load RPT data
function [Q_static, Q_ocv, success] = loadRPTData(rptFolder, ch, cycle)
    rptFile = fullfile(rptFolder, sprintf('Ch%s_RPT_%dcyc.csv', ch, cycle));
    if ~isfile(rptFile)
        Q_static = NaN;
        Q_ocv = NaN;
        success = false;
        return;
    end
    T_rpt = readmatrix(rptFile);
    
    % Static Capacity (StepIdx=3, CycleIdx=2)
    mask_static = (T_rpt(:,2) == 3) & (T_rpt(:,4) == 2);
    Q_static = T_rpt(mask_static, 9);
    Q_static = Q_static(end);
    
    % OCV Capacity (StepIdx=10, CycleIdx=2)
    mask_ocv = (T_rpt(:,2) == 10) & (T_rpt(:,4) == 2);
    Q_ocv = T_rpt(mask_ocv, 9);
    Q_ocv = Q_ocv(end);
    success = true;
end

% Function to load Aging data
function [Q_aging, maxCycle, actualMaxCycle, success] = loadAgingData(agingFolder, ch, startCycle, endCycle)
    agingSubFolder = fullfile(agingFolder, sprintf('%dto%dcyc', startCycle, endCycle));
    agingPattern = sprintf('Ch%s*%dto%dcyc*.csv', ch, startCycle, endCycle);
    agingFiles = dir(fullfile(agingSubFolder, agingPattern));
    
    if isempty(agingFiles)
        Q_aging = [];
        maxCycle = 0;
        actualMaxCycle = 0;
        success = false;
        return;
    end
    
    agingFile = fullfile(agingSubFolder, agingFiles(1).name);
    T_aging = readmatrix(agingFile);
    
    % Find actual cycles in data
    cycles = unique(T_aging(:,4));
    cycles = cycles(cycles >= 1);
    actualMaxCycle = max(cycles);  % Actual cycle number in data (absolute)
    actualMinCycle = min(cycles);  % Actual minimum cycle number in data
    
    % Extract capacity for each cycle (using actual cycle numbers from data)
    % Store capacities in a map first
    Q_aging_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for i = 1:length(cycles)
        cycleNum = cycles(i);
        mask = (T_aging(:,2) == 3) & (T_aging(:,4) == cycleNum);
        Q_tmp = T_aging(mask, 9);
        if ~isempty(Q_tmp)
            Q_aging_map(cycleNum) = Q_tmp(end);
        end
    end
    
    % Create capacity array: use actual cycle numbers as indices
    % For plotting, we'll use a simple array indexed by cycle number
    maxCycle = actualMaxCycle;  % Use actual max cycle number
    
    % Create array with capacity for each cycle
    Q_aging = NaN(1, actualMaxCycle);
    for cyc = actualMinCycle:actualMaxCycle
        if isKey(Q_aging_map, cyc)
            Q_aging(cyc) = Q_aging_map(cyc);
        end
    end
    
    % Trim to only valid portion (remove leading NaNs)
    firstValidIdx = find(~isnan(Q_aging), 1);
    if ~isempty(firstValidIdx)
        Q_aging = Q_aging(firstValidIdx:end);
    else
        Q_aging = [];
    end
    success = true;
end

% =========================================================================
% Process Each Channel
% =========================================================================

% Initialize structure to store capacity data for all channels
allChannelsCapacity = struct();

for chIdx = 1:length(chList)
    ch = chList{chIdx};
    fprintf('Processing Channel %s...\n', ch);
    
    % Initialize structure for this channel
    channelCapacity = struct();
    channelCapacity.cycles = [];
    channelCapacity.capacity = [];  % Static Capacity
    channelCapacity.capacity_ocv = [];  % OCV Capacity
    
% =========================================================================
% Load All Data Dynamically
% =========================================================================

% Load RPT data
rptData = struct();
for i = 1:length(rptCycles)
    cycle = rptCycles(i);
    [Q_static, Q_ocv, success] = loadRPTData(rptFolder, ch, cycle);
    if success
        rptData.(sprintf('cycle_%d', cycle)) = struct('static', Q_static, 'ocv', Q_ocv);
    else
        fprintf('Warning: RPT%d 파일 없음 - 건너뜀\n', cycle);
    end
end

% Load Aging data
agingData = struct();
for i = 1:length(agingCycles)-1
    startCycle = agingCycles(i);
    endCycle = agingCycles(i+1);
    [Q_aging, maxCycle, actualMaxCycle, success] = loadAgingData(agingFolder, ch, startCycle, endCycle);
    if success
        agingData.(sprintf('period_%dto%d', startCycle, endCycle)) = struct('capacity', Q_aging, 'maxCycle', maxCycle, 'actualMaxCycle', actualMaxCycle);
    else
        fprintf('Warning: Aging %dto%dcyc 파일 없음 - 건너뜀\n', startCycle, endCycle);
    end
end

% =========================================================================
% Generate Plot Data
% =========================================================================

% Calculate X-axis positions
x_pos = 1;
x_vals = [];
y_vals = [];
plotInfo = struct();

% Add RPT data points
for i = 1:length(rptCycles)
    cycle = rptCycles(i);
    rptKey = sprintf('cycle_%d', cycle);
    
    % Collect RPT Static and OCV Capacity data
    if isfield(rptData, rptKey)
        channelCapacity.cycles = [channelCapacity.cycles, cycle];
        channelCapacity.capacity = [channelCapacity.capacity, rptData.(rptKey).static];
        channelCapacity.capacity_ocv = [channelCapacity.capacity_ocv, rptData.(rptKey).ocv];
    end
    
    % Static
    if isfield(rptData, rptKey)
        x_vals = [x_vals, x_pos];
        y_vals = [y_vals, rptData.(rptKey).static];
        plotInfo.(sprintf('rpt_%d_static', cycle)) = struct('x', x_pos, 'y', rptData.(rptKey).static, 'type', 'RPT', 'cycle', cycle, 'test', 'Static');
        x_pos = x_pos + 1;
    end
    
    % OCV
    if isfield(rptData, rptKey)
        x_vals = [x_vals, x_pos];
        y_vals = [y_vals, rptData.(rptKey).ocv];
        plotInfo.(sprintf('rpt_%d_ocv', cycle)) = struct('x', x_pos, 'y', rptData.(rptKey).ocv, 'type', 'RPT', 'cycle', cycle, 'test', 'OCV');
        x_pos = x_pos + 1;
    end
    
    % Add aging data between RPTs
    if i < length(rptCycles)
        startCycle = cycle;
        endCycle = rptCycles(i+1);
        agingKey = sprintf('period_%dto%d', startCycle, endCycle);
        
        if isfield(agingData, agingKey)
            agingStart = x_pos;
            agingEnd = x_pos + 7;  % 8 units for aging data
            agingX = linspace(agingStart, agingEnd, length(agingData.(agingKey).capacity));
            
            x_vals = [x_vals, agingX];
            y_vals = [y_vals, agingData.(agingKey).capacity];
            % Use actualMaxCycle from data as the actual end cycle
            actualEndCycle = agingData.(agingKey).actualMaxCycle;
            plotInfo.(sprintf('aging_%dto%d', startCycle, endCycle)) = struct('x', agingX, 'y', agingData.(agingKey).capacity, 'type', 'Aging', 'startCycle', startCycle, 'endCycle', endCycle, 'actualEndCycle', actualEndCycle, 'maxCycle', agingData.(agingKey).maxCycle);
            
            x_pos = agingEnd + 1;
        end
    end
end

% Add aging data after last RPT (if exists)
if length(agingCycles) > length(rptCycles)
    lastRPT = rptCycles(end);
    for i = length(rptCycles):length(agingCycles)-1
        startCycle = agingCycles(i);
        endCycle = agingCycles(i+1);
        agingKey = sprintf('period_%dto%d', startCycle, endCycle);
        
        if isfield(agingData, agingKey)
            agingStart = x_pos;
            agingEnd = x_pos + 7;  % 8 units for aging data
            agingX = linspace(agingStart, agingEnd, length(agingData.(agingKey).capacity));
            
            x_vals = [x_vals, agingX];
            y_vals = [y_vals, agingData.(agingKey).capacity];
            % Use actualMaxCycle from data as the actual end cycle
            actualEndCycle = agingData.(agingKey).actualMaxCycle;
            plotInfo.(sprintf('aging_%dto%d', startCycle, endCycle)) = struct('x', agingX, 'y', agingData.(agingKey).capacity, 'type', 'Aging', 'startCycle', startCycle, 'endCycle', endCycle, 'actualEndCycle', actualEndCycle, 'maxCycle', agingData.(agingKey).maxCycle);
            
            x_pos = agingEnd + 1;
        end
    end
end

% =========================================================================
% Create Plot
% =========================================================================

fig = figure('Name', sprintf('Channel %s - Discharged Capacity Retention', ch), 'Color', 'w');
hold on;

% Get initial RPT values (RPT0) for SOH calculation
if isfield(rptData, 'cycle_0')
    initialStatic = rptData.cycle_0.static;
    initialOCV = rptData.cycle_0.ocv;
else
    initialStatic = NaN;
    initialOCV = NaN;
    fprintf('Warning: RPT0 데이터 없음 - SOH 계산 불가\n');
end

% Plot RPT points
rptFields = fieldnames(plotInfo);
for i = 1:length(rptFields)
    field = rptFields{i};
    if strcmp(plotInfo.(field).type, 'RPT')
        if strcmp(plotInfo.(field).test, 'Static')
            plot(plotInfo.(field).x, plotInfo.(field).y, 'o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
                'MarkerFaceColor', color_static, 'MarkerEdgeColor', color_static, ...
                'DisplayName', sprintf('RPT%d Static', plotInfo.(field).cycle), 'Color', color_static);
            % Calculate SOH and display as percentage
            if ~isnan(initialStatic)
                soh = (plotInfo.(field).y / initialStatic) * 100;
                text(plotInfo.(field).x, plotInfo.(field).y+1, sprintf('%.2f (%.1f%%)', plotInfo.(field).y, soh), ...
                    'HorizontalAlignment', 'center', 'FontSize', 25, 'Color', color_static, 'FontWeight', 'bold');
            else
                text(plotInfo.(field).x, plotInfo.(field).y+1, sprintf('%.2f', plotInfo.(field).y), ...
                    'HorizontalAlignment', 'center', 'FontSize', 25, 'Color', color_static, 'FontWeight', 'bold');
            end
        else
            plot(plotInfo.(field).x, plotInfo.(field).y, 'o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
                'MarkerFaceColor', color_ocv, 'MarkerEdgeColor', color_ocv, ...
                'DisplayName', sprintf('RPT%d OCV', plotInfo.(field).cycle), 'Color', color_ocv);
            % Calculate SOH and display as percentage
            if ~isnan(initialOCV)
                soh = (plotInfo.(field).y / initialOCV) * 100;
                text(plotInfo.(field).x, plotInfo.(field).y+1, sprintf('%.2f (%.1f%%)', plotInfo.(field).y, soh), ...
                    'HorizontalAlignment', 'center', 'FontSize', 25, 'Color', color_ocv, 'FontWeight', 'bold');
            else
                text(plotInfo.(field).x, plotInfo.(field).y+1, sprintf('%.2f', plotInfo.(field).y), ...
                    'HorizontalAlignment', 'center', 'FontSize', 25, 'Color', color_ocv, 'FontWeight', 'bold');
            end
        end
    end
end

% Plot Aging data
for i = 1:length(rptFields)
    field = rptFields{i};
    if strcmp(plotInfo.(field).type, 'Aging')
        % Use endCycle for legend display
        plot(plotInfo.(field).x, plotInfo.(field).y, 'o', 'LineWidth', 1.2, 'MarkerSize', 3, ...
            'DisplayName', sprintf('Aging %dto%d', plotInfo.(field).startCycle + 1, plotInfo.(field).endCycle), 'Color', color_aging);
        
        % Add text for first and last points with same color (below the points)
        validIdx = ~isnan(plotInfo.(field).y);
        if any(validIdx)
            firstIdx = find(validIdx, 1);
            lastIdx = find(validIdx, 1, 'last');
            
            text(plotInfo.(field).x(firstIdx), plotInfo.(field).y(firstIdx)-1, sprintf('%.2f', plotInfo.(field).y(firstIdx)), ...
                'HorizontalAlignment', 'center', 'FontSize', 15, 'Color', color_aging, 'FontWeight', 'bold');
            text(plotInfo.(field).x(lastIdx), plotInfo.(field).y(lastIdx)-1, sprintf('%.2f', plotInfo.(field).y(lastIdx)), ...
                'HorizontalAlignment', 'center', 'FontSize', 15, 'Color', color_aging, 'FontWeight', 'bold');
        end
    end
end

% =========================================================================
% Format Plot
% =========================================================================

% Set axis limits
xlim([0, x_pos]);
ylim([40 66]);

% Create x-axis ticks and labels
xticks_vals = [];
xticklabels_vals = {};

% Add RPT ticks
for i = 1:length(rptCycles)
    cycle = rptCycles(i);
    rptKey = sprintf('cycle_%d', cycle);
    
    % Static tick
    staticField = sprintf('rpt_%d_static', cycle);
    if isfield(plotInfo, staticField)
        xticks_vals = [xticks_vals, plotInfo.(staticField).x];
        xticklabels_vals{end+1} = sprintf('RPT%d Static', cycle);
    end
    
    % OCV tick
    ocvField = sprintf('rpt_%d_ocv', cycle);
    if isfield(plotInfo, ocvField)
        xticks_vals = [xticks_vals, plotInfo.(ocvField).x];
        xticklabels_vals{end+1} = sprintf('RPT%d OCV', cycle);
    end
    
    % Add aging ticks (first and last)
    if i < length(rptCycles)
        startCycle = cycle;
        endCycle = rptCycles(i+1);
        agingField = sprintf('aging_%dto%d', startCycle, endCycle);
        
        if isfield(plotInfo, agingField)
            validIdx = ~isnan(plotInfo.(agingField).y);
            if any(validIdx)
                firstIdx = find(validIdx, 1);
                lastIdx = find(validIdx, 1, 'last');
                
                xticks_vals = [xticks_vals, plotInfo.(agingField).x(firstIdx), plotInfo.(agingField).x(lastIdx)];
                xticklabels_vals{end+1} = sprintf('Aging %d', startCycle + 1);
                xticklabels_vals{end+1} = sprintf('Aging %d', endCycle);
            end
        end
    end
end

% Add aging ticks after last RPT
if length(agingCycles) > length(rptCycles)
    for i = length(rptCycles):length(agingCycles)-1
        startCycle = agingCycles(i);
        endCycle = agingCycles(i+1);
        agingField = sprintf('aging_%dto%d', startCycle, endCycle);
        
        if isfield(plotInfo, agingField)
            validIdx = ~isnan(plotInfo.(agingField).y);
            if any(validIdx)
                firstIdx = find(validIdx, 1);
                lastIdx = find(validIdx, 1, 'last');
                
                xticks_vals = [xticks_vals, plotInfo.(agingField).x(firstIdx), plotInfo.(agingField).x(lastIdx)];
                xticklabels_vals{end+1} = sprintf('Aging %d', startCycle + 1);
                xticklabels_vals{end+1} = sprintf('Aging %d', endCycle);
            end
        end
    end
end

xticks(xticks_vals);
xticklabels(xticklabels_vals);
xtickangle(45);  % Rotate x-axis labels 45 degrees

xlabel('Test Type / Cycle Index', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Discharged Capacity [Ah]', 'FontSize', 20, 'FontWeight', 'bold');
title(sprintf('Channel %s - Capacity Retention', ch), 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'southeast');
grid on;

% =========================================================================
% Save
% =========================================================================

    savefig(fig, fullfile(saveFolder, sprintf('Ch%s_Capacity_Retention.fig', ch)));
    
    fprintf('Plot saved: Ch%s_Capacity_Retention.fig\n', ch);
    close(fig);  % Close figure to free memory
    
    % Store channel capacity data
    allChannelsCapacity.(sprintf('Ch%s', ch)) = channelCapacity;
end

% =========================================================================
% Save Capacity Data to MAT File
% =========================================================================

save(fullfile(saveFolder, 'Capacity_Data_Static.mat'), 'allChannelsCapacity');
fprintf('Capacity data saved to: Capacity_Data_Static.mat\n');

fprintf('All channels processed!\n');
