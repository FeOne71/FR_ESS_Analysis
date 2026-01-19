%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 02b_DriveCycle_BoxplotVisualization.m
% 박스플롯 시각화 (Capacity_C3 vs 저항값)
% 
% 목적: 
% - DriveCycle_Summary_Table.mat를 로드
% - X축: Capacity_C3 (용량값), Y축: 저항값
% - 스캐터 플롯 대신 박스플롯으로 시각화
% - DC별 subplot, SOC별 세분화 (SOC 90%, 70%, 50%)
% - 모든 저항값 (R_1s, R_3s, R_5s, R_10s, R_30s, R_60s)
% - 충전/방전 이벤트 구분
%
% 입력:
% - DriveCycle_Summary_Table.mat (02번 스크립트 출력)
%
% 출력:
% - figures/Boxplots/Boxplot_*.fig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('=== Drive Cycle Boxplot Visualization ===\n');

%% Configuration - User Settings
% =========================================================================
inputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
outputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
% =========================================================================

%% Load Summary Table
fprintf('\n=== Loading Summary Table ===\n');
summaryTablePath = fullfile(inputDir, 'DriveCycle_Summary_Table.mat');
if ~exist(summaryTablePath, 'file')
    fprintf('ERROR: DriveCycle_Summary_Table.mat not found!\n');
    fprintf('Expected path: %s\n', summaryTablePath);
    fprintf('Please run DriveCycle_DataAggregation_02.m first to generate the summary table.\n');
    return;
end

load(summaryTablePath, 'summaryTable');
fprintf('Loaded summary table: %d rows\n', height(summaryTable));

% Debug: Check table structure
fprintf('\n=== Table Structure Check ===\n');
fprintf('Table columns: %d\n', width(summaryTable));
if ismember('EventType', summaryTable.Properties.VariableNames)
    fprintf('Unique Event Types: %s\n', strjoin(unique(summaryTable.EventType), ', '));
else
    fprintf('WARNING: EventType column not found in table!\n');
end

%% Create figures directory
figuresDir = fullfile(outputDir, 'figures', 'Boxplots');
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

%% Generate Boxplot Visualizations (Capacity_C3 vs 저항값)
fprintf('\n=== Generating Boxplot Visualizations (Capacity_C3 vs Resistance) ===\n');

% Check for Capacity_C3 column
if ~ismember('Capacity_C3', summaryTable.Properties.VariableNames)
    fprintf('ERROR: Capacity_C3 column not found in summary table!\n');
    fprintf('Available columns: %s\n', strjoin(summaryTable.Properties.VariableNames, ', '));
    return;
end

if ~isempty(summaryTable) && height(summaryTable) > 0
    % Get unique values
    uniqueDCs = unique(summaryTable.DC_Profile);
    uniqueSOCs = unique(summaryTable.SOC(~isnan(summaryTable.SOC)));
    uniqueSOCs = sort(uniqueSOCs);
    uniqueEventTypes = unique(summaryTable.EventType);
    
    % Resistance time intervals
    resistanceVars = {'R_1s', 'R_3s', 'R_5s', 'R_10s', 'R_30s', 'R_60s'};
    resistanceLabels = {'R 1s', 'R 3s', 'R 5s', 'R 10s', 'R 30s', 'R 60s'};
    
    % Filter available resistance variables
    availableResistanceVars = {};
    availableResistanceLabels = {};
    for r_idx = 1:length(resistanceVars)
        if ismember(resistanceVars{r_idx}, summaryTable.Properties.VariableNames)
            availableResistanceVars{end+1} = resistanceVars{r_idx};
            availableResistanceLabels{end+1} = resistanceLabels{r_idx};
        end
    end
    
    fprintf('Available resistance variables: %d\n', length(availableResistanceVars));
    fprintf('Unique DC Profiles: %s\n', strjoin(uniqueDCs, ', '));
    fprintf('Unique SOC levels: %s\n', mat2str(uniqueSOCs'));
    fprintf('Unique Event Types: %s\n', strjoin(uniqueEventTypes, ', '));
    
    % Capacity_C3 range for binning
    allCapacity = summaryTable.Capacity_C3(~isnan(summaryTable.Capacity_C3));
    if ~isempty(allCapacity)
        minCapacity = min(allCapacity);
        maxCapacity = max(allCapacity);
        fprintf('Capacity_C3 range: [%.2f, %.2f] Ah\n', minCapacity, maxCapacity);
    else
        fprintf('ERROR: No valid Capacity_C3 data found!\n');
        return;
    end
    
    % Configuration: Number of capacity bins
    nCapacityBins = 10;  % Capacity_C3를 10개 구간으로 나눔
    
    % Target SOC levels to visualize
    targetSOCs = [90, 70, 50];
    
    % Create boxplots for each resistance variable
    for r_idx = 1:length(availableResistanceVars)
        rVar = availableResistanceVars{r_idx};
        rLabel = availableResistanceLabels{r_idx};
        
        fprintf('\n--- Creating boxplots for %s ---\n', rLabel);
        
        % Create figure for each event type
        for et_idx = 1:length(uniqueEventTypes)
            eventType = uniqueEventTypes{et_idx};
            
            % Filter data for this event type
            eventData = summaryTable(strcmp(summaryTable.EventType, eventType), :);
            
            if height(eventData) == 0
                fprintf('  No data for %s events. Skipping.\n', eventType);
                continue;
            end
            
            % Create figure for each SOC level
            for soc_idx = 1:length(targetSOCs)
                targetSOC = targetSOCs(soc_idx);
                
                % Filter data for this SOC level (allow small tolerance for floating point)
                socMask = abs(eventData.SOC - targetSOC) < 1.0;  % Within 1% tolerance
                socData = eventData(socMask, :);
                
                if height(socData) == 0
                    fprintf('  No data for %s events at SOC %d%%. Skipping.\n', eventType, targetSOC);
                    continue;
                end
                
                fprintf('  Processing SOC %d%%: %d events\n', targetSOC, height(socData));
                
                % Create figure: DC별 subplot (2 rows x 4 columns)
                fig = figure('Name', sprintf('Boxplot: Capacity_C3 vs %s - %s Events - SOC%d', rLabel, eventType, targetSOC), ...
                            'Position', [100 + r_idx*30 + et_idx*20 + soc_idx*10, 100, 1600, 900], 'Visible', 'on');
                
                nRows = 2;
                nCols = 4;
                
                % Plot each DC Profile
                for dc_idx = 1:length(uniqueDCs)
                    dcProfile = uniqueDCs{dc_idx};
                    
                    subplot(nRows, nCols, dc_idx);
                    hold on;
                    
                    % Filter data for this DC, EventType, and SOC
                    dcMask = strcmp(socData.DC_Profile, dcProfile) & ...
                             ~isnan(socData.Capacity_C3) & ...
                             ~isnan(socData.(rVar));
                    
                    dcData = socData(dcMask, :);
                
                if height(dcData) == 0
                    text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'middle', 'FontSize', 14, ...
                         'Color', 'red', 'Units', 'normalized');
                    title(sprintf('%s - %s - SOC%d%%', dcProfile, eventType, targetSOC), 'FontSize', 11, 'FontWeight', 'bold');
                    continue;
                end
                
                % Create capacity bins for this DC
                dcCapacity = dcData.Capacity_C3;
                dcResistance = dcData.(rVar);
                
                % Determine capacity bin edges (use data range for this DC)
                dcMinCap = min(dcCapacity);
                dcMaxCap = max(dcCapacity);
                capacityBinEdges = linspace(dcMinCap, dcMaxCap, nCapacityBins + 1);
                capacityBinCenters = (capacityBinEdges(1:end-1) + capacityBinEdges(2:end)) / 2;
                
                % Collect data for each capacity bin (단순히 Capacity 구간별로만)
                boxplotData = {};
                boxplotGroups = {};
                boxplotPositions = [];
                
                % For each capacity bin
                for bin_idx = 1:nCapacityBins
                    binMin = capacityBinEdges(bin_idx);
                    binMax = capacityBinEdges(bin_idx + 1);
                    
                    % Last bin includes the upper edge
                    if bin_idx == nCapacityBins
                        binMask = dcData.Capacity_C3 >= binMin & dcData.Capacity_C3 <= binMax;
                    else
                        binMask = dcData.Capacity_C3 >= binMin & dcData.Capacity_C3 < binMax;
                    end
                    
                    binResistance = dcData.(rVar)(binMask);
                    binResistance = binResistance(~isnan(binResistance));
                    
                    if ~isempty(binResistance) && length(binResistance) >= 1
                        boxplotData{end+1} = binResistance;
                        % Group label: Capacity bin center
                        binCenter = capacityBinCenters(bin_idx);
                        boxplotGroups{end+1} = sprintf('%.1f', binCenter);
                        boxplotPositions(end+1) = bin_idx;
                    end
                end
                
                % Create boxplot
                if ~isempty(boxplotData)
                    % Prepare data for boxplot
                    allValues = [];
                    allGroups = [];
                    
                    for b_idx = 1:length(boxplotData)
                        values = boxplotData{b_idx};
                        groupLabel = boxplotGroups{b_idx};
                        allValues = [allValues; values(:)];
                        allGroups = [allGroups; repmat({groupLabel}, length(values), 1)];
                    end
                    
                    % Create boxplot
                    boxplot(allValues, allGroups, 'Positions', boxplotPositions, ...
                           'Colors', 'b', 'Symbol', 'r+', 'OutlierSize', 4);
                    
                    % Set labels
                    xlabel('Capacity_C3 (Ah)', 'FontSize', 10, 'FontWeight', 'bold');
                    ylabel(sprintf('%s (mΩ)', rLabel), 'FontSize', 10, 'FontWeight', 'bold');
                    title(sprintf('%s - %s - SOC%d%%', dcProfile, eventType, targetSOC), 'FontSize', 11, 'FontWeight', 'bold');
                    grid on;
                    
                    % Reverse x-axis so capacity decreases from left to right
                    % (This makes degradation trend more intuitive: left to right = degradation)
                    set(gca, 'XDir', 'reverse');
                    
                    % Add sample count annotation
                    totalSamples = sum(cellfun(@length, boxplotData));
                    text(0.02, 0.98, sprintf('N=%d', totalSamples), ...
                         'Units', 'normalized', 'VerticalAlignment', 'top', ...
                         'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');
                else
                    % No data - show message
                    text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'middle', 'FontSize', 14, ...
                         'Color', 'red', 'Units', 'normalized');
                    title(sprintf('%s - %s - SOC%d%%', dcProfile, eventType, targetSOC), 'FontSize', 11, 'FontWeight', 'bold');
                end
                end  % End of DC Profile loop
                
                % Overall title
                sgtitle(sprintf('Capacity_C3 vs %s (%s Events) - SOC %d%%\nBoxplot by DC Profile', ...
                               rLabel, eventType, targetSOC), 'FontSize', 14, 'FontWeight', 'bold');
                
                % Save figure
                savePath = fullfile(figuresDir, sprintf('Boxplot_Capacity_vs_%s_%s_SOC%d.fig', rVar, eventType, targetSOC));
                saveas(fig, savePath);
                fprintf('  Saved: Boxplot_Capacity_vs_%s_%s_SOC%d.fig\n', rVar, eventType, targetSOC);
                close(fig);
            end  % End of SOC loop
        end  % End of Event Type loop
    end  % End of Resistance Variable loop
    
    fprintf('\nBoxplot visualizations complete.\n');
    fprintf('All figures saved to: %s\n', figuresDir);
else
    fprintf('WARNING: Summary table is empty. Cannot generate boxplots.\n');
end

fprintf('\n=== Boxplot Visualization Complete ===\n');

