%% Field R1s SOC and Voltage Analysis
% This script creates two figures:
% 1. SOC bin (55-69%) vs R1s by year (2021-2025) using boxplot
% 2. Voltage (3.2-4.0V) vs R1s by year (2021-2025)

clear; clc; close all;

% Add daboxplot to path
addpath('frank-pk-DataViz-3.2.3.0/daboxplot');

%% Configuration
fontSize = 12;
saveDir = 'SOC_Voltage_Analysis_Results';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Load OldData and NewData results
fprintf('Loading OldData and NewData results...\n');

% OldData years
oldDataYears = {'2021', '2022', '2023'};
% NewData years  
newDataYears = {'2023', '2024', '2025'};

% Initialize data storage
allSOCData = [];
allVoltageData = [];
allR1sData = [];
allYearData = [];
allSourceData = [];

%% Load OldData results
fprintf('Loading OldData results...\n');
for i = 1:length(oldDataYears)
    year = oldDataYears{i};
    fileName = sprintf('R1s_Results_%s.mat', year);
    filePath = fullfile('R1s_Results_OldData', fileName);
    
    if exist(filePath, 'file')
        fprintf('Loading: %s\n', fileName);
        load(filePath, 'R1s_Results');
        
        % Extract SOC bin data
        dayNames = fieldnames(R1s_Results);
        for d = 1:length(dayNames)
            dayKey = dayNames{d};
            if isfield(R1s_Results.(dayKey), 'Rack01')
                rackData = R1s_Results.(dayKey).Rack01;
                
                % Extract SOC target data from raw dI, dV, SOC data
                if isfield(rackData, 'SOC_targets') && isfield(rackData, 'dI') && isfield(rackData, 'dV') && isfield(rackData, 'SOC')
                    socTargets = rackData.SOC_targets;
                    dI_data = rackData.dI;
                    dV_data = rackData.dV;
                    soc_data = rackData.SOC;
                    
                    % Debug: Check data size
                    if d == 1  % Only for first day
                        fprintf('SOC_targets: ');
                        fprintf('%d ', socTargets);
                        fprintf('\n');
                        fprintf('Raw data size - dI: %d, dV: %d, SOC: %d\n', length(dI_data), length(dV_data), length(soc_data));
                    end
                    
                    for t = 1:length(socTargets)
                        socTarget = socTargets(t);
                        soc_low = socTarget;      % 5% range
                        soc_high = socTarget + 4;
                        
                        % Find data points in this SOC target range
                        targetIdx = (soc_data >= soc_low) & (soc_data <= soc_high);
                        dI_target = dI_data(targetIdx);
                        dV_target = dV_data(targetIdx);
                        
                        % Debug: Check SOC range data count
                        if d == 1 && t == 1  % Only for first day, first target
                            fprintf('SOC %d%%-%d%% range: %d data points\n', socTarget, soc_high, length(dI_target));
                        end
                        
                        % Calculate R1s for this SOC target if enough data points
                        if length(dI_target) >= 5
                            mdl = fitlm(dI_target, dV_target, 'Intercept', false);
                            r1s_value = mdl.Coefficients.Estimate(1);
                            
                            if r1s_value > 0  % Only positive R1s values
                                allSOCData = [allSOCData; socTarget];
                                allR1sData = [allR1sData; r1s_value];
                                allYearData = [allYearData; str2double(year)];
                                allSourceData = [allSourceData; {'OldData'}];
                            end
                        end
                    end
                end
                
                % Debug: Check available fields
                if d == 1  % Only for first day to avoid spam
                    fprintf('Available fields in rackData: ');
                    fieldNames = fieldnames(rackData);
                    for f = 1:length(fieldNames)
                        fprintf('%s ', fieldNames{f});
                    end
                    fprintf('\n');
                end
                
                % Extract voltage data (if available)
                if isfield(rackData, 'V') && isfield(rackData, 'R1s')
                    V_data = rackData.V;
                    R1s_data = rackData.R1s;
                    
                    % Filter voltage range 3.2-4.0V
                    voltageValidIdx = V_data >= 3.2 & V_data <= 4.0;
                    if any(voltageValidIdx)
                        allVoltageData = [allVoltageData; V_data(voltageValidIdx)];
                        allR1sData = [allR1sData; R1s_data(voltageValidIdx)];
                        allYearData = [allYearData; repmat(str2double(year), sum(voltageValidIdx), 1)];
                        allSourceData = [allSourceData; repmat({'OldData'}, sum(voltageValidIdx), 1)];
                    end
                end
            end
        end
    else
        fprintf('File not found: %s\n', fileName);
    end
end

%% Load NewData results
fprintf('Loading NewData results...\n');
for i = 1:length(newDataYears)
    year = newDataYears{i};
    fileName = sprintf('R1s_Results_%s.mat', year);
    filePath = fullfile('R1s_Results_NewData', fileName);
    
    if exist(filePath, 'file')
        fprintf('Loading: %s\n', fileName);
        load(filePath, 'R1s_Results');
        
        % Extract SOC bin data
        dayNames = fieldnames(R1s_Results);
        for d = 1:length(dayNames)
            dayKey = dayNames{d};
            if isfield(R1s_Results.(dayKey), 'Rack01')
                rackData = R1s_Results.(dayKey).Rack01;
                
                % Extract SOC target data from raw dI, dV, SOC data
                if isfield(rackData, 'SOC_targets') && isfield(rackData, 'dI') && isfield(rackData, 'dV') && isfield(rackData, 'SOC')
                    socTargets = rackData.SOC_targets;
                    dI_data = rackData.dI;
                    dV_data = rackData.dV;
                    soc_data = rackData.SOC;
                    
                    % Debug: Check data size for NewData
                    if d == 1  % Only for first day
                        fprintf('NewData - SOC_targets: ');
                        fprintf('%d ', socTargets);
                        fprintf('\n');
                        fprintf('NewData - Raw data size - dI: %d, dV: %d, SOC: %d\n', length(dI_data), length(dV_data), length(soc_data));
                    end
                    
                    for t = 1:length(socTargets)
                        socTarget = socTargets(t);
                        soc_low = socTarget;      % 5% range
                        soc_high = socTarget + 4;
                        
                        % Find data points in this SOC target range
                        targetIdx = (soc_data >= soc_low) & (soc_data <= soc_high);
                        dI_target = dI_data(targetIdx);
                        dV_target = dV_data(targetIdx);
                        
                        % Debug: Check SOC range data count for NewData
                        if d == 1 && t == 1  % Only for first day, first target
                            fprintf('NewData - SOC %d%% range (%d-%d): %d data points\n', socTarget, soc_low, soc_high, length(dI_target));
                        end
                        
                        % Calculate R1s for this SOC target if enough data points
                        if length(dI_target) >= 5
                            mdl = fitlm(dI_target, dV_target, 'Intercept', false);
                            r1s_value = mdl.Coefficients.Estimate(1);
                            
                            if r1s_value > 0  % Only positive R1s values
                                allSOCData = [allSOCData; socTarget];
                                allR1sData = [allR1sData; r1s_value];
                                allYearData = [allYearData; str2double(year)];
                                allSourceData = [allSourceData; {'NewData'}];
                            end
                        end
                    end
                end
                
                % Extract voltage data (if available)
                if isfield(rackData, 'V') && isfield(rackData, 'R1s')
                    V_data = rackData.V;
                    R1s_data = rackData.R1s;
                    
                    % Filter voltage range 3.2-4.0V
                    voltageValidIdx = V_data >= 3.2 & V_data <= 4.0;
                    if any(voltageValidIdx)
                        allVoltageData = [allVoltageData; V_data(voltageValidIdx)];
                        allR1sData = [allR1sData; R1s_data(voltageValidIdx)];
                        allYearData = [allYearData; repmat(str2double(year), sum(voltageValidIdx), 1)];
                        allSourceData = [allSourceData; repmat({'NewData'}, sum(voltageValidIdx), 1)];
                    end
                end
            end
        end
    else
        fprintf('File not found: %s\n', fileName);
    end
end

%% Separate SOC and Voltage data
% SOC data
socValidIdx = ~isnan(allSOCData);
socData = allSOCData(socValidIdx);
socR1sData = allR1sData(socValidIdx);
socYearData = allYearData(socValidIdx);
socSourceData = allSourceData(socValidIdx);

% Voltage data
voltageValidIdx = ~isnan(allVoltageData);
voltageData = allVoltageData(voltageValidIdx);
voltageR1sData = allR1sData(voltageValidIdx);
voltageYearData = allYearData(voltageValidIdx);
voltageSourceData = allSourceData(voltageValidIdx);

fprintf('\nData Summary:\n');
fprintf('SOC data points: %d\n', length(socData));
fprintf('Voltage data points: %d\n', length(voltageData));

% Debug: Check R1s data range
if ~isempty(socData)
    fprintf('\nSOC R1s Data Debug:\n');
    fprintf('Min R1s: %.6f mΩ\n', min(socR1sData)*1000);
    fprintf('Max R1s: %.6f mΩ\n', max(socR1sData)*1000);
    fprintf('Mean R1s: %.6f mΩ\n', mean(socR1sData)*1000);
    fprintf('Negative R1s count: %d\n', sum(socR1sData < 0));
    
    if sum(socR1sData < 0) > 0
        fprintf('WARNING: Found negative R1s values!\n');
        negativeIdx = socR1sData < 0;
        fprintf('Negative R1s values: ');
        fprintf('%.6f ', socR1sData(negativeIdx)*1000);
        fprintf('mΩ\n');
    end
end

if ~isempty(voltageData)
    fprintf('\nVoltage R1s Data Debug:\n');
    fprintf('Min R1s: %.6f mΩ\n', min(voltageR1sData)*1000);
    fprintf('Max R1s: %.6f mΩ\n', max(voltageR1sData)*1000);
    fprintf('Mean R1s: %.6f mΩ\n', mean(voltageR1sData)*1000);
    fprintf('Negative R1s count: %d\n', sum(voltageR1sData < 0));
    
    if sum(voltageR1sData < 0) > 0
        fprintf('WARNING: Found negative R1s values!\n');
        negativeIdx = voltageR1sData < 0;
        fprintf('Negative R1s values: ');
        fprintf('%.6f ', voltageR1sData(negativeIdx)*1000);
        fprintf('mΩ\n');
    end
end

%% Figure 1: SOC target R1s by year (subplot format)
fprintf('Creating Figure 1: SOC target R1s by year...\n');

if ~isempty(socData)
    figure('Position', [100, 100, 1400, 1000]);
    
    % Define SOC targets (5% bins)
    socTargets = [55, 60, 65];  % 55-59%, 60-64%, 65-69%
    years = unique(socYearData);
    colors = lines(length(years));
    
    % Create 1x3 subplot layout
    for t = 1:length(socTargets)
        subplot(1, 3, t);
        
        socTarget = socTargets(t);
        
        % Find data for this SOC target
        targetIdx = socData == socTarget;
        targetSOC = socData(targetIdx);
        targetR1s = socR1sData(targetIdx);
        targetYears = socYearData(targetIdx);
        
        if ~isempty(targetR1s)
            % Group by year for boxplot
            yearGroups = [];
            yearLabels = [];
            yearPositions = [];
            
            for i = 1:length(years)
                year = years(i);
                yearIdx = targetYears == year;
                
                if any(yearIdx)
                    yearR1s = targetR1s(yearIdx) * 1000; % Convert to mΩ
                    yearGroups = [yearGroups; yearR1s];
                    yearLabels = [yearLabels; repmat(year, length(yearR1s), 1)];
                    yearPositions = [yearPositions; year];
                end
            end
            
            if ~isempty(yearGroups)
                % Simple boxplot with vector data
                uniqueYears = unique(yearLabels);
                
                % Create simple boxplot
                boxplot(yearGroups, yearLabels, 'Positions', uniqueYears, ...
                    'Widths', 0.6, 'Colors', 'b', 'Symbol', 'ko');
                
                % Add sample size labels and p-value
                for i = 1:length(uniqueYears)
                    year = uniqueYears(i);
                    yearIdx = yearLabels == year;
                    n_count = sum(yearIdx);
                    if n_count > 0
                        text(year, max(yearGroups(yearIdx)) + 0.01, sprintf('n=%d', n_count), ...
                            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                            'FontSize', 12);
                    end
                end
                
                % Calculate p-value for year comparison (ANOVA)
                if length(uniqueYears) >= 2
                    yearData = {};
                    for i = 1:length(uniqueYears)
                        year = uniqueYears(i);
                        yearIdx = yearLabels == year;
                        yearData{i} = yearGroups(yearIdx);
                        fprintf('Debug - Year %d: %d data points, range %.3f to %.3f\n', ...
                            year, length(yearData{i}), min(yearData{i}), max(yearData{i}));
                    end
                    
                    % Check for valid data
                    validYears = [];
                    for i = 1:length(yearData)
                        if length(yearData{i}) > 0 && ~all(isnan(yearData{i}))
                            validYears = [validYears, i];
                        end
                    end
                    
                    if length(validYears) >= 2
                        % Prepare data for ANOVA using cell array format
                        validData = {};
                        for i = 1:length(validYears)
                            idx = validYears(i);
                            validData{i} = yearData{idx};
                        end
                        
                        % Perform ANOVA test using vector format with group labels
                        allData = [];
                        groupLabels = [];
                        
                        for i = 1:length(validData)
                            data = validData{i};
                            allData = [allData; data(:)];
                            groupLabels = [groupLabels; repmat(i, length(data), 1)];
                        end
                        
                        [p_value, ~, ~] = anova1(allData, groupLabels, 'off');
                        
                        % Add p-value to subtitle
                        subtitle(sprintf('ANOVA p-value = %.3f', p_value), 'FontSize', fontSize-2);
                        
                        fprintf('SOC %d%%-%d%%: ANOVA p-value = %.10f\n', socTarget, socTarget+4, p_value);
                        fprintf('  Actual p-value: %e\n', p_value);
                    else
                        subtitle('ANOVA p-value = N/A (insufficient data)', 'FontSize', fontSize-2);
                        fprintf('SOC %d%%-%d%%: Insufficient valid data for ANOVA\n', socTarget, socTarget+4);
                    end
                end
            end
        end
        
        xlabel('Year', 'FontSize', fontSize);
        ylabel('R1s (mΩ)', 'FontSize', fontSize);
        title(sprintf('SOC %d%%-%d%%', socTarget, socTarget+4), 'FontSize', fontSize);
        grid on;
        
        % Set axis limits
        xlim([2020.5, 2025.5]);
        ylim([0, 1.7]);
        xticks(years);
    end
    
    % Save figure
    saveas(gcf, fullfile(saveDir, 'SOC_Target_R1s_by_Year.fig'));
    fprintf('Saved Figure 1: SOC target R1s by year\n');
else
    fprintf('No SOC data available for Figure 1\n');
end

%% Figure 2: Voltage vs R1s by year
fprintf('Creating Figure 2: Voltage vs R1s by year...\n');

if ~isempty(voltageData)
    figure('Position', [100, 100, 1200, 800]);
    
    % Define years and colors
    years = unique(voltageYearData);
    colors = lines(length(years));
    
    hold on;
    
    % Plot data for each year
    for i = 1:length(years)
        year = years(i);
        yearIdx = voltageYearData == year;
        
        yearVoltage = voltageData(yearIdx);
        yearR1s = voltageR1sData(yearIdx);
        
        % Create voltage bins
        voltageBins = 3.2:0.05:4.0;
        meanR1s = zeros(size(voltageBins));
        stdR1s = zeros(size(voltageBins));
        binCounts = zeros(size(voltageBins));
        
        for j = 1:length(voltageBins)
            if j < length(voltageBins)
                binIdx = yearVoltage >= voltageBins(j) & yearVoltage < voltageBins(j+1);
            else
                binIdx = yearVoltage >= voltageBins(j);
            end
            
            if any(binIdx)
                meanR1s(j) = mean(yearR1s(binIdx));
                stdR1s(j) = std(yearR1s(binIdx));
                binCounts(j) = sum(binIdx);
            else
                meanR1s(j) = NaN;
                stdR1s(j) = NaN;
                binCounts(j) = 0;
            end
        end
        
        % Remove NaN values
        validIdx = ~isnan(meanR1s) & binCounts > 0;
        validVoltage = voltageBins(validIdx);
        validMeanR1s = meanR1s(validIdx);
        validStdR1s = stdR1s(validIdx);
        
        % Plot with error bars
        if any(validIdx)
            errorbar(validVoltage, validMeanR1s*1000, validStdR1s*1000, 'o-', ...
                'Color', colors(i,:), 'LineWidth', 2, 'MarkerSize', 8, ...
                'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('%d', year));
        end
    end
    
    xlabel('Voltage (V)', 'FontSize', fontSize);
    ylabel('R1s (mΩ)', 'FontSize', fontSize);
    title('Voltage vs R1s by Year', 'FontSize', fontSize+2);
    legend('Location', 'best', 'FontSize', fontSize);
    grid on;
    
    % Set axis limits
    xlim([3.2, 4.0]);
    
    % Save figure
    saveas(gcf, fullfile(saveDir, 'Voltage_vs_R1s_by_Year.fig'));
    fprintf('Saved Figure 2: Voltage vs R1s by year\n');
else
    fprintf('No voltage data available for Figure 2\n');
end

%% Create summary statistics
fprintf('\n=== Summary Statistics ===\n');

% SOC data summary
if ~isempty(socData)
    fprintf('\nSOC Data Summary:\n');
    for i = 1:length(years)
        year = years(i);
        yearIdx = socYearData == year;
        yearSOC = socData(yearIdx);
        yearR1s = socR1sData(yearIdx);
        
        fprintf('Year %d: %d data points, Mean R1s: %.3f mΩ, Std: %.3f mΩ\n', ...
            year, length(yearR1s), mean(yearR1s)*1000, std(yearR1s)*1000);
    end
end

% Voltage data summary
if ~isempty(voltageData)
    fprintf('\nVoltage Data Summary:\n');
    voltageYears = unique(voltageYearData);
    for i = 1:length(voltageYears)
        year = voltageYears(i);
        yearIdx = voltageYearData == year;
        yearVoltage = voltageData(yearIdx);
        yearR1s = voltageR1sData(yearIdx);
        
        fprintf('Year %d: %d data points, Mean R1s: %.3f mΩ, Std: %.3f mΩ\n', ...
            year, length(yearR1s), mean(yearR1s)*1000, std(yearR1s)*1000);
    end
end

fprintf('\nAnalysis complete!\n');
