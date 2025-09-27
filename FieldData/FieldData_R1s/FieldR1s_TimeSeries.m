%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldR1s_TimeSeries.m
% ESS Rack01 R1s Time Series Analysis
% Combine OldData and NewData results for time series visualization
% x-axis: YYYY:MM:DD, y-axis: Daily R1s (mΩ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

%% Paths
oldDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_R1s\R1s_Results_OldData';
newDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_R1s\R1s_Results_NewData';
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_R1s\TimeSeries_Results';

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameters
rackName = 'Rack01';
fontSize = 12;
lineWidth = 2;

%% Load OldData Results
fprintf('Loading OldData results...\n');
oldDataFiles = dir(fullfile(oldDataDir, 'R1s_Results_*.mat'));
oldDataResults = struct();

for i = 1:length(oldDataFiles)
    filePath = fullfile(oldDataDir, oldDataFiles(i).name);
    fprintf('Loading: %s\n', oldDataFiles(i).name);
    
    S = load(filePath);
    if isfield(S, 'R1s_Results')
        % Convert filename to valid field name (remove .mat extension)
        fieldName = oldDataFiles(i).name(1:end-4);
        oldDataResults.(fieldName) = S.R1s_Results;
    end
end

%% Load NewData Results
fprintf('Loading NewData results...\n');
newDataFiles = dir(fullfile(newDataDir, 'R1s_Results_*.mat'));
newDataResults = struct();

for i = 1:length(newDataFiles)
    filePath = fullfile(newDataDir, newDataFiles(i).name);
    fprintf('Loading: %s\n', newDataFiles(i).name);
    
    S = load(filePath);
    if isfield(S, 'R1s_Results')
        % Convert filename to valid field name (remove .mat extension)
        fieldName = newDataFiles(i).name(1:end-4);
        newDataResults.(fieldName) = S.R1s_Results;
    end
end

%% Extract and Combine Data
fprintf('Extracting and combining data...\n');

% Initialize arrays
allDates = [];
allR1s = [];
allR2 = [];
allN = [];
allSource = []; % 'OldData' or 'NewData'
allBatteryTemp = []; % Battery temperature
allAmbientTemp = []; % Ambient temperature

% Process OldData
oldFileNames = fieldnames(oldDataResults);
for i = 1:length(oldFileNames)
    fileName = oldFileNames{i};
    yearData = oldDataResults.(fileName);
    
    dayNames = fieldnames(yearData);
    for j = 1:length(dayNames)
        dayKey = dayNames{j};
        if isfield(yearData.(dayKey), rackName)
            % Extract date from dayKey (e.g., 'Raw_20210601' -> '2021-06-01')
            if length(dayKey) >= 12 && strcmp(dayKey(1:4), 'Raw_')
                dateStr = dayKey(5:12); % '20210601'
                year = str2double(dateStr(1:4));
                month = str2double(dateStr(5:6));
                day = str2double(dateStr(7:8));
                
                date = datetime(year, month, day);
                r1s = yearData.(dayKey).(rackName).R1s;
                r2 = yearData.(dayKey).(rackName).R_squared;
                n = yearData.(dayKey).(rackName).n_points;
                
                % Extract temperature data
                if isfield(yearData.(dayKey).(rackName), 'BatteryTemp')
                    batteryTemp = mean(yearData.(dayKey).(rackName).BatteryTemp);
                else
                    batteryTemp = NaN;
                end
                
                if isfield(yearData.(dayKey).(rackName), 'AmbientTemp')
                    ambientTemp = mean(yearData.(dayKey).(rackName).AmbientTemp);
                else
                    ambientTemp = NaN;
                end
                
                allDates = [allDates; date];
                allR1s = [allR1s; r1s];
                allR2 = [allR2; r2];
                allN = [allN; n];
                allSource = [allSource; {'OldData'}];
                allBatteryTemp = [allBatteryTemp; batteryTemp];
                allAmbientTemp = [allAmbientTemp; ambientTemp];
            end
        end
    end
end

% Process NewData
newFileNames = fieldnames(newDataResults);
for i = 1:length(newFileNames)
    fileName = newFileNames{i};
    yearData = newDataResults.(fileName);
    
    dayNames = fieldnames(yearData);
    for j = 1:length(dayNames)
        dayKey = dayNames{j};
        if isfield(yearData.(dayKey), rackName)
            % Extract date from dayKey (e.g., 'Raw_20230101' -> '2023-01-01')
            if length(dayKey) >= 12 && strcmp(dayKey(1:4), 'Raw_')
                dateStr = dayKey(5:12); % '20230101'
                year = str2double(dateStr(1:4));
                month = str2double(dateStr(5:6));
                day = str2double(dateStr(7:8));
                
                date = datetime(year, month, day);
                r1s = yearData.(dayKey).(rackName).R1s;
                r2 = yearData.(dayKey).(rackName).R_squared;
                n = yearData.(dayKey).(rackName).n_points;
                
                % Extract temperature data
                if isfield(yearData.(dayKey).(rackName), 'BatteryTemp')
                    batteryTemp = mean(yearData.(dayKey).(rackName).BatteryTemp);
                else
                    batteryTemp = NaN;
                end
                
                % NewData has no ambient temperature
                ambientTemp = NaN;
                
                allDates = [allDates; date];
                allR1s = [allR1s; r1s];
                allR2 = [allR2; r2];
                allN = [allN; n];
                allSource = [allSource; {'NewData'}];
                allBatteryTemp = [allBatteryTemp; batteryTemp];
                allAmbientTemp = [allAmbientTemp; ambientTemp];
            end
        end
    end
end

%% Sort by Date
fprintf('Sorting data by date...\n');
[allDates, sortIdx] = sort(allDates);
allR1s = allR1s(sortIdx);
allR2 = allR2(sortIdx);
allN = allN(sortIdx);
allSource = allSource(sortIdx);

% Create sequential x-axis (1, 2, 3, ...) instead of actual dates
xAxis = (1:length(allDates))';

%% Create Time Series Plot
fprintf('Creating time series plot...\n');

figure('Position', [100, 100, 1400, 800]);

% Separate OldData and NewData for different colors
oldDataIdx = strcmp(allSource, 'OldData');
newDataIdx = strcmp(allSource, 'NewData');

% Left y-axis: R1s with battery temperature color
yyaxis left;

% Plot all data as connected line with battery temperature color
scatter(xAxis, allR1s*1000, 50, allBatteryTemp, 'filled');
hold on;
plot(xAxis, allR1s*1000, '-', 'LineWidth', 1, 'Color', '#b0b0b0');

ylabel('R1s (mΩ)', 'FontSize', fontSize);
ylim([0.5 1]);

% Add trend line using linear regression
if length(allR1s) > 1
    mdl = fitlm(xAxis, allR1s*1000, 'Intercept', true);
    trendLine = mdl.Fitted;
    plot(xAxis, trendLine, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Trend (R²=%.3f)', mdl.Rsquared.Ordinary));
end

% Right y-axis: Ambient Temperature (only if data exists)
ambientTempValidIdx = ~isnan(allAmbientTemp);
if any(ambientTempValidIdx)
    yyaxis right;
    plot(xAxis(ambientTempValidIdx), allAmbientTemp(ambientTempValidIdx), 'b-', 'LineWidth', 1, 'DisplayName', 'Ambient Temp');
    ylabel('Ambient Temperature (°C)', 'FontSize', fontSize);
    ylim([20 40]);
end

% Colorbar for battery temperature (tie to left axis and fix color limits)
yyaxis left;
ax = gca;
ax.CLim = [20 35]; % Fixed color mapping range
colormap(ax, flipud(autumn));
c = colorbar(ax, 'eastoutside');
c.Label.String = 'Battery Temperature (°C)';
c.Label.FontSize = fontSize;
c.Limits = [20 35]; % Fixed colorbar display range

% Formatting
xlabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
title('ESS Rack01 Daily R1s Time Series with Temperature', 'FontSize', fontSize+2, 'Color', 'k');
grid on;
legend('Location', 'best', 'FontSize', fontSize);

% Set all axis colors to black
set(gca, 'XColor', 'k', 'YColor', 'k');
yyaxis left;
set(gca, 'YColor', 'k');
yyaxis right;
set(gca, 'YColor', 'k');

% Format x-axis - show only one label per month
% Find unique months and their first occurrence
yearVals = zeros(size(allDates));
monthVals = zeros(size(allDates));
for i = 1:length(allDates)
    yearVals(i) = allDates(i).Year;
    monthVals(i) = allDates(i).Month;
end
[uniqueMonths, uniqueIdx] = unique(yearVals*100 + monthVals, 'stable');
xticks(xAxis(uniqueIdx));
xticklabels(datestr(allDates(uniqueIdx), 'yyyy-mm'));
xtickangle(45);

% Add statistics text
if ~isempty(allR1s)
    meanR1s = mean(allR1s)*1000;
    stdR1s = std(allR1s)*1000;
    minR1s = min(allR1s)*1000;
    maxR1s = max(allR1s)*1000;
    
    statsText = sprintf('Statistics:\nMean: %.3f mΩ\nStd: %.3f mΩ\nMin: %.3f mΩ\nMax: %.3f mΩ\nTotal Days: %d', ...
        meanR1s, stdR1s, minR1s, maxR1s, length(allR1s));
    
    text(0.02, 0.98, statsText, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', fontSize-2, ...
        'BackgroundColor', 'white', 'EdgeColor', 'black');
end

%% Save Results
% Save figure
saveas(gcf, fullfile(saveDir, 'R1s_TimeSeries.fig'));

% Save as mat file
save(fullfile(saveDir, 'R1s_TimeSeries_Data.mat'), 'allDates', 'allR1s', 'allR2', 'allN', 'allSource');

%% Second Figure: June Data Only
fprintf('Creating June data plot...\n');

% Extract June data only
juneIdx = monthVals == 6;
juneDates = allDates(juneIdx);
juneR1s = allR1s(juneIdx);
juneSource = allSource(juneIdx);
% Create sequential x-axis for June data (1, 2, 3, ...)
juneXAxis = (1:length(juneDates))';

if any(juneIdx)
    figure('Position', [100, 100, 1400, 600]);
    
    % Separate OldData and NewData for June
    juneOldIdx = strcmp(juneSource, 'OldData');
    juneNewIdx = strcmp(juneSource, 'NewData');
    
    % Plot all June data as connected line first
    plot(juneXAxis, juneR1s*1000, 'k-', 'LineWidth', 1, 'DisplayName', 'June R1s Trend');
    hold on;
    
    % Plot June data markers
    if any(juneOldIdx)
        plot(juneXAxis(juneOldIdx), juneR1s(juneOldIdx)*1000, 'bo-', ...
            'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'OldData (June)');
    end
    
    if any(juneNewIdx)
        plot(juneXAxis(juneNewIdx), juneR1s(juneNewIdx)*1000, 'ro-', ...
            'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'NewData (June)');
    end
    
    % Formatting
    xlabel('Date (YYYY-MM)', 'FontSize', fontSize);
    ylabel('R1s (mΩ)', 'FontSize', fontSize);
    title('ESS Rack01 June R1s Comparison', 'FontSize', fontSize+2);
    grid on;
    legend('Location', 'best', 'FontSize', fontSize);
    
    % Format x-axis for June data - same as Figure 1 (sequential with one label per month)
    juneYearVals = zeros(size(juneDates));
    juneMonthVals = zeros(size(juneDates));
    for i = 1:length(juneDates)
        juneYearVals(i) = juneDates(i).Year;
        juneMonthVals(i) = juneDates(i).Month;
    end
    [juneUniqueMonths, juneUniqueIdx] = unique(juneYearVals*100 + juneMonthVals, 'stable');
    xticks(juneXAxis(juneUniqueIdx));
    xticklabels(datestr(juneDates(juneUniqueIdx), 'yyyy-mm'));
    xtickangle(45);
    
    % Save June figure
    saveas(gcf, fullfile(saveDir, 'R1s_June_Comparison.fig'));
end


%% Third Figure: R1s Time Series with Temperature
fprintf('Creating R1s time series with temperature...\n');

% Remove outliers using IQR method (all data)
Q1 = quantile(allR1s, 0.25);
Q3 = quantile(allR1s, 0.75);
IQR = Q3 - Q1;
lowerBound = Q1 - 1.5 * IQR;
upperBound = Q3 + 1.5 * IQR;



% Filter data (outlier removal)
outlierValidIdx = allR1s >= lowerBound & allR1s <= upperBound;
filteredDates = allDates(outlierValidIdx);
filteredR1s = allR1s(outlierValidIdx);
filteredBatteryTemp = allBatteryTemp(outlierValidIdx);
filteredAmbientTemp = allAmbientTemp(outlierValidIdx);
filteredSource = allSource(outlierValidIdx);

% Use the same x-axis as Figure 1 (1 to total data length)
% Find the indices of filtered data in the original x-axis
[~, filteredIndices] = ismember(filteredDates, allDates);
xAxisTemp = filteredIndices;

% Remove NaN values for battery temperature data only
battTempValidIdx = ~isnan(filteredBatteryTemp);
xAxisTemp = xAxisTemp(battTempValidIdx);
allR1sTemp = filteredR1s(battTempValidIdx);
allBatteryTempClean = filteredBatteryTemp(battTempValidIdx);
allAmbientTempClean = filteredAmbientTemp(battTempValidIdx);
allSourceTemp = filteredSource(battTempValidIdx);

if length(allR1sTemp) > 0
    figure('Position', [100, 100, 1400, 800]);
    
    % Separate OldData and NewData for different colors
    oldDataIdx = strcmp(allSourceTemp, 'OldData');
    newDataIdx = strcmp(allSourceTemp, 'NewData');
    
    % Left y-axis: R1s with battery temperature color
    yyaxis left;
    
        % Plot all data as connected line with battery temperature color
        scatter(xAxisTemp, allR1sTemp*1000, 50, allBatteryTempClean, 'filled');
        hold on;
        plot(xAxisTemp, allR1sTemp*1000, '-', 'LineWidth', 1, 'Color', '#b0b0b0');
    
    ylabel('R1s (mΩ)', 'FontSize', fontSize);
    ylim([0.5 1]);
    
    % Add trend line using linear regression
    if length(allR1sTemp) > 1
        mdl = fitlm(xAxisTemp, allR1sTemp*1000, 'Intercept', true);
        trendLine = mdl.Fitted;
        plot(xAxisTemp, trendLine, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Trend (R²=%.3f)', mdl.Rsquared.Ordinary));
    end
    
    % Right y-axis: Ambient Temperature (only if data exists)
    ambientTempValidIdx = ~isnan(allAmbientTempClean);
    if any(ambientTempValidIdx)
        yyaxis right;
        plot(xAxisTemp(ambientTempValidIdx), allAmbientTempClean(ambientTempValidIdx), 'b-', 'LineWidth', 1, 'DisplayName', 'Ambient Temp');
        ylabel('Ambient Temperature (°C)', 'FontSize', fontSize);
        ylim([20 40]);
    end
    
    % Colorbar for battery temperature (tie to left axis and fix color limits)
    yyaxis left;
    ax = gca;
    ax.CLim = [20 35];
    colormap(ax, flipud(autumn));
    c = colorbar(ax, 'eastoutside');
    c.Label.String = 'Battery Temperature (°C)';
    c.Label.FontSize = fontSize;
    c.Limits = [20 35];
    
    % Formatting
    xlabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
    title('ESS Rack01 Daily R1s Time Series with Temperature (Outliers Removed)', 'FontSize', fontSize+2, 'Color', 'k');
    grid on;
    legend('Location', 'best', 'FontSize', fontSize);
    
    % Set all axis colors to black
    set(gca, 'XColor', 'k', 'YColor', 'k');
    yyaxis left;
    set(gca, 'YColor', 'k');
    yyaxis right;
    set(gca, 'YColor', 'k');
    
    % Format x-axis - show only one label per month (same as Figure 1)
    xticks(xAxis(uniqueIdx));
    xticklabels(datestr(allDates(uniqueIdx), 'yyyy-mm'));
    xtickangle(45);
    
    % Add statistics text
    meanR1sTemp = mean(allR1sTemp)*1000;
    stdR1sTemp = std(allR1sTemp)*1000;
    minR1sTemp = min(allR1sTemp)*1000;
    maxR1sTemp = max(allR1sTemp)*1000;
    
    statsTextTemp = sprintf('Statistics:\nMean: %.3f mΩ\nStd: %.3f mΩ\nMin: %.3f mΩ\nMax: %.3f mΩ\nValid Days: %d', ...
        meanR1sTemp, stdR1sTemp, minR1sTemp, maxR1sTemp, length(allR1sTemp));
    
    text(0.02, 0.98, statsTextTemp, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', fontSize-2, ...
        'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % Save Figure with temperature
    saveas(gcf, fullfile(saveDir, 'R1s_TimeSeries_WithTemperature.fig'));
    fprintf('Saved R1s time series with temperature\n');
end

%% Fourth Figure: 3D Plot (Figure 3 without Ambient Temp)
fprintf('Creating 3D plot...\n');

% Use the same data as Figure 3 (filtered data without ambient temp)
if length(allR1sTemp) > 0
    figure('Position', [100, 100, 1200, 800]);
    
    % Create 3D scatter plot with same color scheme as Figure 3
    scatter3(allBatteryTempClean, xAxisTemp, allR1sTemp*1000, 50, allBatteryTempClean, 'filled');
    hold on;
    
    % Add 3D line connecting the points
    % plot3(allBatteryTempClean, xAxisTemp, allR1sTemp*1000, '-', 'LineWidth', 1, 'Color', '#b0b0b0');
    
    % Apply same colormap as Figure 3
    colormap(flipud(autumn));
    
    % Formatting
    xlabel('Battery Temperature (°C)', 'FontSize', fontSize, 'Color', 'k');
    ylabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
    zlabel('R1s (mΩ)', 'FontSize', fontSize, 'Color', 'k');
    title('3D Plot: Battery Temperature vs Date vs R1s (Outliers Removed)', 'FontSize', fontSize+2, 'Color', 'k');
    
    % Set axis limits and ticks for 3D plot
    xlim([min(allBatteryTempClean), max(allBatteryTempClean)]);
    ylim([min(xAxisTemp), max(xAxisTemp)]);
    zlim([0.5 1]);
    
    % Set y-axis (Date) ticks and labels
    set(gca, 'YTick', xAxis(uniqueIdx));
    set(gca, 'YTickLabel', datestr(allDates(uniqueIdx), 'yyyy-mm'));
    
    % Add colorbar for battery temperature (match Figures 1 & 3)
    ax3d = gca;
    ax3d.CLim = [20 35];
    colormap(ax3d, flipud(autumn));
    c = colorbar(ax3d, 'eastoutside');
    c.Label.String = 'Battery Temperature (°C)';
    c.Label.FontSize = fontSize;
    c.Limits = [20 35];
    
    % Add grid
    grid on;
    
    % Set all axis colors to black (same as Figure 3)
    set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
    
    % Set view angle for better visualization
    view(45, 30);
    
    % Save 3D figure
    saveas(gcf, fullfile(saveDir, 'R1s_3D_Plot.fig'));
    fprintf('Saved 3D plot\n');
end

%% Sixth Figure: 3D Plot (No Outlier Removal)
fprintf('Creating 3D plot without outlier removal...\n');

% Use all data (no outlier removal)
if length(allR1s) > 0
    % Extract month and year information
    yearVals3D = arrayfun(@(x) x.Year, allDates);
    monthVals3D = arrayfun(@(x) x.Month, allDates);
    
    % Remove NaN values for battery temperature data only
    battTempValidIdxAll = ~isnan(allBatteryTemp);
    xAxisAll = xAxis(battTempValidIdxAll);
    allR1sAll = allR1s(battTempValidIdxAll);
    allBatteryTempAll = allBatteryTemp(battTempValidIdxAll);
    allAmbientTempAll = allAmbientTemp(battTempValidIdxAll);
    
    figure('Position', [100, 100, 1200, 800]);
    
    % Create 3D scatter plot with same color scheme
    scatter3(allBatteryTempAll, xAxisAll, allR1sAll*1000, 50, allBatteryTempAll, 'filled');
    
    % Apply same colormap
    colormap(flipud(autumn));
    
    % Formatting
    xlabel('Battery Temperature (°C)', 'FontSize', fontSize, 'Color', 'k');
    ylabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
    zlabel('R1s (mΩ)', 'FontSize', fontSize, 'Color', 'k');
    title('3D Plot: Battery Temperature vs Date vs R1s (No Outlier Removal)', 'FontSize', fontSize+2, 'Color', 'k');
    
    % Set axis limits and ticks for 3D plot
    xlim([min(allBatteryTempAll), max(allBatteryTempAll)]);
    ylim([min(xAxisAll), max(xAxisAll)]);
    zlim([0.5 1]);
    
    % Set y-axis (Date) ticks and labels
    set(gca, 'YTick', xAxis(uniqueIdx));
    set(gca, 'YTickLabel', datestr(allDates(uniqueIdx), 'yyyy-mm'));
    
    % Add colorbar for battery temperature (match Figures 1 & 3)
    ax3d = gca;
    ax3d.CLim = [20 35];
    colormap(ax3d, flipud(autumn));
    c = colorbar(ax3d, 'eastoutside');
    c.Label.String = 'Battery Temperature (°C)';
    c.Label.FontSize = fontSize;
    c.Limits = [20 35];
    
    % Add grid
    grid on;
    
    % Set all axis colors to black
    set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
    
    % Set view angle for better visualization
    view(45, 30);
    
    % Save 3D figure
    saveas(gcf, fullfile(saveDir, 'R1s_3D_Plot_NoOutlierRemoval.fig'));
    fprintf('Saved 3D plot (no outlier removal)\n');
end

%% Fifth Figure: 3-Sigma Outlier Removal
fprintf('Creating 3-sigma outlier removal plot...\n');

% Calculate 3-sigma bounds
meanR1s = mean(allR1s);
stdR1s = std(allR1s);
sigma3LowerBound = meanR1s - 3 * stdR1s;
sigma3UpperBound = meanR1s + 3 * stdR1s;

% Filter data using 3-sigma method
sigma3ValidIdx = allR1s >= sigma3LowerBound & allR1s <= sigma3UpperBound;
sigma3FilteredDates = allDates(sigma3ValidIdx);
sigma3FilteredR1s = allR1s(sigma3ValidIdx);
sigma3FilteredBatteryTemp = allBatteryTemp(sigma3ValidIdx);
sigma3FilteredAmbientTemp = allAmbientTemp(sigma3ValidIdx);
sigma3FilteredSource = allSource(sigma3ValidIdx);

% Use the same x-axis as Figure 1
[~, sigma3FilteredIndices] = ismember(sigma3FilteredDates, allDates);
sigma3XAxisTemp = sigma3FilteredIndices;

% Remove NaN values for battery temperature data only
sigma3BattTempValidIdx = ~isnan(sigma3FilteredBatteryTemp);
sigma3XAxisTemp = sigma3XAxisTemp(sigma3BattTempValidIdx);
sigma3R1sTemp = sigma3FilteredR1s(sigma3BattTempValidIdx);
sigma3BatteryTempClean = sigma3FilteredBatteryTemp(sigma3BattTempValidIdx);
sigma3AmbientTempClean = sigma3FilteredAmbientTemp(sigma3BattTempValidIdx);
sigma3SourceTemp = sigma3FilteredSource(sigma3BattTempValidIdx);

if length(sigma3R1sTemp) > 0
    figure('Position', [100, 100, 1400, 800]);
    
    % Separate OldData and NewData for different colors
    sigma3OldDataIdx = strcmp(sigma3SourceTemp, 'OldData');
    sigma3NewDataIdx = strcmp(sigma3SourceTemp, 'NewData');
    
    % Left y-axis: R1s with battery temperature color
    yyaxis left;
    
        % Plot all data as connected line with battery temperature color
        scatter(sigma3XAxisTemp, sigma3R1sTemp*1000, 50, sigma3BatteryTempClean, 'filled');
        hold on;
        plot(sigma3XAxisTemp, sigma3R1sTemp*1000, '-', 'LineWidth', 1, 'Color', '#b0b0b0');
    
    ylabel('R1s (mΩ)', 'FontSize', fontSize);
    ylim([0.5 1]);
    
    % Add trend line using linear regression
    if length(sigma3R1sTemp) > 1
        mdl = fitlm(sigma3XAxisTemp, sigma3R1sTemp*1000, 'Intercept', true);
        trendLine = mdl.Fitted;
        plot(sigma3XAxisTemp, trendLine, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Trend (R²=%.3f)', mdl.Rsquared.Ordinary));
    end
    
    % Right y-axis: Ambient Temperature (only if data exists)
    sigma3AmbientTempValidIdx = ~isnan(sigma3AmbientTempClean);
    if any(sigma3AmbientTempValidIdx)
        yyaxis right;
        plot(sigma3XAxisTemp(sigma3AmbientTempValidIdx), sigma3AmbientTempClean(sigma3AmbientTempValidIdx), 'b-', 'LineWidth', 2, 'DisplayName', 'Ambient Temp');
        ylabel('Ambient Temperature (°C)', 'FontSize', fontSize);
    end
    
    % Colorbar for battery temperature
    c = colorbar('eastoutside');
    c.Label.String = 'Battery Temperature (°C)';
    c.Label.FontSize = fontSize;
    colormap(flipud(autumn));
    
    % Formatting
    xlabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
    title('ESS Rack01 Daily R1s Time Series with Temperature (3-Sigma Outliers Removed)', 'FontSize', fontSize+2, 'Color', 'k');
    grid on;
    legend('Location', 'best', 'FontSize', fontSize);
    
    % Set all axis colors to black
    set(gca, 'XColor', 'k', 'YColor', 'k');
    yyaxis left;
    set(gca, 'YColor', 'k');
    yyaxis right;
    set(gca, 'YColor', 'k');
    
    % Format x-axis - show only one label per month (same as Figure 1)
    xticks(xAxis(uniqueIdx));
    xticklabels(datestr(allDates(uniqueIdx), 'yyyy-mm'));
    xtickangle(45);
    
    % Add statistics text
    sigma3MeanR1s = mean(sigma3R1sTemp)*1000;
    sigma3StdR1s = std(sigma3R1sTemp)*1000;
    sigma3MinR1s = min(sigma3R1sTemp)*1000;
    sigma3MaxR1s = max(sigma3R1sTemp)*1000;
    
    statsTextSigma3 = sprintf('Statistics (3-Sigma):\nMean: %.3f mΩ\nStd: %.3f mΩ\nMin: %.3f mΩ\nMax: %.3f mΩ\nValid Days: %d\nBounds: %.3f-%.3f mΩ', ...
        sigma3MeanR1s, sigma3StdR1s, sigma3MinR1s, sigma3MaxR1s, length(sigma3R1sTemp), ...
        sigma3LowerBound*1000, sigma3UpperBound*1000);
    
    text(0.02, 0.98, statsTextSigma3, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', fontSize-2, ...
        'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % Save Figure with 3-sigma outlier removal
    saveas(gcf, fullfile(saveDir, 'R1s_TimeSeries_3Sigma.fig'));
    fprintf('Saved 3-sigma outlier removal plot\n');
end

%% Create Statistics Tables
fprintf('Creating statistics tables...\n');

% Extract year and month information
yearVals = arrayfun(@(x) x.Year, allDates);
monthVals = arrayfun(@(x) x.Month, allDates);

% Table 1: Before outlier removal (daily data)
tableBeforeOutlier = table(allDates, allR1s*1000, allR2, allN, allSource, allBatteryTemp, allAmbientTemp, ...
    yearVals, monthVals, ...
    'VariableNames', {'Date', 'R1s_mOhm', 'R_squared', 'N_points', 'Source', 'BatteryTemp_C', 'AmbientTemp_C', 'Year', 'Month'});

% Table 2: After outlier removal (daily data)
if length(allR1sTemp) > 0
    yearValsFiltered = arrayfun(@(x) x.Year, filteredDates(battTempValidIdx));
    monthValsFiltered = arrayfun(@(x) x.Month, filteredDates(battTempValidIdx));
    
    % Get the corresponding R2 and N values for filtered data
    filteredR2 = allR2(outlierValidIdx);
    filteredN = allN(outlierValidIdx);
    
    tableAfterOutlier = table(filteredDates(battTempValidIdx), allR1sTemp*1000, ...
        filteredR2(battTempValidIdx), filteredN(battTempValidIdx), ...
        allSourceTemp, allBatteryTempClean, allAmbientTempClean, ...
        yearValsFiltered, monthValsFiltered, ...
        'VariableNames', {'Date', 'R1s_mOhm', 'R_squared', 'N_points', 'Source', 'BatteryTemp_C', 'AmbientTemp_C', 'Year', 'Month'});
else
    tableAfterOutlier = table();
end

% Table 3: Monthly statistics (before outlier removal)
[uniqueYearMonth, ~, idx] = unique([tableBeforeOutlier.Year, tableBeforeOutlier.Month], 'rows');
monthlyStats = table();
for i = 1:size(uniqueYearMonth, 1)
    mask = idx == i;
    data = tableBeforeOutlier(mask, :);
    
    newRow = table(uniqueYearMonth(i,1), uniqueYearMonth(i,2), ...
        mean(data.R1s_mOhm), std(data.R1s_mOhm), min(data.R1s_mOhm), max(data.R1s_mOhm), height(data), ...
        mean(data.R_squared), std(data.R_squared), min(data.R_squared), max(data.R_squared), ...
        mean(data.BatteryTemp_C, 'omitnan'), std(data.BatteryTemp_C, 'omitnan'), min(data.BatteryTemp_C), max(data.BatteryTemp_C), ...
        mean(data.AmbientTemp_C, 'omitnan'), std(data.AmbientTemp_C, 'omitnan'), min(data.AmbientTemp_C), max(data.AmbientTemp_C), ...
        'VariableNames', {'Year', 'Month', 'R1s_mean', 'R1s_std', 'R1s_min', 'R1s_max', 'R1s_count', ...
        'R2_mean', 'R2_std', 'R2_min', 'R2_max', ...
        'BattTemp_mean', 'BattTemp_std', 'BattTemp_min', 'BattTemp_max', ...
        'AmbTemp_mean', 'AmbTemp_std', 'AmbTemp_min', 'AmbTemp_max'});
    monthlyStats = [monthlyStats; newRow];
end

% Table 4: Yearly statistics (before outlier removal)
[uniqueYears, ~, idx] = unique(tableBeforeOutlier.Year);
yearlyStats = table();
for i = 1:length(uniqueYears)
    mask = idx == i;
    data = tableBeforeOutlier(mask, :);
    
    newRow = table(uniqueYears(i), ...
        mean(data.R1s_mOhm), std(data.R1s_mOhm), min(data.R1s_mOhm), max(data.R1s_mOhm), height(data), ...
        mean(data.R_squared), std(data.R_squared), min(data.R_squared), max(data.R_squared), ...
        mean(data.BatteryTemp_C, 'omitnan'), std(data.BatteryTemp_C, 'omitnan'), min(data.BatteryTemp_C), max(data.BatteryTemp_C), ...
        mean(data.AmbientTemp_C, 'omitnan'), std(data.AmbientTemp_C, 'omitnan'), min(data.AmbientTemp_C), max(data.AmbientTemp_C), ...
        'VariableNames', {'Year', 'R1s_mean', 'R1s_std', 'R1s_min', 'R1s_max', 'R1s_count', ...
        'R2_mean', 'R2_std', 'R2_min', 'R2_max', ...
        'BattTemp_mean', 'BattTemp_std', 'BattTemp_min', 'BattTemp_max', ...
        'AmbTemp_mean', 'AmbTemp_std', 'AmbTemp_min', 'AmbTemp_max'});
    yearlyStats = [yearlyStats; newRow];
end

% Table 5: Monthly statistics (after outlier removal)
if ~isempty(tableAfterOutlier)
    [uniqueYearMonth, ~, idx] = unique([tableAfterOutlier.Year, tableAfterOutlier.Month], 'rows');
    monthlyStatsFiltered = table();
    for i = 1:size(uniqueYearMonth, 1)
        mask = idx == i;
        data = tableAfterOutlier(mask, :);
        
        newRow = table(uniqueYearMonth(i,1), uniqueYearMonth(i,2), ...
            mean(data.R1s_mOhm), std(data.R1s_mOhm), min(data.R1s_mOhm), max(data.R1s_mOhm), height(data), ...
            mean(data.R_squared), std(data.R_squared), min(data.R_squared), max(data.R_squared), ...
            mean(data.BatteryTemp_C, 'omitnan'), std(data.BatteryTemp_C, 'omitnan'), min(data.BatteryTemp_C), max(data.BatteryTemp_C), ...
            mean(data.AmbientTemp_C, 'omitnan'), std(data.AmbientTemp_C, 'omitnan'), min(data.AmbientTemp_C), max(data.AmbientTemp_C), ...
            'VariableNames', {'Year', 'Month', 'R1s_mean', 'R1s_std', 'R1s_min', 'R1s_max', 'R1s_count', ...
            'R2_mean', 'R2_std', 'R2_min', 'R2_max', ...
            'BattTemp_mean', 'BattTemp_std', 'BattTemp_min', 'BattTemp_max', ...
            'AmbTemp_mean', 'AmbTemp_std', 'AmbTemp_min', 'AmbTemp_max'});
        monthlyStatsFiltered = [monthlyStatsFiltered; newRow];
    end
    
    % Table 6: Yearly statistics (after outlier removal)
    [uniqueYears, ~, idx] = unique(tableAfterOutlier.Year);
    yearlyStatsFiltered = table();
    for i = 1:length(uniqueYears)
        mask = idx == i;
        data = tableAfterOutlier(mask, :);
        
        newRow = table(uniqueYears(i), ...
            mean(data.R1s_mOhm), std(data.R1s_mOhm), min(data.R1s_mOhm), max(data.R1s_mOhm), height(data), ...
            mean(data.R_squared), std(data.R_squared), min(data.R_squared), max(data.R_squared), ...
            mean(data.BatteryTemp_C, 'omitnan'), std(data.BatteryTemp_C, 'omitnan'), min(data.BatteryTemp_C), max(data.BatteryTemp_C), ...
            mean(data.AmbientTemp_C, 'omitnan'), std(data.AmbientTemp_C, 'omitnan'), min(data.AmbientTemp_C), max(data.AmbientTemp_C), ...
            'VariableNames', {'Year', 'R1s_mean', 'R1s_std', 'R1s_min', 'R1s_max', 'R1s_count', ...
            'R2_mean', 'R2_std', 'R2_min', 'R2_max', ...
            'BattTemp_mean', 'BattTemp_std', 'BattTemp_min', 'BattTemp_max', ...
            'AmbTemp_mean', 'AmbTemp_std', 'AmbTemp_min', 'AmbTemp_max'});
        yearlyStatsFiltered = [yearlyStatsFiltered; newRow];
    end
else
    monthlyStatsFiltered = table();
    yearlyStatsFiltered = table();
end

% Save tables as mat files
save(fullfile(saveDir, 'R1s_Statistics_BeforeOutlierRemoval.mat'), 'tableBeforeOutlier', 'monthlyStats', 'yearlyStats');
if ~isempty(tableAfterOutlier)
    save(fullfile(saveDir, 'R1s_Statistics_AfterOutlierRemoval.mat'), 'tableAfterOutlier', 'monthlyStatsFiltered', 'yearlyStatsFiltered');
end

%% Print Summary
fprintf('\n=== Time Series Summary ===\n');
fprintf('Total data points: %d\n', length(allR1s));
if any(oldDataIdx)
    fprintf('OldData points: %d\n', sum(oldDataIdx));
end
if any(newDataIdx)
    fprintf('NewData points: %d\n', sum(newDataIdx));
end
fprintf('Date range: %s to %s\n', datestr(min(allDates), 'yyyy-mm-dd'), datestr(max(allDates), 'yyyy-mm-dd'));
fprintf('Mean R1s: %.4f mΩ\n', mean(allR1s)*1000);
fprintf('Std R1s: %.4f mΩ\n', std(allR1s)*1000);
fprintf('Mean R²: %.3f\n', mean(allR2));

if length(allR1sTemp) > 0
    fprintf('\n=== After Outlier Removal ===\n');
    fprintf('Valid data points: %d\n', length(allR1sTemp));
    fprintf('Mean R1s: %.4f mΩ\n', mean(allR1sTemp)*1000);
    fprintf('Std R1s: %.4f mΩ\n', std(allR1sTemp)*1000);
    fprintf('Outliers removed: %d (%.1f%%)\n', length(allR1s)-length(allR1sTemp), ...
        (length(allR1s)-length(allR1sTemp))/length(allR1s)*100);
end

fprintf('\nResults saved to: %s\n', saveDir);
fprintf('Time series analysis completed!\n');
