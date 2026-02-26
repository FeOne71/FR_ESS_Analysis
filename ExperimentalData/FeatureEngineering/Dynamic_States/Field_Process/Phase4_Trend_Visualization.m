% Phase4_Trend_Visualization.m
% Visualizes the long-term trends of the 3 Key Labels: 
% Energy Efficiency, Standard IR, and SOH_BMS over time
% Saves plots as .fig and .pdf to the Field_Process/Results folder

clear; clc; close all;
warning('off', 'all');

%% Configuration
fieldProcessDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Field_Process';
resultsDir = fullfile(fieldProcessDir, 'Results');

inFile = fullfile(resultsDir, 'Daily_Features_Labels.mat');
if ~exist(inFile, 'file')
    error('Daily_Features_Labels.mat not found. Run Phase 3 first.');
end
load(inFile, 'DailyTable');

% Format dates for plotting
% Date strings are 'yyyyMMdd'
dates_dt = datetime(DailyTable.Date, 'InputFormat', 'yyyyMMdd');
[dates_dt, sortIdx] = sort(dates_dt);
DailyTable = DailyTable(sortIdx, :);

%% Visualization: Trend Plot (3 Subplots for 3 Labels)
fig = figure('Name', 'Long-term Trend Analysis', 'Position', [150, 150, 1000, 800], 'Color', 'w');

% Plot 1: Energy Efficiency
subplot(3,1,1);
plot(dates_dt, DailyTable.Label_Energy_Eff, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'Color', [0, 0.4470, 0.7410]);
title('Trend of Energy Efficiency (%)');
ylabel('Efficiency (%)');
grid on;
datetick('x', 'yyyy-mm-dd', 'keepticks');
xlim([min(dates_dt)-days(15), max(dates_dt)+days(15)]);

% Plot 2: Standard IR (5 levels)
subplot(3,1,2);
hold on;
plot(dates_dt, DailyTable.Label_Std_IR_1s, '-s', 'LineWidth', 1, 'MarkerFaceColor', 'r', 'Color', [0.8500, 0.3250, 0.0980]);
plot(dates_dt, DailyTable.Label_Std_IR_3s, '-o', 'LineWidth', 1, 'MarkerFaceColor', 'm', 'Color', [1, 0, 1]);
plot(dates_dt, DailyTable.Label_Std_IR_5s, '-^', 'LineWidth', 1, 'MarkerFaceColor', 'k', 'Color', [0, 0, 0]);
plot(dates_dt, DailyTable.Label_Std_IR_10s, '-d', 'LineWidth', 1, 'MarkerFaceColor', 'c', 'Color', [0, 1, 1]);
plot(dates_dt, DailyTable.Label_Std_IR_30s, '-v', 'LineWidth', 1, 'MarkerFaceColor', 'b', 'Color', [0, 0, 1]);
hold off;
title('Trend of Standard Internal Resistance (5 levels: 1s, 3s, 5s, 10s, 30s)');
ylabel('IR (\Omega)');
legend('R_{1s}', 'R_{3s}', 'R_{5s}', 'R_{10s}', 'R_{30s}', 'Location', 'best');
grid on;
datetick('x', 'yyyy-mm-dd', 'keepticks');
xlim([min(dates_dt)-days(15), max(dates_dt)+days(15)]);

% Plot 3: SOH BMS
subplot(3,1,3);
plot(dates_dt, DailyTable.Label_SOH_BMS, '-^', 'LineWidth', 1.5, 'MarkerFaceColor', 'g', 'Color', [0.4660, 0.6740, 0.1880]);
title('Trend of SOH (BMS Recording)');
ylabel('SOH (%)');
xlabel('Date');
grid on;
datetick('x', 'yyyy-mm-dd', 'keepticks');
xlim([min(dates_dt)-days(15), max(dates_dt)+days(15)]);

sgtitle('Long-term Degradation Trends in Field Data');

%% Save Figures
% User explicitly requested to save as .fig and NOT as .png
figName = fullfile(resultsDir, 'Trend_Analysis_Results');

saveas(fig, [figName, '.fig']);
% Print to PDF for reporting
% print(fig, [figName, '.pdf'], '-dpdf', '-bestfit');

fprintf('Trend visualization complete.\n');
% fprintf('Saved figures to %s as .fig and .pdf\n', resultsDir);
