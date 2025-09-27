%% Field R1s Arrhenius Temperature Correction Analysis
% This script performs Arrhenius temperature correction on R1s data
% and creates the same figure format as the original SOC analysis
% 
% Arrhenius equation: R1s(T) = R1s_ref × exp(Ea/R × (1/T - 1/Tref))
% where Tref = 25°C = 298.15K

clear; clc; close all;

%% Configuration
fontSize = 12;
saveDir = 'Arrhenius_Analysis_Results';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Arrhenius parameters
Tref = 25 + 273.15;  % Reference temperature: 25°C in Kelvin
R_gas = 8.314;       % Gas constant (J/mol/K)

%% Load OldData and NewData results
fprintf('Loading OldData and NewData results for Arrhenius analysis...\n');

% OldData years
oldDataYears = {'2021', '2022', '2023'};
% NewData years  
newDataYears = {'2023', '2024', '2025'};

% Initialize data storage for daily R1s data
allDates = [];
allR1s = [];
allBatteryTemp = [];
allSource = []; % 'OldData' or 'NewData'

%% Load OldData results
fprintf('Loading OldData results...\n');
for i = 1:length(oldDataYears)
    year = oldDataYears{i};
    fileName = sprintf('R1s_Results_%s.mat', year);
    filePath = fullfile('R1s_Results_OldData', fileName);
    
    if exist(filePath, 'file')
        fprintf('Loading: %s\n', fileName);
        load(filePath, 'R1s_Results');
        
        % Extract daily R1s data
        dayNames = fieldnames(R1s_Results);
        for d = 1:length(dayNames)
            dayKey = dayNames{d};
            if isfield(R1s_Results.(dayKey), 'Rack01')
                rackData = R1s_Results.(dayKey).Rack01;
                
                % Extract date from dayKey (e.g., 'Raw_20210601' -> '2021-06-01')
                if length(dayKey) >= 12 && strcmp(dayKey(1:4), 'Raw_')
                    dateStr = dayKey(5:12); % '20210601'
                    year_num = str2double(dateStr(1:4));
                    month = str2double(dateStr(5:6));
                    day = str2double(dateStr(7:8));
                    
                    date = datetime(year_num, month, day);
                    
                    % Extract R1s and temperature data
                    if isfield(rackData, 'R1s') && isfield(rackData, 'BatteryTemp')
                        r1s = rackData.R1s;
                        batteryTemp = mean(rackData.BatteryTemp);
                        
                        if r1s > 0 && ~isnan(batteryTemp)  % Only positive R1s and valid temperature
                            allDates = [allDates; date];
                            allR1s = [allR1s; r1s];
                            allBatteryTemp = [allBatteryTemp; batteryTemp];
                            allSource = [allSource; {'OldData'}];
                        end
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
        
        % Extract daily R1s data
        dayNames = fieldnames(R1s_Results);
        for d = 1:length(dayNames)
            dayKey = dayNames{d};
            if isfield(R1s_Results.(dayKey), 'Rack01')
                rackData = R1s_Results.(dayKey).Rack01;
                
                % Extract date from dayKey (e.g., 'Raw_20230101' -> '2023-01-01')
                if length(dayKey) >= 12 && strcmp(dayKey(1:4), 'Raw_')
                    dateStr = dayKey(5:12); % '20230101'
                    year_num = str2double(dateStr(1:4));
                    month = str2double(dateStr(5:6));
                    day = str2double(dateStr(7:8));
                    
                    date = datetime(year_num, month, day);
                    
                    % Extract R1s and temperature data
                    if isfield(rackData, 'R1s') && isfield(rackData, 'BatteryTemp')
                        r1s = rackData.R1s;
                        batteryTemp = mean(rackData.BatteryTemp);
                        
                        if r1s > 0 && ~isnan(batteryTemp)  % Only positive R1s and valid temperature
                            allDates = [allDates; date];
                            allR1s = [allR1s; r1s];
                            allBatteryTemp = [allBatteryTemp; batteryTemp];
                            allSource = [allSource; {'NewData'}];
                        end
                    end
                end
            end
        end
    else
        fprintf('File not found: %s\n', fileName);
    end
end

%% Sort by Date
fprintf('Sorting data by date...\n');
[allDates, sortIdx] = sort(allDates);
allR1s = allR1s(sortIdx);
allBatteryTemp = allBatteryTemp(sortIdx);
allSource = allSource(sortIdx);

% Create sequential x-axis (1, 2, 3, ...) instead of actual dates
xAxis = (1:length(allDates))';

fprintf('\nData Summary:\n');
fprintf('Total data points: %d\n', length(allDates));
fprintf('Date range: %s to %s\n', datestr(min(allDates), 'yyyy-mm-dd'), datestr(max(allDates), 'yyyy-mm-dd'));
fprintf('Temperature range: %.1f°C to %.1f°C\n', min(allBatteryTemp), max(allBatteryTemp));
fprintf('R1s range: %.4f to %.4f mΩ\n', min(allR1s)*1000, max(allR1s)*1000);

%% Arrhenius Analysis
fprintf('\n=== Arrhenius Analysis ===\n');

% Convert temperature to Kelvin
T_abs = allBatteryTemp + 273.15;  % Convert to Kelvin
x_arrhenius = 1./T_abs - 1/Tref;   % (1/T - 1/Tref)

% Method 1: Exponential form - R1s vs (1/T - 1/Tref)
% Y-axis intercept gives R1s_Tref (R1s at 25°C)
p_exp = polyfit(x_arrhenius, allR1s, 1);
R1s_Tref_exp = p_exp(2);  % Y-axis intercept (R1s at 25°C)

% Method 2: Logarithmic form - log(R1s) vs (1/T - 1/Tref)  
% Y-axis intercept gives log(R1s_Tref)
y_log = log(allR1s);
p_log = polyfit(x_arrhenius, y_log, 1);
R1s_Tref_log = exp(p_log(2));  % Y-axis intercept converted back

fprintf('Method 1 (Exponential): R1s at 25°C = %.6f mΩ\n', R1s_Tref_exp*1000);
fprintf('Method 2 (Logarithmic): R1s at 25°C = %.6f mΩ\n', R1s_Tref_log*1000);
fprintf('Difference between methods: %.6f mΩ\n', abs(R1s_Tref_exp - R1s_Tref_log)*1000);

%% Temperature Correction
fprintf('\n=== Temperature Correction ===\n');

% Calculate temperature correction factors for each measurement
% Method 1: Using exponential form regression
% R1s_corrected = R1s_measured * (R1s_Tref / R1s_predicted_at_T)
% where R1s_predicted_at_T = p_exp(1) * (1/T - 1/Tref) + p_exp(2)

% Predict R1s at each temperature using the regression line
r1s_predicted_at_T = polyval(p_exp, x_arrhenius);

% Apply temperature correction: scale each measurement to 25°C
r1s_corrected_exp = allR1s .* (R1s_Tref_exp ./ r1s_predicted_at_T);

% Method 2: Using logarithmic form regression
% R1s_corrected = R1s_measured * (R1s_Tref / R1s_predicted_at_T)
% where R1s_predicted_at_T = exp(p_log(1) * (1/T - 1/Tref) + p_log(2))

% Predict R1s at each temperature using the logarithmic regression
r1s_predicted_at_T_log = exp(polyval(p_log, x_arrhenius));

% Apply temperature correction: scale each measurement to 25°C
r1s_corrected_log = allR1s .* (R1s_Tref_log ./ r1s_predicted_at_T_log);

fprintf('Original R1s mean: %.6f mΩ\n', mean(allR1s)*1000);
fprintf('Corrected R1s mean (Method 1): %.6f mΩ\n', mean(r1s_corrected_exp)*1000);
fprintf('Corrected R1s mean (Method 2): %.6f mΩ\n', mean(r1s_corrected_log)*1000);
fprintf('Temperature correction factor range (Method 1): %.4f to %.4f\n', ...
    min(R1s_Tref_exp ./ r1s_predicted_at_T), max(R1s_Tref_exp ./ r1s_predicted_at_T));
fprintf('Temperature correction factor range (Method 2): %.4f to %.4f\n', ...
    min(R1s_Tref_log ./ r1s_predicted_at_T_log), max(R1s_Tref_log ./ r1s_predicted_at_T_log));

%% Detailed Example of Temperature Correction
fprintf('\n=== Detailed Temperature Correction Example ===\n');
fprintf('Let''s take the first 5 data points as examples:\n\n');

for i = 1:min(5, length(allR1s))
    T_measured = allBatteryTemp(i);
    T_abs = T_measured + 273.15;  % Convert to Kelvin
    x_arrhenius_i = 1/T_abs - 1/Tref;  % (1/T - 1/Tref)
    
    % Method 1: Exponential form
    r1s_predicted_exp_i = polyval(p_exp, x_arrhenius_i);
    correction_factor_exp = R1s_Tref_exp / r1s_predicted_exp_i;
    r1s_corrected_exp_i = allR1s(i) * correction_factor_exp;
    
    % Method 2: Logarithmic form
    r1s_predicted_log_i = exp(polyval(p_log, x_arrhenius_i));
    correction_factor_log = R1s_Tref_log / r1s_predicted_log_i;
    r1s_corrected_log_i = allR1s(i) * correction_factor_log;
    
    fprintf('Data Point %d:\n', i);
    fprintf('  Date: %s\n', datestr(allDates(i), 'yyyy-mm-dd'));
    fprintf('  Measured Temperature: %.1f°C (%.1fK)\n', T_measured, T_abs);
    fprintf('  Measured R1s: %.4f mΩ\n', allR1s(i)*1000);
    fprintf('  Arrhenius x-value (1/T-1/Tref): %.6f K⁻¹\n', x_arrhenius_i);
    fprintf('  Method 1 - Predicted R1s at T: %.4f mΩ\n', r1s_predicted_exp_i*1000);
    fprintf('  Method 1 - Correction Factor: %.4f\n', correction_factor_exp);
    fprintf('  Method 1 - Corrected R1s: %.4f mΩ\n', r1s_corrected_exp_i*1000);
    fprintf('  Method 2 - Predicted R1s at T: %.4f mΩ\n', r1s_predicted_log_i*1000);
    fprintf('  Method 2 - Correction Factor: %.4f\n', correction_factor_log);
    fprintf('  Method 2 - Corrected R1s: %.4f mΩ\n', r1s_corrected_log_i*1000);
    fprintf('  Difference (M1-M2): %.4f mΩ\n\n', (r1s_corrected_exp_i - r1s_corrected_log_i)*1000);
end

fprintf('=== How Arrhenius Plot Works ===\n');
fprintf('1. We plot R1s vs (1/T - 1/Tref) where Tref = 25°C = 298.15K\n');
fprintf('2. The linear regression gives us: R1s = slope × (1/T - 1/Tref) + intercept\n');
fprintf('3. The Y-axis intercept (when 1/T - 1/Tref = 0) gives us R1s at 25°C\n');
fprintf('4. For temperature correction, we use the ratio: R1s_corrected = R1s_measured × (R1s_25°C / R1s_predicted_at_T)\n');
fprintf('5. This scales each measurement to what it would be at 25°C\n\n');

%% Figure 1: Original R1s Time Series (TimeSeries format)
fprintf('\nCreating Figure 1: Original R1s Time Series...\n');

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
    plot(xAxis, trendLine, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Trend (R²=%.4f)', mdl.Rsquared.Ordinary));
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
c.Ticks = 20:1:35; % Set ticks from 20 to 35 with 1°C intervals

% Formatting
xlabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
title('ESS Rack01 Daily R1s Time Series (Original)', 'FontSize', fontSize+2, 'Color', 'k');
grid on;
legend('Location', 'best', 'FontSize', fontSize);

% Set all axis colors to black
set(gca, 'XColor', 'k', 'YColor', 'k');
yyaxis left;
set(gca, 'YColor', 'k');

% Format x-axis - show only one label per month
% Find unique months and their first occurrence
dateYears = zeros(size(allDates));
dateMonths = zeros(size(allDates));
for i = 1:length(allDates)
    dateYears(i) = allDates(i).Year;
    dateMonths(i) = allDates(i).Month;
end

uniqueMonths = unique(dateYears*100 + dateMonths);
monthLabels = {};
monthPositions = [];

for i = 1:length(uniqueMonths)
    monthVal = uniqueMonths(i);
    yearVal = floor(monthVal/100);
    monthVal = mod(monthVal, 100);
    
    % Find first occurrence of this month
    monthIdx = find(dateYears == yearVal & dateMonths == monthVal, 1);
    if ~isempty(monthIdx)
        monthLabels{end+1} = sprintf('%d-%02d', yearVal, monthVal);
        monthPositions(end+1) = monthIdx;
    end
end

% Set x-axis labels
xticks(monthPositions);
xticklabels(monthLabels);
xtickangle(45);

% Save figure
saveas(gcf, fullfile(saveDir, 'Original_R1s_TimeSeries.fig'));
fprintf('Saved Figure 1: Original R1s Time Series\n');

%% Figure 2: Temperature Corrected R1s Time Series (Method 1)
fprintf('\nCreating Figure 2: Temperature Corrected R1s Time Series (Method 1)...\n');

figure('Position', [100, 100, 1400, 800]);

% Left y-axis: Temperature corrected R1s with battery temperature color
yyaxis left;

% Plot all data as connected line with battery temperature color
scatter(xAxis, r1s_corrected_exp*1000, 50, allBatteryTemp, 'filled');
hold on;
plot(xAxis, r1s_corrected_exp*1000, '-', 'LineWidth', 1, 'Color', '#b0b0b0');

ylabel('R1s (mΩ) - Temperature Corrected to 25°C', 'FontSize', fontSize);
ylim([0.5 1]);

% Add trend line using linear regression
if length(r1s_corrected_exp) > 1
    mdl = fitlm(xAxis, r1s_corrected_exp*1000, 'Intercept', true);
    trendLine = mdl.Fitted;
    plot(xAxis, trendLine, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Trend (R²=%.4f)', mdl.Rsquared.Ordinary));
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
c.Ticks = 20:1:35; % Set ticks from 20 to 35 with 1°C intervals

% Formatting
xlabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
title('ESS Rack01 Daily R1s Time Series (Temperature Corrected - Method 1)', 'FontSize', fontSize+2, 'Color', 'k');
grid on;
legend('Location', 'best', 'FontSize', fontSize);

% Set all axis colors to black
set(gca, 'XColor', 'k', 'YColor', 'k');
yyaxis left;
set(gca, 'YColor', 'k');

% Format x-axis - show only one label per month
xticks(monthPositions);
xticklabels(monthLabels);
xtickangle(45);

% Save figure
saveas(gcf, fullfile(saveDir, 'Temperature_Corrected_R1s_TimeSeries_Method1.fig'));
fprintf('Saved Figure 2: Temperature Corrected R1s Time Series (Method 1)\n');

%% Figure 3: Temperature Corrected R1s Time Series (Method 2)
fprintf('\nCreating Figure 3: Temperature Corrected R1s Time Series (Method 2)...\n');

figure('Position', [100, 100, 1400, 800]);

% Left y-axis: Temperature corrected R1s with battery temperature color
yyaxis left;

% Plot all data as connected line with battery temperature color
scatter(xAxis, r1s_corrected_log*1000, 50, allBatteryTemp, 'filled');
hold on;
plot(xAxis, r1s_corrected_log*1000, '-', 'LineWidth', 1, 'Color', '#b0b0b0');

ylabel('R1s (mΩ) - Temperature Corrected to 25°C', 'FontSize', fontSize);
ylim([0.5 1]);

% Add trend line using linear regression
if length(r1s_corrected_log) > 1
    mdl = fitlm(xAxis, r1s_corrected_log*1000, 'Intercept', true);
    trendLine = mdl.Fitted;
    plot(xAxis, trendLine, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Trend (R²=%.4f)', mdl.Rsquared.Ordinary));
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
c.Ticks = 20:1:35; % Set ticks from 20 to 35 with 1°C intervals

% Formatting
xlabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
title('ESS Rack01 Daily R1s Time Series (Temperature Corrected - Method 2)', 'FontSize', fontSize+2, 'Color', 'k');
grid on;
legend('Location', 'best', 'FontSize', fontSize);

% Set all axis colors to black
set(gca, 'XColor', 'k', 'YColor', 'k');
yyaxis left;
set(gca, 'YColor', 'k');

% Format x-axis - show only one label per month
xticks(monthPositions);
xticklabels(monthLabels);
xtickangle(45);

% Save figure
saveas(gcf, fullfile(saveDir, 'Temperature_Corrected_R1s_TimeSeries_Method2.fig'));
fprintf('Saved Figure 3: Temperature Corrected R1s Time Series (Method 2)\n');

%% Figure 4: Arrhenius Plots
fprintf('\nCreating Figure 4: Arrhenius Plots...\n');

figure('Position', [100, 100, 1400, 600]);

% Method 1: Exponential form - R1s vs (1/T - 1/Tref)
subplot(1, 2, 1);
scatter(x_arrhenius, allR1s*1000, 50, dateYears, 'filled');
hold on;

% Plot regression line
x_fit = linspace(min(x_arrhenius), max(x_arrhenius), 100);
y_fit_exp = polyval(p_exp, x_fit) * 1000;
plot(x_fit, y_fit_exp, 'r-', 'LineWidth', 2);

% Mark Y-axis intercept (R1s at 25°C)
plot(0, R1s_Tref_exp*1000, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red', ...
    'DisplayName', sprintf('R1s at 25°C = %.4f mΩ', R1s_Tref_exp*1000));

xlabel('1/T - 1/T_{ref} (K^{-1})', 'FontSize', fontSize);
ylabel('R1s (mΩ)', 'FontSize', fontSize);
title('Arrhenius Plot - Exponential Form', 'FontSize', fontSize);
colorbar;
colormap(jet);
c = colorbar;
c.Label.String = 'Year';
legend('Location', 'best');
grid on;

% Method 2: Logarithmic form - log(R1s) vs (1/T - 1/Tref)
subplot(1, 2, 2);
scatter(x_arrhenius, y_log, 50, dateYears, 'filled');
hold on;

% Plot regression line
y_fit_log = polyval(p_log, x_fit);
plot(x_fit, y_fit_log, 'r-', 'LineWidth', 2);

% Mark Y-axis intercept (log(R1s) at 25°C)
plot(0, p_log(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red', ...
    'DisplayName', sprintf('log(R1s) at 25°C = %.4f', p_log(2)));

xlabel('1/T - 1/T_{ref} (K^{-1})', 'FontSize', fontSize);
ylabel('log(R1s)', 'FontSize', fontSize);
title('Arrhenius Plot - Logarithmic Form', 'FontSize', fontSize);
colorbar;
colormap(jet);
c = colorbar;
c.Label.String = 'Year';
legend('Location', 'best');
grid on;

% Save figure
saveas(gcf, fullfile(saveDir, 'Arrhenius_Plots.fig'));
fprintf('Saved Figure 4: Arrhenius Plots\n');

%% Figure 5: 3D Plot - Temperature Corrected R1s (Method 1)
fprintf('\nCreating Figure 5: 3D Plot - Temperature Corrected R1s (Method 1)...\n');

figure('Position', [100, 100, 1200, 800]);

% Create 3D scatter plot with same color scheme as TimeSeries
% X-axis: Battery Temperature, Y-axis: Sequential Date Index, Z-axis: R1s
scatter3(allBatteryTemp, xAxis, r1s_corrected_exp*1000, 50, allBatteryTemp, 'filled');
hold on;

% Apply same colormap as TimeSeries
colormap(flipud(autumn));

% Formatting
xlabel('Battery Temperature (°C)', 'FontSize', fontSize, 'Color', 'k');
ylabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
zlabel('R1s (mΩ) - Temperature Corrected to 25°C', 'FontSize', fontSize, 'Color', 'k');
title('3D Plot: Battery Temperature vs Date vs R1s (Temperature Corrected - Method 1)', 'FontSize', fontSize+2, 'Color', 'k');

% Set axis limits and ticks for 3D plot
xlim([min(allBatteryTemp), max(allBatteryTemp)]);
ylim([min(xAxis), max(xAxis)]);
zlim([0.5 1]);

% Set y-axis (Date) ticks and labels - same as TimeSeries
set(gca, 'YTick', monthPositions);
set(gca, 'YTickLabel', monthLabels);

% Add colorbar for battery temperature (match TimeSeries)
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
saveas(gcf, fullfile(saveDir, 'R1s_3D_Plot_Temperature_Corrected_Method1.fig'));
fprintf('Saved Figure 5: 3D Plot - Temperature Corrected R1s (Method 1)\n');

%% Figure 6: 3D Plot - Temperature Corrected R1s (Method 2)
fprintf('\nCreating Figure 6: 3D Plot - Temperature Corrected R1s (Method 2)...\n');

figure('Position', [100, 100, 1200, 800]);

% Create 3D scatter plot with same color scheme as TimeSeries
scatter3(allBatteryTemp, xAxis, r1s_corrected_log*1000, 50, allBatteryTemp, 'filled');
hold on;

% Apply same colormap as TimeSeries
colormap(flipud(autumn));

% Formatting
xlabel('Battery Temperature (°C)', 'FontSize', fontSize, 'Color', 'k');
ylabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
zlabel('R1s (mΩ) - Temperature Corrected to 25°C', 'FontSize', fontSize, 'Color', 'k');
title('3D Plot: Battery Temperature vs Date vs R1s (Temperature Corrected - Method 2)', 'FontSize', fontSize+2, 'Color', 'k');

% Set axis limits and ticks for 3D plot
xlim([min(allBatteryTemp), max(allBatteryTemp)]);
ylim([min(xAxis), max(xAxis)]);
zlim([0.5 1]);

% Set y-axis (Date) ticks and labels - same as TimeSeries
set(gca, 'YTick', monthPositions);
set(gca, 'YTickLabel', monthLabels);

% Add colorbar for battery temperature (match TimeSeries)
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
saveas(gcf, fullfile(saveDir, 'R1s_3D_Plot_Temperature_Corrected_Method2.fig'));
fprintf('Saved Figure 6: 3D Plot - Temperature Corrected R1s (Method 2)\n');

%% Temperature Correction Impact Analysis
fprintf('\n=== Temperature Correction Impact Analysis ===\n');

% 1. Trend Analysis: Original vs Temperature Corrected R1s
fprintf('\n1. TREND ANALYSIS (Original vs Temperature Corrected):\n');
fprintf('====================================================\n');

% Original R1s trend
mdl_original = fitlm(xAxis, allR1s*1000);
slope_original = mdl_original.Coefficients.Estimate(2);
r2_original = mdl_original.Rsquared.Ordinary;
pval_original = mdl_original.Coefficients.pValue(2);

% Temperature corrected R1s trend (Method 1)
mdl_corrected_exp = fitlm(xAxis, r1s_corrected_exp*1000);
slope_corrected_exp = mdl_corrected_exp.Coefficients.Estimate(2);
r2_corrected_exp = mdl_corrected_exp.Rsquared.Ordinary;
pval_corrected_exp = mdl_corrected_exp.Coefficients.pValue(2);

% Temperature corrected R1s trend (Method 2)
mdl_corrected_log = fitlm(xAxis, r1s_corrected_log*1000);
slope_corrected_log = mdl_corrected_log.Coefficients.Estimate(2);
r2_corrected_log = mdl_corrected_log.Rsquared.Ordinary;
pval_corrected_log = mdl_corrected_log.Coefficients.pValue(2);

fprintf('Original R1s Trend:\n');
fprintf('  Slope: %.6f mΩ/day (%.4f mΩ/year)\n', slope_original, slope_original*365);
fprintf('  R²: %.4f\n', r2_original);
fprintf('  P-value: %.2e\n', pval_original);
if slope_original > 0
    fprintf('  Trend: INCREASING\n');
else
    fprintf('  Trend: DECREASING\n');
end

fprintf('\nTemperature Corrected R1s Trend (Method 1):\n');
fprintf('  Slope: %.6f mΩ/day (%.4f mΩ/year)\n', slope_corrected_exp, slope_corrected_exp*365);
fprintf('  R²: %.4f\n', r2_corrected_exp);
fprintf('  P-value: %.2e\n', pval_corrected_exp);
if slope_corrected_exp > 0
    fprintf('  Trend: INCREASING\n');
else
    fprintf('  Trend: DECREASING\n');
end

fprintf('\nTemperature Corrected R1s Trend (Method 2):\n');
fprintf('  Slope: %.6f mΩ/day (%.4f mΩ/year)\n', slope_corrected_log, slope_corrected_log*365);
fprintf('  R²: %.4f\n', r2_corrected_log);
fprintf('  P-value: %.2e\n', pval_corrected_log);
if slope_corrected_log > 0
    fprintf('  Trend: INCREASING\n');
else
    fprintf('  Trend: DECREASING\n');
end

% Calculate slope change percentage
slope_change_exp = (slope_corrected_exp - slope_original) / abs(slope_original) * 100;
slope_change_log = (slope_corrected_log - slope_original) / abs(slope_original) * 100;

fprintf('\nSlope Change After Temperature Correction:\n');
fprintf('  Method 1: %.1f%% change\n', slope_change_exp);
fprintf('  Method 2: %.1f%% change\n', slope_change_log);

% 2. Year-by-Year Change Rate Analysis
fprintf('\n\n2. YEAR-BY-YEAR CHANGE RATE ANALYSIS:\n');
fprintf('=====================================\n');

% Calculate year-by-year statistics
dateYears = zeros(size(allDates));
for i = 1:length(allDates)
    dateYears(i) = allDates(i).Year;
end
years = unique(dateYears);
years = sort(years);

fprintf('Year\tOriginal\tCorrected(M1)\tCorrected(M2)\tTemp_Effect\n');
fprintf('----\t--------\t-------------\t-------------\t----------\n');

yearly_original = zeros(size(years));
yearly_corrected_exp = zeros(size(years));
yearly_corrected_log = zeros(size(years));

for i = 1:length(years)
    yearVal = years(i);
    yearIdx = dateYears == yearVal;
    
    yearly_original(i) = mean(allR1s(yearIdx))*1000;
    yearly_corrected_exp(i) = mean(r1s_corrected_exp(yearIdx))*1000;
    yearly_corrected_log(i) = mean(r1s_corrected_log(yearIdx))*1000;
    
    temp_effect = yearly_original(i) - yearly_corrected_exp(i);
    
    fprintf('%d\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n', ...
        yearVal, yearly_original(i), yearly_corrected_exp(i), yearly_corrected_log(i), temp_effect);
end

% Calculate year-to-year change rates
fprintf('\nYear-to-Year Change Rate (mΩ/year):\n');
fprintf('Period\t\tOriginal\tCorrected(M1)\tCorrected(M2)\n');
fprintf('------\t\t--------\t-------------\t-------------\n');

for i = 2:length(years)
    period = sprintf('%d-%d', years(i-1), years(i));
    change_original = yearly_original(i) - yearly_original(i-1);
    change_corrected_exp = yearly_corrected_exp(i) - yearly_corrected_exp(i-1);
    change_corrected_log = yearly_corrected_log(i) - yearly_corrected_log(i-1);
    
    fprintf('%s\t\t%.4f\t\t%.4f\t\t%.4f\n', ...
        period, change_original, change_corrected_exp, change_corrected_log);
end

% Overall change rate (2021 to 2025)
if length(years) >= 2
    total_years = years(end) - years(1);
    total_change_original = yearly_original(end) - yearly_original(1);
    total_change_corrected_exp = yearly_corrected_exp(end) - yearly_corrected_exp(1);
    total_change_corrected_log = yearly_corrected_log(end) - yearly_corrected_log(1);
    
    annual_rate_original = total_change_original / total_years;
    annual_rate_corrected_exp = total_change_corrected_exp / total_years;
    annual_rate_corrected_log = total_change_corrected_log / total_years;
    
    fprintf('\nOverall Annual Change Rate (%d-%d):\n', years(1), years(end));
    fprintf('Original: %.4f mΩ/year\n', annual_rate_original);
    fprintf('Corrected (M1): %.4f mΩ/year\n', annual_rate_corrected_exp);
    fprintf('Corrected (M2): %.4f mΩ/year\n', annual_rate_corrected_log);
end

% 3. Conclusions
fprintf('\n\n3. CONCLUSIONS:\n');
fprintf('==============\n');

fprintf('\nA. Temperature Correction Impact:\n');
fprintf('   - Standard deviation reduced by %.1f%% (%.4f → %.4f mΩ)\n', ...
    (std(allR1s)*1000 - std(r1s_corrected_exp)*1000) / (std(allR1s)*1000) * 100, ...
    std(allR1s)*1000, std(r1s_corrected_exp)*1000);
fprintf('   - Temperature effect range: %.4f to %.4f mΩ\n', ...
    min(yearly_original - yearly_corrected_exp), max(yearly_original - yearly_corrected_exp));

fprintf('\nB. Trend Analysis Results:\n');
if slope_original > 0 && slope_corrected_exp > 0
    fprintf('   - R1s INCREASES over time in both original and corrected data\n');
    if abs(slope_corrected_exp) < abs(slope_original)
        fprintf('   - Temperature correction REDUCES the increasing trend\n');
    else
        fprintf('   - Temperature correction AMPLIFIES the increasing trend\n');
    end
elseif slope_original > 0 && slope_corrected_exp <= 0
    fprintf('   - Original data shows INCREASING trend\n');
    fprintf('   - Temperature correction REVERSES the trend to DECREASING\n');
elseif slope_original <= 0 && slope_corrected_exp > 0
    fprintf('   - Original data shows DECREASING trend\n');
    fprintf('   - Temperature correction REVERSES the trend to INCREASING\n');
else
    fprintf('   - Both original and corrected data show DECREASING trend\n');
end

fprintf('\nC. Battery Health Assessment:\n');
if slope_corrected_exp > 0
    fprintf('   - Temperature-corrected R1s shows INCREASING trend: %.4f mΩ/year\n', slope_corrected_exp*365);
    fprintf('   - This indicates PROGRESSIVE BATTERY DEGRADATION\n');
    fprintf('   - The increase is %.1f%% of the trend after removing temperature effects\n', ...
        abs(slope_corrected_exp/slope_original)*100);
else
    fprintf('   - Temperature-corrected R1s shows DECREASING or STABLE trend\n');
    fprintf('   - This suggests MINIMAL or NO PROGRESSIVE DEGRADATION\n');
end

fprintf('\nD. Temperature Effect Significance:\n');
temp_effect_magnitude = abs(mean(yearly_original - yearly_corrected_exp));
if temp_effect_magnitude > 0.01
    fprintf('   - Temperature has SIGNIFICANT effect on R1s (%.4f mΩ average)\n', temp_effect_magnitude);
    fprintf('   - Temperature correction is ESSENTIAL for accurate degradation analysis\n');
else
    fprintf('   - Temperature has MINIMAL effect on R1s (%.4f mΩ average)\n', temp_effect_magnitude);
    fprintf('   - Temperature correction has LIMITED impact on trend analysis\n');
end

%% Create summary statistics
fprintf('\n=== Summary Statistics (Temperature Corrected) ===\n');

% Daily R1s data summary
fprintf('\nDaily R1s Data Summary:\n');
fprintf('Original R1s - Mean: %.4f mΩ, Std: %.4f mΩ, Range: %.4f to %.4f mΩ\n', ...
    mean(allR1s)*1000, std(allR1s)*1000, min(allR1s)*1000, max(allR1s)*1000);

fprintf('Temperature Corrected R1s (Method 1) - Mean: %.4f mΩ, Std: %.4f mΩ\n', ...
    mean(r1s_corrected_exp)*1000, std(r1s_corrected_exp)*1000);

fprintf('Temperature Corrected R1s (Method 2) - Mean: %.4f mΩ, Std: %.4f mΩ\n', ...
    mean(r1s_corrected_log)*1000, std(r1s_corrected_log)*1000);

% Year-by-year summary
dateYears = zeros(size(allDates));
for i = 1:length(allDates)
    dateYears(i) = allDates(i).Year;
end
years = unique(dateYears);
fprintf('\nYear-by-Year Summary:\n');
for i = 1:length(years)
    yearVal = years(i);
    yearIdx = dateYears == yearVal;
    yearR1s_orig = allR1s(yearIdx);
    yearR1s_corr1 = r1s_corrected_exp(yearIdx);
    yearR1s_corr2 = r1s_corrected_log(yearIdx);
    
    fprintf('Year %d: %d data points\n', yearVal, sum(yearIdx));
    fprintf('  Original R1s: Mean=%.4f mΩ, Std=%.4f mΩ\n', mean(yearR1s_orig)*1000, std(yearR1s_orig)*1000);
    fprintf('  Corrected R1s (M1): Mean=%.4f mΩ, Std=%.4f mΩ\n', mean(yearR1s_corr1)*1000, std(yearR1s_corr1)*1000);
    fprintf('  Corrected R1s (M2): Mean=%.4f mΩ, Std=%.4f mΩ\n', mean(yearR1s_corr2)*1000, std(yearR1s_corr2)*1000);
end

%% Save results
fprintf('\nSaving analysis results...\n');

% Save corrected data
arrhenius_results = struct();
arrhenius_results.allDates = allDates;
arrhenius_results.allR1s = allR1s;
arrhenius_results.allBatteryTemp = allBatteryTemp;
arrhenius_results.allSource = allSource;
arrhenius_results.r1s_corrected_exp = r1s_corrected_exp;
arrhenius_results.r1s_corrected_log = r1s_corrected_log;
arrhenius_results.R1s_Tref_exp = R1s_Tref_exp;
arrhenius_results.R1s_Tref_log = R1s_Tref_log;
arrhenius_results.Tref = Tref;
arrhenius_results.x_arrhenius = x_arrhenius;
arrhenius_results.y_log = y_log;

save(fullfile(saveDir, 'Arrhenius_Analysis_Results.mat'), 'arrhenius_results');

fprintf('\nArrhenius analysis complete!\n');
fprintf('Results saved to: %s\n', saveDir);
