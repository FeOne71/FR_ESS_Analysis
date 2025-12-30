%% Field R1s Arrhenius Temperature Correction Analysis - Version 3
% This script performs Arrhenius temperature correction on 1-second R1s materials
% with statistical outlier removal (IQR method) and creates selected figures
% 
% Process:
% 1. Load 1-second R1s materials from OldData_New and NewData_New
% 2. Apply statistical outlier removal using IQR method
% 3. Create master Arrhenius model from filtered data
% 4. Apply temperature correction to each 1-second data point
% 5. Generate analysis and visualizations (Fig1,2,3,6 only)

clear; clc; close all;

%% Configuration
fontSize = 12;
saveDir = 'Arrhenius_Analysis_Results_ver03';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Arrhenius parameters
Tref = 25 + 273.15;  % Reference temperature: 25°C in Kelvin

%% Load OldData_New and NewData_New results
fprintf('Loading OldData_New and NewData_New results for Arrhenius analysis...\n');

% OldData years
oldDataYears = {'2021', '2022', '2023'};
% NewData years  
newDataYears = {'2023', '2024', '2025'};

% Initialize data storage for master model creation
all_R_inst_master = [];
all_T_inst_master = [];

% Initialize debugging counters
debug_stats = struct();
debug_stats.files_loaded = 0;
debug_stats.files_failed = 0;
debug_stats.total_days_processed = 0;
debug_stats.total_data_points = 0;
debug_stats.start_time = tic;

%% Load OldData_New results for master model
fprintf('Loading OldData_New results for master model...\n');
for i = 1:length(oldDataYears)
    year = oldDataYears{i};
    fileName = sprintf('R1s_Results_%s.mat', year);
    filePath = fullfile('R1s_Results_OldData_ver02', fileName);
    
    if exist(filePath, 'file')
        try
            load(filePath, 'R1s_Results');
            debug_stats.files_loaded = debug_stats.files_loaded + 1;
            
            % Extract 1-second R1s materials
            dayNames = fieldnames(R1s_Results);
            year_data_points = 0;
            for d = 1:length(dayNames)
                dayKey = dayNames{d};
                if isfield(R1s_Results.(dayKey), 'Rack01')
                    rackData = R1s_Results.(dayKey).Rack01;
                    
                    % Extract R1s materials
                    if isfield(rackData, 'R_inst') && isfield(rackData, 'T_inst')
                        if ~isempty(rackData.R_inst) && ~isempty(rackData.T_inst)
                            % Data validation - only remove NaN values, keep all other data
                            valid_idx = ~isnan(rackData.R_inst) & ~isnan(rackData.T_inst);
                            
                            if sum(valid_idx) > 0
                                all_R_inst_master = [all_R_inst_master; rackData.R_inst(valid_idx)];
                                all_T_inst_master = [all_T_inst_master; rackData.T_inst(valid_idx)];
                                year_data_points = year_data_points + sum(valid_idx);
                                debug_stats.total_days_processed = debug_stats.total_days_processed + 1;
                            end
                        end
                    end
                end
            end
            
            debug_stats.total_data_points = debug_stats.total_data_points + year_data_points;
            
        catch ME
            fprintf('ERROR loading %s: %s\n', fileName, ME.message);
            debug_stats.files_failed = debug_stats.files_failed + 1;
        end
    else
        fprintf('WARNING: File not found: %s\n', fileName);
        debug_stats.files_failed = debug_stats.files_failed + 1;
    end
end

%% Load NewData_New results for master model
fprintf('Loading NewData_New results for master model...\n');
for i = 1:length(newDataYears)
    year = newDataYears{i};
    fileName = sprintf('R1s_Results_%s.mat', year);
    filePath = fullfile('R1s_Results_NewData_ver02', fileName);
    
    if exist(filePath, 'file')
        try
            load(filePath, 'R1s_Results');
            debug_stats.files_loaded = debug_stats.files_loaded + 1;
            
            % Extract 1-second R1s materials
            dayNames = fieldnames(R1s_Results);
            year_data_points = 0;
            for d = 1:length(dayNames)
                dayKey = dayNames{d};
                if isfield(R1s_Results.(dayKey), 'Rack01')
                    rackData = R1s_Results.(dayKey).Rack01;
                    
                    % Extract R1s materials
                    if isfield(rackData, 'R_inst') && isfield(rackData, 'T_inst')
                        if ~isempty(rackData.R_inst) && ~isempty(rackData.T_inst)
                            % Data validation - only remove NaN values, keep all other data
                            valid_idx = ~isnan(rackData.R_inst) & ~isnan(rackData.T_inst);
                            
                            if sum(valid_idx) > 0
                                all_R_inst_master = [all_R_inst_master; rackData.R_inst(valid_idx)];
                                all_T_inst_master = [all_T_inst_master; rackData.T_inst(valid_idx)];
                                year_data_points = year_data_points + sum(valid_idx);
                                debug_stats.total_days_processed = debug_stats.total_days_processed + 1;
                            end
                        end
                    end
                end
            end
            
            debug_stats.total_data_points = debug_stats.total_data_points + year_data_points;
            
        catch ME
            fprintf('ERROR loading %s: %s\n', fileName, ME.message);
            debug_stats.files_failed = debug_stats.files_failed + 1;
        end
    else
        fprintf('WARNING: File not found: %s\n', fileName);
        debug_stats.files_failed = debug_stats.files_failed + 1;
    end
end

%% Data Loading Summary
fprintf('\n=== DATA LOADING SUMMARY ===\n');
fprintf('Files loaded: %d, Failed: %d, Data points: %d, Loading time: %.1f sec\n', ...
    debug_stats.files_loaded, debug_stats.files_failed, debug_stats.total_data_points, toc(debug_stats.start_time));

if isempty(all_R_inst_master)
    error('No valid R1s data found! Check file paths and data structure.');
end

fprintf('Data range - R1s: %.4f to %.4f mΩ, Temperature: %.1f to %.1f°C\n', ...
    min(all_R_inst_master)*1000, max(all_R_inst_master)*1000, min(all_T_inst_master), max(all_T_inst_master));

%% Statistical Outlier Removal using Standard Z-Score Method
fprintf('\n=== Statistical Outlier Removal (Standard Z-Score Method) ===\n');
outlier_start_time = tic;

% Calculate Z-Score for R1s values
mean_val = mean(all_R_inst_master);
std_val = std(all_R_inst_master);
z_scores = abs((all_R_inst_master - mean_val) / std_val);

% Set Z-Score threshold (3-sigma rule)
z_threshold = 3.0;
outlier_idx = z_scores > z_threshold;

% Count outliers
n_outliers = sum(outlier_idx);
n_total = length(all_R_inst_master);
outlier_percentage = n_outliers / n_total * 100;

% Calculate bounds for display
lower_bound = mean_val - z_threshold * std_val;
upper_bound = mean_val + z_threshold * std_val;

fprintf('Standard Z-Score Method Results:\n');
fprintf('  Mean: %.6f mΩ, Std: %.6f mΩ, Z-threshold: %.1f\n', mean_val*1000, std_val*1000, z_threshold);
fprintf('  Lower bound: %.6f mΩ, Upper bound: %.6f mΩ\n', lower_bound*1000, upper_bound*1000);
fprintf('  Outliers removed: %d out of %d (%.1f%%)\n', n_outliers, n_total, outlier_percentage);

% Remove outliers
all_R_inst_master = all_R_inst_master(~outlier_idx);
all_T_inst_master = all_T_inst_master(~outlier_idx);

fprintf('Filtered data range - R1s: %.4f to %.4f mΩ, Temperature: %.1f to %.1f°C\n', ...
    min(all_R_inst_master)*1000, max(all_R_inst_master)*1000, min(all_T_inst_master), max(all_T_inst_master));
fprintf('Outlier removal completed in %.1f seconds\n', toc(outlier_start_time));

%% Master Arrhenius Model Creation
fprintf('\n=== Master Arrhenius Model Creation ===\n');
arrhenius_start_time = tic;

% Convert temperature to Kelvin
T_abs = all_T_inst_master + 273.15;  % Convert to Kelvin
x_arrhenius = 1./T_abs - 1/Tref;   % (1/T - 1/Tref)

% Method 1: Exponential form - R1s vs (1/T - 1/Tref)
try
    p_master = polyfit(x_arrhenius, all_R_inst_master, 1);
    R1s_Tref_master = p_master(2);  % Y-axis intercept (R1s at 25°C)
catch ME
    fprintf('ERROR in Method 1 fit: %s\n', ME.message);
    p_master = [0, 0];
    R1s_Tref_master = 0;
end

% Method 2: Logarithmic form - log(R1s) vs (1/T - 1/Tref) - COMMENTED OUT
% try
%     y_log = log(all_R_inst_master);
%     p_log_master = polyfit(x_arrhenius, y_log, 1);
%     R1s_Tref_log_master = exp(p_log_master(2));  % Y-axis intercept converted back
% catch ME
%     fprintf('ERROR in Method 2 fit: %s\n', ME.message);
%     p_log_master = [0, 0];
%     R1s_Tref_log_master = 0;
% end

fprintf('Arrhenius Model Results:\n');
fprintf('  Method 1 (Exponential): Slope = %.6f, R1s at 25°C = %.6f mΩ\n', p_master(1), R1s_Tref_master*1000);
% fprintf('  Method 2 (Logarithmic): Slope = %.6f, R1s at 25°C = %.6f mΩ\n', p_log_master(1), R1s_Tref_log_master*1000);
% fprintf('  Difference: %.6f mΩ (%.1f%%)\n', ...
%     abs(R1s_Tref_master - R1s_Tref_log_master)*1000, ...
%     abs(R1s_Tref_master - R1s_Tref_log_master)/mean([R1s_Tref_master, R1s_Tref_log_master])*100);

%% Master Arrhenius Regression Statistics
[~, ~, ~, ~, stats_exp] = regress(all_R_inst_master, [x_arrhenius, ones(size(x_arrhenius))]);
% [~, ~, ~, ~, stats_log] = regress(y_log, [x_arrhenius, ones(size(x_arrhenius))]);

fprintf('Regression Statistics:\n');
fprintf('  Method 1 - R²: %.4f, p-value: %.2e\n', stats_exp(1), stats_exp(3));
% fprintf('  Method 2 - R²: %.4f, p-value: %.2e\n', stats_log(1), stats_log(3));

%% 1-second Temperature Correction and Daily Aggregation
fprintf('\n=== 1-second Temperature Correction and Daily Aggregation ===\n');
correction_start_time = tic;

% Initialize final analysis data (1-second level)
final_dates = [];
final_R1s_uncorrected = [];
final_R1s_corrected_exp = [];
final_R1s_corrected_log = [];
final_Temp_inst = [];
final_source = [];

% Process OldData_New results
fprintf('Processing OldData_New results for final analysis...\n');
total_years = length(oldDataYears) + length(newDataYears);
processed_years = 0;

for i = 1:length(oldDataYears)
    year = oldDataYears{i};
    fileName = sprintf('R1s_Results_%s.mat', year);
    filePath = fullfile('R1s_Results_OldData_ver02', fileName);
    
    if exist(filePath, 'file')
        try
            load(filePath, 'R1s_Results');
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
                        
                        % Extract R1s materials
                        if isfield(rackData, 'R_inst') && isfield(rackData, 'T_inst')
                            if ~isempty(rackData.R_inst) && ~isempty(rackData.T_inst)
                                R_inst = rackData.R_inst;
                                T_inst = rackData.T_inst;
                                
                                % Apply same Z-Score filtering to individual data
                                valid_idx = ~isnan(R_inst) & ~isnan(T_inst);
                                R_inst = R_inst(valid_idx);
                                T_inst = T_inst(valid_idx);
                                
                                % Apply Z-Score outlier removal
                                if ~isempty(R_inst) && length(R_inst) > 1
                                    z_scores_inst = abs((R_inst - mean_val) / std_val);
                                    outlier_idx = z_scores_inst > z_threshold;
                                    R_inst = R_inst(~outlier_idx);
                                    T_inst = T_inst(~outlier_idx);
                                end
                                
                                if ~isempty(R_inst)
                                    % 1-second temperature correction (Method 1)
                                    x_arrhenius_inst = (1 ./ (T_inst + 273.15)) - (1 / Tref);
                                    R_predicted_inst = polyval(p_master, x_arrhenius_inst);
                                    correction_factor_exp = R1s_Tref_master ./ R_predicted_inst;
                                    R_inst_corrected_exp = R_inst .* correction_factor_exp;
                                    
                                    % 1-second temperature correction (Method 2) - COMMENTED OUT
                                    % R_predicted_inst_log = exp(polyval(p_log_master, x_arrhenius_inst));
                                    % correction_factor_log = R1s_Tref_log_master ./ R_predicted_inst_log;
                                    % R_inst_corrected_log = R_inst .* correction_factor_log;
                                    
                                    % Store 1-second data (not daily aggregated)
                                    n_points = length(R_inst);
                                    date_vector = repmat(date, n_points, 1);
                                    source_vector = repmat({'OldData'}, n_points, 1);
                                    
                                    final_dates = [final_dates; date_vector];
                                    final_R1s_uncorrected = [final_R1s_uncorrected; R_inst];
                                    final_R1s_corrected_exp = [final_R1s_corrected_exp; R_inst_corrected_exp];
                                    % final_R1s_corrected_log = [final_R1s_corrected_log; R_inst_corrected_log];
                                    final_Temp_inst = [final_Temp_inst; T_inst];
                                    final_source = [final_source; source_vector];
                                end
                            end
                        end
                    end
                end
            end
            processed_years = processed_years + 1;
            
        catch ME
            fprintf('  ERROR processing %s: %s\n', fileName, ME.message);
        end
    else
        fprintf('  WARNING: File not found: %s\n', fileName);
    end
end

% Process NewData_New results
fprintf('Processing NewData_New results for final analysis...\n');
for i = 1:length(newDataYears)
    year = newDataYears{i};
    fileName = sprintf('R1s_Results_%s.mat', year);
    filePath = fullfile('R1s_Results_NewData_ver02', fileName);
    
    if exist(filePath, 'file')
        try
            load(filePath, 'R1s_Results');
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
                        
                        % Extract R1s materials
                        if isfield(rackData, 'R_inst') && isfield(rackData, 'T_inst')
                            if ~isempty(rackData.R_inst) && ~isempty(rackData.T_inst)
                                R_inst = rackData.R_inst;
                                T_inst = rackData.T_inst;
                                
                                % Apply same Z-Score filtering to individual data
                                valid_idx = ~isnan(R_inst) & ~isnan(T_inst);
                                R_inst = R_inst(valid_idx);
                                T_inst = T_inst(valid_idx);
                                
                                % Apply Z-Score outlier removal
                                if ~isempty(R_inst) && length(R_inst) > 1
                                    z_scores_inst = abs((R_inst - mean_val) / std_val);
                                    outlier_idx = z_scores_inst > z_threshold;
                                    R_inst = R_inst(~outlier_idx);
                                    T_inst = T_inst(~outlier_idx);
                                end
                                
                                if ~isempty(R_inst)
                                    % 1-second temperature correction (Method 1)
                                    x_arrhenius_inst = (1 ./ (T_inst + 273.15)) - (1 / Tref);
                                    R_predicted_inst = polyval(p_master, x_arrhenius_inst);
                                    correction_factor_exp = R1s_Tref_master ./ R_predicted_inst;
                                    R_inst_corrected_exp = R_inst .* correction_factor_exp;
                                    
                                    % 1-second temperature correction (Method 2) - COMMENTED OUT
                                    % R_predicted_inst_log = exp(polyval(p_log_master, x_arrhenius_inst));
                                    % correction_factor_log = R1s_Tref_log_master ./ R_predicted_inst_log;
                                    % R_inst_corrected_log = R_inst .* correction_factor_log;
                                    
                                    % Store 1-second data (not daily aggregated)
                                    n_points = length(R_inst);
                                    date_vector = repmat(date, n_points, 1);
                                    source_vector = repmat({'NewData'}, n_points, 1);
                                    
                                    final_dates = [final_dates; date_vector];
                                    final_R1s_uncorrected = [final_R1s_uncorrected; R_inst];
                                    final_R1s_corrected_exp = [final_R1s_corrected_exp; R_inst_corrected_exp];
                                    % final_R1s_corrected_log = [final_R1s_corrected_log; R_inst_corrected_log];
                                    final_Temp_inst = [final_Temp_inst; T_inst];
                                    final_source = [final_source; source_vector];
                                end
                            end
                        end
                    end
                end
            end
            processed_years = processed_years + 1;
            
        catch ME
            fprintf('  ERROR processing %s: %s\n', fileName, ME.message);
        end
    else
        fprintf('  WARNING: File not found: %s\n', fileName);
    end
end

%% Sort by Date and Final Data Summary
[final_dates, sortIdx] = sort(final_dates);
final_R1s_uncorrected = final_R1s_uncorrected(sortIdx);
final_R1s_corrected_exp = final_R1s_corrected_exp(sortIdx);
% final_R1s_corrected_log = final_R1s_corrected_log(sortIdx);
final_Temp_inst = final_Temp_inst(sortIdx);
final_source = final_source(sortIdx);

% Create sequential x-axis (1, 2, 3, ...) instead of actual dates
xAxis = (1:length(final_dates))';

fprintf('\n=== FINAL ANALYSIS DATA SUMMARY ===\n');
fprintf('Processing completed in %.1f seconds, Data points: %d\n', toc(correction_start_time), length(final_dates));
fprintf('Date range: %s to %s\n', datestr(min(final_dates), 'yyyy-mm-dd'), datestr(max(final_dates), 'yyyy-mm-dd'));

% Data source breakdown
oldData_count = sum(strcmp(final_source, 'OldData'));
newData_count = sum(strcmp(final_source, 'NewData'));
fprintf('Data source: OldData %d (%.1f%%), NewData %d (%.1f%%)\n', ...
    oldData_count, oldData_count/length(final_source)*100, newData_count, newData_count/length(final_source)*100);

%% Temperature Correction Impact Analysis
fprintf('\n=== Temperature Correction Impact Analysis ===\n');

% Original R1s trend
mdl_original = fitlm(xAxis, final_R1s_uncorrected*1000);
slope_original = mdl_original.Coefficients.Estimate(2);
r2_original = mdl_original.Rsquared.Ordinary;
pval_original = mdl_original.Coefficients.pValue(2);

% Temperature corrected R1s trend (Method 1)
mdl_corrected_exp = fitlm(xAxis, final_R1s_corrected_exp*1000);
slope_corrected_exp = mdl_corrected_exp.Coefficients.Estimate(2);
r2_corrected_exp = mdl_corrected_exp.Rsquared.Ordinary;
pval_corrected_exp = mdl_corrected_exp.Coefficients.pValue(2);

% Temperature corrected R1s trend (Method 2) - COMMENTED OUT
% mdl_corrected_log = fitlm(xAxis, final_R1s_corrected_log*1000);
% slope_corrected_log = mdl_corrected_log.Coefficients.Estimate(2);
% r2_corrected_log = mdl_corrected_log.Rsquared.Ordinary;
% pval_corrected_log = mdl_corrected_log.Coefficients.pValue(2);

% Calculate slope change percentage
slope_change_exp = (slope_corrected_exp - slope_original) / abs(slope_original) * 100;
% slope_change_log = (slope_corrected_log - slope_original) / abs(slope_original) * 100;

fprintf('Trend Analysis:\n');
fprintf('  Original: %.6f mΩ/day (R²=%.4f)\n', slope_original, r2_original);
fprintf('  Corrected M1: %.6f mΩ/day (R²=%.4f, %.1f%% change)\n', slope_corrected_exp, r2_corrected_exp, slope_change_exp);
% fprintf('  Corrected M2: %.6f mΩ/day (R²=%.4f, %.1f%% change)\n', slope_corrected_log, r2_corrected_log, slope_change_log);

%% Figure 1: Original R1s Time Series (TimeSeries format)
fprintf('\nCreating figures...\n');

figure('Position', [100, 100, 1400, 800]);

% Plot all data as connected line with battery temperature color
scatter(xAxis, final_R1s_uncorrected*1000, 50, final_Temp_inst, 'filled');
hold on;
plot(xAxis, final_R1s_uncorrected*1000, '-', 'LineWidth', 1, 'Color', '#b0b0b0');

ylabel('R1s (mΩ)', 'FontSize', fontSize);

% Add trend line using linear regression
if length(final_R1s_uncorrected) > 1
    mdl = fitlm(xAxis, final_R1s_uncorrected*1000, 'Intercept', true);
    trendLine = mdl.Fitted;
    plot(xAxis, trendLine, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Trend (R²=%.4f)', mdl.Rsquared.Ordinary));
end

% Colorbar for battery temperature
ax = gca;
colormap(ax, flipud(autumn));
c = colorbar(ax, 'eastoutside');
c.Label.String = 'Battery Temperature (°C)';
c.Label.FontSize = fontSize;
% Set colorbar limits to actual data range
c.Limits = [min(final_Temp_inst), max(final_Temp_inst)];
c.Ticks = round(min(final_Temp_inst)):1:round(max(final_Temp_inst));
% Set colorbar limits to actual data range
c.Limits = [min(final_Temp_inst), max(final_Temp_inst)];
c.Ticks = round(min(final_Temp_inst)):1:round(max(final_Temp_inst));

% Formatting
xlabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
title('ESS Rack01 Daily R1s Time Series (Original - 1-second Materials, Z-Score Filtered)', 'FontSize', fontSize+2, 'Color', 'k');
grid on;
legend('Location', 'best', 'FontSize', fontSize);

% Set all axis colors to black
set(gca, 'XColor', 'k', 'YColor', 'k');

% Format x-axis - show only one label per month
% Find unique months and their first occurrence
dateYears = zeros(size(final_dates));
dateMonths = zeros(size(final_dates));
for i = 1:length(final_dates)
    dateYears(i) = final_dates(i).Year;
    dateMonths(i) = final_dates(i).Month;
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
saveas(gcf, fullfile(saveDir, 'Original_R1s_TimeSeries_ZScore_Filtered.fig'));

%% Figure 2: Temperature Corrected R1s Time Series (Method 1)

figure('Position', [100, 100, 1400, 800]);

% Plot all data as connected line with battery temperature color
scatter(xAxis, final_R1s_corrected_exp*1000, 50, final_Temp_inst, 'filled');
hold on;
plot(xAxis, final_R1s_corrected_exp*1000, '-', 'LineWidth', 1, 'Color', '#b0b0b0');

ylabel('R1s (mΩ) - Temperature Corrected to 25°C', 'FontSize', fontSize);

% Add trend line using linear regression
if length(final_R1s_corrected_exp) > 1
    mdl = fitlm(xAxis, final_R1s_corrected_exp*1000, 'Intercept', true);
    trendLine = mdl.Fitted;
    plot(xAxis, trendLine, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Trend (R²=%.4f)', mdl.Rsquared.Ordinary));
end

% Colorbar for battery temperature
ax = gca;
colormap(ax, flipud(autumn));
c = colorbar(ax, 'eastoutside');
c.Label.String = 'Battery Temperature (°C)';
c.Label.FontSize = fontSize;
% Set colorbar limits to actual data range
c.Limits = [min(final_Temp_inst), max(final_Temp_inst)];
c.Ticks = round(min(final_Temp_inst)):1:round(max(final_Temp_inst));
% Set colorbar limits to actual data range
c.Limits = [min(final_Temp_inst), max(final_Temp_inst)];
c.Ticks = round(min(final_Temp_inst)):1:round(max(final_Temp_inst));

% Formatting
xlabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
title('ESS Rack01 Daily R1s Time Series (Temperature Corrected - Method 1, Z-Score Filtered)', 'FontSize', fontSize+2, 'Color', 'k');
grid on;
legend('Location', 'best', 'FontSize', fontSize);

% Set all axis colors to black
set(gca, 'XColor', 'k', 'YColor', 'k');

% Format x-axis - show only one label per month
xticks(monthPositions);
xticklabels(monthLabels);
xtickangle(45);

% Save figure
saveas(gcf, fullfile(saveDir, 'Temperature_Corrected_R1s_TimeSeries_Method1_ZScore_Filtered.fig'));

%% Figure 3: Temperature Corrected R1s Time Series (Method 2) - COMMENTED OUT

% figure('Position', [100, 100, 1400, 800]);
% 
% % Plot all data as connected line with battery temperature color
% scatter(xAxis, final_R1s_corrected_log*1000, 50, final_Temp_inst, 'filled');
% hold on;
% plot(xAxis, final_R1s_corrected_log*1000, '-', 'LineWidth', 1, 'Color', '#b0b0b0');
% 
% ylabel('R1s (mΩ) - Temperature Corrected to 25°C', 'FontSize', fontSize);
% 
% % Add trend line using linear regression
% if length(final_R1s_corrected_log) > 1
%     mdl = fitlm(xAxis, final_R1s_corrected_log*1000, 'Intercept', true);
%     trendLine = mdl.Fitted;
%     plot(xAxis, trendLine, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Trend (R²=%.4f)', mdl.Rsquared.Ordinary));
% end
% 
% % Colorbar for battery temperature
% ax = gca;
% colormap(ax, flipud(autumn));
% c = colorbar(ax, 'eastoutside');
% c.Label.String = 'Battery Temperature (°C)';
% c.Label.FontSize = fontSize;
% Set colorbar limits to actual data range
c.Limits = [min(final_Temp_inst), max(final_Temp_inst)];
c.Ticks = round(min(final_Temp_inst)):1:round(max(final_Temp_inst));
% % Set colorbar limits to actual data range
c.Limits = [min(final_Temp_inst), max(final_Temp_inst)];
c.Ticks = round(min(final_Temp_inst)):1:round(max(final_Temp_inst));
% 
% % Formatting
% xlabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
% title('ESS Rack01 Daily R1s Time Series (Temperature Corrected - Method 2, IQR Filtered)', 'FontSize', fontSize+2, 'Color', 'k');
% grid on;
% legend('Location', 'best', 'FontSize', fontSize);
% 
% % Set all axis colors to black
% set(gca, 'XColor', 'k', 'YColor', 'k');
% 
% % Format x-axis - show only one label per month
% xticks(monthPositions);
% xticklabels(monthLabels);
% xtickangle(45);
% 
% % Save figure
% saveas(gcf, fullfile(saveDir, 'Temperature_Corrected_R1s_TimeSeries_Method2_IQR_Filtered.fig'));

%% Figure 6: 3D Plot - Temperature Corrected R1s (Method 2)

figure('Position', [100, 100, 1200, 800]);

% Create 3D scatter plot with same color scheme as TimeSeries
scatter3(final_Temp_inst, xAxis, final_R1s_corrected_exp*1000, 50, final_Temp_inst, 'filled');
hold on;

% Add 3D trend line
% Fit 3D trend: R1s = a*Temp + b*Date + c
X_matrix = [final_Temp_inst, xAxis, ones(length(final_Temp_inst), 1)];
coefficients = X_matrix \ (final_R1s_corrected_exp*1000);

% Create trend line points
temp_trend = linspace(min(final_Temp_inst), max(final_Temp_inst), 50);
date_trend = linspace(min(xAxis), max(xAxis), 50);
R1s_trend_line = coefficients(1)*temp_trend + coefficients(2)*date_trend + coefficients(3);

% Plot trend line
plot3(temp_trend, date_trend, R1s_trend_line, 'b-', 'LineWidth', 2);


% Apply same colormap as TimeSeries
colormap(flipud(autumn));

% Formatting
xlabel('Battery Temperature (°C)', 'FontSize', fontSize, 'Color', 'k');
ylabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
zlabel('R1s (mΩ) - Temperature Corrected to 25°C', 'FontSize', fontSize, 'Color', 'k');
title('3D Plot: Battery Temperature vs Date vs R1s (Temperature Corrected - Method 1, Z-Score Filtered)', 'FontSize', fontSize+2, 'Color', 'k');

% Set y-axis (Date) ticks and labels - same as TimeSeries
set(gca, 'YTick', monthPositions);
set(gca, 'YTickLabel', monthLabels);

% Add colorbar for battery temperature (match TimeSeries)
ax3d = gca;
colormap(ax3d, flipud(autumn));
c = colorbar(ax3d, 'eastoutside');
c.Label.String = 'Battery Temperature (°C)';
c.Label.FontSize = fontSize;
% Set colorbar limits to actual data range
c.Limits = [min(final_Temp_inst), max(final_Temp_inst)];
c.Ticks = round(min(final_Temp_inst)):1:round(max(final_Temp_inst));

% Add grid
grid on;

% Add legend with 3-sigma values
legend('Location', 'best', 'FontSize', fontSize-2);

% Set all axis colors to black
set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');

% Set view angle for better visualization
view(45, 30);

% Save 3D figure
saveas(gcf, fullfile(saveDir, 'R1s_3D_Plot_Temperature_Corrected_Method1_ZScore_Filtered.fig'));

%% Figure 7: Daily Average R1s 3D Plot (Method 1)

% Calculate daily averages from 1-second data
unique_dates = unique(final_dates);
daily_avg_dates = [];
daily_avg_R1s_uncorrected = [];
daily_avg_R1s_corrected = [];
daily_avg_Temp = [];
daily_avg_source = [];

for i = 1:length(unique_dates)
    date = unique_dates(i);
    date_idx = final_dates == date;
    
    if sum(date_idx) > 0
        daily_avg_dates = [daily_avg_dates; date];
        daily_avg_R1s_uncorrected = [daily_avg_R1s_uncorrected; mean(final_R1s_uncorrected(date_idx))];
        daily_avg_R1s_corrected = [daily_avg_R1s_corrected; mean(final_R1s_corrected_exp(date_idx))];
        daily_avg_Temp = [daily_avg_Temp; mean(final_Temp_inst(date_idx))];
        
        % Determine source (OldData or NewData)
        source_values = final_source(date_idx);
        if any(strcmp(source_values, 'OldData'))
            daily_avg_source = [daily_avg_source; {'OldData'}];
        else
            daily_avg_source = [daily_avg_source; {'NewData'}];
        end
    end
end

% Create sequential x-axis for daily data
daily_xAxis = (1:length(daily_avg_dates))';

figure('Position', [100, 100, 1200, 800]);

% Create 3D scatter plot with same color scheme as TimeSeries
% X-axis: Battery Temperature, Y-axis: Sequential Date Index, Z-axis: Daily Average R1s
scatter3(daily_avg_Temp, daily_xAxis, daily_avg_R1s_corrected*1000, 50, daily_avg_Temp, 'filled');
hold on;

% Add 3D trend line for daily averages
% Fit 3D trend: R1s = a*Temp + b*Date + c
X_matrix_daily = [daily_avg_Temp, daily_xAxis, ones(length(daily_avg_Temp), 1)];
coefficients_daily = X_matrix_daily \ (daily_avg_R1s_corrected*1000);

% Create trend line points
temp_trend_daily = linspace(min(daily_avg_Temp), max(daily_avg_Temp), 30);
date_trend_daily = linspace(min(daily_xAxis), max(daily_xAxis), 30);
R1s_trend_line_daily = coefficients_daily(1)*temp_trend_daily + coefficients_daily(2)*date_trend_daily + coefficients_daily(3);

% Plot trend line
plot3(temp_trend_daily, date_trend_daily, R1s_trend_line_daily, 'r-', 'LineWidth', 2);


% Apply same colormap as TimeSeries
colormap(flipud(autumn));

% Formatting
xlabel('Battery Temperature (°C)', 'FontSize', fontSize, 'Color', 'k');
ylabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
zlabel('Daily Average R1s (mΩ) - Temperature Corrected to 25°C', 'FontSize', fontSize, 'Color', 'k');
title('3D Plot: Daily Average R1s (Temperature Corrected - Method 1, IQR Filtered)', 'FontSize', fontSize+2, 'Color', 'k');

% Set axis limits and ticks for 3D plot
xlim([min(daily_avg_Temp), max(daily_avg_Temp)]);
ylim([min(daily_xAxis), max(daily_xAxis)]);
zlim([min(daily_avg_R1s_corrected)*1000, max(daily_avg_R1s_corrected)*1000]);

% Set y-axis (Date) ticks and labels - same as TimeSeries
% Find unique months and their first occurrence for daily data
dateYears = zeros(size(daily_avg_dates));
dateMonths = zeros(size(daily_avg_dates));
for i = 1:length(daily_avg_dates)
    dateYears(i) = daily_avg_dates(i).Year;
    dateMonths(i) = daily_avg_dates(i).Month;
end

uniqueMonths_daily = unique(dateYears*100 + dateMonths);
monthLabels_daily = {};
monthPositions_daily = [];

for i = 1:length(uniqueMonths_daily)
    monthVal = uniqueMonths_daily(i);
    yearVal = floor(monthVal/100);
    monthVal = mod(monthVal, 100);
    
    % Find first occurrence of this month
    monthIdx = find(dateYears == yearVal & dateMonths == monthVal, 1);
    if ~isempty(monthIdx)
        monthLabels_daily{end+1} = sprintf('%d-%02d', yearVal, monthVal);
        monthPositions_daily(end+1) = monthIdx;
    end
end

% Set y-axis (Date) ticks and labels
set(gca, 'YTick', monthPositions_daily);
set(gca, 'YTickLabel', monthLabels_daily);

% Add colorbar for battery temperature (match TimeSeries)
ax3d = gca;
colormap(ax3d, flipud(autumn));
c = colorbar(ax3d, 'eastoutside');
c.Label.String = 'Battery Temperature (°C)';
c.Label.FontSize = fontSize;
% Set colorbar limits to actual daily average temperature range
c.Limits = [min(daily_avg_Temp), max(daily_avg_Temp)];
c.Ticks = round(min(daily_avg_Temp)):1:round(max(daily_avg_Temp));

% Add grid
grid on;

% Add legend with 3-sigma values
legend('Location', 'best', 'FontSize', fontSize-2);

% Set all axis colors to black
set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');

% Set view angle for better visualization
view(45, 30);

% Save 3D figure
saveas(gcf, fullfile(saveDir, 'R1s_3D_Plot_Daily_Average_Method1_ZScore_Filtered.fig'));

%% Daily Average R1s Validity Analysis
fprintf('\n=== R1s (day-avg) vs. time 유효성 판단 ===\n');

% Calculate daily averages from 1-second data (reuse existing calculation)
% Note: daily_avg_dates, daily_avg_R1s_uncorrected, etc. are already calculated above

% Calculate daily_n_points for validity analysis first
daily_n_points = [];
for i = 1:length(daily_avg_dates)
    date = daily_avg_dates(i);
    date_idx = final_dates == date;
    daily_n_points = [daily_n_points; sum(date_idx)];
end

% Daily average validity analysis
fprintf('Daily Average R1s Validity Analysis (IQR Filtered):\n');
fprintf('  Total days: %d\n', length(daily_avg_dates));
fprintf('  Date range: %s to %s\n', datestr(min(daily_avg_dates), 'yyyy-mm-dd'), datestr(max(daily_avg_dates), 'yyyy-mm-dd'));

% Data points per day analysis
fprintf('\nData Points per Day Analysis:\n');
fprintf('  Mean points/day: %.1f\n', mean(daily_n_points));
fprintf('  Std points/day: %.1f\n', std(daily_n_points));
fprintf('  Min points/day: %d\n', min(daily_n_points));
fprintf('  Max points/day: %d\n', max(daily_n_points));

% Days with insufficient data
insufficient_days = daily_n_points < 10; % Less than 10 points per day
fprintf('  Days with <10 points: %d (%.1f%%)\n', sum(insufficient_days), sum(insufficient_days)/length(daily_n_points)*100);

% Daily R1s range analysis
fprintf('\nDaily R1s Range Analysis (IQR Filtered):\n');
fprintf('  Original - Mean: %.4f mΩ, Std: %.4f mΩ, Range: %.4f to %.4f mΩ\n', ...
    mean(daily_avg_R1s_uncorrected)*1000, std(daily_avg_R1s_uncorrected)*1000, ...
    min(daily_avg_R1s_uncorrected)*1000, max(daily_avg_R1s_uncorrected)*1000);
fprintf('  Corrected M1 - Mean: %.4f mΩ, Std: %.4f mΩ, Range: %.4f to %.4f mΩ\n', ...
    mean(daily_avg_R1s_corrected)*1000, std(daily_avg_R1s_corrected)*1000, ...
    min(daily_avg_R1s_corrected)*1000, max(daily_avg_R1s_corrected)*1000);

% Temperature range analysis
fprintf('\nDaily Temperature Range Analysis:\n');
fprintf('  Mean: %.1f°C, Std: %.1f°C, Range: %.1f to %.1f°C\n', ...
    mean(daily_avg_Temp), std(daily_avg_Temp), min(daily_avg_Temp), max(daily_avg_Temp));

% Year-by-year daily analysis
dateYears = zeros(size(daily_avg_dates));
for i = 1:length(daily_avg_dates)
    dateYears(i) = daily_avg_dates(i).Year;
end
years = unique(dateYears);

fprintf('\nYear-by-Year Daily Analysis (IQR Filtered):\n');
for i = 1:length(years)
    yearVal = years(i);
    yearIdx = dateYears == yearVal;
    yearR1s_orig = daily_avg_R1s_uncorrected(yearIdx);
    yearR1s_corr = daily_avg_R1s_corrected(yearIdx);
    yearTemp = daily_avg_Temp(yearIdx);
    yearPoints = daily_n_points(yearIdx);
    
    fprintf('Year %d: %d days, %.1f±%.1f points/day, %.1f±%.1f°C\n', ...
        yearVal, sum(yearIdx), mean(yearPoints), std(yearPoints), mean(yearTemp), std(yearTemp));
    fprintf('  Original: %.4f±%.4f mΩ, Corrected: %.4f±%.4f mΩ\n', ...
        mean(yearR1s_orig)*1000, std(yearR1s_orig)*1000, ...
        mean(yearR1s_corr)*1000, std(yearR1s_corr)*1000);
end

% Validity assessment
fprintf('\n=== 유효성 판단 결과 (IQR Filtered) ===\n');
valid_days = daily_n_points >= 10; % At least 10 points per day
valid_percentage = sum(valid_days) / length(daily_n_points) * 100;

fprintf('데이터 품질 평가:\n');
fprintf('  유효한 일수: %d/%d (%.1f%%)\n', sum(valid_days), length(daily_n_points), valid_percentage);
fprintf('  일평균 데이터 포인트: %.1f개\n', mean(daily_n_points));
fprintf('  온도 범위: %.1f~%.1f°C (적절함)\n', min(daily_avg_Temp), max(daily_avg_Temp));
fprintf('  IQR 필터링 적용: 이상치 제거됨\n');

if valid_percentage >= 80
    fprintf('  결론: 일평균 R1s 데이터 유효성 양호 (≥80%%)\n');
elseif valid_percentage >= 60
    fprintf('  결론: 일평균 R1s 데이터 유효성 보통 (60-80%%)\n');
else
    fprintf('  결론: 일평균 R1s 데이터 유효성 부족 (<60%%)\n');
end

%% Summary Statistics
fprintf('\n=== Summary Statistics ===\n');
fprintf('Original R1s - Mean: %.4f mΩ, Std: %.4f mΩ\n', mean(final_R1s_uncorrected)*1000, std(final_R1s_uncorrected)*1000);
fprintf('Corrected M1 - Mean: %.4f mΩ, Std: %.4f mΩ\n', mean(final_R1s_corrected_exp)*1000, std(final_R1s_corrected_exp)*1000);
% fprintf('Corrected M2 - Mean: %.4f mΩ, Std: %.4f mΩ\n', mean(final_R1s_corrected_log)*1000, std(final_R1s_corrected_log)*1000);

% Year-by-year summary
dateYears = zeros(size(final_dates));
for i = 1:length(final_dates)
    dateYears(i) = final_dates(i).Year;
end
years = unique(dateYears);
fprintf('\nYear-by-Year Summary:\n');
for i = 1:length(years)
    yearVal = years(i);
    yearIdx = dateYears == yearVal;
    yearR1s_orig = final_R1s_uncorrected(yearIdx);
    yearR1s_corr1 = final_R1s_corrected_exp(yearIdx);
    % yearR1s_corr2 = final_R1s_corrected_log(yearIdx);
    
    fprintf('Year %d: %d points, Original: %.4f±%.4f mΩ, M1: %.4f±%.4f mΩ\n', ...
        yearVal, sum(yearIdx), mean(yearR1s_orig)*1000, std(yearR1s_orig)*1000, ...
        mean(yearR1s_corr1)*1000, std(yearR1s_corr1)*1000);
end

%% Save results

% Save corrected data
arrhenius_results_ver03 = struct();
arrhenius_results_ver03.final_dates = final_dates;
arrhenius_results_ver03.final_R1s_uncorrected = final_R1s_uncorrected;
arrhenius_results_ver03.final_R1s_corrected_exp = final_R1s_corrected_exp;
% arrhenius_results_ver03.final_R1s_corrected_log = final_R1s_corrected_log;
arrhenius_results_ver03.final_Temp_inst = final_Temp_inst;
arrhenius_results_ver03.final_source = final_source;
arrhenius_results_ver03.R1s_Tref_master = R1s_Tref_master;
% arrhenius_results_ver03.R1s_Tref_log_master = R1s_Tref_log_master;
arrhenius_results_ver03.Tref = Tref;
arrhenius_results_ver03.p_master = p_master;
% arrhenius_results_ver03.p_log_master = p_log_master;
arrhenius_results_ver03.all_R_inst_master = all_R_inst_master;
arrhenius_results_ver03.all_T_inst_master = all_T_inst_master;
arrhenius_results_ver03.IQR_bounds = [lower_bound, upper_bound];
arrhenius_results_ver03.outlier_stats = [n_outliers, n_total, outlier_percentage];

save(fullfile(saveDir, 'Arrhenius_Analysis_Results_ZScore_Filtered.mat'), 'arrhenius_results_ver03');

fprintf('\nArrhenius analysis with Z-Score filtering complete! Results saved to: %s\n', saveDir);
