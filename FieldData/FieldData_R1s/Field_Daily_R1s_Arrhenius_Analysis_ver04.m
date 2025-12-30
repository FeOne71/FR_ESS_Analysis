%% Field R1s Arrhenius Temperature Correction Analysis - Version 4
% This script performs year-by-year Arrhenius temperature correction on 1-second R1s materials
% with individual models for each year and creates comprehensive analysis
% 
% Process:
% 1. Load 1-second R1s materials from OldData_New and NewData_New
% 2. Group data by year
% 3. For each year: (a) Remove outliers, (b) Create log-linear Arrhenius model, (c) Validate model
% 4. Apply year-specific temperature correction to each 1-second data point
% 5. Generate comprehensive analysis and visualizations

clear; clc; close all;

%% Configuration
fontSize = 12;
saveDir = 'Arrhenius_Analysis_Results_ver04';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Arrhenius parameters
Tref = 25 + 273.15;  % Reference temperature: 25°C in Kelvin

%% Load OldData_New and NewData_New results
fprintf('Loading OldData_New and NewData_New results for year-by-year Arrhenius analysis...\n');

% OldData years
oldDataYears = {'2021', '2022', '2023'};
% NewData years  
newDataYears = {'2023', '2024', '2025'};

% Initialize data storage for year-by-year analysis
yearly_data = struct();
all_years = [];

% Initialize debugging counters
debug_stats = struct();
debug_stats.files_loaded = 0;
debug_stats.files_failed = 0;
debug_stats.total_days_processed = 0;
debug_stats.total_data_points = 0;
debug_stats.start_time = tic;

%% Load OldData_New results
fprintf('Loading OldData_New results...\n');
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
            year_R_inst = [];
            year_T_inst = [];
            year_dates = [];
            
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
                                % Data validation - only remove NaN values
                                valid_idx = ~isnan(rackData.R_inst) & ~isnan(rackData.T_inst);
                                
                                if sum(valid_idx) > 0
                                    year_R_inst = [year_R_inst; rackData.R_inst(valid_idx)];
                                    year_T_inst = [year_T_inst; rackData.T_inst(valid_idx)];
                                    year_dates = [year_dates; repmat(date, sum(valid_idx), 1)];
                                    year_data_points = year_data_points + sum(valid_idx);
                                    debug_stats.total_days_processed = debug_stats.total_days_processed + 1;
                                end
                            end
                        end
                    end
                end
            end
            
            % Store yearly data
            if ~isempty(year_R_inst)
                yearly_data.(['year_' year]).R_inst = year_R_inst;
                yearly_data.(['year_' year]).T_inst = year_T_inst;
                yearly_data.(['year_' year]).dates = year_dates;
                yearly_data.(['year_' year]).source = 'OldData';
                yearly_data.(['year_' year]).year_str = year;
                all_years = [all_years, {year}];
                debug_stats.total_data_points = debug_stats.total_data_points + year_data_points;
            end
            
        catch ME
            fprintf('ERROR loading %s: %s\n', fileName, ME.message);
            debug_stats.files_failed = debug_stats.files_failed + 1;
        end
    else
        fprintf('WARNING: File not found: %s\n', fileName);
        debug_stats.files_failed = debug_stats.files_failed + 1;
    end
end

%% Load NewData_New results
fprintf('Loading NewData_New results...\n');
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
            year_R_inst = [];
            year_T_inst = [];
            year_dates = [];
            
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
                                % Data validation - only remove NaN values
                                valid_idx = ~isnan(rackData.R_inst) & ~isnan(rackData.T_inst);
                                
                                if sum(valid_idx) > 0
                                    year_R_inst = [year_R_inst; rackData.R_inst(valid_idx)];
                                    year_T_inst = [year_T_inst; rackData.T_inst(valid_idx)];
                                    year_dates = [year_dates; repmat(date, sum(valid_idx), 1)];
                                    year_data_points = year_data_points + sum(valid_idx);
                                    debug_stats.total_days_processed = debug_stats.total_days_processed + 1;
                                end
                            end
                        end
                    end
                end
            end
            
            % Store yearly data (merge with existing if year already exists)
            if ~isempty(year_R_inst)
                field_name = ['year_' year];
                if isfield(yearly_data, field_name)
                    % Merge with existing data
                    yearly_data.(field_name).R_inst = [yearly_data.(field_name).R_inst; year_R_inst];
                    yearly_data.(field_name).T_inst = [yearly_data.(field_name).T_inst; year_T_inst];
                    yearly_data.(field_name).dates = [yearly_data.(field_name).dates; year_dates];
                    yearly_data.(field_name).source = 'Mixed'; % Both OldData and NewData
                else
                    % New year
                    yearly_data.(field_name).R_inst = year_R_inst;
                    yearly_data.(field_name).T_inst = year_T_inst;
                    yearly_data.(field_name).dates = year_dates;
                    yearly_data.(field_name).source = 'NewData';
                    yearly_data.(field_name).year_str = year;
                    all_years = [all_years, {year}];
                end
                debug_stats.total_data_points = debug_stats.total_data_points + year_data_points;
            end
            
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

if isempty(all_years)
    error('No valid R1s data found! Check file paths and data structure.');
end

% Sort years
all_years = sort(all_years);
fprintf('Years with data: %s\n', strjoin(all_years, ', '));

%% Year-by-Year Arrhenius Model Creation
fprintf('\n=== Year-by-Year Arrhenius Model Creation ===\n');
arrhenius_start_time = tic;

% Initialize storage for yearly models
yearly_models = struct();
yearly_outlier_stats = struct();

for i = 1:length(all_years)
    year = all_years{i};
    fprintf('\n--- Processing Year %s ---\n', year);
    
    % Get yearly data
    field_name = ['year_' year];
    R_inst = yearly_data.(field_name).R_inst;
    T_inst = yearly_data.(field_name).T_inst;
    dates = yearly_data.(field_name).dates;
    
    fprintf('  Data points: %d, Date range: %s to %s\n', ...
        length(R_inst), datestr(min(dates), 'yyyy-mm-dd'), datestr(max(dates), 'yyyy-mm-dd'));
    fprintf('  R1s range: %.4f to %.4f mΩ, Temperature range: %.1f to %.1f°C\n', ...
        min(R_inst)*1000, max(R_inst)*1000, min(T_inst), max(T_inst));
    
    %% (a) Outlier Removal for this year
    fprintf('  Applying outlier removal...\n');
    
    % Calculate Z-Score for this year's data
    mean_val = mean(R_inst);
    std_val = std(R_inst);
    z_scores = abs((R_inst - mean_val) / std_val);
    
    % Set Z-Score threshold (3-sigma rule)
    z_threshold = 3.0;
    outlier_idx = z_scores > z_threshold;
    
    % Count outliers
    n_outliers = sum(outlier_idx);
    n_total = length(R_inst);
    outlier_percentage = n_outliers / n_total * 100;
    
    % Calculate bounds for display
    lower_bound = mean_val - z_threshold * std_val;
    upper_bound = mean_val + z_threshold * std_val;
    
    fprintf('    Mean: %.6f mΩ, Std: %.6f mΩ, Z-threshold: %.1f\n', mean_val*1000, std_val*1000, z_threshold);
    fprintf('    Lower bound: %.6f mΩ, Upper bound: %.6f mΩ\n', lower_bound*1000, upper_bound*1000);
    fprintf('    Outliers removed: %d out of %d (%.1f%%)\n', n_outliers, n_total, outlier_percentage);
    
    % Remove outliers
    R_inst_filtered = R_inst(~outlier_idx);
    T_inst_filtered = T_inst(~outlier_idx);
    dates_filtered = dates(~outlier_idx);
    
    fprintf('    Filtered data: %d points, R1s range: %.4f to %.4f mΩ\n', ...
        length(R_inst_filtered), min(R_inst_filtered)*1000, max(R_inst_filtered)*1000);
    
    % Store outlier statistics
    yearly_outlier_stats.(field_name) = struct();
    yearly_outlier_stats.(field_name).n_outliers = n_outliers;
    yearly_outlier_stats.(field_name).n_total = n_total;
    yearly_outlier_stats.(field_name).outlier_percentage = outlier_percentage;
    yearly_outlier_stats.(field_name).mean_val = mean_val;
    yearly_outlier_stats.(field_name).std_val = std_val;
    yearly_outlier_stats.(field_name).bounds = [lower_bound, upper_bound];
    
    %% (b) Create Log-Linear Arrhenius Model for this year
    fprintf('  Creating log-linear Arrhenius model...\n');
    
    if length(R_inst_filtered) < 10
        fprintf('    WARNING: Insufficient data points (%d) for reliable model\n', length(R_inst_filtered));
        continue;
    end
    
    % Convert temperature to Kelvin
    T_abs = T_inst_filtered + 273.15;  % Convert to Kelvin
    x_arrhenius = 1./T_abs - 1/Tref;   % (1/T - 1/Tref)
    
    % Log-linear form: ln(R1s) vs (1/T - 1/Tref)
    try
        y_log = log(R_inst_filtered);
        p_log = polyfit(x_arrhenius, y_log, 1);
        R1s_Tref_log = exp(p_log(2));  % Y-axis intercept converted back
        
        % Calculate regression statistics
        [~, ~, ~, ~, stats_log] = regress(y_log, [x_arrhenius, ones(size(x_arrhenius))]);
        
        fprintf('    Log-linear model: Slope = %.6f, R1s at 25°C = %.6f mΩ\n', p_log(1), R1s_Tref_log*1000);
        fprintf('    R²: %.4f, p-value: %.2e\n', stats_log(1), stats_log(3));
        
        % Store model
        yearly_models.(field_name) = struct();
        yearly_models.(field_name).p_log = p_log;
        yearly_models.(field_name).R1s_Tref_log = R1s_Tref_log;
        yearly_models.(field_name).stats_log = stats_log;
        yearly_models.(field_name).x_arrhenius = x_arrhenius;
        yearly_models.(field_name).y_log = y_log;
        yearly_models.(field_name).R_inst_filtered = R_inst_filtered;
        yearly_models.(field_name).T_inst_filtered = T_inst_filtered;
        yearly_models.(field_name).dates_filtered = dates_filtered;
        
    catch ME
        fprintf('    ERROR in log-linear fit: %s\n', ME.message);
        continue;
    end
    
    %% (c) Model Validation Visualization for this year
    fprintf('  Creating model validation plot...\n');
    
    figure('Position', [100, 100, 1200, 800]);
    
    % Subplot 1: Log-linear Arrhenius plot
    subplot(2,2,1);
    scatter(x_arrhenius, y_log, 50, T_inst_filtered, 'filled');
    hold on;
    
    % Plot fitted line
    x_fit = linspace(min(x_arrhenius), max(x_arrhenius), 100);
    y_fit = polyval(p_log, x_fit);
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    
    xlabel('1/T - 1/T_{ref} (K^{-1})', 'FontSize', fontSize);
    ylabel('ln(R1s)', 'FontSize', fontSize);
    title(sprintf('Year %s: Log-Linear Arrhenius Model (R²=%.4f)', year, stats_log(1)), 'FontSize', fontSize+2);
    colorbar;
    colormap(flipud(autumn));
    grid on;
    
    % Subplot 2: Temperature vs R1s (original)
    subplot(2,2,2);
    scatter(T_inst_filtered, R_inst_filtered*1000, 50, T_inst_filtered, 'filled');
    xlabel('Temperature (°C)', 'FontSize', fontSize);
    ylabel('R1s (mΩ)', 'FontSize', fontSize);
    title(sprintf('Year %s: Temperature vs R1s (Filtered)', year), 'FontSize', fontSize+2);
    colorbar;
    colormap(flipud(autumn));
    grid on;
    
    % Subplot 3: Model prediction vs actual
    subplot(2,2,3);
    R_predicted = exp(polyval(p_log, x_arrhenius));
    scatter(R_inst_filtered*1000, R_predicted*1000, 50, T_inst_filtered, 'filled');
    hold on;
    plot([min(R_inst_filtered)*1000, max(R_inst_filtered)*1000], ...
         [min(R_inst_filtered)*1000, max(R_inst_filtered)*1000], 'r--', 'LineWidth', 2);
    xlabel('Actual R1s (mΩ)', 'FontSize', fontSize);
    ylabel('Predicted R1s (mΩ)', 'FontSize', fontSize);
    title(sprintf('Year %s: Model Validation', year), 'FontSize', fontSize+2);
    colorbar;
    colormap(flipud(autumn));
    grid on;
    
    % Subplot 4: Residuals
    subplot(2,2,4);
    residuals = R_inst_filtered - R_predicted;
    scatter(T_inst_filtered, residuals*1000, 50, T_inst_filtered, 'filled');
    xlabel('Temperature (°C)', 'FontSize', fontSize);
    ylabel('Residuals (mΩ)', 'FontSize', fontSize);
    title(sprintf('Year %s: Residuals', year), 'FontSize', fontSize+2);
    colorbar;
    colormap(flipud(autumn));
    grid on;
    
    % Save validation figure
    saveas(gcf, fullfile(saveDir, sprintf('Year_%s_Model_Validation.fig', year)));
    
    fprintf('    Model validation plot saved\n');
end

fprintf('\nYear-by-year model creation completed in %.1f seconds\n', toc(arrhenius_start_time));

%% Year-Specific Temperature Correction and Final Data Assembly
fprintf('\n=== Year-Specific Temperature Correction and Final Data Assembly ===\n');
correction_start_time = tic;

% Initialize final analysis data (1-second level)
final_dates = [];
final_R1s_uncorrected = [];
final_R1s_corrected = [];
final_Temp_inst = [];
final_source = [];
final_year = [];

% Process each year with its specific model
for i = 1:length(all_years)
    year = all_years{i};
    field_name = ['year_' year];
    
    if ~isfield(yearly_models, field_name)
        fprintf('  WARNING: No model available for year %s, skipping\n', year);
        continue;
    end
    
    fprintf('  Processing year %s with year-specific model...\n', year);
    
    % Get yearly data and model
    R_inst = yearly_data.(field_name).R_inst;
    T_inst = yearly_data.(field_name).T_inst;
    dates = yearly_data.(field_name).dates;
    source = yearly_data.(field_name).source;
    
    p_log = yearly_models.(field_name).p_log;
    R1s_Tref_log = yearly_models.(field_name).R1s_Tref_log;
    
    % Apply same outlier removal as used for model creation
    mean_val = yearly_outlier_stats.(field_name).mean_val;
    std_val = yearly_outlier_stats.(field_name).std_val;
    z_threshold = 3.0;
    z_scores = abs((R_inst - mean_val) / std_val);
    outlier_idx = z_scores > z_threshold;
    
    R_inst_filtered = R_inst(~outlier_idx);
    T_inst_filtered = T_inst(~outlier_idx);
    dates_filtered = dates(~outlier_idx);
    
    if ~isempty(R_inst_filtered)
        % Apply year-specific temperature correction
        x_arrhenius_inst = (1 ./ (T_inst_filtered + 273.15)) - (1 / Tref);
        R_predicted_inst = exp(polyval(p_log, x_arrhenius_inst));
        correction_factor = R1s_Tref_log ./ R_predicted_inst;
        R_inst_corrected = R_inst_filtered .* correction_factor;
        
        % Store corrected data
        n_points = length(R_inst_filtered);
        year_vector = repmat({year}, n_points, 1);
        source_vector = repmat({source}, n_points, 1);
        
        final_dates = [final_dates; dates_filtered];
        final_R1s_uncorrected = [final_R1s_uncorrected; R_inst_filtered];
        final_R1s_corrected = [final_R1s_corrected; R_inst_corrected];
        final_Temp_inst = [final_Temp_inst; T_inst_filtered];
        final_source = [final_source; source_vector];
        final_year = [final_year; year_vector];
        
        fprintf('    Corrected %d data points\n', n_points);
    end
end

%% Sort by Date and Final Data Summary
[final_dates, sortIdx] = sort(final_dates);
final_R1s_uncorrected = final_R1s_uncorrected(sortIdx);
final_R1s_corrected = final_R1s_corrected(sortIdx);
final_Temp_inst = final_Temp_inst(sortIdx);
final_source = final_source(sortIdx);
final_year = final_year(sortIdx);

% Create sequential x-axis (1, 2, 3, ...) instead of actual dates
xAxis = (1:length(final_dates))';

fprintf('\n=== FINAL ANALYSIS DATA SUMMARY ===\n');
fprintf('Processing completed in %.1f seconds, Data points: %d\n', toc(correction_start_time), length(final_dates));
fprintf('Date range: %s to %s\n', datestr(min(final_dates), 'yyyy-mm-dd'), datestr(max(final_dates), 'yyyy-mm-dd'));

% Data source breakdown
oldData_count = sum(strcmp(final_source, 'OldData'));
newData_count = sum(strcmp(final_source, 'NewData'));
mixed_count = sum(strcmp(final_source, 'Mixed'));
fprintf('Data source: OldData %d (%.1f%%), NewData %d (%.1f%%), Mixed %d (%.1f%%)\n', ...
    oldData_count, oldData_count/length(final_source)*100, ...
    newData_count, newData_count/length(final_source)*100, ...
    mixed_count, mixed_count/length(final_source)*100);

% Year breakdown
for i = 1:length(all_years)
    year = all_years{i};
    year_count = sum(strcmp(final_year, year));
    if year_count > 0
        fprintf('Year %s: %d points (%.1f%%)\n', year, year_count, year_count/length(final_source)*100);
    end
end

%% Temperature Correction Impact Analysis

fprintf('\n=== Year-by-Year Degradation Trend Analysis ===\n');
for i = 1:length(all_years)
    year = all_years{i};
    year_idx = strcmp(final_year, year); % 해당 연도의 데이터 인덱스

    if sum(year_idx) > 10 % 데이터가 충분할 경우에만
        year_xAxis = xAxis(year_idx);
        year_R1s_corrected = final_R1s_corrected(year_idx);

        % 해당 연도 데이터만으로 추세선 분석
        mdl_year = fitlm(year_xAxis, year_R1s_corrected * 1000);
        slope_year = mdl_year.Coefficients.Estimate(2);
        r2_year = mdl_year.Rsquared.Ordinary;
        
        fprintf('Year %s Trend: Slope = %.6f mΩ/day, R² = %.4f\n', year, slope_year, r2_year);
    end
end

%% Figure 1: Original R1s Time Series
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
c.Limits = [min(final_Temp_inst), max(final_Temp_inst)];
c.Ticks = round(min(final_Temp_inst)):1:round(max(final_Temp_inst));

% Formatting
xlabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
title('ESS Rack01 Daily R1s Time Series (Original - Year-by-Year Filtered)', 'FontSize', fontSize+2, 'Color', 'k');
grid on;
legend('Location', 'best', 'FontSize', fontSize);

% Set all axis colors to black
set(gca, 'XColor', 'k', 'YColor', 'k');

% Format x-axis - show only one label per month
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
saveas(gcf, fullfile(saveDir, 'Original_R1s_TimeSeries_YearByYear.fig'));

%% Figure 2: Temperature Corrected R1s Time Series

figure('Position', [100, 100, 1400, 800]);

% Plot all data as connected line with battery temperature color
scatter(xAxis, final_R1s_corrected*1000, 50, final_Temp_inst, 'filled');
hold on;
plot(xAxis, final_R1s_corrected*1000, '-', 'LineWidth', 1, 'Color', '#b0b0b0');

ylabel('R1s (mΩ) - Temperature Corrected to 25°C', 'FontSize', fontSize);

% Add trend line using linear regression
if length(final_R1s_corrected) > 1
    mdl = fitlm(xAxis, final_R1s_corrected*1000, 'Intercept', true);
    trendLine = mdl.Fitted;
    plot(xAxis, trendLine, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Trend (R²=%.4f)', mdl.Rsquared.Ordinary));
end

% Colorbar for battery temperature
ax = gca;
colormap(ax, flipud(autumn));
c = colorbar(ax, 'eastoutside');
c.Label.String = 'Battery Temperature (°C)';
c.Label.FontSize = fontSize;
c.Limits = [min(final_Temp_inst), max(final_Temp_inst)];
c.Ticks = round(min(final_Temp_inst)):1:round(max(final_Temp_inst));

% Formatting
xlabel('Date (YYYY:MM:DD)', 'FontSize', fontSize, 'Color', 'k');
title('ESS Rack01 Daily R1s Time Series (Year-by-Year Temperature Corrected)', 'FontSize', fontSize+2, 'Color', 'k');
grid on;
legend('Location', 'best', 'FontSize', fontSize);

% Set all axis colors to black
set(gca, 'XColor', 'k', 'YColor', 'k');

% Format x-axis - show only one label per month
xticks(monthPositions);
xticklabels(monthLabels);
xtickangle(45);

% Save figure
saveas(gcf, fullfile(saveDir, 'Temperature_Corrected_R1s_TimeSeries_YearByYear.fig'));

%% Figure 3: 3D Plot - Temperature Corrected R1s

figure('Position', [100, 100, 1200, 800]);

% Create 3D scatter plot with same color scheme as TimeSeries
scatter3(final_Temp_inst, xAxis, final_R1s_corrected*1000, 50, final_Temp_inst, 'filled');
hold on;

% Add 3D trend line
% Fit 3D trend: R1s = a*Temp + b*Date + c
X_matrix = [final_Temp_inst, xAxis, ones(length(final_Temp_inst), 1)];
coefficients = X_matrix \ (final_R1s_corrected*1000);

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
title('3D Plot: Battery Temperature vs Date vs R1s (Year-by-Year Corrected)', 'FontSize', fontSize+2, 'Color', 'k');

% Set y-axis (Date) ticks and labels - same as TimeSeries
set(gca, 'YTick', monthPositions);
set(gca, 'YTickLabel', monthLabels);

% Add colorbar for battery temperature (match TimeSeries)
ax3d = gca;
colormap(ax3d, flipud(autumn));
c = colorbar(ax3d, 'eastoutside');
c.Label.String = 'Battery Temperature (°C)';
c.Label.FontSize = fontSize;
c.Limits = [min(final_Temp_inst), max(final_Temp_inst)];
c.Ticks = round(min(final_Temp_inst)):1:round(max(final_Temp_inst));

% Add grid
grid on;

% Add legend
legend('Location', 'best', 'FontSize', fontSize-2);

% Set all axis colors to black
set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');

% Set view angle for better visualization
view(45, 30);

% Save 3D figure
saveas(gcf, fullfile(saveDir, 'R1s_3D_Plot_Temperature_Corrected_YearByYear.fig'));

%% Figure 4: Daily Average R1s 3D Plot

% Calculate daily averages from 1-second data
unique_dates = unique(final_dates);
daily_avg_dates = [];
daily_avg_R1s_uncorrected = [];
daily_avg_R1s_corrected = [];
daily_avg_Temp = [];
daily_avg_source = [];
daily_avg_year = [];

for i = 1:length(unique_dates)
    date = unique_dates(i);
    date_idx = final_dates == date;
    
    if sum(date_idx) > 0
        daily_avg_dates = [daily_avg_dates; date];
        daily_avg_R1s_uncorrected = [daily_avg_R1s_uncorrected; mean(final_R1s_uncorrected(date_idx))];
        daily_avg_R1s_corrected = [daily_avg_R1s_corrected; mean(final_R1s_corrected(date_idx))];
        daily_avg_Temp = [daily_avg_Temp; mean(final_Temp_inst(date_idx))];
        
        % Determine source (OldData, NewData, or Mixed)
        source_values = final_source(date_idx);
        if any(strcmp(source_values, 'OldData')) && any(strcmp(source_values, 'NewData'))
            daily_avg_source = [daily_avg_source; {'Mixed'}];
        elseif any(strcmp(source_values, 'OldData'))
            daily_avg_source = [daily_avg_source; {'OldData'}];
        else
            daily_avg_source = [daily_avg_source; {'NewData'}];
        end
        
        % Determine year
        year_values = final_year(date_idx);
        daily_avg_year = [daily_avg_year; year_values(1)];
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
title('3D Plot: Daily Average R1s (Year-by-Year Temperature Corrected)', 'FontSize', fontSize+2, 'Color', 'k');

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
c.Limits = [min(daily_avg_Temp), max(daily_avg_Temp)];
c.Ticks = round(min(daily_avg_Temp)):1:round(max(daily_avg_Temp));

% Add grid
grid on;

% Add legend
legend('Location', 'best', 'FontSize', fontSize-2);

% Set all axis colors to black
set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');

% Set view angle for better visualization
view(45, 30);

% Save 3D figure
saveas(gcf, fullfile(saveDir, 'R1s_3D_Plot_Daily_Average_YearByYear.fig'));

%% Year-by-Year Model Summary
fprintf('\n=== Year-by-Year Model Summary ===\n');
for i = 1:length(all_years)
    year = all_years{i};
    field_name = ['year_' year];
    if isfield(yearly_models, field_name)
        model = yearly_models.(field_name);
        outlier_stats = yearly_outlier_stats.(field_name);
        
        fprintf('Year %s:\n', year);
        fprintf('  Data points: %d (%.1f%% outliers removed)\n', ...
            outlier_stats.n_total - outlier_stats.n_outliers, outlier_stats.outlier_percentage);
        fprintf('  Model: R1s at 25°C = %.6f mΩ, Slope = %.6f\n', ...
            model.R1s_Tref_log*1000, model.p_log(1));
        fprintf('  R² = %.4f, p-value = %.2e\n', model.stats_log(1), model.stats_log(3));
    end
end

%% Summary Statistics
fprintf('\n=== Summary Statistics ===\n');
fprintf('Original R1s - Mean: %.4f mΩ, Std: %.4f mΩ\n', mean(final_R1s_uncorrected)*1000, std(final_R1s_uncorrected)*1000);
fprintf('Corrected R1s - Mean: %.4f mΩ, Std: %.4f mΩ\n', mean(final_R1s_corrected)*1000, std(final_R1s_corrected)*1000);

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
    yearR1s_corr = final_R1s_corrected(yearIdx);
    
    fprintf('Year %d: %d points, Original: %.4f±%.4f mΩ, Corrected: %.4f±%.4f mΩ\n', ...
        yearVal, sum(yearIdx), mean(yearR1s_orig)*1000, std(yearR1s_orig)*1000, ...
        mean(yearR1s_corr)*1000, std(yearR1s_corr)*1000);
end

%% Save results
fprintf('\nSaving results...\n');

% Save corrected data
arrhenius_results_ver04 = struct();
arrhenius_results_ver04.final_dates = final_dates;
arrhenius_results_ver04.final_R1s_uncorrected = final_R1s_uncorrected;
arrhenius_results_ver04.final_R1s_corrected = final_R1s_corrected;
arrhenius_results_ver04.final_Temp_inst = final_Temp_inst;
arrhenius_results_ver04.final_source = final_source;
arrhenius_results_ver04.final_year = final_year;
arrhenius_results_ver04.yearly_models = yearly_models;
arrhenius_results_ver04.yearly_outlier_stats = yearly_outlier_stats;
arrhenius_results_ver04.Tref = Tref;
arrhenius_results_ver04.all_years = all_years;

save(fullfile(saveDir, 'Arrhenius_Analysis_Results_YearByYear.mat'), 'arrhenius_results_ver04');

fprintf('\nYear-by-year Arrhenius analysis complete! Results saved to: %s\n', saveDir);
