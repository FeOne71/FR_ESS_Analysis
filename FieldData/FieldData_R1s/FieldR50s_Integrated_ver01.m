%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldR50s_Integrated_ver01.m
% ESS Rack01 R50s for Integrated Data Processing with 50-second Sampling
% Combines OldData and NewData processing with 30-second moving average
% and 50-second interval sampling for R50s calculation
% 
% Key Features:
% 1. Integrated processing of both OldData and NewData
% 2. 30-second moving average for V and I signals
% 3. 50-second interval sampling for dV and dI calculation
% 4. Daily R50s calculation using linear fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

%% Configuration
% Data directories
oldDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old';
newDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New';

% Save directory
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_R1s','R50s_Results_Integrated_ver01');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameters
% R50s calculation parameters
Cnom = 128;
Cnom_cell = 64;                   % Rack nominal Capacity (Ah)
idle_thr = Cnom_cell*0.05;         % Idle threshold [charge, discharge] (A)

% Moving average parameters
window_size = 30;                 % 30-second moving average window
sampling_interval = 50;          % 50-second sampling interval for R50s

% Topology
Ns = 17*14;    % 238s
Np = 2;        % 2p

% Visualization parameters
Fontsize = 12;
LineWidth = 2;

%% Data Processing Configuration
% OldData years
oldDataYears = {'2021', '2022', '2023'};
% NewData years  
newDataYears = {'2023', '2024', '2025'};
% Rack names
rackNames_all = {'Rack01'};

%% Process OldData
fprintf('=== Processing OldData for R50s ===\n');
for year_idx = 1:length(oldDataYears)
    year = oldDataYears{year_idx};
    fprintf('\nProcessing OldData year: %s\n', year);
    
    % Initialize R50s results structure for current year
    R50s_Results = struct();
    yearPath = fullfile(oldDataDir, year);
    
    fprintf('Year path: %s\n', yearPath);
    fprintf('Year path exists: %d\n', exist(yearPath, 'dir'));
    
    if ~exist(yearPath, 'dir')
        fprintf('Directory not found, skipping year %s\n', year);
        continue;
    end
    
    monthDirs = dir(fullfile(yearPath, '20*'));
    fprintf('Found %d month directories\n', length(monthDirs));

    for m = 1:length(monthDirs)
        if ~monthDirs(m).isdir, continue; end
        monthPath = fullfile(yearPath, monthDirs(m).name);
        fprintf('Processing month: %s\n', monthDirs(m).name);
        matFiles = dir(fullfile(monthPath, '*.mat'));
        fprintf('Found %d mat files in %s\n', length(matFiles), monthDirs(m).name);
        
        % Sort mat files by name (date order)
        [~, idx] = sort({matFiles.name});
        matFiles = matFiles(idx);

        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            fprintf('Loading file: %s\n', matFiles(f).name);
            
            try
                load(matFilePath);
                
                for rack_idx = 1:length(rackNames_all)
                    rackName = rackNames_all{rack_idx};
                    fprintf('  %s: ', rackName);

                    if ~isfield(Raw, rackName)
                        fprintf('No data found\n');
                        continue;
                    end

                    rackData = Raw.(rackName);
                    
                    % Extract signals
                    t = rackData.Time;
                    I = rackData.DCCurrent_A;
                    V = rackData.AverageCV_V;
                    soc = rackData.SOCPct;
                    T_batt = rackData.AverageMT_degC;
                    P = rackData.DCPower_kW;
                    
                    % Convert to cell units
                    I = I / Np;  % Convert rack current to cell current

                    % Convert time data to seconds-from-start
                    t = datetime(t);
                    t = seconds(t - t(1));
                    
                    % Apply moving average for noise reduction
                    fprintf('Applying %d-second moving average...\n', window_size);
                    I_smooth = smoothdata(I, 'movmean', window_size);
                    V_smooth = smoothdata(V, 'movmean', window_size);
                    
                    % Sample at 50-second intervals for R50s calculation
                    fprintf('Sampling at %d-second intervals...\n', sampling_interval);
                    sample_indices = 1:sampling_interval:length(I_smooth);
                    if sample_indices(end) ~= length(I_smooth)
                        sample_indices = [sample_indices, length(I_smooth)];
                    end
                    
                    I_sampled = I_smooth(sample_indices);
                    V_sampled = V_smooth(sample_indices);
                    T_sampled = T_batt(sample_indices);
                    
                    % Get SOC for sampled data
                    soc_sampled = soc(sample_indices);
                    soc_inst = soc_sampled(1:end-1); % SOC for dI/dV points
                    
                    % Calculate dV and dI from sampled signals
                    dI = diff(I_sampled);
                    dV = diff(V_sampled);
                    T_inst = T_sampled(1:end-1); % Temperature not smoothed
                    
                    % Debug: Check dV and dI ranges
                    fprintf('    dI range: %.6f to %.6f A\n', min(dI), max(dI));
                    fprintf('    dV range: %.6f to %.6f V\n', min(dV), max(dV));
                    
                    % Filter valid data points (NaN, dI=0, and dV=0 removal)
                    valid_idx = ~isnan(dI) & ~isnan(dV) & dI ~= 0 & dV ~= 0;
                    dI_valid = dI(valid_idx);
                    dV_valid = dV(valid_idx);
                    T_inst_valid = T_inst(valid_idx);
                    soc_inst_valid = soc_inst(valid_idx);
                    
                    % Calculate R50s using linear fitting: dV = a * dI (y-intercept = 0)
                    if length(dI_valid) >= 10 % Need sufficient data points for reliable fitting
                        % Linear fitting without y-intercept: dV = a * dI
                        p = polyfit(dI_valid, dV_valid, 1);
                        R50s_daily = p(1); % Fitted slope 'a'
                        
                        % Calculate R-squared for model validation
                        dV_predicted = R50s_daily * dI_valid;
                        ss_res = sum((dV_valid - dV_predicted).^2);
                        ss_tot = sum((dV_valid - mean(dV_valid)).^2);
                        r_squared = 1 - (ss_res / ss_tot);
                        
                        % Store results
                        dayKey = matFiles(f).name(1:end-4); % Remove .mat extension
                        R50s_Results.(dayKey).(rackName).R50s_daily = R50s_daily; % Daily R50s (slope)
                        R50s_Results.(dayKey).(rackName).dI_materials = dI_valid; % Store dI materials
                        R50s_Results.(dayKey).(rackName).dV_materials = dV_valid; % Store dV materials
                        R50s_Results.(dayKey).(rackName).T_inst = T_inst_valid;
                        R50s_Results.(dayKey).(rackName).SOC_inst = soc_inst_valid; % Store SOC
                        R50s_Results.(dayKey).(rackName).SOC_daily = mean(soc_inst_valid); % Daily average SOC
                        R50s_Results.(dayKey).(rackName).n_points = length(dI_valid);
                        R50s_Results.(dayKey).(rackName).r_squared = r_squared;
                        R50s_Results.(dayKey).(rackName).source = 'OldData';
                        R50s_Results.(dayKey).(rackName).window_size = window_size;
                        R50s_Results.(dayKey).(rackName).sampling_interval = sampling_interval;
                        
                        fprintf('R50s_daily = %.6f Ω (%.3f mΩ), R² = %.4f, n_points = %d, T range: %.1f-%.1f°C\n', ...
                            R50s_daily, R50s_daily*1000, r_squared, length(dI_valid), min(T_inst_valid), max(T_inst_valid));
                    else
                        fprintf('Insufficient data points (%d) for reliable R50s calculation\n', length(dI_valid));
                    end
                end
                
            catch ME
                fprintf('ERROR processing %s: %s\n', matFiles(f).name, ME.message);
                continue;
            end
        end % End mat file loop
    end % End month loop
    
    % Save results for current year
    if ~isempty(fieldnames(R50s_Results))
        saveFileName = sprintf('R50s_Results_OldData_%s.mat', year);
        save(fullfile(saveDir, saveFileName), 'R50s_Results');
        fprintf('Saved OldData results to: %s\n', saveFileName);
        
        % Create summary statistics
        createYearlySummary(R50s_Results, year, 'OldData', rackNames_all);
    else
        fprintf('No valid OldData found for %s\n', year);
    end
end % End OldData year loop

%% Process NewData
fprintf('\n=== Processing NewData for R50s ===\n');
for year_idx = 1:length(newDataYears)
    year = newDataYears{year_idx};
    fprintf('\nProcessing NewData year: %s\n', year);
    
    % Initialize R50s results structure for current year
    R50s_Results = struct();
    yearPath = fullfile(newDataDir, year);
    
    fprintf('Year path: %s\n', yearPath);
    fprintf('Year path exists: %d\n', exist(yearPath, 'dir'));
    
    if ~exist(yearPath, 'dir')
        fprintf('Directory not found, skipping year %s\n', year);
        continue;
    end
    
    monthDirs = dir(fullfile(yearPath, '20*'));
    fprintf('Found %d month directories\n', length(monthDirs));

    for m = 1:length(monthDirs)
        if ~monthDirs(m).isdir, continue; end
        monthPath = fullfile(yearPath, monthDirs(m).name);
        fprintf('Processing month: %s\n', monthDirs(m).name);
        matFiles = dir(fullfile(monthPath, '*.mat'));
        fprintf('Found %d mat files in %s\n', length(matFiles), monthDirs(m).name);
        
        % Sort mat files by name (date order)
        [~, idx] = sort({matFiles.name});
        matFiles = matFiles(idx);

        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            fprintf('Loading file: %s\n', matFiles(f).name);
            
            try
                S = load(matFilePath);
                
                % Determine the data struct variable
                if isfield(S,'Raw')
                    Raw = S.Raw;
                else
                    % Find first struct containing expected fields
                    vars = fieldnames(S);
                    Raw = [];
                    for vi = 1:numel(vars)
                        val = S.(vars{vi});
                        if isstruct(val)
                            fieldNames = fieldnames(val);
                            if any(strcmpi(fieldNames,'CVavg')) && any(strcmpi(fieldNames,'DCCurrent'))
                                Raw = val; break;
                            end
                        end
                    end
                    if isempty(Raw)
                        fprintf('No suitable struct found in %s\n', matFiles(f).name);
                        continue;
                    end
                end
                
                for rack_idx = 1:length(rackNames_all)
                    rackName = rackNames_all{rack_idx};
                    fprintf('  %s: ', rackName);

                    % New data may be stored directly as a struct (no Rack01 field)
                    if isfield(Raw, rackName)
                        D = Raw.(rackName);
                    else
                        D = Raw; % direct struct (single rack)
                    end

                    % Time (debug and handle different data types)
                    if isfield(D,'Time')
                        t = D.Time;
                    elseif isfield(D,'Date_Time')
                        t = D.Date_Time;
                    else
                        fprintf('No Time/Date_Time found\n');
                        continue;
                    end
                    
                    % Convert time data to seconds-from-start
                    if isdatetime(t)
                        % Already datetime
                        t = seconds(t - t(1));
                    elseif isduration(t)
                        % Duration data
                        t = seconds(t);
                    elseif isnumeric(t)
                        % Numeric data (seconds or other)
                        t = t - t(1);
                    else
                        % Try to convert to datetime
                        try
                            t = datetime(t);
                            t = seconds(t - t(1));
                        catch
                            fprintf('Cannot convert time data to datetime\n');
                            continue;
                        end
                    end

                    % Signals
                    I = D.DCCurrent(:);
                    V = D.CVavg(:);
                    soc = D.SOC_BMS(:);
                    P = D.DCPower(:);
                    T_batt = D.MTavg(:);
                    
                    % Convert to cell units
                    I = I / Np;  % Convert to cell current
                    
                    % Apply moving average for noise reduction
                    fprintf('Applying %d-second moving average...\n', window_size);
                    I_smooth = smoothdata(I, 'movmean', window_size);
                    V_smooth = smoothdata(V, 'movmean', window_size);
                    
                    % Sample at 50-second intervals for R50s calculation
                    fprintf('Sampling at %d-second intervals...\n', sampling_interval);
                    sample_indices = 1:sampling_interval:length(I_smooth);
                    if sample_indices(end) ~= length(I_smooth)
                        sample_indices = [sample_indices, length(I_smooth)];
                    end
                    
                    I_sampled = I_smooth(sample_indices);
                    V_sampled = V_smooth(sample_indices);
                    T_sampled = T_batt(sample_indices);
                    
                    % Calculate dV and dI from sampled signals
                    dI = diff(I_sampled);
                    dV = diff(V_sampled);
                    T_inst = T_sampled(1:end-1); % Temperature not smoothed
                    
                    % Debug: Check dV and dI ranges
                    fprintf('    dI range: %.6f to %.6f A\n', min(dI), max(dI));
                    fprintf('    dV range: %.6f to %.6f V\n', min(dV), max(dV));
                    
                    % Filter valid data points (NaN, dI=0, and dV=0 removal)
                    valid_idx = ~isnan(dI) & ~isnan(dV) & dI ~= 0 & dV ~= 0;
                    dI_valid = dI(valid_idx);
                    dV_valid = dV(valid_idx);
                    T_inst_valid = T_inst(valid_idx);
                    
                    % Calculate R50s using linear fitting: dV = a * dI (y-intercept = 0)
                    if length(dI_valid) >= 10 % Need sufficient data points for reliable fitting
                        % Linear fitting without y-intercept: dV = a * dI
                        p = polyfit(dI_valid, dV_valid, 1);
                        R50s_daily = p(1); % Fitted slope 'a'
                        
                        % Calculate R-squared for model validation
                        dV_predicted = R50s_daily * dI_valid;
                        ss_res = sum((dV_valid - dV_predicted).^2);
                        ss_tot = sum((dV_valid - mean(dV_valid)).^2);
                        r_squared = 1 - (ss_res / ss_tot);
                        
                        % Store results
                        dayKey = matFiles(f).name(1:end-4); % Remove .mat extension
                        R50s_Results.(dayKey).(rackName).R50s_daily = R50s_daily; % Daily R50s (slope)
                        R50s_Results.(dayKey).(rackName).dI_materials = dI_valid; % Store dI materials
                        R50s_Results.(dayKey).(rackName).dV_materials = dV_valid; % Store dV materials
                        R50s_Results.(dayKey).(rackName).T_inst = T_inst_valid;
                        R50s_Results.(dayKey).(rackName).SOC_inst = soc_inst_valid; % Store SOC
                        R50s_Results.(dayKey).(rackName).SOC_daily = mean(soc_inst_valid); % Daily average SOC
                        R50s_Results.(dayKey).(rackName).n_points = length(dI_valid);
                        R50s_Results.(dayKey).(rackName).r_squared = r_squared;
                        R50s_Results.(dayKey).(rackName).source = 'NewData';
                        R50s_Results.(dayKey).(rackName).window_size = window_size;
                        R50s_Results.(dayKey).(rackName).sampling_interval = sampling_interval;
                        
                        fprintf('R50s_daily = %.6f Ω (%.3f mΩ), R² = %.4f, n_points = %d, T range: %.1f-%.1f°C\n', ...
                            R50s_daily, R50s_daily*1000, r_squared, length(dI_valid), min(T_inst_valid), max(T_inst_valid));
                    else
                        fprintf('Insufficient data points (%d) for reliable R50s calculation\n', length(dI_valid));
                    end
                end
                
            catch ME
                fprintf('ERROR processing %s: %s\n', matFiles(f).name, ME.message);
                continue;
            end
        end % End mat file loop
    end % End month loop
    
    % Save results for current year
    if ~isempty(fieldnames(R50s_Results))
        saveFileName = sprintf('R50s_Results_NewData_%s.mat', year);
        save(fullfile(saveDir, saveFileName), 'R50s_Results');
        fprintf('Saved NewData results to: %s\n', saveFileName);
        
        % Create summary statistics
        createYearlySummary(R50s_Results, year, 'NewData', rackNames_all);
    else
        fprintf('No valid NewData found for %s\n', year);
    end
end % End NewData year loop

%% Create ln(R50s) vs 1/T Visualization by Year
fprintf('\n=== Creating ln(R50s) vs 1/T Visualization by Year ===\n');
createR50sVsTemperaturePlot(saveDir);

%% Year-by-Year Temperature Correction and Final Visualization
fprintf('\n=== Year-by-Year Temperature Correction and Final Visualization ===\n');
createYearlyCorrectedVisualization(saveDir);

%% SOC-Controlled Visualization
fprintf('\n=== Creating SOC-Controlled Visualizations ===\n');
createSOCControlledVisualization(saveDir);

fprintf('\nIntegrated R50s processing completed!\n');
fprintf('Results saved to: %s\n', saveDir);

%% Helper Functions

function createYearlySummary(R50s_Results, year, dataType, rackNames_all)
    % Create summary statistics for a year
    allR50s_daily = [];
    allR_squared = [];
    allN = [];
    
    dayNames = fieldnames(R50s_Results);
    for d = 1:length(dayNames)
        if isfield(R50s_Results.(dayNames{d}), rackNames_all{1})
            allR50s_daily = [allR50s_daily; R50s_Results.(dayNames{d}).(rackNames_all{1}).R50s_daily];
            allR_squared = [allR_squared; R50s_Results.(dayNames{d}).(rackNames_all{1}).r_squared];
            allN = [allN; R50s_Results.(dayNames{d}).(rackNames_all{1}).n_points];
        end
    end
    
    if ~isempty(allR50s_daily)
        fprintf('\n=== %s %s Summary ===\n', year, dataType);
        fprintf('Total days processed: %d\n', length(dayNames));
        fprintf('Mean R50s_daily: %.4f mΩ\n', mean(allR50s_daily)*1000);
        fprintf('Std R50s_daily: %.4f mΩ\n', std(allR50s_daily)*1000);
        fprintf('Mean R²: %.4f\n', mean(allR_squared));
        fprintf('Mean materials per day: %.1f\n', mean(allN));
    end
end

function createR50sVsTemperaturePlot(saveDir)
    % Create R50s vs 1/T visualization by year
    
    % Find all saved result files
    resultFiles = dir(fullfile(saveDir, 'R50s_Results_*.mat'));
    
    if isempty(resultFiles)
        fprintf('No result files found in %s\n', saveDir);
        return;
    end
    
    % Collect data by year
    yearly_data = struct();
    
    for i = 1:length(resultFiles)
        fileName = resultFiles(i).name;
        filePath = fullfile(saveDir, fileName);
        
        try
            load(filePath, 'R50s_Results');
            
            % Extract year from filename
            if contains(fileName, 'OldData_')
                year = fileName(strfind(fileName, 'OldData_') + 8:end-4);
                dataType = 'OldData';
            elseif contains(fileName, 'NewData_')
                year = fileName(strfind(fileName, 'NewData_') + 8:end-4);
                dataType = 'NewData';
            else
                continue;
            end
            
            % Extract R50s and temperature data
            dayNames = fieldnames(R50s_Results);
            year_R50s = [];
            year_T = [];
            
            for d = 1:length(dayNames)
                if isfield(R50s_Results.(dayNames{d}), 'Rack01')
                    rackData = R50s_Results.(dayNames{d}).Rack01;
                    if isfield(rackData, 'R50s_daily') && isfield(rackData, 'T_inst')
                        year_R50s = [year_R50s; rackData.R50s_daily];
                        year_T = [year_T; mean(rackData.T_inst)]; % Use mean temperature for the day
                    end
                end
            end
            
            % Store yearly data
            field_name = ['year_' year];
            if isfield(yearly_data, field_name)
                % Merge with existing data
                yearly_data.(field_name).R50s = [yearly_data.(field_name).R50s; year_R50s];
                yearly_data.(field_name).T = [yearly_data.(field_name).T; year_T];
                yearly_data.(field_name).source = 'Mixed';
            else
                % New year
                yearly_data.(field_name).R50s = year_R50s;
                yearly_data.(field_name).T = year_T;
                yearly_data.(field_name).source = dataType;
            end
            
        catch ME
            fprintf('Error loading %s: %s\n', fileName, ME.message);
        end
    end
    
    % Create visualization
    if isempty(fieldnames(yearly_data))
        fprintf('No valid yearly data found for visualization\n');
        return;
    end
    
    % Get all years and sort them
    yearFields = fieldnames(yearly_data);
    years = {};
    for i = 1:length(yearFields)
        years{end+1} = yearFields{i}(6:end); % Remove 'year_' prefix
    end
    years = sort(years);
    
    % Apply year-by-year outlier removal
    fprintf('\n=== Year-by-Year Outlier Removal ===\n');
    yearly_outlier_stats = struct();
    z_threshold = 2.5; % 2.5-sigma rule (more strict than 3-sigma)
    
    for i = 1:length(years)
        year = years{i};
        field_name = ['year_' year];
        
        R50s_data = yearly_data.(field_name).R50s;
        T_data = yearly_data.(field_name).T;
        
        if length(R50s_data) > 1
            % Calculate Z-Score for this year's data
            mean_val = mean(R50s_data);
            std_val = std(R50s_data);
            z_scores = abs((R50s_data - mean_val) / std_val);
            
            % Identify outliers
            outlier_idx = z_scores > z_threshold;
            
            % Count outliers
            n_outliers = sum(outlier_idx);
            n_total = length(R50s_data);
            outlier_percentage = n_outliers / n_total * 100;
            
            % Calculate bounds for display
            lower_bound = mean_val - z_threshold * std_val;
            upper_bound = mean_val + z_threshold * std_val;
            
            % Debug: Show data range and outlier details
            fprintf('Year %s: Mean: %.6f mΩ, Std: %.6f mΩ, Z-threshold: %.1f\n', ...
                year, mean_val*1000, std_val*1000, z_threshold);
            fprintf('  Data range: %.4f to %.4f mΩ\n', min(R50s_data)*1000, max(R50s_data)*1000);
            fprintf('  Lower bound: %.6f mΩ, Upper bound: %.6f mΩ\n', lower_bound*1000, upper_bound*1000);
            
            % Show which values are outliers
            if n_outliers > 0
                outlier_values = R50s_data(outlier_idx);
                outlier_z_scores = z_scores(outlier_idx);
                fprintf('  Outlier values: ');
                for j = 1:min(10, length(outlier_values))
                    fprintf('%.4f mΩ (Z=%.2f) ', outlier_values(j)*1000, outlier_z_scores(j));
                end
                if length(outlier_values) > 10
                    fprintf('... (%d more)', length(outlier_values) - 10);
                end
                fprintf('\n');
            end
            
            fprintf('  Outliers removed: %d out of %d (%.1f%%)\n', n_outliers, n_total, outlier_percentage);
            
            % Remove outliers
            R50s_filtered = R50s_data(~outlier_idx);
            T_filtered = T_data(~outlier_idx);
            
            fprintf('  Filtered data: %d points, R50s range: %.4f to %.4f mΩ\n', ...
                length(R50s_filtered), min(R50s_filtered)*1000, max(R50s_filtered)*1000);
            
            % Store filtered data
            yearly_data.(field_name).R50s = R50s_filtered;
            yearly_data.(field_name).T = T_filtered;
            
            % Store outlier statistics
            yearly_outlier_stats.(field_name) = struct();
            yearly_outlier_stats.(field_name).n_outliers = n_outliers;
            yearly_outlier_stats.(field_name).n_total = n_total;
            yearly_outlier_stats.(field_name).outlier_percentage = outlier_percentage;
            yearly_outlier_stats.(field_name).mean_val = mean_val;
            yearly_outlier_stats.(field_name).std_val = std_val;
            yearly_outlier_stats.(field_name).bounds = [lower_bound, upper_bound];
        end
    end
    
    % Create figure
    figure('Position', [100, 100, 1400, 1000]);
    
    % Calculate subplot layout
    nYears = length(years);
    nCols = min(3, nYears);
    nRows = ceil(nYears / nCols);
    
    % Calculate global axis limits for consistent scaling (using filtered data)
    all_inv_T = [];
    all_ln_R50s = [];
    for i = 1:nYears
        year = years{i};
        field_name = ['year_' year];
        R50s_data = yearly_data.(field_name).R50s;
        T_data = yearly_data.(field_name).T;
        
        if length(R50s_data) > 1
            T_K = T_data + 273.15;
            inv_T = 1 ./ T_K;
            ln_R50s = log(R50s_data * 1000);
            
            all_inv_T = [all_inv_T; inv_T];
            all_ln_R50s = [all_ln_R50s; ln_R50s];
        end
    end
    
    % Set global axis limits with small margins
    x_margin = (max(all_inv_T) - min(all_inv_T)) * 0.05;
    y_margin = (max(all_ln_R50s) - min(all_ln_R50s)) * 0.05;
    x_lim = [min(all_inv_T) - x_margin, max(all_inv_T) + x_margin];
    y_lim = [min(all_ln_R50s) - y_margin, max(all_ln_R50s) + y_margin];
    
    for i = 1:nYears
        year = years{i};
        field_name = ['year_' year];
        
        subplot(nRows, nCols, i);
        
        R50s_data = yearly_data.(field_name).R50s;
        T_data = yearly_data.(field_name).T;
        source = yearly_data.(field_name).source;
        
        % Convert temperature to 1/T (K^-1)
        T_K = T_data + 273.15; % Convert to Kelvin
        inv_T = 1 ./ T_K;
        
        % Convert R50s to ln(R50s)
        ln_R50s = log(R50s_data * 1000); % Convert to mΩ and take natural log
        
        % Create scatter plot
        scatter(inv_T, ln_R50s, 50, T_data, 'filled');
        
        % Add linear trend line
        if length(R50s_data) > 1
            mdl = fitlm(inv_T, ln_R50s);
            slope = mdl.Coefficients.Estimate(2);
            intercept = mdl.Coefficients.Estimate(1);
            r_squared = mdl.Rsquared.Ordinary;
            
            x_trend = linspace(min(inv_T), max(inv_T), 100);
            y_trend = slope * x_trend + intercept;
            hold on;
            plot(x_trend, y_trend, 'r-', 'LineWidth', 2);
            
            % Add trend line info with temperature distribution
            temp_min = min(T_data);
            temp_max = max(T_data);
            temp_mean = mean(T_data);
            temp_std = std(T_data);
            
            text(0.05, 0.95, sprintf('Slope: %.2f K\nIntercept: %.2f\nR²: %.3f\nn = %d days\nTemp: %.1f±%.1f°C\nRange: %.1f-%.1f°C', ...
                slope, intercept, r_squared, length(R50s_data), temp_mean, temp_std, temp_min, temp_max), ...
                'Units', 'normalized', 'VerticalAlignment', 'top', ...
                'BackgroundColor', 'white', 'EdgeColor', 'black');
        end
        
        % Formatting
        xlabel('1/T (K^{-1})', 'FontSize', 10);
        ylabel('ln(R50s) (ln(mΩ))', 'FontSize', 10);
        title(sprintf('Year %s', year), 'FontSize', 12);
        colorbar;
        colormap(flipud(autumn));
        caxis([min(T_data), max(T_data)]);
        grid on;
        
        % Set consistent axis limits for all subplots
        xlim(x_lim);
        ylim(y_lim);
    end
    
    % Add overall title
    sgtitle('ln(R50s) vs 1/T by Year (Arrhenius Relationship)', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save figure
    saveas(gcf, fullfile(saveDir, 'R50s_vs_1T_by_Year.fig'));    
    fprintf('ln(R50s) vs 1/T visualization saved to: %s\n', saveDir);
    
    % Print yearly slope, intercept, R², and temperature distribution summary
    fprintf('\n=== Yearly Arrhenius Analysis Summary ===\n');
    fprintf('Year\t\tSlope (K)\tIntercept\tR²\t\tData Points\tTemp Range (°C)\tTemp Mean±Std (°C)\n');
    fprintf('----\t\t---------\t---------\t---\t\t-----------\t---------------\t------------------\n');
    
    for i = 1:nYears
        year = years{i};
        field_name = ['year_' year];
        
        R50s_data = yearly_data.(field_name).R50s;
        T_data = yearly_data.(field_name).T;
        
        if length(R50s_data) > 1
            % Convert temperature to 1/T (K^-1)
            T_K = T_data + 273.15;
            inv_T = 1 ./ T_K;
            
            % Convert R50s to ln(R50s)
            ln_R50s = log(R50s_data * 1000);
            
            % Calculate slope, intercept, and R² using MATLAB functions
            mdl = fitlm(inv_T, ln_R50s);
            slope = mdl.Coefficients.Estimate(2);
            intercept = mdl.Coefficients.Estimate(1);
            r_squared = mdl.Rsquared.Ordinary;
            
            % Calculate temperature statistics
            temp_min = min(T_data);
            temp_max = max(T_data);
            temp_mean = mean(T_data);
            temp_std = std(T_data);
            
            fprintf('%s\t\t%.2f\t\t%.2f\t\t%.4f\t\t%d\t\t%.1f-%.1f\t\t%.1f±%.1f\n', ...
                year, slope, intercept, r_squared, length(R50s_data), temp_min, temp_max, temp_mean, temp_std);
        else
            fprintf('%s\t\tN/A\t\tN/A\t\tN/A\t\t%d\t\tN/A\t\tN/A\n', year, length(R50s_data));
        end
    end
    fprintf('\n');
    
    % Print outlier removal summary
    fprintf('=== Outlier Removal Summary ===\n');
    for i = 1:nYears
        year = years{i};
        field_name = ['year_' year];
        if isfield(yearly_outlier_stats, field_name)
            stats = yearly_outlier_stats.(field_name);
            fprintf('Year %s: Removed %d outliers (%.1f%%) from %d total points\n', ...
                year, stats.n_outliers, stats.outlier_percentage, stats.n_total);
        end
    end
    fprintf('\n');
end

function createYearlyCorrectedVisualization(saveDir)
    % Create year-by-year temperature corrected visualization
    
    % Find all saved result files
    resultFiles = dir(fullfile(saveDir, 'R50s_Results_*.mat'));
    
    if isempty(resultFiles)
        fprintf('No result files found in %s\n', saveDir);
        return;
    end
    
    % Collect data by year
    yearly_data = struct();
    
    for i = 1:length(resultFiles)
        fileName = resultFiles(i).name;
        filePath = fullfile(saveDir, fileName);
        
        try
            load(filePath, 'R50s_Results');
            
            % Extract year from filename
            if contains(fileName, 'OldData_')
                year = fileName(strfind(fileName, 'OldData_') + 8:end-4);
                dataType = 'OldData';
            elseif contains(fileName, 'NewData_')
                year = fileName(strfind(fileName, 'NewData_') + 8:end-4);
                dataType = 'NewData';
            else
                continue;
            end
            
            % Extract R50s, temperature, and SOC data
            dayNames = fieldnames(R50s_Results);
            year_R50s = [];
            year_T = [];
            year_SOC = [];
            year_dates = [];
            
            for d = 1:length(dayNames)
                if isfield(R50s_Results.(dayNames{d}), 'Rack01')
                    rackData = R50s_Results.(dayNames{d}).Rack01;
                    if isfield(rackData, 'R50s_daily') && isfield(rackData, 'T_inst')
                        year_R50s = [year_R50s; rackData.R50s_daily];
                        year_T = [year_T; mean(rackData.T_inst)]; % Use mean temperature for the day
                        if isfield(rackData, 'SOC_daily')
                            year_SOC = [year_SOC; rackData.SOC_daily]; % Use daily average SOC
                        else
                            year_SOC = [year_SOC; NaN]; % If SOC not available
                        end
                        
                        % Extract date from dayKey
                        dayKey = dayNames{d};
                        if length(dayKey) >= 12 && strcmp(dayKey(1:4), 'Raw_')
                            dateStr = dayKey(5:12); % '20210601'
                            year_num = str2double(dateStr(1:4));
                            month = str2double(dateStr(5:6));
                            day = str2double(dateStr(7:8));
                            date = datetime(year_num, month, day);
                            year_dates = [year_dates; date];
                        elseif length(dayKey) >= 8
                            % Try to parse as YYYYMMDD format
                            try
                                year_num = str2double(dayKey(1:4));
                                month = str2double(dayKey(5:6));
                                day = str2double(dayKey(7:8));
                                if ~isnan(year_num) && ~isnan(month) && ~isnan(day)
                                    date = datetime(year_num, month, day);
                                    year_dates = [year_dates; date];
                                end
                            catch
                                % If parsing fails, use a default date
                                year_dates = [year_dates; datetime(str2double(year), 1, 1)];
                            end
                        else
                            % If no date can be extracted, use a default date
                            year_dates = [year_dates; datetime(str2double(year), 1, 1)];
                        end
                    end
                end
            end
            
            % Store yearly data
            field_name = ['year_' year];
            if isfield(yearly_data, field_name)
                % Merge with existing data
                yearly_data.(field_name).R50s = [yearly_data.(field_name).R50s; year_R50s];
                yearly_data.(field_name).T = [yearly_data.(field_name).T; year_T];
                yearly_data.(field_name).SOC = [yearly_data.(field_name).SOC; year_SOC];
                yearly_data.(field_name).dates = [yearly_data.(field_name).dates; year_dates];
                yearly_data.(field_name).source = 'Mixed';
            else
                % New year
                yearly_data.(field_name).R50s = year_R50s;
                yearly_data.(field_name).T = year_T;
                yearly_data.(field_name).SOC = year_SOC;
                yearly_data.(field_name).dates = year_dates;
                yearly_data.(field_name).source = dataType;
            end
            
        catch ME
            fprintf('Error loading %s: %s\n', fileName, ME.message);
        end
    end
    
    % Create yearly Arrhenius models
    if isempty(fieldnames(yearly_data))
        fprintf('No valid yearly data found for visualization\n');
        return;
    end
    
    % Get all years and sort them
    yearFields = fieldnames(yearly_data);
    years = {};
    for i = 1:length(yearFields)
        years{end+1} = yearFields{i}(6:end); % Remove 'year_' prefix
    end
    years = sort(years);
    
    % Apply year-by-year outlier removal
    fprintf('\n=== Year-by-Year Outlier Removal ===\n');
    yearly_outlier_stats = struct();
    z_threshold = 2.5; % 2.5-sigma rule (more strict than 3-sigma)
    
    for i = 1:length(years)
        year = years{i};
        field_name = ['year_' year];
        
        R50s_data = yearly_data.(field_name).R50s;
        T_data = yearly_data.(field_name).T;
        dates = yearly_data.(field_name).dates;
        
        if length(R50s_data) > 1
            % Calculate Z-Score for this year's data
            mean_val = mean(R50s_data);
            std_val = std(R50s_data);
            z_scores = abs((R50s_data - mean_val) / std_val);
            
            % Identify outliers
            outlier_idx = z_scores > z_threshold;
            
            % Count outliers
            n_outliers = sum(outlier_idx);
            n_total = length(R50s_data);
            outlier_percentage = n_outliers / n_total * 100;
            
            % Calculate bounds for display
            lower_bound = mean_val - z_threshold * std_val;
            upper_bound = mean_val + z_threshold * std_val;
            
            % Debug: Show data range and outlier details
            fprintf('Year %s: Mean: %.6f mΩ, Std: %.6f mΩ, Z-threshold: %.1f\n', ...
                year, mean_val*1000, std_val*1000, z_threshold);
            fprintf('  Data range: %.4f to %.4f mΩ\n', min(R50s_data)*1000, max(R50s_data)*1000);
            fprintf('  Lower bound: %.6f mΩ, Upper bound: %.6f mΩ\n', lower_bound*1000, upper_bound*1000);
            
            % Show which values are outliers
            if n_outliers > 0
                outlier_values = R50s_data(outlier_idx);
                outlier_z_scores = z_scores(outlier_idx);
                outlier_dates = dates(outlier_idx);
                fprintf('  Outlier values: ');
                for j = 1:min(10, length(outlier_values))
                    fprintf('%.4f mΩ (Z=%.2f, %s) ', outlier_values(j)*1000, outlier_z_scores(j), datestr(outlier_dates(j), 'yyyy-mm-dd'));
                end
                if length(outlier_values) > 10
                    fprintf('... (%d more)', length(outlier_values) - 10);
                end
                fprintf('\n');
            end
            
            fprintf('  Outliers removed: %d out of %d (%.1f%%)\n', n_outliers, n_total, outlier_percentage);
            
            % Remove outliers
            R50s_filtered = R50s_data(~outlier_idx);
            T_filtered = T_data(~outlier_idx);
            dates_filtered = dates(~outlier_idx);
            
            fprintf('  Filtered data: %d points, R50s range: %.4f to %.4f mΩ\n', ...
                length(R50s_filtered), min(R50s_filtered)*1000, max(R50s_filtered)*1000);
            
            % Store filtered data
            yearly_data.(field_name).R50s = R50s_filtered;
            yearly_data.(field_name).T = T_filtered;
            yearly_data.(field_name).dates = dates_filtered;
            
            % Store outlier statistics
            yearly_outlier_stats.(field_name) = struct();
            yearly_outlier_stats.(field_name).n_outliers = n_outliers;
            yearly_outlier_stats.(field_name).n_total = n_total;
            yearly_outlier_stats.(field_name).outlier_percentage = outlier_percentage;
            yearly_outlier_stats.(field_name).mean_val = mean_val;
            yearly_outlier_stats.(field_name).std_val = std_val;
            yearly_outlier_stats.(field_name).bounds = [lower_bound, upper_bound];
        end
    end
    
    % Create yearly models
    yearly_models = struct();
    Tref = 25 + 273.15; % Reference temperature: 25°C in Kelvin
    
    fprintf('\nCreating yearly Arrhenius models (after outlier removal)...\n');
    for i = 1:length(years)
        year = years{i};
        field_name = ['year_' year];
        
        R50s_data = yearly_data.(field_name).R50s;
        T_data = yearly_data.(field_name).T;
        
        if length(R50s_data) > 1
            % Convert temperature to 1/T (K^-1)
            T_K = T_data + 273.15;
            inv_T = 1 ./ T_K;
            
            % Convert R50s to ln(R50s)
            ln_R50s = log(R50s_data * 1000);
            
            % Calculate slope and intercept using MATLAB functions
            mdl = fitlm(inv_T, ln_R50s);
            slope = mdl.Coefficients.Estimate(2);
            intercept = mdl.Coefficients.Estimate(1);
            
            % Store model
            yearly_models.(field_name).slope = slope;
            yearly_models.(field_name).intercept = intercept;
            
            fprintf('Year %s: Slope = %.2f K, Intercept = %.2f\n', year, slope, intercept);
        end
    end
    
    % Collect all data for final visualization
    all_dates = [];
    all_R50s_uncorrected = [];
    all_R50s_corrected = [];
    all_T = [];
    all_SOC = [];
    all_years = [];
    
    for i = 1:length(years)
        year = years{i};
        field_name = ['year_' year];
        
        if isfield(yearly_models, field_name)
            R50s_data = yearly_data.(field_name).R50s;
            T_data = yearly_data.(field_name).T;
            SOC_data = yearly_data.(field_name).SOC;
            dates = yearly_data.(field_name).dates;
            
            % Apply year-specific temperature correction
            model = yearly_models.(field_name);
            slope = model.slope;
            intercept = model.intercept;
            
            R50s_corrected = zeros(size(R50s_data));
            for j = 1:length(R50s_data)
                % Current temperature in Kelvin
                T_i_kelvin = T_data(j) + 273.15;
                
                % Correction formula
                ln_R_ref25 = slope * (1/298.15) + intercept;
                ln_R_model_i = slope * (1/T_i_kelvin) + intercept;
                Factor_i = exp(ln_R_ref25) / exp(ln_R_model_i);
                
                % Apply correction
                R50s_corrected(j) = R50s_data(j) * Factor_i;
            end
            
            % Store corrected data
            all_dates = [all_dates; dates];
            all_R50s_uncorrected = [all_R50s_uncorrected; R50s_data];
            all_R50s_corrected = [all_R50s_corrected; R50s_corrected];
            all_T = [all_T; T_data];
            all_SOC = [all_SOC; SOC_data];
            all_years = [all_years; repmat({year}, length(R50s_data), 1)];
        end
    end
    
    % Sort by date
    [all_dates, sortIdx] = sort(all_dates);
    all_R50s_uncorrected = all_R50s_uncorrected(sortIdx);
    all_R50s_corrected = all_R50s_corrected(sortIdx);
    all_T = all_T(sortIdx);
    all_SOC = all_SOC(sortIdx);
    all_years = all_years(sortIdx);
    
    % Save temperature corrected data for future use (after sorting)
    corrected_data = struct();
    corrected_data.dates = all_dates;
    corrected_data.R50s_uncorrected = all_R50s_uncorrected;
    corrected_data.R50s_corrected = all_R50s_corrected;
    corrected_data.T = all_T;
    corrected_data.SOC = all_SOC;
    corrected_data.years = all_years;
    corrected_data.yearly_models = yearly_models;
    corrected_data.yearly_outlier_stats = yearly_outlier_stats;
    
    save(fullfile(saveDir, 'Temperature_Corrected_R50s_Data.mat'), 'corrected_data', '-v7.3');
    fprintf('Temperature corrected data saved to: Temperature_Corrected_R50s_Data.mat\n');
    
    % Create sequential x-axis
    xAxis = (1:length(all_dates))';
    
    % Create final visualization
    figure('Position', [100, 100, 1400, 800]);
    
    % Plot all data as scatter points with battery temperature color
    scatter(xAxis, all_R50s_corrected*1000, 50, all_T, 'filled', 'HandleVisibility', 'off');
    hold on;
    
    % Add year division lines
    unique_years = unique(all_years);
    for i = 1:length(unique_years)
        year = unique_years{i};
        year_idx = strcmp(all_years, year);
        if sum(year_idx) > 0
            first_idx = find(year_idx, 1);
            xline(first_idx, '--', 'Color', 'k', 'LineWidth', 1, 'Alpha', 0.7, 'HandleVisibility', 'off');
        end
    end
    
    % Add trend line (green thick dashed line)
    if length(all_R50s_corrected) > 1
        mdl = fitlm(xAxis, all_R50s_corrected*1000, 'Intercept', true);
        trendLine = mdl.Fitted;
        plot(xAxis, trendLine, 'g--', 'LineWidth', 3, 'DisplayName', sprintf('Trend (R²=%.4f)', mdl.Rsquared.Ordinary));
    end
    
    % Colorbar for battery temperature
    ax = gca;
    colormap(ax, flipud(autumn));
    c = colorbar(ax, 'eastoutside');
    c.Label.String = 'Battery Temperature (°C)';
    c.Label.FontSize = 12;
    c.Label.FontWeight = 'bold';
    c.Limits = [min(all_T), max(all_T)];
    c.Ticks = round(min(all_T)):1:round(max(all_T));
    
    % Formatting with bold fonts and thick axes
    xlabel('Date (YYYY:MM:DD)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('R50s (mΩ) - Temperature Corrected to 25°C', 'FontSize', 12, 'FontWeight', 'bold');
    title('ESS Rack01 Daily R50s Time Series (Year-by-Year Temperature Corrected)', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % Set thick axes
    ax.LineWidth = 2;
    ax.FontWeight = 'bold';
    
    % Show only trend legend
    legend('Location', 'best', 'FontSize', 12);
    
    % Format x-axis - show only one label per month
    dateYears = zeros(size(all_dates));
    dateMonths = zeros(size(all_dates));
    for i = 1:length(all_dates)
        dateYears(i) = all_dates(i).Year;
        dateMonths(i) = all_dates(i).Month;
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
    saveas(gcf, fullfile(saveDir, 'Yearly_Corrected_R50s_TimeSeries.fig'));
    fprintf('Year-by-year corrected visualization saved to: %s\n', saveDir);
    
    % Print summary
    fprintf('\n=== Year-by-Year Correction Summary ===\n');
    fprintf('Total data points: %d\n', length(all_R50s_corrected));
    fprintf('Date range: %s to %s\n', datestr(min(all_dates), 'yyyy-mm-dd'), datestr(max(all_dates), 'yyyy-mm-dd'));
    fprintf('Original R50s - Mean: %.4f mΩ, Std: %.4f mΩ, Range: %.4f-%.4f mΩ\n', ...
        mean(all_R50s_uncorrected)*1000, std(all_R50s_uncorrected)*1000, ...
        min(all_R50s_uncorrected)*1000, max(all_R50s_uncorrected)*1000);
    fprintf('Corrected R50s - Mean: %.4f mΩ, Std: %.4f mΩ, Range: %.4f-%.4f mΩ\n', ...
        mean(all_R50s_corrected)*1000, std(all_R50s_corrected)*1000, ...
        min(all_R50s_corrected)*1000, max(all_R50s_corrected)*1000);
    
    % Check correlation with temperature (should be low after correction)
    corr_orig = corrcoef(all_R50s_uncorrected*1000, all_T);
    corr_corr = corrcoef(all_R50s_corrected*1000, all_T);
    fprintf('Correlation with Temperature:\n');
    fprintf('  Original R50s vs Temp: %.4f\n', corr_orig(1,2));
    fprintf('  Corrected R50s vs Temp: %.4f (should be close to 0)\n', corr_corr(1,2));
    
    % Year-by-year summary with detailed statistics
    fprintf('\n=== Year-by-Year Detailed Statistics ===\n');
    fprintf('Year\tPoints\tOriginal (mΩ)\t\t\tCorrected (mΩ)\t\t\tTemp Range\tCorr vs Temp\n');
    fprintf('\t\tMean±Std\tRange\t\tMean±Std\tRange\t\t(°C)\t\tOriginal\tCorrected\n');
    fprintf('----\t------\t----------------\t\t----------------\t\t--------\t--------\t--------\n');
    
    for i = 1:length(years)
        year = years{i};
        year_idx = strcmp(all_years, year);
        if sum(year_idx) > 0
            year_R50s_orig = all_R50s_uncorrected(year_idx);
            year_R50s_corr = all_R50s_corrected(year_idx);
            year_T = all_T(year_idx);
            
            % Calculate correlation with temperature
            corr_orig_year = corrcoef(year_R50s_orig*1000, year_T);
            corr_corr_year = corrcoef(year_R50s_corr*1000, year_T);
            
            fprintf('%s\t%d\t%.4f±%.4f\t%.4f-%.4f\t%.4f±%.4f\t%.4f-%.4f\t%.1f-%.1f\t%.4f\t%.4f\n', ...
                year, sum(year_idx), ...
                mean(year_R50s_orig)*1000, std(year_R50s_orig)*1000, ...
                min(year_R50s_orig)*1000, max(year_R50s_orig)*1000, ...
                mean(year_R50s_corr)*1000, std(year_R50s_corr)*1000, ...
                min(year_R50s_corr)*1000, max(year_R50s_corr)*1000, ...
                min(year_T), max(year_T), ...
                corr_orig_year(1,2), corr_corr_year(1,2));
        end
    end
    
    % Calculate correction factors
    fprintf('\n=== Temperature Correction Factor Analysis ===\n');
    all_correction_factors = [];
    for i = 1:length(years)
        year = years{i};
        field_name = ['year_' year];
        
        if isfield(yearly_models, field_name)
            R50s_data = yearly_data.(field_name).R50s;
            T_data = yearly_data.(field_name).T;
            model = yearly_models.(field_name);
            slope = model.slope;
            intercept = model.intercept;
            
            year_factors = zeros(size(R50s_data));
            for j = 1:length(R50s_data)
                T_i_kelvin = T_data(j) + 273.15;
                ln_R_ref25 = slope * (1/298.15) + intercept;
                ln_R_model_i = slope * (1/T_i_kelvin) + intercept;
                year_factors(j) = exp(ln_R_ref25) / exp(ln_R_model_i);
            end
            all_correction_factors = [all_correction_factors; year_factors];
            
            fprintf('Year %s: Correction factor range: %.4f to %.4f (Mean: %.4f)\n', ...
                year, min(year_factors), max(year_factors), mean(year_factors));
        end
    end
    fprintf('Overall correction factor range: %.4f to %.4f (Mean: %.4f)\n', ...
        min(all_correction_factors), max(all_correction_factors), mean(all_correction_factors));
    
    % Print outlier removal summary
    fprintf('\n=== Outlier Removal Summary ===\n');
    for i = 1:length(years)
        year = years{i};
        field_name = ['year_' year];
        if isfield(yearly_outlier_stats, field_name)
            stats = yearly_outlier_stats.(field_name);
            fprintf('Year %s: Removed %d outliers (%.1f%%) from %d total points\n', ...
                year, stats.n_outliers, stats.outlier_percentage, stats.n_total);
        end
    end
    fprintf('\n');
end

function createSOCControlledVisualization(saveDir)
    % Create SOC-controlled visualizations by grouping SOC into bins
    
    % Load temperature corrected data
    dataFile = fullfile(saveDir, 'Temperature_Corrected_R50s_Data.mat');
    if ~exist(dataFile, 'file')
        fprintf('Temperature corrected data file not found: %s\n', dataFile);
        return;
    end
    
    load(dataFile, 'corrected_data');
    
    dates = corrected_data.dates;
    R50s_uncorrected = corrected_data.R50s_uncorrected;
    R50s_corrected = corrected_data.R50s_corrected;
    SOC = corrected_data.SOC;
    T = corrected_data.T;
    years = corrected_data.years;
    
    % Filter out NaN SOC values
    valid_soc_idx = ~isnan(SOC);
    dates = dates(valid_soc_idx);
    R50s_uncorrected = R50s_uncorrected(valid_soc_idx);
    R50s_corrected = R50s_corrected(valid_soc_idx);
    SOC = SOC(valid_soc_idx);
    T = T(valid_soc_idx);
    years = years(valid_soc_idx);
    
    fprintf('Total data points with valid SOC: %d\n', length(SOC));
    fprintf('SOC range: %.1f%% to %.1f%%\n', min(SOC), max(SOC));
    
    % Define SOC groups
    soc_groups = {
        struct('name', 'SOC_50_55', 'min', 50, 'max', 55, 'label', 'SOC 50-55%');
        struct('name', 'SOC_55_55', 'min', 55, 'max', 60, 'label', 'SOC 55-60%');
        struct('name', 'SOC_60_65', 'min', 60, 'max', 65, 'label', 'SOC 60-65%');
        struct('name', 'SOC_65_70', 'min', 65, 'max', 70, 'label', 'SOC 65-70%');
    };
    
    % Create visualization for each SOC group
    for g = 1:length(soc_groups)
        group = soc_groups{g};
        
        % Filter data for this SOC group
        soc_idx = SOC >= group.min & SOC < group.max;
        
        if sum(soc_idx) == 0
            fprintf('No data found for %s\n', group.label);
            continue;
        end
        
        group_dates = dates(soc_idx);
        group_R50s_uncorrected = R50s_uncorrected(soc_idx);
        group_R50s_corrected = R50s_corrected(soc_idx);
        group_SOC = SOC(soc_idx);
        group_T = T(soc_idx);
        group_years = years(soc_idx);
        
        fprintf('\n=== %s ===\n', group.label);
        fprintf('Initial data points: %d\n', length(group_R50s_corrected));
        fprintf('Date range: %s to %s\n', datestr(min(group_dates), 'yyyy-mm-dd'), datestr(max(group_dates), 'yyyy-mm-dd'));
        fprintf('Initial R50s (uncorrected) range: %.4f to %.4f mΩ\n', min(group_R50s_uncorrected)*1000, max(group_R50s_uncorrected)*1000);
        fprintf('Initial R50s (corrected) range: %.4f to %.4f mΩ\n', min(group_R50s_corrected)*1000, max(group_R50s_corrected)*1000);
        fprintf('SOC range: %.1f%% to %.1f%%\n', min(group_SOC), max(group_SOC));
        
        % Apply year-by-year outlier removal within this SOC group (using corrected data)
        fprintf('\n  === Year-by-Year Outlier Removal (within %s) ===\n', group.label);
        z_threshold = 2.5; % 2.5-sigma rule
        
        unique_years_in_group = unique(group_years);
        all_outlier_idx = false(size(group_R50s_corrected));
        
        for y = 1:length(unique_years_in_group)
            year = unique_years_in_group{y};
            year_idx = strcmp(group_years, year);
            
            if sum(year_idx) > 1 % Need at least 2 points for outlier detection
                year_R50s = group_R50s_corrected(year_idx);
                
                % Calculate Z-Score for this year's data within this SOC group
                mean_val = mean(year_R50s);
                std_val = std(year_R50s);
                z_scores = abs((year_R50s - mean_val) / std_val);
                
                % Identify outliers
                year_outlier_idx = z_scores > z_threshold;
                
                % Count outliers
                n_outliers = sum(year_outlier_idx);
                n_total = length(year_R50s);
                
                if n_outliers > 0
                    fprintf('    Year %s: Mean: %.6f mΩ, Std: %.6f mΩ\n', ...
                        year, mean_val*1000, std_val*1000);
                    fprintf('      Outliers removed: %d out of %d (%.1f%%)\n', ...
                        n_outliers, n_total, n_outliers/n_total*100);
                    
                    % Mark outliers in the full group array
                    year_indices = find(year_idx);
                    all_outlier_idx(year_indices(year_outlier_idx)) = true;
                end
            end
        end
        
        % Remove outliers
        if sum(all_outlier_idx) > 0
            fprintf('  Total outliers removed: %d out of %d (%.1f%%)\n', ...
                sum(all_outlier_idx), length(group_R50s_corrected), sum(all_outlier_idx)/length(group_R50s_corrected)*100);
            
            group_dates = group_dates(~all_outlier_idx);
            group_R50s_uncorrected = group_R50s_uncorrected(~all_outlier_idx);
            group_R50s_corrected = group_R50s_corrected(~all_outlier_idx);
            group_SOC = group_SOC(~all_outlier_idx);
            group_T = group_T(~all_outlier_idx);
            group_years = group_years(~all_outlier_idx);
            
            fprintf('  Filtered data points: %d\n', length(group_R50s_corrected));
            fprintf('  Filtered R50s (corrected) range: %.4f to %.4f mΩ\n', min(group_R50s_corrected)*1000, max(group_R50s_corrected)*1000);
        else
            fprintf('  No outliers found.\n');
        end
        
        % Sort by date
        [group_dates, sortIdx] = sort(group_dates);
        group_R50s_uncorrected = group_R50s_uncorrected(sortIdx);
        group_R50s_corrected = group_R50s_corrected(sortIdx);
        group_SOC = group_SOC(sortIdx);
        group_T = group_T(sortIdx);
        group_years = group_years(sortIdx);
        
        % Create sequential x-axis
        xAxis = (1:length(group_dates))';
        
        % Create figure with two subplots: uncorrected and corrected
        figure('Position', [100, 100, 1400, 1000]);
        
        % Subplot 1: Uncorrected R50s
        subplot(2, 1, 1);
        scatter(xAxis, group_R50s_uncorrected*1000, 50, group_T, 'filled', 'HandleVisibility', 'off');
        hold on;
        
        % Add year division lines
        unique_years = unique(group_years);
        for i = 1:length(unique_years)
            year = unique_years{i};
            year_idx = strcmp(group_years, year);
            if sum(year_idx) > 0
                first_idx = find(year_idx, 1);
                xline(first_idx, '--', 'Color', 'k', 'LineWidth', 1, 'Alpha', 0.7, 'HandleVisibility', 'off');
            end
        end
        
        % Add trend line
        if length(group_R50s_uncorrected) > 1
            mdl = fitlm(xAxis, group_R50s_uncorrected*1000, 'Intercept', true);
            trendLine = mdl.Fitted;
            plot(xAxis, trendLine, 'g--', 'LineWidth', 3, 'DisplayName', sprintf('Trend (R²=%.4f)', mdl.Rsquared.Ordinary));
        end
        
        % Colorbar for battery temperature
        ax1 = gca;
        colormap(ax1, flipud(autumn));
        c1 = colorbar(ax1, 'eastoutside');
        c1.Label.String = 'Battery Temperature (°C)';
        c1.Label.FontSize = 12;
        c1.Label.FontWeight = 'bold';
        
        % Set colorbar limits safely
        if length(group_T) > 0
            temp_min = min(group_T);
            temp_max = max(group_T);
            if temp_min == temp_max
                c1.Limits = [temp_min - 0.5, temp_max + 0.5];
                c1.Ticks = round(temp_min);
            else
                c1.Limits = [temp_min, temp_max];
                c1.Ticks = round(temp_min):1:round(temp_max);
            end
        end
        
        % Formatting
        xlabel('Date (YYYY:MM:DD)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('R50s (mΩ) - Uncorrected', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('ESS Rack01 Daily R50s Time Series (%s, Uncorrected)', group.label), ...
            'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        ax1.LineWidth = 2;
        ax1.FontWeight = 'bold';
        legend('Location', 'best', 'FontSize', 12);
        
        % Format x-axis
        dateYears = zeros(size(group_dates));
        dateMonths = zeros(size(group_dates));
        for i = 1:length(group_dates)
            dateYears(i) = group_dates(i).Year;
            dateMonths(i) = group_dates(i).Month;
        end
        
        uniqueMonths = unique(dateYears*100 + dateMonths);
        monthLabels = {};
        monthPositions = [];
        
        for i = 1:length(uniqueMonths)
            monthVal = uniqueMonths(i);
            yearVal = floor(monthVal/100);
            monthVal = mod(monthVal, 100);
            
            monthIdx = find(dateYears == yearVal & dateMonths == monthVal, 1);
            if ~isempty(monthIdx)
                monthLabels{end+1} = sprintf('%d-%02d', yearVal, monthVal);
                monthPositions(end+1) = monthIdx;
            end
        end
        
        xticks(monthPositions);
        xticklabels(monthLabels);
        xtickangle(45);
        
        % Subplot 2: Corrected R50s
        subplot(2, 1, 2);
        scatter(xAxis, group_R50s_corrected*1000, 50, group_T, 'filled', 'HandleVisibility', 'off');
        hold on;
        
        % Add year division lines
        for i = 1:length(unique_years)
            year = unique_years{i};
            year_idx = strcmp(group_years, year);
            if sum(year_idx) > 0
                first_idx = find(year_idx, 1);
                xline(first_idx, '--', 'Color', 'k', 'LineWidth', 1, 'Alpha', 0.7, 'HandleVisibility', 'off');
            end
        end
        
        % Add trend line
        if length(group_R50s_corrected) > 1
            mdl = fitlm(xAxis, group_R50s_corrected*1000, 'Intercept', true);
            trendLine = mdl.Fitted;
            plot(xAxis, trendLine, 'g--', 'LineWidth', 3, 'DisplayName', sprintf('Trend (R²=%.4f)', mdl.Rsquared.Ordinary));
        end
        
        % Colorbar for battery temperature
        ax2 = gca;
        colormap(ax2, flipud(autumn));
        c2 = colorbar(ax2, 'eastoutside');
        c2.Label.String = 'Battery Temperature (°C)';
        c2.Label.FontSize = 12;
        c2.Label.FontWeight = 'bold';
        
        % Set colorbar limits safely
        if length(group_T) > 0
            temp_min = min(group_T);
            temp_max = max(group_T);
            if temp_min == temp_max
                c2.Limits = [temp_min - 0.5, temp_max + 0.5];
                c2.Ticks = round(temp_min);
            else
                c2.Limits = [temp_min, temp_max];
                c2.Ticks = round(temp_min):1:round(temp_max);
            end
        end
        
        % Formatting
        xlabel('Date (YYYY:MM:DD)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('R50s (mΩ) - Temperature Corrected to 25°C', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('ESS Rack01 Daily R50s Time Series (%s, Temperature Corrected)', group.label), ...
            'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        ax2.LineWidth = 2;
        ax2.FontWeight = 'bold';
        legend('Location', 'best', 'FontSize', 12);
        
        % Format x-axis
        xticks(monthPositions);
        xticklabels(monthLabels);
        xtickangle(45);
        
        % Save figure
        saveFileName = sprintf('Yearly_Corrected_R50s_TimeSeries_%s.fig', group.name);
        saveas(gcf, fullfile(saveDir, saveFileName));
        fprintf('SOC-controlled visualization saved to: %s\n', saveFileName);
        
        % Print statistics
        fprintf('Statistics for %s:\n', group.label);
        fprintf('  Uncorrected R50s - Mean: %.4f mΩ, Std: %.4f mΩ\n', mean(group_R50s_uncorrected)*1000, std(group_R50s_uncorrected)*1000);
        fprintf('  Corrected R50s - Mean: %.4f mΩ, Std: %.4f mΩ\n', mean(group_R50s_corrected)*1000, std(group_R50s_corrected)*1000);
        fprintf('  Mean SOC: %.1f%%, Std: %.1f%%\n', mean(group_SOC), std(group_SOC));
        
        % Year-by-year statistics
        fprintf('  Year-by-year statistics (Corrected):\n');
        for i = 1:length(unique_years)
            year = unique_years{i};
            year_idx = strcmp(group_years, year);
            if sum(year_idx) > 0
                year_R50s_corr = group_R50s_corrected(year_idx);
                year_R50s_uncorr = group_R50s_uncorrected(year_idx);
                fprintf('    Year %s: %d points, Uncorrected: %.4f±%.4f mΩ, Corrected: %.4f±%.4f mΩ\n', ...
                    year, sum(year_idx), mean(year_R50s_uncorr)*1000, std(year_R50s_uncorr)*1000, ...
                    mean(year_R50s_corr)*1000, std(year_R50s_corr)*1000);
            end
        end
    end
    
    fprintf('\nSOC-controlled visualizations completed!\n');
end
