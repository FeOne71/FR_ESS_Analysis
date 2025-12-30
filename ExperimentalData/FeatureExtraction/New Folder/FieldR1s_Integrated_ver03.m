%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldR1s_Integrated_ver03.m
% ESS Rack01 R1s for Integrated Data Processing with Moving Average
% Combines OldData and NewData processing with 30-second moving average
% for noise reduction while preserving signal integrity
% 
% Key Features:
% 1. Integrated processing of both OldData and NewData
% 2. 30-second moving average for V and I signals
% 3. Noise removal while preserving meaningful signal changes
% 4. Daily median R1s calculation for robust statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

%% Configuration
% Data directories
oldDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old';
newDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New';

% Save directory
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_R1s','R1s_Results_Integrated_ver03');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameters
% R1s calculation parameters
Cnom = 128;
Cnom_cell = 64;                   % Rack nominal Capacity (Ah)
idle_thr = Cnom_cell*0.05;         % Idle threshold [charge, discharge] (A)

% Moving average parameters
window_size = 30;                 % 10-second moving average window

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
fprintf('=== Processing OldData ===\n');
for year_idx = 1:length(oldDataYears)
    year = oldDataYears{year_idx};
    fprintf('\nProcessing OldData year: %s\n', year);
    
    % Initialize R1s results structure for current year
    R1s_Results = struct();
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
                    
                    % Calculate dV and dI from smoothed signals
                    dI = diff(I_smooth);
                    dV = diff(V_smooth);
                    T_inst = T_batt(1:end-1); % Temperature not smoothed
                    
                    % Debug: Check dV and dI ranges
                    fprintf('    dI range: %.6f to %.6f A\n', min(dI), max(dI));
                    fprintf('    dV range: %.6f to %.6f V\n', min(dV), max(dV));
                    
                    % Filter valid data points (NaN, dI=0, and dV=0 removal)
                    valid_idx = ~isnan(dI) & ~isnan(dV) & dI ~= 0 & dV ~= 0;
                    dI_valid = dI(valid_idx);
                    dV_valid = dV(valid_idx);
                    T_inst_valid = T_inst(valid_idx);
                    
                    % Calculate R1s using linear fitting: dV = a * dI (y-intercept = 0)
                    if length(dI_valid) >= 10 % Need sufficient data points for reliable fitting
                        % Linear fitting without y-intercept: dV = a * dI
                        p = polyfit(dI_valid, dV_valid, 1);
                        R1s_daily = p(1); % Fitted slope 'a'
                        
                        % Calculate R-squared for model validation
                        dV_predicted = R1s_daily * dI_valid;
                        ss_res = sum((dV_valid - dV_predicted).^2);
                        ss_tot = sum((dV_valid - mean(dV_valid)).^2);
                        r_squared = 1 - (ss_res / ss_tot);
                        
                        % Store results
                        dayKey = matFiles(f).name(1:end-4); % Remove .mat extension
                        R1s_Results.(dayKey).(rackName).R1s_daily = R1s_daily; % Daily R1s (slope)
                        R1s_Results.(dayKey).(rackName).dI_materials = dI_valid; % Store dI materials
                        R1s_Results.(dayKey).(rackName).dV_materials = dV_valid; % Store dV materials
                        R1s_Results.(dayKey).(rackName).T_inst = T_inst_valid;
                        R1s_Results.(dayKey).(rackName).n_points = length(dI_valid);
                        R1s_Results.(dayKey).(rackName).r_squared = r_squared;
                        R1s_Results.(dayKey).(rackName).source = 'OldData';
                        R1s_Results.(dayKey).(rackName).window_size = window_size;
                        
                        fprintf('R1s_daily = %.6f Ω (%.3f mΩ), R² = %.4f, n_points = %d, T range: %.1f-%.1f°C\n', ...
                            R1s_daily, R1s_daily*1000, r_squared, length(dI_valid), min(T_inst_valid), max(T_inst_valid));
                    else
                        fprintf('Insufficient data points (%d) for reliable R1s calculation\n', length(dI_valid));
                    end
                end
                
            catch ME
                fprintf('ERROR processing %s: %s\n', matFiles(f).name, ME.message);
                continue;
            end
        end % End mat file loop
    end % End month loop
    
    % Save results for current year
    if ~isempty(fieldnames(R1s_Results))
        saveFileName = sprintf('R1s_Results_OldData_%s.mat', year);
        save(fullfile(saveDir, saveFileName), 'R1s_Results');
        fprintf('Saved OldData results to: %s\n', saveFileName);
        
        % Create summary statistics
        createYearlySummary(R1s_Results, year, 'OldData', rackNames_all);
    else
        fprintf('No valid OldData found for %s\n', year);
    end
end % End OldData year loop

%% Process NewData
fprintf('\n=== Processing NewData ===\n');
for year_idx = 1:length(newDataYears)
    year = newDataYears{year_idx};
    fprintf('\nProcessing NewData year: %s\n', year);
    
    % Initialize R1s results structure for current year
    R1s_Results = struct();
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
                    
                    % Calculate dV and dI from smoothed signals
                    dI = diff(I_smooth);
                    dV = diff(V_smooth);
                    T_inst = T_batt(1:end-1); % Temperature not smoothed
                    
                    % Debug: Check dV and dI ranges
                    fprintf('    dI range: %.6f to %.6f A\n', min(dI), max(dI));
                    fprintf('    dV range: %.6f to %.6f V\n', min(dV), max(dV));
                    
                    % Filter valid data points (NaN, dI=0, and dV=0 removal)
                    valid_idx = ~isnan(dI) & ~isnan(dV) & dI ~= 0 & dV ~= 0;
                    dI_valid = dI(valid_idx);
                    dV_valid = dV(valid_idx);
                    T_inst_valid = T_inst(valid_idx);
                    
                    % Calculate R1s using linear fitting: dV = a * dI (y-intercept = 0)
                    if length(dI_valid) >= 10 % Need sufficient data points for reliable fitting
                        % Linear fitting without y-intercept: dV = a * dI
                        p = polyfit(dI_valid, dV_valid, 1);
                        R1s_daily = p(1); % Fitted slope 'a'
                        
                        % Calculate R-squared for model validation
                        dV_predicted = R1s_daily * dI_valid;
                        ss_res = sum((dV_valid - dV_predicted).^2);
                        ss_tot = sum((dV_valid - mean(dV_valid)).^2);
                        r_squared = 1 - (ss_res / ss_tot);
                        
                        % Store results
                        dayKey = matFiles(f).name(1:end-4); % Remove .mat extension
                        R1s_Results.(dayKey).(rackName).R1s_daily = R1s_daily; % Daily R1s (slope)
                        R1s_Results.(dayKey).(rackName).dI_materials = dI_valid; % Store dI materials
                        R1s_Results.(dayKey).(rackName).dV_materials = dV_valid; % Store dV materials
                        R1s_Results.(dayKey).(rackName).T_inst = T_inst_valid;
                        R1s_Results.(dayKey).(rackName).n_points = length(dI_valid);
                        R1s_Results.(dayKey).(rackName).r_squared = r_squared;
                        R1s_Results.(dayKey).(rackName).source = 'NewData';
                        R1s_Results.(dayKey).(rackName).window_size = window_size;
                        
                        fprintf('R1s_daily = %.6f Ω (%.3f mΩ), R² = %.4f, n_points = %d, T range: %.1f-%.1f°C\n', ...
                            R1s_daily, R1s_daily*1000, r_squared, length(dI_valid), min(T_inst_valid), max(T_inst_valid));
                    else
                        fprintf('Insufficient data points (%d) for reliable R1s calculation\n', length(dI_valid));
                    end
                end
                
            catch ME
                fprintf('ERROR processing %s: %s\n', matFiles(f).name, ME.message);
                continue;
            end
        end % End mat file loop
    end % End month loop
    
    % Save results for current year
    if ~isempty(fieldnames(R1s_Results))
        saveFileName = sprintf('R1s_Results_NewData_%s.mat', year);
        save(fullfile(saveDir, saveFileName), 'R1s_Results');
        fprintf('Saved NewData results to: %s\n', saveFileName);
        
        % Create summary statistics
        createYearlySummary(R1s_Results, year, 'NewData', rackNames_all);
    else
        fprintf('No valid NewData found for %s\n', year);
    end
end % End NewData year loop

%% Create ln(R1s) vs 1/T Visualization by Year
fprintf('\n=== Creating ln(R1s) vs 1/T Visualization by Year ===\n');
createR1sVsTemperaturePlot(saveDir);

%% Year-by-Year Temperature Correction and Final Visualization
fprintf('\n=== Year-by-Year Temperature Correction and Final Visualization ===\n');
createYearlyCorrectedVisualization(saveDir);

fprintf('\nIntegrated R1s processing with moving average completed!\n');
fprintf('Results saved to: %s\n', saveDir);

%% Helper Functions

function createYearlySummary(R1s_Results, year, dataType, rackNames_all)
    % Create summary statistics for a year
    allR1s_daily = [];
    allR_squared = [];
    allN = [];
    
    dayNames = fieldnames(R1s_Results);
    for d = 1:length(dayNames)
        if isfield(R1s_Results.(dayNames{d}), rackNames_all{1})
            allR1s_daily = [allR1s_daily; R1s_Results.(dayNames{d}).(rackNames_all{1}).R1s_daily];
            allR_squared = [allR_squared; R1s_Results.(dayNames{d}).(rackNames_all{1}).r_squared];
            allN = [allN; R1s_Results.(dayNames{d}).(rackNames_all{1}).n_points];
        end
    end
    
    if ~isempty(allR1s_daily)
        fprintf('\n=== %s %s Summary ===\n', year, dataType);
        fprintf('Total days processed: %d\n', length(dayNames));
        fprintf('Mean R1s_daily: %.4f mΩ\n', mean(allR1s_daily)*1000);
        fprintf('Std R1s_daily: %.4f mΩ\n', std(allR1s_daily)*1000);
        fprintf('Mean R²: %.4f\n', mean(allR_squared));
        fprintf('Mean materials per day: %.1f\n', mean(allN));
    end
end

function createR1sVsTemperaturePlot(saveDir)
    % Create R1s vs 1/T visualization by year
    
    % Find all saved result files
    resultFiles = dir(fullfile(saveDir, 'R1s_Results_*.mat'));
    
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
            load(filePath, 'R1s_Results');
            
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
            
            % Extract R1s and temperature data
            dayNames = fieldnames(R1s_Results);
            year_R1s = [];
            year_T = [];
            
            for d = 1:length(dayNames)
                if isfield(R1s_Results.(dayNames{d}), 'Rack01')
                    rackData = R1s_Results.(dayNames{d}).Rack01;
                    if isfield(rackData, 'R1s_daily') && isfield(rackData, 'T_inst')
                        year_R1s = [year_R1s; rackData.R1s_daily];
                        year_T = [year_T; mean(rackData.T_inst)]; % Use mean temperature for the day
                    end
                end
            end
            
            % Store yearly data
            field_name = ['year_' year];
            if isfield(yearly_data, field_name)
                % Merge with existing data
                yearly_data.(field_name).R1s = [yearly_data.(field_name).R1s; year_R1s];
                yearly_data.(field_name).T = [yearly_data.(field_name).T; year_T];
                yearly_data.(field_name).source = 'Mixed';
            else
                % New year
                yearly_data.(field_name).R1s = year_R1s;
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
    
    % Create figure
    figure('Position', [100, 100, 1400, 1000]);
    
    % Calculate subplot layout
    nYears = length(years);
    nCols = min(3, nYears);
    nRows = ceil(nYears / nCols);
    
    % Calculate global axis limits for consistent scaling
    all_inv_T = [];
    all_ln_R1s = [];
    for i = 1:nYears
        year = years{i};
        field_name = ['year_' year];
        R1s_data = yearly_data.(field_name).R1s;
        T_data = yearly_data.(field_name).T;
        
        if length(R1s_data) > 1
            T_K = T_data + 273.15;
            inv_T = 1 ./ T_K;
            ln_R1s = log(R1s_data * 1000);
            
            all_inv_T = [all_inv_T; inv_T];
            all_ln_R1s = [all_ln_R1s; ln_R1s];
        end
    end
    
    % Set global axis limits with small margins
    x_margin = (max(all_inv_T) - min(all_inv_T)) * 0.05;
    y_margin = (max(all_ln_R1s) - min(all_ln_R1s)) * 0.05;
    x_lim = [min(all_inv_T) - x_margin, max(all_inv_T) + x_margin];
    y_lim = [min(all_ln_R1s) - y_margin, max(all_ln_R1s) + y_margin];
    
    for i = 1:nYears
        year = years{i};
        field_name = ['year_' year];
        
        subplot(nRows, nCols, i);
        
        R1s_data = yearly_data.(field_name).R1s;
        T_data = yearly_data.(field_name).T;
        source = yearly_data.(field_name).source;
        
        % Convert temperature to 1/T (K^-1)
        T_K = T_data + 273.15; % Convert to Kelvin
        inv_T = 1 ./ T_K;
        
        % Convert R1s to ln(R1s)
        ln_R1s = log(R1s_data * 1000); % Convert to mΩ and take natural log
        
        % Create scatter plot
        scatter(inv_T, ln_R1s, 50, T_data, 'filled');
        
        % Add linear trend line
        if length(R1s_data) > 1
            mdl = fitlm(inv_T, ln_R1s);
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
                slope, intercept, r_squared, length(R1s_data), temp_mean, temp_std, temp_min, temp_max), ...
                'Units', 'normalized', 'VerticalAlignment', 'top', ...
                'BackgroundColor', 'white', 'EdgeColor', 'black');
        end
        
        % Formatting
        xlabel('1/T (K^{-1})', 'FontSize', 10);
        ylabel('ln(R1s) (ln(mΩ))', 'FontSize', 10);
        title(sprintf('Year %s', year), 'FontSize', 12);
        colorbar;
        colormap(flipud(autumn));
        caxis([min(T_data), max(T_data)]);
        grid on;
        
        % Set consistent axis limits for all subplots
        xlim(x_lim);
        ylim(y_lim);
        
        % Data info is now included in the trend line info above
    end
    
    % Add overall title
    sgtitle('ln(R1s) vs 1/T by Year (Arrhenius Relationship)', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save figure
    saveas(gcf, fullfile(saveDir, 'R1s_vs_1T_by_Year.fig'));    
    fprintf('ln(R1s) vs 1/T visualization saved to: %s\n', saveDir);
    
    % Print yearly slope, intercept, R², and temperature distribution summary
    fprintf('\n=== Yearly Arrhenius Analysis Summary ===\n');
    fprintf('Year\t\tSlope (K)\tIntercept\tR²\t\tData Points\tTemp Range (°C)\tTemp Mean±Std (°C)\n');
    fprintf('----\t\t---------\t---------\t---\t\t-----------\t---------------\t------------------\n');
    
    for i = 1:nYears
        year = years{i};
        field_name = ['year_' year];
        
        R1s_data = yearly_data.(field_name).R1s;
        T_data = yearly_data.(field_name).T;
        
        if length(R1s_data) > 1
            % Convert temperature to 1/T (K^-1)
            T_K = T_data + 273.15;
            inv_T = 1 ./ T_K;
            
            % Convert R1s to ln(R1s)
            ln_R1s = log(R1s_data * 1000);
            
            % Calculate slope, intercept, and R² using MATLAB functions
            mdl = fitlm(inv_T, ln_R1s);
            slope = mdl.Coefficients.Estimate(2);
            intercept = mdl.Coefficients.Estimate(1);
            r_squared = mdl.Rsquared.Ordinary;
            
            % Calculate temperature statistics
            temp_min = min(T_data);
            temp_max = max(T_data);
            temp_mean = mean(T_data);
            temp_std = std(T_data);
            
            fprintf('%s\t\t%.2f\t\t%.2f\t\t%.4f\t\t%d\t\t%.1f-%.1f\t\t%.1f±%.1f\n', ...
                year, slope, intercept, r_squared, length(R1s_data), temp_min, temp_max, temp_mean, temp_std);
        else
            fprintf('%s\t\tN/A\t\tN/A\t\tN/A\t\t%d\t\tN/A\t\tN/A\n', year, length(R1s_data));
        end
    end
    fprintf('\n');
end

function createYearlyCorrectedVisualization(saveDir)
    % Create year-by-year temperature corrected visualization
    
    % Find all saved result files
    resultFiles = dir(fullfile(saveDir, 'R1s_Results_*.mat'));
    
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
            load(filePath, 'R1s_Results');
            
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
            
            % Extract R1s and temperature data
            dayNames = fieldnames(R1s_Results);
            year_R1s = [];
            year_T = [];
            year_dates = [];
            
            for d = 1:length(dayNames)
                if isfield(R1s_Results.(dayNames{d}), 'Rack01')
                    rackData = R1s_Results.(dayNames{d}).Rack01;
                    if isfield(rackData, 'R1s_daily') && isfield(rackData, 'T_inst')
                        year_R1s = [year_R1s; rackData.R1s_daily];
                        year_T = [year_T; mean(rackData.T_inst)]; % Use mean temperature for the day
                        
                        % Extract date from dayKey
                        dayKey = dayNames{d};
                        if length(dayKey) >= 12 && strcmp(dayKey(1:4), 'Raw_')
                            dateStr = dayKey(5:12); % '20210601'
                            year_num = str2double(dateStr(1:4));
                            month = str2double(dateStr(5:6));
                            day = str2double(dateStr(7:8));
                            date = datetime(year_num, month, day);
                            year_dates = [year_dates; date];
                        end
                    end
                end
            end
            
            % Store yearly data
            field_name = ['year_' year];
            if isfield(yearly_data, field_name)
                % Merge with existing data
                yearly_data.(field_name).R1s = [yearly_data.(field_name).R1s; year_R1s];
                yearly_data.(field_name).T = [yearly_data.(field_name).T; year_T];
                yearly_data.(field_name).dates = [yearly_data.(field_name).dates; year_dates];
                yearly_data.(field_name).source = 'Mixed';
            else
                % New year
                yearly_data.(field_name).R1s = year_R1s;
                yearly_data.(field_name).T = year_T;
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
    
    % Create yearly models
    yearly_models = struct();
    Tref = 25 + 273.15; % Reference temperature: 25°C in Kelvin
    
    fprintf('Creating yearly Arrhenius models...\n');
    for i = 1:length(years)
        year = years{i};
        field_name = ['year_' year];
        
        R1s_data = yearly_data.(field_name).R1s;
        T_data = yearly_data.(field_name).T;
        
        if length(R1s_data) > 1
            % Convert temperature to 1/T (K^-1)
            T_K = T_data + 273.15;
            inv_T = 1 ./ T_K;
            
            % Convert R1s to ln(R1s)
            ln_R1s = log(R1s_data * 1000);
            
            % Calculate slope and intercept using MATLAB functions
            mdl = fitlm(inv_T, ln_R1s);
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
    all_R1s_uncorrected = [];
    all_R1s_corrected = [];
    all_T = [];
    all_years = [];
    
    for i = 1:length(years)
        year = years{i};
        field_name = ['year_' year];
        
        if isfield(yearly_models, field_name)
            R1s_data = yearly_data.(field_name).R1s;
            T_data = yearly_data.(field_name).T;
            dates = yearly_data.(field_name).dates;
            
            % Apply year-specific temperature correction
            model = yearly_models.(field_name);
            slope = model.slope;
            intercept = model.intercept;
            
            R1s_corrected = zeros(size(R1s_data));
            for j = 1:length(R1s_data)
                % Current temperature in Kelvin
                T_i_kelvin = T_data(j) + 273.15;
                
                % Correction formula
                ln_R_ref25 = slope * (1/298.15) + intercept;
                ln_R_model_i = slope * (1/T_i_kelvin) + intercept;
                Factor_i = exp(ln_R_ref25) / exp(ln_R_model_i);
                
                % Apply correction
                R1s_corrected(j) = R1s_data(j) * Factor_i;
            end
            
            % Store corrected data
            all_dates = [all_dates; dates];
            all_R1s_uncorrected = [all_R1s_uncorrected; R1s_data];
            all_R1s_corrected = [all_R1s_corrected; R1s_corrected];
            all_T = [all_T; T_data];
            all_years = [all_years; repmat({year}, length(R1s_data), 1)];
        end
    end
    
    % Sort by date
    [all_dates, sortIdx] = sort(all_dates);
    all_R1s_uncorrected = all_R1s_uncorrected(sortIdx);
    all_R1s_corrected = all_R1s_corrected(sortIdx);
    all_T = all_T(sortIdx);
    all_years = all_years(sortIdx);
    
    % Create sequential x-axis
    xAxis = (1:length(all_dates))';
    
    % Create final visualization
    figure('Position', [100, 100, 1400, 800]);
    
    % Plot all data as scatter points with battery temperature color
    scatter(xAxis, all_R1s_corrected*1000, 50, all_T, 'filled', 'HandleVisibility', 'off');
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
    if length(all_R1s_corrected) > 1
        mdl = fitlm(xAxis, all_R1s_corrected*1000, 'Intercept', true);
        trendLine = mdl.Fitted;
        plot(xAxis, trendLine, 'g--', 'LineWidth', 3, 'DisplayName', sprintf('Trend (R²=%.4f)', mdl.Rsquared.Ordinary));
    end
    
    % Colorbar for battery temperature (ver04 style)
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
    ylabel('R1s (mΩ) - Temperature Corrected to 25°C', 'FontSize', 12, 'FontWeight', 'bold');
    title('ESS Rack01 Daily R1s Time Series (Year-by-Year Temperature Corrected)', 'FontSize', 14, 'FontWeight', 'bold');
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
    saveas(gcf, fullfile(saveDir, 'Yearly_Corrected_R1s_TimeSeries.fig'));
    fprintf('Year-by-year corrected visualization saved to: %s\n', saveDir);
    
    % Print summary
    fprintf('\n=== Year-by-Year Correction Summary ===\n');
    fprintf('Total data points: %d\n', length(all_R1s_corrected));
    fprintf('Date range: %s to %s\n', datestr(min(all_dates), 'yyyy-mm-dd'), datestr(max(all_dates), 'yyyy-mm-dd'));
    fprintf('Original R1s - Mean: %.4f mΩ, Std: %.4f mΩ\n', mean(all_R1s_uncorrected)*1000, std(all_R1s_uncorrected)*1000);
    fprintf('Corrected R1s - Mean: %.4f mΩ, Std: %.4f mΩ\n', mean(all_R1s_corrected)*1000, std(all_R1s_corrected)*1000);
    
    % Debug: Find extreme values
    fprintf('\n=== Extreme Values Analysis ===\n');
    [sorted_orig, orig_idx] = sort(all_R1s_uncorrected*1000, 'descend');
    [sorted_corr, corr_idx] = sort(all_R1s_corrected*1000, 'descend');
    
    fprintf('Top 10 Original R1s values:\n');
    for i = 1:min(10, length(sorted_orig))
        idx = orig_idx(i);
        fprintf('  %.4f mΩ on %s (Year: %s, Temp: %.1f°C)\n', ...
            sorted_orig(i), datestr(all_dates(idx), 'yyyy-mm-dd'), ...
            all_years{idx}, all_T(idx));
    end
    
    fprintf('\nTop 10 Corrected R1s values:\n');
    for i = 1:min(10, length(sorted_corr))
        idx = corr_idx(i);
        fprintf('  %.4f mΩ on %s (Year: %s, Temp: %.1f°C)\n', ...
            sorted_corr(i), datestr(all_dates(idx), 'yyyy-mm-dd'), ...
            all_years{idx}, all_T(idx));
    end
    
    % Year 2023 specific analysis
    year_2023_idx = strcmp(all_years, '2023');
    if sum(year_2023_idx) > 0
        year_2023_orig = all_R1s_uncorrected(year_2023_idx)*1000;
        year_2023_corr = all_R1s_corrected(year_2023_idx)*1000;
        year_2023_dates = all_dates(year_2023_idx);
        year_2023_temp = all_T(year_2023_idx);
        
        fprintf('\n=== Year 2023 Extreme Values ===\n');
        [sorted_2023_orig, idx_2023_orig] = sort(year_2023_orig, 'descend');
        fprintf('Top 5 Original R1s values in 2023:\n');
        for i = 1:min(5, length(sorted_2023_orig))
            idx = idx_2023_orig(i);
            fprintf('  %.4f mΩ on %s (Temp: %.1f°C)\n', ...
                sorted_2023_orig(i), datestr(year_2023_dates(idx), 'yyyy-mm-dd'), year_2023_temp(idx));
        end
    end
    
    % Year-by-year summary
    for i = 1:length(years)
        year = years{i};
        year_idx = strcmp(all_years, year);
        if sum(year_idx) > 0
            year_R1s_orig = all_R1s_uncorrected(year_idx);
            year_R1s_corr = all_R1s_corrected(year_idx);
            fprintf('Year %s: %d points, Original: %.4f±%.4f mΩ, Corrected: %.4f±%.4f mΩ\n', ...
                year, sum(year_idx), mean(year_R1s_orig)*1000, std(year_R1s_orig)*1000, ...
                mean(year_R1s_corr)*1000, std(year_R1s_corr)*1000);
        end
    end
end


