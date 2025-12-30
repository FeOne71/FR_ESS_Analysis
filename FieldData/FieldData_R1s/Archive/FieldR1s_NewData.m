%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldR1s_NewData.m
% ESS Rack01 R1s for New Data
% y=ax
% x=dI, y=dV, a=R1s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New';
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_R1s','R1s_Results_NewData');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameters
% R1s calculation parameters
Cnom = 128;
Cnom_cell = 64;                   % Rack nominal Capacity (Ah)
idle_thr = Cnom_cell*0.05;         % Idle threshold [charge, discharge] (A)
% dt will be calculated based on data characteristics

% Topology
Ns = 17*14;    % 238s
Np = 2;        % 2p

% Visualization parameters
Fontsize = 12;
LineWidth = 2;

%% Load New Data
yearList = {'2025'}; % Process
rackNames_all = {'Rack01'};

%% Process each year
for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    fprintf('Processing year: %s\n', year);

    % Initialize R1s results structure for current year
    R1s_Results = struct();
    yearPath = fullfile(dataDir, year);
    
    fprintf('Year path: %s\n', yearPath);
    fprintf('Year path exists: %d\n', exist(yearPath, 'dir'));
    
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
        
        % Initialize monthly table
        monthly_table = table();
        month_name = monthDirs(m).name;

        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            fprintf('Loading file: %s\n', matFiles(f).name);
            
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

            rackNames = rackNames_all;

            for rack_idx = 1:length(rackNames)
                rackName = rackNames{rack_idx};
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
                
                % Debug: check data type
                fprintf('Time data type: %s, size: %s\n', class(t), mat2str(size(t)));
                if length(t) > 0
                    fprintf('First time value: %s\n', string(t(1)));
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
                
                % Additional signals
                soc = D.SOC_BMS(:);
                P = D.DCPower(:);
                T_batt = D.MTavg(:);
                % T_amb = D.AmbientTemp(:);

    % To cell units
                I = I / Np;  % Convert to cell current
    
                % Calculate dV/dI for R1s estimation
                fprintf('Calculating R1s for %s...\n', rackName);
                
                % Calculate dI and dV using continuous segments only
                dt = diff(t);
                
                % Find continuous segments (remove large gaps)
                expected_interval = mode(dt);
                fprintf('Expected interval = %.1f seconds\n', expected_interval);
                
                gap_indices = find(dt > expected_interval * 1.5); % 1.5배 이상 간격은 구간 분리
                
                % Define segment boundaries
                if isempty(gap_indices)
                    % No gaps found, use entire data as one segment
                    segment_starts = [1];
                    segment_ends = [length(t)];
                else
                    segment_starts = [1; gap_indices + 1];
                    segment_ends = [gap_indices; length(t)];
                end
                
                % Calculate dI and dV for each continuous segment
                dI = [];
                dV = [];
                soc_valid = [];
                P_valid = [];
                T_batt_valid = [];
                % T_amb_valid = [];
                
                for i = 1:length(segment_starts)
                    start_idx = segment_starts(i);
                    end_idx = segment_ends(i);
                    
                    if end_idx > start_idx % Need at least 2 data points
                        dI = [dI; diff(I(start_idx:end_idx))];
                        dV = [dV; diff(V(start_idx:end_idx))];
                        % Store corresponding SOC, Power, Battery Temp, Ambient Temp values
                        soc_valid = [soc_valid; soc(start_idx:end_idx-1)]; % diff reduces length by 1
                        P_valid = [P_valid; P(start_idx:end_idx-1)];
                        T_batt_valid = [T_batt_valid; T_batt(start_idx:end_idx-1)];
                        % T_amb_valid = [T_amb_valid; T_amb(start_idx:end_idx-1)];
                    end
                end
                
                % Remove cases where dI or dV is NaN
                validIdx = ~isnan(dI) & ~isnan(dV);
                dI = dI(validIdx);
                dV = dV(validIdx);
                soc_valid = soc_valid(validIdx);
                P_valid = P_valid(validIdx);
                T_batt_valid = T_batt_valid(validIdx);
                % T_amb_valid = T_amb_valid(validIdx);
                
                % Check if we have enough data for calculation
                if length(dI) < 5
                    fprintf('Insufficient data for R1s calculation\n');
                    continue;
                end
                
                % Linear regression: y = ax (dV = R1s * dI) using fitlm
                mdl = fitlm(dI, dV, 'Intercept', false);
                R1s = mdl.Coefficients.Estimate(1); % Slope (no intercept)
                R_squared = mdl.Rsquared.Ordinary;
                
                % Store results in structure
                dayKey = matFiles(f).name(1:end-4); % Remove .mat extension
                R1s_Results.(dayKey).(rackName).R1s = R1s;
                R1s_Results.(dayKey).(rackName).R_squared = R_squared;
                R1s_Results.(dayKey).(rackName).n_points = length(dI);
                R1s_Results.(dayKey).(rackName).dI = dI;
                R1s_Results.(dayKey).(rackName).dV = dV;
                R1s_Results.(dayKey).(rackName).time = t(2:end);
                R1s_Results.(dayKey).(rackName).dt = expected_interval;
                
                % Store additional data (only for dI/dV calculation points)
                R1s_Results.(dayKey).(rackName).SOC = soc_valid;
                R1s_Results.(dayKey).(rackName).Power = P_valid;
                R1s_Results.(dayKey).(rackName).BatteryTemp = T_batt_valid;
                % R1s_Results.(dayKey).(rackName).AmbientTemp = T_amb_valid;
                
                % Add to monthly table
                new_row = table({dayKey}, R1s, R_squared, length(dI), expected_interval, {dI}, {dV}, {t(2:end)}, ...
                    'VariableNames', {'Date', 'R1s', 'R_squared', 'N_points', 'dt', 'dI', 'dV', 'time'});
                monthly_table = [monthly_table; new_row];
                
                fprintf('R1s = %.4f mOhm, R² = %.3f, n = %d\n', R1s*1000, R_squared, length(dI));
                
                % SOC target analysis (±1% intervals)
                SOC_targets = [45, 50, 55, 60, 65, 70]; % 6 target SOC values
                R1s_target = NaN(size(SOC_targets));
                R2_target = NaN(size(SOC_targets));
                n_points_target = zeros(size(SOC_targets));
                
                fprintf('Creating %d SOC targets (±1%% intervals)\n', length(SOC_targets));
                
                for target_idx = 1:length(SOC_targets)
                    soc_target = SOC_targets(target_idx);
                    soc_low = soc_target - 1;  % ±1% range
                    soc_high = soc_target + 1;
                    
                    % Find data points in this SOC target range
                    target_idx_data = (soc_valid >= soc_low) & (soc_valid < soc_high);
                    dI_target = dI(target_idx_data);
                    dV_target = dV(target_idx_data);
                    
                    n_points_target(target_idx) = length(dI_target);
                    
                    % Calculate R1s for this SOC target if enough data points
                    if length(dI_target) >= 5
                        mdl_target = fitlm(dI_target, dV_target, 'Intercept', false);
                        R1s_target(target_idx) = mdl_target.Coefficients.Estimate(1);
                        R2_target(target_idx) = mdl_target.Rsquared.Ordinary;
                    end
                end
                
                % Store SOC target results
                R1s_Results.(dayKey).(rackName).SOC_targets = SOC_targets;
                R1s_Results.(dayKey).(rackName).R1s_target = R1s_target;
                R1s_Results.(dayKey).(rackName).R2_target = R2_target;
                R1s_Results.(dayKey).(rackName).n_points_target = n_points_target;
                    
            end % End rack loop
        end % End mat file loop
        
        % Save monthly table and create visualizations (same as OldData)
        if ~isempty(monthly_table)
            tableFileName = sprintf('Monthly_Table_%s.mat', month_name);
            save(fullfile(saveDir, tableFileName), 'monthly_table');
            fprintf('Saved monthly table: %s\n', tableFileName);
            
            % Create subplot visualization for all days in the month
            nDays = height(monthly_table);
            if nDays > 0
                % Calculate subplot layout
                nCols = min(6, nDays); % Max 6 columns
                nRows = ceil(nDays / nCols);
                
                figure('Position', [100, 100, 1200, 800]);
                sgtitle(sprintf('%s - All Days R1s Analysis', month_name), 'FontSize', 16);
                
                for day = 1:nDays
                    subplot(nRows, nCols, day);
                    
                    % Get data for current day
                    dI_day = monthly_table.dI{day};
                    dV_day = monthly_table.dV{day};
                    R1s_day = monthly_table.R1s(day);
                    R2_day = monthly_table.R_squared(day);
                    date_day = monthly_table.Date{day};
                    
                    % Plot
                    scatter(dI_day, dV_day, 10, 'b', 'filled');
                    hold on;
                    
                    % Calculate fitted line across full axis range
                    x_range = xlim;
                    dV_fit_full = R1s_day * x_range;
                    plot(x_range, dV_fit_full, 'r-', 'LineWidth', 1.5);
                    
                    % Set symmetric axis limits around origin based on data range
                    xlim_current = xlim;
                    ylim_current = ylim;
                    x_max = max(abs(xlim_current));
                    y_max = max(abs(ylim_current));
                    xlim([-x_max, x_max]);
                    ylim([-y_max, y_max]);
                    
                    % Draw axes through origin
                    plot([0 0], ylim, 'k-', 'LineWidth', 0.5);
                    plot(xlim, [0 0], 'k-', 'LineWidth', 0.5);
                    
                    title(sprintf('%s\nR1s=%.3fmΩ, R²=%.3f', date_day(5:end), R1s_day*1000, R2_day), 'FontSize', 8);
                    xlabel('ΔI (A)', 'FontSize', 8);
                    ylabel('ΔV (V)', 'FontSize', 8);
                    grid on;
                end
                
                % Save subplot figure
                subplotFileName = sprintf('Monthly_Subplot_%s.fig', month_name);
                saveas(gcf, fullfile(saveDir, subplotFileName));
                % close(gcf);
                fprintf('Saved subplot visualization: %s\n', subplotFileName);
                
                % Create figure 2: Daily R1s trend
                figure(2);
                clf;
                
                % Extract day numbers and R1s values
                day_numbers = [];
                R1s_values = [];
                valid_days = [];
                
                for day = 1:nDays
                    date_str = monthly_table.Date{day};
                    day_num = str2double(date_str(end-1:end)); % Extract day number
                    if ~isnan(monthly_table.R1s(day)) % Only valid R1s values
                        day_numbers = [day_numbers; day_num];
                        R1s_values = [R1s_values; monthly_table.R1s(day)];
                        valid_days = [valid_days; day];
                    end
                end
                
                % Plot R1s trend
                plot(day_numbers, R1s_values*1000, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
                xlabel('Day of Month', 'FontSize', 12);
                ylabel('R1s (mΩ)', 'FontSize', 12);
                title(sprintf('%s - Daily R1s Trend', month_name), 'FontSize', 14);
                grid on;
                
                % Add data point labels
                for i = 1:length(day_numbers)
                    text(day_numbers(i), R1s_values(i)*1000, sprintf('%.3f', R1s_values(i)*1000), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
                end
                
                % Save R1s trend figure
                trendFileName = sprintf('Daily_R1s_Trend_%s.fig', month_name);
                saveas(gcf, fullfile(saveDir, trendFileName));
                % close(gcf);
                fprintf('Saved R1s trend visualization: %s\n', trendFileName);
                
                % Create SOC target figures: Each SOC target gets its own figure
                SOC_targets = [45, 50, 55, 60, 65, 70]; % 6 target SOC values
                
                for target_idx = 1:length(SOC_targets)
                    soc_target = SOC_targets(target_idx);
                    
                    % Create figure for this SOC target
                    figure('Position', [100, 100, 1400, 800]);
                    
                    % Calculate subplot layout for days
                    nCols = min(6, nDays); % Max 6 columns
                    nRows = ceil(nDays / nCols);
                    
                    sgtitle(sprintf('%s - SOC %d%%-%d%% Target Analysis', month_name, soc_target-1, soc_target+1), 'FontSize', 16);
                    
                    for day = 1:nDays
                        subplot(nRows, nCols, day);
                        
                        date_str = monthly_table.Date{day};
                        dayKey = date_str;
                        
                        if isfield(R1s_Results, dayKey) && isfield(R1s_Results.(dayKey), rackNames_all{1})
                            rackData = R1s_Results.(dayKey).(rackNames_all{1});
                            
                            if isfield(rackData, 'SOC_targets') && isfield(rackData, 'SOC') && isfield(rackData, 'dI') && isfield(rackData, 'dV')
                                % Get data for this day
                                soc_data = rackData.SOC;
                                dI_data = rackData.dI;
                                dV_data = rackData.dV;
                                
                                % Find data points in this SOC target range
                                soc_low = soc_target - 1;
                                soc_high = soc_target + 1;
                                target_idx_data = (soc_data >= soc_low) & (soc_data < soc_high);
                                
                                dI_target = dI_data(target_idx_data);
                                dV_target = dV_data(target_idx_data);
                                
                                % Plot if we have data
                                if length(dI_target) >= 5
                                    % Scatter plot
                                    scatter(dI_target, dV_target, 10, 'b', 'filled');
                                    hold on;
                                    
                                    % Fitting line
                                    mdl = fitlm(dI_target, dV_target, 'Intercept', false);
                                    R1s_value = mdl.Coefficients.Estimate(1);
                                    
                                    % Plot fitting line across full axis range
                                    x_range = xlim;
                                    dV_fit_full = R1s_value * x_range;
                                    plot(x_range, dV_fit_full, 'r-', 'LineWidth', 1.5);
                                    
                                    R2_value = mdl.Rsquared.Ordinary;
                                    
                                    title(sprintf('%s: SOC %d%%\nR1s=%.3fmΩ, R²=%.3f, n=%d', ...
                                        date_str(5:end), soc_target, R1s_value*1000, R2_value, length(dI_target)), 'FontSize', 8);
                                    
                                    % Set symmetric axis limits around origin based on data range
                                    xlim_current = xlim;
                                    ylim_current = ylim;
                                    x_max = max(abs(xlim_current));
                                    y_max = max(abs(ylim_current));
                                    xlim([-x_max, x_max]);
                                    ylim([-y_max, y_max]);
                                    
                                    % Draw axes through origin
                                    plot([0 0], ylim, 'k-', 'LineWidth', 0.5);
                                    plot(xlim, [0 0], 'k-', 'LineWidth', 0.5);
                                    
                                    xlabel('ΔI (A)', 'FontSize', 8);
                                    ylabel('ΔV (V)', 'FontSize', 8);
                                    grid on;
                                else
                                    title(sprintf('%s: SOC %d%%\nNo data', date_str(5:end), soc_target), 'FontSize', 8);
                                    xlabel('ΔI (A)', 'FontSize', 8);
                                    ylabel('ΔV (V)', 'FontSize', 8);
                                end
                            else
                                title(sprintf('%s: SOC %d%%\nNo data', date_str(5:end), soc_target), 'FontSize', 8);
                                xlabel('ΔI (A)', 'FontSize', 8);
                                ylabel('ΔV (V)', 'FontSize', 8);
                            end
                        else
                            title(sprintf('Day %d: SOC %d%%\nNo data', day, soc_target), 'FontSize', 8);
                            xlabel('ΔI (A)', 'FontSize', 8);
                            ylabel('ΔV (V)', 'FontSize', 8);
                        end
                    end
                    
                    % Save SOC target figure
                    socFileName = sprintf('SOC_%d_Target_%s.fig', soc_target, month_name);
                    saveas(gcf, fullfile(saveDir, socFileName));
                    % close(gcf);
                    fprintf('Saved SOC %d%% target analysis: %s\n', soc_target, socFileName);
                end
            end
        end
        
    end % End month loop
    
    % Save results for current year
    if ~isempty(fieldnames(R1s_Results))
        saveFileName = sprintf('R1s_Results_%s.mat', year);
        save(fullfile(saveDir, saveFileName), 'R1s_Results');
        fprintf('Saved results to: %s\n', saveFileName);
        
        % Create summary statistics
        allR1s = [];
        allR2 = [];
        allN = [];
        
        dayNames = fieldnames(R1s_Results);
        for d = 1:length(dayNames)
            if isfield(R1s_Results.(dayNames{d}), rackNames_all{1})
                allR1s = [allR1s; R1s_Results.(dayNames{d}).(rackNames_all{1}).R1s];
                allR2 = [allR2; R1s_Results.(dayNames{d}).(rackNames_all{1}).R_squared];
                allN = [allN; R1s_Results.(dayNames{d}).(rackNames_all{1}).n_points];
            end
        end
        
        if ~isempty(allR1s)
            fprintf('\n=== %s Summary ===\n', year);
            fprintf('Total days processed: %d\n', length(allR1s));
            fprintf('Mean R1s: %.4f mΩ\n', mean(allR1s)*1000);
            fprintf('Std R1s: %.4f mΩ\n', std(allR1s)*1000);
            fprintf('Mean R²: %.3f\n', mean(allR2));
            fprintf('Mean data points per day: %.1f\n', mean(allN));
        end
    else
        fprintf('No valid data found for %s\n', year);
    end
    
end % End year loop

fprintf('\nR1s calculation completed!\n');