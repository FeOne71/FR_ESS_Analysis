%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Qmax Calculation - Find OCV1 and OCV2 points
% Pattern: Rest(1hr) -> Charge/Discharge -> Rest(1hr)
% Qmax = integral(i*dt) / (SOC2 - SOC1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory and file paths
dataDir  = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\raw2mat_ver04';
ocvFuncPath = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\Summary\OCV_interpolation_functions.mat';
saveDir  = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\Dataset_Overview';

if ~exist(saveDir, 'dir')
    mkdir(saveDir); 
end

%% Load OCV interpolation function
load(ocvFuncPath);
fprintf('Loaded OCV interpolation function\n');

%% Initialize arrays to store Qmax data
qmax_data = struct();
qmax_data.ocv1_soc = [];
qmax_data.ocv1_voltage = [];
qmax_data.ocv2_soc = [];
qmax_data.ocv2_voltage = [];
qmax_data.charge_type = []; % 'charge' or 'discharge'
qmax_data.integral_idt = [];
qmax_data.soc_diff = [];
qmax_data.qmax = [];
qmax_data.year = [];
qmax_data.time = [];

yearList = {'2021', '2022', '2023'};

fprintf('=====Searching for Qmax Patterns=====\n');

for y = 1:length(yearList)
    year = yearList{y};
    yearPath = fullfile(dataDir, year);
    monthDirs = dir(fullfile(yearPath, '20*'));

    for m = 1:length(monthDirs)
        monthPath = fullfile(yearPath, monthDirs(m).name);
        matFiles = dir(fullfile(monthPath, 'Raw_*.mat'));

        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            [~, name, ~] = fileparts(matFiles(f).name);
            dateStr = extractAfter(name, 'Raw_');
            fieldName = ['date_', dateStr];

            fprintf('Processing: %s\n', matFiles(f).name);
            
            load(matFilePath); 
            t          = Raw.Time;
            I          = Raw.DCCurrent_A;
            V          = Raw.Total_AverageCVSum_V;
            soc        = Raw.Total_AverageSOC;
            bsc_charge = Raw.Charge;  % string array
            
            % Convert datetime to seconds if needed
            if isdatetime(t)
                t_seconds = seconds(t - t(1));  % Convert to seconds from start
            else
                t_seconds = t;  % Already in seconds
            end
            
            % Find Qmax patterns: Rest -> Charge/Discharge -> Rest
            qmax_patterns = findQmaxPatterns(t_seconds, bsc_charge, I, V, soc, OCV_integrated_RPT0cyc.OCV_SOC_func);
            
            % Process each pattern
            for i = 1:length(qmax_patterns)
                pattern = qmax_patterns{i};
                if ~isempty(pattern)
                    % Extract OCV1 and OCV2 values
                    ocv1_idx = pattern.ocv1_idx;
                    ocv2_idx = pattern.ocv2_idx;
                    charge_start = pattern.charge_start;
                    charge_end = pattern.charge_end;
                    charge_type = pattern.charge_type;
                    
                    % Get OCV values using interpolation function
                    ocv1_soc = soc(ocv1_idx);
                    ocv2_soc = soc(ocv2_idx);
                    
                    % Use original SOC values directly with interpolation function
                    ocv1_voltage = OCV_integrated_RPT0cyc.OCV_SOC_func(ocv1_soc);
                    ocv2_voltage = OCV_integrated_RPT0cyc.OCV_SOC_func(ocv2_soc);
                    
                    % Debug: Show SOC and voltage values
                    fprintf('    OCV1: SOC %.3f, Voltage = %.3f V\n', ocv1_soc, ocv1_voltage);
                    fprintf('    OCV2: SOC %.3f, Voltage = %.3f V\n', ocv2_soc, ocv2_voltage);
                    
                    % Calculate integral of current over time
                    charge_indices = charge_start:charge_end;
                    integral_idt = trapz(t_seconds(charge_indices), I(charge_indices));
                    
                    % Calculate SOC difference
                    soc_diff = ocv2_soc - ocv1_soc;
                    
                    % Calculate Qmax
                    if abs(soc_diff) > 0.1  % Avoid division by very small numbers
                        qmax = integral_idt / soc_diff;
                        
                        % Store data
                        qmax_data.ocv1_soc = [qmax_data.ocv1_soc; ocv1_soc];
                        qmax_data.ocv1_voltage = [qmax_data.ocv1_voltage; ocv1_voltage];
                        qmax_data.ocv2_soc = [qmax_data.ocv2_soc; ocv2_soc];
                        qmax_data.ocv2_voltage = [qmax_data.ocv2_voltage; ocv2_voltage];
                        qmax_data.charge_type = [qmax_data.charge_type; {charge_type}];
                        qmax_data.integral_idt = [qmax_data.integral_idt; integral_idt];
                        qmax_data.soc_diff = [qmax_data.soc_diff; soc_diff];
                        qmax_data.qmax = [qmax_data.qmax; qmax];
                        qmax_data.year = [qmax_data.year; y];
                        qmax_data.time = [qmax_data.time; t_seconds(ocv1_idx)];
                        
                        fprintf('  Found pattern: %s, Qmax = %.2f Ah, SOC diff = %.3f\n', ...
                                charge_type, qmax, soc_diff);
                    end
                end
            end
        end
    end
end

%% Process and save Qmax data
fprintf('\n=====Processing Qmax Data=====\n');
fprintf('Total Qmax patterns found: %d\n', length(qmax_data.qmax));

if length(qmax_data.qmax) > 0
    % Remove any invalid data points
    valid_idx = ~isnan(qmax_data.qmax) & ~isinf(qmax_data.qmax);
    qmax_data.ocv1_soc = qmax_data.ocv1_soc(valid_idx);
    qmax_data.ocv1_voltage = qmax_data.ocv1_voltage(valid_idx);
    qmax_data.ocv2_soc = qmax_data.ocv2_soc(valid_idx);
    qmax_data.ocv2_voltage = qmax_data.ocv2_voltage(valid_idx);
    qmax_data.charge_type = qmax_data.charge_type(valid_idx);
    qmax_data.integral_idt = qmax_data.integral_idt(valid_idx);
    qmax_data.soc_diff = qmax_data.soc_diff(valid_idx);
    qmax_data.qmax = qmax_data.qmax(valid_idx);
    qmax_data.year = qmax_data.year(valid_idx);
    qmax_data.time = qmax_data.time(valid_idx);
    
    % Save Qmax data
    save(fullfile(saveDir, 'Qmax_Data.mat'), 'qmax_data');
    
    %% Plot Qmax results
    figure('Position', [100, 100, 1200, 800]);
    
    % Define colors for each year
    colors = [32 133 78; 7 115 194; 205 83 76] / 255;
    year_names = {'2021', '2022', '2023'};
    
    % Qmax vs SOC difference
    subplot(2,2,1);
    hold on;
    for y = 1:length(yearList)
        year_idx = qmax_data.year == y;
        if sum(year_idx) > 0
            scatter(qmax_data.soc_diff(year_idx), qmax_data.qmax(year_idx), 50, colors(y,:), 'filled', 'DisplayName', year_names{y});
        end
    end
    xlabel('SOC Difference [-]');
    ylabel('Qmax [Ah]');
    title('Qmax vs SOC Difference');
    legend('show', 'Location', 'best');
    grid on;
    hold off;
    
    % Qmax distribution by year
    subplot(2,2,2);
    hold on;
    for y = 1:length(yearList)
        year_idx = qmax_data.year == y;
        if sum(year_idx) > 0
            histogram(qmax_data.qmax(year_idx), 15, 'FaceColor', colors(y,:), 'FaceAlpha', 0.7, 'DisplayName', year_names{y});
        end
    end
    xlabel('Qmax [Ah]');
    ylabel('Count');
    title('Qmax Distribution by Year');
    legend('show', 'Location', 'best');
    grid on;
    hold off;
    
    % OCV1 vs OCV2
    subplot(2,2,3);
    hold on;
    for y = 1:length(yearList)
        year_idx = qmax_data.year == y;
        if sum(year_idx) > 0
            scatter(qmax_data.ocv1_soc(year_idx), qmax_data.ocv2_soc(year_idx), 50, colors(y,:), 'filled', 'DisplayName', year_names{y});
        end
    end
    xlabel('OCV1 SOC [-]');
    ylabel('OCV2 SOC [-]');
    title('OCV1 vs OCV2 SOC');
    legend('show', 'Location', 'best');
    grid on;
    hold off;
    
    % Charge type distribution
    subplot(2,2,4);
    charge_types = unique(qmax_data.charge_type);
    charge_counts = zeros(length(charge_types), 1);
    for i = 1:length(charge_types)
        charge_counts(i) = sum(strcmp(qmax_data.charge_type, charge_types{i}));
    end
    bar(charge_counts);
    set(gca, 'XTickLabel', charge_types);
    xlabel('Charge Type');
    ylabel('Count');
    title('Charge Type Distribution');
    grid on;
    
    % Save figure
    saveas(gcf, fullfile(saveDir, 'Qmax_Results.fig'));
    
    fprintf('Qmax data saved to: %s\n', saveDir);
    fprintf('Total valid Qmax points: %d\n', length(qmax_data.qmax));
    
    % Show statistics
    fprintf('\nQmax Statistics:\n');
    fprintf('Mean Qmax: %.2f Ah\n', mean(qmax_data.qmax));
    fprintf('Std Qmax: %.2f Ah\n', std(qmax_data.qmax));
    fprintf('Min Qmax: %.2f Ah\n', min(qmax_data.qmax));
    fprintf('Max Qmax: %.2f Ah\n', max(qmax_data.qmax));
    
else
    fprintf('No valid Qmax patterns found!\n');
end

%% Function to find Qmax patterns
function qmax_patterns = findQmaxPatterns(t, bsc_charge, I, V, soc, ocv_func)
    qmax_patterns = {};
    
    % Find all periods
    idle_idx = strcmp(bsc_charge, 'Idle');
    charge_idx = strcmp(bsc_charge, 'Charge');
    discharge_idx = strcmp(bsc_charge, 'Discharge');
    
    % Get start/end indices for each type
    idle_start = find(diff([0; idle_idx]) == 1);
    idle_end = find(diff([idle_idx; 0]) == -1);
    charge_start = find(diff([0; charge_idx]) == 1);
    charge_end = find(diff([charge_idx; 0]) == -1);
    discharge_start = find(diff([0; discharge_idx]) == 1);
    discharge_end = find(diff([discharge_idx; 0]) == -1);
    
    fprintf('Found: %d idle, %d charge, %d discharge periods\n', ...
            length(idle_start), length(charge_start), length(discharge_start));
    
    % Combine all periods and sort by start time
    all_periods = {};
    period_count = 0;
    
    % Add idle periods
    for i = 1:length(idle_start)
        duration = t(idle_end(i)) - t(idle_start(i));
        if duration >= 600  % Only periods >= 600s
            period_count = period_count + 1;
            all_periods{period_count} = [idle_start(i), idle_end(i), 0, duration];
        end
    end
    
    % Add charge periods (no time limit)
    for i = 1:length(charge_start)
        duration = t(charge_end(i)) - t(charge_start(i));
        period_count = period_count + 1;
        all_periods{period_count} = [charge_start(i), charge_end(i), 1, duration];
    end
    
    % Add discharge periods (no time limit)
    for i = 1:length(discharge_start)
        duration = t(discharge_end(i)) - t(discharge_start(i));
        period_count = period_count + 1;
        all_periods{period_count} = [discharge_start(i), discharge_end(i), 2, duration];
    end
    
    % Check if we have any periods
    if isempty(all_periods)
        fprintf('  No valid periods found\n');
        return;
    end
    
    fprintf('  Total periods found: %d\n', length(all_periods));
    
    % Sort by start time
    start_times = cellfun(@(x) x(1), all_periods);
    [~, sort_idx] = sort(start_times);
    all_periods = all_periods(sort_idx);
    
    % Find patterns: Idle -> Charge/Discharge -> Idle
    if length(all_periods) < 3
        fprintf('  Not enough periods to form patterns (need at least 3, found %d)\n', length(all_periods));
        return;
    end
    
    fprintf('  Looking for patterns in %d periods\n', length(all_periods));
    
    for i = 1:length(all_periods)-2
        fprintf('  Checking pattern at index %d\n', i);
        if all_periods{i}(3) == 0 && all_periods{i+1}(3) > 0 && all_periods{i+2}(3) == 0
            % Found pattern!
            rest1_start = all_periods{i}(1);
            rest1_end = all_periods{i}(2);
            charge_start_idx = all_periods{i+1}(1);
            charge_end_idx = all_periods{i+1}(2);
            rest2_start = all_periods{i+2}(1);
            rest2_end = all_periods{i+2}(2);
            
            % Check if charge/discharge period is between the two rests
            if charge_start_idx > rest1_end && rest2_start > charge_end_idx
                pattern = struct();
                pattern.ocv1_idx = rest1_end;
                pattern.ocv2_idx = rest2_end;
                pattern.charge_start = charge_start_idx;
                pattern.charge_end = charge_end_idx;
                
                if all_periods{i+1}(3) == 1
                    pattern.charge_type = 'charge';
                    fprintf('  Found: Rest(%.1fmin) -> Charge(%.1fmin) -> Rest(%.1fmin)\n', ...
                            all_periods{i}(4)/60, all_periods{i+1}(4)/60, all_periods{i+2}(4)/60);
                else
                    pattern.charge_type = 'discharge';
                    fprintf('  Found: Rest(%.1fmin) -> Discharge(%.1fmin) -> Rest(%.1fmin)\n', ...
                            all_periods{i}(4)/60, all_periods{i+1}(4)/60, all_periods{i+2}(4)/60);
                end
                
                qmax_patterns{end+1} = pattern;
            end
        end
    end
    
    fprintf('  Valid Qmax patterns: %d\n', length(qmax_patterns));
end 