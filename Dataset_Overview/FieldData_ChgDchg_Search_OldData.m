%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Charge/Discharge Event Search with Rest Periods
% 2021-2023 June data
% Find SOC variation >= 10% with 10+ min rest periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Parameters
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old';
yearList = {'2021', '2022', '2023'}; 
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Dataset_Overview\ChgDchg_Search_OldData';

C_nom_cell = 64;
Np = 2;
soc_threshold = 10;  % SOC variation threshold (%)
rest_duration_min = 3;  % Minimum rest period (minutes)
idle_threshold = C_nom_cell * 0.05;  % Idle current threshold (A)

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

fprintf('Searching for charge/discharge events with SOC variation >= %d%% and rest periods >= %d minutes...\n', soc_threshold, rest_duration_min);

%% Main processing
results = struct();
result_count = 0;

for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    year_path = fullfile(dataDir, year);
    
    if ~exist(year_path, 'dir')
        continue;
    end
    
    % Get months up to June only
    month_folders = dir(fullfile(year_path, '202*'));
    month_folders = month_folders([month_folders.isdir]);
    
    valid_months = {};
    for i = 1:length(month_folders)
        month_name = month_folders(i).name;
        if length(month_name) >= 6
            month_num = str2double(month_name(5:6));
            if month_num <= 6
                valid_months{end+1} = month_name;
            end
        end
    end
    
    for month_idx = 1:length(valid_months)
        month_folder = valid_months{month_idx};
        month_path = fullfile(year_path, month_folder);
        
        mat_files = dir(fullfile(month_path, 'Raw_*.mat'));
        
        for file_idx = 1:length(mat_files)
            file_name = mat_files(file_idx).name;
            file_path = fullfile(month_path, file_name);
            
            try
                fprintf('Processing: %s\n', file_path);
                data = load(file_path);
                
                % Extract data
                event_data = data.Raw;
                rack_data = event_data.Rack01;
                
                soc_data = rack_data.SOCPct;
                current_data = rack_data.DCCurrent_A / Np;  % Convert to cell level
                time_strings = rack_data.Time;
                voltage_data = rack_data.AverageCV_V;
                power_data = rack_data.DCPower_kW;
                temp_data = rack_data.AverageMT_degC;
                
                % Filter valid data
                valid_mask = ~isnan(soc_data) & ~isinf(soc_data) & ~isnan(current_data) & ~isinf(current_data);
                valid_soc = soc_data(valid_mask);
                valid_current = current_data(valid_mask);
                valid_time = datetime(time_strings(valid_mask), 'Format', 'yyyy-MM-dd HH:mm:ss');
                valid_voltage = voltage_data(valid_mask);
                valid_power = power_data(valid_mask);
                valid_temp = temp_data(valid_mask);
                
                if length(valid_soc) > 10
                    % Find charge/discharge events
                    events = find_events(valid_current, valid_soc, valid_time, soc_threshold, rest_duration_min, idle_threshold);
                    
                    % Save results
                    for seg_idx = 1:length(events)
                        event = events(seg_idx);
                        result_count = result_count + 1;
                        
                        results(result_count).file_name = file_name;
                        results(result_count).year = year;
                        results(result_count).month = month_folder;
                        results(result_count).event_type = event.type;
                        results(result_count).soc1 = event.soc1;
                        results(result_count).soc2 = event.soc2;
                        results(result_count).delta_soc = event.delta_soc;
                        results(result_count).rest1_duration_min = event.rest1_duration_min;
                        results(result_count).rest2_duration_min = event.rest2_duration_min;
                        results(result_count).t1 = event.t1;
                        results(result_count).t2 = event.t2;
                        
                        % Store data for plotting
                        start_idx = event.start_idx;
                        end_idx = event.end_idx;
                        results(result_count).soc_data = valid_soc(start_idx:end_idx);
                        results(result_count).current_data = valid_current(start_idx:end_idx);
                        results(result_count).voltage_data = valid_voltage(start_idx:end_idx);
                        results(result_count).power_data = valid_power(start_idx:end_idx);
                        results(result_count).temp_data = valid_temp(start_idx:end_idx);
                        results(result_count).time_data = valid_time(start_idx:end_idx);
                        
                        % Store full day data for context plotting
                        results(result_count).full_soc_data = valid_soc;
                        results(result_count).full_current_data = valid_current;
                        results(result_count).full_voltage_data = valid_voltage;
                        results(result_count).full_power_data = valid_power;
                        results(result_count).full_temp_data = valid_temp;
                        results(result_count).full_time_data = valid_time;
                        results(result_count).segment_start_idx = start_idx;
                        results(result_count).segment_end_idx = end_idx;
                        
                        fprintf('  -> %s event: SOC1=%.2f%%, SOC2=%.2f%%, ΔSOC=%.2f%%, Rest1=%.1fmin, Rest2=%.1fmin\n', ...
                            event.type, event.soc1, event.soc2, event.delta_soc, event.rest1_duration_min, event.rest2_duration_min);
                    end
                end
                
            catch ME
                fprintf('  -> Error: %s\n', ME.message);
            end
        end
                    end
                end
                
%% Save and display results
fprintf('\n=== Analysis Complete ===\n');
fprintf('Found %d events\n', result_count);

if result_count > 0
    result_table = struct2table(results, 'AsArray', true);
    result_file = fullfile(saveDir, sprintf('Events_SOC%d%%_Rest%dmin_Results.mat', soc_threshold, rest_duration_min));
    save(result_file, 'results', 'result_table', 'soc_threshold', 'rest_duration_min', 'idle_threshold', 'C_nom_cell', 'Np');
    
    fprintf('Results saved to: %s\n', result_file);
    
    % Display results
    fprintf('\n=== All Results ===\n');
    [~, sorted_idx] = sort([results.delta_soc], 'descend');
    for i = 1:length(sorted_idx)
        idx = sorted_idx(i);
        fprintf('%d. %s: %s - SOC1=%.2f%%, SOC2=%.2f%%, ΔSOC=%.2f%%, Rest1=%.1fmin, Rest2=%.1fmin\n', ...
            i, results(idx).file_name, results(idx).event_type, results(idx).soc1, ...
            results(idx).soc2, results(idx).delta_soc, results(idx).rest1_duration_min, results(idx).rest2_duration_min);
    end
    
    %% Plot results
    fprintf('\n=== Plotting Events ===\n');
    
    for plot_idx = 1:length(sorted_idx)
        idx = sorted_idx(plot_idx);
        
        figure('Position', [100, 100, 1200, 800]);
        
        % Get full day data and segment data
        full_time = results(idx).full_time_data;
        full_soc = results(idx).full_soc_data;
        full_current = results(idx).full_current_data;
        full_voltage = results(idx).full_voltage_data;
        full_power = results(idx).full_power_data;
        
        segment_time = results(idx).time_data;
        segment_soc = results(idx).soc_data;
        segment_current = results(idx).current_data;
        segment_voltage = results(idx).voltage_data;
        segment_power = results(idx).power_data;
        
        % Subplot 1: SOC and Current
        subplot(3, 1, 1);
        yyaxis left;
        plot(full_time, full_soc, 'b-', 'LineWidth', 0.5, 'Color', [0.7 0.7 1]);
        hold on;
        plot(segment_time, segment_soc, 'b-', 'LineWidth', 2);
        ylabel('SOC (%)');
        ylim([0, 100]);
        grid on;
        
        yyaxis right;
        plot(full_time, full_current, 'r-', 'LineWidth', 0.5, 'Color', [1 0.7 0.7]);
        hold on;
        plot(segment_time, segment_current, 'r-', 'LineWidth', 2);
        ylabel('Current (A)');
        title(sprintf('%s: %s Event - SOC & Current', results(idx).file_name, results(idx).event_type));
        xlabel('Time');
        legend('SOC (full)', 'SOC (segment)', 'Current (full)', 'Current (segment)', 'Location', 'best');
        
        % Subplot 2: Current and Voltage
        subplot(3, 1, 2);
        yyaxis left;
        plot(full_time, full_current, 'r-', 'LineWidth', 0.5, 'Color', [1 0.7 0.7]);
        hold on;
        plot(segment_time, segment_current, 'r-', 'LineWidth', 2);
        ylabel('Current (A)');
        grid on;
        
        yyaxis right;
        plot(full_time, full_voltage, 'g-', 'LineWidth', 0.5, 'Color', [0.7 1 0.7]);
        hold on;
        plot(segment_time, segment_voltage, 'g-', 'LineWidth', 2);
        ylabel('Voltage (V)');
        title(sprintf('Current & Voltage - SOC1=%.2f%%, SOC2=%.2f%%, ΔSOC=%.2f%%', ...
            results(idx).soc1, results(idx).soc2, results(idx).delta_soc));
        xlabel('Time');
        legend('Current (full)', 'Current (segment)', 'Voltage (full)', 'Voltage (segment)', 'Location', 'best');
        
        % Subplot 3: Power and Voltage
        subplot(3, 1, 3);
        yyaxis left;
        plot(full_time, full_power, 'm-', 'LineWidth', 0.5, 'Color', [1 0.7 1]);
        hold on;
        plot(segment_time, segment_power, 'm-', 'LineWidth', 2);
        ylabel('Power (kW)');
        grid on;
        
        yyaxis right;
        plot(full_time, full_voltage, 'g-', 'LineWidth', 0.5, 'Color', [0.7 1 0.7]);
        hold on;
        plot(segment_time, segment_voltage, 'g-', 'LineWidth', 2);
        ylabel('Voltage (V)');
        
        xlabel('Time');
        title(sprintf('Power & Voltage - Rest1=%.1fmin, Rest2=%.1fmin', ...
            results(idx).rest1_duration_min, results(idx).rest2_duration_min));
        legend('Power (full)', 'Power (segment)', 'Voltage (full)', 'Voltage (segment)', 'Location', 'best');
        
        % Save plot
        plot_file = fullfile(saveDir, sprintf('Event_%d_%s_%s.fig', ...
            plot_idx, results(idx).file_name(1:end-4), results(idx).event_type));
        savefig(plot_file);
        fprintf('Plot %d saved: %s\n', plot_idx, plot_file);
    end
    
    fprintf('All plots saved to: %s\n', saveDir);
    
else
    fprintf('No events found.\n');
end

%% Simple function to find events
function events = find_events(current_data, soc_data, time_data, soc_threshold, rest_duration_min, idle_threshold)
    events = struct('type', {}, 'start_idx', {}, 'end_idx', {}, 'soc1', {}, 'soc2', {}, 'delta_soc', {}, ...
                   'rest1_duration_min', {}, 'rest2_duration_min', {}, 't1', {}, 't2', {});
    
    n = length(current_data);
    if n < 10
        return;
    end
    
    % Find all rest periods (|I| < idle_threshold)
    rest_mask = abs(current_data) < idle_threshold;
    
    % Find rest period boundaries
    rest_starts = [];
    rest_ends = [];
    
    in_rest = false;
    for i = 1:n
        if rest_mask(i) && ~in_rest
            rest_starts(end+1) = i;
            in_rest = true;
        elseif ~rest_mask(i) && in_rest
            rest_ends(end+1) = i-1;
            in_rest = false;
        end
    end
    
    if in_rest
        rest_ends(end+1) = n;
    end
    
    % Check each rest period for potential events
    for i = 1:length(rest_starts)
        rest1_start = rest_starts(i);
        rest1_end = rest_ends(i);
        
        % Check if rest1 is long enough
        rest1_duration = time_data(rest1_end) - time_data(rest1_start);
        rest1_duration_min = minutes(rest1_duration);
        
        if rest1_duration_min >= rest_duration_min
            % Look for charge/discharge after rest1
            activity_start = rest1_end + 1;
            activity_end = find_activity_end(current_data, activity_start, idle_threshold);
            
            if activity_end > activity_start
                % Look for rest2 after activity
                rest2_start = find_next_rest(current_data, activity_end, idle_threshold);
                
                if rest2_start > activity_end
                    rest2_end = find_rest_end(current_data, rest2_start, idle_threshold);
                    
                    if rest2_end > rest2_start
                        % Check if rest2 is long enough
                        rest2_duration = time_data(rest2_end) - time_data(rest2_start);
                        rest2_duration_min = minutes(rest2_duration);
                        
                        if rest2_duration_min >= rest_duration_min
                            % Get SOC values
                            soc1 = soc_data(rest1_end);
                            soc2 = soc_data(rest2_end);
                            delta_soc = abs(soc2 - soc1);
                            
                            % Check SOC variation
                            if delta_soc >= soc_threshold
                                % Determine event type
                                if soc2 > soc1
                                    event_type = 'charging';
                                else
                                    event_type = 'discharging';
                                end
                                
                                % Create event
                                event = struct();
                                event.type = event_type;
                                event.start_idx = rest1_start;
                                event.end_idx = rest2_end;
                                event.soc1 = soc1;
                                event.soc2 = soc2;
                                event.delta_soc = delta_soc;
                                event.rest1_duration_min = rest1_duration_min;
                                event.rest2_duration_min = rest2_duration_min;
                                event.t1 = time_data(rest1_end);
                                event.t2 = time_data(rest2_end);
                                
                                events = [events; event];
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Helper functions
function activity_end = find_activity_end(current_data, start_idx, idle_threshold)
    n = length(current_data);
    activity_end = start_idx;
    
    for i = start_idx:n
        if abs(current_data(i)) >= idle_threshold
            activity_end = i;
        else
            break;
        end
    end
end

function rest_start = find_next_rest(current_data, start_idx, idle_threshold)
    n = length(current_data);
    rest_start = n + 1;
    
    for i = start_idx:n
        if abs(current_data(i)) < idle_threshold
            rest_start = i;
            break;
        end
    end
end

function rest_end = find_rest_end(current_data, start_idx, idle_threshold)
    n = length(current_data);
    rest_end = start_idx;
    
    for i = start_idx:n
        if abs(current_data(i)) < idle_threshold
            rest_end = i;
        else
            break;
        end
    end
end