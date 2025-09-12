%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New Kimje Data Rated Measurement Profile Search
% 202309- 
% 202401-
% 202501-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New';
yearList = {'2024'}; 
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\Dataset_Overview\FieldData_RPT_Check');

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameters
C_nom_rack = 128;
soc_variation_threshold = 8;  % SOC variation range threshold for charge/discharge cycle (%)
current_variation_threshold = 10;  % Current variation threshold for charge/discharge cycle (A)

%% Find continuous charge/discharge cycles
fprintf('Searching for continuous charge/discharge cycles with SOC variation >= %d%% and current variation >= %dA...\n', soc_variation_threshold, current_variation_threshold);

% Structure to store results
results = struct();
result_count = 0;

% Iterate through years
for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    year_path = fullfile(dataDir, year);
    
    if ~exist(year_path, 'dir')
        fprintf('Year folder does not exist: %s\n', year_path);
        continue;
    end
    
    % Find month folders
    month_folders = dir(fullfile(year_path, '202*'));
    month_folders = month_folders([month_folders.isdir]);
    
    for month_idx = 1:length(month_folders)
        month_folder = month_folders(month_idx).name;
        month_path = fullfile(year_path, month_folder);
        
        % Find daily mat files
        mat_files = dir(fullfile(month_path, 'Raw_*.mat'));
        
        for file_idx = 1:length(mat_files)
            file_name = mat_files(file_idx).name;
            file_path = fullfile(month_path, file_name);
            
            try
                % Load mat file
                fprintf('Processing: %s\n', file_path);
                data = load(file_path);
                
                % Extract data from Raw
                event_data = data.Raw;
                
                % Get SOC_BMS and current data
                soc_data = event_data.SOC_BMS;
                current_data = event_data.DCCurrent;
                
                % Filter valid data (remove NaN, Inf)
                valid_mask = ~isnan(soc_data) & ~isinf(soc_data) & ~isnan(current_data) & ~isinf(current_data);
                valid_soc = soc_data(valid_mask);
                valid_current = current_data(valid_mask);
                valid_time = event_data.Date_Time(valid_mask);
                
                if length(valid_soc) > 10  % Need at least 10 data points
                    % Find continuous charging segments (positive current, increasing SOC)
                    charging_segments = find_continuous_segments(valid_current, valid_soc, 'charging', soc_variation_threshold, current_variation_threshold);
                    
                    % Find continuous discharging segments (negative current, decreasing SOC)
                    discharging_segments = find_continuous_segments(valid_current, valid_soc, 'discharging', soc_variation_threshold, current_variation_threshold);
                    
                    % Combine all segments
                    all_segments = [charging_segments; discharging_segments];
                    
                    % Save results for each valid segment
                    for seg_idx = 1:length(all_segments)
                        segment = all_segments(seg_idx);
                        
                        result_count = result_count + 1;
                        
                        % Save results
                        results(result_count).file_path = file_path;
                        results(result_count).file_name = file_name;
                        results(result_count).year = year;
                        results(result_count).month = month_folder;
                        results(result_count).cycle_type = segment.type;
                        results(result_count).soc_range = segment.soc_range;
                        results(result_count).soc_min = segment.soc_min;
                        results(result_count).soc_max = segment.soc_max;
                        results(result_count).current_range = segment.current_range;
                        results(result_count).current_min = segment.current_min;
                        results(result_count).current_max = segment.current_max;
                        results(result_count).data_points = segment.length;
                        results(result_count).start_idx = segment.start_idx;
                        results(result_count).end_idx = segment.end_idx;
                        
                                                 % Store segment data for plotting
                         results(result_count).soc_data = segment.soc_data;
                         results(result_count).current_data = segment.current_data;
                         
                         % Store voltage data for the segment (matching current data length)
                         segment_voltage = event_data.CVavg(valid_mask);
                         results(result_count).voltage_data = segment_voltage(segment.start_idx:segment.end_idx);
                         
                         % Store temperature data for the segment (matching current data length)
                         segment_temp = event_data.MTavg(valid_mask);
                         results(result_count).temp_data = segment_temp(segment.start_idx:segment.end_idx);
                         
                         % Store power data for the segment (matching current data length)
                         segment_power = event_data.DCPower(valid_mask);
                         results(result_count).power_data = segment_power(segment.start_idx:segment.end_idx);
                         
                         results(result_count).time_data = valid_time(segment.start_idx:segment.end_idx);
                        
                        % Store original filtered data for extended plotting
                        results(result_count).original_time_data = valid_time;
                        results(result_count).original_soc_data = valid_soc;
                        results(result_count).original_current_data = valid_current;
                        
                        % Store original filtered voltage, power, and temperature data
                        results(result_count).original_voltage_data = event_data.CVavg(valid_mask);
                        results(result_count).original_power_data = event_data.DCPower(valid_mask);
                        results(result_count).original_temp_data = event_data.MTavg(valid_mask);
                        
                        fprintf('  -> %s segment found: SOC range %.2f%% (%.2f%% ~ %.2f%%), Current range %.2fA (%.2fA ~ %.2fA), Length: %d points\n', ...
                            segment.type, segment.soc_range, segment.soc_min, segment.soc_max, segment.current_range, segment.current_min, segment.current_max, segment.length);
                    end
                end
                
            catch ME
                fprintf('  -> File load error: %s\n', ME.message);
            end
        end
    end
end

%% Summary and save results
fprintf('\n=== Analysis Complete ===\n');
fprintf('Found %d continuous charge/discharge segments (SOC >= %d%%, Current >= %dA).\n', ...
    result_count, soc_variation_threshold, current_variation_threshold);

if result_count > 0
    % Convert results to table
    result_table = struct2table(results, 'AsArray', true);
    
    % Save results
    result_file = fullfile(saveDir, sprintf('Continuous_Cycles_SOC%d%%_Current%dA_Results.mat', soc_variation_threshold, current_variation_threshold));
    save(result_file, 'results', 'result_table', 'soc_variation_threshold', 'current_variation_threshold');
    
    % Also save as CSV file
    csv_file = fullfile(saveDir, sprintf('Continuous_Cycles_SOC%d%%_Current%dA_Results.csv', soc_variation_threshold, current_variation_threshold));
    writetable(result_table, csv_file);
    
    fprintf('Results saved to:\n');
    fprintf('  - MAT file: %s\n', result_file);
    fprintf('  - CSV file: %s\n', csv_file);
    
    % Display top 10 results
    fprintf('\n=== Top 10 Results ===\n');
    [~, sorted_idx] = sort([results.soc_range], 'descend');
    for i = 1:min(10, length(sorted_idx))
        idx = sorted_idx(i);
        fprintf('%d. %s: %s - SOC range %.2f%% (%.2f%% ~ %.2f%%), Current range %.2fA (%.2fA ~ %.2fA), Length: %d points\n', ...
            i, results(idx).file_name, results(idx).cycle_type, results(idx).soc_range, ...
            results(idx).soc_min, results(idx).soc_max, results(idx).current_range, results(idx).current_min, results(idx).current_max, results(idx).data_points);
    end
    
         %% Plot continuous charge/discharge segments
     fprintf('\n=== Plotting Continuous Charge/Discharge Segments ===\n');
     
     % Plot top 3 results (2 subplots each)
     num_plots = min(3, length(sorted_idx));
     
     for plot_idx = 1:num_plots
         idx = sorted_idx(plot_idx);
         
         % Create figure for each result
         figure('Position', [100, 100, 1400, 600]);
         
         % Get segment data
         soc_data = results(idx).soc_data;
         current_data = results(idx).current_data;
         voltage_data = results(idx).voltage_data;
         power_data = results(idx).power_data;
         time_data = results(idx).time_data;
         start_idx = results(idx).start_idx;
         end_idx = results(idx).end_idx;
         

         
         % Extract data for the segment (include 1 point before start)
         plot_start_idx = max(1, start_idx - 1);
         
         % Extract extended data from the filtered data
         time_data_extended = time_data;
         soc_data_extended = soc_data;
         current_data_extended = current_data;
         
         % If we can include one point before, do so
         if plot_start_idx < start_idx
             % Get the original filtered data
             original_filtered_time = results(idx).original_time_data;
             original_filtered_soc = results(idx).original_soc_data;
             original_filtered_current = results(idx).original_current_data;
             
             % Add the previous point to the data
             time_data_extended = [original_filtered_time(plot_start_idx); time_data];
             soc_data_extended = [original_filtered_soc(plot_start_idx); soc_data];
             current_data_extended = [original_filtered_current(plot_start_idx); current_data];
         end
         
         % Extract voltage data for the segment
         if plot_start_idx < start_idx
             % Get the original filtered voltage data
             original_filtered_voltage = results(idx).original_voltage_data;
             % Include one point before
             segment_voltage = [original_filtered_voltage(plot_start_idx); voltage_data];
         else
             segment_voltage = voltage_data;
         end
         
         % Extract power data for the segment
         if plot_start_idx < start_idx
             % Get the original filtered power data
             original_filtered_power = results(idx).original_power_data;
             % Include one point before
             segment_power = [original_filtered_power(plot_start_idx); power_data];
         else
             segment_power = power_data;
         end
         
         % Subplot 1: Voltage and Current vs Time
         subplot(1, 2, 1);
         
         % Plot voltage on left y-axis
         yyaxis left;
         plot(time_data_extended, segment_voltage, 'g-', 'LineWidth', 1.5);
         ylabel('Voltage (V)');
         grid on;
         
         % Plot current on right y-axis
         yyaxis right;
         plot(time_data_extended, current_data_extended, 'r-', 'LineWidth', 1.5);
         ylabel('Current (A)');
         
         % Set title and labels
         title(sprintf('%s: %s Segment - Voltage & Current vs Time', results(idx).file_name, results(idx).cycle_type));
         xlabel('Time (HH:MM:SS)');
         
         % Format x-axis to show time in HH:MM:SS format
         if isa(time_data_extended, 'duration')
             % Set ticks at regular intervals
             num_ticks = min(8, length(time_data_extended));
             tick_indices = round(linspace(1, length(time_data_extended), num_ticks));
             xticks(time_data_extended(tick_indices));
             xticklabels(string(time_data_extended(tick_indices), 'hh:mm:ss'));
         elseif isa(time_data_extended, 'datetime')
             datetick('x', 'HH:MM:SS');
         else
             % For numeric data, use regular ticks
             num_ticks = min(8, length(time_data_extended));
             tick_indices = round(linspace(1, length(time_data_extended), num_ticks));
             xticks(tick_indices);
             xticklabels(arrayfun(@(x) sprintf('%.0f', x), time_data_extended(tick_indices), 'UniformOutput', false));
         end
         
         % Add legend
         legend('Voltage', 'Current', 'Location', 'best');
         
         % Subplot 2: Q vs SOC
         subplot(1, 2, 2);
         
         % Calculate Q = I*dt (charge/discharge capacity in Ah)
         dt = diff(time_data_extended); % Time differences
         dt_seconds = seconds(dt); % Convert to seconds
         
         % Calculate Q for each interval (A*s/3600 = Ah)
         Q_intervals = abs(current_data_extended(1:end-1)) .* dt_seconds / 3600;
         
         % Cumulative Q (charge/discharge capacity in Ah)
         Q_cumulative = cumsum(Q_intervals);
         
         % Plot Voltage vs Q (excluding the last voltage point since Q has one less point)
         plot(Q_cumulative, segment_voltage(1:end-1), 'b-', 'LineWidth', 1.5);
         ylabel('Voltage (V)');
         xlabel('Q (Ah)');
         grid on;
         
         % Set title and labels
         title(sprintf('Voltage vs Q - SOC: %.2f%% (%.2f%%~%.2f%%), Current: %.2fA (%.2fA~%.2fA)', ...
             results(idx).soc_range, results(idx).soc_min, results(idx).soc_max, ...
             results(idx).current_range, results(idx).current_min, results(idx).current_max));
         
         % Add legend
         legend('Voltage', 'Location', 'best');
         
         % Save individual plot
         plot_file = fullfile(saveDir, sprintf('Segment_%d_%s_%s_SOC%d%%_Current%dA.fig', ...
             plot_idx, results(idx).file_name(1:end-4), results(idx).cycle_type, ...
             soc_variation_threshold, current_variation_threshold));
         savefig(plot_file);
         fprintf('Plot %d saved to: %s\n', plot_idx, plot_file);
     end
    
         fprintf('All individual plots saved to: %s\n', saveDir);
    
else
    fprintf('No data found matching the criteria.\n');
end

%% Function to find continuous segments
function segments = find_continuous_segments(current_data, soc_data, type, soc_threshold, current_threshold)
    segments = struct('type', {}, 'start_idx', {}, 'end_idx', {}, 'soc_data', {}, 'current_data', {}, ...
                     'time_data', {}, 'soc_range', {}, 'soc_min', {}, 'soc_max', {}, ...
                     'current_range', {}, 'current_min', {}, 'current_max', {}, 'length', {});
    
    n = length(current_data);
    if n < 10
        return;
    end
    
    % Define conditions based on type
    if strcmp(type, 'charging')
        % Charging: positive current, SOC increases or stays the same
        current_condition = @(curr) curr > 0;
        soc_condition = @(soc_prev, soc_curr) soc_curr >= soc_prev;
    else % discharging
        % Discharging: negative current, SOC decreases or stays the same
        current_condition = @(curr) curr < 0;
        soc_condition = @(soc_prev, soc_curr) soc_curr <= soc_prev;
    end
    
    % Find continuous segments
    start_idx = 1;
    current_segment_start = -1;
    
    for i = 1:n
        curr = current_data(i);
        
        % Check if current meets the condition
        if current_condition(curr)
            if current_segment_start == -1
                current_segment_start = i;
            end
        else
            % Current condition violated, check if we have a valid segment
            if current_segment_start ~= -1
                end_idx = i - 1;
                
                % Check SOC condition for the segment
                segment_soc = soc_data(current_segment_start:end_idx);
                segment_current = current_data(current_segment_start:end_idx);
                
                % Verify SOC condition throughout the segment
                soc_valid = true;
                for j = 2:length(segment_soc)
                    if ~soc_condition(segment_soc(j-1), segment_soc(j))
                        soc_valid = false;
                        break;
                    end
                end
                
                if soc_valid && length(segment_soc) >= 10
                    % Calculate variations
                    soc_range = max(segment_soc) - min(segment_soc);
                    current_range = max(segment_current) - min(segment_current);
                    
                    % Check if segment meets thresholds
                    if soc_range >= soc_threshold && current_range >= current_threshold
                        % Create segment struct
                        segment = struct();
                        segment.type = type;
                        segment.start_idx = current_segment_start;
                        segment.end_idx = end_idx;
                        segment.soc_data = segment_soc;
                        segment.current_data = segment_current;
                        segment.time_data = []; % Will be filled by caller
                        segment.soc_range = soc_range;
                        segment.soc_min = min(segment_soc);
                        segment.soc_max = max(segment_soc);
                        segment.current_range = current_range;
                        segment.current_min = min(segment_current);
                        segment.current_max = max(segment_current);
                        segment.length = length(segment_soc);
                        
                        segments = [segments; segment];
                    end
                end
                
                current_segment_start = -1;
            end
        end
    end
    
    % Check for segment at the end
    if current_segment_start ~= -1
        end_idx = n;
        segment_soc = soc_data(current_segment_start:end_idx);
        segment_current = current_data(current_segment_start:end_idx);
        
        % Verify SOC condition throughout the segment
        soc_valid = true;
        for j = 2:length(segment_soc)
            if ~soc_condition(segment_soc(j-1), segment_soc(j))
                soc_valid = false;
                break;
            end
        end
        
        if soc_valid && length(segment_soc) >= 10
            % Calculate variations
            soc_range = max(segment_soc) - min(segment_soc);
            current_range = max(segment_current) - min(segment_current);
            
            % Check if segment meets thresholds
            if soc_range >= soc_threshold && current_range >= current_threshold
                % Create segment struct
                segment = struct();
                segment.type = type;
                segment.start_idx = current_segment_start;
                segment.end_idx = end_idx;
                segment.soc_data = segment_soc;
                segment.current_data = segment_current;
                segment.time_data = []; % Will be filled by caller
                segment.soc_range = soc_range;
                segment.soc_min = min(segment_soc);
                segment.soc_max = max(segment_soc);
                segment.current_range = current_range;
                segment.current_min = min(segment_current);
                segment.current_max = max(segment_current);
                segment.length = length(segment_soc);
                
                segments = [segments; segment];
            end
        end
    end
end

