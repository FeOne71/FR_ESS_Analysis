%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New Kimje Data Rated Measurement Profile Search
% 202106 
% 202206
% 202306
% For charging    : I>0, I<0 break
% For discharging : I<0, I>0 break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
yearList = {'2023'}; 
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\Dataset_Overview\FieldData_RPT_OldDataCheck');

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameters
C_nom_rack = 128;
soc_variation_threshold = 10;  % SOC variation range threshold for charge/discharge cycle (%)
current_variation_threshold = 50;  % Current variation threshold for charge/discharge cycle (A)
min_duration_minutes = 10;  % Minimum duration for continuous segments (minutes)

%% Find continuous charge/discharge cycles
fprintf('Searching for continuous charge/discharge cycles with SOC variation >= %d%%, current variation >= %dA, and duration >= %d minutes...\n', soc_variation_threshold, current_variation_threshold, min_duration_minutes);

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
    
    fprintf('Found %d month folders in %s\n', length(month_folders), year_path);
    
    for month_idx = 1:length(month_folders)
        month_folder = month_folders(month_idx).name;
        month_path = fullfile(year_path, month_folder);
        
        fprintf('Processing month: %s\n', month_folder);
        
        % Find daily mat files
        mat_files = dir(fullfile(month_path, 'Raw_*.mat'));
        fprintf('Found %d Raw_*.mat files in %s\n', length(mat_files), month_path);
        
        for file_idx = 1:length(mat_files)
            file_name = mat_files(file_idx).name;
            file_path = fullfile(month_path, file_name);
            
            try
                % Load mat file
                fprintf('Processing: %s\n', file_path);
                data = load(file_path);
                
                % Debug: Check what variables are in the loaded data
                fprintf('  -> Variables in data: ');
                data_fields = fieldnames(data);
                for i = 1:length(data_fields)
                    fprintf('%s ', data_fields{i});
                end
                fprintf('\n');
                
                % Extract data from Raw_Rack
                if isfield(data, 'Raw_Rack')
                    event_data = data.Raw_Rack;
                    fprintf('  -> Raw_Rack found\n');
                else
                    fprintf('  -> Raw_Rack not found, trying Raw\n');
                    if isfield(data, 'Raw')
                        event_data = data.Raw;
                    else
                        fprintf('  -> Neither Raw_Rack nor Raw found\n');
                        continue;
                    end
                end
                
                % Check if Rack01 exists
                if isfield(event_data, 'Rack01')
                    fprintf('  -> Rack01 found\n');
                    rack_data = event_data.Rack01;
                else
                    fprintf('  -> Rack01 not found, using event_data directly\n');
                    rack_data = event_data;
                end
                
                % Get SOC and current data
                if isfield(rack_data, 'SOCPct')
                    soc_data = rack_data.SOCPct;
                    fprintf('  -> SOCPct found, size: %dx%d\n', size(soc_data));
                else
                    fprintf('  -> SOCPct not found\n');
                    continue;
                end
                
                if isfield(rack_data, 'DCCurrent_A')
                    current_data = rack_data.DCCurrent_A;
                    fprintf('  -> DCCurrent_A found, size: %dx%d\n', size(current_data));
                else
                    fprintf('  -> DCCurrent_A not found\n');
                    continue;
                end
                
                % Filter valid data (remove NaN, Inf)
                valid_mask = ~isnan(soc_data) & ~isinf(soc_data) & ~isnan(current_data) & ~isinf(current_data);
                valid_soc = soc_data(valid_mask);
                valid_current = current_data(valid_mask);
                if isfield(rack_data, 'Time')
                    time_strings = rack_data.Time(valid_mask);
                    fprintf('  -> Time found, size: %dx%d\n', size(rack_data.Time));
                    
                    % Convert string time to datetime
                    try
                        valid_time = datetime(time_strings, 'Format', 'yyyy-MM-dd HH:mm:ss');
                        fprintf('  -> Time converted to datetime successfully\n');
                    catch
                        fprintf('  -> Time conversion failed, using numeric indices\n');
                        valid_time = (1:length(time_strings))';
                    end
                else
                    fprintf('  -> Time not found\n');
                    continue;
                end
                
                fprintf('  -> Valid data points: %d\n', length(valid_soc));
                fprintf('  -> SOC range: %.2f%% ~ %.2f%%\n', min(valid_soc), max(valid_soc));
                fprintf('  -> Current range: %.2fA ~ %.2fA\n', min(valid_current), max(valid_current));
                
                if length(valid_soc) > 10  % Need at least 10 data points
                    % Find continuous charging segments (I >= 0 and dI >= 0)
                    charging_segments = find_continuous_segments(valid_current, valid_soc, valid_time, 'charging', soc_variation_threshold, current_variation_threshold, min_duration_minutes);
                    
                    % Find continuous discharging segments (I <= 0 and dI <= 0)
                    discharging_segments = find_continuous_segments(valid_current, valid_soc, valid_time, 'discharging', soc_variation_threshold, current_variation_threshold, min_duration_minutes);
                    
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
                         if isfield(rack_data, 'AverageCV_V')
                             segment_voltage = rack_data.AverageCV_V(valid_mask);
                             results(result_count).voltage_data = segment_voltage(segment.start_idx:segment.end_idx);
                         else
                             results(result_count).voltage_data = [];
                         end
                         
                         % Store power data for the segment (matching current data length)
                         if isfield(rack_data, 'DCPower_kW')
                             segment_power = rack_data.DCPower_kW(valid_mask);
                             results(result_count).power_data = segment_power(segment.start_idx:segment.end_idx);
                         else
                             results(result_count).power_data = [];
                         end
                         
                         % Store temperature data for the segment (matching current data length)
                         if isfield(rack_data, 'MTavg')
                             segment_temp = rack_data.MTavg(valid_mask);
                             results(result_count).temp_data = segment_temp(segment.start_idx:segment.end_idx);
                         else
                             results(result_count).temp_data = [];
                         end
                         
                         results(result_count).time_data = valid_time(segment.start_idx:segment.end_idx);
                        
                                                 % Store original filtered data for extended plotting
                         results(result_count).original_time_data = valid_time;
                         results(result_count).original_soc_data = valid_soc;
                         results(result_count).original_current_data = valid_current;
                         
                         % Store original filtered voltage, power, and temperature data
                         if isfield(rack_data, 'AverageCV_V')
                             results(result_count).original_voltage_data = rack_data.AverageCV_V(valid_mask);
                         else
                             results(result_count).original_voltage_data = [];
                         end
                         
                         if isfield(rack_data, 'DCPower_kW')
                             results(result_count).original_power_data = rack_data.DCPower_kW(valid_mask);
                         else
                             results(result_count).original_power_data = [];
                         end
                         
                         if isfield(rack_data, 'MTavg')
                             results(result_count).original_temp_data = rack_data.MTavg(valid_mask);
                         else
                             results(result_count).original_temp_data = [];
                         end
                        
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
fprintf('Found %d continuous charge/discharge segments (SOC >= %d%%, Current >= %dA, Duration >= %d min).\n', ...
    result_count, soc_variation_threshold, current_variation_threshold, min_duration_minutes);

if result_count > 0
    % Convert results to table
    result_table = struct2table(results, 'AsArray', true);
    
    % Save results
    result_file = fullfile(saveDir, sprintf('Continuous_Cycles_SOC%d%%_Current%dA_Duration%dmin_Results.mat', soc_variation_threshold, current_variation_threshold, min_duration_minutes));
    save(result_file, 'results', 'result_table', 'soc_variation_threshold', 'current_variation_threshold', 'min_duration_minutes');
    
         fprintf('Results saved to:\n');
     fprintf('  - MAT file: %s\n', result_file);
    
         % Display all results
     fprintf('\n=== All Results ===\n');
     [~, sorted_idx] = sort([results.soc_range], 'descend');
     for i = 1:length(sorted_idx)
         idx = sorted_idx(i);
         fprintf('%d. %s: %s - SOC range %.2f%% (%.2f%% ~ %.2f%%), Current range %.2fA (%.2fA ~ %.2fA), Length: %d points\n', ...
             i, results(idx).file_name, results(idx).cycle_type, results(idx).soc_range, ...
             results(idx).soc_min, results(idx).soc_max, results(idx).current_range, results(idx).current_min, results(idx).current_max, results(idx).data_points);
     end
    
         %% Plot continuous charge/discharge segments
     fprintf('\n=== Plotting Continuous Charge/Discharge Segments ===\n');
     
     % Plot all results (2 subplots each)
     num_plots = length(sorted_idx);
     
     for plot_idx = 1:num_plots
         idx = sorted_idx(plot_idx);
         
                   % Create figure for each result
          figure('Position', [100, 100, 1200, 1000]);
         
                              % Get segment data
           soc_data = results(idx).soc_data;
           current_data = results(idx).current_data;
           voltage_data = results(idx).voltage_data;
           power_data = results(idx).power_data;
           time_data = results(idx).time_data;
           start_idx = results(idx).start_idx;
           end_idx = results(idx).end_idx;
          
          % Check if data is valid for plotting
          if isempty(soc_data) || isempty(current_data) || isempty(time_data)
              fprintf('  -> Skipping plot for %s: missing data\n', results(idx).file_name);
              continue;
          end
          
          % Check if data lengths match
          if length(soc_data) ~= length(time_data) || length(current_data) ~= length(time_data)
              fprintf('  -> Skipping plot for %s: data length mismatch\n', results(idx).file_name);
              continue;
          end
          
                     % Extract data for the segment (include 1 point before start)
           % Use the filtered data that was used for segment detection
           % Include 1 point before the segment start (within filtered data)
           plot_start_idx = max(1, start_idx - 1);
           plot_end_idx = end_idx;
           
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
            if ~isempty(voltage_data) && length(voltage_data) == length(time_data)
                if plot_start_idx < start_idx
                    % Get the original filtered voltage data
                    original_filtered_voltage = results(idx).original_voltage_data;
                    if ~isempty(original_filtered_voltage) && plot_start_idx <= length(original_filtered_voltage)
                        % Include one point before
                        segment_voltage = [original_filtered_voltage(plot_start_idx); voltage_data];
                    else
                        segment_voltage = voltage_data;
                    end
                else
                    segment_voltage = voltage_data;
                end
            else
                segment_voltage = [];
            end
         
                              % Subplot 1: SOC and Current
           subplot(3, 1, 1);
         
                   % Plot SOC on left y-axis
          yyaxis left;
          plot(time_data_extended, soc_data_extended, 'b-', 'LineWidth', 1.5, 'Color', '#20854E');
          ylabel('SOC (%)');
          ylim([0, 100]);
          grid on;
          
          % Plot current on right y-axis
          yyaxis right;
          plot(time_data_extended, current_data_extended, 'r-', 'LineWidth', 1.5, 'Color', '#CD534C');
         ylabel('Current (A)');
         
         % Set title and labels
         title(sprintf('%s: %s Segment - SOC & Current', results(idx).file_name, results(idx).cycle_type));
         xlabel('Time (HH:MM:SS)');
         
         % Format x-axis to show time in HH:MM:SS format
         if isa(time_data_extended, 'datetime')
             % For datetime, use datetick
             datetick('x', 'HH:MM:SS');
         elseif isa(time_data_extended, 'duration')
             % For duration, set ticks manually
             num_ticks = min(8, length(time_data_extended));
             tick_indices = round(linspace(1, length(time_data_extended), num_ticks));
             xticks(time_data_extended(tick_indices));
             xticklabels(string(time_data_extended(tick_indices), 'hh:mm:ss'));
         else
             % For numeric data, use regular ticks
             num_ticks = min(8, length(time_data_extended));
             tick_indices = round(linspace(1, length(time_data_extended), num_ticks));
             xticks(tick_indices);
             xticklabels(arrayfun(@(x) sprintf('%.0f', x), time_data_extended(tick_indices), 'UniformOutput', false));
         end
         
         % Add legend
         legend('SOC', 'Current', 'Location', 'best');
         
                              % Subplot 2: Current and Voltage
           subplot(3, 1, 2);
         
                   % Plot current on left y-axis
          yyaxis left;
          plot(time_data_extended, current_data_extended, 'r-', 'LineWidth', 1.5, 'Color', '#CD534C');
          ylabel('Current (A)');
          grid on;
          
          % Plot voltage on right y-axis
          yyaxis right;
          if ~isempty(segment_voltage) && length(segment_voltage) == length(time_data_extended)
              plot(time_data_extended, segment_voltage, 'g-', 'LineWidth', 1.5, 'Color', '#0073C2');
              ylabel('Voltage (V)');
          else
              ylabel('Voltage (V) - No data');
          end
         
         % Set title and labels
         title(sprintf('Current & Voltage - SOC: %.2f%% (%.2f%%~%.2f%%), Current: %.2fA (%.2fA~%.2fA)', ...
             results(idx).soc_range, results(idx).soc_min, results(idx).soc_max, ...
             results(idx).current_range, results(idx).current_min, results(idx).current_max));
         xlabel('Time (HH:MM:SS)');
         
         % Format x-axis to show time in HH:MM:SS format
         if isa(time_data_extended, 'datetime')
             % For datetime, use datetick
             datetick('x', 'HH:MM:SS');
         elseif isa(time_data_extended, 'duration')
             % For duration, set ticks manually
             num_ticks = min(8, length(time_data_extended));
             tick_indices = round(linspace(1, length(time_data_extended), num_ticks));
             xticks(time_data_extended(tick_indices));
             xticklabels(string(time_data_extended(tick_indices), 'hh:mm:ss'));
         else
             % For numeric data, use regular ticks
             num_ticks = min(8, length(time_data_extended));
             tick_indices = round(linspace(1, length(time_data_extended), num_ticks));
             xticks(tick_indices);
             xticklabels(arrayfun(@(x) sprintf('%.0f', x), time_data_extended(tick_indices), 'UniformOutput', false));
         end
         
                              % Add legend
           if ~isempty(segment_voltage) && length(segment_voltage) == length(time_data_extended)
               legend('Current', 'Voltage', 'Location', 'best');
           else
               legend('Current', 'Location', 'best');
           end
           
           % Subplot 3: Power vs Time
           subplot(3, 1, 3);
           
                       % Extract power data for the segment
            if ~isempty(power_data) && length(power_data) == length(time_data)
                if plot_start_idx < start_idx
                    % Get the original filtered power data
                    original_filtered_power = results(idx).original_power_data;
                    if ~isempty(original_filtered_power) && plot_start_idx <= length(original_filtered_power)
                        % Include one point before
                        segment_power = [original_filtered_power(plot_start_idx); power_data];
                    else
                        segment_power = power_data;
                    end
                else
                    segment_power = power_data;
                end
            else
                segment_power = [];
            end
           
                       % Plot power vs time
            if ~isempty(segment_power) && length(segment_power) == length(time_data_extended)
                plot(time_data_extended, segment_power, 'm-', 'LineWidth', 1.5, 'Color', '#7E2F8E');
                ylabel('Power (kW)');
                grid on;
            else
                ylabel('Power (kW) - No data');
            end
           
           % Set title and labels
           title(sprintf('Power vs Time - SOC: %.2f%% (%.2f%%~%.2f%%), Current: %.2fA (%.2fA~%.2fA)', ...
               results(idx).soc_range, results(idx).soc_min, results(idx).soc_max, ...
               results(idx).current_range, results(idx).current_min, results(idx).current_max));
           xlabel('Time (HH:MM:SS)');
           
           % Format x-axis to show time in HH:MM:SS format
           if isa(time_data_extended, 'datetime')
               % For datetime, use datetick
               datetick('x', 'HH:MM:SS');
           elseif isa(time_data_extended, 'duration')
               % For duration, set ticks manually
               num_ticks = min(8, length(time_data_extended));
               tick_indices = round(linspace(1, length(time_data_extended), num_ticks));
               xticks(time_data_extended(tick_indices));
               xticklabels(string(time_data_extended(tick_indices), 'hh:mm:ss'));
           else
               % For numeric data, use regular ticks
               num_ticks = min(8, length(time_data_extended));
               tick_indices = round(linspace(1, length(time_data_extended), num_ticks));
               xticks(tick_indices);
               xticklabels(arrayfun(@(x) sprintf('%.0f', x), time_data_extended(tick_indices), 'UniformOutput', false));
           end
           
                       % Add legend
            if ~isempty(segment_power) && length(segment_power) == length(time_data_extended)
                legend('Power (kW)', 'Location', 'best');
            else
                % Don't add legend if no power data
            end
         
         % Save individual plot
         plot_file = fullfile(saveDir, sprintf('Segment_%d_%s_%s_SOC%d%%_Current%dA_Duration%dmin.fig', ...
             plot_idx, results(idx).file_name(1:end-4), results(idx).cycle_type, ...
             soc_variation_threshold, current_variation_threshold, min_duration_minutes));
         savefig(plot_file);
         fprintf('Plot %d saved to: %s\n', plot_idx, plot_file);
     end
    
         fprintf('All individual plots saved to: %s\n', saveDir);
    
else
    fprintf('No data found matching the criteria.\n');
end

%% Function to find continuous segments
function segments = find_continuous_segments(current_data, soc_data, time_data, type, soc_threshold, current_threshold, min_duration_minutes)
    segments = struct('type', {}, 'start_idx', {}, 'end_idx', {}, 'soc_data', {}, 'current_data', {}, ...
                     'time_data', {}, 'soc_range', {}, 'soc_min', {}, 'soc_max', {}, ...
                     'current_range', {}, 'current_min', {}, 'current_max', {}, 'length', {});
    
    n = length(current_data);
    if n < 10
        return;
    end
    

    
    % Define conditions based on type
    if strcmp(type, 'charging')
        % Charging: I > 0 (positive current)
        current_condition = @(curr) curr > 0;
    else % discharging
        % Discharging: I < 0 (negative current)
        current_condition = @(curr) curr < 0;
    end
    
    % Find continuous segments
    current_segment_start = -1;
    
    for i = 1:n
        curr = current_data(i);
        
        % Check if current meets the condition
        if current_condition(curr)
            % Start new segment if not already started
            if current_segment_start == -1
                current_segment_start = i;
            end
        else
            % Current condition violated (sign changed), check if we have a valid segment
            if current_segment_start ~= -1
                end_idx = i - 1;
                
                % Check segment duration and other conditions
                segment_duration = time_data(end_idx) - time_data(current_segment_start);
                segment_duration_minutes = minutes(segment_duration);
                
                if segment_duration_minutes >= min_duration_minutes
                    % Check SOC condition for the segment
                    segment_soc = soc_data(current_segment_start:end_idx);
                    segment_current = current_data(current_segment_start:end_idx);
                    
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
                        segment.time_data = time_data(current_segment_start:end_idx);
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
        
        % Check segment duration and other conditions
        segment_duration = time_data(end_idx) - time_data(current_segment_start);
        segment_duration_minutes = minutes(segment_duration);
        
        if segment_duration_minutes >= min_duration_minutes
            % Check SOC condition for the segment
            segment_soc = soc_data(current_segment_start:end_idx);
            segment_current = current_data(current_segment_start:end_idx);
            
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
                segment.time_data = time_data(current_segment_start:end_idx);
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

