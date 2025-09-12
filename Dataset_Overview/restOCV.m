%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESS Field Data rest OCV reconstruction
% BSC_Charge -> "IDLE"
% rest duration > 0.5 hr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory
dataDir  = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\Rack_raw2mat\New';
yearList = {'2023', '2024', '2025'}; 
saveDir  = fullfile('G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\Dataset_Overview');

if ~exist(saveDir, 'dir')
    mkdir(saveDir); 
end

fprintf('=====Raw File Traversal=====\n');

% Initialize arrays to store OCV data
all_ocv_soc = [];
all_ocv_voltage = [];
all_ocv_time = [];
all_ocv_year = [];

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
            V          = Raw.CVavg;
            soc        = Raw.SOC_BMS;
            status     = Raw.Charge;  % string array
            
            % Convert datetime to seconds if needed
            if isdatetime(t)
                t_seconds = seconds(t - t(1));  % Convert to seconds from start
                fprintf('  Loaded data - Time: %d points, Charge states: %s\n', ...
                        length(t), strjoin(unique(status), ', '));
                fprintf('  Time range: %.1f to %.1f seconds\n', min(t_seconds), max(t_seconds));
            else
                t_seconds = t;  % Already in seconds
                fprintf('  Loaded data - Time: %d points, Charge states: %s\n', ...
                        length(t), strjoin(unique(status), ', '));
                fprintf('  Time range: %.1f to %.1f seconds\n', min(t_seconds), max(t_seconds));
            end
                
                % Find rest periods (IDLE status with duration > 1 hour)
                rest_periods = findRestPeriods(t_seconds, status, V, soc);
            
                            % Extract final values from each rest period
                for i = 1:length(rest_periods)
                    period = rest_periods{i};
                    if ~isempty(period)
                        % Get the last value from the rest period
                        final_idx = period(end);
                        all_ocv_soc = [all_ocv_soc; soc(final_idx)];
                        all_ocv_voltage = [all_ocv_voltage; V(final_idx)];
                        all_ocv_time = [all_ocv_time; t_seconds(final_idx)];
                        all_ocv_year = [all_ocv_year; y]; % Store year index
                    end
                end
        end
    end
end

%% Process and save OCV data
fprintf('\n=====Processing OCV Data=====\n');
fprintf('Total rest periods found: %d\n', length(all_ocv_soc));

% Remove any invalid data points
valid_idx = ~isnan(all_ocv_soc) & ~isnan(all_ocv_voltage);
all_ocv_soc = all_ocv_soc(valid_idx);
all_ocv_voltage = all_ocv_voltage(valid_idx);
all_ocv_time = all_ocv_time(valid_idx);
all_ocv_year = all_ocv_year(valid_idx);

% Save OCV data
ocv_data = struct();
ocv_data.soc = all_ocv_soc;
ocv_data.voltage = all_ocv_voltage;
ocv_data.time = all_ocv_time;
ocv_data.year = all_ocv_year;
ocv_data.year_list = yearList;

save(fullfile(saveDir, 'OCV_Rest_Data.mat'), 'ocv_data');

%% Plot OCV curve by year
figure('Position', [100, 100, 1000, 800]);

% Define colors for each year
colors = [32 133 78; 7 115 194; 205 83 76] / 255;
year_names = {'2021', '2022', '2023'};

% Main OCV plot
subplot(2,1,1);
hold on;
for y = 1:length(yearList)
    year_idx = all_ocv_year == y;
    if sum(year_idx) > 0
        scatter(all_ocv_soc(year_idx), all_ocv_voltage(year_idx), 30, colors(y,:), 'filled', 'DisplayName', year_names{y});
    end
end
xlabel('SOC [%]');
xlim([0 100]);
ylabel('Voltage [V]');
ylim([820 1000]);
title('OCV Curve from Rest Periods (>0.5hr) by Year');
legend('show', 'Location', 'southeast');
grid on;
hold off;

% SOC distribution by year
subplot(2,1,2);
hold on;
for y = 1:length(yearList)
    year_idx = all_ocv_year == y;
    if sum(year_idx) > 0
        histogram(all_ocv_soc(year_idx), 15, 'FaceColor', colors(y,:), 'FaceAlpha', 0.7, 'DisplayName', year_names{y});
    end
end
xlabel('SOC [%]');
xlim([0 100]);
ylabel('Count');
title('SOC Distribution by Year');
legend('show', 'Location', 'southeast');
grid on;
hold off;

% Save figure
saveas(gcf, fullfile(saveDir, 'OCV_Rest_Curve.fig'));

%% Create averaged OCV curve by year
figure('Position', [100, 100, 1000, 800]);

% Debug: Check if we have data for plotting
fprintf('\n=====Debugging Second Figure=====\n');
fprintf('Total data points: %d\n', length(all_ocv_soc));
for y = 1:length(yearList)
    year_count = sum(all_ocv_year == y);
    fprintf('Year %s: %d points\n', year_names{y}, year_count);
end

hold on;
for y = 1:length(yearList)
    year_idx = all_ocv_year == y;
    if sum(year_idx) > 0
        year_soc = all_ocv_soc(year_idx);
        year_voltage = all_ocv_voltage(year_idx);
        
        % Debug: Check year data
        fprintf('Year %s: SOC range %.3f to %.3f, Voltage range %.3f to %.3f\n', ...
                year_names{y}, min(year_soc), max(year_soc), min(year_voltage), max(year_voltage));
        
        % Group by similar SOC values and average OCV
        unique_soc = unique(round(year_soc * 100) / 100); % Round to 2 decimal places
        avg_ocv_by_soc = zeros(length(unique_soc), 1);
        
        fprintf('  Unique SOC values: %d\n', length(unique_soc));
        
        for i = 1:length(unique_soc)
            soc_target = unique_soc(i);
            % Find all points within 0.1 of this SOC value
            similar_soc_idx = abs(year_soc - soc_target) < 0.1;
            if sum(similar_soc_idx) > 0
                avg_ocv_by_soc(i) = mean(year_voltage(similar_soc_idx));
                fprintf('    SOC %.3f: %d points, avg voltage = %.3f\n', ...
                        soc_target, sum(similar_soc_idx), avg_ocv_by_soc(i));
            end
        end
        
        % Sort by SOC for proper line plotting
        [sorted_soc, sort_idx] = sort(unique_soc);
        sorted_voltage = avg_ocv_by_soc(sort_idx);
        
        fprintf('  Sorted SOC points: %d\n', length(sorted_soc));
        
        % Plot averaged OCV curve with gap detection (5% threshold)
        if length(sorted_soc) > 1
            % Find gaps larger than 5% (relative to SOC value)
            gaps = diff(sorted_soc) > (sorted_soc(1:end-1) * 0.05);
            gap_indices = find(gaps);
            
            % Plot segments between gaps
            start_idx = 1;
            for gap_idx = gap_indices'
                % Plot segment up to this gap
                plot(sorted_soc(start_idx:gap_idx), sorted_voltage(start_idx:gap_idx), 'o-', 'MarkerSize', 2,...
                     'Color', colors(y,:), 'LineWidth', 2, 'DisplayName', year_names{y}, 'HandleVisibility', 'off');
                start_idx = gap_idx + 1;
            end
            
            % Plot final segment
            if start_idx <= length(sorted_soc)
                plot(sorted_soc(start_idx:end), sorted_voltage(start_idx:end), 'o-', 'MarkerSize', 2,...
                     'Color', colors(y,:), 'LineWidth', 2, 'DisplayName', year_names{y});
            end
        else
            % Single point - just plot it
            plot(sorted_soc, sorted_voltage, 'o-', 'MarkerSize', 2,...
                 'Color', colors(y,:), 'LineWidth', 2, 'DisplayName', year_names{y});
        end
    end
end

xlabel('SOC [%]');
xlim([0 100]);
ylabel('Voltage [V]');
ylim([820 1000]);
title('Averaged OCV Curve by Year ');
legend('show', 'Location', 'southeast');
grid on;
hold off;

% Save averaged OCV figure
saveas(gcf, fullfile(saveDir, 'OCV_Averaged_Curve.fig'));

fprintf('OCV data saved to: %s\n', saveDir);
fprintf('Total valid OCV points: %d\n', length(all_ocv_soc));

%% Function to find rest periods
function rest_periods = findRestPeriods(t, bsc_charge, V, soc)
    rest_periods = {};
    
    % Debug: Check data validity
    fprintf('Data check - Time: %d points, V: %d points, SOC: %d points\n', ...
            length(t), length(V), length(soc));
    fprintf('V range: %.3f to %.3f V, SOC range: %.3f to %.3f\n', ...
            min(V), max(V), min(soc), max(soc));
    
    % Find IDLE periods
    idle_idx = strcmp(bsc_charge, 'Idle');
    
    if sum(idle_idx) == 0
        fprintf('No IDLE periods found in this file\n');
        return;
    end
    
    % Find continuous IDLE periods
    idle_start = find(diff([0; idle_idx]) == 1);
    idle_end = find(diff([idle_idx; 0]) == -1);
    
    fprintf('Found %d IDLE periods\n', length(idle_start));
    
    for i = 1:length(idle_start)
        start_idx = idle_start(i);
        end_idx = idle_end(i);
        
        % Calculate duration (assuming t is in seconds)
        duration = t(end_idx) - t(start_idx);
        
        % Debug: Show all periods regardless of duration
        fprintf('  Period %d: duration = %.1f minutes (%.1f hours)\n', ...
                i, duration/60, duration/3600);
        
        % Check if duration is greater than 30 minutes (1800 seconds)
        if duration >= 3600
            % Get the period indices
            period_idx = start_idx:end_idx;
            
            % Debug: Check V and SOC values at the end of this period
            final_idx = period_idx(end);
            fprintf('    -> VALID: final V = %.3f V, final SOC = %.3f\n', ...
                    V(final_idx), soc(final_idx));
            
            rest_periods{end+1} = period_idx;
        else
            fprintf('    -> TOO SHORT: %.1f minutes < 30 minutes\n', duration/60);
        end
    end
    
    fprintf('  Valid periods (>0.5hr): %d\n', length(rest_periods));
end