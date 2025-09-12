%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESS Field Data rest OCV reconstruction (New Rack)
% DCCurrent < abs(idle_thresh = 1A)
% rest duration > 30 minutes (0.5 hours)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory
dataDir  = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\Rack_raw2mat\New';
yearList = {'2024'}; 
saveDir  = fullfile('G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\Dataset_Overview\Results');

if ~exist(saveDir, 'dir')
    mkdir(saveDir); 
end

fprintf('=====Raw File Traversal=====\n');

% Initialize arrays to store OCV data
all_ocv_soc = [];
all_ocv_voltage = [];
all_ocv_time = [];
all_ocv_year = [];
all_ocv_month = [];

for y = 1:length(yearList)
    year = yearList{y};
    yearPath = fullfile(dataDir, year);
    
    % Find existing month directories
    monthDirs = dir(fullfile(yearPath, '2024*'));
    monthList = {monthDirs.name};
    
    fprintf('Found months: %s\n', strjoin(monthList, ', '));
    
    % Sort months chronologically
    monthList = sort(monthList);
    
    for m = 1:length(monthList)
        monthPath = fullfile(yearPath, monthList{m});
        matFiles = dir(fullfile(monthPath, 'Raw_*.mat'));

        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            [~, name, ~] = fileparts(matFiles(f).name);
            dateStr = extractAfter(name, 'Raw_');
            fieldName = ['date_', dateStr];

            fprintf('Processing: %s\n', matFiles(f).name);
            
            load(matFilePath); 
            t          = Raw.Date_Time_seconds;
            V          = Raw.CVavg;
            I          = Raw.DCCurrent;
            soc        = Raw.SOC_BMS;
            
            % Convert datetime to seconds if needed
            if isdatetime(t)
                t_seconds = seconds(t - t(1));  % Convert to seconds from start
                fprintf('  Loaded data - Time: %d points\n', length(t));
                fprintf('  Time range: %.1f to %.1f seconds\n', min(t_seconds), max(t_seconds));
            else
                t_seconds = t;  % Already in seconds
                fprintf('  Loaded data - Time: %d points\n', length(t));
                fprintf('  Time range: %.1f to %.1f seconds\n', min(t_seconds), max(t_seconds));
            end
                
            % Find rest periods (DCCurrent < abs(1A) with duration > 30 minutes)
            rest_periods = findRestPeriods(t_seconds, I, V, soc);
            
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
                    % Extract actual month from folder name (e.g., '202405' -> 5)
                    actual_month = str2double(monthList{m}(5:6));
                    all_ocv_month = [all_ocv_month; actual_month]; % Store actual month number
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
all_ocv_month = all_ocv_month(valid_idx);

% Save OCV data
ocv_data = struct();
ocv_data.soc = all_ocv_soc;
ocv_data.voltage = all_ocv_voltage;
ocv_data.time = all_ocv_time;
ocv_data.year = all_ocv_year;
ocv_data.month = all_ocv_month;
ocv_data.year_list = yearList;

save(fullfile(saveDir, 'NewRack_OCV_Rest_Data.mat'), 'ocv_data');

%% Plot OCV curve by year
figure('Position', [100, 100, 1000, 800]);

% Define colors for each month (12 colors)
all_colors = [
    32 133 78;   % 1월 - 녹색
    7 115 194;   % 2월 - 파란색
    205 83 76;   % 3월 - 빨간색
    255 193 7;   % 4월 - 노란색
    156 39 176;  % 5월 - 보라색
    255 87 34;   % 6월 - 주황색
    0 150 136;   % 7월 - 청록색
    233 30 99;   % 8월 - 분홍색
    121 85 72;   % 9월 - 갈색
    63 81 181;   % 10월 - 인디고
    255 152 0;   % 11월 - 주황색
    76 175 80    % 12월 - 연두색
] / 255;

all_month_names = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
                   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

% Main OCV plot
subplot(2,1,1);
hold on;
for month = 1:max(all_ocv_month)
    month_idx = all_ocv_month == month;
    if sum(month_idx) > 0
        scatter(all_ocv_soc(month_idx), all_ocv_voltage(month_idx), 30, all_colors(month,:), 'filled', 'DisplayName', all_month_names{month});
    end
end
xlabel('SOC [%]');
xlim([0 100]);
ylabel('Voltage [V]');
ylim([2.5 4.5]);
title('OCV Curve from Rest Periods (>30min, |I|<1A) - 2024 by Month');
legend('show', 'Location', 'southeast');
grid on;
hold off;

% SOC distribution by month
subplot(2,1,2);
hold on;
for month = 1:max(all_ocv_month)
    month_idx = all_ocv_month == month;
    if sum(month_idx) > 0
        histogram(all_ocv_soc(month_idx), 15, 'FaceColor', all_colors(month,:), 'FaceAlpha', 0.7, 'DisplayName', all_month_names{month});
    end
end
xlabel('SOC [%]');
xlim([0 100]);
ylabel('Count');
title('SOC Distribution - 2024 by Month');
legend('show', 'Location', 'southeast');
grid on;
hold off;

% Save figure
saveas(gcf, fullfile(saveDir, 'NewRack_OCV_Rest_Curve.fig'));

%% Create averaged OCV curve by year
figure('Position', [100, 100, 1000, 800]);

% Debug: Check if we have data for plotting
fprintf('\n=====Debugging Second Figure=====\n');
fprintf('Total data points: %d\n', length(all_ocv_soc));
fprintf('Available months: %s\n', strjoin(all_month_names(1:max(all_ocv_month)), ', '));

hold on;
for month = 1:max(all_ocv_month)
    month_idx = all_ocv_month == month;
    if sum(month_idx) > 0
        month_soc = all_ocv_soc(month_idx);
        month_voltage = all_ocv_voltage(month_idx);
        
        % Debug: Check month data
        fprintf('Month %s: %d points, SOC range %.3f to %.3f, Voltage range %.3f to %.3f\n', ...
            all_month_names{month}, sum(month_idx), min(month_soc), max(month_soc), min(month_voltage), max(month_voltage));

        % Group by similar SOC values and average OCV
        unique_soc = unique(round(month_soc * 100) / 100); % Round to 2 decimal places
        avg_ocv_by_soc = zeros(length(unique_soc), 1);

        fprintf('  Unique SOC values: %d\n', length(unique_soc));

        for i = 1:length(unique_soc)
            soc_target = unique_soc(i);
            % Find all points within 0.1 of this SOC value
            similar_soc_idx = abs(month_soc - soc_target) < 0.1;
            if sum(similar_soc_idx) > 0
                avg_ocv_by_soc(i) = mean(month_voltage(similar_soc_idx));
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
                    'Color', all_colors(month,:), 'LineWidth', 2, 'DisplayName', all_month_names{month}, 'HandleVisibility', 'off');
                start_idx = gap_idx + 1;
            end

            % Plot final segment
            if start_idx <= length(sorted_soc)
                plot(sorted_soc(start_idx:end), sorted_voltage(start_idx:end), 'o-', 'MarkerSize', 2,...
                    'Color', all_colors(month,:), 'LineWidth', 2, 'DisplayName', all_month_names{month});
            end
        else
            % Single point - just plot it
            plot(sorted_soc, sorted_voltage, 'o-', 'MarkerSize', 2,...
                'Color', all_colors(month,:), 'LineWidth', 2, 'DisplayName', all_month_names{month});
        end
    end
end


xlabel('SOC [%]');
xlim([0 100]);
ylabel('Voltage [V]');
ylim([2.5 4.5]);
title('Averaged OCV Curve - 2024 by Month (>30min, |I|<1A)');
legend('show', 'Location', 'southeast');
grid on;
hold off;

% Save averaged OCV figure
saveas(gcf, fullfile(saveDir, 'NewRack_OCV_Averaged_Curve.fig'));

fprintf('OCV data saved to: %s\n', saveDir);
fprintf('Total valid OCV points: %d\n', length(all_ocv_soc));

%% Function to find rest periods
function rest_periods = findRestPeriods(t, DCCurrent, V, soc)
    rest_periods = {};
    idle_thresh = 1.0; % 1A threshold
    
    % Debug: Check data validity
    fprintf('Data check - Time: %d points, V: %d points, SOC: %d points\n', ...
            length(t), length(V), length(soc));
    fprintf('V range: %.3f to %.3f V, SOC range: %.3f to %.3f\n', ...
            min(V), max(V), min(soc), max(soc));
    fprintf('Current range: %.3f to %.3f A\n', min(DCCurrent), max(DCCurrent));
    
    % Find idle periods (|DCCurrent| < idle_thresh)
    idle_idx = abs(DCCurrent) < idle_thresh;
    
    if sum(idle_idx) == 0
        fprintf('No idle periods (|I|<%.1fA) found in this file\n', idle_thresh);
        return;
    end
    
    % Find continuous idle periods
    idle_start = find(diff([0; idle_idx]) == 1);
    idle_end = find(diff([idle_idx; 0]) == -1);
    
    fprintf('Found %d idle periods\n', length(idle_start));
    
    for i = 1:length(idle_start)
        start_idx = idle_start(i);
        end_idx = idle_end(i);
        
        % Calculate duration (assuming t is in seconds)
        duration = t(end_idx) - t(start_idx);
        
        % Debug: Show all periods regardless of duration
        fprintf('  Period %d: duration = %.1f minutes (%.1f hours)\n', ...
                i, duration/60, duration/3600);
        
        % Check if duration is greater than 30 minutes (1800 seconds)
        if duration >= 1800
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
    
    fprintf('  Valid periods (>30min): %d\n', length(rest_periods));
end 