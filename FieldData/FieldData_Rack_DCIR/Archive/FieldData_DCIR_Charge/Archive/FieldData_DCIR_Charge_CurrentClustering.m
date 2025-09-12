%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field Data DCIR Charge Current Clustering
% ESS charging events and resistance analysis (current-based clustering)
% Step 1: Folder traversal
% Step 2: BSC_Charge based event detection
% Step 3: DCIR calculation (time-based + difference values)
% Step 4: Current-based clustering and visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory
dataDir  = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\raw2mat_ver04';
yearList = {'2021'}; 
saveDir  = fullfile('G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\FieldData_DCIR_Charge\Results');

if ~exist(saveDir, 'dir')
    mkdir(saveDir); 
end

%% Variables
C_nom        = 1024;          % Ah
min_charge_duration = 30;     % [s] - Charging duration (10 seconds or more)
max_P_std    = 10;            % Max power standard deviation [kW] 3.7V * 14S * 17S = 880.6V \ Power = 880.6V * 1024A = 902.9kW \ 902.9 * 0.005 = 4.5 around 5kW
max_I_std    = C_nom *0.01;  % Max current standard deviation [A] 5.12A
dt_sec       = 1;             % Sampling time [sec]

% Current range classification settings
current_bins = [50, 100, 150, 200, 250, 300];  % [A] - Charging current range boundaries
current_labels = {'range_50_100A', 'range_100_150A', 'range_150_200A', 'range_200_250A', 'range_250_300A'};  % Range labels

%% Step 1: Raw Files Traversal
fprintf('=====Step 1: Loading Field Data=====\n');
eventStruct = struct();
eventCount = 0;

total_events    = 0;
filtered_events = 0;

for y = 1:length(yearList)
    year = yearList{y};
    yearPath = fullfile(dataDir, year);
    monthDirs = dir(fullfile(yearPath, '20*'));

    for m = 1:length(monthDirs)
        if ~monthDirs(m).isdir, continue; end
        monthPath = fullfile(yearPath, monthDirs(m).name);
        matFiles = dir(fullfile(monthPath, 'Raw_*.mat'));

        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            [~, name, ~] = fileparts(matFiles(f).name);
            dateStr = extractAfter(name, 'Raw_');
            fieldName = ['date_' dateStr];
            
            load(matFilePath); 
            t          = Raw.Time;
            I          = Raw.DCCurrent_A;
            V          = Raw.Total_AverageCVSum_V;
            soc        = Raw.Total_AverageSOC;
            T_batt     = Raw.Total_AverageMT_degC;
            bsc_charge = Raw.Charge;  % string array
            status     = Raw.Status;
            dc_power   = Raw.DCPower_kW;
            
%% Step 2: Searching Idle to Charging Transition
            fprintf('=====Step 2: Searching Idle to Charging Transition =====\n');
            % 1) Define "Idle" and "Charging" status
            is_idle = strcmp(bsc_charge, 'Idle');
            is_charging = strcmp(bsc_charge, 'Charging');
            
            % Detect "Idle" -> "Charging"
            idle_to_charge = find(is_idle(1:end-1) & is_charging(2:end));
            total_events = total_events + length(idle_to_charge);
            fprintf('File: %s - Found %d potential events\n', matFiles(f).name, length(idle_to_charge));
            
            fprintf('  Total Idle points: %d\n', sum(is_idle));
            fprintf('  Total Charging points: %d\n', sum(is_charging));
            fprintf('  Idle→Charging transitions: %d\n', length(idle_to_charge));

            step2_count = 0;
            step3_count = 0;
            step4_count = 0;

%% Step 3: Events Filtering              
            fprintf('=====Step 3: Filtering Events =====\n');
            for i = 1:length(idle_to_charge)
                idx1 = idle_to_charge(i);  % Last point of "Idle"
                start_charge_idx = idx1 + 1;  % Start point of "Charging"
                
                % 2) Find end point of "Charging" period
                charge_end_idx = start_charge_idx;
                while charge_end_idx <= length(bsc_charge) && strcmp(bsc_charge(charge_end_idx), 'Charging')
                    charge_end_idx = charge_end_idx + 1;
                end
                charge_end_idx = charge_end_idx - 1;  % Last point of "Charging"
                
                % Complete event period (from Idle start to Charging end)
                start_idx = idx1;
                end_idx = charge_end_idx;
                
                step2_count = step2_count + 1;
                
                % Condition #1 Check charging duration (10 seconds or more)
                charge_duration = charge_end_idx - start_charge_idx + 1;
                if charge_duration < min_charge_duration
                    fprintf('  Event %d filtered: Charge duration too short (%ds)\n', i, charge_duration);
                    continue;
                end
                
                step3_count = step3_count + 1;
                
                t_seg   = t(start_idx:end_idx);
                I_seg   = I(start_idx:end_idx);
                V_seg   = V(start_idx:end_idx);
                soc_seg = soc(start_idx:end_idx);
                T_seg   = T_batt(start_idx:end_idx);
                bsc_seg = bsc_charge(start_idx:end_idx);
                P_seg   = dc_power(start_idx:end_idx);
                
                duration_sec = seconds(t_seg(end) - t_seg(1));
                deltaI = max(I_seg) - min(I_seg);
                I_max = max(abs(I_seg));
                
                fprintf('Event %d: Duration=%.1fs, Charge_Duration=%ds, deltaI=%.2fA, I_max=%.2fA\n', i, duration_sec, charge_duration, deltaI, I_max);
                
                % idx2 is the last point of "Charging"
                idx2 = charge_end_idx - start_idx + 1;  % Relative position within segment
                
                % Start point of Charging period (relative position within segment)
                charge_start_in_segment = start_charge_idx - start_idx + 1;

                % Condition #3: Output stability (check only Charging period)
                power_std = std(P_seg(3:idx2));
                current_std = std(I_seg(3:idx2));
                
                fprintf('  Event %d: power_std=%.1fkW (limit=%.1fkW)\n', i, power_std, max_P_std);

                if power_std >= max_P_std || current_std >= max_I_std
                    fprintf('    → FILTERED: Power & Current instability (std=%.1fkW > %.1fkW) | (std=%.1fA > %.1fA)\n', power_std, max_P_std, current_std, max_I_std);
                    continue;    
                end
                step4_count = step4_count + 1;

%% Step 4: DCIR Calculation
                fprintf('=====Step 4: DCIR Calculation =====\n');
                % Save Events as Struct
                eventCount = eventCount + 1;
                evtName = sprintf('event%d', eventCount);
                eventStruct.(evtName).start_idx = start_idx;
                eventStruct.(evtName).end_idx = end_idx;
                eventStruct.(evtName).idx1 = idx1;  % Last point of Idle
                eventStruct.(evtName).idx2 = idx2;  % Last point of Charging (relative position within segment)
                eventStruct.(evtName).charge_duration = charge_duration;

                eventStruct.(evtName).t_seq   = t_seg;
                eventStruct.(evtName).I_seq   = I_seg;
                eventStruct.(evtName).V_seq   = V_seg;
                eventStruct.(evtName).T_seq   = T_seg;
                eventStruct.(evtName).soc_seq = soc_seg;
                eventStruct.(evtName).bsc_seq = bsc_seg;
                eventStruct.(evtName).P_seq   = P_seg;

                eventStruct.(evtName).t_seg_ridIdle = t_seg(1:idx2);
                eventStruct.(evtName).I_seg_ridIdle = I_seg(1:idx2);
                eventStruct.(evtName).V_seg_ridIdle = V_seg(1:idx2);
                eventStruct.(evtName).T_seg_ridIdle = T_seg(1:idx2);
                eventStruct.(evtName).soc_seq_ridIdle = soc_seg(1:idx2);
                eventStruct.(evtName).P_seg_ridIdle = P_seg(1:idx2);

                % DCIR calculation: time-based
                dt_list = [1, 3, 5, 10, 30, 50];
                for d = 1:length(dt_list)
                    dt_sec = dt_list(d);
                    idx_dt = 1 + dt_sec;

                    if idx_dt <= length(I_seg)
                        V1 = V_seg(1);
                        V2 = V_seg(idx_dt);
                        I1 = I_seg(1);
                        I2 = I_seg(idx_dt);
                        dV = V2 - V1;
                        dI = I2 - I1;

                        if (dI) > 0 && dV > 0
                            dcir_val = (dV / dI) * 1000; % [mOhm]
                        else
                            dcir_val = NaN;
                        end
                    else
                        V1 = NaN; V2 = NaN; I1 = NaN; I2 = NaN;
                        dcir_val = NaN;
                    end

                    fieldName = sprintf('DCIR_%ds', dt_sec);

                    eventStruct.(evtName).(fieldName).val = dcir_val;
                    eventStruct.(evtName).(fieldName).V1 = V1;
                    eventStruct.(evtName).(fieldName).V2 = V2;
                    eventStruct.(evtName).(fieldName).I1 = I1;
                    eventStruct.(evtName).(fieldName).I2 = I2;
                    eventStruct.(evtName).(fieldName).dV = dV;
                    eventStruct.(evtName).(fieldName).dI = dI;
                end              

                % DCIR calculation based on Charging end point
                if idx2 <= length(I_seg)
                    V1_r = V_seg(1);
                    V2_r = V_seg(idx2);
                    I1_r = I_seg(1);
                    I2_r = I_seg(idx2);
                    dV_ramp = V2_r - V1_r;
                    dI_ramp = I2_r - I1_r;

                    if (dV_ramp) > 0 && (dI_ramp) > 0
                        dcir_ramp = (dV_ramp / dI_ramp) * 1000; % [mOhm]
                    else
                        dcir_ramp = NaN;
                    end
                else
                    V1_r = NaN; V2_r = NaN; I1_r = NaN; I2_r = NaN;
                    dcir_ramp = NaN;
                end

                eventStruct.(evtName).DCIR_ramp.val = dcir_ramp;
                eventStruct.(evtName).DCIR_ramp.V1  = V1_r;
                eventStruct.(evtName).DCIR_ramp.V2  = V2_r;
                eventStruct.(evtName).DCIR_ramp.I1  = I1_r;
                eventStruct.(evtName).DCIR_ramp.I2  = I2_r;
                eventStruct.(evtName).DCIR_ramp.dV_ramp  = dV_ramp;
                eventStruct.(evtName).DCIR_ramp.dI_ramp  = dI_ramp;

%% Step 5: DCIR Difference Calculation
                fprintf('=====Step 5: DCIR Difference Calculation Ts-1s =====\n');
                % Additional DCIR difference calculation
                if isfield(eventStruct.(evtName), 'DCIR_5s') && isfield(eventStruct.(evtName), 'DCIR_1s')
                    if ~isnan(eventStruct.(evtName).DCIR_5s.val) && ~isnan(eventStruct.(evtName).DCIR_1s.val)
                        eventStruct.(evtName).DCIR_diff_5s_1s = eventStruct.(evtName).DCIR_5s.val - eventStruct.(evtName).DCIR_1s.val;
                    else
                        eventStruct.(evtName).DCIR_diff_5s_1s = NaN;
                    end
                end
                
                if isfield(eventStruct.(evtName), 'DCIR_10s') && isfield(eventStruct.(evtName), 'DCIR_1s')
                    if ~isnan(eventStruct.(evtName).DCIR_10s.val) && ~isnan(eventStruct.(evtName).DCIR_1s.val)
                        eventStruct.(evtName).DCIR_diff_10s_1s = eventStruct.(evtName).DCIR_10s.val - eventStruct.(evtName).DCIR_1s.val;
                    else
                        eventStruct.(evtName).DCIR_diff_10s_1s = NaN;
                    end
                end
                
                if isfield(eventStruct.(evtName), 'DCIR_30s') && isfield(eventStruct.(evtName), 'DCIR_1s')
                    if ~isnan(eventStruct.(evtName).DCIR_30s.val) && ~isnan(eventStruct.(evtName).DCIR_1s.val)
                        eventStruct.(evtName).DCIR_diff_30s_1s = eventStruct.(evtName).DCIR_30s.val - eventStruct.(evtName).DCIR_1s.val;
                    else
                        eventStruct.(evtName).DCIR_diff_30s_1s = NaN;
                    end
                end
                
                if isfield(eventStruct.(evtName), 'DCIR_50s') && isfield(eventStruct.(evtName), 'DCIR_1s')
                    if ~isnan(eventStruct.(evtName).DCIR_50s.val) && ~isnan(eventStruct.(evtName).DCIR_1s.val)
                        eventStruct.(evtName).DCIR_diff_50s_1s = eventStruct.(evtName).DCIR_50s.val - eventStruct.(evtName).DCIR_1s.val;
                    else
                        eventStruct.(evtName).DCIR_diff_50s_1s = NaN;
                    end
                end

                filtered_events = filtered_events + 1;

                % Must add event_number field
                eventStruct.(evtName).event_number = eventCount;
            end
            
            % File-by-file filtering result summary
            fprintf('  Step summary:\n');
            fprintf('    After charging end detection: %d\n', step2_count);
            fprintf('    After duration filter (≥%ds): %d\n', min_charge_duration, step3_count);
            fprintf('    After power stability check: %d\n', step4_count);
            fprintf('    Final valid events from this file: %d\n', step4_count);
        end
    end
end

save(fullfile(saveDir, 'all_chg_events_current_clustering_2021.mat'), 'eventStruct');
fprintf('Step 1 complete. Events saved.\n');
fprintf('Total potential events detected: %d\n', total_events);
fprintf('Events that passed all filters: %d\n', filtered_events);
fprintf('Final stored events: %d\n', eventCount);
fprintf('\n');

%% Step 6: Current-based Clustering and Visualization
fprintf('=====Step 6: Current-based Clustering =====\n');
% Event classification by current range
current_clusters = struct();
for i = 1:length(current_labels)
    current_clusters.(current_labels{i}) = [];
end

eventNames = fieldnames(eventStruct);
fprintf('Step 4: Current-based clustering analysis...\n');

for i = 1:length(eventNames)
    evt = eventStruct.(eventNames{i});
    
    % Calculate average charging current
    if isfield(evt, 'I_seg_ridIdle') && ~isempty(evt.I_seg_ridIdle)
        avg_current = mean(evt.I_seg_ridIdle(3:end));  % Average current in charging period
        
        % Determine current range
        current_bin_idx = find(avg_current >= current_bins, 1, 'last');
        if ~isempty(current_bin_idx) && current_bin_idx < length(current_bins)
            current_range = current_labels{current_bin_idx};
            % Save event_number
            current_clusters.(current_range) = [current_clusters.(current_range), evt.event_number];
            
            % Must save as current_label
            eventStruct.(eventNames{i}).avg_current = avg_current;
            eventStruct.(eventNames{i}).current_label = current_range;
        end
    end
end

% Statistics output by range
fprintf('\nCurrent-based clustering results:\n');
for i = 1:length(current_labels)
    cluster_size = length(current_clusters.(current_labels{i}));
    fprintf('  %s: %d events\n', current_labels{i}, cluster_size);
end

% Added after Step 6
% Convert eventStruct to current label → event number → variables structure

eventStruct_new = struct();
original_eventNames = fieldnames(eventStruct);
for i = 1:length(original_eventNames)
    evt = eventStruct.(original_eventNames{i});
    if ~isfield(evt, 'current_label')
        continue;
    end
    label = evt.current_label;      % e.g., 'range_50_100A'
    eventN = sprintf('event%d', evt.event_number); % e.g., 'event13'
    eventStruct_new.(label).(eventN) = evt;
end

eventStruct = eventStruct_new;

%% Step 7: Current-based Visualization
fprintf('=====Step 7: Current-based Visualization =====\n');

% Display labels for plots (with proper formatting)
display_labels = {'50-100A', '100-150A', '150-200A', '200-250A', '250-300A'};

% 1. Current vs Time and Voltage vs Time graphs by Current label
figure('Name', 'Current vs Time and Voltage vs Time by Current Range', 'Position', [50, 50, 1200, 800]);

for i = 1:length(current_labels)
    label = current_labels{i};
    cluster_events = current_clusters.(label);
    if isempty(cluster_events)
        continue;
    end
    
    num_events_to_show = length(cluster_events);
    subplot(2, 3, i); hold on;
    
    for j = 1:num_events_to_show
        eventN = sprintf('event%d', cluster_events(j));
        evt = eventStruct.(label).(eventN);
        if isfield(evt, 't_seg_ridIdle') && isfield(evt, 'I_seg_ridIdle')
            t_data = evt.t_seg_ridIdle;
            I_data = evt.I_seg_ridIdle;
            t_seconds = seconds(t_data - t_data(1));
            plot(t_seconds, I_data, 'LineWidth', 1.5);
        end
    end
    title(sprintf('Current vs Time: %s (%d events)', display_labels{i}, length(cluster_events)));
    xlabel('Time [s]');
    ylabel('Current [A]');
    grid on;
    % legend('Location', 'best');
end

sgtitle('Current vs Time by Current Range (2021)');
saveas(gcf, fullfile(saveDir, 'fig_2021_Current_vs_Time_by_range.fig'));

% Voltage vs Time graph
figure('Name', 'Voltage vs Time by Current Range', 'Position', [100, 100, 1200, 800]);

for i = 1:length(current_labels)
    label = current_labels{i};
    cluster_events = current_clusters.(label);
    if isempty(cluster_events)
        continue;
    end
    num_events_to_show = length(cluster_events);
    subplot(2, 3, i); hold on;
    for j = 1:num_events_to_show
        eventN = sprintf('event%d', cluster_events(j));
        evt = eventStruct.(label).(eventN);
        if isfield(evt, 't_seg_ridIdle') && isfield(evt, 'V_seg_ridIdle')
            t_data = evt.t_seg_ridIdle;
            V_data = evt.V_seg_ridIdle;
            t_seconds = seconds(t_data - t_data(1));
            plot(t_seconds, V_data, 'LineWidth', 1.5);
        end
    end
    title(sprintf('Voltage vs Time: %s (%d events)', display_labels{i}, length(cluster_events)));
    xlabel('Time [s]');
    ylabel('Voltage [V]');
    grid on;
    % legend('Location', 'best');
end

sgtitle('Voltage vs Time by Current Range (2021)');
saveas(gcf, fullfile(saveDir, 'fig_2021_Voltage_vs_Time_by_range.fig'));

% DCIR fields and labels
% Replace all existing DCIR analysis visualization parts

dcir_fields = {'DCIR_1s', 'DCIR_3s', 'DCIR_5s', 'DCIR_10s', 'DCIR_30s', 'DCIR_50s'};
dcir_labels = {'1s', '3s', '5s', '10s', '30s', '50s'};
dcir_diff_fields = {'DCIR_diff_3s_1s', 'DCIR_diff_5s_1s', 'DCIR_diff_10s_1s', 'DCIR_diff_30s_1s', 'DCIR_diff_50s_1s'};
dcir_diff_labels = {'3s-1s', '5s-1s', '10s-1s', '30s-1s', '50s-1s'};

for i = 1:length(current_labels)
    label = current_labels{i};
    cluster_events = current_clusters.(label);
    if isempty(cluster_events)
        continue;
    end
    figure('Name', sprintf('DCIR & Diff - %s', display_labels{i}), 'Position', [200, 200, 1400, 800]);
    % 6+5=11 subplot, 3 rows 4 columns, last cell empty
    for k = 1:length(dcir_fields)
        subplot(3, 4, k); hold on;
        vals = nan(1, length(cluster_events));
        event_numbers = nan(1, length(cluster_events));
        for j = 1:length(cluster_events)
            eventN = sprintf('event%d', cluster_events(j));
            evt = eventStruct.(label).(eventN);
            if isfield(evt, dcir_fields{k}) && isfield(evt.(dcir_fields{k}), 'val')
                vals(j) = evt.(dcir_fields{k}).val;
                event_numbers(j) = cluster_events(j);  % Use actual event number
            end
        end
        m = nanmean(vals);
        s = nanstd(vals);
        % Outlier detection
        outlier_idx = find(abs(vals - m) > s);
        non_outlier_idx = find(abs(vals - m) <= s);
        
        % Plot non-outlier points with circles
        if ~isempty(non_outlier_idx)
            plot(event_numbers(non_outlier_idx), vals(non_outlier_idx), 'o', 'Color', '#0072BD', 'MarkerSize', 5, 'DisplayName', dcir_labels{k});
        end
        
        % Plot outlier points with X markers
        if ~isempty(outlier_idx)
            plot(event_numbers(outlier_idx), vals(outlier_idx), 'x', 'Color', '#D95319', 'MarkerSize', 5, 'DisplayName', 'Outliers');
        end
        
        yline(m, '--g', 'LineWidth', 2, 'DisplayName', sprintf('Mean: %.2f', m));
        yline(m+s, ':r', 'LineWidth', 2, 'DisplayName', sprintf('+1σ: %.2f', s));
        yline(m-s, ':r', 'LineWidth', 2, 'DisplayName', sprintf('-1σ: %.2f', s));
        title(sprintf('DCIR %s (%d/%d)', dcir_labels{k}, length(outlier_idx), length(cluster_events)));
        xlabel('Event Number');
        ylim([0 40]);
        ylabel('[mΩ]');
        grid on;
        legend('Location', 'northeast');
    end
    for k = 1:length(dcir_diff_fields)
        subplot(3, 4, 6+k); hold on;
        vals = nan(1, length(cluster_events));
        event_numbers = nan(1, length(cluster_events));
        for j = 1:length(cluster_events)
            eventN = sprintf('event%d', cluster_events(j));
            evt = eventStruct.(label).(eventN);
            if isfield(evt, dcir_diff_fields{k})
                vals(j) = evt.(dcir_diff_fields{k});
                event_numbers(j) = cluster_events(j);  % Use actual event number
            end
        end
        m = nanmean(vals);
        s = nanstd(vals);
        % Outlier detection
        outlier_idx = find(abs(vals - m) > s);
        non_outlier_idx = find(abs(vals - m) <= s);
        
        % Plot non-outlier points with circles
        if ~isempty(non_outlier_idx)
            plot(event_numbers(non_outlier_idx), vals(non_outlier_idx), 'o', 'Color', '#0072BD', 'MarkerSize', 5, 'DisplayName', dcir_diff_labels{k});
        end
        
        % Plot outlier points with X markers
        if ~isempty(outlier_idx)
            plot(event_numbers(outlier_idx), vals(outlier_idx), 'x', 'Color', '#D95319', 'MarkerSize', 5, 'DisplayName', 'Outliers');
        end
        
        yline(m, '--g', 'LineWidth', 2, 'DisplayName', sprintf('Mean: %.2f', m));
        yline(m+s, ':r', 'LineWidth', 2, 'DisplayName', sprintf('+1σ: %.2f', s));
        yline(m-s, ':r', 'LineWidth', 2, 'DisplayName', sprintf('-1σ: %.2f', s));
        title(sprintf('DCIR Diff %s (%d/%d)', dcir_diff_labels{k}, length(outlier_idx), length(cluster_events)));
        xlabel('Event Number');
        ylabel('[mΩ]');
        ylim([0 40]);
        grid on;
        legend('Location', 'northeast');
    end
    sgtitle(sprintf('DCIR & Diff by Event (%s, %d events)', display_labels{i}, length(cluster_events)));
    saveas(gcf, fullfile(saveDir, sprintf('fig_2021_DCIR_Diff_vs_Event_%s.fig', label)));
end

fprintf('\n=== Event Preprocessing Completed ===\n'); 
