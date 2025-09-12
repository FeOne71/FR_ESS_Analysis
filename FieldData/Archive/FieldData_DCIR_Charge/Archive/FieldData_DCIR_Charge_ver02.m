%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field Data DCIR Charge Current Clustering (Automatic K-means)
% ESS charging events and resistance analysis (automatic current-based clustering)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory
dataDir  = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\raw2mat_ver04';
yearList = {'2022'}; 
saveDir  = fullfile('G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\FieldData_DCIR_Charge\AutoResults');

if ~exist(saveDir, 'dir')
    mkdir(saveDir); 
end

%% Variables
C_nom        = 1024;          % Ah
min_charge_duration = 30;     % [s] - Charging duration (30 seconds or more)
max_P_std    = 9;            % Max power standard deviation [kW]
max_I_std    = C_nom *0.01;   % Max current standard deviation [A] 10.24A

% Automatic clustering settings
min_cluster_size = 3;          % Minimum number of events per cluster
max_clusters = 5;              % Maximum number of clusters to try

%% Step 1: Raw Files Traversal
fprintf('Loading field data...\n');
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
            load(matFilePath); 
            
            t          = Raw.Time;
            I          = Raw.DCCurrent_A;
            V          = Raw.Total_AverageCVSum_V;
            soc        = Raw.Total_AverageSOC;
            T_batt     = Raw.Total_AverageMT_degC;
            bsc_charge = Raw.Charge;
            dc_power   = Raw.DCPower_kW;
            
%% Step 2: Searching Idle to Charging Transition
            is_idle = strcmp(bsc_charge, 'Idle');
            is_charging = strcmp(bsc_charge, 'Charging');
            
            idle_to_charge = find(is_idle(1:end-1) & is_charging(2:end));
            total_events = total_events + length(idle_to_charge);

%% Step 3: Events Filtering              
            for i = 1:length(idle_to_charge)
                idx1 = idle_to_charge(i);
                start_charge_idx = idx1 + 1;
                
                % Find end point of "Charging" period
                charge_end_idx = start_charge_idx;
                while charge_end_idx <= length(bsc_charge) && strcmp(bsc_charge(charge_end_idx), 'Charging')
                    charge_end_idx = charge_end_idx + 1;
                end
                charge_end_idx = charge_end_idx - 1;
                
                start_idx = idx1;
                end_idx = charge_end_idx;
                
                % Check charging duration
                charge_duration = charge_end_idx - start_charge_idx + 1;
                if charge_duration < min_charge_duration
                    continue;
                end
                
                t_seg   = t(start_idx:end_idx);
                I_seg   = I(start_idx:end_idx);
                V_seg   = V(start_idx:end_idx);
                soc_seg = soc(start_idx:end_idx);
                T_seg   = T_batt(start_idx:end_idx);
                bsc_seg = bsc_charge(start_idx:end_idx);
                P_seg   = dc_power(start_idx:end_idx);
                
                idx2 = charge_end_idx - start_idx + 1;

                % Check output stability
                power_std = std(P_seg(3:idx2));
                current_std = std(I_seg(3:idx2));

                if power_std >= max_P_std || current_std >= max_I_std
                    continue;    
                end

%% Step 4: DCIR Calculation
                eventCount = eventCount + 1;
                evtName = sprintf('event%d', eventCount);
                eventStruct.(evtName).start_idx = start_idx;
                eventStruct.(evtName).end_idx = end_idx;
                eventStruct.(evtName).idx1 = idx1;
                eventStruct.(evtName).idx2 = idx2;
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
                            dcir_val = (dV / dI) * 1000;
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
                        dcir_ramp = (dV_ramp / dI_ramp) * 1000;
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
                eventStruct.(evtName).event_number = eventCount;
            end
        end
    end
end

save(fullfile(saveDir, 'all_chg_events_current_clustering_auto_2022.mat'), 'eventStruct');
fprintf('Events saved. Total: %d, Filtered: %d\n', total_events, filtered_events);

%% Step 6: Automatic K-means Clustering
fprintf('Performing automatic K-means clustering...\n');

eventNames = fieldnames(eventStruct);
avg_currents = [];
valid_event_indices = [];

for i = 1:length(eventNames)
    evt = eventStruct.(eventNames{i});
    
    if isfield(evt, 'I_seg_ridIdle') && ~isempty(evt.I_seg_ridIdle)
        avg_current = mean(evt.I_seg_ridIdle(3:end));
        avg_currents = [avg_currents, avg_current];
        valid_event_indices = [valid_event_indices, i];
        eventStruct.(eventNames{i}).avg_current = avg_current;
    end
end

if length(avg_currents) >= min_cluster_size
    silhouette_scores = [];
    cluster_results = {};
    
    for k = 2:min(max_clusters, length(avg_currents)-1)
        [idx, centroids] = kmeans(avg_currents', k, 'Replicates', 10, 'Distance', 'sqeuclidean');
        silhouette_avg = mean(silhouette(avg_currents', idx));
        silhouette_scores = [silhouette_scores, silhouette_avg];
        cluster_results{k} = {idx, centroids};
    end
    
    [max_score, best_k_idx] = max(silhouette_scores);
    best_k = best_k_idx + 1;
    
    similar_scores = find(silhouette_scores >= max_score * 0.95);
    if length(similar_scores) > 1
        best_k = similar_scores(end) + 1;
    end
    
    best_idx = cluster_results{best_k}{1};
    best_centroids = cluster_results{best_k}{2};
    
    current_labels = {};
    for i = 1:best_k
        centroid_val = round(best_centroids(i));
        label = sprintf('cluster_%dA', centroid_val);
        current_labels{i} = label;
    end
    
    current_clusters = struct();
    for i = 1:length(current_labels)
        current_clusters.(current_labels{i}) = [];
    end
    
    for i = 1:length(valid_event_indices)
        event_idx = valid_event_indices(i);
        cluster_idx = best_idx(i);
        event_name = eventNames{event_idx};
        evt = eventStruct.(event_name);
        
        current_label = current_labels{cluster_idx};
        eventStruct.(event_name).current_label = current_label;
        current_clusters.(current_label) = [current_clusters.(current_label), evt.event_number];
    end
    
else
    current_labels = {'cluster_all'};
    current_clusters = struct();
    current_clusters.cluster_all = [];
    
    for i = 1:length(valid_event_indices)
        event_idx = valid_event_indices(i);
        event_name = eventNames{event_idx};
        evt = eventStruct.(event_name);
        eventStruct.(event_name).current_label = 'cluster_all';
        current_clusters.cluster_all = [current_clusters.cluster_all, evt.event_number];
    end
end

eventStruct_new = struct();
original_eventNames = fieldnames(eventStruct);
for i = 1:length(original_eventNames)
    evt = eventStruct.(original_eventNames{i});
    if ~isfield(evt, 'current_label')
        continue;
    end
    label = evt.current_label;
    eventN = sprintf('event%d', evt.event_number);
    eventStruct_new.(label).(eventN) = evt;
end

eventStruct = eventStruct_new;

%% Step 7: Visualization
fprintf('Generating visualizations...\n');

display_labels = {};
for i = 1:length(current_labels)
    label = current_labels{i};
    if contains(label, 'cluster_')
        current_val = regexp(label, 'cluster_(\d+)A', 'tokens');
        if ~isempty(current_val)
            current_num = str2double(current_val{1}{1});
            c_rate = current_num / C_nom;
            display_labels{i} = sprintf('%.2fC', c_rate);
        else
            display_labels{i} = label;
        end
    else
        display_labels{i} = label;
    end
end

% Current vs Time graphs
figure('Name', 'Current vs Time by Current Range', 'Position', [50, 50, 1200, 800]);

for i = 1:length(current_labels)
    label = current_labels{i};
    cluster_events = current_clusters.(label);
    if isempty(cluster_events)
        continue;
    end
    
    subplot(2, 3, i); hold on;
    
    for j = 1:length(cluster_events)
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
end

sgtitle('Current vs Time by Current Range (2022)');
saveas(gcf, fullfile(saveDir, 'fig_2022_Current_vs_Time_by_range_auto.fig'));

% Voltage vs Time graph
figure('Name', 'Voltage vs Time by Current Range', 'Position', [100, 100, 1200, 800]);

for i = 1:length(current_labels)
    label = current_labels{i};
    cluster_events = current_clusters.(label);
    if isempty(cluster_events)
        continue;
    end
    subplot(2, 3, i); hold on;
    for j = 1:length(cluster_events)
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
end

sgtitle('Voltage vs Time by Current Range (2022)');
saveas(gcf, fullfile(saveDir, 'fig_2022_Voltage_vs_Time_by_range_auto.fig'));

% DCIR Analysis
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
    
    display_idx = find(strcmp(current_labels, label));
    if isempty(display_idx)
        display_label = label;
    else
        display_label = display_labels{display_idx};
    end
    
    figure('Name', sprintf('DCIR & Diff - %s', display_label), 'Position', [200, 200, 1400, 800]);
    
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
        outlier_idx = find(abs(vals - m) > s);
        non_outlier_idx = find(abs(vals - m) <= s);
        
        if ~isempty(non_outlier_idx)
            plot(event_numbers(non_outlier_idx), vals(non_outlier_idx), 'o', 'Color', '#0072BD', 'MarkerSize', 5, 'DisplayName', dcir_labels{k});
        end
        
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
        outlier_idx = find(abs(vals - m) > s);
        non_outlier_idx = find(abs(vals - m) <= s);
        
        if ~isempty(non_outlier_idx)
            plot(event_numbers(non_outlier_idx), vals(non_outlier_idx), 'o', 'Color', '#0072BD', 'MarkerSize', 5, 'DisplayName', dcir_diff_labels{k});
        end
        
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
    sgtitle(sprintf('DCIR & Diff by Event (%s, %d events)', display_label, length(cluster_events)));
    saveas(gcf, fullfile(saveDir, sprintf('fig_2022_DCIR_Diff_vs_Event_%s_auto.fig', label)));
end

fprintf('Analysis completed successfully.\n'); 
fprintf('\n=== Event Preprocessing Completed ===\n'); 