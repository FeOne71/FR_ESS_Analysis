%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field Data DCIR Charge Current Clustering (Automatic K-means)
% ESS charging events and resistance analysis (automatic current-based clustering)
% Multi-year processing: 2021, 2022, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory
dataDir  = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
yearList = {'2021', '2022', '2023'}; 
saveDir  = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR\AutoResults_Charge');
if ~exist(saveDir, 'dir')
    mkdir(saveDir); 
end

%% Variables
C_nom        = 1024;          % Ah
min_charge_duration = 300;     % [s] - Charging duration (30 seconds or more)
max_P_std    = 5;            % Max power standard deviation [kW]
max_I_std    = C_nom *0.02;   % Max current standard deviation [A] 10.24A

% Automatic clustering settings
min_cluster_size = 3;          % Minimum number of events per cluster
max_clusters = 5;              % Maximum number of clusters to try

%% Initialize global event structure
global_eventStruct = struct();

%% Process each year
for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    fprintf('Processing year: %s\n', year);
    
    % Initialize year-specific variables
    eventStruct = struct();
    eventCount = 0;
    total_events = 0;
    filtered_events = 0;
    
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

                if current_std >= max_I_std || power_std >= max_P_std
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
    
    fprintf('Year %s: Total: %d, Filtered: %d\n', year, total_events, filtered_events);

%% Step 6: Automatic K-means Clustering for each year
    fprintf('Performing automatic K-means clustering for %s...\n', year);
    
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
            c_rate = centroid_val / C_nom;
            % Convert to C-rate with underscore as decimal point
            if c_rate >= 1
                label = sprintf('cluster_%.0fC', c_rate);
            else
                label = sprintf('cluster_0_%02dC', round(c_rate * 100));
            end
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

    % Convert to cluster-based structure and add to global structure
    eventStruct_new = struct();
    original_eventNames = fieldnames(eventStruct);
    
    % Initialize cluster-specific event counters
    cluster_event_counters = struct();
    
    for i = 1:length(original_eventNames)
        evt = eventStruct.(original_eventNames{i});
        if ~isfield(evt, 'current_label')
            continue;
        end
        label = evt.current_label;
        
        % Initialize counter for this cluster if not exists
        if ~isfield(cluster_event_counters, label)
            cluster_event_counters.(label) = 0;
        end
        
        % Increment cluster-specific counter
        cluster_event_counters.(label) = cluster_event_counters.(label) + 1;
        eventN = sprintf('event%d', cluster_event_counters.(label));
        
        % Add to global structure with year information
        if ~isfield(global_eventStruct, label)
            global_eventStruct.(label) = struct();
        end
        % Use year as string to avoid invalid field name error
        year_str = sprintf('year_%s', year);
        if ~isfield(global_eventStruct.(label), year_str)
            global_eventStruct.(label).(year_str) = struct();
        end
        global_eventStruct.(label).(year_str).(eventN) = evt;
    end
    
    fprintf('Year %s clustering completed.\n', year);
end

%% Save combined results
save(fullfile(saveDir, 'all_chg_events_current_clustering_all_years.mat'), 'global_eventStruct');
fprintf('All years processing completed. Results saved to AutoResults folder.\n');

fprintf('Analysis completed successfully.\n'); 
fprintf('\n=== Event Preprocessing Completed ===\n');
fprintf('Starting visualization...\n');

%% Run visualization script
DCIR_All_years;

%% Run outlier removal and conditional plotting script
fprintf('\n=== Starting Outlier Removal and Conditional Plotting ===\n');
DCIR_Outlier_Removal_and_Conditional_Plots; 