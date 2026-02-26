
% RPT Feature Dataset Generation (Micro-Segmentation)
% Charge: 3.7-3.95V (step 0.05), Discharge: 3.75-3.87V (step 0.03)
% Purpose: Generate Separate Charge/Discharge Datasets for ML (Feature -> SOH)

clear; clc; close all;
warning off;

%% 1. Configuration
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing';
crateFile = fullfile(baseDir, 'Crate_integrated\Crate_integrated.mat');
% Unified output: Feature_Dataset/Dataset/ (mat files)
resultRoot = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\Feature_Dataset';
outputDir = fullfile(resultRoot, 'Dataset');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

% Segment definition: Charge vs Discharge (separate datasets)
% Charge: 3.7 - 3.95 V, step 0.05 -> 5 segments
win_min_chg = 3.7;
win_max_chg = 3.95;
win_step_chg = 0.05;
segments_chg = win_min_chg:win_step_chg:win_max_chg;

% Discharge: 3.75 - 3.87 V, step 0.03 -> 4 segments
win_min_dchg = 3.75;
win_max_dchg = 3.87;
win_step_dchg = 0.03;
segments_dchg = win_min_dchg:win_step_dchg:win_max_dchg; 

% C-rates to include
target_crates = [0.1, 0.5, 1, 2, 3];
crate_keys = {'c01', 'c05', 'c1', 'c2', 'c3'};
channels = {'Ch09','Ch10','Ch11','Ch12','Ch13','Ch14','Ch15','Ch16'}; 

fprintf('Loading Data: %s...\n', crateFile);
load(crateFile, 'Crate_data');

%% 2. Feature Extraction
fprintf('\n=== Charge: %.2f-%.2fV (step %.2f) | Discharge: %.2f-%.2fV (step %.2f) ===\n', ...
    win_min_chg, win_max_chg, win_step_chg, win_min_dchg, win_max_dchg, win_step_dchg);

% Dataset Containers (separate for Charge / Discharge)
X_data_chg = [];
Y_data_chg = [];
Meta_data_chg = table();
summary_time_matrix_chg = [];

X_data_dchg = [];
Y_data_dchg = [];
Meta_data_dchg = table();
summary_time_matrix_dchg = [];

count_chg = 0;
count_dchg = 0;

for chIdx = 1:length(channels)
    chName = channels{chIdx};
    if ~isfield(Crate_data, chName), continue; end
    cycles = fieldnames(Crate_data.(chName));
    
    for cycIdx = 1:length(cycles)
        cycleName = cycles{cycIdx};
        cycNum = str2double(strrep(cycleName, 'cyc', ''));
        
        for k = 1:length(crate_keys)
            cKey = crate_keys{k};
            crateVal = target_crates(k);
            
            if isfield(Crate_data.(chName).(cycleName), cKey)
            % Process both Charge and Discharge
            stepTypes = {'charge', 'discharge'};
            
            for sType = 1:length(stepTypes)
                typeKey = stepTypes{sType};
                isCharge = strcmp(typeKey, 'charge');
                
                if isfield(Crate_data.(chName).(cycleName).(cKey), typeKey)
                    rawData = Crate_data.(chName).(cycleName).(cKey).(typeKey);
                else
                    continue;
                end
                
                if isempty(rawData.V) || isempty(rawData.Q), continue; end
                
                V = rawData.V;
                Q = rawData.Q;
                
                if isfield(rawData, 'I') && ~isempty(rawData.I)
                    I_avg = mean(abs(rawData.I));
                else
                    I_avg = NaN;
                end
                
                % Segment definition per type
                if isCharge
                    segments = segments_chg;
                    win_min = win_min_chg;
                    win_max = win_max_chg;
                else
                    segments = segments_dchg;
                    win_min = win_min_dchg;
                    win_max = win_max_dchg;
                end
                
                % === Preprocessing: Standardize V-Q to Common Grid ===
                v_global_min = 3.0;
                v_global_max = 4.2;
                n_grid_points = 1000;
                v_common = linspace(v_global_min, v_global_max, n_grid_points);
                
                [V_u, uid] = unique(V);
                Q_u = Q(uid);
                
                try
                    Q_common = interp1(V_u, Q_u, v_common, 'linear');
                catch
                    continue;
                end
                
                % === Feature 1: Voltage Segment Capacity & Time ===
                n_seg = length(segments) - 1;
                feat_vec_cap = zeros(1, n_seg);
                feat_vec_time = zeros(1, n_seg);
                
                for s = 1:n_seg
                    v_start = segments(s);
                    v_end = segments(s+1);
                    
                    idx_start = find(v_common >= v_start, 1, 'first');
                    idx_end = find(v_common <= v_end, 1, 'last');
                    
                    if ~isempty(idx_start) && ~isempty(idx_end) && (idx_end > idx_start)
                        dq = abs(Q_common(idx_end) - Q_common(idx_start));
                        feat_vec_cap(s) = dq;
                        if ~isnan(I_avg) && I_avg > 0
                            dt = (dq / I_avg) * 3600;
                        else
                            dt = NaN;
                        end
                        feat_vec_time(s) = dt;
                    else
                        feat_vec_cap(s) = NaN;
                        feat_vec_time(s) = NaN;
                    end
                end
                
                % === Feature 2: Differential Capacity (window per type) ===
                dv_step = v_common(2) - v_common(1);
                dq_common = gradient(Q_common);
                dqdv_common = abs(dq_common / dv_step);
                
                win_indices = v_common >= win_min & v_common <= win_max;
                dqdv_window = dqdv_common(win_indices);
                v_window = v_common(win_indices);
                
                if sum(~isnan(dqdv_window)) > 5
                    [peak_height, max_idx] = max(dqdv_window);
                    peak_volt = v_window(max_idx);
                    valid_mask = ~isnan(dqdv_window);
                    peak_area = trapz(v_window(valid_mask), dqdv_window(valid_mask));
                else
                    peak_height = NaN;
                    peak_volt = NaN;
                    peak_area = NaN;
                end
                
                soh = max(Q) - min(Q);
                total_Q = max(Q) - min(Q);
                if ~isnan(I_avg)
                    total_cc_time = (total_Q / I_avg) * 3600;
                else
                    total_cc_time = NaN;
                end
                
                if ~any(isnan(feat_vec_cap)) && ~isnan(peak_height) && ~isnan(total_cc_time)
                    row = [feat_vec_cap, feat_vec_time, peak_height, peak_volt, peak_area, crateVal];
                    
                    if isCharge
                        X_data_chg = [X_data_chg; row];
                        Y_data_chg = [Y_data_chg; soh];
                        Meta_table_row = cell2table({chName, cycleName, crateVal, soh}, 'VariableNames', {'Channel', 'Cycle', 'Crate', 'SOH'});
                        Meta_data_chg = [Meta_data_chg; Meta_table_row];
                        summary_time_matrix_chg = [summary_time_matrix_chg; [feat_vec_time, total_cc_time, crateVal]];
                        count_chg = count_chg + 1;
                    else
                        X_data_dchg = [X_data_dchg; row];
                        Y_data_dchg = [Y_data_dchg; soh];
                        Meta_table_row = cell2table({chName, cycleName, crateVal, soh}, 'VariableNames', {'Channel', 'Cycle', 'Crate', 'SOH'});
                        Meta_data_dchg = [Meta_data_dchg; Meta_table_row];
                        summary_time_matrix_dchg = [summary_time_matrix_dchg; [feat_vec_time, total_cc_time, crateVal]];
                        count_dchg = count_dchg + 1;
                    end
                end
            end % End stepTypes loop
            end
        end
    end
end

%% 3. Generate Statistics Table (Per C-rate, per dataset)
fprintf('Extraction Complete: Charge %d samples, Discharge %d samples.\n', count_chg, count_dchg);
fprintf('Charge Feature Count: %d (5 Cap + 5 Time + 3 Peak + Crate)\n', size(X_data_chg, 2));
fprintf('Discharge Feature Count: %d (4 Cap + 4 Time + 3 Peak + Crate)\n', size(X_data_dchg, 2));

% Summary: Charge
if ~isempty(summary_time_matrix_chg)
    num_segs = length(segments_chg) - 1;
    unique_crates = unique(summary_time_matrix_chg(:, end));
    fprintf('\n====== Summary: Charge Steps (%.2f-%.2fV) ======\n', win_min_chg, win_max_chg);
    for uc = 1:length(unique_crates)
        target_c = unique_crates(uc);
        c_mask = summary_time_matrix_chg(:, end) == target_c;
        sub_data = summary_time_matrix_chg(c_mask, :);
        if isempty(sub_data), continue; end
        mean_times = mean(sub_data(:, 1:num_segs));
        median_times = median(sub_data(:, 1:num_segs));
        total_times = sub_data(:, num_segs+1);
        valid_totals = total_times; valid_totals(valid_totals==0) = NaN;
        props = sub_data(:, 1:num_segs) ./ valid_totals * 100;
        mean_props = mean(props, 'omitnan'); median_props = median(props, 'omitnan');
        fprintf('\n--- %.1f C Charge (N=%d) ---\n', target_c, sum(c_mask));
        fprintf('%-14s | %-12s | %-12s | %-12s | %-12s\n', 'Segment', 'Avg Time(s)', 'Avg Prop(%)', 'Med Time(s)', 'Med Prop(%)');
        fprintf('Total CC | %-12.1f | - | %-12.1f | -\n', mean(total_times), median(total_times));
        for i = 1:num_segs
            fprintf('Seg %d (%.2fV) | %-12.1f | %-12.1f | %-12.1f | %-12.1f\n', i, segments_chg(i), mean_times(i), mean_props(i), median_times(i), median_props(i));
        end
        fprintf('================================================================================\n');
    end
end

% Summary: Discharge
if ~isempty(summary_time_matrix_dchg)
    num_segs = length(segments_dchg) - 1;
    unique_crates = unique(summary_time_matrix_dchg(:, end));
    fprintf('\n====== Summary: Discharge Steps (%.2f-%.2fV) ======\n', win_min_dchg, win_max_dchg);
    for uc = 1:length(unique_crates)
        target_c = unique_crates(uc);
        c_mask = summary_time_matrix_dchg(:, end) == target_c;
        sub_data = summary_time_matrix_dchg(c_mask, :);
        if isempty(sub_data), continue; end
        mean_times = mean(sub_data(:, 1:num_segs));
        median_times = median(sub_data(:, 1:num_segs));
        total_times = sub_data(:, num_segs+1);
        valid_totals = total_times; valid_totals(valid_totals==0) = NaN;
        props = sub_data(:, 1:num_segs) ./ valid_totals * 100;
        mean_props = mean(props, 'omitnan'); median_props = median(props, 'omitnan');
        fprintf('\n--- %.1f C Discharge (N=%d) ---\n', target_c, sum(c_mask));
        fprintf('%-14s | %-12s | %-12s | %-12s | %-12s\n', 'Segment', 'Avg Time(s)', 'Avg Prop(%)', 'Med Time(s)', 'Med Prop(%)');
        fprintf('Total CC | %-12.1f | - | %-12.1f | -\n', mean(total_times), median(total_times));
        for i = 1:num_segs
            fprintf('Seg %d (%.2fV) | %-12.1f | %-12.1f | %-12.1f | %-12.1f\n', i, segments_dchg(i), mean_times(i), mean_props(i), median_times(i), median_props(i));
        end
        fprintf('================================================================================\n');
    end
end

%% 4. Save to MAT (separate Charge / Discharge)
% Charge
feature_names_chg = arrayfun(@(i) sprintf('Cap_%.2f_%.2fV', segments_chg(i), segments_chg(i+1)), 1:length(segments_chg)-1, 'UniformOutput', false);
time_names_chg = arrayfun(@(i) sprintf('Time_%.2f_%.2fV', segments_chg(i), segments_chg(i+1)), 1:length(segments_chg)-1, 'UniformOutput', false);
feature_names_chg = [feature_names_chg, time_names_chg, {'Peak_Height', 'Peak_Volt', 'Peak_Area', 'Crate'}];

savePathChg = fullfile(outputDir, 'Feature_Dataset_Charge.mat');
save(savePathChg, 'X_data_chg', 'Y_data_chg', 'Meta_data_chg', 'feature_names_chg', 'segments_chg');
fprintf('\nSaved: %s\n', savePathChg);

% Discharge
feature_names_dchg = arrayfun(@(i) sprintf('Cap_%.2f_%.2fV', segments_dchg(i), segments_dchg(i+1)), 1:length(segments_dchg)-1, 'UniformOutput', false);
time_names_dchg = arrayfun(@(i) sprintf('Time_%.2f_%.2fV', segments_dchg(i), segments_dchg(i+1)), 1:length(segments_dchg)-1, 'UniformOutput', false);
feature_names_dchg = [feature_names_dchg, time_names_dchg, {'Peak_Height', 'Peak_Volt', 'Peak_Area', 'Crate'}];

savePathDchg = fullfile(outputDir, 'Feature_Dataset_Discharge.mat');
save(savePathDchg, 'X_data_dchg', 'Y_data_dchg', 'Meta_data_dchg', 'feature_names_dchg', 'segments_dchg');
fprintf('Saved: %s\n', savePathDchg);
