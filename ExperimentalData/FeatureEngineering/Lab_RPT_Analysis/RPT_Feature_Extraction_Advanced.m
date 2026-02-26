%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Feature Extraction Advanced
% Master Ruler Creation (Time-Balanced)
% Feature Extraction (dQ, dQ/dV with Savitzky-Golay)
% Label Generation (SOH, LLI, LAM)
% Matrix Integration & Standardization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
warning off;
%% ========================================================================
% Section 1: Configuration & Data Loading
% ========================================================================
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_static_mat = fullfile(baseDir, 'Capacity_Trend_Figures', 'Capacity_Data_Static.mat');
path_rpt_vq_mat = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\Dataset';

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
fprintf('Loading Data...\n');
load(path_static_mat, 'allChannelsCapacity');
load(path_rpt_vq_mat, 'RPT_VQ_grid');
fprintf('Data Loaded.\n');
win_chg_min = 3.7;
win_chg_max = 3.95;
win_dch_max = 3.88;
win_dch_min = 3.75;
channels = fieldnames(allChannelsCapacity);
%% ========================================================================
% Section 2: Phase 1 - Master Ruler Creation
% ========================================================================
fprintf('\nPhase 1: Generating Global Master Ruler (Static Capacity)...\n');
MasterRulers = struct();
num_segments = 5;

% === Step 1: Collect Static discharge time grids from all 8 channels ===
% Static 방전 V-t 데이터 하나로 충전/방전 모두 결정
standard_V_grid = RPT_VQ_grid.cyc0.(channels{1}).Static.V_grid;

% 충전 윈도우 (3.700 ~ 3.950)
mask_chg = standard_V_grid >= win_chg_min & standard_V_grid <= win_chg_max;
V_standard_chg = standard_V_grid(mask_chg);

% 방전 윈도우 (3.750 ~ 3.880)
mask_dch = standard_V_grid >= win_dch_min & standard_V_grid <= win_dch_max;
V_standard_dch = standard_V_grid(mask_dch);

all_T_grids_chg = zeros(length(channels), length(V_standard_chg));
all_T_grids_dch = zeros(length(channels), length(V_standard_dch));
valid_channels = 0;

for i = 1:length(channels)
    ch = channels{i};
    
    if ~isfield(RPT_VQ_grid.cyc0, ch) || ~isfield(RPT_VQ_grid.cyc0.(ch), 'Static')
        warning('Channel %s: Static data not found. Skipping.', ch);
        continue;
    end
    
    data_s = RPT_VQ_grid.cyc0.(ch).Static;
    V_s = data_s.V_grid;
    t_s = data_s.t;
    valid_channels = valid_channels + 1;
    
    % 충전 윈도우 구간
    mask_c = V_s >= win_chg_min & V_s <= win_chg_max;
    V_c = V_s(mask_c);
    t_c = t_s(mask_c);
    if ~isempty(t_c) && ~isempty(V_c)
        t_c_interp = interp1(V_c, t_c, V_standard_chg, 'linear');
        all_T_grids_chg(valid_channels, :) = t_c_interp';
    end
    
    % 방전 윈도우 구간
    mask_d = V_s >= win_dch_min & V_s <= win_dch_max;
    V_d = V_s(mask_d);
    t_d = t_s(mask_d);
    if ~isempty(t_d) && ~isempty(V_d)
        t_d_interp = interp1(V_d, t_d, V_standard_dch, 'linear');
        all_T_grids_dch(valid_channels, :) = t_d_interp';
    end
end

all_T_grids_chg = all_T_grids_chg(1:valid_channels, :);
all_T_grids_dch = all_T_grids_dch(1:valid_channels, :);
fprintf('  Static data collected from %d channels.\n', valid_channels);

% === Step 2: Average time grids and create global voltage boundaries ===
avg_T_chg = mean(all_T_grids_chg, 1, 'omitnan');
avg_T_dch = mean(all_T_grids_dch, 1, 'omitnan');
avg_T_chg = fillmissing(avg_T_chg, 'linear', 'EndValues', 'nearest');
avg_T_dch = fillmissing(avg_T_dch, 'linear', 'EndValues', 'nearest');

% Ensure monotonic
[avg_T_chg, uid_c] = unique(avg_T_chg, 'stable');
V_chg_std = V_standard_chg(uid_c);
[avg_T_dch, uid_d] = unique(avg_T_dch, 'stable');
V_dch_std = V_standard_dch(uid_d);

% 시간 균등 분할 → 전압 경계
T_start_chg = min(avg_T_chg); T_end_chg = max(avg_T_chg);
target_Ts_chg = linspace(T_start_chg, T_end_chg, num_segments + 1);
Global_V_bounds_chg = interp1(avg_T_chg, V_chg_std, target_Ts_chg, 'linear');

T_start_dch = min(avg_T_dch); T_end_dch = max(avg_T_dch);
target_Ts_dch = linspace(T_start_dch, T_end_dch, num_segments + 1);
Global_V_bounds_dch = interp1(avg_T_dch, V_dch_std, target_Ts_dch, 'linear');

fprintf('  Global voltage boundaries (from Static discharge V-t):\n');
fprintf('  Charge  (%.3f~%.3fV): ', win_chg_min, win_chg_max);
fprintf('%.4f ', Global_V_bounds_chg); fprintf('\n');
fprintf('  Discharge (%.3f~%.3fV): ', win_dch_min, win_dch_max);
fprintf('%.4f ', Global_V_bounds_dch); fprintf('\n');

% === Step 3: Store global ruler for all channels ===
for i = 1:length(channels)
    ch = channels{i};
    
    if ~isfield(RPT_VQ_grid, 'cyc0') || ~isfield(RPT_VQ_grid.cyc0, ch)
        continue;
    end
    
    % All channels use the same global ruler
    % Charge V_bounds → 오름차순 (충전 V_grid가 오름차순이므로)
    MasterRulers.(ch).V_bounds_chg = sort(Global_V_bounds_chg, 'ascend');
    MasterRulers.(ch).V_bounds_dch = Global_V_bounds_dch;
    
    % Store OCV data
    if isfield(RPT_VQ_grid.cyc0.(ch), 'OCV_charge')
        MasterRulers.(ch).Fresh_OCV_Charge = RPT_VQ_grid.cyc0.(ch).OCV_charge;
    else
        MasterRulers.(ch).Fresh_OCV_Charge = [];
    end
    
    fprintf('  %s: Global Ruler assigned.\n', ch);
end
save(fullfile(saveDir, 'MasterRulers.mat'), 'MasterRulers');

%% ========================================================================
% Section 3: Phase 2 - Feature Extraction Loop
% ========================================================================
fprintf('\nPhase 2: Extracting Features...\n');
data_CellID = {};
data_Cycle = [];
data_Crate = {}; 
data_CrateNum = []; 
data_Features = [];
data_Labels = [];
target_crates = {'c01', 'c05', 'c1', 'c2', 'c3'};
target_crates_val = [0.1, 0.5, 1, 2, 3];
cnt = 0;
for i = 1:length(channels)
    ch = channels{i};
    ch_key = ch;
    
    
    if ~isfield(MasterRulers, ch)
        continue;
    end
    
    % Step 3.1: Get Rated Capacity for this channel (Q_rated from Cycle 0)
    % Try to find Cycle 0 capacity from allChannelsCapacity
    if isfield(allChannelsCapacity, ch)
        static_cycles_ch = allChannelsCapacity.(ch).cycles;
        idx_fresh = find(static_cycles_ch == 0, 1);
        if ~isempty(idx_fresh)
             Q_rated_ch = max(allChannelsCapacity.(ch).Q{idx_fresh});
        else
             % Fallback: Try to use MasterRulers Fresh OCV Q max if available
             if isfield(MasterRulers.(ch), 'Fresh_OCV_Charge')
                 Q_rated_ch = max(MasterRulers.(ch).Fresh_OCV_Charge.Q);
             else
                 Q_rated_ch = 55.6; % Hard fallback
             end
        end
    else
        Q_rated_ch = 55.6;
    end
    
    V_ruler_chg = MasterRulers.(ch).V_bounds_chg;
    V_ruler_dch = MasterRulers.(ch).V_bounds_dch;
    cyc_fields = fieldnames(RPT_VQ_grid);
    
    for c = 1:length(cyc_fields)
        cyc_key = cyc_fields{c};
        cyc_num = sscanf(cyc_key, 'cyc%d');
        
        if ~isfield(RPT_VQ_grid.(cyc_key), ch_key)
            continue;
        end
        
        ch_data = RPT_VQ_grid.(cyc_key).(ch_key);
        
        for r = 1:length(target_crates)
            crate_label = target_crates{r};
            crate_val = target_crates_val(r);
            f_chg = [crate_label '_charge'];
            f_dch = [crate_label '_discharge'];
            
            if ~isfield(ch_data, f_chg) || ~isfield(ch_data, f_dch)
                continue;
            end
            
            [fea_chg_dQ, fea_chg_diff] = extract_features_half(ch_data.(f_chg), ...
                V_ruler_chg, 'charge', win_chg_min, win_chg_max);
            [fea_dch_dQ, fea_dch_diff] = extract_features_half(ch_data.(f_dch), ...
                V_ruler_dch, 'discharge', win_dch_min, win_dch_max);
            
            row_features = [fea_chg_dQ, fea_dch_dQ, fea_chg_diff, fea_dch_diff];
            
            static_cycles = allChannelsCapacity.(ch).cycles;
            idx_s = find(static_cycles == cyc_num, 1);
            if ~isempty(idx_s)
                 q_vec = allChannelsCapacity.(ch).Q{idx_s};
                 lbl_soh_cap = max(q_vec); 
            else
                lbl_soh_cap = NaN;
            end
            
            if isfield(ch_data, 'OCV_charge') && isfield(MasterRulers.(ch), 'Fresh_OCV_Charge')
                [lbl_lli, lbl_lam] = analyze_ica_aging(ch_data.OCV_charge, ...
                    MasterRulers.(ch).Fresh_OCV_Charge, Q_rated_ch);
            else
                lbl_lli = NaN;
                lbl_lam = NaN;
            end
            
            cnt = cnt + 1;
            data_CellID{cnt,1} = ch;
            data_Cycle(cnt,1) = cyc_num;
            data_Crate{cnt,1} = crate_label;
            data_CrateNum(cnt,1) = crate_val;
            data_Features(cnt,:) = row_features;
            data_Labels(cnt,:) = [lbl_soh_cap, lbl_lli, lbl_lam];
        end
    end
end
%% ========================================================================
% Section 4: Matrix Integration & Standardization
% ========================================================================
fprintf('\nPhase 4: Integration & Standardization...\n');
FeatureTable = table(data_CellID, data_Cycle, data_Crate, data_CrateNum, ...
    data_Features, data_Labels, 'VariableNames', ...
    {'CellID', 'Cycle', 'CrateLabel', 'CrateNum', 'X_Features', 'Y_Labels'});
% Calculate total number of features dynamically
num_features = 2*num_segments + 4;  % charge_dQ + discharge_dQ + 4 dQ/dV features
mean_X = zeros(length(target_crates), num_features);
std_X  = zeros(length(target_crates), num_features);
FeatureTable.X_Normalized = FeatureTable.X_Features;
for r = 1:length(target_crates)
    c_val = target_crates_val(r);
    idx = FeatureTable.CrateNum == c_val;
    X_sub = FeatureTable.X_Features(idx, :);
    
    if isempty(X_sub)
        continue;
    end
    
    mu = mean(X_sub, 1, 'omitnan');
    sigma = std(X_sub, 0, 1, 'omitnan');
    sigma(sigma==0) = 1; 
    X_norm = (X_sub - mu) ./ sigma;
    FeatureTable.X_Normalized(idx, :) = X_norm;
    mean_X(r,:) = mu;
    std_X(r,:) = sigma;
end
savePath = fullfile(saveDir, 'Feature_Matrix_Final.mat');
save(savePath, 'FeatureTable', 'mean_X', 'std_X');
fprintf('Saved Final Feature Matrix to: %s\n', savePath);

%% ========================================================================
% Section 5: Export to Excel with Feature Names
% ========================================================================
fprintf('\nExporting to Excel with Feature Names...\n');

% Define feature names dynamically based on num_segments
num_features = 2*num_segments + 4;
feature_names = cell(1, num_features);

% Charge dQ features (dynamic segments)
for seg = 1:num_segments
    feature_names{seg} = sprintf('Chg_dQ_Seg%d', seg);
end

% Discharge dQ features (dynamic segments)
for seg = 1:num_segments
    feature_names{num_segments + seg} = sprintf('Dch_dQ_Seg%d', seg);
end

% Charge dQ/dV features (peak height and area)
feature_names{2*num_segments + 1} = 'Chg_dQdV_PkHeight';
feature_names{2*num_segments + 2} = 'Chg_dQdV_PkArea';

% Discharge dQ/dV features (peak height and area)
feature_names{2*num_segments + 3} = 'Dch_dQdV_PkHeight';
feature_names{2*num_segments + 4} = 'Dch_dQdV_PkArea';

% Create expanded table with individual feature columns
ExportTable = table();
ExportTable.CellID = data_CellID;
ExportTable.Cycle = data_Cycle;
ExportTable.CrateLabel = data_Crate;
ExportTable.CrateNum = data_CrateNum;

% Add individual feature columns with names
for i = 1:num_features
    ExportTable.(feature_names{i}) = data_Features(:, i);
end

% Add normalized features
for i = 1:num_features
    col_name = sprintf('%s_Norm', feature_names{i});
    ExportTable.(col_name) = FeatureTable.X_Normalized(:, i);
end

% Add label columns
ExportTable.SOH_Capacity = data_Labels(:, 1);
ExportTable.LLI = data_Labels(:, 2);
ExportTable.LAM = data_Labels(:, 3);

% Save to Excel
excelPath = fullfile(saveDir, 'Feature_Matrix_Final.xlsx');
writetable(ExportTable, excelPath, 'Sheet', 'Features');
fprintf('Saved Feature Table to Excel: %s\n', excelPath);

%% ========================================================================
% Section 6: Visualization - Segment-wise Q-V Curves and dQ/dV Features
% ========================================================================
fprintf('\nGenerating Visualizations...\n');

% Target channels for visualization
vis_channels = {'Ch09', 'Ch10', 'Ch11', 'Ch12', 'Ch13', 'Ch14', 'Ch15', 'Ch16'};
figDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis', 'Feature_Extraction_Advanced');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

% Process each C-rate
for r = 1:length(target_crates)
    crate_label = target_crates{r};
    crate_val = target_crates_val(r);
    
    fprintf('  Plotting C-rate: %s (%.1fC)\n', crate_label, crate_val);
    
    % === CHARGE SEGMENT VISUALIZATION ===
    % Determine subplot layout based on num_segments
    if num_segments <= 5
        nrows = 1; ncols = num_segments;
    elseif num_segments <= 10
        nrows = 2; ncols = 5;
    elseif num_segments <= 15
        nrows = 3; ncols = 5;
    else
        nrows = 4; ncols = 5;
    end
    
    fig_chg = figure('Position', [50, 50, 1600, 1000], 'Visible', 'off', 'Name', sprintf('Charge Segments - %s', crate_label));
    
    % Get all cycle fields
    cyc_fields = fieldnames(RPT_VQ_grid);
    cycle_colors = lines(length(cyc_fields)); % Different color for each cycle
    
    for seg_idx = 1:num_segments
        subplot(nrows, ncols, seg_idx);
        hold on;
        
        % Get first channel's ruler for voltage range in title
        first_ch = vis_channels{1};
        if isfield(MasterRulers, first_ch)
            V_ruler_chg = MasterRulers.(first_ch).V_bounds_chg;
            V_start_title = V_ruler_chg(seg_idx);
            V_end_title = V_ruler_chg(seg_idx + 1);
        end
        
        % Track which cycles have been plotted for legend
        cycles_plotted = false(length(cyc_fields), 1);
        
        % Plot each channel
        for ch_idx = 1:length(vis_channels)
            ch = vis_channels{ch_idx};
            
            if ~isfield(MasterRulers, ch)
                continue;
            end
            
            V_ruler_chg = MasterRulers.(ch).V_bounds_chg;
            f_chg = [crate_label '_charge'];
            
            % Loop through all cycles
            for cyc_idx = 1:length(cyc_fields)
                cyc_key = cyc_fields{cyc_idx};
                
                if isfield(RPT_VQ_grid, cyc_key) && isfield(RPT_VQ_grid.(cyc_key), ch) && ...
                   isfield(RPT_VQ_grid.(cyc_key).(ch), f_chg)
                    
                    data_struct = RPT_VQ_grid.(cyc_key).(ch).(f_chg);
                    V_grid = data_struct.V_grid;
                    Q_val = data_struct.Q;
                    
                    % Extract segment data with forced boundaries
                    V_start = V_ruler_chg(seg_idx);
                    V_end = V_ruler_chg(seg_idx + 1);
                    
                    % Get grid data within segment (approximate)
                    mask = V_grid >= V_start & V_grid <= V_end;
                    V_seg = V_grid(mask);
                    Q_seg = Q_val(mask);
                    
                    if ~isempty(V_seg) && ~isempty(Q_seg)
                        % Interpolate Q at exact boundary voltages
                        [V_unique, uid] = unique(V_grid);
                        Q_unique = Q_val(uid);                        
                        Q_start = interp1(V_unique, Q_unique, V_start, 'linear');
                        Q_end   = interp1(V_unique, Q_unique, V_end, 'linear');

                        % Remove points too close to boundaries (within 0.0005V)
                        eps_v = 0.0005;
                        mask_interior = (V_seg > V_start + eps_v) & (V_seg < V_end - eps_v);
                        V_seg_clean = V_seg(mask_interior);
                        Q_seg_clean = Q_seg(mask_interior);
                        
                        % Force insert exact boundary points
                        V_seg_forced = [V_start; V_seg_clean; V_end];
                        Q_seg_forced = [Q_start; Q_seg_clean; Q_end];
                        
                        % Sort only (no unique to preserve exact boundaries)
                        [V_seg_forced, sort_idx] = sort(V_seg_forced, 'ascend');
                        Q_seg_forced = Q_seg_forced(sort_idx);
                        
                        % Normalize Q to start from 0
                        Q_seg_norm = Q_seg_forced - Q_start;
                        
                        % Only add to legend if this cycle hasn't been plotted yet
                        if ~cycles_plotted(cyc_idx)
                            plot(Q_seg_norm, V_seg_forced, 'LineWidth', 1.2, 'Color', cycle_colors(cyc_idx, :),'DisplayName', cyc_key);
                            cycles_plotted(cyc_idx) = true;
                        else
                            plot(Q_seg_norm, V_seg_forced, 'LineWidth', 1.2, 'Color', cycle_colors(cyc_idx, :),'HandleVisibility', 'off');
                        end
                    end
                end
            end
        end
        
        xlabel('Capacity (Ah)', 'FontSize', 10);
        ylabel('Voltage (V)', 'FontSize', 10);
        title(sprintf('Segment %d (%.2f-%.2fV)', seg_idx, V_start_title, V_end_title), ...
            'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        if seg_idx == 1
            legend('Location', 'best', 'FontSize', 6);
        end
        hold off;
    end
    
    sgtitle(sprintf('Charge Q-V Segments - All Cycles - %s (%.1fC)', crate_label, crate_val), ...
        'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig_chg, fullfile(figDir, sprintf('Charge_Segments_%s.fig', crate_label)));
    close(fig_chg);
    
    % === DISCHARGE SEGMENT VISUALIZATION ===
    fig_dch = figure('Position', [50, 50, 1600, 1000], 'Visible', 'off', 'Name', sprintf('Discharge Segments - %s', crate_label));
    
    for seg_idx = 1:num_segments
        subplot(nrows, ncols, seg_idx);
        hold on;
        
        % Get first channel's ruler for voltage range in title
        first_ch = vis_channels{1};
        if isfield(MasterRulers, first_ch)
            V_ruler_dch = MasterRulers.(first_ch).V_bounds_dch;
            V_start_title = V_ruler_dch(seg_idx);
            V_end_title = V_ruler_dch(seg_idx + 1);
        end
        
        % Track which cycles have been plotted for legend
        cycles_plotted = false(length(cyc_fields), 1);
        
        % Plot each channel
        for ch_idx = 1:length(vis_channels)
            ch = vis_channels{ch_idx};
            
            if ~isfield(MasterRulers, ch)
                continue;
            end
            
            V_ruler_dch = MasterRulers.(ch).V_bounds_dch;
            f_dch = [crate_label '_discharge'];
            
            % Loop through all cycles
            for cyc_idx = 1:length(cyc_fields)
                cyc_key = cyc_fields{cyc_idx};
                
                if isfield(RPT_VQ_grid, cyc_key) && isfield(RPT_VQ_grid.(cyc_key), ch) && ...
                   isfield(RPT_VQ_grid.(cyc_key).(ch), f_dch)
                    
                    data_struct = RPT_VQ_grid.(cyc_key).(ch).(f_dch);
                    V_grid = data_struct.V_grid;
                    Q_val = data_struct.Q;
                    
                    % Extract segment data with forced boundaries
                    V_start = V_ruler_dch(seg_idx);
                    V_end = V_ruler_dch(seg_idx + 1);
                    
                    % Get grid data within segment (approximate)
                    mask = V_grid <= V_start & V_grid >= V_end;
                    V_seg = V_grid(mask);
                    Q_seg = Q_val(mask);
                    
                    if ~isempty(V_seg) && ~isempty(Q_seg)
                        % Interpolate Q at exact boundary voltages
                        [V_unique, uid] = unique(V_grid);
                        Q_unique = Q_val(uid);
                        Q_start = interp1(V_unique, Q_unique, V_start, 'linear');
                        Q_end = interp1(V_unique, Q_unique, V_end, 'linear');
                        
                        % Remove points too close to boundaries (within 0.0005V)
                        eps_v = 0.0005;
                        mask_interior = (V_seg < V_start - eps_v) & (V_seg > V_end + eps_v);
                        V_seg_clean = V_seg(mask_interior);
                        Q_seg_clean = Q_seg(mask_interior);
                        
                        % Force insert exact boundary points
                        V_seg_forced = [V_start; V_seg_clean; V_end];
                        Q_seg_forced = [Q_start; Q_seg_clean; Q_end];
                        
                        % Sort only (no unique to preserve exact boundaries, descending for discharge)
                        [V_seg_forced, sort_idx] = sort(V_seg_forced, 'descend');
                        Q_seg_forced = Q_seg_forced(sort_idx);
                        
                        % Normalize Q to start from 0
                        Q_seg_norm = Q_seg_forced - Q_start;
                        
                        % Only add to legend if this cycle hasn't been plotted yet
                        if ~cycles_plotted(cyc_idx)
                            plot(Q_seg_norm, V_seg_forced, 'LineWidth', 1.2, 'Color', cycle_colors(cyc_idx, :), ...
                                'DisplayName', cyc_key);
                            cycles_plotted(cyc_idx) = true;
                        else
                            plot(Q_seg_norm, V_seg_forced, 'LineWidth', 1.2, 'Color', cycle_colors(cyc_idx, :), ...
                                'HandleVisibility', 'off');
                        end
                    end
                end
            end
        end
        
        xlabel('Capacity (Ah)', 'FontSize', 10);
        ylabel('Voltage (V)', 'FontSize', 10);
        title(sprintf('Segment %d (%.2f-%.2fV)', seg_idx, V_start_title, V_end_title), ...
            'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        if seg_idx == 1
            legend('Location', 'best', 'FontSize', 6);
        end
        hold off;
    end
    
    sgtitle(sprintf('Discharge Q-V Segments - All Cycles - %s (%.1fC)', crate_label, crate_val), ...
        'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig_dch, fullfile(figDir, sprintf('Discharge_Segments_%s.fig', crate_label)));
    close(fig_dch);
    
    % === CHARGE dQ/dV SEGMENT VISUALIZATION ===
    fig_dqdv_chg = figure('Position', [50, 50, 1600, 1000], 'Visible', 'off', 'Name', sprintf('dQdV Charge Segments - %s', crate_label));
    
    for seg_idx = 1:num_segments
        subplot(nrows, ncols, seg_idx);
        hold on;
        
        % Get first channel's ruler for voltage range in title
        first_ch = vis_channels{1};
        if isfield(MasterRulers, first_ch)
            V_ruler_chg = MasterRulers.(first_ch).V_bounds_chg;
            V_start_title = V_ruler_chg(seg_idx);
            V_end_title = V_ruler_chg(seg_idx + 1);
        end
        
        % Track which cycles have been plotted for legend
        cycles_plotted = false(length(cyc_fields), 1);
        
        % Plot each channel
        for ch_idx = 1:length(vis_channels)
            ch = vis_channels{ch_idx};
            
            if ~isfield(MasterRulers, ch)
                continue;
            end
            
            V_ruler_chg = MasterRulers.(ch).V_bounds_chg;
            f_chg = [crate_label '_charge'];
            
            % Loop through all cycles
            for cyc_idx = 1:length(cyc_fields)
                cyc_key = cyc_fields{cyc_idx};
                
                if isfield(RPT_VQ_grid, cyc_key) && isfield(RPT_VQ_grid.(cyc_key), ch) && ...
                   isfield(RPT_VQ_grid.(cyc_key).(ch), f_chg)
                    
                    data_struct = RPT_VQ_grid.(cyc_key).(ch).(f_chg);
                    V_grid = data_struct.V_grid;
                    Q_val = data_struct.Q;
                    
                    [V_u, uid] = unique(V_grid);
                    Q_u = Q_val(uid);
                    
                    if numel(V_u) > 10
                        dV = gradient(V_u);
                        dQ = gradient(Q_u);
                        dV(dV==0) = NaN;
                        dQdV_raw = dQ ./ dV;
                        
                        % Apply Savitzky-Golay filter
                        if length(dQdV_raw) > 51
                            try
                                dQdV_filt = sgolayfilt(dQdV_raw, 3, 51);
                            catch
                                dQdV_filt = dQdV_raw;
                            end
                        else
                            dQdV_filt = dQdV_raw;
                        end
                        
                        % Extract segment data
                        V_start = V_ruler_chg(seg_idx);
                        V_end = V_ruler_chg(seg_idx + 1);
                        mask = V_u >= V_start & V_u <= V_end;
                        
                        if sum(mask) > 0
                            % Only add to legend if this cycle hasn't been plotted yet
                            if ~cycles_plotted(cyc_idx)
                                plot(V_u(mask), dQdV_filt(mask), 'LineWidth', 1.2, 'Color', cycle_colors(cyc_idx, :), ...
                                    'DisplayName', cyc_key);
                                cycles_plotted(cyc_idx) = true;
                            else
                                plot(V_u(mask), dQdV_filt(mask), 'LineWidth', 1.2, 'Color', cycle_colors(cyc_idx, :), ...
                                    'HandleVisibility', 'off');
                            end
                        end
                    end
                end
            end
        end
        
        xlabel('Voltage (V)', 'FontSize', 10);
        ylabel('dQ/dV (Ah/V)', 'FontSize', 10);
        title(sprintf('Segment %d (%.2f-%.2fV)', seg_idx, V_start_title, V_end_title), ...
            'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        if seg_idx == 1
            legend('Location', 'best', 'FontSize', 6);
        end
        hold off;
    end
    
    sgtitle(sprintf('Charge dQ/dV Segments - All Cycles - %s (%.1fC)', crate_label, crate_val), ...
        'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig_dqdv_chg, fullfile(figDir, sprintf('dQdV_Charge_Segments_%s.fig', crate_label)));
    close(fig_dqdv_chg);
    
    % === DISCHARGE dQ/dV SEGMENT VISUALIZATION ===
    fig_dqdv_dch = figure('Position', [50, 50, 1600, 1000], 'Visible', 'off', 'Name', sprintf('dQdV Discharge Segments - %s', crate_label));
    
    for seg_idx = 1:num_segments
        subplot(nrows, ncols, seg_idx);
        hold on;
        
        % Get first channel's ruler for voltage range in title
        first_ch = vis_channels{1};
        if isfield(MasterRulers, first_ch)
            V_ruler_dch = MasterRulers.(first_ch).V_bounds_dch;
            V_start_title = V_ruler_dch(seg_idx);
            V_end_title = V_ruler_dch(seg_idx + 1);
        end
        
        % Track which cycles have been plotted for legend
        cycles_plotted = false(length(cyc_fields), 1);
        
        % Plot each channel
        for ch_idx = 1:length(vis_channels)
            ch = vis_channels{ch_idx};
            
            if ~isfield(MasterRulers, ch)
                continue;
            end
            
            V_ruler_dch = MasterRulers.(ch).V_bounds_dch;
            f_dch = [crate_label '_discharge'];
            
            % Loop through all cycles
            for cyc_idx = 1:length(cyc_fields)
                cyc_key = cyc_fields{cyc_idx};
                
                if isfield(RPT_VQ_grid, cyc_key) && isfield(RPT_VQ_grid.(cyc_key), ch) && ...
                   isfield(RPT_VQ_grid.(cyc_key).(ch), f_dch)
                    
                    data_struct = RPT_VQ_grid.(cyc_key).(ch).(f_dch);
                    V_grid = data_struct.V_grid;
                    Q_val = data_struct.Q;
                    
                    [V_u, uid] = unique(V_grid);
                    Q_u = Q_val(uid);
                    
                    if numel(V_u) > 10
                        dV = gradient(V_u);
                        dQ = gradient(Q_u);
                        dV(dV==0) = NaN;
                        dQdV_raw = dQ ./ dV;
                        
                        % Apply Savitzky-Golay filter
                        if length(dQdV_raw) > 51
                            try
                                dQdV_filt = sgolayfilt(dQdV_raw, 3, 51);
                            catch
                                dQdV_filt = dQdV_raw;
                            end
                        else
                            dQdV_filt = dQdV_raw;
                        end
                        
                        % Extract segment data
                        V_start = V_ruler_dch(seg_idx);
                        V_end = V_ruler_dch(seg_idx + 1);
                        mask = V_u <= V_start & V_u >= V_end;
                        
                        if sum(mask) > 0
                            % Only add to legend if this cycle hasn't been plotted yet
                            if ~cycles_plotted(cyc_idx)
                                plot(V_u(mask), dQdV_filt(mask), 'LineWidth', 1.2, 'Color', cycle_colors(cyc_idx, :), ...
                                    'DisplayName', cyc_key);
                                cycles_plotted(cyc_idx) = true;
                            else
                                plot(V_u(mask), dQdV_filt(mask), 'LineWidth', 1.2, 'Color', cycle_colors(cyc_idx, :), ...
                                    'HandleVisibility', 'off');
                            end
                        end
                    end
                end
            end
        end
        
        xlabel('Voltage (V)', 'FontSize', 10);
        ylabel('dQ/dV (Ah/V)', 'FontSize', 10);
        title(sprintf('Segment %d (%.2f-%.2fV)', seg_idx, V_start_title, V_end_title), ...
            'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        if seg_idx == 1
            legend('Location', 'best', 'FontSize', 6);
        end
        hold off;
    end
    
    sgtitle(sprintf('Discharge dQ/dV Segments - All Cycles - %s (%.1fC)', crate_label, crate_val), ...
        'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig_dqdv_dch, fullfile(figDir, sprintf('dQdV_Discharge_Segments_%s.fig', crate_label)));
    close(fig_dqdv_dch);
end

%% ========================================================================
% Section 7: C-rate Comparison Visualization (dQ/dV vs V)
% ========================================================================
fprintf('\nGenerating C-rate Comparison Visualizations...\n');

% Get all cycle fields
cyc_fields = fieldnames(RPT_VQ_grid);
cycle_colors = lines(length(cyc_fields));

% === CHARGE dQ/dV Comparison Across C-rates ===
fig_crate_chg = figure('Position', [50, 50, 1800, 1000], 'Visible', 'off', 'Name', 'dQdV Charge - C-rate Comparison');

% Calculate global Y-limit for Charge
% Column 11: Chg PkH
if ~isempty(data_Features)
    y_max_chg = max(data_Features(:, 11), [], 'omitnan') * 1.1; 
else
    y_max_chg = 10; % Fallback
end
if isempty(y_max_chg) || isnan(y_max_chg) || y_max_chg <= 0
    y_max_chg = 10;
end



% Calculate global Y-limit for Discharge
% Column 13: Dch PkH
if ~isempty(data_Features)
    Y_LIMIT_DCH = max(data_Features(:, 13), [], 'omitnan') * 1.1; 
else
    Y_LIMIT_DCH = 10; % Fallback
end

% Strict validation
if isempty(Y_LIMIT_DCH) || isnan(Y_LIMIT_DCH) || ~isreal(Y_LIMIT_DCH) || Y_LIMIT_DCH <= 0
    Y_LIMIT_DCH = 10;
else
    Y_LIMIT_DCH = double(Y_LIMIT_DCH(1)); % Force scalar
end
fprintf('DEBUG: Y_LIMIT_DCH = %f\n', Y_LIMIT_DCH);

for r = 1:length(target_crates)
    crate_label = target_crates{r};
    crate_val = target_crates_val(r);
    
    subplot(2, 3, r);
    hold on;
    
    % Track which cycles have been plotted for legend
    cycles_plotted = false(length(cyc_fields), 1);
    
    % --- [Modified] Visualize Feature Extraction Window (Shaded Area) ---
    fill([win_chg_min, win_chg_max, win_chg_max, win_chg_min], ...
         [0, 0, y_max_chg, y_max_chg], ...
         [0.8, 0.8, 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Plot all channels and all cycles for this C-rate
    for i = 1:length(vis_channels)
        ch = vis_channels{i};
        
        % [Modified] Restrict to Channel 9 only
        if ~strcmpi(ch, 'ch09')
            continue;
        end
        
        if ~isfield(MasterRulers, ch)
            continue;
        end
        
        f_chg = [crate_label '_charge'];
        
        % Loop through all cycles
        for cyc_idx = 1:length(cyc_fields)
            cyc_key = cyc_fields{cyc_idx};
            
            if isfield(RPT_VQ_grid, cyc_key) && isfield(RPT_VQ_grid.(cyc_key), ch) && ...
               isfield(RPT_VQ_grid.(cyc_key).(ch), f_chg)
                
                data_struct = RPT_VQ_grid.(cyc_key).(ch).(f_chg);
                V_grid = data_struct.V_grid;
                Q_val = data_struct.Q;
                
                % Remove duplicates
                [V_u, uid] = unique(V_grid);
                Q_u = Q_val(uid);
                
                if numel(V_u) > 10
                    % Calculate dQ/dV
                    dV = gradient(V_u);
                    dQ = gradient(Q_u);
                    dV(dV==0) = NaN;
                    dQdV_raw = dQ ./ dV;
                    
                    % Apply Moving Average filter (consistent with feature extraction)
                    window_size = 21;
                    dQdV_filt = movmean(dQdV_raw, window_size);
                    
                    % Extract window region
                    % [Modified] Show Full Range (User Request)
                    % mask = V_u >= win_chg_min & V_u <= win_chg_max;
                    mask = true(size(V_u));
                    
                    if sum(mask) > 0
                        % Plot Line
                        if ~cycles_plotted(cyc_idx)
                            plot(V_u(mask), dQdV_filt(mask), 'LineWidth', 1.5, 'Color', cycle_colors(cyc_idx, :), ...
                                'DisplayName', cyc_key);
                            cycles_plotted(cyc_idx) = true;
                        else
                            plot(V_u(mask), dQdV_filt(mask), 'LineWidth', 1.5, 'Color', cycle_colors(cyc_idx, :), ...
                                'HandleVisibility', 'off');
                        end
                        
                        % [Modified] Mark Peak Position
                        [max_val, max_idx] = max(dQdV_filt(mask));
                        V_mask = V_u(mask); 
                        V_peak = V_mask(max_idx);
                        plot(V_peak, max_val, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
                             'MarkerFaceColor', cycle_colors(cyc_idx, :), 'HandleVisibility', 'off');
                    end
                end
            end
        end
    end

    
    xlabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('dQ/dV (Ah/V)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('%s (%.1fC)', crate_label, crate_val), 'FontSize', 13, 'FontWeight', 'bold');
    grid on;
    % [Modified] Show Full Voltage Range
    xlim([3.0, 4.2]);
    ylim([0, y_max_chg]);
    
    if r == 1
        legend('Location', 'best', 'FontSize', 8);
    end
    hold off;
end

sgtitle('Charge dQ/dV - C-rate Comparison (All Cycles)', 'FontSize', 16, 'FontWeight', 'bold');
saveas(fig_crate_chg, fullfile(figDir, 'dQdV_Charge_Crate_Comparison.fig'));
close(fig_crate_chg);

% === DISCHARGE dQ/dV Comparison Across C-rates ===
fig_crate_dch = figure('Position', [50, 50, 1800, 1000], 'Visible', 'off', 'Name', 'dQdV Discharge - C-rate Comparison');

% Calculate global Y-limit for Discharge
% Column 13: Dch PkH
if ~isempty(data_Features)
    y_max_dch = max(abs(data_Features(:, 13))) * 1.2; % [Modified] Use Abs Max for Discharge
else
    y_max_dch = 10; % Fallback
end

for r = 1:length(target_crates)
    crate_label = target_crates{r};
    crate_val = target_crates_val(r);
    
    subplot(2, 3, r);
    hold on;
    
    % Track which cycles have been plotted for legend
    cycles_plotted = false(length(cyc_fields), 1);
    
    % --- [Modified] Visualize Feature Extraction Window (Shaded Area) ---
    % Moved to after plotting to determine correct Y-limit dynamically

    % Plot all channels and all cycles for this C-rate
    for i = 1:length(vis_channels)
        ch = vis_channels{i};
        
        % [Modified] Restrict to Channel 9 only
        if ~strcmpi(ch, 'ch09')
            continue;
        end
        
        if ~isfield(MasterRulers, ch)
            continue;
        end
        
        f_dch = [crate_label '_discharge'];
        
        % Loop through all cycles
        for cyc_idx = 1:length(cyc_fields)
            cyc_key = cyc_fields{cyc_idx};
            
            if isfield(RPT_VQ_grid, cyc_key) && isfield(RPT_VQ_grid.(cyc_key), ch) && ...
               isfield(RPT_VQ_grid.(cyc_key).(ch), f_dch)
                

                data_struct = RPT_VQ_grid.(cyc_key).(ch).(f_dch);
                V_grid = data_struct.V_grid;
                Q_val = data_struct.Q;
                
                % Remove duplicates
                [V_u, uid] = unique(V_grid);
                Q_u = Q_val(uid);
                
                if numel(V_u) > 10
                    % Calculate dQ/dV
                    dV = gradient(V_u);
                    dQ = gradient(Q_u);
                    dV(dV==0) = NaN;
                    dQdV_raw = dQ ./ dV;
                    
                    % Apply Moving Average filter (consistent, Window=21)
                    window_size = 21;
                    dQdV_filt = abs(movmean(dQdV_raw, window_size)); % [Modified] Plot Absolute Value for Discharge
                    
                    % Extract window region
                    % [Modified] Show Full Range (User Request)
                    % mask = V_u <= win_dch_max & V_u >= win_dch_min;
                    mask = true(size(V_u));
                    
                    if sum(mask) > 0
                        % Plot Line
                        if ~cycles_plotted(cyc_idx)
                            plot(V_u(mask), dQdV_filt(mask), 'LineWidth', 1.5, 'Color', cycle_colors(cyc_idx, :), ...
                                'DisplayName', cyc_key);
                            cycles_plotted(cyc_idx) = true;
                        else
                            plot(V_u(mask), dQdV_filt(mask), 'LineWidth', 1.5, 'Color', cycle_colors(cyc_idx, :), ...
                                'HandleVisibility', 'off');
                        end
                        
                        % [Modified] Mark Peak Position
                        [max_val, max_idx] = max(dQdV_filt(mask));
                        V_mask = V_u(mask); 
                        V_peak = V_mask(max_idx);
                        plot(V_peak, max_val, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
                             'MarkerFaceColor', cycle_colors(cyc_idx, :), 'HandleVisibility', 'off');
                    end
                end
            end
        end
    end
    
    xlabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('dQ/dV (Ah/V)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('%s (%.1fC)', crate_label, crate_val), 'FontSize', 13, 'FontWeight', 'bold');
    grid on;
    % [Modified] Show Full Voltage Range
    xlim([3.0, 4.2]);
    ylim auto; % Let axes scale automatically to show the curve
    
    % --- [Modified] Draw Shaded Window Dynamically ---
    % Get current Y-axis limits (determined by the data)
    yl = ylim; 
    y_top = yl(2);
    
    % Draw the shading extending to the top of the graph
    h_fill = fill([win_dch_min, win_dch_max, win_dch_max, win_dch_min], ...
                  [0, 0, y_top, y_top], ...
                  [0.8, 0.8, 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                  
    % Move shading to the background (behind the plot lines)
    uistack(h_fill, 'bottom');
    
    if r == 1
        legend('Location', 'best', 'FontSize', 8);
    end
    hold off;
end

sgtitle('Discharge dQ/dV - C-rate Comparison (All Cycles)', 'FontSize', 16, 'FontWeight', 'bold');
saveas(fig_crate_dch, fullfile(figDir, 'dQdV_Discharge_Crate_Comparison.fig'));
close(fig_crate_dch);

fprintf('C-rate comparison visualizations saved.\n');
fprintf('All visualizations saved to: %s\n', figDir);
%% ========================================================================
% Helper Functions
% ========================================================================
function [features_dQ, features_diff] = extract_features_half(data_struct, V_ruler, mode, minV, maxV)
    V_grid = data_struct.V_grid;
    Q_val  = data_struct.Q;
    [V_u, uid] = unique(V_grid); 
    Q_u = Q_val(uid);
    
    if numel(V_u) < 2
         features_dQ = nan(1, length(V_ruler)-1);
         features_diff = [NaN, NaN];
         return;
    end
    Q_ruler = interp1(V_u, Q_u, V_ruler, 'linear');
    features_dQ = abs(diff(Q_ruler));
    
    dV = gradient(V_u);
    dQ = gradient(Q_u);
    dV(dV==0) = NaN;
    dQdV_raw = dQ ./ dV;
    
    % --- Smoothing Method: Moving Average (Window=21, ~0.021V) ---
    % As per documentation Section 2.2 (Features)
    
    window_size = 21;
    dQdV_filt = movmean(dQdV_raw, window_size);
    
    if strcmp(mode, 'charge')
        mask = V_u >= minV & V_u <= maxV;
    else
        mask = V_u <= maxV & V_u >= minV;
    end
    
    dQdV_win = dQdV_filt(mask);
    V_win = V_u(mask);
    
    if isempty(dQdV_win)
        features_diff = [NaN, NaN];
    else
        pk_height = max(dQdV_win);
        pk_area = trapz(V_win, dQdV_win);
        features_diff = [pk_height, abs(pk_area)];
    end
end
function [lli, lam] = analyze_ica_aging(curr_ocv_struct, fresh_ocv_struct, Q_rated)
    if nargin < 3, Q_rated = 55.6; end
    
    if isempty(curr_ocv_struct) || isempty(fresh_ocv_struct)
        lli = NaN;
        lam = NaN;
        return;
    end
    roi_min = 3.40;
    roi_max = 4.00; 
    [peak_V_fresh, peak_H_fresh] = get_main_peak(fresh_ocv_struct, roi_min, roi_max);
    [peak_V_aged, peak_H_aged] = get_main_peak(curr_ocv_struct, roi_min, roi_max);
    
    if isnan(peak_V_fresh) || isnan(peak_V_aged)
        lli = NaN;
        lam = NaN;
    else
        % LLI (%) = |dV_peak| * H_peak_fresh / Q_rated * 100
        lli = (abs(peak_V_aged - peak_V_fresh) * peak_H_fresh / Q_rated) * 100;
        
        % LAM (%) = (1 - Area_peak_aged / Area_peak_fresh) * 100
        % Calculate peak area (approximate by H * fixed width or integration)
        % Since we don't have full integration here, use peak height ratio as proxy for area ratio in consistent peak shape
        % OR better: Integrate dQ/dV around the peak
        
        % Method 1: Use Peak Height Ratio as Area Ratio Proxy (assuming peak width constant)
        % lam = (1 - peak_H_aged / peak_H_fresh) * 100;
        
        % Method 2: Integrate area (more accurate based on slide formula)
        area_fresh = trapz(fresh_ocv_struct.V_grid, abs(movmean(gradient(fresh_ocv_struct.Q)./gradient(fresh_ocv_struct.V_grid), 21)));
        area_aged  = trapz(curr_ocv_struct.V_grid,  abs(movmean(gradient(curr_ocv_struct.Q) ./gradient(curr_ocv_struct.V_grid),  21)));
        
        % Calculate integration around peak window only (3.4~4.0V) to be precise
        mask_f = fresh_ocv_struct.V_grid >= 3.4 & fresh_ocv_struct.V_grid <= 4.0;
        mask_a = curr_ocv_struct.V_grid >= 3.4 & curr_ocv_struct.V_grid <= 4.0;
        
        dQdV_f = gradient(fresh_ocv_struct.Q)./gradient(fresh_ocv_struct.V_grid); dQdV_f(isinf(dQdV_f))=NaN; dQdV_f=fillmissing(dQdV_f,'linear');
        dQdV_a = gradient(curr_ocv_struct.Q)./gradient(curr_ocv_struct.V_grid); dQdV_a(isinf(dQdV_a))=NaN; dQdV_a=fillmissing(dQdV_a,'linear');
        
        area_peak_fresh = trapz(fresh_ocv_struct.V_grid(mask_f), abs(movmean(dQdV_f(mask_f), 21)));
        area_peak_aged  = trapz(curr_ocv_struct.V_grid(mask_a),  abs(movmean(dQdV_a(mask_a), 21)));

        lam = (1 - area_peak_aged / area_peak_fresh) * 100;
    end
end
function [pk_V, pk_H] = get_main_peak(data_struct, minV, maxV)
    % Extracts main peak from OCV dQ/dV curve for LLI/LAM labeling
    % Data Source: OCV Charge Profile
    % Preprocessing: 0.001V Resampling -> dQ/dV -> Moving Average (Window=21)
    
    V = data_struct.V_grid;
    Q = data_struct.Q;
    
    if numel(V) < 10
        pk_V = NaN;
        pk_H = NaN;
        return;
    end
    
    [V_u, uid] = unique(V);
    Q_u = Q(uid);
    dV = gradient(V_u);
    dQ = gradient(Q_u);
    dV(dV==0) = NaN;
    dQdV = dQ ./ dV;
    dQdV(isinf(dQdV)) = NaN;
    dQdV = fillmissing(dQdV, 'linear');
    
    % --- Smoothing Method: Moving Average (Window=21, ~0.021V) ---
    % As per documentation Section 3 (Labels)
    window_size = 21;
    dQdV_filt = movmean(dQdV, window_size);
    
    mask = V_u >= minV & V_u <= maxV;
    V_roi = V_u(mask);
    dH_roi = dQdV_filt(mask);
    
    if isempty(V_roi)
        pk_V = NaN;
        pk_H = NaN;
        return;
    end
    
    [pks, locs] = findpeaks(dH_roi);
    
    if isempty(pks)
        [pk_H, idx] = max(dH_roi); 
        pk_V = V_roi(idx);
    else
        [pk_H, idx] = max(pks);
        pk_V = V_roi(locs(idx));
    end
end
