%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Label Extractor (ver02 Final)
% - 5 Diagnostic Labels for Battery State Estimation
% - Labels:
%   [1] SOH         : (Q_actual / Q_BOL) × 100
%   [2] LLI         : Lithium Loss Index (ICA peak shift, 0.1C OCV)
%   [3] LAM         : Active Material Loss (ICA area shrinkage, 0.1C OCV)
%   [4] SOP_dch_10s : mean P10s across discharge segments (HPPC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ========================================================================
% Section 1: Configuration & Data Loading
% ========================================================================
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_static_mat = fullfile(baseDir, 'Capacity_Trend_Figures', 'Capacity_Data_Static.mat');
path_ocv_avg_mat = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'OCV_integrated.mat');
path_rpt_vq_mat = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');
path_master_ruler = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', ...
    'ExperimentalData', 'FeatureEngineering', 'MasterRulers_v3.mat');

% Output Configuration
currentScriptPath = mfilename('fullpath');
[ver02Dir, ~, ~] = fileparts(currentScriptPath);
saveDir = fullfile(ver02Dir, 'RPT_FeatureLabelExtractor'); % Saving to same dir as features

if ~exist(saveDir, 'dir'), mkdir(saveDir); end

fprintf('Loading Label Data...\n');
load(path_static_mat, 'allChannelsCapacity');
load(path_ocv_avg_mat, 'OCV_data');
load(path_rpt_vq_mat, 'RPT_VQ_grid');
load(path_master_ruler, 'MasterRulers');

channels = fieldnames(allChannelsCapacity);

% Target evaluation points (matching C-rates used in Feature extraction)
target_crates = {'c01', 'c05', 'c1', 'c2', 'c3'};
target_crates_val = [0.1, 0.5, 1, 2, 3];

% Constants for SoP
V_min_sop = 3.0; 
V_max_sop = 4.2; 
soc_grid = linspace(0, 1, length(OCV_data.OCV_integrated_0.V_avg_SOC_charge));

% Initialize SOP Data Struct & Flattened column names
SOP_Data = struct();
num_segments = 5;
sop_col_names = cell(1, 40);
sop_prefixes = {'P2s_chg', 'P10s_chg', 'P30s_chg', 'P60s_chg', 'P2s_dchg', 'P10s_dchg', 'P30s_dchg', 'P60s_dchg'};
col_idx = 1;
for sp = 1:length(sop_prefixes)
    for seg = 1:num_segments
        sop_col_names{col_idx} = sprintf('%s_Seg%d', sop_prefixes{sp}, seg);
        col_idx = col_idx + 1;
    end
end

%% ========================================================================
% Section 2: Label Extraction Loop
% ========================================================================
fprintf('\nExtracting SOH, LLI, LAM, and SOP Labels...\n');
y_CellID = {};
y_Cycle = [];
y_CrateLabel = {};
y_CrateNum = [];
y_Labels = []; % [SOH, LLI, LAM, SOP_flat(40)]

cnt = 0;
for i = 1:length(channels)
    ch = channels{i};
    if ~isfield(allChannelsCapacity, ch) || ~isfield(MasterRulers, ch), continue; end
    cap_data = allChannelsCapacity.(ch);
    
    VR_chg = MasterRulers.(ch).V_bounds_chg;
    VR_dch = MasterRulers.(ch).V_bounds_dch;
    
    if ~isfield(cap_data, 'cycles') || ~isfield(cap_data, 'Q')
        continue;
    end
    
    idx_cyc0 = find(cap_data.cycles == 0, 1);
    if isempty(idx_cyc0)
        idx_cyc0 = 1; % Fallback first recorded
    end
    
    Q_0 = NaN;
    try
        item = cap_data.Q{1, idx_cyc0}; % Static Capacity (low-rate full discharge)
        if ~isempty(item), Q_0 = max(item); end
    catch
        % ignore
    end
    if isnan(Q_0), warning('Q_0 not found for %s, skipping channel.', ch); continue; end
    
    for c = 1:length(cap_data.cycles)
        cyc_num = cap_data.cycles(c);
        cyc_str = sprintf('cyc%d', cyc_num);
        
        % --- 1. Extract SOP for this cycle ---
        current_sop_flat = nan(1, 40);
        SOP_Data.(ch).(cyc_str).P2s_chg = nan(1, num_segments);
        SOP_Data.(ch).(cyc_str).P10s_chg = nan(1, num_segments);
        SOP_Data.(ch).(cyc_str).P30s_chg = nan(1, num_segments);
        SOP_Data.(ch).(cyc_str).P60s_chg = nan(1, num_segments);
        SOP_Data.(ch).(cyc_str).P2s_dchg = nan(1, num_segments);
        SOP_Data.(ch).(cyc_str).P10s_dchg = nan(1, num_segments);
        SOP_Data.(ch).(cyc_str).P30s_dchg = nan(1, num_segments);
        SOP_Data.(ch).(cyc_str).P60s_dchg = nan(1, num_segments);

        parsed_file = fullfile(baseDir, 'RPT', 'Postprocessing', 'Parsed', sprintf('RPT%d_%s_parsed.mat', cyc_num, ch));
        if exist(parsed_file, 'file')
            parsed_data = load(parsed_file, 'pdata');
            if isfield(parsed_data, 'pdata')
                pdata = parsed_data.pdata;
                I_1C = 64.0; tolC = 0.05;
                
                % Determine Cycle-Specific OCV
                ocv_field_name = sprintf('OCV_integrated_%d', cyc_num);
                if ~isfield(OCV_data, ocv_field_name)
                    ocv_field_name = 'OCV_integrated_0'; 
                end
                curr_ocv_struct = OCV_data.(ocv_field_name);
                curr_V_ocv_chg_avg = curr_ocv_struct.V_avg_SOC_charge;
                curr_V_ocv_dch_avg = curr_ocv_struct.V_avg_SOC_discharge;
                
                % Map 5 voltage feature segments to SOC centers using current Cycle's OCV
                [V_u_chg, uid_chg] = unique(curr_V_ocv_chg_avg);
                soc_u_chg = soc_grid(uid_chg);
                soc_bnd_chg = interp1(V_u_chg, soc_u_chg, VR_chg, 'linear');
                soc_mid_chg = (soc_bnd_chg(1:end-1) + soc_bnd_chg(2:end)) / 2;
                SOP_Data.(ch).(cyc_str).soc_mid_chg = soc_mid_chg;
                
                [V_u_dch, uid_dch] = unique(curr_V_ocv_dch_avg);
                soc_u_dch = soc_grid(uid_dch);
                soc_bnd_dch = interp1(V_u_dch, soc_u_dch, VR_dch, 'linear');
                soc_mid_dch = (soc_bnd_dch(1:end-1) + soc_bnd_dch(2:end)) / 2;
                SOP_Data.(ch).(cyc_str).soc_mid_dch = soc_mid_dch;
                
                
                charge_n1C_idx = [];
                disch_n1C_idx = [];
                for k = 1:length(pdata)
                    if isempty(pdata(k).I), continue; end
                    avg_I = mean(pdata(k).I);
                    dur = pdata(k).t(end) - pdata(k).t(1);
                    if avg_I >= I_1C*(1-tolC) && avg_I <= I_1C*(1+tolC) && dur >= 50 && dur <= 70
                        charge_n1C_idx(end+1) = k;
                    end
                    if avg_I <= -I_1C*(1-tolC) && avg_I >= -I_1C*(1+tolC) && dur >= 50 && dur <= 70
                        disch_n1C_idx(end+1) = k;
                    end
                end
                
                % Charge SOP
                if ~isempty(charge_n1C_idx)
                    valid_idx = max(1, charge_n1C_idx(1)-1) : charge_n1C_idx(end);
                    [t_all, I_all, V_all, SOC_all] = extract_basic_window(pdata, valid_idx, 0.0, Q_0);
                    [meas_soc, P_res] = calc_sop_from_window(t_all, I_all, V_all, SOC_all, I_1C, tolC, V_max_sop, V_min_sop, true, soc_grid, curr_V_ocv_chg_avg);
                    if ~isempty(meas_soc) && length(meas_soc) > 1
                        [meas_soc_u, idx_u] = unique(meas_soc);
                        P_res_u = P_res(idx_u, :);
                        
                        % Ensure strictly monotonically increasing for interp1
                        if length(meas_soc_u) > 1 && meas_soc_u(1) > meas_soc_u(end)
                            meas_soc_u = flipud(meas_soc_u);
                            P_res_u = flipud(P_res_u);
                        end
                        
                        SOP_Data.(ch).(cyc_str).P2s_chg = interp1(meas_soc_u, P_res_u(:,1), soc_mid_chg, 'linear');
                        SOP_Data.(ch).(cyc_str).P10s_chg = interp1(meas_soc_u, P_res_u(:,2), soc_mid_chg, 'linear');
                        SOP_Data.(ch).(cyc_str).P30s_chg = interp1(meas_soc_u, P_res_u(:,3), soc_mid_chg, 'linear');
                        SOP_Data.(ch).(cyc_str).P60s_chg = interp1(meas_soc_u, P_res_u(:,4), soc_mid_chg, 'linear');
                    end
                end
                
                % Discharge SOP
                if ~isempty(disch_n1C_idx)
                    valid_idx = max(1, disch_n1C_idx(1)-1) : disch_n1C_idx(end);
                    [t_all, I_all, V_all, SOC_all] = extract_basic_window(pdata, valid_idx, 1.0, Q_0);
                    [meas_soc, P_res] = calc_sop_from_window(t_all, I_all, V_all, SOC_all, I_1C, tolC, V_max_sop, V_min_sop, false, soc_grid, curr_V_ocv_dch_avg);
                    if ~isempty(meas_soc) && length(meas_soc) > 1
                        [meas_soc_u, idx_u] = unique(meas_soc);
                        P_res_u = P_res(idx_u, :);
                        
                        % Ensure strictly monotonically increasing for interp1
                        if length(meas_soc_u) > 1 && meas_soc_u(1) > meas_soc_u(end)
                            meas_soc_u = flipud(meas_soc_u);
                            P_res_u = flipud(P_res_u);
                        end
                        
                        SOP_Data.(ch).(cyc_str).P2s_dchg = interp1(meas_soc_u, P_res_u(:,1), soc_mid_dch, 'linear');
                        SOP_Data.(ch).(cyc_str).P10s_dchg = interp1(meas_soc_u, P_res_u(:,2), soc_mid_dch, 'linear');
                        SOP_Data.(ch).(cyc_str).P30s_dchg = interp1(meas_soc_u, P_res_u(:,3), soc_mid_dch, 'linear');
                        SOP_Data.(ch).(cyc_str).P60s_dchg = interp1(meas_soc_u, P_res_u(:,4), soc_mid_dch, 'linear');
                    end
                end
                
                % Flatten for table
                current_sop_flat = [SOP_Data.(ch).(cyc_str).P2s_chg, ...
                                    SOP_Data.(ch).(cyc_str).P10s_chg, ...
                                    SOP_Data.(ch).(cyc_str).P30s_chg, ...
                                    SOP_Data.(ch).(cyc_str).P60s_chg, ...
                                    SOP_Data.(ch).(cyc_str).P2s_dchg, ...
                                    SOP_Data.(ch).(cyc_str).P10s_dchg, ...
                                    SOP_Data.(ch).(cyc_str).P30s_dchg, ...
                                    SOP_Data.(ch).(cyc_str).P60s_dchg];
            end
        end

        % --- 2. Extract Labels for each C-rate ---
        for r = 1:length(target_crates)
            crate_label = target_crates{r};
            crate_val = target_crates_val(r);
            
            % SOH: always uses Static capacity (row=1), C-rate independent
            Q_static_bol     = get_max_q(cap_data.Q, 1, idx_cyc0);  % BOL Static (row=1, cyc0)
            Q_static_current = get_max_q(cap_data.Q, 1, c);          % Current Static (row=1)
            if isnan(Q_static_current) || isnan(Q_static_bol)
                SOH = NaN;
            else
                SOH = (Q_static_current / Q_static_bol) * 100;
            end
            
            fresh_str = 'cyc0';
            if isfield(RPT_VQ_grid, cyc_str) && isfield(RPT_VQ_grid.(cyc_str), ch) && isfield(RPT_VQ_grid.(cyc_str).(ch), 'OCV_charge') && ...
               isfield(RPT_VQ_grid, fresh_str) && isfield(RPT_VQ_grid.(fresh_str), ch) && isfield(RPT_VQ_grid.(fresh_str).(ch), 'OCV_charge')
                curr_ocv_chg = RPT_VQ_grid.(cyc_str).(ch).OCV_charge;
                fresh_ocv_chg = RPT_VQ_grid.(fresh_str).(ch).OCV_charge;
                [LLI, LAM] = analyze_ica_aging(curr_ocv_chg, fresh_ocv_chg, Q_0);
            else
                LLI = NaN; LAM = NaN;
            end
            
            if LLI < 0, LLI = 0; end
            if LAM < 0, LAM = 0; end
            
            % SOP: P10s max across 5 discharge segments (peak power capability)
            SOP_dch_10s = max(SOP_Data.(ch).(cyc_str).P10s_dchg);
            
            cnt = cnt + 1;
            y_CellID{cnt,1} = ch;
            y_Cycle(cnt,1) = cyc_num;
            y_CrateLabel{cnt,1} = crate_label;
            y_CrateNum(cnt,1) = crate_val;
            y_Labels(cnt,:) = [SOH, LLI, LAM, SOP_dch_10s];
        end
    end
end

% Finalize Table — 4 labels only (segment detail in SOP_Data struct)
label_names = {'SOH', 'LLI', 'LAM', 'SOP_dch_10s'};

LabelTable_ver02 = table(y_CellID, y_Cycle, y_CrateLabel, y_CrateNum, ...
    y_Labels(:,1), y_Labels(:,2), y_Labels(:,3), y_Labels(:,4), ...
    'VariableNames', [{'CellID', 'Cycle', 'CrateLabel', 'CrateNum'}, label_names]);

save(fullfile(saveDir, 'Label_Matrix_ver02.mat'), 'LabelTable_ver02', 'label_names');
try
    writetable(LabelTable_ver02, fullfile(saveDir, 'Label_Matrix_ver02.xlsx'));
catch ME
    warning('Could not write to Excel file. (.mat file saved successfully)');
end
fprintf('ver02 Labels Extracted and Saved to: %s\\Label_Matrix_ver02.mat\n', saveDir);
fprintf('  Labels: %s\n', strjoin(label_names, ', '));

% Save Struct
save(fullfile(saveDir, 'Label_Struct_ver02.mat'), 'SOP_Data');

%% ========================================================================
% Helper Functions
% ========================================================================
function val = get_max_q(Q_cell, r, c)
    val = NaN;
    try
        item = Q_cell{r, c};
        if ~isempty(item)
            val = max(item);
        end
    catch
        % ignore
    end
end

function [lli, lam] = analyze_ica_aging(curr_ocv_struct, fresh_ocv_struct, Q_rated)
    if nargin < 3, Q_rated = Q_0; end  % Must be passed from caller
    
    if isempty(curr_ocv_struct) || isempty(fresh_ocv_struct)
        lli = NaN;
        lam = NaN;
        return;
    end
    roi_min = 3.40;
    roi_max = 4.00; 
    
    % Get the FRESH peak first (unbiased, absolute max is fine initially)
    [peak_V_fresh, peak_H_fresh] = get_main_peak(fresh_ocv_struct, roi_min, roi_max, NaN);
    
    % Get the AGED peak, TRACKING the fresh peak voltage to prevent U-shaped drift
    [peak_V_aged, peak_H_aged] = get_main_peak(curr_ocv_struct, roi_min, roi_max, peak_V_fresh);
    
    if isnan(peak_V_fresh) || isnan(peak_V_aged)
        lli = NaN;
        lam = NaN;
    else
        % LLI (%): Horizontal peak shift (Dubarry et al.)
        % = delta_V * (dQ/dV)|_fresh / Q_rated * 100
        lli = (abs(peak_V_aged - peak_V_fresh) * peak_H_fresh / Q_rated) * 100;
        
        % LAM (%): Peak Height Ratio
        % Vertical shrinkage of main dQ/dV peak = loss of active material
        lam = (1 - peak_H_aged / peak_H_fresh) * 100;
    end
end

function [pk_V, pk_H] = get_main_peak(data_struct, minV, maxV, target_V)
    % target_V: Optional initial tracking point. If NaN, finds absolute max.
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
        % No clear peaks, fallback to absolute max in ROI
        [pk_H, idx] = max(dH_roi); 
        pk_V = V_roi(idx);
    else
        % We have peaks. Decide if we strictly track 'target_V' or just take 'max'
        if nargin >= 4 && ~isnan(target_V)
            % TRACKING MODE: Find the peak closest to target_V
            peak_voltages = V_roi(locs);
            [~, closest_idx] = min(abs(peak_voltages - target_V));
            pk_V = peak_voltages(closest_idx);
            pk_H = pks(closest_idx);
        else
            % FRESH MODE: Simply pick the highest peak (Absolute Max)
            [pk_H, idx] = max(pks);
            pk_V = V_roi(locs(idx));
        end
    end
end

function [t, I, V, SOC] = extract_basic_window(pdata, valid_idx, initial_soc, C_nom)
    pdata_filtered = pdata(valid_idx);
    I = vertcat(pdata_filtered.I);
    V = vertcat(pdata_filtered.V);
    if isfield(pdata, 'steptime_double') && ~isempty(pdata(1).steptime_double)
        steptime_concat = vertcat(pdata_filtered.steptime_double);
        steptime_concat(1) = 0;
        t = zeros(size(steptime_concat));
        for i=2:length(t)
            dt_val = steptime_concat(i) - steptime_concat(i-1);
            if dt_val < 0, t(i) = t(i-1) + 0.1; else, t(i) = t(i-1) + dt_val; end
        end
    else
        t = vertcat(pdata_filtered.t);
        t = t - t(1);
    end
    SOC = initial_soc + cumtrapz(t, I) / (3600 * C_nom);
end

function [meas_soc, P_res] = calc_sop_from_window(t, I, V, SOC, I_1C, tolC, V_max, V_min, is_charge_target, soc_grid, ocv_avg)
    meas_soc = []; P_res = [];
    threshold = abs(I_1C) * 0.5;
    is_active = abs(I) >= threshold;
    edges = find(diff([0; is_active]) == 1);
    
    for k = 1:length(edges)
        idx_active = edges(k);
        idx_start = idx_active;
        while idx_start > 1 && abs(I(idx_start)) > 1.0
            idx_start = idx_start - 1;
        end
        idx_rest = max(1, idx_start);
        
        idx_2s = find(t >= t(idx_active) + 2, 1, 'first');
        idx_10s = find(t >= t(idx_active) + 10, 1, 'first');
        idx_30s = find(t >= t(idx_active) + 30, 1, 'first');
        idx_60s = find(t >= t(idx_active) + 59, 1, 'first');
        
        if ~isempty(idx_60s) && ~isempty(idx_30s) && ~isempty(idx_10s) && ~isempty(idx_2s)
            is_charge_actual = I(idx_60s) > 0;
            if is_charge_actual == is_charge_target
                if (is_charge_actual && I(idx_60s) >= I_1C*(1-tolC)) || (~is_charge_actual && I(idx_60s) <= -I_1C*(1-tolC))
                    V_start = V(idx_rest);
                    I_start = I(idx_rest);
                    
                    R2 = abs((V(idx_2s) - V_start) / (I(idx_2s) - I_start));
                    R10 = abs((V(idx_10s) - V_start) / (I(idx_10s) - I_start));
                    R30 = abs((V(idx_30s) - V_start) / (I(idx_30s) - I_start));
                    R60 = abs((V(idx_60s) - V_start) / (I(idx_60s) - I_start));
                    
                    soc_pulse = SOC(idx_rest);
                    if soc_pulse > 1.0, soc_pulse = 1.0; elseif soc_pulse < 0.0, soc_pulse = 0.0; end
                    
                    ocv_pulse = interp1(soc_grid, ocv_avg, soc_pulse, 'linear');
                    
                    if is_charge_actual
                        P2 = V_max * (V_max - ocv_pulse) / R2;
                        P10 = V_max * (V_max - ocv_pulse) / R10;
                        P30 = V_max * (V_max - ocv_pulse) / R30;
                        P60 = V_max * (V_max - ocv_pulse) / R60;
                    else
                        P2 = V_min * (ocv_pulse - V_min) / R2;
                        P10 = V_min * (ocv_pulse - V_min) / R10;
                        P30 = V_min * (ocv_pulse - V_min) / R30;
                        P60 = V_min * (ocv_pulse - V_min) / R60;
                    end
                    meas_soc(end+1, 1) = soc_pulse;
                    P_res(end+1, :) = [P2, P10, P30, P60];
                end
            end
        end
    end
end

RPT_Visualize_ver02