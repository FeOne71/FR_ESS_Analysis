%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT_FeatureExtractor_v5_Charge.m
% - NCM Equi-capacity 11 Segments (Charge Only)
% - No Discharge
% - Pooled Data approach (No ECM correction by default, just raw V)
% - Extracts 11 dQ segments, Peak H/A/Pos, C_eff, T_avg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
warning off;

%% ========================================================================
% Section 1: Configuration & Data Loading
% ========================================================================
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_rpt_vq_mat  = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');
path_static_mat  = fullfile(baseDir, 'Capacity_Trend_Figures', 'Capacity_Data_Static.mat');

currentScriptPath = mfilename('fullpath');
[scriptDir, ~, ~] = fileparts(currentScriptPath);
saveDir = fullfile(scriptDir, 'RPT_FeatureLabelExtractor_v5_Charge');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

fprintf('Loading Data...\n');
load(path_rpt_vq_mat,  'RPT_VQ_grid');
load(path_static_mat,  'allChannelsCapacity');

channels       = fieldnames(allChannelsCapacity);
target_crates  = {'c01', 'c05', 'c1', 'c2', 'c3'};
target_crates_val = [0.1, 0.5, 1, 2, 3];

% NEW: Equi-capacity NCM boundaries (Phase 5)
VR_chg = [ 3.40, 3.52, 3.59, 3.64, 3.67, 3.70, 3.75, 3.84, 3.92, 4.01, 4.10, 4.20 ];
num_segments   = length(VR_chg) - 1; % 11
ma_window      = 21;

%% ========================================================================
% Section 2: Feature Extraction Loop
% ========================================================================
fprintf('\nExtracting v5 Charge-Only Features (Equi-Capacity Segments)...\n');
data_CellID    = {};
data_Cycle     = [];
data_CrateLabel = {};
data_CrateNum  = [];
data_X         = [];

cnt = 0;
for i = 1:length(channels)
    ch = channels{i};
    
    cap_data  = allChannelsCapacity.(ch);
    idx_cyc0  = find(cap_data.cycles == 0, 1);
    if isempty(idx_cyc0), idx_cyc0 = 1; end

    Q_0 = NaN;
    item = cap_data.Q{1, idx_cyc0};
    if ~isempty(item), Q_0 = max(item); end
    if isnan(Q_0), warning('Q_0 not found for %s, skipping.', ch); continue; end

    cyc_fields = fieldnames(RPT_VQ_grid);
    for c = 1:length(cyc_fields)
        cyc_key = cyc_fields{c};
        cyc_num = sscanf(cyc_key, 'cyc%d');
        if ~isfield(RPT_VQ_grid.(cyc_key), ch), continue; end

        ch_data = RPT_VQ_grid.(cyc_key).(ch);

        for r = 1:length(target_crates)
            crate_label = target_crates{r};
            crate_val   = target_crates_val(r);
            f_chg = [crate_label '_charge'];
            
            if ~isfield(ch_data, f_chg), continue; end

            s_chg = ch_data.(f_chg);
            
            % Use Raw Voltage for extraction (no ECM correction to match CP field nature)
            V_chg = s_chg.V_raw(:);
            Q_chg = s_chg.Q_raw(:);
            I_chg = s_chg.I_raw(:);
            
            s_chg_corr = build_vgrid(V_chg, Q_chg, true);

            %% --- [1] dQ_seg: 11 features ---
            dQ_chg = nan(1, num_segments);
            if length(s_chg_corr.V_grid) > 1
                q_min_c = min(s_chg_corr.Q); q_max_c = max(s_chg_corr.Q);
                q_grid_c = q_min_c : 0.001 : q_max_c; % High resolution Q grid
                
                [Q_u_c, uid_c] = unique(s_chg_corr.Q);
                V_u_c = s_chg_corr.V_grid(uid_c);
                V_curve_c = interp1(Q_u_c, V_u_c, q_grid_c, 'linear', 'extrap');
                
                QR_chg = nan(1, length(VR_chg));
                for b = 1:length(VR_chg)
                    [~, idx_v] = min(abs(V_curve_c - VR_chg(b)));
                    QR_chg(b) = q_grid_c(idx_v);
                end
                dQ_chg = abs(diff(QR_chg));
            end

            %% --- [2] Peak Features: 3 features ---
            [PkH_chg, PkA_chg, PkPos_chg] = extract_peak_features(s_chg_corr, ma_window);

            %% --- [3] C_eff ---
            C_eff_chg = mean(abs(I_chg)) / Q_0;

            %% --- [4] T_avg ---
            T_chg_avg = NaN;
            if isfield(s_chg,'T1_raw') && ~isempty(s_chg.T1_raw)
                T_chg_avg = mean(s_chg.T1_raw(:), 'omitnan');
            end

            %% --- Collect row ---
            cnt = cnt + 1;
            data_CellID{cnt,1}     = ch;
            data_Cycle(cnt,1)      = cyc_num;
            data_CrateLabel{cnt,1} = crate_label;
            data_CrateNum(cnt,1)   = crate_val;
            data_X(cnt,:) = [dQ_chg, PkH_chg, PkA_chg, PkPos_chg, C_eff_chg, T_chg_avg];
        end
    end
    fprintf('  Done: %s\n', ch);
end

%% ========================================================================
% Section 3: Save
% ========================================================================
feature_names = cell(1, num_segments);
for k = 1:num_segments
    feature_names{k} = sprintf('dQ_chg_S%d', k);
end

feature_names = [feature_names, {'Peak_H_chg', 'Peak_A_chg', 'Peak_Pos_chg', 'C_eff_chg', 'T_chg_avg'}];

FeatureTable_v5_Charge = table(data_CellID, data_Cycle, data_CrateLabel, data_CrateNum, ...
    'VariableNames', {'CellID','Cycle','CrateLabel','CrateNum'});
for f = 1:length(feature_names)
    FeatureTable_v5_Charge.(feature_names{f}) = data_X(:, f);
end

save(fullfile(saveDir, 'Feature_Matrix_v5_Charge.mat'), 'FeatureTable_v5_Charge', 'feature_names');
fprintf('\n=== v5 Charge-Only Feature Extraction Complete ===\n');
fprintf('  Total samples : %d\n', cnt);
fprintf('  Saved to      : %s\n', saveDir);

%% ========================================================================
% Helper Functions
% ========================================================================
function s_out = build_vgrid(V_raw, Q_raw, is_charge)
    Q_raw = cumtrapz(1:length(Q_raw), abs(diff([0; Q_raw])));  % make monotonic Q
    Q_raw = Q_raw(:);

    [V_u, u] = unique(V_raw(:), 'stable');
    Q_u = Q_raw(u);
    mono = true(size(V_u));
    if is_charge
        for ii = 2:length(V_u), if V_u(ii) <= V_u(ii-1), mono(ii) = false; end, end
    else
        for ii = 2:length(V_u), if V_u(ii) >= V_u(ii-1), mono(ii) = false; end, end
    end
    V_u = V_u(mono); Q_u = Q_u(mono);
    if ~is_charge, V_u = flipud(V_u); Q_u = flipud(Q_u); end  % ascending for interp

    Vg = (min(V_u):0.001:max(V_u))';
    Qg = interp1(V_u, Q_u, Vg, 'linear');
    good = ~isnan(Qg);
    Vg = Vg(good); Qg = Qg(good);

    s_out.V_grid = Vg;
    s_out.Q = Qg;
end

function [pk_height, pk_area, pk_pos] = extract_peak_features(s, ma_win)
    pk_height = NaN; pk_area = NaN; pk_pos = NaN;
    V_grid = s.V_grid; Q = s.Q;
    if numel(V_grid) < 10, return; end

    [V_u, uid] = unique(V_grid);
    Q_u = Q(uid);
    if V_u(1) > V_u(end), V_u = flipud(V_u); Q_u = flipud(Q_u); end

    dV = gradient(V_u); dQ = gradient(Q_u);
    dV(dV==0) = NaN;
    dQdV = dQ ./ dV; dQdV(isinf(dQdV)|isnan(dQdV)) = 0;
    dQdV_filt = movmean(dQdV, ma_win);

    pk_area = abs(trapz(V_u, dQdV_filt));
    [pks, locs] = findpeaks(abs(dQdV_filt), 'SortStr', 'descend', 'NPeaks', 1);
    if ~isempty(pks)
        pk_height = pks(1);
        pk_pos = V_u(locs(1));
    end
end
