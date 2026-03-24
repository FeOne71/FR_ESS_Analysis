%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Feature Extractor (ver02 Final)
% - 21 Physics/Shape Features for Battery Diagnosis
% - Data Source: V_grid / Q (0.001V pre-resampled grid)
% - Peak extraction: findpeaks on full dQ/dV curve (no window)
% - Features:
%   [1-5]   dQ_chg_S1~5    : Charge segmented capacity
%   [6-10]  dQ_dch_S1~5    : Discharge segmented capacity
%   [11]    Peak_H_chg     : dominant peak |dQ/dV| in charge (findpeaks)
%   [12]    Peak_H_dch     : dominant peak |dQ/dV| in discharge (findpeaks)
%   [13]    Peak_A_chg     : Total |trapz(dQ/dV)| in charge
%   [14]    Peak_A_dch     : Total |trapz(dQ/dV)| in discharge
%   [15]    Peak_Pos_chg   : V at dominant peak in charge
%   [16]    Peak_Pos_dch   : V at dominant peak in discharge
%   [17]    Energy_dch     : Discharge energy (Wh) within segment range
%   [18]    C_eff_chg      : I_chg_avg / Q_0  (charge rate)
%   [19]    C_eff_dch      : I_dch_avg / Q_0  (discharge rate)
%   [20]    T_chg_avg      : Average temperature during charge
%   [21]    T_dch_avg      : Average temperature during discharge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
warning off;

%% ========================================================================
% Section 1: Configuration & Data Loading
% ========================================================================
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_rpt_vq_mat = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');
path_static_mat = fullfile(baseDir, 'Capacity_Trend_Figures', 'Capacity_Data_Static.mat');
path_master_ruler = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', ...
    'ExperimentalData', 'FeatureEngineering', 'MasterRulers_v3.mat');

currentScriptPath = mfilename('fullpath');
[ver02Dir, ~, ~] = fileparts(currentScriptPath);
saveDir = fullfile(ver02Dir, 'RPT_FeatureLabelExtractor');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

fprintf('Loading Data...\n');
load(path_rpt_vq_mat, 'RPT_VQ_grid');
load(path_static_mat, 'allChannelsCapacity');
load(path_master_ruler, 'MasterRulers');

channels      = fieldnames(allChannelsCapacity);
target_crates = {'c01', 'c05', 'c1', 'c2', 'c3'};
target_crates_val = [0.1, 0.5, 1, 2, 3];
num_segments  = 5;

% Peak window = segment boundary full range (phase transition tracking)
% Charge: [VR_chg(1), VR_chg(end)]
% Discharge: [min(VR_dch), max(VR_dch)]
ma_window = 21; % MovMean window for dQ/dV smoothing

%% ========================================================================
% Section 2: Feature Extraction Loop
% ========================================================================
fprintf('\nExtracting 21 Physics Features...\n');
data_CellID   = {};
data_Cycle    = [];
data_CrateLabel = {};
data_CrateNum = [];
data_X        = []; % [dQ_chg(5), dQ_dch(5), PkH(2), PkA(2), PkPos(2), Energy_dch(1), C_eff, T_avg(2)]

cnt = 0;
for i = 1:length(channels)
    ch = channels{i};
    if ~isfield(MasterRulers, ch), continue; end

    cap_data  = allChannelsCapacity.(ch);
    idx_cyc0  = find(cap_data.cycles == 0, 1);
    if isempty(idx_cyc0), idx_cyc0 = 1; end

    Q_0 = NaN;
    item = cap_data.Q{1, idx_cyc0};
    if ~isempty(item), Q_0 = max(item); end
    if isnan(Q_0), warning('Q_0 not found for %s, skipping.', ch); continue; end

    VR_chg = MasterRulers.(ch).V_bounds_chg;
    VR_dch = MasterRulers.(ch).V_bounds_dch;

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
            f_dch = [crate_label '_discharge'];

            if ~isfield(ch_data, f_chg) || ~isfield(ch_data, f_dch), continue; end

            s_chg = ch_data.(f_chg);
            s_dch = ch_data.(f_dch);

            %% --- [1] dQ_seg: 10 features (V_grid / Q, 0.001V grid) ---
            % X-axis = Q (strictly monotonic), Y-axis = V
            dQ_chg = nan(1, num_segments);
            dQ_dch = nan(1, num_segments);
            if length(s_chg.Q) > 1
                q_min_c = min(s_chg.Q); q_max_c = max(s_chg.Q);
                q_grid_c = q_min_c : 0.01 : q_max_c;
                [Q_u_c, uid_c] = unique(s_chg.Q);
                V_u_c = s_chg.V_grid(uid_c);
                V_curve_c = interp1(Q_u_c, V_u_c, q_grid_c, 'linear', 'extrap');
                QR_chg = nan(1, length(VR_chg));
                for b = 1:length(VR_chg)
                    [~, idx_v] = min(abs(V_curve_c - VR_chg(b)));
                    QR_chg(b) = q_grid_c(idx_v);
                end
                dQ_chg = abs(diff(QR_chg));
            end
            if length(s_dch.Q) > 1
                q_min_d = min(s_dch.Q); q_max_d = max(s_dch.Q);
                q_grid_d = q_min_d : 0.01 : q_max_d;
                [Q_u_d, uid_d] = unique(s_dch.Q);
                V_u_d = s_dch.V_grid(uid_d);
                V_curve_d = interp1(Q_u_d, V_u_d, q_grid_d, 'linear', 'extrap');
                QR_dch = nan(1, length(VR_dch));
                for b = 1:length(VR_dch)
                    [~, idx_v] = min(abs(V_curve_d - VR_dch(b)));
                    QR_dch(b) = q_grid_d(idx_v);
                end
                dQ_dch = abs(diff(QR_dch));
            end

            %% --- [2] Peak Features: 6 features (V_grid / Q, 0.001V grid) ---
            % No window — findpeaks on full dQ/dV curve (C-rate adaptive)
            [PkH_chg, PkA_chg, PkPos_chg] = extract_peak_features(s_chg, ma_window);
            [PkH_dch, PkA_dch, PkPos_dch] = extract_peak_features(s_dch, ma_window);

            %% --- [3] Energy (Wh) Features: 1 feature (trapz(Q_raw, V_raw)) ---
            % Energy window = full segment boundary range
            Energy_dch = NaN;
            idx_d = s_dch.V_raw >= min(VR_dch) & s_dch.V_raw <= max(VR_dch);
            if sum(idx_d) > 5
                Energy_dch = abs(trapz(s_dch.Q_raw(idx_d), s_dch.V_raw(idx_d)));
            end

            %% --- [4] C_eff_chg / C_eff_dch: 2 features (I_raw) ---
            C_eff_chg = mean(abs(s_chg.I_raw)) / Q_0;
            C_eff_dch = mean(abs(s_dch.I_raw)) / Q_0;

            %% --- [5] T_chg_avg / T_dch_avg: 2 features (T1_raw) ---
            T_chg_avg = NaN;
            T_dch_avg = NaN;
            if isfield(s_chg, 'T1_raw') && ~isempty(s_chg.T1_raw)
                T_chg_avg = mean(s_chg.T1_raw(:), 'omitnan');
            end
            if isfield(s_dch, 'T1_raw') && ~isempty(s_dch.T1_raw)
                T_dch_avg = mean(s_dch.T1_raw(:), 'omitnan');
            end

            %% --- Collect row (20 features) ---
            cnt = cnt + 1;
            data_CellID{cnt,1}   = ch;
            data_Cycle(cnt,1)    = cyc_num;
            data_CrateLabel{cnt,1} = crate_label;
            data_CrateNum(cnt,1) = crate_val;
            data_X(cnt,:) = [dQ_chg, dQ_dch, ...
                             PkH_chg, PkH_dch, PkA_chg, PkA_dch, PkPos_chg, PkPos_dch, ...
                             Energy_dch, C_eff_chg, C_eff_dch, T_chg_avg, T_dch_avg];
        end
    end
end

%% ========================================================================
% Section 3: Finalize & Save
% ========================================================================
feature_names = {'dQ_chg_S1','dQ_chg_S2','dQ_chg_S3','dQ_chg_S4','dQ_chg_S5', ...
                 'dQ_dch_S1','dQ_dch_S2','dQ_dch_S3','dQ_dch_S4','dQ_dch_S5', ...
                 'Peak_H_chg','Peak_H_dch', ...
                 'Peak_A_chg','Peak_A_dch', ...
                 'Peak_Pos_chg','Peak_Pos_dch', ...
                 'Energy_dch', ...
                 'C_eff_chg','C_eff_dch','T_chg_avg','T_dch_avg'};

FeatureTable_ver02 = table(data_CellID, data_Cycle, data_CrateLabel, data_CrateNum, ...
    'VariableNames', {'CellID','Cycle','CrateLabel','CrateNum'});

for f = 1:length(feature_names)
    FeatureTable_ver02.(feature_names{f}) = data_X(:, f);
end

save(fullfile(saveDir, 'Feature_Matrix_ver02.mat'), 'FeatureTable_ver02', 'feature_names');
try
    writetable(FeatureTable_ver02(:, [{'CellID','Cycle','CrateLabel','CrateNum'}, feature_names]), ...
        fullfile(saveDir, 'Feature_Matrix_ver02.xlsx'));
catch
    warning('Could not write Excel file. (.mat saved successfully)');
end

fprintf('\n=== Feature Extraction Complete ===\n');
fprintf('  Total samples : %d\n', cnt);
fprintf('  Feature width : %d columns\n', size(data_X, 2));
fprintf('  Saved to      : %s\n', saveDir);

%% ========================================================================
% Helper: Extract Peak Features from V_grid / Q data (0.001V grid)
% ========================================================================
function [pk_height, pk_area, pk_pos] = extract_peak_features(data_struct, ma_win)
    pk_height = NaN; pk_area = NaN; pk_pos = NaN;

    V_grid = data_struct.V_grid;
    Q      = data_struct.Q;
    if numel(V_grid) < 10, return; end

    [V_u, uid] = unique(V_grid);
    Q_u = Q(uid);
    if V_u(1) > V_u(end)   % discharge: sort ascending
        V_u = flipud(V_u); Q_u = flipud(Q_u);
    end

    % dQ/dV (smoothed, movmean)
    dV_all = gradient(V_u);
    dQ_all = gradient(Q_u);
    dV_all(dV_all == 0) = NaN;
    dQdV_raw = dQ_all ./ dV_all;
    dQdV_raw(isinf(dQdV_raw) | isnan(dQdV_raw)) = 0;
    dQdV_filt = movmean(dQdV_raw, ma_win);

    % Peak_A: trapz over full curve
    pk_area = abs(trapz(V_u, dQdV_filt));

    % Peak_H & Peak_Pos: dominant peak via findpeaks (C-rate adaptive)
    dQdV_abs = abs(dQdV_filt);
    [pks, locs] = findpeaks(dQdV_abs, 'SortStr', 'descend', 'NPeaks', 1);
    if ~isempty(pks)
        pk_height = pks(1);
        pk_pos = V_u(locs(1));
    end
end

% RPT_LabelExtractor   % (commented out: run separately if needed)
% RPT_Visualize_ver02