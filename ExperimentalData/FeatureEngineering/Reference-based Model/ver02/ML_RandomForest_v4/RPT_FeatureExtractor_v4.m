%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT_FeatureExtractor_v4.m
% - 2RC ECM 보정 적용 버전
% - V_raw에서 과전압(I×R₀ + V_RC1 + V_RC2) 제거 → OCV_est
% - OCV_est 기반 dQ/dV 계산
% - 동일 21개 피처 추출
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
warning off;

%% ========================================================================
% Section 0: ECM Parameters Load
% ========================================================================
ecm_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4\ECM_2RC_AllChannels_cyc0.mat';
fprintf('Loading ECM 2RC parameters...\n');
ecm_data = load(ecm_path, 'All_ECM');
All_ECM = ecm_data.All_ECM;

%% ========================================================================
% Section 1: Configuration & Data Loading
% ========================================================================
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_rpt_vq_mat  = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');
path_static_mat  = fullfile(baseDir, 'Capacity_Trend_Figures', 'Capacity_Data_Static.mat');
path_master_ruler = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', ...
    'ExperimentalData', 'FeatureEngineering', 'MasterRulers_v3.mat');

currentScriptPath = mfilename('fullpath');
[scriptDir, ~, ~] = fileparts(currentScriptPath);
saveDir = fullfile(scriptDir, 'RPT_FeatureLabelExtractor_v4');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

fprintf('Loading Data...\n');
load(path_rpt_vq_mat,  'RPT_VQ_grid');
load(path_static_mat,  'allChannelsCapacity');
load(path_master_ruler, 'MasterRulers');

channels       = fieldnames(allChannelsCapacity);
target_crates  = {'c01', 'c05', 'c1', 'c2', 'c3'};
target_crates_val = [0.1, 0.5, 1, 2, 3];
num_segments   = 5;
ma_window      = 21;

%% ========================================================================
% Section 2: Feature Extraction Loop
% ========================================================================
fprintf('\nExtracting v4 ECM-Corrected Features...\n');
data_CellID    = {};
data_Cycle     = [];
data_CrateLabel = {};
data_CrateNum  = [];
data_X         = [];

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

    % Load ECM params for this channel (cyc0) - charge and discharge separately
    if isfield(All_ECM, ch) && isfield(All_ECM.(ch), 'discharge')
        ecm_dch = All_ECM.(ch).discharge;
        ecm_chg = All_ECM.(ch).charge;
    else
        warning('No ECM params for %s, skipping.', ch);
        continue;
    end

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

            %% --- ECM Correction ---
            V_ocv_chg = apply_ecm_correction(s_chg.V_raw(:), s_chg.I_raw(:), s_chg.t_raw(:), ecm_chg, Q_0, true);
            V_ocv_dch = apply_ecm_correction(s_dch.V_raw(:), s_dch.I_raw(:), s_dch.t_raw(:), ecm_dch, Q_0, false);

            % Build corrected V_grid using V-grid method on OCV_est
            s_chg_corr = build_vgrid(V_ocv_chg, s_chg.Q_raw(:), true);
            s_dch_corr = build_vgrid(V_ocv_dch, s_dch.Q_raw(:), false);

            %% --- [1] dQ_seg: 10 features ---
            dQ_chg = nan(1, num_segments);
            dQ_dch = nan(1, num_segments);
            if length(s_chg_corr.V_grid) > 1
                q_min_c = min(s_chg_corr.Q); q_max_c = max(s_chg_corr.Q);
                q_grid_c = q_min_c : 0.01 : q_max_c;
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
            if length(s_dch_corr.V_grid) > 1
                q_min_d = min(s_dch_corr.Q); q_max_d = max(s_dch_corr.Q);
                q_grid_d = q_min_d : 0.01 : q_max_d;
                [Q_u_d, uid_d] = unique(s_dch_corr.Q);
                V_u_d = s_dch_corr.V_grid(uid_d);
                V_curve_d = interp1(Q_u_d, V_u_d, q_grid_d, 'linear', 'extrap');
                QR_dch = nan(1, length(VR_dch));
                for b = 1:length(VR_dch)
                    [~, idx_v] = min(abs(V_curve_d - VR_dch(b)));
                    QR_dch(b) = q_grid_d(idx_v);
                end
                dQ_dch = abs(diff(QR_dch));
            end

            %% --- [2] Peak Features: 6 features ---
            [PkH_chg, PkA_chg, PkPos_chg] = extract_peak_features(s_chg_corr, ma_window);
            [PkH_dch, PkA_dch, PkPos_dch] = extract_peak_features(s_dch_corr, ma_window);

            %% --- [3] Energy: 1 feature (from original V_raw) ---
            Energy_dch = NaN;
            idx_d = s_dch.V_raw >= min(VR_dch) & s_dch.V_raw <= max(VR_dch);
            if sum(idx_d) > 5
                Energy_dch = abs(trapz(s_dch.Q_raw(idx_d), s_dch.V_raw(idx_d)));
            end

            %% --- [4] C_eff ---
            C_eff_chg = mean(abs(s_chg.I_raw)) / Q_0;
            C_eff_dch = mean(abs(s_dch.I_raw)) / Q_0;

            %% --- [5] T_avg ---
            T_chg_avg = NaN; T_dch_avg = NaN;
            if isfield(s_chg,'T1_raw') && ~isempty(s_chg.T1_raw)
                T_chg_avg = mean(s_chg.T1_raw(:), 'omitnan');
            end
            if isfield(s_dch,'T1_raw') && ~isempty(s_dch.T1_raw)
                T_dch_avg = mean(s_dch.T1_raw(:), 'omitnan');
            end

            %% --- Collect row ---
            cnt = cnt + 1;
            data_CellID{cnt,1}     = ch;
            data_Cycle(cnt,1)      = cyc_num;
            data_CrateLabel{cnt,1} = crate_label;
            data_CrateNum(cnt,1)   = crate_val;
            data_X(cnt,:) = [dQ_chg, dQ_dch, ...
                             PkH_chg, PkH_dch, PkA_chg, PkA_dch, PkPos_chg, PkPos_dch, ...
                             Energy_dch, C_eff_chg, C_eff_dch, T_chg_avg, T_dch_avg];
        end
    end
    fprintf('  Done: %s\n', ch);
end

%% ========================================================================
% Section 3: Save
% ========================================================================
feature_names = {'dQ_chg_S1','dQ_chg_S2','dQ_chg_S3','dQ_chg_S4','dQ_chg_S5', ...
                 'dQ_dch_S1','dQ_dch_S2','dQ_dch_S3','dQ_dch_S4','dQ_dch_S5', ...
                 'Peak_H_chg','Peak_H_dch', ...
                 'Peak_A_chg','Peak_A_dch', ...
                 'Peak_Pos_chg','Peak_Pos_dch', ...
                 'Energy_dch', ...
                 'C_eff_chg','C_eff_dch','T_chg_avg','T_dch_avg'};

FeatureTable_v4 = table(data_CellID, data_Cycle, data_CrateLabel, data_CrateNum, ...
    'VariableNames', {'CellID','Cycle','CrateLabel','CrateNum'});
for f = 1:length(feature_names)
    FeatureTable_v4.(feature_names{f}) = data_X(:, f);
end

save(fullfile(saveDir, 'Feature_Matrix_v4.mat'), 'FeatureTable_v4', 'feature_names');
fprintf('\n=== v4 ECM Feature Extraction Complete ===\n');
fprintf('  Total samples : %d\n', cnt);
fprintf('  Saved to      : %s\n', saveDir);

%% ========================================================================
% Helper: ECM Correction (Steady-State)
% CC 충방전에서 RC 정상상태 가정
% V_ocv = V_terminal - I(t) × R_total(SOC)
% R_total = R₀ + R₁ + R₂  (SOC 10~95% 범위만 사용)
% ========================================================================
function V_ocv = apply_ecm_correction(V_raw, I_raw, t_raw, ecm_ch, Q_0, is_charge)
    if isa(t_raw, 'duration'), t_raw = seconds(t_raw); end
    t_raw = double(t_raw(:)); V_raw = double(V_raw(:)); I_raw = double(I_raw(:));

    % SOC via Coulomb counting
    I_abs = abs(I_raw);
    Q_cum = cumtrapz(t_raw, I_abs) / 3600;
    if is_charge, SOC = Q_cum / Q_0 * 100;
    else, SOC = 100 - Q_cum / Q_0 * 100; end
    SOC = max(0, min(100, SOC));

    % Build R_total LUT (SOC 10~95% only — 극한 SOC 이상치 제거)
    [soc_s, si] = sort(ecm_ch.SOC);
    R_total_s = (ecm_ch.R0(si) + ecm_ch.R1(si) + ecm_ch.R2(si)) / 1e3;  % Ω
    valid = soc_s >= 10 & soc_s <= 95;
    soc_s     = soc_s(valid);
    R_total_s = R_total_s(valid);

    % Clamp SOC & interpolate
    SOC_c   = max(min(soc_s), min(max(soc_s), SOC));
    R_total = interp1(soc_s, R_total_s, SOC_c, 'linear');

    V_ocv = V_raw - I_raw .* R_total;
end


%% Helper: build V_grid from OCV_est and Q (0.001V)
function s_out = build_vgrid(V_ocv, Q_raw, is_charge)
    Q_raw = cumtrapz(1:length(Q_raw), abs(diff([0; Q_raw])));  % make monotonic Q
    Q_raw = Q_raw(:);

    [V_u, u] = unique(V_ocv(:), 'stable');
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

%% Helper: Extract Peak Features
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
