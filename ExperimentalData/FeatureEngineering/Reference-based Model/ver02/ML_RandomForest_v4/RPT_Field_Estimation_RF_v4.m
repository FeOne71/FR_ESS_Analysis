% RPT_Field_Estimation_RF_v4.m
% =====================================================================
% 충전/방전 분리 RF 모델을 사용한 필드 SOH/LLI/LAM 추정
%
% [흐름]
%   1. 학습된 충전/방전 모델 로드
%   2. 필드 데이터에서 충전/방전 세그먼트 추출
%   3. 충전 세그먼트 → 충전 피처 9개 → C_eff_chg 기준 NormStats 보간 → RF_chg
%   4. 방전 세그먼트 → 방전 피처 10개 → C_eff_dch 기준 NormStats 보간 → RF_dch
%   5. 결과 결합: 둘 다 있으면 CV RMSE 기반 가중평균, 한쪽만 있으면 그 결과 사용
% =====================================================================

clear; clc; close all;
warning on;

%% ========================================================================
% Section 1: 모델 로드
%% ========================================================================
fprintf('=== Section 1: Loading Split RF Models (v4) ===\n');
baseDir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
                   'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
v4Dir   = fullfile(baseDir, 'ML_RandomForest_v4');
modelFile = fullfile(v4Dir, 'Result_RF_v4.mat');
path_master_ruler = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', ...
    'ExperimentalData', 'FeatureEngineering', 'MasterRulers_v3.mat');

if ~exist(modelFile, 'file'), error('RF v4 Model not found: %s', modelFile); end
load(modelFile, 'Results_v4');
load(path_master_ruler, 'MasterRulers');

label_names = Results_v4.label_names;

%% ========================================================================
% Section 2: 필드 데이터 로드
%% ========================================================================
fprintf('\n=== Section 2: Loading Field Data (2021~2025) ===\n');
rawDir   = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'Rack_raw2mat');
dataFiles = {
    fullfile(rawDir, 'Old', '2021', '202106', 'Raw_20210603 .mat'), 'Y2021', 'old',  datetime(2021,6,3);
    fullfile(rawDir, 'New', '2023', '202310', 'Raw_20231016.mat'), 'Y2023', 'new',  datetime(2023,10,16);
    fullfile(rawDir, 'New', '2024', '202409', 'Raw_20240909.mat'), 'Y2024', 'new',  datetime(2024,9,9);
    fullfile(rawDir, 'New', '2025', '202507', 'Raw_20250711.mat'), 'Y2025', 'new',  datetime(2024,9,9);
};

min_charge_secs    = [600, 300, 300, 300];
min_discharge_secs = [300, 150, 300, 300];

dt = 1;
fns = fieldnames(MasterRulers);
VR_chg = MasterRulers.(fns{1}).V_bounds_chg;
VR_dch = MasterRulers.(fns{1}).V_bounds_dch;
Np = 2; C_cell_Ah = 64; thr_A = C_cell_Ah * 0.05;
ma_window = 60; Q_0 = 64;

FieldData = struct();
for k = 1:size(dataFiles, 1)
    fpath      = dataFiles{k, 1};
    year_label = dataFiles{k, 2};
    dataType   = dataFiles{k, 3};
    base_date  = dataFiles{k, 4};
    min_chg_sec = min_charge_secs(k);
    min_dch_sec = min_discharge_secs(k);

    if ~exist(fpath, 'file'), continue; end
    fprintf('  Processing %s...\n', year_label);

    S = load(fpath);
    if strcmp(dataType, 'old')
        D = S.Raw.Rack01;
        if isfield(D, 'Time'), t = datetime(D.Time);
        elseif isfield(D, 'Date_Time')
            if isduration(D.Date_Time), t = base_date + D.Date_Time;
            else, t = datetime(D.Date_Time); end
        else, continue; end
        I_rack = D.DCCurrent_A(:);
        V_avg  = D.AverageCV_V(:);
        if isfield(D, 'AverageMT_degC'), T_avg = D.AverageMT_degC(:);
        else, T_avg = nan(size(I_rack)); end
        if isfield(D, 'SOHPct'), raw_soh_bms = D.SOHPct(:);
        else, raw_soh_bms = nan(size(I_rack)); end
    else
        D = S.Raw;
        if isduration(D.Date_Time), t = base_date + D.Date_Time;
        else, t = datetime(D.Date_Time); end
        I_rack = D.DCCurrent(:);
        V_avg  = D.CVavg(:);
        if isfield(D, 'MTavg'), T_avg = D.MTavg(:);
        else, T_avg = nan(size(I_rack)); end
        if isfield(D, 'SOH_BMS'), raw_soh_bms = D.SOH_BMS(:);
        else, raw_soh_bms = nan(size(I_rack)); end
    end

    valid_soh = raw_soh_bms;
    valid_soh(valid_soh <= 0) = NaN;
    soh_bms_day = round(median(valid_soh, 'omitnan'));

    tsec   = seconds(t - t(1));
    I_cell = I_rack / Np;

    % --- 충전/방전 세그먼트 추출 (기존 로직 동일) ---
    chg_mask = I_cell > thr_A;
    dch_mask = I_cell < -thr_A;

    % Longest charge segment
    chg_s = NaN; chg_e = NaN;
    chg_diff = diff([0; chg_mask(:); 0]);
    chg_starts = find(chg_diff == 1);
    chg_ends   = find(chg_diff == -1) - 1;
    if ~isempty(chg_starts)
        chg_lens = chg_ends - chg_starts + 1;
        valid_chg_segs = find(chg_lens * dt >= min_chg_sec);
        if ~isempty(valid_chg_segs)
            [~, mx] = max(chg_lens(valid_chg_segs));
            idx = valid_chg_segs(mx);
            chg_s = chg_starts(idx); chg_e = chg_ends(idx);
        end
    end

    % Longest discharge segment
    dch_s = NaN; dch_e = NaN;
    dch_diff = diff([0; dch_mask(:); 0]);
    dch_starts = find(dch_diff == 1);
    dch_ends   = find(dch_diff == -1) - 1;
    if ~isempty(dch_starts)
        dch_lens = dch_ends - dch_starts + 1;
        valid_dch_segs = find(dch_lens * dt >= min_dch_sec);
        if ~isempty(valid_dch_segs)
            [~, mx] = max(dch_lens(valid_dch_segs));
            idx = valid_dch_segs(mx);
            dch_s = dch_starts(idx); dch_e = dch_ends(idx);
        end
    end

    has_chg = ~isnan(chg_s);
    has_dch = ~isnan(dch_s);

    if has_chg
        FieldData.(year_label).Chg.Time = tsec(chg_s:chg_e);
        FieldData.(year_label).Chg.V = V_avg(chg_s:chg_e);
        FieldData.(year_label).Chg.I = I_cell(chg_s:chg_e);
        FieldData.(year_label).Chg.T = T_avg(chg_s:chg_e);
    end
    if has_dch
        FieldData.(year_label).Dch.Time = tsec(dch_s:dch_e);
        FieldData.(year_label).Dch.V = V_avg(dch_s:dch_e);
        FieldData.(year_label).Dch.I = I_cell(dch_s:dch_e);
        FieldData.(year_label).Dch.T = T_avg(dch_s:dch_e);
    end
    FieldData.(year_label).SOH_BMS = soh_bms_day;
end

%% ========================================================================
% Section 3: 독립 피처 추출 및 추정
%% ========================================================================
fprintf('\n=== Section 3: Split Feature Extraction & Estimation ===\n');
years = fieldnames(FieldData);
Results_Ev = struct();

fallback_priority = [3, 4, 2, 5, 1];
num_segs = length(VR_chg) - 1;

for k = 1:length(years)
    yr = years{k};
    has_chg = isfield(FieldData.(yr), 'Chg');
    has_dch = isfield(FieldData.(yr), 'Dch');

    if ~has_chg && ~has_dch
        fprintf('  %s: No segments found. Skipped.\n', yr);
        continue;
    end

    % ================================================================
    % 3A: 충전 모델 추정
    % ================================================================
    if has_chg
        t_c = FieldData.(yr).Chg.Time(:);
        v_c = FieldData.(yr).Chg.V(:);
        i_c = FieldData.(yr).Chg.I(:);

        i_c_sm = movmean(abs(i_c), ma_window);
        q_c = cumtrapz(t_c, i_c_sm) / 3600;
        v_c_sm = movmean(v_c, ma_window);
        q_c_sm = movmean(q_c, ma_window);

        [v_c_uniq, uc] = unique(v_c_sm, 'stable');
        q_c_uniq = q_c_sm(uc);
        mono_c = true(size(v_c_uniq));
        for ii = 2:length(v_c_uniq)
            if v_c_uniq(ii) <= v_c_uniq(ii-1), mono_c(ii) = false; end
        end
        v_c_uniq = v_c_uniq(mono_c); q_c_uniq = q_c_uniq(mono_c);

        % dQ_chg segments
        dQ_chg_raw = nan(1, num_segs);
        V_min_c = min(v_c); V_max_c = max(v_c);
        valid_chg = false(1, num_segs);
        for s = 1:num_segs
            valid_chg(s) = (VR_chg(s) >= V_min_c - 0.02) && (VR_chg(s+1) <= V_max_c + 0.02);
        end

        if length(v_c_uniq) > 1
            try
                v_grid_c = (VR_chg(1)-0.05):0.001:(VR_chg(end)+0.05);
                Q_interp_c = interp1(v_c_uniq, q_c_uniq, v_grid_c, 'linear', 'extrap');
                Q_sm_c = movmean(Q_interp_c, 30);
                QR_chg = nan(1, length(VR_chg));
                for b = 1:length(VR_chg)
                    [~, mi] = min(abs(v_grid_c - VR_chg(b)));
                    QR_chg(b) = Q_sm_c(mi);
                end
                dQ_chg_raw = abs(diff(QR_chg));
                dQ_chg_raw(~valid_chg) = NaN;
            catch
            end
        end
        dQ_chg = segment_fallback(dQ_chg_raw, valid_chg, fallback_priority);

        % Peak features (charge)
        [PkH_chg, PkA_chg, PkPos_chg] = field_extract_peak(v_c_uniq, q_c_uniq);

        % C_eff_chg
        C_eff_chg = mean(abs(i_c)) / Q_0;

        % Assemble charge feature vector (9 features)
        X_chg_raw = [dQ_chg, PkH_chg, PkA_chg, PkPos_chg, C_eff_chg];
        X_chg_raw(isnan(X_chg_raw)) = 0;

        % Normalize: C_eff_chg 기준 NormStats 보간
        NS_chg = Results_v4.chg.NormStats;
        X_chg_norm = field_normalize(X_chg_raw, C_eff_chg, NS_chg);

        % Predict
        for i = 1:length(label_names)
            lbl = label_names{i};
            Results_Ev.(yr).chg.(lbl) = predict(Results_v4.chg.(lbl).Model, X_chg_norm);
        end
        Results_Ev.(yr).chg.C_eff = C_eff_chg;
        Results_Ev.(yr).chg.X_raw = X_chg_raw;

        fprintf('  %s [CHG] C_eff=%.3fC', yr, C_eff_chg);
        for i = 1:length(label_names)
            fprintf(' | %s=%.2f', label_names{i}, Results_Ev.(yr).chg.(label_names{i}));
        end
        fprintf('\n');
    end

    % ================================================================
    % 3B: 방전 모델 추정
    % ================================================================
    if has_dch
        t_d = FieldData.(yr).Dch.Time(:);
        v_d = FieldData.(yr).Dch.V(:);
        i_d = FieldData.(yr).Dch.I(:);

        i_d_sm = movmean(abs(i_d), ma_window);
        q_d = cumtrapz(t_d, i_d_sm) / 3600;
        v_d_sm = movmean(v_d, ma_window);
        q_d_sm = movmean(q_d, ma_window);

        [v_d_uniq, ud] = unique(v_d_sm, 'stable');
        q_d_uniq = q_d_sm(ud);
        mono_d = true(size(v_d_uniq));
        for ii = 2:length(v_d_uniq)
            if v_d_uniq(ii) >= v_d_uniq(ii-1), mono_d(ii) = false; end
        end
        v_d_uniq = v_d_uniq(mono_d); q_d_uniq = q_d_uniq(mono_d);

        % dQ_dch segments
        dQ_dch_raw = nan(1, num_segs);
        V_min_d = min(v_d); V_max_d = max(v_d);
        valid_dch = false(1, num_segs);
        for s = 1:num_segs
            v_hi_d = VR_dch(s); v_lo_d = VR_dch(s+1);
            valid_dch(s) = (v_hi_d <= V_max_d + 0.02) && (v_lo_d >= V_min_d - 0.02);
        end

        if length(v_d_uniq) > 1
            try
                v_grid_d = (VR_dch(1)+0.05):-0.001:(VR_dch(end)-0.05);
                Q_interp_d = interp1(v_d_uniq, q_d_uniq, v_grid_d, 'linear', 'extrap');
                Q_sm_d = movmean(Q_interp_d, 30);
                QR_dch = nan(1, length(VR_dch));
                for b = 1:length(VR_dch)
                    [~, mi] = min(abs(v_grid_d - VR_dch(b)));
                    QR_dch(b) = Q_sm_d(mi);
                end
                dQ_dch_raw = abs(diff(QR_dch));
                dQ_dch_raw(~valid_dch) = NaN;
            catch
            end
        end
        dQ_dch = segment_fallback(dQ_dch_raw, valid_dch, fallback_priority);

        % Peak features (discharge)
        [PkH_dch, PkA_dch, PkPos_dch] = field_extract_peak(v_d_uniq, q_d_uniq);

        % Energy_dch
        Energy_dch = NaN;
        idx_e = v_d >= min(VR_dch) & v_d <= max(VR_dch);
        if sum(idx_e) > 5
            q_e = cumtrapz(t_d(idx_e), abs(i_d(idx_e))) / 3600;
            Energy_dch = abs(trapz(q_e, v_d(idx_e)));
        end
        if isnan(Energy_dch), Energy_dch = 0; end

        % C_eff_dch
        C_eff_dch = mean(abs(i_d)) / Q_0;

        % Assemble discharge feature vector (10 features)
        X_dch_raw = [dQ_dch, PkH_dch, PkA_dch, PkPos_dch, Energy_dch, C_eff_dch];
        X_dch_raw(isnan(X_dch_raw)) = 0;

        % Normalize: C_eff_dch 기준 NormStats 보간
        NS_dch = Results_v4.dch.NormStats;
        X_dch_norm = field_normalize(X_dch_raw, C_eff_dch, NS_dch);

        % Predict
        for i = 1:length(label_names)
            lbl = label_names{i};
            Results_Ev.(yr).dch.(lbl) = predict(Results_v4.dch.(lbl).Model, X_dch_norm);
        end
        Results_Ev.(yr).dch.C_eff = C_eff_dch;
        Results_Ev.(yr).dch.X_raw = X_dch_raw;

        fprintf('  %s [DCH] C_eff=%.3fC', yr, C_eff_dch);
        for i = 1:length(label_names)
            fprintf(' | %s=%.2f', label_names{i}, Results_Ev.(yr).dch.(label_names{i}));
        end
        fprintf('\n');
    end

    % ================================================================
    % 3C: 결과 결합 (CV RMSE 기반 가중 평균)
    % ================================================================
    for i = 1:length(label_names)
        lbl = label_names{i};
        val_chg = NaN; val_dch = NaN;
        if has_chg && isfield(Results_Ev.(yr), 'chg')
            val_chg = Results_Ev.(yr).chg.(lbl);
        end
        if has_dch && isfield(Results_Ev.(yr), 'dch')
            val_dch = Results_Ev.(yr).dch.(lbl);
        end

        if ~isnan(val_chg) && ~isnan(val_dch)
            % 가중 평균: RMSE 역수 비례 가중치
            rmse_chg = Results_v4.chg.(lbl).RMSE;
            rmse_dch = Results_v4.dch.(lbl).RMSE;
            w_chg = (1/rmse_chg) / (1/rmse_chg + 1/rmse_dch);
            w_dch = 1 - w_chg;
            Results_Ev.(yr).combined.(lbl) = w_chg * val_chg + w_dch * val_dch;
        elseif ~isnan(val_chg)
            Results_Ev.(yr).combined.(lbl) = val_chg;
        elseif ~isnan(val_dch)
            Results_Ev.(yr).combined.(lbl) = val_dch;
        else
            Results_Ev.(yr).combined.(lbl) = NaN;
        end
    end

    Results_Ev.(yr).SOH_BMS = FieldData.(yr).SOH_BMS;

    fprintf('  %s [COMBINED]', yr);
    for i = 1:length(label_names)
        fprintf(' | %s=%.2f', label_names{i}, Results_Ev.(yr).combined.(label_names{i}));
    end
    fprintf(' (BMS SOH=%d)\n', FieldData.(yr).SOH_BMS);
end

%% ========================================================================
% Section 4: 시각화 — 연도별 Trajectory
%% ========================================================================
fprintf('\n=== Section 4: Visualization ===\n');

fig = figure('Name', 'v4 Field Trajectory', 'Position', [100, 150, 1200, 500], 'Visible', 'off');

for lbl_i = 1:length(label_names)
    lname = label_names{lbl_i};
    subplot(1, 3, lbl_i);

    yr_list = fieldnames(Results_Ev);
    x_pos = 1:length(yr_list);
    val_chg = nan(1, length(yr_list));
    val_dch = nan(1, length(yr_list));
    val_comb = nan(1, length(yr_list));
    val_bms = nan(1, length(yr_list));

    for yi = 1:length(yr_list)
        yr = yr_list{yi};
        if isfield(Results_Ev.(yr), 'chg') && isfield(Results_Ev.(yr).chg, lname)
            val_chg(yi) = Results_Ev.(yr).chg.(lname);
        end
        if isfield(Results_Ev.(yr), 'dch') && isfield(Results_Ev.(yr).dch, lname)
            val_dch(yi) = Results_Ev.(yr).dch.(lname);
        end
        val_comb(yi) = Results_Ev.(yr).combined.(lname);
        if strcmp(lname, 'SOH')
            val_bms(yi) = Results_Ev.(yr).SOH_BMS;
        end
    end

    plot(x_pos, val_chg, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Charge Model');
    hold on;
    plot(x_pos, val_dch, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Discharge Model');
    plot(x_pos, val_comb, 'k-^', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Combined');
    if strcmp(lname, 'SOH')
        plot(x_pos, val_bms, 'g--d', 'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'BMS SOH');
    end

    xticks(x_pos); xticklabels(yr_list); xtickangle(30);
    ylabel(lname); title(lname, 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 8); grid on;
end
sgtitle('[v4] Field Estimation — Split Charge/Discharge Models', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig, fullfile(v4Dir, 'RF_v4_Field_Trajectory.fig'));
close(fig);

fprintf(' >> Saved: RF_v4_Field_Trajectory.fig\n');

%% ========================================================================
% Helper Functions
%% ========================================================================

function dQ_out = segment_fallback(dQ_raw, valid_mask, priority)
% Fill NaN segments using neighbor interpolation
    dQ_out = dQ_raw;
    nan_segs = find(isnan(dQ_out));
    valid_segs = find(~isnan(dQ_out));
    if isempty(valid_segs), return; end
    for s = nan_segs
        % Use nearest valid neighbor
        [~, near_idx] = min(abs(valid_segs - s));
        dQ_out(s) = dQ_out(valid_segs(near_idx));
    end
end

function [pk_height, pk_area, pk_pos] = field_extract_peak(V_uniq, Q_uniq)
% Extract dQ/dV peak features from monotonic V-Q data
    pk_height = NaN; pk_area = NaN; pk_pos = NaN;
    if numel(V_uniq) < 10, return; end

    V_u = V_uniq(:);
    Q_u = Q_uniq(:);
    if V_u(1) > V_u(end)
        V_u = flipud(V_u); Q_u = flipud(Q_u);
    end

    dV = gradient(V_u);
    dQ = gradient(Q_u);
    dV(dV == 0) = NaN;
    dQdV = dQ ./ dV;
    dQdV(isinf(dQdV) | isnan(dQdV)) = 0;
    dQdV_filt = movmean(dQdV, 21);

    pk_area = abs(trapz(V_u, dQdV_filt));

    dQdV_abs = abs(dQdV_filt);
    [pks, locs] = findpeaks(dQdV_abs, 'SortStr', 'descend', 'NPeaks', 1);
    if ~isempty(pks)
        pk_height = pks(1);
        pk_pos = V_u(locs(1));
    end
end

function X_norm = field_normalize(X_raw, C_eff, NS)
% 필드 데이터 정규화: C_eff 기준 가장 가까운 2개 C-rate 그룹 가중 보간
    crate_vals = NS.crate_vals_num;
    crate_groups = NS.crate_groups;
    norm_idx = NS.norm_idx;

    [~, sc] = sort(abs(crate_vals - C_eff));
    cr1 = crate_groups{sc(1)};
    cr2 = crate_groups{sc(2)};
    d1 = abs(C_eff - crate_vals(sc(1)));
    d2 = abs(C_eff - crate_vals(sc(2)));
    if (d1 + d2) < eps
        w1 = 1; w2 = 0;
    else
        w1 = d2 / (d1 + d2);
        w2 = d1 / (d1 + d2);
    end

    mu  = w1 * NS.(cr1).mu  + w2 * NS.(cr2).mu;
    sig = w1 * NS.(cr1).sigma + w2 * NS.(cr2).sigma;

    X_norm = X_raw;
    X_norm(norm_idx) = (X_raw(norm_idx) - mu) ./ (sig + eps);
    % C_eff: raw 유지 (마지막 피처)
end
