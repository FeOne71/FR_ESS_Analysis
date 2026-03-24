% RPT_Field_Estimation_RF.m
% 1. 학습 완료된 100% Data 기반 Random Forest 실생성 모델 로드
% 2. 필드 데이터(Field Data, 2021~2025) 로드 및 20개 물리 피처 추출
% 3. 학습된 RF 모델을 사용하여 필드 데이터의 SOH, LLI, LAM, SOP 추정 및 시각화

clear; clc; close all;
warning on;

%% ========================================================================
% Section 1: 모델 로드 (Random Forest)
% ========================================================================
fprintf('=== Section 1: Loading Trained RF Model ===\n');
baseDir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
rfDir   = fullfile(baseDir, 'ML_RandomForest');
modelFile = fullfile(rfDir, 'Result_RandomForest.mat');
path_master_ruler = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', ...
    'ExperimentalData', 'FeatureEngineering', 'MasterRulers_v3.mat');

if ~exist(modelFile, 'file'), error('RF Model not found.'); end
load(modelFile, 'Results_RF');
load(path_master_ruler, 'MasterRulers');

label_names = Results_RF.label_names;

%% ========================================================================
% Section 2: 필드 데이터 로드 (KIMJ ESS 2021~2025)
% Segment logic: identical to Visualize_VIT_byYear.m / FieldQmax_dQdV.m
% ========================================================================
fprintf('\n=== Section 2: Loading Field Data (2021~2025) ===\n');
rawDir   = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'Rack_raw2mat');
dataFiles = {
    fullfile(rawDir, 'Old', '2021', '202106', 'Raw_20210603.mat'), 'Y2021', 'old',  datetime(2021,6,3);
    fullfile(rawDir, 'New', '2023', '202310', 'Raw_20231016.mat'), 'Y2023', 'new',  datetime(2023,10,16);
    % fullfile(rawDir, 'New', '2024', '202409', 'Raw_20240909.mat'), 'Y2024_Auto', 'new',  datetime(2024,9,9);
    fullfile(rawDir, 'New', '2024', '202409', 'Raw_20240909.mat'), 'Y2024', 'new',  datetime(2024,9,9);
    fullfile(rawDir, 'New', '2025', '202507', 'Raw_20250711.mat'), 'Y2025', 'new',  datetime(2024,9,9);
    % Y2025 excluded: voltage range too narrow
};

% Per-year minimum segment durations (seconds) — same as reference script
min_charge_secs    = [600, 300, 300  300];   % Y2021, Y2023, Y2024_Auto, Y2024_Manual, Y2025
min_discharge_secs = [300, 150, 300, 300];

dt = 1;

% Global Master Ruler Definition
fns = fieldnames(MasterRulers);
VR_chg = MasterRulers.(fns{1}).V_bounds_chg;
VR_dch = MasterRulers.(fns{1}).V_bounds_dch;

Np = 2; 
C_cell_Ah = 64; 
thr_A = C_cell_Ah * 0.05;

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

    % BMS SOH: 0은 센서 이상값 → NaN 처리 후 하루 대표값 = 중앙값
    valid_soh = raw_soh_bms;
    valid_soh(valid_soh <= 0) = NaN;
    soh_bms_day = round(median(valid_soh, 'omitnan'));

    tsec   = seconds(t - t(1));
    I_cell = I_rack / Np;

    % Store full profile
    FieldData.(year_label).Full.V    = V_avg;
    FieldData.(year_label).Full.I    = I_cell;
    FieldData.(year_label).Full.Time = tsec;
    FieldData.(year_label).SOH_BMS  = soh_bms_day;

    % Segment detection (same logic as reference)
    isChg   = I_cell >  thr_A;
    isDch   = I_cell < -thr_A;
    chgSegs = local_find_segments(isChg);
    dchSegs = local_find_segments(isDch);

    % Filter by minimum duration
    if ~isempty(chgSegs)
        dur_c = chgSegs(:,2) - chgSegs(:,1) + 1;
        chgSegs = chgSegs(dur_c >= ceil(min_chg_sec/dt), :);
    end
    if ~isempty(dchSegs)
        dur_d = dchSegs(:,2) - dchSegs(:,1) + 1;
        dchSegs = dchSegs(dur_d >= ceil(min_dch_sec/dt), :);
    end

    % -------------------------------------------------------
    % Y2024_Manual: Use specific segments Chg02 & Dchg03 by time
    % -------------------------------------------------------
    if strcmp(year_label, 'Y2024_Manual')
        % Charge: target 12:55 ~ 13:13
        target_chg_start = datetime(2024, 9, 9, 12, 55, 0);
        target_chg_end   = datetime(2024, 9, 9, 13, 13, 0);
        chg_seg_idx = [];
        for si = 1:size(chgSegs, 1)
            if abs(t(chgSegs(si,1)) - target_chg_start) < minutes(5) && ...
               abs(t(chgSegs(si,2)) - target_chg_end)   < minutes(5)
                chg_seg_idx = si; break;
            end
        end
        if isempty(chg_seg_idx), chgSegs = [];
        else, chgSegs = chgSegs(chg_seg_idx, :); end

        % Discharge: starting ~14:14, ending at 14:25:49
        target_dch_start = datetime(2024, 9, 9, 14, 14, 0);
        target_rest_end  = datetime(2024, 9, 9, 14, 25, 49);
        dch_seg_idx = [];
        for si = 1:size(dchSegs, 1)
            if abs(t(dchSegs(si,1)) - target_dch_start) < minutes(1)
                dch_seg_idx = si; break;
            end
        end
        if ~isempty(dch_seg_idx)
            dch_end = dchSegs(dch_seg_idx, 2);
            if t(dch_end) > target_rest_end
                idx_limit = find(t <= target_rest_end, 1, 'last');
                if ~isempty(idx_limit) && idx_limit >= dchSegs(dch_seg_idx, 1)
                    dchSegs(dch_seg_idx, 2) = idx_limit;
                end
            end
            dchSegs = dchSegs(dch_seg_idx, :);
        else
            dchSegs = [];
        end
    end

    % Pick the LONGEST surviving charge segment

    chg_s = NaN; chg_e = NaN;
    if ~isempty(chgSegs)
        [~, best] = max(chgSegs(:,2) - chgSegs(:,1));
        chg_s = chgSegs(best, 1);  chg_e = chgSegs(best, 2);
        FieldData.(year_label).Chg.Time = tsec(chg_s:chg_e);
        FieldData.(year_label).Chg.V = V_avg(chg_s:chg_e);
        FieldData.(year_label).Chg.I = I_cell(chg_s:chg_e);
        FieldData.(year_label).Chg.T = T_avg(chg_s:chg_e);
        FieldData.(year_label).Chg.Q = cumtrapz(abs(I_cell(chg_s:chg_e))) / 3600;
    end

    dch_s = NaN; dch_e = NaN;
    if ~isempty(dchSegs)
        [~, best] = max(dchSegs(:,2) - dchSegs(:,1));
        dch_s = dchSegs(best, 1);  dch_e = dchSegs(best, 2);
        FieldData.(year_label).Dch.Time = tsec(dch_s:dch_e);
        FieldData.(year_label).Dch.V = V_avg(dch_s:dch_e);
        FieldData.(year_label).Dch.I = I_cell(dch_s:dch_e);
        FieldData.(year_label).Dch.T = T_avg(dch_s:dch_e);
        FieldData.(year_label).Dch.Q = cumtrapz(abs(I_cell(dch_s:dch_e))) / 3600;
    end

    % Segment Visualization
    fig_seg = figure('Name', sprintf('[%s] Field Selection', year_label), ...
        'Position', [100, 100, 1000, 600], 'Visible', 'on');
    sgtitle(sprintf('KIMJ ESS %s - Daily Profile & Selected Segments', year_label), ...
        'FontSize', 14, 'FontWeight', 'bold');
    t_hr = tsec / 3600;
    ax1 = subplot(2,1,1); hold on; grid on; box on;
    plot(t_hr, V_avg, 'k-', 'LineWidth', 1);
    ylims1 = [min(V_avg)-0.05, max(V_avg)+0.05];
    if ~isnan(chg_s), patch(t_hr([chg_s chg_e chg_e chg_s]), ylims1([1 1 2 2]), 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); end
    if ~isnan(dch_s), patch(t_hr([dch_s dch_e dch_e dch_s]), ylims1([1 1 2 2]), 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); end
    ylim(ylims1); ylabel('Voltage (V)'); title('Voltage Profile');
    ax2 = subplot(2,1,2); hold on; grid on; box on;
    plot(t_hr, I_cell, 'k-', 'LineWidth', 1);
    ylims2 = [min(I_cell)-1, max(I_cell)+1];
    if ~isnan(chg_s), patch(t_hr([chg_s chg_e chg_e chg_s]), ylims2([1 1 2 2]), 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); end
    if ~isnan(dch_s), patch(t_hr([dch_s dch_e dch_e dch_s]), ylims2([1 1 2 2]), 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); end
    ylim(ylims2); ylabel('Current (A)'); xlabel('Time (Hours)'); title('Current Profile');
    legend('Raw Data', 'Extracted Charge', 'Extracted Discharge', 'Location', 'best');
    linkaxes([ax1, ax2], 'x'); xlim([0, 24]);
    saveas(fig_seg, fullfile(rfDir, sprintf('Segment_Selection_%s.fig', year_label)));
    close(fig_seg);
end

%% ========================================================================
% Section 3: 피처 추출 및 상태 추정 (Approach A: Raw Features + Merged RF)
% ========================================================================
fprintf('\n=== Section 3: Feature Extraction & RF Estimation ===\n');
years = fieldnames(FieldData);
Results_Ev = struct();

% Segment Fallback priority (중앙-안정 세그먼트 우선)
fallback_priority = [3, 4, 2, 5, 1];


% Load MasterRulers_v3 for voltage segments
mr_path = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering', 'MasterRulers_v3.mat');
mr_data = load(mr_path, 'MasterRulers');
mr_ch = fieldnames(mr_data.MasterRulers); mr_ch = mr_ch{1};
VR_chg = mr_data.MasterRulers.(mr_ch).V_bounds_chg;  % 6 voltage boundaries (ascending)
VR_dch = mr_data.MasterRulers.(mr_ch).V_bounds_dch;  % 6 voltage boundaries (descending)
ma_window = 60; Q_0 = 64;

for k = 1:length(years)
    yr = years{k};
    if ~isfield(FieldData.(yr), 'Chg') || ~isfield(FieldData.(yr), 'Dch')
        fprintf('  %s: Missing Chg/Dch Segment. Skipped.\n', yr); continue;
    end

    % -- 1. Load raw data --
    t_c = FieldData.(yr).Chg.Time(:); v_c = FieldData.(yr).Chg.V(:); i_c = FieldData.(yr).Chg.I(:);
    t_d = FieldData.(yr).Dch.Time(:); v_d = FieldData.(yr).Dch.V(:); i_d = FieldData.(yr).Dch.I(:);
    i_c_sm = movmean(abs(i_c), ma_window);
    i_d_sm = movmean(abs(i_d), ma_window);
    q_c = cumtrapz(t_c, i_c_sm) / 3600;
    q_d = cumtrapz(t_d, i_d_sm) / 3600;

    % -- 2. Smooth V and Q --
    v_c_sm = movmean(v_c, ma_window); q_c_sm = movmean(q_c, ma_window);
    v_d_sm = movmean(v_d, ma_window); q_d_sm = movmean(q_d, ma_window);

    % -- 3. Remove duplicate V (preserve time order) --
    [v_c_uniq, uc_idx] = unique(v_c_sm, 'stable'); q_c_uniq = q_c_sm(uc_idx);
    [v_d_uniq, ud_idx] = unique(v_d_sm, 'stable'); q_d_uniq = q_d_sm(ud_idx);

    % -- 4. Remove non-monotonic points --
    mono_c = true(size(v_c_uniq));
    for ii = 2:length(v_c_uniq)
        if v_c_uniq(ii) <= v_c_uniq(ii-1), mono_c(ii) = false; end
    end
    v_c_uniq = v_c_uniq(mono_c); q_c_uniq = q_c_uniq(mono_c);

    mono_d = true(size(v_d_uniq));
    for ii = 2:length(v_d_uniq)
        if v_d_uniq(ii) >= v_d_uniq(ii-1), mono_d(ii) = false; end
    end
    v_d_uniq = v_d_uniq(mono_d); q_d_uniq = q_d_uniq(mono_d);

    % ================================================================
    % SECTION 3A: Segment Validity Check + Fallback 
    % ================================================================
    % Validity check: use RAW voltage range (MasterRuler was built from raw field V)
    % Using smoothed v_c_uniq would truncate edges and cause false invalids
    V_min_field_chg = min(v_c); V_max_field_chg = max(v_c);
    V_min_field_dch = min(v_d); V_max_field_dch = max(v_d);


    % Valid mask: segment [VR(k), VR(k+1)] must be fully within field V range
    num_segs = length(VR_chg) - 1;  % = 5
    valid_chg = false(1, num_segs);
    valid_dch = false(1, num_segs);
    for s = 1:num_segs
        v_lo_c = VR_chg(s);   v_hi_c = VR_chg(s+1);
        valid_chg(s) = (v_lo_c >= V_min_field_chg - 0.02) && (v_hi_c <= V_max_field_chg + 0.02);
        v_hi_d = VR_dch(s);   v_lo_d = VR_dch(s+1);  % descending order
        valid_dch(s) = (v_hi_d <= V_max_field_dch + 0.02) && (v_lo_d >= V_min_field_dch - 0.02);
    end

    % Compute raw dQ for valid segments using 0.001V grid & moving average
    dQ_chg_raw = nan(1, num_segs);
    dQ_dch_raw = nan(1, num_segs);
    
    if length(v_c_uniq) > 1
        try
            v_grid_c = (VR_chg(1) - 0.05) : 0.001 : (VR_chg(end) + 0.05);
            Q_interp_c = interp1(v_c_uniq, q_c_uniq, v_grid_c, 'linear', 'extrap');
            Q_sm_c = movmean(Q_interp_c, 30);
            
            QR_chg = nan(1, length(VR_chg));
            for b = 1:length(VR_chg)
                [~, min_idx] = min(abs(v_grid_c - VR_chg(b)));
                QR_chg(b) = Q_sm_c(min_idx);
            end
            dQ_chg_raw = abs(diff(QR_chg));
            dQ_chg_raw(~valid_chg) = NaN;  % invalidate out-of-range segs
        catch ME
            fprintf('  %s Chg interpolation failed: %s\n', yr, ME.message);
        end
    end
    
    if length(v_d_uniq) > 1
        try
            v_grid_d = (VR_dch(1) + 0.05) : -0.001 : (VR_dch(end) - 0.05);
            Q_interp_d = interp1(v_d_uniq, q_d_uniq, v_grid_d, 'linear', 'extrap');
            Q_sm_d = movmean(Q_interp_d, 30);
            
            QR_dch = nan(1, length(VR_dch));
            for b = 1:length(VR_dch)
                [~, min_idx] = min(abs(v_grid_d - VR_dch(b)));
                QR_dch(b) = Q_sm_d(min_idx);
            end
            dQ_dch_raw = abs(diff(QR_dch));
            dQ_dch_raw(~valid_dch) = NaN;
        catch ME
            fprintf('  %s Dch interpolation failed: %s\n', yr, ME.message);
        end
    end

    % Fallback: fill NaN segs using linear interpolation of neighbors
    dQ_chg = segment_fallback(dQ_chg_raw, valid_chg, fallback_priority);
    dQ_dch = segment_fallback(dQ_dch_raw, valid_dch, fallback_priority);

    % Report
    fprintf('  %s | Valid Chg Segs: [%s]  Dch Segs: [%s]\n', yr, ...
        strjoin(arrayfun(@(x) num2str(x), valid_chg, 'UniformOutput', false), ' '), ...
        strjoin(arrayfun(@(x) num2str(x), valid_dch, 'UniformOutput', false), ' '));

    % -- 6. dQ/dV Peak Extraction --
    PkH_chg = NaN; PkA_chg = NaN; PkPos_chg = NaN;
    PkH_dch = NaN; PkA_dch = NaN; PkPos_dch = NaN;
    [PkH_chg, PkA_chg, PkPos_chg, V_chg_curve, dQdV_chg_curve] = field_extract_peak(v_c_uniq, q_c_uniq);
    [PkH_dch, PkA_dch, PkPos_dch, V_dch_curve, dQdV_dch_curve] = field_extract_peak(v_d_uniq, q_d_uniq);

    % -- 7. Energy (Wh) --
    Energy_dch = NaN;
    idx_d = FieldData.(yr).Dch.V >= min(VR_dch) & FieldData.(yr).Dch.V <= max(VR_dch);
    if sum(idx_d) > 5
        Energy_dch = abs(trapz(FieldData.(yr).Dch.Q(idx_d), FieldData.(yr).Dch.V(idx_d)));
    end

    % -- 8. C_eff --
    C_eff_chg = mean(abs(FieldData.(yr).Chg.I)) / Q_0;
    C_eff_dch = mean(abs(FieldData.(yr).Dch.I)) / Q_0;
    T_chg_avg = mean(FieldData.(yr).Chg.T, 'omitnan');
    T_dch_avg = mean(FieldData.(yr).Dch.T, 'omitnan');

    % ================================================================
    % SECTION 3B: Split Normalization + Predict
    %   Charge features → C_eff_chg 기준 NormStats 그룹 보간 정규화
    %   Discharge features → C_eff_dch 기준 NormStats 그룹 보간 정규화
    %   C_eff features → raw 유지 (정규화 없음)
    %
    %   [v3] NormStats 구조 변경 반영:
    %     기존: NormStats.(cr).mu(idx) / .sigma(idx)
    %     신규: NormStats.(cr).mu_chg / .sigma_chg  (충전 피처용)
    %           NormStats.(cr).mu_dch / .sigma_dch  (방전 피처용)
    % ================================================================
    NormStats      = Results_RF.NormStats;
    crate_groups   = NormStats.crate_groups;
    crate_vals_num = NormStats.crate_vals_num;

    % Feature indices (19 features, T_avg 제외 기준)
    chg_idx  = [1:5, 11, 13, 15];        % dQ_chg×5, PkH_chg, PkA_chg, PkPos_chg (8개)
    dch_idx  = [6:10, 12, 14, 16, 17];   % dQ_dch×5, PkH_dch, PkA_dch, PkPos_dch, Energy_dch (9개)
    ceff_idx = [18, 19];                  % C_eff_chg, C_eff_dch (정규화 없음)

    X_raw = [dQ_chg, dQ_dch, PkH_chg, PkH_dch, PkA_chg, PkA_dch, PkPos_chg, PkPos_dch, Energy_dch, C_eff_chg, C_eff_dch];

    % NaN check
    if any(isnan(X_raw))
        n_nan = sum(isnan(X_raw));
        fprintf('  %s: %d NaN features remain after fallback. Replacing with 0.\n', yr, n_nan);
        X_raw(isnan(X_raw)) = 0;
    end

    % --- 충전 피처: C_eff_chg 기준 가장 가까운 2개 그룹 가중 보간 ---
    [~, sc] = sort(abs(crate_vals_num - C_eff_chg));
    cr1_c = crate_groups{sc(1)}; cr2_c = crate_groups{sc(2)};
    d1c = abs(C_eff_chg - crate_vals_num(sc(1)));
    d2c = abs(C_eff_chg - crate_vals_num(sc(2)));
    if (d1c+d2c) < eps, w1c=1; w2c=0; else, w1c=d2c/(d1c+d2c); w2c=d1c/(d1c+d2c); end
    mu_chg  = w1c * NormStats.(cr1_c).mu_chg    + w2c * NormStats.(cr2_c).mu_chg;
    sig_chg = w1c * NormStats.(cr1_c).sigma_chg + w2c * NormStats.(cr2_c).sigma_chg;

    % --- 방전 피처: C_eff_dch 기준 가장 가까운 2개 그룹 가중 보간 ---
    [~, sd] = sort(abs(crate_vals_num - C_eff_dch));
    cr1_d = crate_groups{sd(1)}; cr2_d = crate_groups{sd(2)};
    d1d = abs(C_eff_dch - crate_vals_num(sd(1)));
    d2d = abs(C_eff_dch - crate_vals_num(sd(2)));
    if (d1d+d2d) < eps, w1d=1; w2d=0; else, w1d=d2d/(d1d+d2d); w2d=d1d/(d1d+d2d); end
    mu_dch  = w1d * NormStats.(cr1_d).mu_dch    + w2d * NormStats.(cr2_d).mu_dch;
    sig_dch = w1d * NormStats.(cr1_d).sigma_dch + w2d * NormStats.(cr2_d).sigma_dch;

    % --- 분리 정규화 적용 ---
    X_norm = X_raw;
    X_norm(chg_idx)  = (X_raw(chg_idx)  - mu_chg)  ./ (sig_chg  + eps);
    X_norm(dch_idx)  = (X_raw(dch_idx)  - mu_dch)  ./ (sig_dch  + eps);
    % ceff_idx: raw 유지 (C-rate 단위 그대로 사용)

    fprintf('  %s | C_eff_chg=%.3fC→%s(w=%.0f%%)  C_eff_dch=%.3fC→%s(w=%.0f%%)\n', ...
        yr, C_eff_chg, cr1_c, w1c*100, C_eff_dch, cr1_d, w1d*100);

    % Predict using Merged RF model (trained on normalized features)
    for i = 1:length(label_names)
        lbl = label_names{i};
        Results_Ev.(yr).Merged.(lbl) = predict(Results_RF.Merged.(lbl).Model, X_norm);
    end

    % Store
    Results_Ev.(yr).X_field  = X_raw;
    Results_Ev.(yr).X_norm   = X_norm;
    Results_Ev.(yr).Segment.valid_chg = valid_chg;
    Results_Ev.(yr).Segment.valid_dch = valid_dch;
    Results_Ev.(yr).Ceff.chg = C_eff_chg;
    Results_Ev.(yr).Ceff.dch = C_eff_dch;
    Results_Ev.(yr).Curves.V_chg = V_chg_curve;
    Results_Ev.(yr).Curves.dQdV_chg = dQdV_chg_curve;
    Results_Ev.(yr).Curves.V_dch = V_dch_curve;
    Results_Ev.(yr).Curves.dQdV_dch = dQdV_dch_curve;

    fprintf('  %s [Merged-Normalized]', yr);
    for i = 1:length(label_names)
        lbl = label_names{i};
        fprintf(' | %s: %.2f', lbl, Results_Ev.(yr).Merged.(lbl));
    end
    fprintf('\n');



    % Feature bar chart (unchanged)
    fig_fea = figure('Name', sprintf('[%s] Extracted Features', yr), 'Position', [150, 150, 1000, 400], 'Visible', 'on');
    bar(X_raw, 'FaceColor', [0.2 0.6 0.5]); grid on;
    title(sprintf('Extracted Physics Features - KIMJ %s', yr), 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Feature Value');
    xticks(1:length(Results_RF.Merged.feature_names)); xticklabels(Results_RF.Merged.feature_names); xtickangle(45);
    saveas(fig_fea, fullfile(rfDir, sprintf('Features_%s.fig', yr)));
    close(fig_fea);
end



%% ========================================================================
% Section 4: 결과 시각화 — Merged vs Method B Comparison
% ========================================================================
yr_strs = fieldnames(Results_Ev);
n_yrs   = length(yr_strs);
subplot_colors = lines(length(label_names));

% --- 4A. Label Estimation Trajectory (Line Plot) ---
fig_rf_field = figure('Name', 'RF Field Estimation', 'Position', [100, 100, 1400, 450]);
sgtitle('RF Field Estimation Trajectory (KIMJ ESS)', 'FontSize', 14, 'FontWeight', 'bold');

% Extract actual year numbers, add offset for Manual to prevent overlap
yr_nums = zeros(n_yrs, 1);
for k = 1:n_yrs
    tmp = regexp(yr_strs{k}, '\d{4}', 'match');
    base_yr = str2double(tmp{1});
    if contains(yr_strs{k}, 'Manual')
        yr_nums(k) = base_yr + 0.4; % offset for visual separation
    else
        yr_nums(k) = base_yr;
    end
end

% X-axis labels (공통)
x_labels = cell(n_yrs, 1);
for k = 1:n_yrs
    yr_name = strrep(yr_strs{k}, 'Y', '');
    yr_name = strrep(yr_name, '_', ' ');
    x_labels{k} = sprintf('%s', yr_name);
end

line_colors = [0.2 0.4 0.8; 0.9 0.3 0.1; 0.1 0.7 0.3];
y_labels    = {'SOH (%)', 'LLI (%)', 'LAM (%)'};

% --- Subplot 1: SOH (BMS vs RF Model) ---
subplot(1, 3, 1); hold on; box on; grid on;

% BMS SOH 추출
soh_bms_vals = nan(n_yrs, 1);
for k = 1:n_yrs
    if isfield(FieldData.(yr_strs{k}), 'SOH_BMS')
        soh_bms_vals(k) = FieldData.(yr_strs{k}).SOH_BMS;
    end
end

% RF 추정 SOH
soh_rf_vals = zeros(n_yrs, 1);
for k = 1:n_yrs
    soh_rf_vals(k) = Results_Ev.(yr_strs{k}).Merged.SOH;
end

% BMS SOH (검정 점선)
plot(yr_nums, soh_bms_vals, 'k--s', 'LineWidth', 2.5, 'MarkerSize', 10, ...
    'MarkerFaceColor', 'k', 'DisplayName', 'BMS SOH');
% RF 추정 SOH (파란 실선)
plot(yr_nums, soh_rf_vals, '-o', 'Color', line_colors(1,:), 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', line_colors(1,:), 'DisplayName', 'RF Estimated SOH');

% Data labels
for k = 1:n_yrs
    if ~isnan(soh_bms_vals(k))
        text(yr_nums(k), soh_bms_vals(k)+0.5, sprintf('%.1f%%', soh_bms_vals(k)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
    end
    text(yr_nums(k), soh_rf_vals(k)-0.8, sprintf('%.1f%%', soh_rf_vals(k)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', ...
        'Color', line_colors(1,:)*0.7);
end

set(gca, 'XTick', yr_nums, 'XTickLabel', x_labels, 'FontSize', 9);
xtickangle(45);
xlabel('Year'); ylabel('SOH (%)');
title('SOH', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
xlim([min(yr_nums)-0.5, max(yr_nums)+0.5]);

% --- Subplot 2: LLI (RF Model only, no field data) ---
subplot(1, 3, 2); hold on; box on; grid on;

lli_rf_vals = zeros(n_yrs, 1);
for k = 1:n_yrs
    lli_rf_vals(k) = Results_Ev.(yr_strs{k}).Merged.LLI;
end

plot(yr_nums, lli_rf_vals, '-o', 'Color', line_colors(2,:), 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', line_colors(2,:), 'DisplayName', 'RF Estimated LLI');

for k = 1:n_yrs
    text(yr_nums(k), lli_rf_vals(k), sprintf('  %.2f', lli_rf_vals(k)), ...
        'VerticalAlignment', 'bottom', 'FontSize', 9, 'FontWeight', 'bold', ...
        'Color', line_colors(2,:)*0.7);
end

set(gca, 'XTick', yr_nums, 'XTickLabel', x_labels, 'FontSize', 9);
xtickangle(45);
xlabel('Year'); ylabel('LLI (%)');
title('LLI', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
xlim([min(yr_nums)-0.5, max(yr_nums)+0.5]);

% --- Subplot 3: LAM (RF Model only, no field data) ---
subplot(1, 3, 3); hold on; box on; grid on;

lam_rf_vals = zeros(n_yrs, 1);
for k = 1:n_yrs
    lam_rf_vals(k) = Results_Ev.(yr_strs{k}).Merged.LAM;
end

plot(yr_nums, lam_rf_vals, '-o', 'Color', line_colors(3,:), 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', line_colors(3,:), 'DisplayName', 'RF Estimated LAM');

for k = 1:n_yrs
    text(yr_nums(k), lam_rf_vals(k), sprintf('  %.2f', lam_rf_vals(k)), ...
        'VerticalAlignment', 'bottom', 'FontSize', 9, 'FontWeight', 'bold', ...
        'Color', line_colors(3,:)*0.7);
end

set(gca, 'XTick', yr_nums, 'XTickLabel', x_labels, 'FontSize', 9);
xtickangle(45);
xlabel('Year'); ylabel('LAM (%)');
title('LAM', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
xlim([min(yr_nums)-0.5, max(yr_nums)+0.5]);

saveas(fig_rf_field, fullfile(rfDir, 'RF_Field_Trajectory.fig'));


% --- 4B. Segment Validity Heatmap (Charge & Discharge) ---
fig_seg = figure('Name', 'Segment Validity', 'Position', [150, 150, 900, 350]);
sgtitle('Field Segment Validity (1=Valid, 0=Fallback-filled)', ...
    'FontSize', 13, 'FontWeight', 'bold');

num_segs = 5;
valid_chg_mat = zeros(n_yrs, num_segs);
valid_dch_mat = zeros(n_yrs, num_segs);
for k = 1:n_yrs
    if isfield(Results_Ev.(yr_strs{k}), 'Segment')
        valid_chg_mat(k,:) = double(Results_Ev.(yr_strs{k}).Segment.valid_chg);
        valid_dch_mat(k,:) = double(Results_Ev.(yr_strs{k}).Segment.valid_dch);
    end
end

subplot(1,2,1);
imagesc(valid_chg_mat');
colormap([0.95 0.3 0.2; 0.2 0.75 0.3]);  % red=invalid, green=valid
clim([0 1]);
set(gca, 'XTick', 1:n_yrs, 'XTickLabel', yr_strs, ...
         'YTick', 1:num_segs, 'YTickLabel', arrayfun(@(x) sprintf('Seg%d',x), 1:num_segs, 'UniformOutput', false));
xlabel('Year'); title('Charge Segments');
for k=1:n_yrs; for s=1:num_segs
    text(k, s, num2str(valid_chg_mat(k,s)), 'HorizontalAlignment','center', ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', 'w');
end; end

subplot(1,2,2);
imagesc(valid_dch_mat');
colormap([0.95 0.3 0.2; 0.2 0.75 0.3]);
clim([0 1]);
set(gca, 'XTick', 1:n_yrs, 'XTickLabel', yr_strs, ...
         'YTick', 1:num_segs, 'YTickLabel', arrayfun(@(x) sprintf('Seg%d',x), 1:num_segs, 'UniformOutput', false));
xlabel('Year'); title('Discharge Segments');
for k=1:n_yrs; for s=1:num_segs
    text(k, s, num2str(valid_dch_mat(k,s)), 'HorizontalAlignment','center', ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', 'w');
end; end

saveas(fig_seg, fullfile(rfDir, 'RF_Segment_Validity.fig'));


%% ========================================================================
% Section 5: dQ/dV Peak Features Visualization
% ========================================================================
fig_peaks = figure('Name', 'dQ/dV Peak Features', 'Position', [150, 150, 1000, 600]);
sgtitle('Field Data dQ/dV Peak Progression (2021-2025)', 'FontSize', 15, 'FontWeight', 'bold');

% X_merged indeces: PkH_chg(11), PkH_dch(12), PkA_chg(13), PkA_dch(14)
peak_labels = {'Peak Height (Charge)', 'Peak Height (Discharge)', 'Peak Area (Charge)', 'Peak Area (Discharge)'};
peak_indices = [11, 12, 13, 14];
colors_pk = lines(4);

for i = 1:4
    subplot(2, 2, i); hold on; box on; grid on;
    y_vals = zeros(1, n_yrs);
    for k = 1:n_yrs
        X_f = Results_Ev.(yr_strs{k}).X_field;
        y_vals(k) = X_f(peak_indices(i));
    end
    
    plot(1:n_yrs, y_vals, '-s', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', colors_pk(i,:), 'Color', colors_pk(i,:));
    for k = 1:n_yrs
        if y_vals(k) ~= 0 % Exclude completely NaN/missing segments
            text(k, y_vals(k), sprintf('%.2f', y_vals(k)), 'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', colors_pk(i,:)*0.8);
        end
    end
    set(gca, 'XTick', 1:n_yrs, 'XTickLabel', yr_strs);
    title(peak_labels{i}); 
    if contains(peak_labels{i}, 'Height'), ylabel('dQ/dV (Ah/V)'); else, ylabel('Capacity (Ah)'); end
    xlim([0.5, n_yrs+0.5]);
    
    % Adjust Y-limits slightly if all values are very close
    if max(y_vals) - min(y_vals) < 0.1 && max(y_vals) > 0
        ylim([mean(y_vals)*0.95, mean(y_vals)*1.05]);
    end
end
saveas(fig_peaks, fullfile(rfDir, 'Peak_Features_Trajectory.fig'));

%% ========================================================================
% Section 6: Overlaid dQ/dV Curve Visualization
% ========================================================================
fig_curves = figure('Name', 'Overlaid dQ/dV Curves', 'Position', [150, 150, 1000, 500]);
sgtitle('Field Data dQ/dV Curves Overlay (2021-2025)', 'FontSize', 15, 'FontWeight', 'bold');

colors_yr = lines(n_yrs);

subplot(1, 2, 1); hold on; box on; grid on;
for k = 1:n_yrs
    yr = yr_strs{k};
    V_c = Results_Ev.(yr).Curves.V_chg;
    dQ_c = Results_Ev.(yr).Curves.dQdV_chg;
    if ~isempty(V_c)
        plot(V_c, dQ_c, '-', 'LineWidth', 2, 'Color', colors_yr(k,:), 'DisplayName', yr);
    end
    % Peak marker
    pk_pos = Results_Ev.(yr).X_field(15); % Peak_Pos_chg index
    pk_h   = Results_Ev.(yr).X_field(11); % Peak_H_chg index
    if ~isnan(pk_pos) && ~isnan(pk_h)
        scatter(pk_pos, pk_h, 80, colors_yr(k,:), 'o', 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
        text(pk_pos, pk_h, sprintf('  %.1f@%.3fV', pk_h, pk_pos), 'FontSize', 8, 'Color', colors_yr(k,:)*0.8, 'HandleVisibility','off');
    end
end
title('Charge dQ/dV'); xlabel('Voltage (V)'); ylabel('dQ/dV (Ah/V)');
legend('Location', 'best');

subplot(1, 2, 2); hold on; box on; grid on;
for k = 1:n_yrs
    yr = yr_strs{k};
    V_d = Results_Ev.(yr).Curves.V_dch;
    dQ_d = Results_Ev.(yr).Curves.dQdV_dch;
    if ~isempty(V_d)
        plot(V_d, dQ_d, '-', 'LineWidth', 2, 'Color', colors_yr(k,:), 'DisplayName', yr);
    end
    % Peak marker
    pk_pos = Results_Ev.(yr).X_field(16); % Peak_Pos_dch index
    pk_h   = Results_Ev.(yr).X_field(12); % Peak_H_dch index
    if ~isnan(pk_pos) && ~isnan(pk_h)
        scatter(pk_pos, pk_h, 80, colors_yr(k,:), 'o', 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
        text(pk_pos, pk_h, sprintf('  %.1f@%.3fV', pk_h, pk_pos), 'FontSize', 8, 'Color', colors_yr(k,:)*0.8, 'HandleVisibility','off');
    end
end
title('Discharge dQ/dV'); xlabel('Voltage (V)'); ylabel('dQ/dV (Ah/V)');
legend('Location', 'best');

saveas(fig_curves, fullfile(rfDir, 'Overlaid_dQdV_Curves.fig'));

fprintf('=== Field Estimation Strategy Executed (RF) ===\n');

%% Helper Function
function [pk_height, pk_area, pk_pos, V_out, dQdV_out] = field_extract_peak(V_uniq, Q_uniq)
    % findpeaks on full dQ/dV curve — matches Lab RPT_FeatureExtractor logic
    pk_height = NaN; pk_area = NaN; pk_pos = NaN; V_out = []; dQdV_out = [];
    if length(V_uniq) <= 5, return; end

    % 1. Full 0.001V grid spanning the ENTIRE extracted segment
    V_grid = (ceil(min(V_uniq)*1000)/1000 : 0.001 : floor(max(V_uniq)*1000)/1000)';

    % 2. Interpolate Q onto 0.001V grid
    Q_grid = interp1(V_uniq, Q_uniq, V_grid, 'linear', 'extrap');

    % 3. dQ/dV with fixed 0.001V denominator
    dV = 0.001;
    dQdV = abs(diff(Q_grid)) ./ dV;
    V_mid = V_grid(1:end-1) + dV/2;

    % 4. Moving average post-processing (window=10)
    dQdV = movmean(dQdV, 10);
    dQdV = max(dQdV, 0);

    % 5. Trim boundary artifacts
    edge_trim = 10;
    if length(dQdV) > 2*edge_trim + 1
        dQdV = dQdV(edge_trim+1:end-edge_trim);
        V_mid = V_mid(edge_trim+1:end-edge_trim);
    end

    % 6. Full curve output for visualization
    V_out = V_mid;
    dQdV_out = dQdV;

    % 7. Peak_A: trapz over full curve
    pk_area = abs(trapz(V_mid, dQdV));

    % 8. Peak_H & Peak_Pos: dominant peak via findpeaks (same as Lab)
    [pks, locs] = findpeaks(dQdV, 'SortStr', 'descend', 'NPeaks', 1);
    if ~isempty(pks)
        pk_height = pks(1);
        pk_pos = V_mid(locs(1));
    end
end

function dQ_out = segment_fallback(dQ_raw, valid_mask, ~)
    % segment_fallback: Fill invalid dQ segments using linear interpolation
    %   of nearest valid neighbours.
    %   dQ_raw    : 1×N raw dQ values (NaN where invalid)
    %   valid_mask: 1×N logical, true = valid
    %   fallback_priority: (unused, kept for API compatibility)
    %
    %   Returns dQ_out with invalid segs filled.  If NO valid segs exist,
    %   dQ_out is all-NaN (caller must handle).

    dQ_out = dQ_raw;  % start with raw (valid values preserved)
    N = length(dQ_raw);
    seg_idx = 1:N;

    valid_idx = seg_idx(valid_mask);
    if isempty(valid_idx)
        return;  % all invalid — leave as NaN
    end

    for s = seg_idx
        if valid_mask(s), continue; end  % already valid, skip

        % Find nearest valid segments to the left and right
        left_valid  = valid_idx(valid_idx < s);
        right_valid = valid_idx(valid_idx > s);

        if ~isempty(left_valid) && ~isempty(right_valid)
            l = left_valid(end);   r = right_valid(1);
            % Linear interpolation between l and r by segment index
            dQ_out(s) = dQ_raw(l) + (dQ_raw(r) - dQ_raw(l)) * (s - l) / (r - l);
        elseif ~isempty(left_valid)
            dQ_out(s) = dQ_raw(left_valid(end));   % extrapolate: copy nearest left
        elseif ~isempty(right_valid)
            dQ_out(s) = dQ_raw(right_valid(1));    % extrapolate: copy nearest right
        end
        % else: remains NaN (handled above by empty valid_idx check)
    end
end

function segs = local_find_segments(mask)

    segs = [];
    n = length(mask);
    i = 1;
    while i <= n
        if mask(i)
            j = i;
            while j < n && mask(j+1), j = j + 1; end
            segs = [segs; i, j];
            i = j + 1;
        else
            i = i + 1;
        end
    end
end
