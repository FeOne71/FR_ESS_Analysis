% RPT_Field_Estimation_SVM.m
% 1. 학습 완료된 100% Data 기반 SVM 실생성 모델 로드
% 2. 필드 데이터(Field Data, 2021~2025) 로드 및 20개 물리 피처 추출
% 3. 학습된 SVM 모델을 사용하여 필드 데이터의 SOH, LLI, LAM, SOP 추정 및 시각화

clear; clc; close all;
warning off;

%% ========================================================================
% Section 1: 모델 로드 (Support Vector Machine)
% ========================================================================
fprintf('=== Section 1: Loading Trained SVM Model ===\n');
baseDir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
svmDir  = fullfile(baseDir, 'ML_SVM');
modelFile = fullfile(svmDir, 'Result_SVM.mat');
rulerPath = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', 'FeatureEngineering', 'MasterRulers_v3.mat');

if ~exist(modelFile, 'file'), error('SVM Model not found.'); end
load(modelFile, 'Results_SVM');
load(rulerPath, 'MasterRulers');

feature_names = Results_SVM.feature_names;
label_names = Results_SVM.label_names;

% Extract the 100% Final Models from the Results struct (identically trained)
SVM_Models = struct();
for i = 1:length(label_names)
    lbl = label_names{i};
    SVM_Models.(lbl) = Results_SVM.(lbl).FinalModel;
end

%% ========================================================================
% Section 2: 필드 데이터 로드 (KIMJ ESS 2021~2025)
% Segment logic: identical to Visualize_VIT_byYear.m / FieldQmax_dQdV.m
% ========================================================================
fprintf('\n=== Section 2: Loading Field Data (2021~2025) ===\n');
rawDir   = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'Rack_raw2mat');
dataFiles = {
    fullfile(rawDir, 'Old', '2021', '202106', 'Raw_20210603.mat'), 'Y2021', 'old',  datetime(2021,6,3);
    fullfile(rawDir, 'New', '2023', '202310', 'Raw_20231016.mat'), 'Y2023', 'new',  datetime(2023,10,16);
    fullfile(rawDir, 'New', '2024', '202409', 'Raw_20240909.mat'), 'Y2024_Auto', 'new',  datetime(2024,9,9);
    fullfile(rawDir, 'New', '2024', '202409', 'Raw_20240909.mat'), 'Y2024_Manual', 'new',  datetime(2024,9,9);
    % Y2025 excluded: voltage range too narrow
};

% Per-year minimum segment durations (seconds)
min_charge_secs    = [600, 300, 300, 300];   % Y2021, Y2023, Y2024_Auto, Y2024_Manual
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
        if isfield(S, 'Raw') && isfield(S.Raw, 'Rack01')
            D = S.Raw.Rack01;
            if isfield(D, 'Time'), t = datetime(D.Time);
            elseif isfield(D, 'Date_Time')
                if isduration(D.Date_Time), t = base_date + D.Date_Time;
                else, t = datetime(D.Date_Time); end
            else, continue; end
            I_rack = D.DCCurrent_A(:);
            V_avg  = D.AverageCV_V(:);
            if isfield(D, 'AverageCT_C'), T_avg = D.AverageCT_C(:);
            elseif isfield(D, 'AverageMT_degC'), T_avg = D.AverageMT_degC(:);
            else, T_avg = nan(size(I_rack)); end
        else, continue; end
    else
        D = S.Raw;
        if isduration(D.Date_Time), t = base_date + D.Date_Time;
        else, t = datetime(D.Date_Time); end
        I_rack = D.DCCurrent(:);
        V_avg  = D.CVavg(:);
        if isfield(D, 'MTavg'), T_avg = D.MTavg(:);
        else, T_avg = nan(size(I_rack)); end
    end
    
    tsec   = seconds(t - t(1));
    I_cell = I_rack / Np;
    
    % Store full profiles for plotting
    FieldData.(year_label).Full.V = V_avg;
    FieldData.(year_label).Full.I = I_cell;
    FieldData.(year_label).Full.Time = tsec;
    
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
    
    % Data Segment Visualization
    fig_seg = figure('Name', sprintf('[%s] Field Selection', year_label), 'Position', [100, 100, 1000, 600], 'Visible', 'off');
    sgtitle(sprintf('KIMJ ESS %s - Daily Profile & Selected Segments', year_label), 'FontSize', 14, 'FontWeight', 'bold');
    
    t_full = FieldData.(year_label).Full.Time / 3600; % hours
    
    % Subplot 1: Voltage
    ax1 = subplot(2, 1, 1); hold on; grid on; box on;
    plot(t_full, FieldData.(year_label).Full.V, 'k-', 'LineWidth', 1);
    ylabel('Voltage (V)'); title('Voltage Profile');
    y_lim_v = ylim;
    if ~isnan(chg_s), patch([chg_s, chg_e, chg_e, chg_s]/3600, [y_lim_v(1), y_lim_v(1), y_lim_v(2), y_lim_v(2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); end
    if ~isnan(dch_s), patch([dch_s, dch_e, dch_e, dch_s]/3600, [y_lim_v(1), y_lim_v(1), y_lim_v(2), y_lim_v(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); end
    
    % Subplot 2: Current
    ax2 = subplot(2, 1, 2); hold on; grid on; box on;
    plot(t_full, FieldData.(year_label).Full.I, 'k-', 'LineWidth', 1);
    ylabel('Current (A)'); xlabel('Time (Hours)'); title('Current Profile');
    y_lim_i = ylim;
    if ~isnan(chg_s), patch([chg_s, chg_e, chg_e, chg_s]/3600, [y_lim_i(1), y_lim_i(1), y_lim_i(2), y_lim_i(2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); end
    if ~isnan(dch_s), patch([dch_s, dch_e, dch_e, dch_s]/3600, [y_lim_i(1), y_lim_i(1), y_lim_i(2), y_lim_i(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); end
    legend('Raw Data', 'Extracted Charge', 'Extracted Discharge', 'Location', 'best');
    linkaxes([ax1, ax2], 'x'); xlim([0, 24]);
    
    saveas(fig_seg, fullfile(svmDir, sprintf('Segment_Selection_%s.fig', year_label)));
    close(fig_seg);
end

%% ========================================================================
% Section 3: 피처 추출 및 상태 추정 (SVM)
% ========================================================================
fprintf('\n=== Section 3: Feature Extraction & SVM Estimation ===\n');
years = fieldnames(FieldData);
Results_Ev = struct();

win_chg = [3.674, 3.976]; win_dch = [3.644, 3.872]; ma_window = 21; Q_0 = 64; 

for k = 1:length(years)
    yr = years{k};
    if ~isfield(FieldData.(yr), 'Chg') || ~isfield(FieldData.(yr), 'Dch')
        fprintf('  %s: Missing Chg/Dch Segment. Skipped.\n', yr); continue;
    end
    [v_c_u, uidc] = unique(FieldData.(yr).Chg.V, 'stable'); 
    q_c_u = FieldData.(yr).Chg.Q(uidc);
    
    dQ_chg = nan(1, 5);
    try 
        % Create standard 0.001V grid based on Master Ruler boundaries
        v_grid_c = (VR_chg(1) - 0.05) : 0.001 : (VR_chg(end) + 0.05);
        
        % Interpolate Q values onto the 0.001V grid (Q as function of V)
        % Warning: V must be monotonically increasing for charge
        [v_c_sorted, sort_idx] = sort(v_c_u);
        q_c_sorted = q_c_u(sort_idx);
        [v_c_strict, unique_idx] = unique(v_c_sorted);
        q_c_strict = q_c_sorted(unique_idx);
        
        Q_interp_c = interp1(v_c_strict, q_c_strict, v_grid_c, 'linear', 'extrap');
        
        % Apply Moving Average (Window 30) on the interpolated Q-V curve
        Q_sm_c = movmean(Q_interp_c, 30);
        
        % Extract dQ at MasterRuler bounds
        QR_chg = nan(1, length(VR_chg));
        for b = 1:length(VR_chg)
            [~, min_idx] = min(abs(v_grid_c - VR_chg(b)));
            QR_chg(b) = Q_sm_c(min_idx);
        end
        dQ_chg = abs(diff(QR_chg)); 
    catch ME
        fprintf('  %s Chg interpolation failed: %s\n', yr, ME.message);
    end
    
    % Discharge segment Processing
    [v_d_u, uidd] = unique(FieldData.(yr).Dch.V, 'stable'); 
    q_d_u = FieldData.(yr).Dch.Q(uidd);
    
    dQ_dch = nan(1, 5);
    try 
        % Create standard 0.001V grid (descending for discharge)
        v_grid_d = (VR_dch(1) + 0.05) : -0.001 : (VR_dch(end) - 0.05);
        
        % Interpolate Q values onto the 0.001V grid
        [v_d_sorted, sort_idx] = sort(v_d_u, 'descend');
        q_d_sorted = q_d_u(sort_idx);
        [v_d_strict, unique_idx] = unique(v_d_sorted, 'stable');
        q_d_strict = q_d_sorted(unique_idx);
        
        Q_interp_d = interp1(v_d_strict, q_d_strict, v_grid_d, 'linear', 'extrap');
        
        % Apply Moving Average (Window 30)
        Q_sm_d = movmean(Q_interp_d, 30);
        
        % Extract dQ at MasterRuler bounds
        QR_dch = nan(1, length(VR_dch));
        for b = 1:length(VR_dch)
            [~, min_idx] = min(abs(v_grid_d - VR_dch(b)));
            QR_dch(b) = Q_sm_d(min_idx);
        end
        dQ_dch = abs(diff(QR_dch));
    catch ME
        fprintf('  %s Dch interpolation failed: %s\n', yr, ME.message);
    end
    
    % -- 2. Peak Extraction
    PkH_chg = NaN; PkA_chg = NaN; PkPos_chg = NaN;
    PkH_dch = NaN; PkA_dch = NaN; PkPos_dch = NaN;
    
    [PkH_chg, PkA_chg, PkPos_chg] = field_extract_peak(v_c_u, q_c_u, win_chg(1), win_chg(2), ma_window);
    [PkH_dch, PkA_dch, PkPos_dch] = field_extract_peak(v_d_u, q_d_u, win_dch(1), win_dch(2), ma_window);
    
    % -- 3. Energy Extraction (Wh)
    Energy_dch = NaN;
    try
        idx_d = FieldData.(yr).Dch.V >= win_dch(1) & FieldData.(yr).Dch.V <= win_dch(2);
        if sum(idx_d) > 5
            Energy_dch = abs(trapz(FieldData.(yr).Dch.Q(idx_d), FieldData.(yr).Dch.V(idx_d)));
        end
    catch, end
    
    % -- 4. C_eff and Temp
    C_eff_chg = mean(abs(FieldData.(yr).Chg.I)) / Q_0;
    C_eff_dch = mean(abs(FieldData.(yr).Dch.I)) / Q_0;
    
    T_chg_avg = mean(FieldData.(yr).Chg.T, 'omitnan');
    T_dch_avg = mean(FieldData.(yr).Dch.T, 'omitnan');
    
    % Store C_eff for later reference
    Results_Ev.(yr).Ceff.chg = C_eff_chg;
    Results_Ev.(yr).Ceff.dch = C_eff_dch;
    
    % ================================================================
    % SECTION 3B: Split Normalization (v2 신규 NormStats 구조 사용)
    %   chg 피처 → C_eff_chg 기준 그룹 mu_chg/sigma_chg
    %   dch 피처 → C_eff_dch 기준 그룹 mu_dch/sigma_dch
    %   ceff 피처 → raw 유지
    % ================================================================
    NormStats      = Results_SVM.NormStats;
    crate_groups   = NormStats.crate_groups;
    crate_vals_num = NormStats.crate_vals_num;

    % Feature indices (19개, T_avg 제외)
    chg_idx  = [1:5, 11, 13, 15];        % dQ_chg×5, PkH_chg, PkA_chg, PkPos_chg (8개)
    dch_idx  = [6:10, 12, 14, 16, 17];   % dQ_dch×5, PkH_dch, PkA_dch, PkPos_dch, Energy_dch (9개)
    ceff_idx = [18, 19];                  % C_eff_chg, C_eff_dch (정규화 없음)

    X_raw = [dQ_chg, dQ_dch, PkH_chg, PkH_dch, PkA_chg, PkA_dch, PkPos_chg, PkPos_dch, Energy_dch, C_eff_chg, C_eff_dch];

    % NaN check
    if any(isnan(X_raw))
        n_nan = sum(isnan(X_raw));
        fprintf('  %s: %d NaN features remain. Replacing with 0.\n', yr, n_nan);
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
    X_norm(chg_idx) = (X_raw(chg_idx) - mu_chg)  ./ (sig_chg + eps);
    X_norm(dch_idx) = (X_raw(dch_idx) - mu_dch)  ./ (sig_dch + eps);
    % ceff_idx: raw 유지

    fprintf('  %s | C_eff_chg=%.3fC→%s(w=%.0f%%)  C_eff_dch=%.3fC→%s(w=%.0f%%)\n', ...
        yr, C_eff_chg, cr1_c, w1c*100, C_eff_dch, cr1_d, w1d*100);

    Results_Ev.(yr).X_field = X_norm;

    % Predict using final SVM models
    for i = 1:length(label_names)
        lbl = label_names{i};
        Results_Ev.(yr).(lbl) = predict(SVM_Models.(lbl), X_norm);
    end
    
    % Visualize the extracted 20 features for this year
    fig_fea = figure('Name', sprintf('[%s] Extracted Features', yr), 'Position', [150, 150, 1000, 400], 'Visible', 'off');
    bar(X_raw, 'FaceColor', [0.8 0.4 0.2]); grid on;
    title(sprintf('Extracted Physics Features - KIMJ %s', yr), 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Feature Value');
    xticks(1:length(feature_names)); xticklabels(feature_names); xtickangle(45);
    saveas(fig_fea, fullfile(svmDir, sprintf('Features_%s.fig', yr)));
    close(fig_fea);
end

fprintf('=== Field Estimation Strategy Executed (SVM) ===\n');

%% ========================================================================
% Section 4: 결과 시각화
% ========================================================================
fig_svm_field = figure('Name', 'SVM Field Estimation', 'Position', [150, 150, 1200, 800]);
sgtitle('SVM (Gaussian) - Field Trajectory (KIMJ ESS, 2021-2025)', 'FontSize', 16, 'FontWeight', 'bold');

yr_strs = fieldnames(Results_Ev);
n_yrs = length(yr_strs);
colors = lines(length(label_names));

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

line_colors = [0.2 0.4 0.8; 0.9 0.3 0.1; 0.1 0.7 0.3];
y_labels    = {'SOH (%)', 'LLI (%)', 'LAM (%)'};

for i = 1:length(label_names)
    lbl = label_names{i};
    subplot(1, 3, i); hold on; box on; grid on;
    
    y_vals = zeros(n_yrs, 1);
    for k = 1:n_yrs
        y_vals(k) = Results_Ev.(yr_strs{k}).(lbl);
    end
    
    plot(yr_nums, y_vals, '-o', 'Color', line_colors(i,:), 'LineWidth', 2, ...
        'MarkerSize', 8, 'MarkerFaceColor', line_colors(i,:));

    % Data point labels
    for k = 1:n_yrs
        text(yr_nums(k), y_vals(k), sprintf('  %.2f', y_vals(k)), ...
            'VerticalAlignment', 'bottom', 'FontSize', 9, 'FontWeight', 'bold', ...
            'Color', line_colors(i,:)*0.7);
    end

    % X-axis labels
    x_labels = cell(n_yrs, 1);
    for k = 1:n_yrs
        yr_name = strrep(yr_strs{k}, 'Y', ''); 
        yr_name = strrep(yr_name, '_', ' '); 
        x_labels{k} = sprintf('%s', yr_name);
    end
    set(gca, 'XTick', yr_nums, 'XTickLabel', x_labels, 'FontSize', 9);
    xtickangle(45);
    
    xlabel('Year'); ylabel(y_labels{i});
    title(lbl, 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
    xlim([min(yr_nums)-0.5, max(yr_nums)+0.5]);
end

saveas(fig_svm_field, fullfile(svmDir, 'SVM_Field_Trajectory.fig'));

%% Helper Function
function [pk_height, pk_area, pk_pos, V_out, dQdV_out] = field_extract_peak(V_u, Q_u, minV, maxV, ma_win)
    pk_height = NaN; pk_area = NaN; pk_pos = NaN; V_out = []; dQdV_out = [];
    try
        % 1. Calculate native physical diffs directly from full raw points
        dV = diff(V_u);
        dQ = diff(Q_u);
        
        if length(dV) > 5
            % 2. Avoid division by zero artifacts and compute raw dQ/dV
            valid_idx = abs(dV) > 0;
            V_mid_raw = V_u(1:end-1);
            V_mid = V_mid_raw(valid_idx) + dV(valid_idx)/2;
            dQdV_raw = abs(dQ(valid_idx)) ./ abs(dV(valid_idx));
            
            % 3. Smooth the derivative natively over the full range
            dQdV_sm_raw = movmean(dQdV_raw, ma_win);
            
            % 4. Create full 0.001V grid for visualization
            V_grid_full = (ceil(min(V_mid)*1000)/1000 : 0.001 : floor(max(V_mid)*1000)/1000)';
            
            % 5. Map the safe, smoothed derivative onto full grid
            [V_u_mid, id_m] = unique(V_mid, 'stable');
            dQdV_sm_u = dQdV_sm_raw(id_m);
            
            [V_u_mid, sort_idx] = sort(V_u_mid);
            dQdV_sm_u = dQdV_sm_u(sort_idx);
            
            if length(V_u_mid) > 1
                dQdV_grid_full = interp1(V_u_mid, dQdV_sm_u, V_grid_full, 'linear', 'extrap');
                
                % V_out and dQdV_out contain the ENTIRE curve for plotting
                V_out = V_grid_full; 
                dQdV_out = dQdV_grid_full;
                
                % 6. Extract Peaks strictly within MasterRuler boundaries
                idx_win = V_grid_full >= minV & V_grid_full <= maxV;
                V_win = V_grid_full(idx_win);
                dQdV_win = dQdV_grid_full(idx_win);
                
                if ~isempty(dQdV_win)
                    [pk_height, max_idx] = max(dQdV_win);
                    pk_pos = V_win(max_idx);
                    pk_area = abs(trapz(V_win, dQdV_win));
                end
            end
        end
    catch
    end
end

function segs = local_find_segments(tf_array)
    % tf_array: logical array (e.g., isChg or isDch)
    % returns: Nx2 matrix of [start_idx, end_idx] for contiguous true blocks
    diff_arr = diff([0; tf_array(:); 0]);
    idx_start = find(diff_arr == 1);
    idx_end   = find(diff_arr == -1) - 1;
    segs = [idx_start, idx_end];
end

