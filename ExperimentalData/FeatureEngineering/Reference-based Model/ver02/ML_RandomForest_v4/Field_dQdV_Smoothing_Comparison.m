% Field_dQdV_Smoothing_Comparison.m
% =====================================================================
% 필드 데이터 dQ/dV 스무딩 비교 시각화
% 각 연도별 1×2 (충전, 방전), 3가지 스무딩 방법 오버레이
% =====================================================================
clear; clc; close all;

rawDir   = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'Rack_raw2mat');
path_master_ruler = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', ...
    'ExperimentalData', 'FeatureEngineering', 'MasterRulers_v3.mat');
load(path_master_ruler, 'MasterRulers');

saveDir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
    'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02', 'ML_RandomForest_v4');

dataFiles = {
    fullfile(rawDir, 'Old', '2021', '202106', 'Raw_20210603.mat'), 'Y2021', 'old',  datetime(2021,6,3);
    fullfile(rawDir, 'New', '2023', '202310', 'Raw_20231016.mat'), 'Y2023', 'new',  datetime(2023,10,16);
    fullfile(rawDir, 'New', '2024', '202409', 'Raw_20240909.mat'), 'Y2024', 'new',  datetime(2024,9,9);
    fullfile(rawDir, 'New', '2025', '202507', 'Raw_20250711.mat'), 'Y2025', 'new',  datetime(2024,9,9);
};

min_charge_secs    = [600, 300, 300, 300];
min_discharge_secs = [300, 150, 300, 300];
Np = 2; C_cell_Ah = 64; thr_A = C_cell_Ah * 0.05;

for k = 1:size(dataFiles, 1)
    fpath      = dataFiles{k, 1};
    year_label = dataFiles{k, 2};
    dataType   = dataFiles{k, 3};
    base_date  = dataFiles{k, 4};
    min_chg_sec = min_charge_secs(k);
    min_dch_sec = min_discharge_secs(k);

    if ~exist(fpath, 'file'), continue; end
    fprintf('Processing %s...\n', year_label);

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
    else
        D = S.Raw;
        if isduration(D.Date_Time), t = base_date + D.Date_Time;
        else, t = datetime(D.Date_Time); end
        I_rack = D.DCCurrent(:);
        V_avg  = D.CVavg(:);
    end

    tsec   = seconds(t - t(1));
    I_cell = I_rack / Np;

    % --- 충전/방전 세그먼트 추출 ---
    chg_mask = I_cell > thr_A;
    dch_mask = I_cell < -thr_A;

    chg_s = NaN; chg_e = NaN;
    chg_diff = diff([0; chg_mask(:); 0]);
    chg_starts = find(chg_diff == 1);
    chg_ends   = find(chg_diff == -1) - 1;
    if ~isempty(chg_starts)
        chg_lens = chg_ends - chg_starts + 1;
        valid_segs = find(chg_lens >= min_chg_sec);
        if ~isempty(valid_segs)
            [~, mx] = max(chg_lens(valid_segs));
            idx = valid_segs(mx);
            chg_s = chg_starts(idx); chg_e = chg_ends(idx);
        end
    end

    dch_s = NaN; dch_e = NaN;
    dch_diff = diff([0; dch_mask(:); 0]);
    dch_starts = find(dch_diff == 1);
    dch_ends   = find(dch_diff == -1) - 1;
    if ~isempty(dch_starts)
        dch_lens = dch_ends - dch_starts + 1;
        valid_segs = find(dch_lens >= min_dch_sec);
        if ~isempty(valid_segs)
            [~, mx] = max(dch_lens(valid_segs));
            idx = valid_segs(mx);
            dch_s = dch_starts(idx); dch_e = dch_ends(idx);
        end
    end

    % === Figure: 1×2 (충전, 방전) — 3 methods 오버레이 ===
    fig = figure('Name', sprintf('dQdV Smoothing - %s', year_label), ...
                 'Position', [50, 100, 1400, 500], 'Visible', 'off');

    % --- Charge ---
    subplot(1, 2, 1);
    if ~isnan(chg_s)
        V_raw = V_avg(chg_s:chg_e);
        I_raw = I_cell(chg_s:chg_e);
        t_seg = tsec(chg_s:chg_e);
        Q_raw = cumtrapz(t_seg, abs(I_raw)) / 3600;

        % Method 1: Raw
        [V1, dQdV1] = compute_dQdV(V_raw, Q_raw, true);
        plot(V1, dQdV1, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'DisplayName', '1) Raw');
        hold on;

        % Method 2: V,Q movmean(60)
        V_sm = movmean(V_raw, 60);
        Q_sm = movmean(Q_raw, 60);
        [V2, dQdV2] = compute_dQdV(V_sm, Q_sm, true);
        plot(V2, dQdV2, 'b-', 'LineWidth', 1, 'DisplayName', '2) V,Q ma(60)');

        % Method 3: V,Q movmean(60) + dQdV movmean(30)
        dQdV3 = movmean(dQdV2, 30);
        plot(V2, dQdV3, 'r-', 'LineWidth', 1.5, 'DisplayName', '3) V,Q ma(60) + dQdV ma(30)');

        ylim_auto([dQdV2; dQdV3]);
    end
    title(sprintf('Charge — %s', year_label), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('V'); ylabel('dQ/dV (Ah/V)');
    legend('Location', 'best', 'FontSize', 9); grid on;

    % --- Discharge ---
    subplot(1, 2, 2);
    if ~isnan(dch_s)
        V_raw = V_avg(dch_s:dch_e);
        I_raw = I_cell(dch_s:dch_e);
        t_seg = tsec(dch_s:dch_e);
        Q_raw = cumtrapz(t_seg, abs(I_raw)) / 3600;

        % Method 1: Raw
        [V1, dQdV1] = compute_dQdV(V_raw, Q_raw, false);
        plot(V1, dQdV1, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'DisplayName', '1) Raw');
        hold on;

        % Method 2: V,Q movmean(60)
        V_sm = movmean(V_raw, 60);
        Q_sm = movmean(Q_raw, 60);
        [V2, dQdV2] = compute_dQdV(V_sm, Q_sm, false);
        plot(V2, dQdV2, 'b-', 'LineWidth', 1, 'DisplayName', '2) V,Q ma(60)');

        % Method 3: V,Q movmean(60) + dQdV movmean(30)
        dQdV3 = movmean(dQdV2, 30);
        plot(V2, dQdV3, 'r-', 'LineWidth', 1.5, 'DisplayName', '3) V,Q ma(60) + dQdV ma(30)');

        ylim_auto([dQdV2; dQdV3]);
    end
    title(sprintf('Discharge — %s', year_label), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('V'); ylabel('dQ/dV (Ah/V)');
    legend('Location', 'best', 'FontSize', 9); grid on;

    sgtitle(sprintf('[%s] dQ/dV Smoothing Comparison (3 Methods Overlaid)', year_label), ...
            'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig, fullfile(saveDir, sprintf('dQdV_Smoothing_%s.fig', year_label)));
    close(fig);
    fprintf('  Saved: dQdV_Smoothing_%s.fig\n', year_label);
end

fprintf('\nDone.\n');

%% Helper
function [V_out, dQdV_out] = compute_dQdV(V_in, Q_in, is_ascending)
    [V_u, uid] = unique(V_in(:), 'stable');
    Q_u = Q_in(uid);
    if is_ascending
        mono = true(size(V_u));
        for ii = 2:length(V_u)
            if V_u(ii) <= V_u(ii-1), mono(ii) = false; end
        end
    else
        mono = true(size(V_u));
        for ii = 2:length(V_u)
            if V_u(ii) >= V_u(ii-1), mono(ii) = false; end
        end
    end
    V_u = V_u(mono); Q_u = Q_u(mono);
    if length(V_u) < 5
        V_out = V_u; dQdV_out = zeros(size(V_u)); return;
    end
    if ~is_ascending
        V_u = flipud(V_u); Q_u = flipud(Q_u);
    end
    dV = gradient(V_u); dQ = gradient(Q_u);
    dV(dV == 0) = NaN;
    dQdV = dQ ./ dV;
    dQdV(isinf(dQdV) | isnan(dQdV)) = 0;
    V_out = V_u; dQdV_out = dQdV;
end

function ylim_auto(data)
    data = data(~isnan(data) & ~isinf(data));
    if isempty(data), return; end
    q1 = prctile(data, 2); q99 = prctile(data, 98);
    margin = (q99 - q1) * 0.15;
    if margin < 0.1, margin = 0.1; end
    ylim([q1 - margin, q99 + margin]);
end
