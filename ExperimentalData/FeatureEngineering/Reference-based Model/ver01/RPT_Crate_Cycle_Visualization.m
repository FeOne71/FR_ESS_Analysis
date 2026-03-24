%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT_Crate_Cycle_Visualization.m
%
% FieldQmax_dQdV.m 스크립트처럼: 충전·방전 각각 별도 figure.
% - Figure 1: Charge 전용 — 모든 사이클·모든 채널 충전 데이터 → V-t, I-t, P-t, Temp-t
% - Figure 2: Discharge 전용 — 모든 사이클·모든 채널 방전 데이터 → V-t, I-t, P-t, Temp-t
% 데이터: RPT_VQ_grid.mat (OCV_integrated). P = V*I 유도. 온도: 세그먼트별 T1_raw [℃] 사용.
% 범례: 사이클만 표시. 색상: cyc0=2021, cyc200=2023, cyc400=2024, cyc600=2025, 나머지 구분.
% 참고: Reference/FieldQmax_dQdV.m (Charge All Years / Discharge All Years 분리)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning on;

%% ========================================================================
% 1. 경로 및 파라미터
% ========================================================================
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_rpt_vq_mat = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');

saveDir = fullfile(baseDir, 'FeatureEngineering', 'Lab_RPT_Analysis', 'Crate_Cycle_Visualizations');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% RPT_VQ_grid 세그먼트: c01_charge, c01_discharge, c05, c1, c2, c3 충방전 전부 사용
crate_prefix_list = {'c01','c05','c1','c2','c3'};
apply_smoothing = true;
smooth_window = 7;

% FieldQmax_dQdV.m 연도 색상: 2021 Blue, 2023 Green, 2024 Orange, 2025 Purple
year_colors_RGB = [0 0 1; 0 0.8 0; 1 0.5 0; 0.8 0 0.8];
cycle_to_year = [0, 200, 400, 600];  % cyc0→2021, cyc200→2023, cyc400→2024, cyc600→2025
% 그 외 사이클용 추가 색상
other_colors = [0.5 0.5 0.5; 0 0.7 0.7; 0.6 0.35 0.2; 0.9 0.5 0.1; 0.4 0.2 0.6; 0.2 0.5 0.5];

fprintf('Loading RPT_VQ_grid from:\n  %s\n', path_rpt_vq_mat);
S = load(path_rpt_vq_mat, 'RPT_VQ_grid');
RPT_VQ_grid = S.RPT_VQ_grid;

cyc_fields = fieldnames(RPT_VQ_grid);
fprintf('Found %d cycle fields: %s\n', numel(cyc_fields), strjoin(cyc_fields', ', '));

%% 사이클별 색상 할당 (cyc0→2021색, cyc200→2023, cyc400→2024, cyc600→2025, 나머지 구분)
cycle_colors = zeros(numel(cyc_fields), 3);
other_idx = 0;
for c_idx = 1:numel(cyc_fields)
    cyc_name = cyc_fields{c_idx};
    num_str = regexp(cyc_name, '\d+', 'match');
    if isempty(num_str), cyc_num = NaN; else, cyc_num = str2double(num_str{1}); end
    [~, year_idx] = ismember(cyc_num, cycle_to_year);
    if year_idx > 0
        cycle_colors(c_idx, :) = year_colors_RGB(year_idx, :);
    else
        other_idx = other_idx + 1;
        idx_other = mod(other_idx - 1, size(other_colors, 1)) + 1;
        cycle_colors(c_idx, :) = other_colors(idx_other, :);
    end
end

%% ========================================================================
% 2. 모든 채널 × 모든 사이클 Charge/Discharge 데이터 수집
% ========================================================================
% 채널 목록: 첫 번째 사이클에서 Ch* 필드만 사용
cyc1 = cyc_fields{1};
all_field_names = fieldnames(RPT_VQ_grid.(cyc1));
channels = all_field_names(strncmp(all_field_names, 'Ch', 2));
fprintf('Using %d channels: %s\n', numel(channels), strjoin(channels', ', '));

Charge_data_list = [];
Discharge_data_list = [];

for c_idx = 1:numel(cyc_fields)
    cyc_name = cyc_fields{c_idx};
    if ~isfield(RPT_VQ_grid, cyc_name), continue; end
    cyc_struct = RPT_VQ_grid.(cyc_name);

    for cr_idx = 1:numel(crate_prefix_list)
        crate_prefix = crate_prefix_list{cr_idx};
        chg_name = [crate_prefix '_charge'];
        dch_name = [crate_prefix '_discharge'];

        for ch_idx = 1:numel(channels)
            ch_name = channels{ch_idx};
            if ~isfield(cyc_struct, ch_name), continue; end
            ch_data = cyc_struct.(ch_name);
            if ~isfield(ch_data, chg_name) || ~isfield(ch_data, dch_name), continue; end

            % Charge
            data_chg = ch_data.(chg_name);
            if ~isfield(data_chg, 'V_grid') || ~isfield(data_chg, 'Q') || ~isfield(data_chg, 't'), continue; end
            V_chg = data_chg.V_grid(:); Q_chg = data_chg.Q(:); t_chg = data_chg.t(:);
            if numel(V_chg) < 3 || numel(t_chg) ~= numel(Q_chg), continue; end

            % Discharge (같은 채널·같은 C-rate에서 둘 다 있을 때만 추가)
            data_dch = ch_data.(dch_name);
            if ~isfield(data_dch, 'V_grid') || ~isfield(data_dch, 'Q') || ~isfield(data_dch, 't'), continue; end
            V_dch = data_dch.V_grid(:); Q_dch = data_dch.Q(:); t_dch = data_dch.t(:);
            if numel(V_dch) < 3 || numel(t_dch) ~= numel(Q_dch), continue; end

            [V_dQdV_chg, dQdV_chg] = local_compute_dQdV(V_chg, Q_chg, apply_smoothing, smooth_window);
            I_chg = local_derive_current(Q_chg, t_chg);
            idx = length(Charge_data_list) + 1;
            Charge_data_list(idx).V_dQdV = V_dQdV_chg; %#ok<AGROW>
            Charge_data_list(idx).dQdV = dQdV_chg;
            Charge_data_list(idx).V_raw = V_chg;
            Charge_data_list(idx).Q_raw = Q_chg;
            Charge_data_list(idx).I_raw = I_chg;
            Charge_data_list(idx).t_raw = t_chg;
            Charge_data_list(idx).P_raw = V_chg .* I_chg;
            [t_temp_chg, T_temp_chg] = local_get_temp_series(data_chg, t_chg);
            Charge_data_list(idx).t_raw_temp = t_temp_chg;
            Charge_data_list(idx).Temp_raw = T_temp_chg;
            Charge_data_list(idx).cycle_name = cyc_name;
            Charge_data_list(idx).crate_name = chg_name; % e.g. 'c05_charge'

            [V_dQdV_dch, dQdV_dch] = local_compute_dQdV(V_dch, Q_dch, apply_smoothing, smooth_window);
            I_dch = local_derive_current(Q_dch, t_dch);
            Discharge_data_list(idx).V_dQdV = V_dQdV_dch; %#ok<AGROW>
            Discharge_data_list(idx).dQdV = dQdV_dch;
            Discharge_data_list(idx).V_raw = V_dch;
            Discharge_data_list(idx).Q_raw = Q_dch;
            Discharge_data_list(idx).I_raw = I_dch;
            Discharge_data_list(idx).t_raw = t_dch;
            Discharge_data_list(idx).P_raw = V_dch .* I_dch;
            [t_temp_dch, T_temp_dch] = local_get_temp_series(data_dch, t_dch);
            Discharge_data_list(idx).t_raw_temp = t_temp_dch;
            Discharge_data_list(idx).Temp_raw = T_temp_dch;
            Discharge_data_list(idx).cycle_name = cyc_name;
            Discharge_data_list(idx).crate_name = dch_name;
        end
    end
end
if isempty(Charge_data_list)
    fprintf('No charge/discharge data found. Check channel and crate.\n');
    return;
end

% T-Time: (사이클×C-rate)별 8채널 평균. 공통 시간 = 8채널 모두 유효한 구간(1번)으로 끝 급상승 방지
Charge_T_avg_list = local_avg_temp_by_cycle_crate(Charge_data_list);
Discharge_T_avg_list = local_avg_temp_by_cycle_crate(Discharge_data_list);

%% ========================================================================
% 3. Field_Visualize_VIT_byYear.m 형식: 충전/방전 각각 1 fig, 2x2 = V-Time, I-Time, Power-Time, T-Time
%    V,I,P: 모든 채널 scatter / T-Time: (사이클×C-rate)별 8채널 평균 1곡선
% ========================================================================
msz = 3;   % scatter marker size (C-rate 5개×사이클×채널 많아서 겹침 방지로 작게)
% ----- Figure 1: Charge 전용 (모든 사이클/채널 취합) -----
fig_charge = figure('Name', 'VIT Charge (all cycles combined)', 'NumberTitle', 'off');
tl_chg = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(tl_chg, 1); hold(ax1, 'on'); grid(ax1, 'on'); xlabel(ax1, 'Time [s]'); ylabel(ax1, 'V [V]'); title(ax1, 'V - Time');
ax2 = nexttile(tl_chg, 2); hold(ax2, 'on'); grid(ax2, 'on'); xlabel(ax2, 'Time [s]'); ylabel(ax2, 'I [A]'); title(ax2, 'I - Time');
ax3 = nexttile(tl_chg, 3); hold(ax3, 'on'); grid(ax3, 'on'); xlabel(ax3, 'Time [s]'); ylabel(ax3, 'Power [W]'); title(ax3, 'Power - Time');
ax4 = nexttile(tl_chg, 4); hold(ax4, 'on'); grid(ax4, 'on'); xlabel(ax4, 'Time [s]'); ylabel(ax4, 'T [°C]'); title(ax4, 'T - Time (8ch avg)');
axChg = [ax1, ax2, ax3, ax4];   % 충전 4개 서브플롯 핸들 (방전과 축 범위 통일용)
cycle_shown_chg = {};

for k = 1:length(Charge_data_list)
    D = Charge_data_list(k);
    cname = D.cycle_name;
    [~, cidx] = ismember(cname, cyc_fields);
    cidx = max(1, cidx);
    col = cycle_colors(cidx, :);
    show_legend = ~ismember(cname, cycle_shown_chg);
    if show_legend, cycle_shown_chg{end+1} = cname; end %#ok<AGROW>
    t_sec = local_time_to_seconds(D.t_raw);
    hv = {'off','on'};
    scatter(ax1, t_sec, D.V_raw, msz, col, 'filled', 'DisplayName', cname, 'HandleVisibility', hv{show_legend+1});
    scatter(ax2, t_sec, D.I_raw, msz, col, 'filled', 'HandleVisibility', 'off');
    scatter(ax3, t_sec, D.P_raw, msz, col, 'filled', 'HandleVisibility', 'off');
end
cycle_shown_T = {};
for g = 1:length(Charge_T_avg_list)
    A = Charge_T_avg_list(g);
    [~, cidx] = ismember(A.cycle_name, cyc_fields);
    cidx = max(1, cidx);
    col = cycle_colors(cidx, :);
    show_leg = ~ismember(A.cycle_name, cycle_shown_T);
    if show_leg, cycle_shown_T{end+1} = A.cycle_name; end %#ok<AGROW>
    hv = {'off','on'};
    plot(ax4, A.t_sec, A.T_avg, '-', 'Color', col, 'LineWidth', 1.5, 'DisplayName', A.cycle_name, 'HandleVisibility', hv{show_leg+1});
end
legend(ax1, 'Location', 'best');
legend(ax4, 'Location', 'best');
title(tl_chg, 'Charge (all cycles)');
uq_crates_chg = unique({Charge_data_list.crate_name});
for cr = 1:numel(uq_crates_chg)
    crate = uq_crates_chg{cr};
    idx = find(strcmp({Charge_data_list.crate_name}, crate));
    t_all = []; V_all = []; I_all = []; P_all = [];
    for i = idx
        D = Charge_data_list(i);
        ts = local_time_to_seconds(D.t_raw);
        t_all = [t_all; ts(:)]; V_all = [V_all; D.V_raw(:)]; I_all = [I_all; D.I_raw(:)]; P_all = [P_all; D.P_raw(:)];
    end
    pre = crate(1:find(crate=='_',1)-1);
    lbl = [strrep(strrep(pre, 'c0', '0.'), 'c', '') 'C'];
    xp = max(t_all) * 1.02;
    text(axChg(1), xp, median(V_all), lbl, 'FontSize', 8);
    text(axChg(2), xp, median(I_all), lbl, 'FontSize', 8);
    text(axChg(3), xp, median(P_all), lbl, 'FontSize', 8);
    idxT = find(strcmp({Charge_T_avg_list.crate_name}, crate));
    if ~isempty(idxT)
        A = Charge_T_avg_list(idxT(1));
        text(axChg(4), max(A.t_sec)*1.02, A.T_avg(end), lbl, 'FontSize', 8);
    end
end

% ----- Figure 2: Discharge 전용 (모든 사이클/채널 취합) -----
fig_discharge = figure('Name', 'VIT Discharge (all cycles combined)', 'NumberTitle', 'off');
tl_dch = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(tl_dch, 1); hold(ax1, 'on'); grid(ax1, 'on'); xlabel(ax1, 'Time [s]'); ylabel(ax1, 'V [V]'); title(ax1, 'V - Time');
ax2 = nexttile(tl_dch, 2); hold(ax2, 'on'); grid(ax2, 'on'); xlabel(ax2, 'Time [s]'); ylabel(ax2, 'I [A]'); title(ax2, 'I - Time');
ax3 = nexttile(tl_dch, 3); hold(ax3, 'on'); grid(ax3, 'on'); xlabel(ax3, 'Time [s]'); ylabel(ax3, 'Power [W]'); title(ax3, 'Power - Time');
ax4 = nexttile(tl_dch, 4); hold(ax4, 'on'); grid(ax4, 'on'); xlabel(ax4, 'Time [s]'); ylabel(ax4, 'T [°C]'); title(ax4, 'T - Time (8ch avg)');
cycle_shown_dch = {};

for k = 1:length(Discharge_data_list)
    D = Discharge_data_list(k);
    cname = D.cycle_name;
    [~, cidx] = ismember(cname, cyc_fields);
    cidx = max(1, cidx);
    col = cycle_colors(cidx, :);
    show_legend = ~ismember(cname, cycle_shown_dch);
    if show_legend, cycle_shown_dch{end+1} = cname; end %#ok<AGROW>
    t_sec = local_time_to_seconds(D.t_raw);
    hv = {'off','on'};
    scatter(ax1, t_sec, D.V_raw, msz, col, 'filled', 'DisplayName', cname, 'HandleVisibility', hv{show_legend+1});
    scatter(ax2, t_sec, D.I_raw, msz, col, 'filled', 'HandleVisibility', 'off');
    scatter(ax3, t_sec, D.P_raw, msz, col, 'filled', 'HandleVisibility', 'off');
end
cycle_shown_T = {};
for g = 1:length(Discharge_T_avg_list)
    A = Discharge_T_avg_list(g);
    [~, cidx] = ismember(A.cycle_name, cyc_fields);
    cidx = max(1, cidx);
    col = cycle_colors(cidx, :);
    show_leg = ~ismember(A.cycle_name, cycle_shown_T);
    if show_leg, cycle_shown_T{end+1} = A.cycle_name; end %#ok<AGROW>
    hv = {'off','on'};
    plot(ax4, A.t_sec, A.T_avg, '-', 'Color', col, 'LineWidth', 1.5, 'DisplayName', A.cycle_name, 'HandleVisibility', hv{show_leg+1});
end
legend(ax1, 'Location', 'best');
legend(ax4, 'Location', 'best');
title(tl_dch, 'Discharge (all cycles)');
axDch = [ax1, ax2, ax3, ax4];
uq_crates_dch = unique({Discharge_data_list.crate_name});
for cr = 1:numel(uq_crates_dch)
    crate = uq_crates_dch{cr};
    idx = find(strcmp({Discharge_data_list.crate_name}, crate));
    t_all = []; V_all = []; I_all = []; P_all = [];
    for i = idx
        D = Discharge_data_list(i);
        ts = local_time_to_seconds(D.t_raw);
        t_all = [t_all; ts(:)]; V_all = [V_all; D.V_raw(:)]; I_all = [I_all; D.I_raw(:)]; P_all = [P_all; D.P_raw(:)];
    end
    pre = crate(1:find(crate=='_',1)-1);
    lbl = [strrep(strrep(pre, 'c0', '0.'), 'c', '') 'C'];
    xp = max(t_all) * 1.02;
    text(axDch(1), xp, median(V_all), lbl, 'FontSize', 8);
    text(axDch(2), xp, median(I_all), lbl, 'FontSize', 8);
    text(axDch(3), xp, median(P_all), lbl, 'FontSize', 8);
    idxT = find(strcmp({Discharge_T_avg_list.crate_name}, crate));
    if ~isempty(idxT)
        A = Discharge_T_avg_list(idxT(1));
        text(axDch(4), max(A.t_sec)*1.02, A.T_avg(end), lbl, 'FontSize', 8);
    end
end

% 4개 서브플롯 x,y축 범위를 충전·방전 동일하게 설정
for tile = 1:4
    x1 = xlim(axChg(tile)); x2 = xlim(axDch(tile));
    y1 = ylim(axChg(tile)); y2 = ylim(axDch(tile));
    xcommon = [min(x1(1), x2(1)), max(x1(2), x2(2))];
    ycommon = [min(y1(1), y2(1)), max(y1(2), y2(2))];
    xlim(axChg(tile), xcommon); ylim(axChg(tile), ycommon);
    xlim(axDch(tile), xcommon); ylim(axDch(tile), ycommon);
end
saveas(fig_charge, fullfile(saveDir, 'RPT_Crate_Charge_AllCycles.fig'));
saveas(fig_discharge, fullfile(saveDir, 'RPT_Crate_Discharge_AllCycles.fig'));

fprintf('\n=== RPT_Crate_Cycle_Visualization complete (Field_Visualize_VIT_byYear format: 2x2 scatter). ===\n');

%% ========================================================================
% Local function: dQ/dV 계산 유틸리티
% ========================================================================
function [V_mid, dQdV] = local_compute_dQdV(V_grid, Q, apply_smoothing, smooth_window)
    % V_grid, Q: 같은 길이 벡터
    V_grid = V_grid(:);
    Q = Q(:);
    
    % NaN 제거
    valid = ~(isnan(V_grid) | isnan(Q));
    V_grid = V_grid(valid);
    Q = Q(valid);
    
    if numel(V_grid) < 3
        V_mid = [];
        dQdV = [];
        return;
    end
    
    % 전압이 증가하도록 정렬 (charge 기준, discharge도 정렬해서 사용)
    [V_sorted, idx] = sort(V_grid);
    Q_sorted = Q(idx);
    
    % smoothing (옵션)
    if apply_smoothing
        try
            Q_sorted = smoothdata(Q_sorted, 'movmean', smooth_window);
        catch
            % 구버전 MATLAB 호환: movmean 옵션 실패 시, simple moving average
            w = smooth_window;
            if mod(w, 2) == 0
                w = w + 1; % 홀수로 강제
            end
            Q_sorted = conv(Q_sorted, ones(w,1)/w, 'same');
        end
    end
    
    dQ = diff(Q_sorted);
    dV = diff(V_sorted);
    
    % 0이나 매우 작은 dV 방지
    dV(abs(dV) < 1e-6) = NaN;
    
    dQdV = dQ ./ dV;
    V_mid = (V_sorted(1:end-1) + V_sorted(2:end)) / 2;
end

%% 전류 유도: I [A] = 3600 * dQ/dt (Q: Ah, t: s). t 감소(방전)도 처리: t 기준 정렬 후 유도하고 원래 순서로 복원
function I = local_derive_current(Q, t)
    Q = Q(:);
    t = t(:);
    if numel(Q) ~= numel(t) || numel(t) < 2
        I = [];
        return;
    end
    % t가 감소하는 경우(방전 그리드)에도 dQ/dt가 나오도록 t 오름차순 정렬 후 계산
    [t_sorted, idx_sort] = sort(t);
    Q_sorted = Q(idx_sort);
    dt = gradient(t_sorted);
    dt(dt <= 0) = NaN;
    dQ = gradient(Q_sorted);
    I_sorted = 3600 * (dQ ./ dt);
    I_sorted(~isfinite(I_sorted)) = NaN;
    % 원래 t 순서로 복원 (plot(t, I) 쌍 맞추기)
    I = zeros(size(t));
    I(idx_sort) = I_sorted;
end

%% 온도 시계열 추출 (시각화용: t_T, T 쌍 그대로 반환. T1_raw는 t_raw에 맞춰 보간하지 않음)
function [t_T, T_out] = local_get_temp_series(seg_struct, t_vec)
    t_vec = t_vec(:);
    n = numel(t_vec);
    t_T = [];
    T_out = [];

    % T1_raw + t_raw 있으면 그대로 (t_raw, T1_raw) 반환
    if isfield(seg_struct, 'T1_raw') && isfield(seg_struct, 't_raw')
        T1r = seg_struct.T1_raw(:);
        tr = seg_struct.t_raw(:);
        if numel(tr) == numel(T1r) && numel(tr) >= 1
            t_T = tr;
            T_out = T1r;
            return;
        end
    end
    % T1 (길이 맞으면 t_vec과 쌍)
    if isfield(seg_struct, 'T1') && numel(seg_struct.T1(:)) == n
        t_T = t_vec;
        T_out = seg_struct.T1(:);
        return;
    end
    % T1_raw만 있고 길이 맞으면 t_vec과 쌍
    if isfield(seg_struct, 'T1_raw')
        T1r = seg_struct.T1_raw(:);
        if numel(T1r) == n
            t_T = t_vec;
            T_out = T1r;
            return;
        end
    end
    % Temp / Temperature (길이 맞으면 t_vec과 쌍)
    if isfield(seg_struct, 'Temp') && numel(seg_struct.Temp(:)) == n
        t_T = t_vec;
        T_out = seg_struct.Temp(:);
        return;
    end
    if isfield(seg_struct, 'Temperature') && numel(seg_struct.Temperature(:)) == n
        t_T = t_vec;
        T_out = seg_struct.Temperature(:);
        return;
    end
end

%% (사이클 × C-rate)별 8채널 평균 T 시계열. 공통 시간 = 8채널 모두 유효한 구간(교집합) → 끝 급상승 방지
function out_list = local_avg_temp_by_cycle_crate(data_list)
    out_list = [];
    if isempty(data_list), return; end
    % (cycle_name, crate_name) 쌍으로 그룹화
    keys = arrayfun(@(k) sprintf('%s|%s', data_list(k).cycle_name, data_list(k).crate_name), 1:numel(data_list), 'UniformOutput', false);
    [uq_keys, ~, ic] = unique(keys);
    n_pts = 300;   % 공통 시간 샘플 수
    for g = 1:length(uq_keys)
        idx = find(ic == g);
        t_all = [];
        T_all = [];
        % 1번: 공통 시간 = 8채널 모두 유효한 구간(교집합) → 끝구간 급상승 방지
        t_min_per_ch = [];
        t_max_per_ch = [];
        for i = 1:numel(idx)
            D = data_list(idx(i));
            t_sec = local_time_to_seconds( D.t_raw_temp );
            if isempty(t_sec) || isempty(D.Temp_raw) || numel(t_sec) ~= numel(D.Temp_raw), continue; end
            t_min_per_ch(end+1) = min(t_sec(:)); %#ok<AGROW>
            t_max_per_ch(end+1) = max(t_sec(:)); %#ok<AGROW>
        end
        if isempty(t_min_per_ch), continue; end
        t_min_common = max(t_min_per_ch);
        t_max_common = min(t_max_per_ch);
        if t_max_common <= t_min_common
            t_common = t_min_common;
            T_interp = NaN;
        else
            t_common = linspace(t_min_common, t_max_common, n_pts)';
            T_mat = NaN(n_pts, numel(idx));
            cnt = 0;
            for i = 1:numel(idx)
                D = data_list(idx(i));
                t_sec = local_time_to_seconds( D.t_raw_temp );
                if isempty(D.Temp_raw) || numel(t_sec) ~= numel(D.Temp_raw), continue; end
                cnt = cnt + 1;
                T_mat(:, cnt) = interp1(t_sec(:), D.Temp_raw(:), t_common, 'linear');
            end
            T_interp = mean(T_mat(:, 1:cnt), 2);
        end
        part = strsplit(uq_keys{g}, '|');
        out_list(g).cycle_name = part{1}; %#ok<AGROW>
        out_list(g).crate_name = part{2}; %#ok<AGROW>
        out_list(g).t_sec = t_common; %#ok<AGROW>
        out_list(g).T_avg = T_interp; %#ok<AGROW>
    end
end

%% 시간을 [s] double로 통일. 항상 "세그먼트 시작=0, 경과시간 증가" (왼쪽=시작, 오른쪽=끝)
%  방전처럼 t가 감소 저장된 경우: t_sec = max(t)-t 로 해서 시작=0, 끝=max-min
function t_sec = local_time_to_seconds(t_vec)
    if isempty(t_vec)
        t_sec = [];
        return;
    end
    t_vec = t_vec(:);
    if isduration(t_vec)
        t_dbl = double(seconds(t_vec));
    else
        t_dbl = double(t_vec);
    end
    t_min = min(t_dbl);
    t_max = max(t_dbl);
    % t가 감소 순이면 첫 점이 세그먼트 시작 → 경과시간 = max(t)-t 로 0부터 증가
    if t_dbl(1) > t_dbl(end)
        t_sec = t_max - t_dbl;
    else
        t_sec = t_dbl - t_min;
    end
    t_sec = double(t_sec);
end
