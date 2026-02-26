%% Visualize_VIT_byYear.m - 연도별 V, I, T, Power 시각화
% FieldQmax_dQdV.m과 동일한 데이터 소스 및 세그먼트 선택 사용
% 생성: 연도별 V-time, I-time, Power-time, T-time (충전/방전), 연도별 V-I
% 저장: Field_VIT 폴더에 .fig

clear; clc; close all;

%% 데이터 소스 (FieldQmax_dQdV.m과 동일)
% 2021: Rack_raw2mat\Old\2021\202106\Raw_20210603.mat
% 2023: Rack_raw2mat\New\2023\202310\Raw_20231016.mat
% 2024: Rack_raw2mat\New\2024\202409\Raw_20240909.mat
% 2025: Rack_raw2mat\New\2025\202507\Raw_20250711.mat
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
dataFile_2021 = fullfile(baseDir, 'Old', '2021', '202106', 'Raw_20210603.mat');
dataFile_2023 = fullfile(baseDir, 'New', '2023', '202310', 'Raw_20231016.mat');
dataFile_2024 = fullfile(baseDir, 'New', '2024', '202409', 'Raw_20240909.mat');
dataFile_2025 = fullfile(baseDir, 'New', '2025', '202507', 'Raw_20250711.mat');

dates = {'2021-06-03', '2023-10-16', '2024-09-09', '2025-07-11'};
dataFiles = {dataFile_2021, dataFile_2023, dataFile_2024, dataFile_2025};
dataTypes = {'old', 'new', 'new', 'new'};
base_dates = {datetime(2021,6,3), datetime(2023,10,16), datetime(2024,9,9), datetime(2025,7,11)};
min_charge_secs = [600, 300, 300, 480];
min_discharge_secs = [300, 150, 300, 300];

Np = 2;
dt = 1;
C_cell_Ah = 64;
thr_A = C_cell_Ah * 0.05;
thr_cell = thr_A / Np;

%% 출력 폴더 (스크립트와 동일한 Field_VIT)
scriptDir = fileparts(mfilename('fullpath'));
saveDir = scriptDir;
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

%% 세그먼트 검출 헬퍼 (FieldQmax_dQdV.m과 동일 로직)
find_segments = @(mask) local_find_segments(mask);

%% 연도별 충전/방전 세그먼트 데이터 수집 (V, I, P, T, t_normalized)
Charge_VIT = struct();
Discharge_VIT = struct();
yearLabels = {'Y2021', 'Y2023', 'Y2024', 'Y2025'};

for date_idx = 1:length(dates)
    date_str = dates{date_idx};
    dataFile = dataFiles{date_idx};
    dataType = dataTypes{date_idx};
    base_date = base_dates{date_idx};
    min_charge_sec = min_charge_secs(date_idx);
    min_discharge_sec = min_discharge_secs(date_idx);
    yearKey = yearLabels{date_idx};

    if ~exist(dataFile, 'file')
        warning('File not found: %s', dataFile);
        continue;
    end

    S = load(dataFile);
    if strcmp(dataType, 'old')
        if ~isfield(S, 'Raw') || ~isfield(S.Raw, 'Rack01')
            warning('Old format: Raw.Rack01 not found in %s', dataFile);
            continue;
        end
        D = S.Raw.Rack01;
        if isfield(D, 'Time')
            t = datetime(D.Time);
        elseif isfield(D, 'Date_Time')
            if isduration(D.Date_Time), t = base_date + D.Date_Time;
            else, t = datetime(D.Date_Time); end
        else
            warning('No Time/Date_Time in %s', dataFile); continue;
        end
        I_rack = D.DCCurrent_A(:);
        Vcell_avg = D.AverageCV_V(:);
        if isfield(D, 'DCPower_kW'), P_rack_kW = D.DCPower_kW(:);
        else, P_rack_kW = zeros(size(I_rack)); end
        if isfield(D, 'AverageMT_degC'), T_degC = D.AverageMT_degC(:);
        else, T_degC = nan(size(I_rack)); end
    else
        D = S.Raw;
        if isduration(D.Date_Time), t = base_date + D.Date_Time;
        else, t = datetime(D.Date_Time); end
        I_rack = D.DCCurrent(:);
        Vcell_avg = D.CVavg(:);
        if isfield(D, 'DCPower')
            pw = D.DCPower(:);
            if max(abs(pw)) < 1e4, P_rack_kW = pw; else, P_rack_kW = pw / 1000; end
        else
            P_rack_kW = zeros(size(I_rack));
        end
        if isfield(D, 'MTavg'), T_degC = D.MTavg(:);
        else, T_degC = nan(size(I_rack)); end
    end

    t0 = t(1);
    tsec = seconds(t - t0);
    I_cell = I_rack / Np;

    isIdle = abs(I_cell) < thr_cell;
    isChg = I_cell > thr_cell;
    isDischg = I_cell < -thr_cell;
    chgSegs = find_segments(isChg);
    dischgSegs = find_segments(isDischg);
    dur = @(seg) seg(:,2) - seg(:,1) + 1;
    chgSegs = chgSegs(dur(chgSegs) >= ceil(min_charge_sec/dt), :);
    dischgSegs = dischgSegs(dur(dischgSegs) >= ceil(min_discharge_sec/dt), :);

    % 2024: Chg02, Dchg03만 사용 (FieldQmax_dQdV.m과 동일)
    if strcmp(date_str, '2024-09-09')
        target_chg_start = datetime(2024, 9, 9, 12, 55, 0);
        target_chg_end = datetime(2024, 9, 9, 13, 13, 0);
        chg_seg_idx = [];
        for k = 1:size(chgSegs,1)
            if abs(t(chgSegs(k,1)) - target_chg_start) < minutes(5) && abs(t(chgSegs(k,2)) - target_chg_end) < minutes(5)
                chg_seg_idx = k; break;
            end
        end
        if isempty(chg_seg_idx), chgSegs = []; else, chgSegs = chgSegs(chg_seg_idx,:); end

        target_dischg_start = datetime(2024, 9, 9, 14, 14, 0);
        target_rest_start = datetime(2024, 9, 9, 14, 25, 49);
        dischg_seg_idx = [];
        for k = 1:size(dischgSegs,1)
            if abs(t(dischgSegs(k,1)) - target_dischg_start) < minutes(1)
                dischg_seg_idx = k; break;
            end
        end
        if ~isempty(dischg_seg_idx)
            dischg_end = dischgSegs(dischg_seg_idx,2);
            if t(dischg_end) > target_rest_start
                idx_limit = find(t <= target_rest_start, 1, 'last');
                if ~isempty(idx_limit) && idx_limit >= dischgSegs(dischg_seg_idx,1)
                    dischgSegs(dischg_seg_idx,2) = idx_limit;
                end
            end
            dischgSegs = dischgSegs(dischg_seg_idx,:);
        else
            dischgSegs = [];
        end
    end

    Charge_VIT.(yearKey) = [];
    for k = 1:size(chgSegs,1)
        s = chgSegs(k,1); e = chgSegs(k,2);
        t_norm = tsec(s:e) - tsec(s);
        Charge_VIT.(yearKey)(end+1).t = t_norm;
        Charge_VIT.(yearKey)(end).V = Vcell_avg(s:e);
        Charge_VIT.(yearKey)(end).I = I_cell(s:e);
        % 셀 단위 전력 [kW]: P_cell = V_cell * I_cell (W) / 1000
        Charge_VIT.(yearKey)(end).P_kW = (Vcell_avg(s:e) .* I_cell(s:e)) / 1000;
        Charge_VIT.(yearKey)(end).T = T_degC(s:e);
        Charge_VIT.(yearKey)(end).date = date_str;
    end

    Discharge_VIT.(yearKey) = [];
    for k = 1:size(dischgSegs,1)
        s = dischgSegs(k,1); e = dischgSegs(k,2);
        t_norm = tsec(s:e) - tsec(s);
        Discharge_VIT.(yearKey)(end+1).t = t_norm;
        Discharge_VIT.(yearKey)(end).V = Vcell_avg(s:e);
        Discharge_VIT.(yearKey)(end).I = I_cell(s:e);
        % 셀 단위 전력 [kW]: P_cell = V_cell * I_cell (W) / 1000
        Discharge_VIT.(yearKey)(end).P_kW = (Vcell_avg(s:e) .* I_cell(s:e)) / 1000;
        Discharge_VIT.(yearKey)(end).T = T_degC(s:e);
        Discharge_VIT.(yearKey)(end).date = date_str;
    end
end

% Scatter marker size (연도별 scatter용)
msz = 10;

%% 1) 연도별 V-time, I-time, Power-time, T-time (충전 + 방전 각각) — scatter
colors = lines(length(yearLabels));
for y = 1:length(yearLabels)
    yk = yearLabels{y};
    fig_chg = figure('Name', ['VIT Charge ', yk], 'NumberTitle', 'off');
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax1 = nexttile(tl, 1); hold(ax1, 'on'); grid(ax1, 'on'); xlabel(ax1, 'Time [s]'); ylabel(ax1, 'V [V]'); title(ax1, 'V - Time');
    ax2 = nexttile(tl, 2); hold(ax2, 'on'); grid(ax2, 'on'); xlabel(ax2, 'Time [s]'); ylabel(ax2, 'I [A]'); title(ax2, 'I - Time');
    ax3 = nexttile(tl, 3); hold(ax3, 'on'); grid(ax3, 'on'); xlabel(ax3, 'Time [s]'); ylabel(ax3, 'Power [W]'); title(ax3, 'Power - Time');
    ax4 = nexttile(tl, 4); hold(ax4, 'on'); grid(ax4, 'on'); xlabel(ax4, 'Time [s]'); ylabel(ax4, 'T [°C]'); title(ax4, 'T - Time');
    if ~isempty(Charge_VIT.(yk))
        all_t = []; all_V = []; all_I = []; all_P = []; all_T = [];
        for seg = 1:length(Charge_VIT.(yk))
            d = Charge_VIT.(yk)(seg);
            scatter(ax1, d.t, d.V, msz, colors(y,:), 'filled');
            scatter(ax2, d.t, d.I, msz, colors(y,:), 'filled');
            scatter(ax3, d.t, d.P_kW*1000, msz, colors(y,:), 'filled');
            scatter(ax4, d.t, d.T, msz, colors(y,:), 'filled');
            all_t = [all_t; d.t(:)];
            all_V = [all_V; d.V(:)];
            all_I = [all_I; d.I(:)];
            all_P = [all_P; d.P_kW(:)];
            all_T = [all_T; d.T(:)];
        end
        max_t = max(all_t);
        max_V = max(all_V, [], 'omitnan');
        max_I = max(all_I, [], 'omitnan');
        max_P = max(all_P, [], 'omitnan');
        max_T = max(all_T, [], 'omitnan');

        if ~(isnan(max_I) || isnan(max_P))
            C_rate = max_I / C_cell_Ah;
            P_nom_W = 3.7 * C_cell_Ah;
            P_rate = (max_P * 1000) / P_nom_W;
            str_CR = sprintf('~%.2f C, ~%.2f P', C_rate, P_rate);
            text(ax1, max_t, max_V, str_CR, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
            text(ax2, max_t, max_I, str_CR, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
            text(ax3, max_t, max_P*1000, str_CR, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
            text(ax4, max_t, max_T, str_CR, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
        end
    end
    title(tl, [yk ' Charge']);
    saveas(fig_chg, fullfile(saveDir, ['VIT_Charge_', yk, '.fig']));
end

for y = 1:length(yearLabels)
    yk = yearLabels{y};
    fig_dchg = figure('Name', ['VIT Discharge ', yk], 'NumberTitle', 'off');
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax1 = nexttile(tl, 1); hold(ax1, 'on'); grid(ax1, 'on'); xlabel(ax1, 'Time [s]'); ylabel(ax1, 'V [V]'); title(ax1, 'V - Time');
    ax2 = nexttile(tl, 2); hold(ax2, 'on'); grid(ax2, 'on'); xlabel(ax2, 'Time [s]'); ylabel(ax2, 'I [A]'); title(ax2, 'I - Time');
    ax3 = nexttile(tl, 3); hold(ax3, 'on'); grid(ax3, 'on'); xlabel(ax3, 'Time [s]'); ylabel(ax3, 'Power [W]'); title(ax3, 'Power - Time');
    ax4 = nexttile(tl, 4); hold(ax4, 'on'); grid(ax4, 'on'); xlabel(ax4, 'Time [s]'); ylabel(ax4, 'T [°C]'); title(ax4, 'T - Time');
    if ~isempty(Discharge_VIT.(yk))
        all_t = []; all_V = []; all_I = []; all_P = []; all_T = [];
        for seg = 1:length(Discharge_VIT.(yk))
            d = Discharge_VIT.(yk)(seg);
            scatter(ax1, d.t, d.V, msz, colors(y,:), 'filled');
            scatter(ax2, d.t, d.I, msz, colors(y,:), 'filled');
            scatter(ax3, d.t, d.P_kW*1000, msz, colors(y,:), 'filled');
            scatter(ax4, d.t, d.T, msz, colors(y,:), 'filled');
            all_t = [all_t; d.t(:)];
            all_V = [all_V; d.V(:)];
            all_I = [all_I; d.I(:)];
            all_P = [all_P; d.P_kW(:)];
            all_T = [all_T; d.T(:)];
        end
        max_t = max(all_t);
        max_V = max(all_V, [], 'omitnan');
        max_I_abs = max(abs(all_I), [], 'omitnan');
        max_P_abs = max(abs(all_P), [], 'omitnan');
        max_T = max(all_T, [], 'omitnan');

        if ~(isnan(max_I_abs) || isnan(max_P_abs))
            C_rate = max_I_abs / C_cell_Ah;
            P_nom_W = 3.7 * C_cell_Ah;
            P_rate = (max_P_abs * 1000) / P_nom_W;
            str_CR = sprintf('~%.2f C, ~%.2f P', C_rate, P_rate);
            text(ax1, max_t, max_V, str_CR, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
            text(ax2, max_t, max(all_I, [], 'omitnan'), str_CR, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
            text(ax3, max_t, max(all_P, [], 'omitnan')*1000, str_CR, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
            text(ax4, max_t, max_T, str_CR, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
        end
    end
    title(tl, [yk ' Discharge']);
    saveas(fig_dchg, fullfile(saveDir, ['VIT_Discharge_', yk, '.fig']));
end

%% 2) 모든 연도 합쳐서 — 선으로 연결 (충전 / 방전 각각)
fig_all_chg = figure('Name', 'VIT Charge (all years combined)', 'NumberTitle', 'off');
tl_all_chg = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(tl_all_chg, 1); hold(ax1, 'on'); grid(ax1, 'on'); xlabel(ax1, 'Time [s]'); ylabel(ax1, 'V [V]'); title(ax1, 'V - Time');
ax2 = nexttile(tl_all_chg, 2); hold(ax2, 'on'); grid(ax2, 'on'); xlabel(ax2, 'Time [s]'); ylabel(ax2, 'I [A]'); title(ax2, 'I - Time');
ax3 = nexttile(tl_all_chg, 3); hold(ax3, 'on'); grid(ax3, 'on'); xlabel(ax3, 'Time [s]'); ylabel(ax3, 'Power [W]'); title(ax3, 'Power - Time');
ax4 = nexttile(tl_all_chg, 4); hold(ax4, 'on'); grid(ax4, 'on'); xlabel(ax4, 'Time [s]'); ylabel(ax4, 'T [°C]'); title(ax4, 'T - Time');
for y = 1:length(yearLabels)
    yk = yearLabels{y};
    for seg = 1:length(Charge_VIT.(yk))
        d = Charge_VIT.(yk)(seg);
        if seg == 1
            scatter(ax1, d.t, d.V, msz, colors(y,:), 'filled', 'DisplayName', yk);
        else
            scatter(ax1, d.t, d.V, msz, colors(y,:), 'filled', 'HandleVisibility', 'off');
        end
        scatter(ax2, d.t, d.I, msz, colors(y,:), 'filled', 'HandleVisibility', 'off');
        scatter(ax3, d.t, d.P_kW*1000, msz, colors(y,:), 'filled', 'HandleVisibility', 'off');
        scatter(ax4, d.t, d.T, msz, colors(y,:), 'filled', 'HandleVisibility', 'off');
    end
end
% 각 연도 데이터 끝에 C, P 텍스트 (합쳐진 그림)
P_nom_W = 3.7 * C_cell_Ah;
for y = 1:length(yearLabels)
    yk = yearLabels{y};
    if isempty(Charge_VIT.(yk)), continue; end
    all_t = []; all_V = []; all_I = []; all_P = []; all_T = [];
    for seg = 1:length(Charge_VIT.(yk))
        d = Charge_VIT.(yk)(seg);
        all_t = [all_t; d.t(:)]; all_V = [all_V; d.V(:)]; all_I = [all_I; d.I(:)]; all_P = [all_P; d.P_kW(:)]; all_T = [all_T; d.T(:)];
    end
    max_t = max(all_t);
    max_V = max(all_V, [], 'omitnan'); max_I = max(all_I, [], 'omitnan');
    max_P = max(all_P, [], 'omitnan'); max_T = max(all_T, [], 'omitnan');
    if isnan(max_I) || isnan(max_P), continue; end
    C_rate = max_I / C_cell_Ah;
    P_rate = (max_P * 1000) / P_nom_W;
    str_CR = sprintf('~%.2f C, ~%.2f P', C_rate, P_rate);
    text(ax1, max_t, max_V, str_CR, 'Color', colors(y,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
    text(ax2, max_t, max_I, str_CR, 'Color', colors(y,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
    text(ax3, max_t, max_P*1000, str_CR, 'Color', colors(y,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
    text(ax4, max_t, max_T, str_CR, 'Color', colors(y,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
end
legend(ax1, 'Location', 'best');
title(tl_all_chg, 'Charge (all years)');
saveas(fig_all_chg, fullfile(saveDir, 'VIT_Charge_allYears_combined.fig'));

fig_all_dchg = figure('Name', 'VIT Discharge (all years combined)', 'NumberTitle', 'off');
tl_all_dchg = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(tl_all_dchg, 1); hold(ax1, 'on'); grid(ax1, 'on'); xlabel(ax1, 'Time [s]'); ylabel(ax1, 'V [V]'); title(ax1, 'V - Time');
ax2 = nexttile(tl_all_dchg, 2); hold(ax2, 'on'); grid(ax2, 'on'); xlabel(ax2, 'Time [s]'); ylabel(ax2, 'I [A]'); title(ax2, 'I - Time');
ax3 = nexttile(tl_all_dchg, 3); hold(ax3, 'on'); grid(ax3, 'on'); xlabel(ax3, 'Time [s]'); ylabel(ax3, 'Power [W]'); title(ax3, 'Power - Time');
ax4 = nexttile(tl_all_dchg, 4); hold(ax4, 'on'); grid(ax4, 'on'); xlabel(ax4, 'Time [s]'); ylabel(ax4, 'T [°C]'); title(ax4, 'T - Time');
for y = 1:length(yearLabels)
    yk = yearLabels{y};
    for seg = 1:length(Discharge_VIT.(yk))
        d = Discharge_VIT.(yk)(seg);
        if seg == 1
            scatter(ax1, d.t, d.V, msz, colors(y,:), 'filled', 'DisplayName', yk);
        else
            scatter(ax1, d.t, d.V, msz, colors(y,:), 'filled', 'HandleVisibility', 'off');
        end
        scatter(ax2, d.t, d.I, msz, colors(y,:), 'filled', 'HandleVisibility', 'off');
        scatter(ax3, d.t, d.P_kW*1000, msz, colors(y,:), 'filled', 'HandleVisibility', 'off');
        scatter(ax4, d.t, d.T, msz, colors(y,:), 'filled', 'HandleVisibility', 'off');
    end
end
% 각 연도 데이터 끝에 C, P 텍스트 (합쳐진 그림)
P_nom_W = 3.7 * C_cell_Ah;
for y = 1:length(yearLabels)
    yk = yearLabels{y};
    if isempty(Discharge_VIT.(yk)), continue; end
    all_t = []; all_V = []; all_I = []; all_P = []; all_T = [];
    for seg = 1:length(Discharge_VIT.(yk))
        d = Discharge_VIT.(yk)(seg);
        all_t = [all_t; d.t(:)]; all_V = [all_V; d.V(:)]; all_I = [all_I; d.I(:)]; all_P = [all_P; d.P_kW(:)]; all_T = [all_T; d.T(:)];
    end
    max_t = max(all_t);
    max_V = max(all_V, [], 'omitnan'); max_I_abs = max(abs(all_I), [], 'omitnan');
    max_P_abs = max(abs(all_P), [], 'omitnan'); max_T = max(all_T, [], 'omitnan');
    if isnan(max_I_abs) || isnan(max_P_abs), continue; end
    C_rate = max_I_abs / C_cell_Ah;
    P_rate = (max_P_abs * 1000) / P_nom_W;
    str_CR = sprintf('~%.2f C, ~%.2f P', C_rate, P_rate);
    text(ax1, max_t, max_V, str_CR, 'Color', colors(y,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
    text(ax2, max_t, max(all_I, [], 'omitnan'), str_CR, 'Color', colors(y,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
    text(ax3, max_t, max(all_P, [], 'omitnan')*1000, str_CR, 'Color', colors(y,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
    text(ax4, max_t, max_T, str_CR, 'Color', colors(y,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
end
legend(ax1, 'Location', 'best');
title(tl_all_dchg, 'Discharge (all years)');
saveas(fig_all_dchg, fullfile(saveDir, 'VIT_Discharge_allYears_combined.fig'));

%% 3) 연도별 V-I 시각화 (충전 / 방전 각각) — scatter
fig_VI_chg = figure('Name', 'V-I Charge (all years)', 'NumberTitle', 'off');
ax_chg = axes(fig_VI_chg); hold(ax_chg, 'on'); grid(ax_chg, 'on');
fig_VI_dchg = figure('Name', 'V-I Discharge (all years)', 'NumberTitle', 'off');
ax_dchg = axes(fig_VI_dchg); hold(ax_dchg, 'on'); grid(ax_dchg, 'on');

for y = 1:length(yearLabels)
    yk = yearLabels{y};
    for seg = 1:length(Charge_VIT.(yk))
        d = Charge_VIT.(yk)(seg);
        scatter(ax_chg, d.I, d.V, msz, colors(y,:), 'filled', 'DisplayName', [yk ' seg' num2str(seg)]);
    end
    for seg = 1:length(Discharge_VIT.(yk))
        d = Discharge_VIT.(yk)(seg);
        scatter(ax_dchg, d.I, d.V, msz, colors(y,:), 'filled', 'DisplayName', [yk ' seg' num2str(seg)]);
    end
end
xlabel(ax_chg, 'Current [A]'); ylabel(ax_chg, 'Voltage [V]'); title(ax_chg, 'V-I (Charge, by year)'); legend(ax_chg, 'Location', 'best');
xlabel(ax_dchg, 'Current [A]'); ylabel(ax_dchg, 'Voltage [V]'); title(ax_dchg, 'V-I (Discharge, by year)'); legend(ax_dchg, 'Location', 'best');

saveas(fig_VI_chg, fullfile(saveDir, 'VIT_VI_Charge_allYears.fig'));
saveas(fig_VI_dchg, fullfile(saveDir, 'VIT_VI_Discharge_allYears.fig'));
fprintf('VIT figures saved to: %s\n', saveDir);

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
