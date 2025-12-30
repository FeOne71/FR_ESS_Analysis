%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldQmax_VoltageAnalysis.m
%
% Rest 전압 강하 분석 및 OCV 안정화 요약
% - Figure 1: Rest Voltage Drop Analysis (2RC 피팅)
% - Figure 2: OCV Stabilization Summary
%
% 이 스크립트는 FieldQmax_Olddata.m에서 분리된 전압 안정화 분석 코드입니다.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% Parameters
Cnom = 128;                         % Rack nominal Capacity (Ah) (참조용)
C_cell_Ah = 64;                     % Cell capacity (Ah)
thr_A = C_cell_Ah * 0.05;                % Idle threshold (A)  Cnom*0.05가 ppt 분석 내용
min_charge_sec = 600;               % 연속 10분 이상 충전 구간
dt = 1;                             % s (고정 가정)

% config (참조: 17 modules * 14S, 2P)
Ns = 17 * 14;                       % 238S (참조)
Np = 2;                             % 2P (참조)

%% Paths
dataFile = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old\2021\202106\Raw_20210607.mat';
ocvFile  = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\OCV_integrated.mat';

% Create FieldQmax folder for results
if ~exist('FieldQmax', 'dir')
    mkdir('FieldQmax');
end

%% Load rack data
S = load(dataFile);                 % expects Raw_Rack
if ~isfield(S, 'Raw')
    error('Raw_Rack not found in %s', dataFile);
end
Raw_Rack = S.Raw;
if ~isfield(Raw_Rack, 'Rack01')
    error('Raw_Rack.Rack01 not found in %s', dataFile);
end
D = Raw_Rack.Rack01;

% Time vector
if isfield(D, 'Time')
    t = datetime(D.Time);
elseif isfield(D, 'Date_Time')
    if isduration(D.Date_Time)
        t = datetime(2021,6,7) + D.Date_Time;
    else
        t = datetime(D.Date_Time);
    end
else
    error('No Time/Date_Time field present.');
end
t0 = t(1);
tsec = seconds(t - t0);

% Signals (rack variables per FieldData_RPT_Old.m)
I_rack    = D.DCCurrent_A(:);       % A (rack)
Vcell_avg = D.AverageCV_V(:);       % V (cell average)
SOC_raw   = D.SOCPct(:);            % Raw SOC from BMS (%)
SOH_raw   = D.SOHPct(:);           % Raw SOH from BMS (%)
SOH_end   = SOH_raw(end);           % Get the last SOH value
% Convert to cell units
I_cell    = I_rack / Np;            % A (cell)

% Masks (cell units)
thr_cell = thr_A / Np;
isIdle = abs(I_cell) < thr_cell;
isChg  = I_cell > thr_cell;

% Helper: find contiguous true segments -> Nx2 [start_idx, end_idx]
find_segments = @(mask) local_find_segments(mask);

idleSegs = find_segments(isIdle);
chgSegs  = find_segments(isChg);

% Filter charge segments by duration >= 600 s (dt=1s)
dur = @(seg) seg(:,2) - seg(:,1) + 1;   % samples
min_samples = ceil(min_charge_sec / dt);
chgSegs = chgSegs(dur(chgSegs) >= min_samples, :);

% Load OCV data and build inverse function (OCV -> SOC)
T = load(ocvFile);
OCV_data = T.OCV_data;

% Build inverse OCV->SOC using avg_ocv_rpt0 and soc_grid
ocv = OCV_data.avg_ocv_rpt0(:);
soc = OCV_data.soc_grid(:);
% ensure monotonic pairs
[ocv_sorted, uq] = unique(ocv, 'stable');
soc_sorted = soc(uq);
% Ensure SOC_from_OCV returns percent [0..100]
if max(soc_sorted) <= 1.5
    soc_sorted = soc_sorted * 100; % fraction -> percent
end
SOC_from_OCV = @(v) interp1(ocv_sorted, soc_sorted, v, 'linear', 'extrap');

%% For each charge segment, compute anchors and Qmax (rack)
Results = struct();
rows = {};
for k = 1:size(chgSegs,1)
    chg_start = chgSegs(k,1);
    chg_end   = chgSegs(k,2);

    % Find idle segments around charge
    prevIdleIdx = find(idleSegs(:,2) < chg_start, 1, 'last');
    nextIdleIdx = find(idleSegs(:,1) > chg_end, 1, 'first');

    if isempty(prevIdleIdx) || isempty(nextIdleIdx)
        % Skip if either anchor missing (no fallback policy)
        continue;
    end

    % Define time points
    befChg = idleSegs(prevIdleIdx,2);  % 충전 전 휴지구간 마지막 시점
    aftChg = idleSegs(nextIdleIdx,2);  % 충전 후 휴지구간 마지막 시점
    
    % SOC at idle end points (for SOC-OCV data)
    V_befChg = Vcell_avg(befChg);  % 충전 전 휴지 종료 시점 전압
    V_aftChg = Vcell_avg(aftChg);  % 충전 후 휴지 종료 시점 전압
    I_befChg = I_cell(befChg);     % 충전 전 휴지 종료 시점 전류
    I_aftChg = I_cell(aftChg);     % 충전 후 휴지 종료 시점 전류

    % SOC at idle end points (OCV-based and BMS)
    SOC1 = SOC_from_OCV(V_befChg);   % 충전 전 휴지 종료 시점 SOC (OCV-based)
    SOC2 = SOC_from_OCV(V_aftChg);   % 충전 후 휴지 종료 시점 SOC (OCV-based)
    SOC1_raw = SOC_raw(befChg);      % Raw SOC at befChg
    SOC2_raw = SOC_raw(aftChg);      % Raw SOC at aftChg

    % integrate cell current using trapz (convert time to seconds)
    t_sec = seconds(t - t(1));  % convert to seconds from start
    Q_Ah_cell = abs(trapz(t_sec(chg_start:chg_end), I_cell(chg_start:chg_end)) / 3600);   % charge-only (cell)
    Ah_total = trapz(t_sec(befChg:aftChg), I_cell(befChg:aftChg)) / 3600; % from befChg to aftChg (cell)
    Ah_to_chg_start = trapz(t_sec(befChg:chg_start), I_cell(befChg:chg_start)) / 3600; % from befChg to chg_start (cell)
    Ah_to_chg_end = trapz(t_sec(befChg:chg_end), I_cell(befChg:chg_end)) / 3600; % from befChg to chg_end (cell)

    % SOC change used for Qmax (from befChg to aftChg - using only rest periods)
    dSOC_q = SOC2 - SOC1;   % percent
    dSOC_q_BMS = SOC2_raw - SOC1_raw;   % percent (BMS)

    % Qmax uses dSOC between rest periods (befChg to aftChg)
    if ~isnan(dSOC_q) && dSOC_q ~= 0
        Qmax_cell_Ah = Q_Ah_cell / (abs(dSOC_q)/100);
    else
        Qmax_cell_Ah = NaN;
    end
    
    % Qmax using BMS SOC
    if ~isnan(dSOC_q_BMS) && dSOC_q_BMS ~= 0
        Qmax_cell_Ah_BMS = Q_Ah_cell / (abs(dSOC_q_BMS)/100);
    else
        Qmax_cell_Ah_BMS = NaN;
    end
    
    % Calculate rest period durations
    % First rest period (before charge)
    rest1_start_idx = idleSegs(prevIdleIdx, 1);
    rest1_duration_sec = befChg - rest1_start_idx + 1;  % seconds
    rest1_duration = duration(0, 0, rest1_duration_sec);
    
    % Second rest period (after charge)
    rest2_start_idx = idleSegs(nextIdleIdx, 1);
    rest2_duration_sec = aftChg - rest2_start_idx + 1;  % seconds
    rest2_duration = duration(0, 0, rest2_duration_sec);

    % store
    Results(k).idx = k; %#ok<SAGROW>
    Results(k).chg_start = chg_start;
    Results(k).chg_end   = chg_end;
    Results(k).befChg = befChg;  % 충전 전 휴지구간 마지막 시점
    Results(k).aftChg = aftChg;  % 충전 후 휴지구간 마지막 시점
    Results(k).V1 = Vcell_avg(befChg); Results(k).V2 = Vcell_avg(aftChg);
    Results(k).I1 = I_cell(befChg); Results(k).I2 = I_cell(aftChg);
    Results(k).SOC1 = SOC1; Results(k).SOC2 = SOC2; Results(k).dSOC = abs(dSOC_q);
    Results(k).SOC1_raw = SOC1_raw; Results(k).SOC2_raw = SOC2_raw;
    Results(k).Q_Ah_cell = Q_Ah_cell; Results(k).Qmax_cell_Ah = Qmax_cell_Ah; Results(k).Qmax_cell_Ah_BMS = Qmax_cell_Ah_BMS;
    Results(k).rest1_duration = rest1_duration; Results(k).rest2_duration = rest2_duration;
end

%% Figure 1: Rest 전압 강하 시각화 (2RC 피팅)
fig1 = figure('Name','Rest Voltage Drop Analysis (2021-06-07)','NumberTitle','off');
tl1 = tiledlayout(fig1, 2, 2, 'TileSpacing','compact', 'Padding','compact');

for k = 1:numel(Results)
    % --- 데이터 준비 ---
    chg_end = Results(k).chg_end;
    aftChg = Results(k).aftChg;
    
    start_time = t(chg_end) - minutes(2);
    end_time = t(aftChg);
    
    time_mask = t >= start_time & t <= end_time;
    time_indices = find(time_mask);
    
    if isempty(time_indices)
        continue;
    end
    
    V_rest = Vcell_avg(time_indices);
    t_rest = t(time_indices);
    
    V_chg_end = Vcell_avg(chg_end);
    V_aftChg = Vcell_avg(aftChg);
    deltaV_rest2 = V_chg_end - V_aftChg;
    
    % Subplot 1: 전압 변화 및 2RC 피팅
    ax1 = nexttile(tl1, 1); hold(ax1,'on'); grid(ax1,'on');
    plot(ax1, t_rest, V_rest, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 1.5, 'DisplayName', 'Measured Voltage');
    
    % --- 2RC 피팅 로직 ---
    chg_end_idx = find(time_indices == chg_end, 1);
    aftChg_idx = find(time_indices == aftChg, 1);
    
    Results(k).fit_V_final = NaN; Results(k).fit_A1 = NaN; Results(k).fit_tau1 = NaN;
    Results(k).fit_A2 = NaN; Results(k).fit_tau2 = NaN; Results(k).fit_R_squared = NaN;
    fit_successful = false;

    if ~isempty(chg_end_idx) && ~isempty(aftChg_idx) && aftChg_idx > chg_end_idx
        t_fit_data = t_rest(chg_end_idx:aftChg_idx);
        V_fit_data = V_rest(chg_end_idx:aftChg_idx);
        t_fit_num = seconds(t_fit_data - t_fit_data(1));
        
        try
            V_final_est = V_fit_data(end);
            total_A_est = V_fit_data(1) - V_final_est;
            A1_est = total_A_est * 0.5;
            A2_est = total_A_est * 0.5;
            tau1_est = t_fit_num(end) * 0.1;
            tau2_est = t_fit_num(end) * 0.6;

            ft = fittype('a + b*exp(-x/c) + d*exp(-x/f)', 'independent', 'x', 'dependent', 'y');
            opts = fitoptions('Method', 'NonlinearLeastSquares');
            opts.Display = 'Off';
            opts.StartPoint = [V_final_est, A1_est, tau1_est, A2_est, tau2_est];
            opts.Lower = [-Inf, 0, 0.1, 0, 0.1];
            opts.Upper = [Inf, Inf, Inf, Inf, Inf];
            
            [fitresult, gof] = fit(t_fit_num, V_fit_data, ft, opts);
            V_fit_result = fitresult(t_fit_num);
            
            plot(ax1, t_fit_data, V_fit_result, '--', 'Color', [0.8 0.4 0.2], 'LineWidth', 2, 'DisplayName', '2RC Fit Curve');
            
            Results(k).fit_V_final = fitresult.a;
            Results(k).fit_A1 = fitresult.b;
            Results(k).fit_tau1 = fitresult.c;
            Results(k).fit_A2 = fitresult.d;
            Results(k).fit_tau2 = fitresult.f;
            Results(k).fit_R_squared = gof.rsquare;
            fit_successful = true;

        catch ME
            fprintf('Fitting failed for segment %d: %s\n', k, ME.message);
        end
    end
    
    xline(ax1, t(chg_end), 'r--', 'LineWidth', 1.5, 'Alpha', 0.7, 'DisplayName', 'Chg End');
    xline(ax1, t(aftChg), 'g--', 'LineWidth', 1.5, 'Alpha', 0.7, 'DisplayName', 'Rest End (t2)');
    plot(ax1, t(aftChg), V_aftChg, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'HandleVisibility', 'off');
    
    if fit_successful
        V_final_predicted = Results(k).fit_V_final;
        plot(ax1, [t(chg_end), t(aftChg)], [V_final_predicted, V_final_predicted], 'k:', 'LineWidth', 2, 'DisplayName', 'Predicted Stable OCV');
        plot(ax1, [t(aftChg), t(aftChg)], [V_aftChg, V_final_predicted], 'm-', 'LineWidth', 2.5, 'DisplayName', 'Stabilization Error');
        error_mV = (V_aftChg - V_final_predicted) * 1000;
        text(ax1, t(aftChg) + seconds(30), V_final_predicted, sprintf('%.1f mV', error_mV), 'Color', 'm', 'FontSize', 12, 'FontWeight', 'bold');
    end

    title(ax1, sprintf('Rest Voltage Analysis - Segment %d (R^2=%.4f)', k, Results(k).fit_R_squared));
    xlabel(ax1, 'Time'); ylabel(ax1, 'V_{cell,avg} [V]');
    legend(ax1, 'Location', 'best');
    
    % Subplot 2: 피팅 오차
    ax2 = nexttile(tl1, 2); hold(ax2,'on'); grid(ax2,'on');
    if fit_successful
        t_fit_data = t_rest(chg_end_idx:aftChg_idx);
        V_fit_data = V_rest(chg_end_idx:aftChg_idx);
        t_fit_num = seconds(t_fit_data - t_fit_data(1));
        
        ft_2rc = fittype('a + b*exp(-x/c) + d*exp(-x/f)', 'independent', 'x', 'dependent', 'y');
        fitobj = cfit(ft_2rc, Results(k).fit_V_final, Results(k).fit_A1, Results(k).fit_tau1, Results(k).fit_A2, Results(k).fit_tau2);
        V_fit_calc = fitobj(t_fit_num);
        
        fitting_error = V_fit_data - V_fit_calc;
        plot(ax2, t_fit_data, fitting_error * 1000, '-', 'Color', [0.8 0.2 0.4], 'LineWidth', 1.5);
        yline(ax2, 0, 'k--');
        rmse = sqrt(mean(fitting_error.^2));
        title(ax2, sprintf('Fitting Error - Segment %d (RMSE=%.2f mV)', k, rmse*1000));
        xlabel(ax2, 'Time'); ylabel(ax2, 'Error [mV] (Actual - Fitted)');
    else
        text(ax2, 0.5, 0.5, 'Fitting Failed', 'Units', 'normalized', 'HorizontalAlignment', 'center');
        title(ax2, sprintf('Fitting Error - Segment %d', k));
    end
    
    % Subplot 3: 전압 변화율 (dV/dt)
    ax3 = nexttile(tl1, 3); hold(ax3,'on'); grid(ax3,'on');
    dt_sec = seconds(diff(t_rest));
    dV = diff(V_rest);
    valid_idx = dt_sec > 1e-6;
    dV_dt = zeros(size(dV));
    dV_dt(valid_idx) = dV(valid_idx) ./ dt_sec(valid_idx);
    t_mid = t_rest(1:end-1) + diff(t_rest)/2;
    plot(ax3, t_mid(valid_idx), dV_dt(valid_idx) * 1000 * 60, '-', 'Color', [0.6 0.2 0.6], 'LineWidth', 1.5);
    yline(ax3, 1, 'k--');
    yline(ax3, -1, 'k--');
    title(ax3, sprintf('Voltage Change Rate - Segment %d', k));
    xlabel(ax3, 'Time'); ylabel(ax3, 'dV/dt [mV/min]');
    
    % Subplot 4: 누적 전압 변화량
    ax4 = nexttile(tl1, 4); hold(ax4,'on'); grid(ax4,'on');
    V_change_cumulative = V_chg_end - V_rest;
    plot(ax4, t_rest, V_change_cumulative * 1000, '-', 'Color', [0.2 0.6 0.4], 'LineWidth', 1.5);
    plot(ax4, t(aftChg), deltaV_rest2 * 1000, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    title(ax4, sprintf('Cumulative Voltage Change from Chg End - Seg %d', k));
    xlabel(ax4, 'Time'); ylabel(ax4, 'Cumulative Drop [mV]');
    
    % --- [중요] 누락되었던 결과 저장 부분 ---
    Results(k).V_chg_end = V_chg_end;
    Results(k).V_aftChg = V_aftChg;
    Results(k).deltaV_rest2 = deltaV_rest2;
    Results(k).rest_duration_min = minutes(end_time - start_time); % 오류가 발생했던 필드
    
    % 안정화 판정
    is_voltage_stable = abs(deltaV_rest2) < 0.01;
    is_rate_stable = all(abs(dV_dt(valid_idx)) < (0.001/60)); 
    Results(k).is_voltage_stable = is_voltage_stable;
    Results(k).is_rate_stable = is_rate_stable;
    Results(k).is_fully_stable = is_voltage_stable && is_rate_stable;
end

% Save Figure 1
saveas(fig1, 'FieldQmax/FieldQmax_RestVoltage_2RC_2021.fig');

%% Figure 2: OCV 안정화 요약
fig2 = figure('Name','OCV Stabilization Summary (2021-06-07)','NumberTitle','off');
tl2 = tiledlayout(fig2, 2, 2, 'TileSpacing','compact', 'Padding','compact');

% 모든 세그먼트의 안정화 결과 요약
segment_ids = [Results.idx];
deltaV_values = [Results.deltaV_rest2];
rest_durations = [Results.rest_duration_min];
is_stable = [Results.is_fully_stable];

% Subplot 1: 전압 변화량 분포
ax1 = nexttile(tl2, 1); hold(ax1,'on'); grid(ax1,'on');
bar(ax1, segment_ids, deltaV_values, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'k');
yline(ax1, 0.01, '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 2, 'Alpha', 0.7);
yline(ax1, -0.01, '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 2, 'Alpha', 0.7);
title(ax1, 'Voltage Change by Segment');
xlabel(ax1, 'Segment ID'); ylabel(ax1, 'ΔV_{rest2} [V]');
legend(ax1, 'ΔV', 'Stability Threshold', 'Location', 'best');

% Subplot 2: 휴지 시간 분포
ax2 = nexttile(tl2, 2); hold(ax2,'on'); grid(ax2,'on');
bar(ax2, segment_ids, rest_durations, 'FaceColor', [0.9 0.6 0.3], 'EdgeColor', 'k');
title(ax2, 'Rest Duration by Segment');
xlabel(ax2, 'Segment ID'); ylabel(ax2, 'Rest Duration [min]');

% Subplot 3: 안정화 상태
ax3 = nexttile(tl2, 3); hold(ax3,'on'); grid(ax3,'on');
colors = zeros(length(is_stable), 3);
colors(is_stable, :) = repmat([0.2 0.8 0.2], sum(is_stable), 1); % Green for stable
colors(~is_stable, :) = repmat([0.8 0.2 0.2], sum(~is_stable), 1); % Red for unstable
bar(ax3, segment_ids, double(is_stable), 'FaceColor', 'flat', 'CData', colors, 'EdgeColor', 'k');
title(ax3, 'Stabilization Status by Segment');
xlabel(ax3, 'Segment ID'); ylabel(ax3, 'Stable (1) / Unstable (0)');
ylim(ax3, [-0.1 1.1]);

% Subplot 4: 전압 변화량 vs 휴지 시간
ax4 = nexttile(tl2, 4); hold(ax4,'on'); grid(ax4,'on');
scatter(ax4, rest_durations, abs(deltaV_values), 100, double(is_stable), 'filled', 'MarkerEdgeColor', 'k');
colorbar(ax4, 'Ticks', [0 1], 'TickLabels', {'Unstable', 'Stable'});
title(ax4, 'Voltage Change vs Rest Duration');
xlabel(ax4, 'Rest Duration [min]'); ylabel(ax4, '|ΔV_{rest2}| [V]');

% Save Figure 2
saveas(fig2, 'FieldQmax/FieldQmax_Stabilization_2021.fig');

% Save results
save('FieldQmax/FieldQmax_VoltageAnalysis_Results_2021.mat', 'Results');

fprintf('\n=== Voltage Analysis Complete ===\n');
fprintf('Figure 1: Rest Voltage Drop Analysis saved\n');
fprintf('Figure 2: OCV Stabilization Summary saved\n');
fprintf('Results saved to FieldQmax_VoltageAnalysis_Results_2021.mat\n');

%% Local function: find contiguous true segments
function segs = local_find_segments(mask)
    mask = mask(:)';
    d = diff([false mask false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    segs = [starts(:) ends(:)];
end
