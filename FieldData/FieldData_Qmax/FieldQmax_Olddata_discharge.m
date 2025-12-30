%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldQmax_discharge_v1.m
%
% Rack01 데이터(2021-06-03)에서 방전 구간(연속 10분, I_rack<-thr_A)만 검출하고
% 셀 기준 변환 Qmax(Ah)를 OCV 기반 SOC로 추정
% - Idle 조건: |I_rack| < thr_A 
% - dt = 1 s
% - SOC2: 방전 종료 이후 첫 rest의 끝 시점 
% - Rack01
%   (1,1) I_rack vs time, (1,2) P_rack_kW vs time, (2,1) 요약 표(uitable)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% Parameters
Cnom = 128;                         % Rack nominal Capacity (Ah) (참조용)
C_cell_Ah = 64;                     % Cell capacity (Ah)
thr_A = C_cell_Ah * 0.05;           % Idle threshold (A)
min_discharge_sec = 300;            % 연속 5분 이상 방전 구간
dt = 1;                             % s (고정 가정)

% config (참조: 17 modules * 14S, 2P)
Ns = 17 * 14;                       % 238S (참조)
Np = 2;                             % 2P (참조)

%% Paths
dataFile_0603 = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old\2021\202106\Raw_20210603.mat';
ocvFile  = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\OCV_integrated.mat';

% Create Qmax/discharge folder for results
if ~exist('Qmax', 'dir')
    mkdir('Qmax');
end
if ~exist('Qmax/discharge', 'dir')
    mkdir('Qmax/discharge');
end

%% Process Multiple Old Data Files
dates = {'0603'};
dataFiles = {dataFile_0603};
base_dates = {datetime(2021,6,3)};

for date_idx = 1:length(dates)
    date_str = dates{date_idx};
    dataFile = dataFiles{date_idx};
    base_date = base_dates{date_idx};
    
    fprintf('\n=== Processing 2021-%s Data (Discharge) ===\n', date_str);

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
        t = base_date + D.Date_Time;
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
isDischg  = I_cell < -thr_cell;

% Helper: find contiguous true segments -> Nx2 [start_idx, end_idx]
find_segments = @(mask) local_find_segments(mask);

idleSegs = find_segments(isIdle);
dischgSegs  = find_segments(isDischg);

% Filter discharge segments by duration >= 600 s (dt=1s)
dur = @(seg) seg(:,2) - seg(:,1) + 1;   % samples
min_samples = ceil(min_discharge_sec / dt);
fprintf('Found %d discharge segments before filtering (min duration: %d sec)\n', size(dischgSegs,1), min_discharge_sec);
dischgSegs = dischgSegs(dur(dischgSegs) >= min_samples, :);
fprintf('Found %d discharge segments after filtering (>= %d sec)\n', size(dischgSegs,1), min_discharge_sec);
if size(dischgSegs,1) > 0
    for k = 1:size(dischgSegs,1)
        fprintf('  Discharge segment %d: idx %d-%d (duration: %.1f min, time: %s to %s)\n', ...
            k, dischgSegs(k,1), dischgSegs(k,2), (dischgSegs(k,2)-dischgSegs(k,1)+1)/60, ...
            datestr(t(dischgSegs(k,1)), 'HH:MM:SS'), datestr(t(dischgSegs(k,2)), 'HH:MM:SS'));
    end
end

% Load OCV data and build inverse function (OCV -> SOC)
T = load(ocvFile);
OCV_data = T.OCV_data;

% Build inverse OCV->SOC using avg_ocv_rpt0 and soc_grid (data already clean)
ocv = OCV_data.avg_ocv_rpt0(:);
soc = OCV_data.soc_grid(:);
SOC_from_OCV = @(v) interp1(ocv, soc, v, 'linear', 'extrap');

%% For each discharge segment, compute anchors and Qmax (rack)
Results = [];
rows = {};
for k = 1:size(dischgSegs,1)
    dischg_start = dischgSegs(k,1);
    dischg_end   = dischgSegs(k,2);

    % Find idle segments around discharge
    nextIdleIdx = find(idleSegs(:,1) > dischg_end, 1, 'first');

    if isempty(nextIdleIdx)
        % Skip if next idle segment missing
        continue;
    end

    % 방전 시작점에서 역방향으로 실제 연속 휴지구간 찾기
    befDischg = dischg_start - 1;  % 방전 시작 직전부터 검색
    
    % 역방향으로 연속 휴지구간 찾기 (제한 없이)
    continuous_idle_start = befDischg;
    for i = befDischg:-1:1  % 데이터 시작까지 역방향 검색
        if abs(I_cell(i)) >= thr_cell
            fprintf('  Debug: Stopped reverse search at idx %d, I_cell=%.6f >= thr_cell=%.4f, time=%s\n', ...
                i, abs(I_cell(i)), thr_cell, datestr(t(i), 'HH:MM:SS'));
            break;
        end
        continuous_idle_start = i;
    end
    
    % 실제 연속 휴지구간이 충분히 긴지 확인 (최소 3분)
    actual_rest_duration_sec = befDischg - continuous_idle_start + 1;
    fprintf('  Debug: befDischg=%d (time=%s), continuous_idle_start=%d (time=%s), rest_duration=%.1f min\n', ...
        befDischg, datestr(t(befDischg), 'HH:MM:SS'), continuous_idle_start, datestr(t(continuous_idle_start), 'HH:MM:SS'), ...
        actual_rest_duration_sec/60);
    if actual_rest_duration_sec < 180  % 3분 미만이면 스킵
        fprintf('Skipping segment %d: rest duration too short (%.1f min)\n', k, actual_rest_duration_sec/60);
        continue;
    end
    
    % Define time points
    % befDischg는 이미 정의됨 (방전 시작 직전)
    
    % 방전 종료 후부터 실제 idle threshold를 만족하는 연속 구간 찾기 (SOC2 시점)
    aftDischg = dischg_end;  % 방전 종료 시점부터 시작
    for i = dischg_end+1:idleSegs(nextIdleIdx,2)
        if abs(I_cell(i)) >= thr_cell  % idle threshold 위반시 중단
            break;
        end
        aftDischg = i;  % idle threshold를 만족하는 마지막 시점
    end
    
        % SOC2 휴지구간 최소 시간 확인 (3분 이상) - 초기 확인용
        soc2_rest_duration_initial = aftDischg - dischg_end;
        if soc2_rest_duration_initial < 180  % 3분 미만이면 스킵
            fprintf('Skipping segment %d: SOC2 rest duration too short (%.1f min)\n', k, soc2_rest_duration_initial/60);
            continue;
        end
    
    % Calculate rest period durations
    % First rest period (before discharge) - 실제 연속 휴지구간 길이
    rest1_start_idx = continuous_idle_start;
    rest1_end_idx = befDischg;
    rest1_duration_sec = actual_rest_duration_sec;
    rest1_duration = duration(0, 0, rest1_duration_sec);
    
    % Second rest period (after discharge) - 전체 휴지구간 길이
    rest2_start_idx = idleSegs(nextIdleIdx, 1);
    rest2_end_idx = idleSegs(nextIdleIdx, 2);
    rest2_duration_sec = rest2_end_idx - rest2_start_idx + 1;  % seconds
    
    % 휴지기2가 2시간 이상이면 SOC2를 2시간 마지막 시점으로 제한
    if rest2_duration_sec >= 7200  % 2시간 = 7200초
        aftDischg_2hour = rest2_start_idx + 7200 - 1;  % 2시간 마지막 시점
        if aftDischg_2hour <= rest2_end_idx
            aftDischg = aftDischg_2hour;
            fprintf('Rest2 duration is %.1f hours (>= 2 hours). Limiting SOC2 to 2-hour mark.\n', rest2_duration_sec/3600);
        end
    end
    
    % SOC2 휴지구간 최종 시간 계산 (aftDischg 수정 후)
    soc2_rest_duration = aftDischg - dischg_end;
    
    % === 디버깅: 실제 연속 휴지구간 정보 ===
    fprintf('\n=== Debugging Segment %d (Date: %s) ===\n', k, date_str);
    fprintf('Discharge: idx %d-%d (%.1f min)\n', dischg_start, dischg_end, (dischg_end-dischg_start+1)/60);
    fprintf('Actual continuous idle before discharge: %.1f min (from idx %d to %d)\n', ...
        actual_rest_duration_sec/60, continuous_idle_start, befDischg);
    fprintf('Idle threshold: %.4f A (cell level)\n', thr_cell);
    fprintf('Current at SOC1: %.6f A\n', I_cell(befDischg));
    fprintf('SOC1 time: %s\n', datestr(t(befDischg), 'yyyy-mm-dd HH:MM:SS'));
    fprintf('SOC2 time: %s\n', datestr(t(aftDischg), 'yyyy-mm-dd HH:MM:SS'));
    fprintf('SOC2 rest duration: %.1f min\n', soc2_rest_duration/60);
    fprintf('=== End Debugging ===\n');
    
    % SOC at idle end points (for SOC-OCV data) - aftDischg 수정 후 계산
    V_befDischg = Vcell_avg(befDischg);  % 방전 전 휴지 종료 시점 전압
    V_aftDischg = Vcell_avg(aftDischg);  % 방전 후 휴지 종료 시점 전압 (2시간 제한 적용)
    I_befDischg = I_cell(befDischg);     % 방전 전 휴지 종료 시점 전류
    I_aftDischg = I_cell(aftDischg);     % 방전 후 휴지 종료 시점 전류 (2시간 제한 적용)

    % SOC at idle end points (OCV-based and BMS)
    SOC1 = SOC_from_OCV(V_befDischg);   % 방전 전 휴지 종료 시점 SOC (OCV-based)
    SOC2 = SOC_from_OCV(V_aftDischg);   % 방전 후 휴지 종료 시점 SOC (OCV-based, 2시간 제한 적용)
    SOC1_raw = SOC_raw(befDischg);      % Raw SOC at befDischg
    SOC2_raw = SOC_raw(aftDischg);      % Raw SOC at aftDischg (2시간 제한 적용)

    % integrate cell current using trapz (convert time to seconds)
    t_sec = seconds(t - t(1));  % convert to seconds from start
    Q_Ah_cell = abs(trapz(t_sec(dischg_start:dischg_end), I_cell(dischg_start:dischg_end)) / 3600);   % discharge-only (cell)
    Ah_total = trapz(t_sec(befDischg:aftDischg), I_cell(befDischg:aftDischg)) / 3600; % from befDischg to aftDischg (cell)
    Ah_to_dischg_start = trapz(t_sec(befDischg:dischg_start), I_cell(befDischg:dischg_start)) / 3600; % from befDischg to dischg_start (cell)
    Ah_to_dischg_end = trapz(t_sec(befDischg:dischg_end), I_cell(befDischg:dischg_end)) / 3600; % from befDischg to dischg_end (cell)

    % SOC change used for Qmax (from befDischg to aftDischg - using only rest periods)
    dSOC_q = SOC2 - SOC1;   % percent
    dSOC_q_BMS = SOC2_raw - SOC1_raw;   % percent (BMS)

    % Qmax uses dSOC between rest periods (befDischg to aftDischg)
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
    
    rest2_duration = duration(0, 0, rest2_duration_sec);
    
    % 실제 사용된 SOC2 휴지 시간 (2시간 제한 적용된 값) - 표에 표시할 값
    rest2_duration_used_sec = soc2_rest_duration;  % aftDischg - dischg_end (2시간 제한 적용)
    rest2_duration_used = duration(0, 0, rest2_duration_used_sec);

    % store
    idx = length(Results) + 1;
    Results(idx).idx = k; %#ok<SAGROW>
    Results(idx).dischg_start = dischg_start;
    Results(idx).dischg_end   = dischg_end;
    Results(idx).befDischg = befDischg;  % 방전 전 휴지구간 마지막 시점
    Results(idx).aftDischg = aftDischg;  % 방전 후 휴지구간 마지막 시점 (2시간 제한 적용)
    Results(idx).V1 = Vcell_avg(befDischg); Results(idx).V2 = Vcell_avg(aftDischg);
    Results(idx).I1 = I_cell(befDischg); Results(idx).I2 = I_cell(aftDischg);
    Results(idx).SOC1 = SOC1; Results(idx).SOC2 = SOC2; Results(idx).dSOC = abs(dSOC_q);
    Results(idx).SOC1_raw = SOC1_raw; Results(idx).SOC2_raw = SOC2_raw;
    Results(idx).Q_Ah_cell = Q_Ah_cell; Results(idx).Qmax_cell_Ah = Qmax_cell_Ah; Results(idx).Qmax_cell_Ah_BMS = Qmax_cell_Ah_BMS;
    Results(idx).rest1_duration = rest1_duration; Results(idx).rest2_duration = rest2_duration_used;  % 실제 사용된 시간 표시
    Results(idx).rest1_start_idx = rest1_start_idx; Results(idx).rest1_end_idx = rest1_end_idx;
    Results(idx).rest2_start_idx = rest2_start_idx; Results(idx).rest2_end_idx = rest2_end_idx;
    % 실제 연속 휴지구간 정보 저장
    Results(idx).continuous_idle_start = continuous_idle_start;
    Results(idx).actual_rest_duration_sec = actual_rest_duration_sec;

    dSOC_row = abs(dSOC_q); % percent, befDischg->aftDischg
    dSOC_row_BMS = abs(dSOC_q_BMS); % percent, befDischg->aftDischg (BMS)
    rows(end+1, :) = {k, datestr(t(dischg_start)), datestr(t(dischg_end)), ...
        char(rest1_duration), char(rest2_duration_used), SOC1, SOC2, dSOC_row, SOC1_raw, SOC2_raw, dSOC_row_BMS, Qmax_cell_Ah, Qmax_cell_Ah_BMS, Q_Ah_cell, Vcell_avg(dischg_start), Vcell_avg(dischg_end), I_cell(dischg_start), I_cell(dischg_end), SOH_end, SOH_end}; %#ok<AGROW>
end

%% Visualization
fig = figure('Name',sprintf('FieldQmax Discharge Rack01 (2021-%s)', date_str),'NumberTitle','off');
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

% (1,1) Current with SOC overlay (yyaxis)
ax1 = nexttile(tl, 1); hold(ax1,'on'); grid(ax1,'on');

% Plot all data
plot(ax1, t, I_cell, '-', 'Color', [0.8 0.3 0.1], 'LineWidth', 2.2);

% Set xlim to show SOC1-1h to SOC2+1h range
if ~isempty(Results) && numel(Results) > 0
    min_time = t(min([Results.befDischg])) - hours(1); % SOC1 - 1 hour (datetime)
    max_time = t(max([Results.aftDischg])) + hours(1); % SOC2 + 1 hour (datetime)
    xlim(ax1, [min_time, max_time]);
end

title(ax1, 'Cell Current I_{cell} [A]'); xlabel(ax1,'Time'); ylabel(ax1,'I_{cell} [A]');
ylim([-200 150]);
% Add vertical dashed lines for SOC1, SOC2 time points
for k = 1:numel(Results)
    befDischg = Results(k).befDischg;
    aftDischg = Results(k).aftDischg;
    
    % Vertical dashed line at SOC1 time point (befDischg)
    xline(ax1, t(befDischg), '--', 'Color', [0 0 1], 'LineWidth', 3, 'Alpha', 0.7);
    
    % Vertical dashed line at SOC2 time point (aftDischg)
    xline(ax1, t(aftDischg), '--', 'Color', [0 0 1], 'LineWidth', 3, 'Alpha', 0.7);
end
% SOC overlay per discharge segment (using only rest period SOC values)
axes(ax1); yyaxis right; ylabel('SOC [%]');
for k = 1:numel(Results)
    dstart = Results(k).dischg_start; dend = Results(k).dischg_end;
    befDischg = Results(k).befDischg; aftDischg = Results(k).aftDischg;
    SOC1 = Results(k).SOC1; SOC2 = Results(k).SOC2;
    
    % Plot horizontal lines for rest period SOC values
    % befDischg to dischg_start: constant SOC1
    plot([t(befDischg) t(dstart)], [SOC1 SOC1], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 4);
    
    % dischg_start to dischg_end: linear interpolation between SOC1 and SOC2
    plot([t(dstart) t(dend)], [SOC1 SOC2], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 4);
    
    % dischg_end to aftDischg: constant SOC2
    plot([t(dend) t(aftDischg)], [SOC2 SOC2], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 4);
    
    % markers at rest periods (blue dots)
    plot(t(befDischg), SOC1, 'o', 'Color', [0 0 1], 'MarkerFaceColor',[0 0 1], 'MarkerSize', 8);
    plot(t(aftDischg), SOC2, 'o', 'Color', [0 0 1], 'MarkerFaceColor',[0 0 1], 'MarkerSize', 8);
end
yyaxis left;
% unify axis colors to black
ax1.XColor = [0 0 0]; if numel(ax1.YAxis)>=1, ax1.YAxis(1).Color = [0 0 0]; end
if numel(ax1.YAxis)>=2, ax1.YAxis(2).Color = [0 0 0]; end

% (2,1) Voltage (avg cell) with SOC overlay (yyaxis)
ax2 = nexttile(tl, 3); hold(ax2,'on'); grid(ax2,'on');

% Plot all data
plot(ax2, t, Vcell_avg, '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 2.2);

% Set xlim to show SOC1-1h to SOC2+1h range
if ~isempty(Results) && numel(Results) > 0
    min_time = t(min([Results.befDischg])) - hours(1); % SOC1 - 1 hour (datetime)
    max_time = t(max([Results.aftDischg])) + hours(1); % SOC2 + 1 hour (datetime)
    xlim(ax2, [min_time, max_time]);
end

title(ax2, 'Average Cell Voltage V_{cell,avg} [V]'); xlabel(ax2,'Time'); ylabel(ax2,'V_{cell,avg} [V]');
ylim([2.7 4.3]);
% Add vertical dashed lines for SOC1, SOC2 time points
for k = 1:numel(Results)
    befDischg = Results(k).befDischg;
    aftDischg = Results(k).aftDischg;
    
    % Vertical dashed line at SOC1 time point (befDischg)
    xline(ax2, t(befDischg), '--', 'Color', [0 0 1], 'LineWidth', 3, 'Alpha', 0.7);
    
    % Vertical dashed line at SOC2 time point (aftDischg)
    xline(ax2, t(aftDischg), '--', 'Color', [0 0 1], 'LineWidth', 3, 'Alpha', 0.7);
end

% SOC overlay per discharge segment (using only rest period SOC values)
axes(ax2); yyaxis right; ylabel('SOC [%]');
for k = 1:numel(Results)
    dstart = Results(k).dischg_start; dend = Results(k).dischg_end;
    befDischg = Results(k).befDischg; aftDischg = Results(k).aftDischg;
    SOC1 = Results(k).SOC1; SOC2 = Results(k).SOC2;
    
    % Plot horizontal lines for rest period SOC values
    % befDischg to dischg_start: constant SOC1
    plot([t(befDischg) t(dstart)], [SOC1 SOC1], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 4);
    
    % dischg_start to dischg_end: linear interpolation between SOC1 and SOC2
    plot([t(dstart) t(dend)], [SOC1 SOC2], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 4);
    
    % dischg_end to aftDischg: constant SOC2
    plot([t(dend) t(aftDischg)], [SOC2 SOC2], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 4);
    
    % markers at rest periods (blue dots)
    plot(t(befDischg), SOC1, 'o', 'Color', [0 0 1], 'MarkerFaceColor',[0 0 1], 'MarkerSize', 8);
    plot(t(aftDischg), SOC2, 'o', 'Color', [0 0 1], 'MarkerFaceColor',[0 0 1], 'MarkerSize', 8);
end
yyaxis left;
% unify axis colors to black
ax2.XColor = [0 0 0]; if numel(ax2.YAxis)>=1, ax2.YAxis(1).Color = [0 0 0]; end
if numel(ax2.YAxis)>=2, ax2.YAxis(2).Color = [0 0 0]; end

% Right column: summary table spanning both rows
axTbl = nexttile(tl, 2, [2 1]);
posTbl = axTbl.Position; delete(axTbl);
varNames = {'Start','End','Rest1 [HH:MM:SS]','Rest2 [HH:MM:SS]','SOC1(EST)[%]','SOC2(EST)[%]','ΔSOC(EST)[%]','SOC1(Raw)[%]','SOC2(Raw)[%]','ΔSOC(Raw)[%]','Qmax(EST) [Ah]','Qmax(BMS) [Ah]','∫I dt [Ah]','OCV1[V]','OCV2[V]','I1[A]','I2[A]','SOH(Raw) [%]','SOH(EST) [%]','SOH(BMS) [%]'};
numSeg = numel(Results);
tblData = cell(numel(varNames), max(numSeg,1));
colNames = cell(1, max(numSeg,1));
if numSeg == 0
    colNames{1} = 'Rack01';
else
    for s = 1:numSeg
        colNames{s}   = sprintf('Rack%02d', s);
        tblData{1,s}  = datestr(t(Results(s).befDischg));
        tblData{2,s}  = datestr(t(Results(s).aftDischg));
        tblData{3,s}  = char(Results(s).rest1_duration);
        tblData{4,s}  = char(Results(s).rest2_duration);
        tblData{5,s}  = sprintf('%.4f', Results(s).SOC1);
        tblData{6,s}  = sprintf('%.4f', Results(s).SOC2);
        tblData{7,s}  = sprintf('%.4f', Results(s).dSOC);
        tblData{8,s}  = sprintf('%.4f', Results(s).SOC1_raw);
        tblData{9,s}  = sprintf('%.4f', Results(s).SOC2_raw);
        tblData{10,s} = sprintf('%.4f', abs(Results(s).SOC2_raw - Results(s).SOC1_raw));
        tblData{11,s} = sprintf('%.4f', Results(s).Qmax_cell_Ah);
        tblData{12,s} = sprintf('%.4f', Results(s).Qmax_cell_Ah_BMS);
        tblData{13,s} = sprintf('%.4f', Results(s).Q_Ah_cell);
        tblData{14,s} = sprintf('%.4f', Results(s).V1);
        tblData{15,s} = sprintf('%.4f', Results(s).V2);
        tblData{16,s} = sprintf('%.4f', Results(s).I1);
        tblData{17,s} = sprintf('%.4f', Results(s).I2);
        tblData{18,s} = sprintf('%.4f', SOH_end);
        tblData{19,s} = sprintf('%.4f', (Results(s).Qmax_cell_Ah/64)*100);
        tblData{20,s} = sprintf('%.4f', (Results(s).Qmax_cell_Ah_BMS/64)*100);
    end
end
uit = uitable(fig, 'Data', tblData, 'ColumnName', colNames, 'RowName', varNames, 'Units','normalized', 'Position', posTbl);
uit.ColumnWidth = {300}; % 모든 컬럼 너비 통일
uit.ColumnFormat = {'char', 'char', 'char', 'char', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'}; % 컬럼 포맷
set(findall(fig,'-property','FontSize'),'FontSize',12);
set(findall(fig,'-property','FontWeight'),'FontWeight','bold');

% Save figure
saveas(fig, sprintf('Qmax/discharge/FieldQmax_Discharge_Figure_2021_%s.fig', date_str));

% Save table data to mat file
% Create table with proper variable names
T = array2table(tblData', 'VariableNames', varNames, 'RowNames', colNames);
save(sprintf('Qmax/discharge/FieldQmax_Discharge_Results_2021_%s.mat', date_str), 'Results', 'T', 'tblData', 'varNames', 'colNames');

%% 전압 필터링 (전역 사용을 위해 미리 정의)
window_size_samples = 300; % 300초 윈도우 (1초 샘플링 가정)
Vcell_avg_filtered = movmean(Vcell_avg, window_size_samples);

%% Figure: 원본 vs 필터링된 전압 비교
fig_voltage_comp = figure('Name', sprintf('Voltage Comparison: Raw vs Filtered Discharge (2021-%s)', date_str), 'NumberTitle', 'off');
tl_voltage_comp = tiledlayout(fig_voltage_comp, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (1) 원본 전압
ax1_comp = nexttile(tl_voltage_comp, 1); hold(ax1_comp, 'on'); grid(ax1_comp, 'on');
plot(ax1_comp, t, Vcell_avg, '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 2.2);

% Set xlim to show SOC1-1h to SOC2+1h range
if ~isempty(Results) && numel(Results) > 0
    min_time = t(min([Results.befDischg])) - hours(1);
    max_time = t(max([Results.aftDischg])) + hours(1);
    xlim(ax1_comp, [min_time, max_time]);
end

title(ax1_comp, 'Raw Voltage V_{cell,avg} [V]');
xlabel(ax1_comp, 'Time');
ylabel(ax1_comp, 'V_{cell,avg} [V]');
ylim(ax1_comp, [2.7 4.3]);

% Add vertical dashed lines for SOC1, SOC2 time points
for k = 1:numel(Results)
    befDischg = Results(k).befDischg;
    aftDischg = Results(k).aftDischg;
    xline(ax1_comp, t(befDischg), '--', 'Color', [0 0 1], 'LineWidth', 3, 'Alpha', 0.7, 'Label', 'SOC1');
    xline(ax1_comp, t(aftDischg), '--', 'Color', [0 0 1], 'LineWidth', 3, 'Alpha', 0.7, 'Label', 'SOC2');
end

% (2) 필터링된 전압
ax2_comp = nexttile(tl_voltage_comp, 2); hold(ax2_comp, 'on'); grid(ax2_comp, 'on');
plot(ax2_comp, t, Vcell_avg_filtered, '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 2.2);

% Set xlim to show SOC1-1h to SOC2+1h range
if ~isempty(Results) && numel(Results) > 0
    min_time = t(min([Results.befDischg])) - hours(1);
    max_time = t(max([Results.aftDischg])) + hours(1);
    xlim(ax2_comp, [min_time, max_time]);
end

title(ax2_comp, sprintf('Filtered Voltage V_{cell,avg} [V] (%ds Moving Average)', window_size_samples));
xlabel(ax2_comp, 'Time');
ylabel(ax2_comp, 'V_{cell,avg} [V]');
ylim(ax2_comp, [2.7 4.3]);

% Add vertical dashed lines for SOC1, SOC2 time points
for k = 1:numel(Results)
    befDischg = Results(k).befDischg;
    aftDischg = Results(k).aftDischg;
    xline(ax2_comp, t(befDischg), '--', 'Color', [0 0 1], 'LineWidth', 3, 'Alpha', 0.7, 'Label', 'SOC1');
    xline(ax2_comp, t(aftDischg), '--', 'Color', [0 0 1], 'LineWidth', 3, 'Alpha', 0.7, 'Label', 'SOC2');
end

% Figure 저장
set(findall(fig_voltage_comp,'-property','FontSize'),'FontSize',12);
set(findall(fig_voltage_comp,'-property','FontWeight'),'FontWeight','bold');
saveas(fig_voltage_comp, sprintf('Qmax/discharge/FieldQmax_Discharge_Voltage_Comparison_2021_%s.fig', date_str));

%% Figure: 휴지기별 전압 분석 (2x2 서브플롯)
% 각 세그먼트마다 휴지기 1과 휴지기 2의 원시/필터링 전압 비교
% dV/dt 계산 (각 날짜별로 새로 계산)
dt_sec = seconds(diff(t));
t_mid = t(1:end-1) + diff(t)/2;
dV_raw = diff(Vcell_avg);
dV_dt_raw = dV_raw ./ dt_sec;
dV_filtered = diff(Vcell_avg_filtered);
dV_dt_filtered = dV_filtered ./ dt_sec;

for seg_idx = 1:numel(Results)
    fig_rest_voltage = figure('Name', sprintf('Rest Period Voltage Analysis Discharge: Segment %d (2021-%s)', seg_idx, date_str), 'NumberTitle', 'off');
    tl_rest = tiledlayout(fig_rest_voltage, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % 휴지기 1 정보
    rest1_start = Results(seg_idx).rest1_start_idx;
    rest1_end = Results(seg_idx).rest1_end_idx;
    befDischg = Results(seg_idx).befDischg;
    
    % 휴지기 2 정보
    rest2_start = Results(seg_idx).rest2_start_idx;
    rest2_end = Results(seg_idx).rest2_end_idx;
    dischg_end = Results(seg_idx).dischg_end;
    aftDischg = Results(seg_idx).aftDischg;  % SOC2 시점 (2시간 제한 적용)
    
    % 휴지기 1 인덱스 계산
    rest1_start_time = t(rest1_start);
    start_time_30min_before = rest1_start_time - minutes(30);
    start_idx_30min = find(t >= start_time_30min_before, 1, 'first');
    if isempty(start_idx_30min) || start_idx_30min < 1
        start_idx_30min = 1;
    end
    indices_rest1 = start_idx_30min:rest1_end;
    
    % 휴지기 2 인덱스 계산 (SOC2 시점까지 포함)
    indices_rest2 = dischg_end:max(rest2_end, aftDischg);
    
    % 모든 데이터에서 y축 범위 계산 (전압)
    all_V_rest1_raw = Vcell_avg(indices_rest1);
    all_V_rest1_filtered = Vcell_avg_filtered(indices_rest1);
    all_V_rest2_raw = Vcell_avg(indices_rest2);
    all_V_rest2_filtered = Vcell_avg_filtered(indices_rest2);
    
    all_V_values = [all_V_rest1_raw(:); all_V_rest1_filtered(:); all_V_rest2_raw(:); all_V_rest2_filtered(:)];
    V_min = min(all_V_values);
    V_max = max(all_V_values);
    V_range = V_max - V_min;
    V_ylim = [V_min - V_range*0.05, V_max + V_range*0.05]; % 5% 여유
    
    % dV/dt 범위 계산
    if length(indices_rest1) > 1 && length(indices_rest2) > 1
        dVdt_indices_rest1 = indices_rest1(1:end-1);
        dVdt_indices_rest2 = indices_rest2(1:end-1);
        valid_rest1_dVdt = (t_mid(dVdt_indices_rest1) >= t(indices_rest1(1))) & (t_mid(dVdt_indices_rest1) <= t(indices_rest1(end)));
        valid_rest2_dVdt = (t_mid(dVdt_indices_rest2) >= t(indices_rest2(1))) & (t_mid(dVdt_indices_rest2) <= t(indices_rest2(end)));
        
        all_dVdt_rest1_raw = dV_dt_raw(dVdt_indices_rest1(valid_rest1_dVdt));
        all_dVdt_rest1_filtered = dV_dt_filtered(dVdt_indices_rest1(valid_rest1_dVdt));
        all_dVdt_rest2_raw = dV_dt_raw(dVdt_indices_rest2(valid_rest2_dVdt));
        all_dVdt_rest2_filtered = dV_dt_filtered(dVdt_indices_rest2(valid_rest2_dVdt));
        
        all_dVdt_values = [all_dVdt_rest1_raw(:); all_dVdt_rest1_filtered(:); all_dVdt_rest2_raw(:); all_dVdt_rest2_filtered(:)] * 1000 * 60; % mV/min
        dVdt_min = min(all_dVdt_values);
        dVdt_max = max(all_dVdt_values);
        dVdt_range = dVdt_max - dVdt_min;
        dVdt_ylim = [dVdt_min - dVdt_range*0.05, dVdt_max + dVdt_range*0.05]; % 5% 여유
    else
        dVdt_ylim = [-10, 10]; % 기본값
    end
    
    % (1,1) 휴지기 1의 30분 전까지 원시 전압 + 원시 dV/dt
    ax1 = nexttile(tl_rest, 1); hold(ax1, 'on'); grid(ax1, 'on');
    yyaxis(ax1, 'left');
    plot(ax1, t(indices_rest1), Vcell_avg(indices_rest1), '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 2.2);
    ylabel(ax1, 'V_{cell,avg} [V]');
    ylim(ax1, V_ylim);
    ax1.YAxis(1).Color = [0.5 0.2 0.6];
    
    yyaxis(ax1, 'right');
    if length(indices_rest1) > 1
        dVdt_indices = indices_rest1(1:end-1);
        valid_dVdt = (t_mid(dVdt_indices) >= t(indices_rest1(1))) & (t_mid(dVdt_indices) <= t(indices_rest1(end)));
        if any(valid_dVdt)
            plot(ax1, t_mid(dVdt_indices(valid_dVdt)), dV_dt_raw(dVdt_indices(valid_dVdt)) * 1000 * 60, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
        end
    end
    ylabel(ax1, 'dV/dt [mV/min]');
    ylim(ax1, dVdt_ylim);
    ax1.YAxis(2).Color = [0.2 0.4 0.8];
    
    xline(ax1, t(rest1_start), '--', 'Color', [0 0.8 0], 'LineWidth', 2, 'Alpha', 0.7, 'Label', 'Rest1 Start');
    xline(ax1, t(befDischg), '--', 'Color', [0 0 1], 'LineWidth', 3, 'Alpha', 0.7, 'Label', 'SOC1');
    title(ax1, sprintf('Rest1: Raw Voltage + dV/dt (30min before to end)'));
    xlabel(ax1, 'Time');
    
    % (1,2) 휴지기 1의 30분 전까지 이동평균 적용 전압 + 필터링 dV/dt
    ax2 = nexttile(tl_rest, 2); hold(ax2, 'on'); grid(ax2, 'on');
    yyaxis(ax2, 'left');
    plot(ax2, t(indices_rest1), Vcell_avg_filtered(indices_rest1), '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 2.2);
    ylabel(ax2, 'V_{cell,avg} [V]');
    ylim(ax2, V_ylim);
    ax2.YAxis(1).Color = [0.5 0.2 0.6];
    
    yyaxis(ax2, 'right');
    if length(indices_rest1) > 1
        dVdt_indices = indices_rest1(1:end-1);
        valid_dVdt = (t_mid(dVdt_indices) >= t(indices_rest1(1))) & (t_mid(dVdt_indices) <= t(indices_rest1(end)));
        if any(valid_dVdt)
            plot(ax2, t_mid(dVdt_indices(valid_dVdt)), dV_dt_filtered(dVdt_indices(valid_dVdt)) * 1000 * 60, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
        end
    end
    ylabel(ax2, 'dV/dt [mV/min]');
    ylim(ax2, dVdt_ylim);
    ax2.YAxis(2).Color = [0.2 0.4 0.8];
    
    xline(ax2, t(rest1_start), '--', 'Color', [0 0.8 0], 'LineWidth', 2, 'Alpha', 0.7, 'Label', 'Rest1 Start');
    xline(ax2, t(befDischg), '--', 'Color', [0 0 1], 'LineWidth', 3, 'Alpha', 0.7, 'Label', 'SOC1');
    title(ax2, sprintf('Rest1: Filtered Voltage + dV/dt (30min before to end)'));
    xlabel(ax2, 'Time');
    
    % (2,1) 방전 이후부터 휴지기2까지의 원시 전압 + 원시 dV/dt
    ax3 = nexttile(tl_rest, 3); hold(ax3, 'on'); grid(ax3, 'on');
    yyaxis(ax3, 'left');
    plot(ax3, t(indices_rest2), Vcell_avg(indices_rest2), '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 2.2);
    ylabel(ax3, 'V_{cell,avg} [V]');
    ylim(ax3, V_ylim);
    ax3.YAxis(1).Color = [0.5 0.2 0.6];
    
    yyaxis(ax3, 'right');
    if length(indices_rest2) > 1
        dVdt_indices = indices_rest2(1:end-1);
        valid_dVdt = (t_mid(dVdt_indices) >= t(indices_rest2(1))) & (t_mid(dVdt_indices) <= t(indices_rest2(end)));
        if any(valid_dVdt)
            plot(ax3, t_mid(dVdt_indices(valid_dVdt)), dV_dt_raw(dVdt_indices(valid_dVdt)) * 1000 * 60, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
        end
    end
    ylabel(ax3, 'dV/dt [mV/min]');
    ylim(ax3, dVdt_ylim);
    ax3.YAxis(2).Color = [0.2 0.4 0.8];
    
    xline(ax3, t(dischg_end), '--', 'Color', [1 0 0], 'LineWidth', 3, 'Alpha', 0.7, 'Label', 'Discharge End');
    xline(ax3, t(rest2_start), '--', 'Color', [0 0.8 0], 'LineWidth', 2, 'Alpha', 0.7, 'Label', 'Rest2 Start');
    xline(ax3, t(aftDischg), '--', 'Color', [0 0 1], 'LineWidth', 3, 'Alpha', 0.7, 'Label', 'SOC2');
    title(ax3, sprintf('Rest2: Raw Voltage + dV/dt (Dischg End to Rest2 End)'));
    xlabel(ax3, 'Time');
    
    % (2,2) 방전 이후부터 휴지기2까지의 이동평균 적용 전압 + 필터링 dV/dt
    ax4 = nexttile(tl_rest, 4); hold(ax4, 'on'); grid(ax4, 'on');
    yyaxis(ax4, 'left');
    plot(ax4, t(indices_rest2), Vcell_avg_filtered(indices_rest2), '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 2.2);
    ylabel(ax4, 'V_{cell,avg} [V]');
    ylim(ax4, V_ylim);
    ax4.YAxis(1).Color = [0.5 0.2 0.6];
    
    yyaxis(ax4, 'right');
    if length(indices_rest2) > 1
        dVdt_indices = indices_rest2(1:end-1);
        valid_dVdt = (t_mid(dVdt_indices) >= t(indices_rest2(1))) & (t_mid(dVdt_indices) <= t(indices_rest2(end)));
        if any(valid_dVdt)
            plot(ax4, t_mid(dVdt_indices(valid_dVdt)), dV_dt_filtered(dVdt_indices(valid_dVdt)) * 1000 * 60, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
        end
    end
    ylabel(ax4, 'dV/dt [mV/min]');
    ylim(ax4, dVdt_ylim);
    ax4.YAxis(2).Color = [0.2 0.4 0.8];
    
    xline(ax4, t(dischg_end), '--', 'Color', [1 0 0], 'LineWidth', 3, 'Alpha', 0.7, 'Label', 'Discharge End');
    xline(ax4, t(rest2_start), '--', 'Color', [0 0.8 0], 'LineWidth', 2, 'Alpha', 0.7, 'Label', 'Rest2 Start');
    xline(ax4, t(aftDischg), '--', 'Color', [0 0 1], 'LineWidth', 3, 'Alpha', 0.7, 'Label', 'SOC2');
    title(ax4, sprintf('Rest2: Filtered Voltage + dV/dt (Dischg End to Rest2 End)'));
    xlabel(ax4, 'Time');
    
    % Figure 저장
    set(findall(fig_rest_voltage,'-property','FontSize'),'FontSize',12);
    set(findall(fig_rest_voltage,'-property','FontWeight'),'FontWeight','bold');
    saveas(fig_rest_voltage, sprintf('Qmax/discharge/FieldQmax_Discharge_Rest_Voltage_Segment%d_2021_%s.fig', seg_idx, date_str));
end

%% SOC2 안정화 비교 분석 (긴 휴지기와 비교)
fprintf('\n=== SOC2 Stabilization Analysis (2021-%s, Discharge) ===\n', date_str);

% 긴 휴지기 찾기 (30분 이상)
min_long_idle_sec = 1800; % 30분
long_idle_mask = (idleSegs(:,2) - idleSegs(:,1)) >= min_long_idle_sec;
long_idleSegs = idleSegs(long_idle_mask, :);

fprintf('Found %d long idle periods (>= 30 minutes)\n', size(long_idleSegs,1));

% 각 Results 세그먼트에 대해 비교 분석
for k = 1:numel(Results)
    SOC2_target = Results(k).SOC2;
    aftDischg_idx = Results(k).aftDischg;
    dischg_end_idx = Results(k).dischg_end;
    
    % 현재 SOC2 시점의 휴지 시간 계산
    current_rest_duration_sec = aftDischg_idx - dischg_end_idx;
    
    % SOC2 구간 dV/dt 분석 (필터링된 전압 사용)
    soc2_indices = dischg_end_idx:aftDischg_idx;
    if length(soc2_indices) > 1
        V_soc2_segment = Vcell_avg_filtered(soc2_indices); % 필터링된 전압 사용
        t_soc2_segment = t(soc2_indices);
        
        % dV/dt 계산
        dt_soc2 = seconds(diff(t_soc2_segment));
        dV_soc2 = diff(V_soc2_segment);
        valid_dt = dt_soc2 > 1e-6;
        
        if any(valid_dt)
            dVdt_soc2 = dV_soc2(valid_dt) ./ dt_soc2(valid_dt); % V/s
            dVdt_soc2_mV_min = dVdt_soc2 * 1000 * 60; % mV/min
            
            % 통계 계산
            mean_dVdt = mean(dVdt_soc2_mV_min);
            std_dVdt = std(dVdt_soc2_mV_min);
            max_dVdt = max(dVdt_soc2_mV_min);
            min_dVdt = min(dVdt_soc2_mV_min);
            
            fprintf('  SOC2 dV/dt Analysis (Segment %d):\n', k);
            fprintf('    Rest Duration: %.1f min (%.0f sec)\n', current_rest_duration_sec/60, current_rest_duration_sec);
            fprintf('    Voltage Range: %.3fV -> %.3fV (ΔV=%.1fmV)\n', ...
                V_soc2_segment(1), V_soc2_segment(end), (V_soc2_segment(end)-V_soc2_segment(1))*1000);
            fprintf('    dV/dt Stats [mV/min]: Mean=%.2f, Std=%.2f, Max=%.2f, Min=%.2f\n', ...
                mean_dVdt, std_dVdt, max_dVdt, min_dVdt);
            
            % 안정화 기준 확인
            stability_threshold = 1.0; % mV/min
            stable_points = abs(dVdt_soc2_mV_min) < stability_threshold;
            stability_ratio = sum(stable_points) / length(dVdt_soc2_mV_min) * 100;
            
            fprintf('    Stability: %.1f%% of points < %.1f mV/min\n', stability_ratio, stability_threshold);
            
            % Results에 저장
            Results(k).soc2_dVdt_mean = mean_dVdt;
            Results(k).soc2_dVdt_std = std_dVdt;
            Results(k).soc2_dVdt_max = max_dVdt;
            Results(k).soc2_dVdt_min = min_dVdt;
            Results(k).soc2_stability_ratio = stability_ratio;
            Results(k).soc2_voltage_start = V_soc2_segment(1);
            Results(k).soc2_voltage_end = V_soc2_segment(end);
        else
            fprintf('  SOC2 dV/dt Analysis: No valid time intervals\n');
        end
    else
        fprintf('  SOC2 dV/dt Analysis: Insufficient data points\n');
    end
    
    % SOC2와 유사한 전압 레벨의 긴 휴지기 찾기
    V_SOC2_target = Vcell_avg(aftDischg_idx);  % SOC2 시점의 전압
    voltage_tolerance = 0.1; % ±100mV 허용 범위 (더 넓게 설정)
    matching_long_idles = [];
    
    fprintf('  SOC2 Target Voltage: %.3fV (SOC=%.1f%%)\n', V_SOC2_target, SOC2_target);
    fprintf('  Voltage Tolerance: ±%.0fmV\n', voltage_tolerance*1000);
    
    for j = 1:size(long_idleSegs,1)
        long_idle_start = long_idleSegs(j,1);
        long_idle_end = long_idleSegs(j,2);
        
        % 긴 휴지기 시작점의 전압
        V_long_start = Vcell_avg(long_idle_start);
        voltage_diff = abs(V_long_start - V_SOC2_target);
        
        fprintf('  Long Idle %d: V=%.3fV, Diff=%.0fmV', j, V_long_start, voltage_diff*1000);
        
        % 전압 매칭 확인
        if voltage_diff <= voltage_tolerance
            SOC_long_start = SOC_from_OCV(V_long_start);  % 결과 저장용
            matching_long_idles(end+1,:) = [j, long_idle_start, long_idle_end, SOC_long_start]; %#ok<SAGROW>
            fprintf(' -> MATCHED!\n');
        else
            fprintf(' -> No match\n');
        end
    end
    
    fprintf('  Total Matches Found: %d\n', size(matching_long_idles,1));

% 긴 휴지기들의 dV/dt 분석도 추가
if size(long_idleSegs,1) > 0
    fprintf('  Long Idle Periods dV/dt Analysis:\n');
    for j = 1:size(long_idleSegs,1)
        long_start = long_idleSegs(j,1);
        long_end = long_idleSegs(j,2);
        long_duration_min = (long_end - long_start + 1) / 60;
        
        % 처음 10분간의 dV/dt 계산 (SOC2와 비교 가능하도록)
        analysis_duration = min(600, long_end - long_start + 1); % 10분 또는 전체 기간
        analysis_indices = long_start:(long_start + analysis_duration - 1);
        
        if length(analysis_indices) > 1
            V_long_segment = Vcell_avg_filtered(analysis_indices); % 필터링된 전압 사용
            t_long_segment = t(analysis_indices);
            
            dt_long = seconds(diff(t_long_segment));
            dV_long = diff(V_long_segment);
            valid_dt_long = dt_long > 1e-6;
            
            if any(valid_dt_long)
                dVdt_long = dV_long(valid_dt_long) ./ dt_long(valid_dt_long);
                dVdt_long_mV_min = dVdt_long * 1000 * 60;
                
                mean_dVdt_long = mean(dVdt_long_mV_min);
                std_dVdt_long = std(dVdt_long_mV_min);
                
                fprintf('    Long Idle %d: Duration=%.1fmin, V=%.3fV, dV/dt=%.2f±%.2f mV/min\n', ...
                    j, long_duration_min, Vcell_avg(long_start), mean_dVdt_long, std_dVdt_long);
            end
        end
    end
end
    
    % 매칭되는 긴 휴지기가 있는 경우 비교 분석
    if ~isempty(matching_long_idles)
        % 가장 가까운 전압 레벨의 긴 휴지기 선택
        voltage_differences = abs(Vcell_avg(matching_long_idles(:,2)) - V_SOC2_target);
        [~, best_match_idx] = min(voltage_differences);
        best_match = matching_long_idles(best_match_idx,:);
        
        long_idle_start = best_match(2);
        long_idle_end = best_match(3);
        SOC_long_start = best_match(4);
        
        % 현재 SOC2 시점 데이터
        V_SOC2_start = Vcell_avg(dischg_end_idx);  % 방전 종료 시점 전압
        V_SOC2_end = Vcell_avg(aftDischg_idx);     % SOC2 시점 전압
        
        % 긴 휴지기에서 동일한 시간 경과 후 전압
        comparison_idx = long_idle_start + current_rest_duration_sec;
        if comparison_idx <= long_idle_end
            V_long_start = Vcell_avg(long_idle_start);
            V_long_comparison = Vcell_avg(comparison_idx);
            
            % 전압 변화량 비교
            dV_SOC2 = V_SOC2_end - V_SOC2_start;
            dV_long = V_long_comparison - V_long_start;
            
            % 안정화 오차 계산
            stabilization_error = abs(dV_SOC2) - abs(dV_long);
            
            % 결과 저장
            Results(k).long_idle_match_found = true;
            Results(k).matching_SOC = SOC_long_start;
            Results(k).current_rest_duration_min = current_rest_duration_sec / 60;
            Results(k).SOC2_voltage_drop = dV_SOC2;
            Results(k).long_idle_voltage_drop = dV_long;
            Results(k).stabilization_error_V = stabilization_error;
            Results(k).stabilization_error_mV = stabilization_error * 1000;
            
            % 보정된 SOC2 추정 (긴 휴지기 패턴 기준)
            if abs(dV_long) > 0.001 % 1mV 이상 변화가 있는 경우만
                V_SOC2_corrected = V_SOC2_start + dV_long;
                SOC2_corrected = SOC_from_OCV(V_SOC2_corrected);
                Results(k).SOC2_corrected = SOC2_corrected;
                Results(k).SOC2_correction_needed = abs(SOC2_corrected - SOC2_target) > 0.5; % 0.5% 이상 차이
            else
                Results(k).SOC2_corrected = SOC2_target;
                Results(k).SOC2_correction_needed = false;
            end
            
            fprintf('Segment %d: SOC2=%.1f%%, Match SOC=%.1f%%, Error=%.1fmV\n', ...
                k, SOC2_target, SOC_long_start, stabilization_error*1000);
        else
            % 긴 휴지기가 비교하기에 너무 짧음
            Results(k).long_idle_match_found = false;
            Results(k).SOC2_corrected = SOC2_target;
            Results(k).SOC2_correction_needed = false;
        end
    else
        % 매칭되는 긴 휴지기 없음
        Results(k).long_idle_match_found = false;
        Results(k).SOC2_corrected = SOC2_target;
        Results(k).SOC2_correction_needed = false;
        fprintf('Segment %d: No matching long idle found for SOC2=%.1f%%\n', k, SOC2_target);
    end
end

%% Figure 2: SOC2 안정화 비교 시각화 (간소화)
fig2 = figure('Name', sprintf('SOC2 Stabilization Comparison Discharge (2021-%s)', date_str), 'NumberTitle', 'off');
tl2 = tiledlayout(fig2, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% (1) 전류 & dV/dt 오버랩
ax1 = nexttile(tl2, 1); hold(ax1, 'on'); grid(ax1, 'on');

% 전류 그래프 (왼쪽 y축) - 원본 전류 사용
yyaxis(ax1, 'left');
plot(ax1, t, I_cell, '-', 'Color', [0.8 0.3 0.1], 'LineWidth', 2.2);
ylabel(ax1, 'I_{cell} [A]');
ax1.YAxis(1).Color = [0.8 0.3 0.1];
ylim(ax1, [-200 150]);

% dV/dt 그래프 (오른쪽 y축) - 원시 전압 사용
yyaxis(ax1, 'right');
% 원시 전압으로 dV/dt 계산
dt_sec = seconds(diff(t));
dV_raw = diff(Vcell_avg); % 원시 전압의 차분 사용
dV_dt_raw = dV_raw ./ dt_sec;
t_mid = t(1:end-1) + diff(t)/2;

% xlim 범위에 맞는 dV/dt만 플롯
if ~isempty(Results) && numel(Results) > 0
    min_time = t(min([Results.befDischg])) - hours(1);
    max_time = t(max([Results.aftDischg])) + hours(1);
    valid_time_indices = (t_mid >= min_time) & (t_mid <= max_time);
    plot(ax1, t_mid(valid_time_indices), dV_dt_raw(valid_time_indices) * 1000 * 60, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
else
    plot(ax1, t_mid, dV_dt_raw * 1000 * 60, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
end
ylabel(ax1, 'Raw dV/dt [mV/min]');
ax1.YAxis(2).Color = [0.2 0.4 0.8];

% SOC1, 방전종료, SOC2 시점 표시
for k = 1:numel(Results)
    befDischg = Results(k).befDischg;
    dischg_end = Results(k).dischg_end;
    aftDischg = Results(k).aftDischg;
    xline(ax1, t(befDischg), '--', 'Color', [0 0 1], 'LineWidth', 2.5, 'Alpha', 0.7, 'Label', 'SOC1');
    xline(ax1, t(dischg_end), '--', 'Color', [1 0 0], 'LineWidth', 2.5, 'Alpha', 0.7, 'Label', 'Discharge End');
    xline(ax1, t(aftDischg), '--', 'Color', [0 0.8 0], 'LineWidth', 2.5, 'Alpha', 0.7, 'Label', 'SOC2');
end

% xlim 설정 (원본 Results 기준)
if ~isempty(Results) && numel(Results) > 0
    min_time = t(min([Results.befDischg])) - hours(1);
    max_time = t(max([Results.aftDischg])) + hours(1);
    xlim(ax1, [min_time, max_time]);
end

title(ax1, 'Current & dV/dt Overview');
xlabel(ax1, 'Time');

% (2) 전압 & dV/dt 오버랩
ax2 = nexttile(tl2, 2); hold(ax2, 'on'); grid(ax2, 'on');

% 전압 그래프 (왼쪽 y축) - 원본 전압 사용
yyaxis(ax2, 'left');
plot(ax2, t, Vcell_avg, '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 2.2);
ylabel(ax2, 'V_{cell,avg} [V]');
ax2.YAxis(1).Color = [0.5 0.2 0.6];
ylim(ax2, [2.7 4.3]);

% dV/dt 그래프 (오른쪽 y축) - 원시 전압 사용
yyaxis(ax2, 'right');
% 이미 계산된 원시 dV/dt 데이터 재사용
% xlim 범위에 맞는 dV/dt만 플롯
if ~isempty(Results) && numel(Results) > 0
    min_time = t(min([Results.befDischg])) - hours(1);
    max_time = t(max([Results.aftDischg])) + hours(1);
    valid_time_indices = (t_mid >= min_time) & (t_mid <= max_time);
    plot(ax2, t_mid(valid_time_indices), dV_dt_raw(valid_time_indices) * 1000 * 60, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
else
    plot(ax2, t_mid, dV_dt_raw * 1000 * 60, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
end
ylabel(ax2, 'Raw dV/dt [mV/min]');
ax2.YAxis(2).Color = [0.2 0.4 0.8];

% SOC1, 방전종료, SOC2 시점 표시
for k = 1:numel(Results)
    befDischg = Results(k).befDischg;
    dischg_end = Results(k).dischg_end;
    aftDischg = Results(k).aftDischg;
    xline(ax2, t(befDischg), '--', 'Color', [0 0 1], 'LineWidth', 2.5, 'Alpha', 0.7, 'Label', 'SOC1');
    xline(ax2, t(dischg_end), '--', 'Color', [1 0 0], 'LineWidth', 2.5, 'Alpha', 0.7, 'Label', 'Discharge End');
    xline(ax2, t(aftDischg), '--', 'Color', [0 0.8 0], 'LineWidth', 2.5, 'Alpha', 0.7, 'Label', 'SOC2');
end

% xlim 설정 (Figure 1과 동일)
if ~isempty(Results) && numel(Results) > 0
    min_time = t(min([Results.befDischg])) - hours(1);
    max_time = t(max([Results.aftDischg])) + hours(1);
    xlim(ax2, [min_time, max_time]);
end

title(ax2, 'Voltage & dV/dt Overview');
xlabel(ax2, 'Time');

% Figure 저장
set(findall(fig2,'-property','FontWeight'),'FontWeight','bold');
saveas(fig2, sprintf('Qmax/discharge/FieldQmax_Discharge_SOC2_Comparison_2021_%s.fig', date_str));

% 업데이트된 Results 저장
save(sprintf('Qmax/discharge/FieldQmax_Discharge_Results_Extended_2021_%s.mat', date_str), 'Results', 'T', 'tblData', 'varNames', 'colNames');

%% Figure 3: 필터링된 전압, 원본 전류 시각화
fig3 = figure('Name', sprintf('Filtered Voltage & Raw Current Overview Discharge (2021-%s)', date_str), 'NumberTitle', 'off');

% 원본 전류로 방전/유휴 구간 재계산 (원본 로직과 동일)
% 원본 전류로 마스크 및 구간 재계산
isIdle_filtered = abs(I_cell) < thr_cell;
isDischg_filtered = I_cell < -thr_cell;
idleSegs_filtered = find_segments(isIdle_filtered);
dischgSegs_filtered = find_segments(isDischg_filtered);

% Filter discharge segments by duration >= 600 s (dt=1s)
dur = @(seg) seg(:,2) - seg(:,1) + 1;   % samples
min_samples = ceil(min_discharge_sec / dt);
fprintf('Found %d filtered discharge segments before filtering (min duration: %d sec)\n', size(dischgSegs_filtered,1), min_discharge_sec);
dischgSegs_filtered = dischgSegs_filtered(dur(dischgSegs_filtered) >= min_samples, :);
fprintf('Found %d filtered discharge segments after filtering (>= %d sec)\n', size(dischgSegs_filtered,1), min_discharge_sec);
if size(dischgSegs_filtered,1) > 0
    for k = 1:size(dischgSegs_filtered,1)
        fprintf('  Filtered discharge segment %d: idx %d-%d (duration: %.1f min, time: %s to %s)\n', ...
            k, dischgSegs_filtered(k,1), dischgSegs_filtered(k,2), (dischgSegs_filtered(k,2)-dischgSegs_filtered(k,1)+1)/60, ...
            datestr(t(dischgSegs_filtered(k,1)), 'HH:MM:SS'), datestr(t(dischgSegs_filtered(k,2)), 'HH:MM:SS'));
    end
end

% 필터링된 전류로 SOC1, SOC2 재계산 (원본 로직과 동일)
Results_filtered = [];
fprintf('Original Results count: %d\n', numel(Results));
fprintf('Filtered discharge segments count: %d\n', size(dischgSegs_filtered, 1));

for k = 1:size(dischgSegs_filtered, 1)
    dischg_start = dischgSegs_filtered(k,1);
    dischg_end   = dischgSegs_filtered(k,2);

    % Find idle segments around discharge (필터링된 구간 사용)
    nextIdleIdx = find(idleSegs_filtered(:,1) > dischg_end, 1, 'first');

    if isempty(nextIdleIdx)
        % Skip if next idle segment missing
        continue;
    end

    % 방전 시작점에서 역방향으로 실제 연속 휴지구간 찾기 (원본 전류 기준)
    befDischg_filtered = dischg_start - 1;  % 방전 시작 직전부터 검색
    
    % 역방향으로 연속 휴지구간 찾기 (제한 없이)
    continuous_idle_start_filtered = befDischg_filtered;
    for i = befDischg_filtered:-1:1  % 데이터 시작까지 역방향 검색
        if abs(I_cell(i)) >= thr_cell
            fprintf('  Debug: Stopped reverse search at idx %d, I_cell=%.6f >= thr_cell=%.4f, time=%s\n', ...
                i, abs(I_cell(i)), thr_cell, datestr(t(i), 'HH:MM:SS'));
            break;
        end
        continuous_idle_start_filtered = i;
    end
    
    % 실제 연속 휴지구간이 충분히 긴지 확인 (최소 3분)
    actual_rest_duration_sec = befDischg_filtered - continuous_idle_start_filtered + 1;
    fprintf('  Debug (Filtered): befDischg_filtered=%d (time=%s), continuous_idle_start_filtered=%d (time=%s), rest_duration=%.1f min\n', ...
        befDischg_filtered, datestr(t(befDischg_filtered), 'HH:MM:SS'), continuous_idle_start_filtered, datestr(t(continuous_idle_start_filtered), 'HH:MM:SS'), ...
        actual_rest_duration_sec/60);
    if actual_rest_duration_sec < 180  % 3분 미만이면 스킵
        fprintf('Skipping segment %d: rest duration too short (%.1f min)\n', k, actual_rest_duration_sec/60);
        continue;
    end
    
    % 방전 종료 후부터 실제 idle threshold를 만족하는 연속 구간 찾기 (SOC2 시점)
    aftDischg_filtered = dischg_end;  % 방전 종료 시점부터 시작
    for i = dischg_end+1:idleSegs_filtered(nextIdleIdx,2)
        if abs(I_cell(i)) >= thr_cell  % idle threshold 위반시 중단
            break;
        end
        aftDischg_filtered = i;  % idle threshold를 만족하는 마지막 시점
    end
    
    % SOC2 휴지구간 최소 시간 확인 (3분 이상)
    soc2_rest_duration_filtered = aftDischg_filtered - dischg_end;
    if soc2_rest_duration_filtered < 120  % 3분 미만이면 스킵
        fprintf('Skipping segment %d: SOC2 rest duration too short (%.1f min)\n', k, soc2_rest_duration_filtered/60);
        continue;
    end
    
    % Results_filtered에 추가
    Results_filtered(end+1).befDischg = befDischg_filtered; %#ok<SAGROW>
    Results_filtered(end).aftDischg = aftDischg_filtered;
    Results_filtered(end).dischg_start = dischg_start;
    Results_filtered(end).dischg_end = dischg_end;
    fprintf('Segment %d: Added to Results_filtered (dischg: %d-%d, SOC1: %d, SOC2: %d)\n', ...
        k, dischg_start, dischg_end, befDischg_filtered, aftDischg_filtered);
end

fprintf('Final Results_filtered count: %d\n', numel(Results_filtered));

tl3 = tiledlayout(fig3, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% (1) 원본 전류 & dV/dt 오버랩
ax1 = nexttile(tl3); hold(ax1, 'on'); grid(ax1, 'on');

% 원본 전류 그래프 (왼쪽 y축)
yyaxis(ax1, 'left');
plot(ax1, t, I_cell, '-', 'Color', [0.8 0.3 0.1], 'LineWidth', 2.2);
ylabel(ax1, 'I_{cell} [A]');
ax1.YAxis(1).Color = [0.8 0.3 0.1];
ylim(ax1, [-200 150]);

% dV/dt 그래프 (오른쪽 y축) - 필터링된 전압으로 dV/dt 계산
yyaxis(ax1, 'right');
% 필터링된 전압으로 dV/dt 계산
dV_filtered = diff(Vcell_avg_filtered);
dV_dt_filtered = dV_filtered ./ dt_sec;
% xlim 범위에 맞는 dV/dt만 플롯
if ~isempty(Results_filtered) && numel(Results_filtered) > 0
    min_time = t(min([Results_filtered.befDischg])) - hours(1);
    max_time = t(max([Results_filtered.aftDischg])) + hours(1);
    valid_time_indices = (t_mid >= min_time) & (t_mid <= max_time);
    plot(ax1, t_mid(valid_time_indices), dV_dt_filtered(valid_time_indices) * 1000 * 60, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
else
    plot(ax1, t_mid, dV_dt_filtered * 1000 * 60, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
end
ylabel(ax1, 'Filtered dV/dt [mV/min]');
ax1.YAxis(2).Color = [0.2 0.4 0.8];

% SOC1, 방전종료, SOC2 시점 표시 (원본 전류 기준)
for k = 1:numel(Results_filtered)
    befDischg = Results_filtered(k).befDischg;
    dischg_end = Results_filtered(k).dischg_end;
    aftDischg = Results_filtered(k).aftDischg;
    xline(ax1, t(befDischg), '--', 'Color', [0 0 1], 'LineWidth', 2.5, 'Alpha', 0.7, 'Label', 'SOC1');
    xline(ax1, t(dischg_end), '--', 'Color', [1 0 0], 'LineWidth', 2.5, 'Alpha', 0.7, 'Label', 'Discharge End');
    xline(ax1, t(aftDischg), '--', 'Color', [0 0.8 0], 'LineWidth', 2.5, 'Alpha', 0.7, 'Label', 'SOC2');
end

% xlim 설정 (필터링된 Results 기준)
if ~isempty(Results_filtered) && numel(Results_filtered) > 0
    min_time = t(min([Results_filtered.befDischg])) - hours(1);
    max_time = t(max([Results_filtered.aftDischg])) + hours(1);
    xlim(ax1, [min_time, max_time]);
end

title(ax1, 'Raw Current & dV/dt Overview');
xlabel(ax1, 'Time');

% (2) 필터링된 전압 & dV/dt 오버랩
ax2 = nexttile(tl3); hold(ax2, 'on'); grid(ax2, 'on');

% 필터링된 전압 그래프 (왼쪽 y축)
yyaxis(ax2, 'left');
plot(ax2, t, Vcell_avg_filtered, '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 2.2);
ylabel(ax2, 'V_{cell,avg} [V] (Filtered)');
ax2.YAxis(1).Color = [0.5 0.2 0.6];
ylim(ax2, [2.7 4.3]);

% dV/dt 그래프 (오른쪽 y축) - 필터링된 전압으로 dV/dt 계산
yyaxis(ax2, 'right');
% 이미 계산된 필터링된 dV/dt 데이터 재사용
% xlim 범위에 맞는 dV/dt만 플롯
if ~isempty(Results_filtered) && numel(Results_filtered) > 0
    min_time = t(min([Results_filtered.befDischg])) - hours(1);
    max_time = t(max([Results_filtered.aftDischg])) + hours(1);
    valid_time_indices = (t_mid >= min_time) & (t_mid <= max_time);
    plot(ax2, t_mid(valid_time_indices), dV_dt_filtered(valid_time_indices) * 1000 * 60, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
else
    plot(ax2, t_mid, dV_dt_filtered * 1000 * 60, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
end
ylabel(ax2, 'Filtered dV/dt [mV/min]');
ax2.YAxis(2).Color = [0.2 0.4 0.8];

% SOC1, 방전종료, SOC2 시점 표시 (원본 전류 기준)
for k = 1:numel(Results_filtered)
    befDischg = Results_filtered(k).befDischg;
    dischg_end = Results_filtered(k).dischg_end;
    aftDischg = Results_filtered(k).aftDischg;
    xline(ax2, t(befDischg), '--', 'Color', [0 0 1], 'LineWidth', 2.5, 'Alpha', 0.7, 'Label', 'SOC1');
    xline(ax2, t(dischg_end), '--', 'Color', [1 0 0], 'LineWidth', 2.5, 'Alpha', 0.7, 'Label', 'Discharge End');
    xline(ax2, t(aftDischg), '--', 'Color', [0 0.8 0], 'LineWidth', 2.5, 'Alpha', 0.7, 'Label', 'SOC2');
end

% xlim 설정 (필터링된 Results 기준)
if ~isempty(Results_filtered) && numel(Results_filtered) > 0
    min_time = t(min([Results_filtered.befDischg])) - hours(1);
    max_time = t(max([Results_filtered.aftDischg])) + hours(1);
    xlim(ax2, [min_time, max_time]);
end

title(ax2, 'Filtered Voltage & dV/dt Overview');
xlabel(ax2, 'Time');

% Figure 저장
set(findall(fig3,'-property','FontWeight'),'FontWeight','bold');
saveas(fig3, sprintf('Qmax/discharge/FieldQmax_Discharge_Filtered_Overview_2021_%s.fig', date_str));

%% Figure 4: 방전 시작 전/후 dV/dt 분석 (필터링된 데이터 기준)
fig4 = figure('Name', sprintf('Pre/Post Discharge dV/dt Analysis (2021-%s)', date_str), 'NumberTitle', 'off');

    % dV/dt 계산 (각 날짜별로 새로 계산)
    dt_sec = seconds(diff(t));
    t_mid = t(1:end-1) + diff(t)/2;
    dV_filtered = diff(Vcell_avg_filtered);
    dV_dt_filtered = dV_filtered ./ dt_sec;

% 각 세그먼트에 대해 서브플롯 생성 (필터링된 Results 사용)
num_segments = numel(Results_filtered);
if num_segments > 0
    % 세그먼트 수에 따라 레이아웃 결정 (각 세그먼트당 2행 서브플롯)
    % 각 세그먼트는 2행씩 차지하므로, 전체 행 수는 세그먼트 수 * 2
    rows = num_segments * 2;  % 각 세그먼트당 2행
    cols = 1;  % 1열
    
    tl4 = tiledlayout(fig4, rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for k = 1:num_segments
        % 각 세그먼트마다 2행 1열 서브플롯 생성
        % 서브플롯1: 방전 시작 30분 전부터 방전 시작까지 (위쪽)
        ax1 = nexttile(tl4, (k-1)*2+1); hold(ax1, 'on'); grid(ax1, 'on');
        
        dischg_start_idx = Results_filtered(k).dischg_start;
        befDischg_idx = Results_filtered(k).befDischg;
        
        % 방전 시작 30분 전부터 방전 시작까지의 인덱스
        dischg_start_time = t(dischg_start_idx);
        start_time = dischg_start_time - minutes(30);
        start_idx = find(t >= start_time, 1, 'first');
        if isempty(start_idx) || start_idx < 1
            start_idx = 1;
        end
        pre_dischg_indices = start_idx:dischg_start_idx;
        
        if length(pre_dischg_indices) > 1
            t_pre = t(pre_dischg_indices);
            
            % 전체 dV/dt에서 해당 구간 찾기
            pre_start_time = t_pre(1);
            pre_end_time = t_pre(end);
            valid_indices_pre = (t_mid >= pre_start_time) & (t_mid <= pre_end_time);
            
                if any(valid_indices_pre)
                t_mid_pre = t_mid(valid_indices_pre);
                dVdt_pre_mV_min = dV_dt_filtered(valid_indices_pre) * 1000 * 60; % mV/min
                
                % 시간을 분 단위로 변환 (시작 시점을 0으로)
                t_mid_pre_min = minutes(t_mid_pre - t_pre(1));
                
                % dV/dt 플롯
                plot(ax1, t_mid_pre_min, dVdt_pre_mV_min, '-o', ...
                    'Color', [0.2 0.4 0.8], 'LineWidth', 2, 'MarkerSize', 2);
                
                % 안정화 기준선 표시
                yline(ax1, 1, 'k--', 'LineWidth', 2, 'Alpha', 0.7, 'Label', '1mV/min', 'LabelHorizontalAlignment', 'right');
                yline(ax1, -1, 'k--', 'LineWidth', 2, 'Alpha', 0.7, 'Label', '-1mV/min', 'LabelHorizontalAlignment', 'right');
                yline(ax1, 0, 'k-', 'LineWidth', 1.5, 'Alpha', 0.5);
                
                % y축 범위 고정
                ylim(ax1, [-10, 10]);
                
                % 방전 시작 시점 표시
                dischg_start_time_min = minutes(t(dischg_start_idx) - t_pre(1));
                xline(ax1, dischg_start_time_min, '--', 'Color', [1 0 0], 'LineWidth', 2, 'Alpha', 0.7, 'Label', 'Discharge Start');
            else
                text(ax1, 0.5, 0.5, 'No Valid dV/dt Data', 'Units', 'normalized', ...
                    'HorizontalAlignment', 'center', 'FontSize', 12);
            end
            
            title(ax1, sprintf('Segment %d: dV/dt (30min before Dischg Start)', k));
            xlabel(ax1, 'Time from Start [min]');
            ylabel(ax1, 'dV/dt [mV/min]');
        else
            text(ax1, 0.5, 0.5, 'Insufficient Data Points', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 12);
            title(ax1, sprintf('Segment %d: No Data', k));
        end
        
        % 서브플롯2: 방전 종료부터 SOC2까지 (아래쪽)
        ax2 = nexttile(tl4, (k-1)*2+2); hold(ax2, 'on'); grid(ax2, 'on');
        
        dischg_end_idx = Results_filtered(k).dischg_end;
        aftDischg_idx = Results_filtered(k).aftDischg;
        
        % 방전 종료 시점부터 SOC2시점까지의 인덱스
        rest_indices = dischg_end_idx:aftDischg_idx;
        
        if length(rest_indices) > 1
            t_rest = t(rest_indices);
            
            % 전체 dV/dt에서 해당 구간 찾기
            rest_start_time = t_rest(1);
            rest_end_time = t_rest(end);
            valid_indices_rest = (t_mid >= rest_start_time) & (t_mid <= rest_end_time);
            
            if any(valid_indices_rest)
                t_mid_rest = t_mid(valid_indices_rest);
                dVdt_rest_mV_min = dV_dt_filtered(valid_indices_rest) * 1000 * 60; % mV/min
                
                % 시간을 분 단위로 변환 (시작 시점을 0으로)
                t_mid_rest_min = minutes(t_mid_rest - t_rest(1));
                
                % dV/dt 플롯
                plot(ax2, t_mid_rest_min, dVdt_rest_mV_min, '-o', ...
                    'Color', [0.2 0.4 0.8], 'LineWidth', 2, 'MarkerSize', 2);
                
                % 안정화 기준선 표시
                yline(ax2, 1, 'k--', 'LineWidth', 2, 'Alpha', 0.7, 'Label', '1mV/min', 'LabelHorizontalAlignment', 'right');
                yline(ax2, -1, 'k--', 'LineWidth', 2, 'Alpha', 0.7, 'Label', '-1mV/min', 'LabelHorizontalAlignment', 'right');
                yline(ax2, 0, 'k-', 'LineWidth', 1.5, 'Alpha', 0.5);
                
                % y축 범위 고정
                ylim(ax2, [-10, 10]);
                
                % 방전종료, SOC2 시점 표시 (분 단위로 변환)
                dischg_end_time_min = minutes(t(dischg_end_idx) - t_rest(1));
                aftDischg_time_min = minutes(t(aftDischg_idx) - t_rest(1));
                xline(ax2, dischg_end_time_min, '--', 'Color', [1 0 0], 'LineWidth', 2, 'Alpha', 0.7, 'Label', 'Discharge End');
                xline(ax2, aftDischg_time_min, '--', 'Color', [0 0.8 0], 'LineWidth', 2, 'Alpha', 0.7, 'Label', 'SOC2');
            else
                text(ax2, 0.5, 0.5, 'No Valid dV/dt Data', 'Units', 'normalized', ...
                    'HorizontalAlignment', 'center', 'FontSize', 12);
            end
            
            title(ax2, sprintf('Segment %d: dV/dt (Dischg End to SOC2)', k));
            xlabel(ax2, 'Time from Discharge End [min]');
            ylabel(ax2, 'dV/dt [mV/min]');
        else
            text(ax2, 0.5, 0.5, 'Insufficient Data Points', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 12);
            title(ax2, sprintf('Segment %d: No Data', k));
        end
    end
else
    % 세그먼트가 없는 경우
    ax = axes(fig4);
    text(ax, 0.5, 0.5, 'No Discharge Segments Found', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 16);
    title(ax, 'dV/dt Analysis');
end

% Figure 저장
set(findall(fig4,'-property','FontWeight'),'FontWeight','bold');
saveas(fig4, sprintf('Qmax/discharge/FieldQmax_Discharge_SOC2_Analysis_2021_%s.fig', date_str));

% 최종 Results 저장
save(sprintf('Qmax/discharge/FieldQmax_Discharge_Results_Final_2021_%s.mat', date_str), 'Results', 'T', 'tblData', 'varNames', 'colNames');

end % for date_idx


%% Local function: find contiguous true segments
function segs = local_find_segments(mask)
    mask = mask(:)';
    d = diff([false mask false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    segs = [starts(:) ends(:)];
end

