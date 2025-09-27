%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldQmax_v1.m
%
% Rack01 데이터(2021-06-07)에서 충전 구간(연속 10분, I_rack>thr_A)만 검출하고
% 랙 기준 Qmax(Ah)를 OCV 기반 SOC로 추정
% - Idle 조건: |I_rack| < thr_A 
% - dt = 1 s
% - SOC2 충전 종료 이후 첫 rest의 끝 시점, Qmax SOC2는 충전 종료 시점 
% - Rack01
%   (1,1) I_rack vs time, (1,2) P_rack_kW vs time, (2,1) 요약 표(uitable)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% Parameters
Cnom = 128;                         % Rack nominal Capacity (Ah) (참조용)
C_cell_Ah = 64;                     % Cell capacity (Ah)
thr_A = Cnom * 0.05;                % Idle threshold (A)
min_charge_sec = 600;               % 연속 10분 이상 충전 구간
dt = 1;                             % s (고정 가정)

% config (참조: 17 modules * 14S, 2P)
Ns = 17 * 14;                       % 238S (참조)
Np = 2;                             % 2P (참조)

%% Paths
dataFile = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\2021\202106\Raw_20210607.mat';
ocvFile  = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\OCV_integrated.mat';

% Create FieldQmax folder for results
if ~exist('FieldQmax', 'dir')
    mkdir('FieldQmax');
end

%% Load rack data
S = load(dataFile);                 % expects Raw_Rack
if ~isfield(S, 'Raw_Rack')
    error('Raw_Rack not found in %s', dataFile);
end
Raw_Rack = S.Raw_Rack;
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

    dSOC_row = abs(dSOC_q); % percent, befChg->aftChg
    dSOC_row_BMS = abs(dSOC_q_BMS); % percent, befChg->aftChg (BMS)
    rows(end+1, :) = {k, datestr(t(chg_start)), datestr(t(chg_end)), ...
        char(rest1_duration), char(rest2_duration), SOC1, SOC2, dSOC_row, SOC1_raw, SOC2_raw, dSOC_row_BMS, Qmax_cell_Ah, Qmax_cell_Ah_BMS, Q_Ah_cell, Vcell_avg(chg_start), Vcell_avg(chg_end), I_cell(chg_start), I_cell(chg_end), SOH_end, SOH_end}; %#ok<AGROW>
end

%% Visualization
fig = figure('Name','FieldQmax Rack01 (2021-06-07)','NumberTitle','off');
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

% (1,1) Current with SOC overlay (yyaxis)
ax1 = nexttile(tl, 1); hold(ax1,'on'); grid(ax1,'on');

% Plot all data
plot(ax1, t, I_cell, '-', 'Color', [0.8 0.3 0.1], 'LineWidth', 1.2);

% Set xlim to show SOC1-1h to SOC2+1h range
if ~isempty(Results)
    min_time = t(min([Results.befChg])) - hours(1); % SOC1 - 1 hour (datetime)
    max_time = t(max([Results.aftChg])) + hours(1); % SOC2 + 1 hour (datetime)
    xlim(ax1, [min_time, max_time]);
end

title(ax1, 'Cell Current I_{cell} [A]'); xlabel(ax1,'Time'); ylabel(ax1,'I_{cell} [A]');
ylim([-200 150]);
% Add vertical dashed lines for SOC1, SOC2 time points
for k = 1:numel(Results)
    befChg = Results(k).befChg;
    aftChg = Results(k).aftChg;
    
    % Vertical dashed line at SOC1 time point (befChg)
    xline(ax1, t(befChg), '--', 'Color', [0 0 1], 'LineWidth', 2, 'Alpha', 0.7);
    
    % Vertical dashed line at SOC2 time point (aftChg)
    xline(ax1, t(aftChg), '--', 'Color', [0 0 1], 'LineWidth', 2, 'Alpha', 0.7);
end
% SOC overlay per charge segment (using only rest period SOC values)
axes(ax1); yyaxis right; ylabel('SOC [%]');
for k = 1:numel(Results)
    cstart = Results(k).chg_start; cend = Results(k).chg_end;
    befChg = Results(k).befChg; aftChg = Results(k).aftChg;
    SOC1 = Results(k).SOC1; SOC2 = Results(k).SOC2;
    
    % Plot horizontal lines for rest period SOC values
    % befChg to chg_start: constant SOC1
    plot([t(befChg) t(cstart)], [SOC1 SOC1], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 3);
    
    % chg_start to chg_end: linear interpolation between SOC1 and SOC2
    plot([t(cstart) t(cend)], [SOC1 SOC2], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 3);
    
    % chg_end to aftChg: constant SOC2
    plot([t(cend) t(aftChg)], [SOC2 SOC2], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 3);
    
    % markers at rest periods (blue dots)
    plot(t(befChg), SOC1, 'o', 'Color', [0 0 1], 'MarkerFaceColor',[0 0 1], 'MarkerSize', 8);
    plot(t(aftChg), SOC2, 'o', 'Color', [0 0 1], 'MarkerFaceColor',[0 0 1], 'MarkerSize', 8);
end
yyaxis left;
% unify axis colors to black
ax1.XColor = [0 0 0]; if numel(ax1.YAxis)>=1, ax1.YAxis(1).Color = [0 0 0]; end
if numel(ax1.YAxis)>=2, ax1.YAxis(2).Color = [0 0 0]; end

% (2,1) Voltage (avg cell) with SOC overlay (yyaxis)
ax2 = nexttile(tl, 3); hold(ax2,'on'); grid(ax2,'on');

% Plot all data
plot(ax2, t, Vcell_avg, '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 1.2);

% Set xlim to show SOC1-1h to SOC2+1h range
if ~isempty(Results)
    min_time = t(min([Results.befChg])) - hours(1); % SOC1 - 1 hour (datetime)
    max_time = t(max([Results.aftChg])) + hours(1); % SOC2 + 1 hour (datetime)
    xlim(ax2, [min_time, max_time]);
end

title(ax2, 'Average Cell Voltage V_{cell,avg} [V]'); xlabel(ax2,'Time'); ylabel(ax2,'V_{cell,avg} [V]');
ylim([2.7 4.3]);
% Add vertical dashed lines for SOC1, SOC2 time points
for k = 1:numel(Results)
    befChg = Results(k).befChg;
    aftChg = Results(k).aftChg;
    
    % Vertical dashed line at SOC1 time point (befChg)
    xline(ax2, t(befChg), '--', 'Color', [0 0 1], 'LineWidth', 2, 'Alpha', 0.7);
    
    % Vertical dashed line at SOC2 time point (aftChg)
    xline(ax2, t(aftChg), '--', 'Color', [0 0 1], 'LineWidth', 2, 'Alpha', 0.7);
end

% SOC overlay per charge segment (using only rest period SOC values)
axes(ax2); yyaxis right; ylabel('SOC [%]');
for k = 1:numel(Results)
    cstart = Results(k).chg_start; cend = Results(k).chg_end;
    befChg = Results(k).befChg; aftChg = Results(k).aftChg;
    SOC1 = Results(k).SOC1; SOC2 = Results(k).SOC2;
    
    % Plot horizontal lines for rest period SOC values
    % befChg to chg_start: constant SOC1
    plot([t(befChg) t(cstart)], [SOC1 SOC1], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 3);
    
    % chg_start to chg_end: linear interpolation between SOC1 and SOC2
    plot([t(cstart) t(cend)], [SOC1 SOC2], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 3);
    
    % chg_end to aftChg: constant SOC2
    plot([t(cend) t(aftChg)], [SOC2 SOC2], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 3);
    
    % markers at rest periods (blue dots)
    plot(t(befChg), SOC1, 'o', 'Color', [0 0 1], 'MarkerFaceColor',[0 0 1], 'MarkerSize', 8);
    plot(t(aftChg), SOC2, 'o', 'Color', [0 0 1], 'MarkerFaceColor',[0 0 1], 'MarkerSize', 8);
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
        tblData{1,s}  = datestr(t(Results(s).befChg));
        tblData{2,s}  = datestr(t(Results(s).aftChg));
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

% Save figure
saveas(fig, 'FieldQmax/FieldQmax_Figure_2021.fig');

% Save table data to mat file
% Create table with proper variable names
T = array2table(tblData', 'VariableNames', varNames, 'RowNames', colNames);
save('FieldQmax/FieldQmax_Results_2021.mat', 'Results', 'T', 'tblData', 'varNames', 'colNames');

%% Local function: find contiguous true segments
function segs = local_find_segments(mask)
    mask = mask(:)';
    d = diff([false mask false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    segs = [starts(:) ends(:)];
end