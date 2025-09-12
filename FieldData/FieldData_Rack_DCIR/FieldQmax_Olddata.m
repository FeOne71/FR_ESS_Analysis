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
thr_A = Cnom * 0.02;                % Idle threshold (A)
min_charge_sec = 600;               % 연속 10분 이상 충전 구간
dt = 1;                             % s (고정 가정)

% config (참조: 17 modules * 14S, 2P)
Ns = 17 * 14;                       % 238S (참조)
Np = 2;                             % 2P (참조)

%% Paths
dataFile = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\2021\202106\Raw_20210607.mat';
ocvFile  = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\RPT_processed_data.mat';

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
T = load(ocvFile);  % expects OCV_data struct
if isfield(T, 'OCV_data')
    OCV_data = T.OCV_data;
elseif isfield(T, 'RPT_data') && isfield(T.RPT_data, 'ocv_function')
    % Fallback: compose struct-like view
    OCV_data = T.RPT_data;  % may contain ocv_function only
else
    error('OCV data not found in %s', ocvFile);
end

% Build inverse OCV->SOC using avg_ocv_rpt0 and soc_grid when available
SOC_from_OCV = [];
if isfield(OCV_data, 'avg_ocv_rpt0') && isfield(OCV_data, 'soc_grid')
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
elseif isfield(OCV_data, 'ocv_function')
    % ocv_function: SOC->OCV, need inverse via grid
    soc_grid = linspace(0,1,101)';
    ocv_vals = OCV_data.ocv_function(soc_grid);
    [ocv_sorted, uq] = unique(ocv_vals, 'stable');
    soc_sorted = soc_grid(uq) * 100; % percent
    SOC_from_OCV = @(v) interp1(ocv_sorted, soc_sorted, v, 'linear', 'extrap');
else
    error('No suitable OCV fields to construct inverse OCV->SOC.');
end

%% For each charge segment, compute anchors and Qmax (rack)
Results = struct();
rows = {};
for k = 1:size(chgSegs,1)
    chg_start = chgSegs(k,1);
    chg_end   = chgSegs(k,2);

    % previous idle (end before chg_start)
    prevIdleIdx = find(idleSegs(:,2) < chg_start, 1, 'last');
    % next idle (start after chg_end)
    nextIdleIdx = find(idleSegs(:,1) > chg_end, 1, 'first');

    if isempty(prevIdleIdx) || isempty(nextIdleIdx)
        % Skip if either anchor missing (no fallback policy)
        continue;
    end

    idle_before_end = idleSegs(prevIdleIdx,2);
    idle_after_end  = idleSegs(nextIdleIdx,2);
    idle_before_start = idleSegs(prevIdleIdx,1);
    idle_after_start  = idleSegs(nextIdleIdx,1);

    % Use instantaneous voltage at rest end points (no averaging)
    V1 = Vcell_avg(idle_before_end);
    V2 = Vcell_avg(idle_after_end);
    I1 = I_cell(idle_before_end);   % cell current at SOC1/OCV1 time
    I2 = I_cell(idle_after_end);    % cell current at SOC2/OCV2 time

    SOC1 = SOC_from_OCV(V1);   % percent
    SOC2 = SOC_from_OCV(V2);   % percent

    % integrate cell current using trapz (convert time to seconds)
    t_sec = seconds(t - t(1));  % convert to seconds from start
    Q_Ah_cell = abs(trapz(t_sec(chg_start:chg_end), I_cell(chg_start:chg_end)) / 3600);   % charge-only (cell)
    Ah_rest2rest = trapz(t_sec(idle_before_end+1:idle_after_end), I_cell(idle_before_end+1:idle_after_end)) / 3600; % from rest end to next rest end (cell)
    Ah_to_chgend = trapz(t_sec(idle_before_end+1:chg_end), I_cell(idle_before_end+1:chg_end)) / 3600;        % from rest end to charge end (cell)

    % SOC change used for Qmax (from chg_start to chg_end)
    if Ah_rest2rest == 0
        dSOC_q = NaN;
    else
        dSOC_q = (SOC2 - SOC1) * (Ah_to_chgend / Ah_rest2rest);   % percent
    end

    % Qmax uses dSOC between chg_start and chg_end
    if ~isnan(dSOC_q) && dSOC_q ~= 0
        Qmax_cell_Ah = Q_Ah_cell / (abs(dSOC_q)/100);
    else
        Qmax_cell_Ah = NaN;
    end

    % store
    Results(k).idx = k; %#ok<SAGROW>
    Results(k).chg_start = chg_start;
    Results(k).chg_end   = chg_end;
    Results(k).idle_before_start = idle_before_start;
    Results(k).idle_before_end   = idle_before_end;
    Results(k).idle_after_start  = idle_after_start;
    Results(k).idle_after_end    = idle_after_end;
    Results(k).V1 = V1; Results(k).V2 = V2;
    Results(k).I1 = I1; Results(k).I2 = I2;
    Results(k).SOC1 = SOC1; Results(k).SOC2 = SOC2; Results(k).dSOC = abs(dSOC_q);
    Results(k).Q_Ah_cell = Q_Ah_cell; Results(k).Qmax_cell_Ah = Qmax_cell_Ah;

    dSOC_row = abs(dSOC_q); % percent, chg_start->chg_end
    rows(end+1, :) = {k, datestr(t(chg_start)), datestr(t(chg_end)), ...
        SOC1, V1, I1, SOC2, V2, I2, Q_Ah_cell, dSOC_row, Qmax_cell_Ah}; %#ok<AGROW>
end

%% Visualization
fig = figure('Name','FieldQmax Rack01 (2021-06-07)','NumberTitle','off');
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

% (1,1) Current with SOC overlay (yyaxis)
ax1 = nexttile(tl, 1); hold(ax1,'on'); grid(ax1,'on');
plot(ax1, t, I_cell, '-', 'Color', [0.8 0.3 0.1], 'LineWidth', 1.2);
title(ax1, 'Cell Current I_{cell} [A]'); xlabel(ax1,'Time'); ylabel(ax1,'I_{cell} [A]');
ylim([-200 150]);
for k = 1:size(chgSegs,1)
    xs = [t(chgSegs(k,1)) t(chgSegs(k,2))];
    area(ax1, xs, [max(I_cell) max(I_cell)], 'FaceColor',[0.9 0.95 1.0], 'EdgeColor','none', 'ShowBaseLine','off');
end
uistack(findobj(ax1,'Type','line'),'top');
% % SOC overlay per charge segment
% axes(ax1); yyaxis right; ylabel('SOC [%]');
% for k = 1:numel(Results)
%     cstart = Results(k).chg_start; cend = Results(k).chg_end;
%     ibs = Results(k).idle_before_start; ibe = Results(k).idle_before_end;
%     ias = Results(k).idle_after_start;  iae = Results(k).idle_after_end;
%     SOC1 = Results(k).SOC1; SOC2 = Results(k).SOC2;
%     % compute SOC across rest-charge-rest
%     Iseg = I_cell(cstart:cend);
%     if isempty(Iseg)
%         continue;
%     end
%     cumAh = cumtrapz(t_sec(cstart:cend), Iseg) / 3600;
%     if cumAh(end) == 0
%         continue;
%     end
%     soc_charge = SOC1 + (SOC2 - SOC1) * (cumAh / cumAh(end));
%     t_idx_full = ibs:iae;
%     soc_full = nan(size(t_idx_full));
%     if ibs < cstart
%         soc_full(1:(cstart-ibs)) = SOC1;
%     end
%     soc_full((cstart-ibs+1):(cend-ibs+1)) = soc_charge;
%     if cend < iae
%         soc_full((cend-ibs+2):end) = SOC2;
%     end
%     plot(t(t_idx_full), soc_full, '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 3);
%     % markers at charge start/end (SOC1 at cstart, SOC_end at cend)
%     plot(t(cstart), SOC1, 'o', 'Color', [0.45 0.7 0.2], 'MarkerFaceColor',[0.45 0.7 0.2]);
%     plot(t(cend), soc_charge(end), 's', 'Color', [0.45 0.7 0.2], 'MarkerFaceColor',[0.45 0.7 0.2]);
% end
yyaxis left;
% unify axis colors to black
ax1.XColor = [0 0 0]; if numel(ax1.YAxis)>=1, ax1.YAxis(1).Color = [0 0 0]; end
if numel(ax1.YAxis)>=2, ax1.YAxis(2).Color = [0 0 0]; end

% (2,1) Voltage (avg cell) with SOC overlay (yyaxis)
ax2 = nexttile(tl, 3); hold(ax2,'on'); grid(ax2,'on');
plot(ax2, t, Vcell_avg, '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 1.2);
title(ax2, 'Average Cell Voltage V_{cell,avg} [V]'); xlabel(ax2,'Time'); ylabel(ax2,'V_{cell,avg} [V]');
ylim([2.7 4.3]);
for k = 1:size(chgSegs,1)
    xs = [t(chgSegs(k,1)) t(chgSegs(k,2))];
    area(ax2, xs, [max(Vcell_avg) max(Vcell_avg)], 'FaceColor',[0.95 0.95 0.9], 'EdgeColor','none', 'ShowBaseLine','off');
end
uistack(findobj(ax2,'Type','line'),'top');
% % SOC overlay per charge segment
% axes(ax2); yyaxis right; ylabel('SOC [%]');
% for k = 1:numel(Results)
%     cstart = Results(k).chg_start; cend = Results(k).chg_end;
%     ibs = Results(k).idle_before_start; ibe = Results(k).idle_before_end;
%     ias = Results(k).idle_after_start;  iae = Results(k).idle_after_end;
%     SOC1 = Results(k).SOC1; SOC2 = Results(k).SOC2;
%     Iseg = I_cell(cstart:cend);
%     if isempty(Iseg)
%         continue;
%     end
%     cumAh = cumtrapz(t_sec(cstart:cend), Iseg) / 3600;
%     if cumAh(end) == 0
%         continue;
%     end
%     soc_charge = SOC1 + (SOC2 - SOC1) * (cumAh / cumAh(end));
%     t_idx_full = ibs:iae;
%     soc_full = nan(size(t_idx_full));
%     if ibs < cstart
%         soc_full(1:(cstart-ibs)) = SOC1;
%     end
%     soc_full((cstart-ibs+1):(cend-ibs+1)) = soc_charge;
%     if cend < iae
%         soc_full((cend-ibs+2):end) = SOC2;
%     end
%     plot(t(t_idx_full), soc_full, '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 3);
%     % charge start/end markers
%     plot(t(cstart), SOC1, 'o', 'Color', [0.45 0.7 0.2], 'MarkerFaceColor',[0.45 0.7 0.2]);
%     plot(t(cend), soc_charge(end), 's', 'Color', [0.45 0.7 0.2], 'MarkerFaceColor',[0.45 0.7 0.2]);
% end
yyaxis left;
% unify axis colors to black
ax2.XColor = [0 0 0]; if numel(ax2.YAxis)>=1, ax2.YAxis(1).Color = [0 0 0]; end
if numel(ax2.YAxis)>=2, ax2.YAxis(2).Color = [0 0 0]; end

% Right column: summary table spanning both rows
axTbl = nexttile(tl, 2, [2 1]);
posTbl = axTbl.Position; delete(axTbl);
varNames = {'Start','End','SOC1[%]','OCV1[V]','I1[A]','SOC2[%]','OCV2[V]','I2[A]','∫I dt [Ah]','ΔSOC[%]','Qmax_cell [Ah]','SOH [%]', 'Raw SOH (%)'};
numSeg = numel(Results);
tblData = cell(numel(varNames), max(numSeg,1));
colNames = cell(1, max(numSeg,1));
if numSeg == 0
    colNames{1} = 'Rack01';
else
    for s = 1:numSeg
        colNames{s} = sprintf('Rack%02d', s);
        tblData{1,s} = datestr(t(Results(s).chg_start));
        tblData{2,s} = datestr(t(Results(s).chg_end));
        tblData{3,s} = Results(s).SOC1;
        tblData{4,s} = Results(s).V1;
        tblData{5,s} = Results(s).I1;
        tblData{6,s} = Results(s).SOC2;
        tblData{7,s} = Results(s).V2;
        tblData{8,s} = Results(s).I2;
        tblData{9,s} = Results(s).Q_Ah_cell;
        tblData{10,s} = Results(s).dSOC;
        tblData{11,s} = Results(s).Qmax_cell_Ah;
        tblData{12,s} = (Results(s).Qmax_cell_Ah/64)*100;
        tblData{13,s} = 98;
    end
end
uit = uitable(fig, 'Data', tblData, 'ColumnName', colNames, 'RowName', varNames, 'Units','normalized', 'Position', posTbl);
uit.ColumnWidth = 'auto';

set(findall(fig,'-property','FontSize'),'FontSize',12);

%% Local function: find contiguous true segments
function segs = local_find_segments(mask)
    mask = mask(:)';
    d = diff([false mask false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    segs = [starts(:) ends(:)];
end


