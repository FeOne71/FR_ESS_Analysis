%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldQmax_v1.m
%
% Rack01 데이터(2023-10-16, 2025-07-11)에서 충전 구간(연속 10분, I_rack>thr_A)만 검출하고
% 랙 기준 Qmax(Ah)를 OCV 기반 SOC추정
% - Idle 조건: |I_rack| < thr_A 
% - dt = 1 s
% - SOC2 충전 종료 이후 첫 rest의 끝 시점, Qmax SOC2는 충전 종료 시점 
% - Rack01
%   (1,1) I_rack vs time, (1,2) P_rack_kW vs time, (2,1) 요약 표(uitable)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% Parameters
Cnom = 128;                 % ref only
C_cell_Ah = 64;             % Cell capacity (Ah)
thr_A = Cnom*0.02;          % idle threshold A
min_charge_sec = 300;      % 2023
% min_charge_sec = 480;       % 2025
dt = 1;                     % s

Ns = 17*14; 
Np = 2;         % config (ref)

%% Paths
dataFile = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New\2023\202310\Raw_20231016.mat';
% dataFile = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New\2025\202507\Raw_20250711.mat';
ocvFile  = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\RPT_processed_data.mat';

%% Load (robust, same as FieldData_RPT_New)
S = load(dataFile);            % loads variables
if isfield(S,'Raw_Rack')
    Raw_Rack = S.Raw_Rack;
else
    % Find first struct containing expected fields
    vars = fieldnames(S);
    Raw_Rack = [];
    for vi = 1:numel(vars)
        val = S.(vars{vi});
        if isstruct(val)
            f = fieldnames(val);
            if any(strcmpi(f,'CVavg')) && any(strcmpi(f,'DCCurrent'))
                Raw_Rack = val; break;
            end
        end
    end
    if isempty(Raw_Rack)
        error('No suitable struct (with CVavg/DCCurrent) found in %s', dataFile);
    end
end
if isfield(Raw_Rack,'Rack01')
    D = Raw_Rack.Rack01;
else
    D = Raw_Rack;
end

% Time (align with RPT New)
if isfield(D,'Date_Time')
    if isduration(D.Date_Time)
        t = datetime(2023,10,16) + D.Date_Time;
        % t = datetime(2025,07,11) + D.Date_Time;
        
    else
        t = datetime(D.Date_Time);
    end
elseif isfield(D,'Time')
    t = datetime(2023,10,16) + duration(string(D.Time),'InputFormat','hh:mm:ss');
    % t = datetime(2025,07,11) + duration(string(D.Time),'InputFormat','hh:mm:ss');
else
    error('No Date_Time/Time found.');
end
t0 = t(1); tsec = seconds(t - t0);

% Signals (New): DCCurrent, CVavg, DCPower
I_rack = D.DCCurrent(:);
Vcell_avg = D.CVavg(:);
P_rack_kW = D.DCPower(:);
SOC_raw = D.SOC_BMS(:);            % Raw SOC from BMS (%)

% SOH data
SOH_raw = D.SOH_BMS(:);
SOH_end = SOH_raw(end);  % Get the last SOH value

% Convert to cell units
I_cell = I_rack / Np;            % A (cell)
thr_cell = thr_A / Np;           % cell idle threshold

% Masks/segments (cell units)
isIdle = abs(I_cell) < thr_cell; isChg = I_cell > thr_cell;
idleSegs = local_find_segments(isIdle);
chgSegs = local_find_segments(isChg);
chgSegs = chgSegs((chgSegs(:,2)-chgSegs(:,1)+1) >= min_charge_sec, :);

%% OCV inverse (V->SOC%)
T = load(ocvFile);
if isfield(T,'OCV_data')
    OCV_data = T.OCV_data;
elseif isfield(T,'RPT_data') && isfield(T.RPT_data,'ocv_function')
    OCV_data = T.RPT_data;
else
    error('OCV data missing');
end
SOC_from_OCV = [];
if isfield(OCV_data,'avg_ocv_rpt0') && isfield(OCV_data,'soc_grid')
    ocv = OCV_data.avg_ocv_rpt0(:); soc = OCV_data.soc_grid(:);
    [ocv_sorted, uq] = unique(ocv,'stable'); soc_sorted = soc(uq);
    if max(soc_sorted) <= 1.5, soc_sorted = soc_sorted*100; end
    SOC_from_OCV = @(v) interp1(ocv_sorted, soc_sorted, v, 'linear','extrap');
elseif isfield(OCV_data,'ocv_function')
    sg = linspace(0,1,101)'; ov = OCV_data.ocv_function(sg);
    [ocv_sorted, uq] = unique(ov,'stable'); soc_sorted = sg(uq)*100;
    SOC_from_OCV = @(v) interp1(ocv_sorted, soc_sorted, v, 'linear','extrap');
else
    error('No OCV fields to invert');
end

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

    SOC_befChg = SOC_from_OCV(V_befChg);   % 충전 전 휴지 종료 시점 SOC (OCV-based)
    SOC_aftChg = SOC_from_OCV(V_aftChg);   % 충전 후 휴지 종료 시점 SOC (OCV-based)
    
    % Raw SOC data at idle points
    SOC_befChg_raw = SOC_raw(befChg);   % Raw SOC at befChg
    SOC_aftChg_raw = SOC_raw(aftChg);   % Raw SOC at aftChg

    % integrate cell current using trapz (convert time to seconds)
    t_sec = seconds(t - t(1));  % convert to seconds from start
    Q_Ah_cell = abs(trapz(t_sec(chg_start:chg_end), I_cell(chg_start:chg_end)) / 3600);   % charge-only (cell)
    Ah_total = trapz(t_sec(befChg:aftChg), I_cell(befChg:aftChg)) / 3600; % from befChg to aftChg (cell)
    Ah_to_chg_start = trapz(t_sec(befChg:chg_start), I_cell(befChg:chg_start)) / 3600; % from befChg to chg_start (cell)
    Ah_to_chg_end = trapz(t_sec(befChg:chg_end), I_cell(befChg:chg_end)) / 3600; % from befChg to chg_end (cell)

    % SOC(t) calculation: SOC(t) = (SOC_2nd_rest - SOC_1st_rest) × [∫I dt(befChg to t) / ∫I dt(befChg to aftChg)] + SOC_1st_rest
    if Ah_total == 0
        SOC1 = NaN; SOC2 = NaN;
    else
        % SOC at charge start/end points
        SOC1 = (SOC_aftChg - SOC_befChg) * (Ah_to_chg_start / Ah_total) + SOC_befChg;  % chg_start 시점
        SOC2 = (SOC_aftChg - SOC_befChg) * (Ah_to_chg_end / Ah_total) + SOC_befChg;  % chg_end 시점
        
        % Debugging output
        fprintf('=== Charge Segment %d ===\n', k);
        fprintf('befChg: %d, aftChg: %d\n', befChg, aftChg);
        fprintf('chg_start: %d, chg_end: %d\n', chg_start, chg_end);
        fprintf('Time befChg: %s\n', datestr(t(befChg)));
        fprintf('Time chg_start: %s\n', datestr(t(chg_start)));
        fprintf('Time chg_end: %s\n', datestr(t(chg_end)));
        fprintf('Time aftChg: %s\n', datestr(t(aftChg)));
        fprintf('V_befChg: %.3f V, V_aftChg: %.3f V\n', V_befChg, V_aftChg);
        fprintf('SOC_befChg: %.3f%%, SOC_aftChg: %.3f%%\n', SOC_befChg, SOC_aftChg);
        fprintf('Ah_total: %.6f Ah\n', Ah_total);
        fprintf('Ah_to_chg_start: %.6f Ah\n', Ah_to_chg_start);
        fprintf('Ah_to_chg_end: %.6f Ah\n', Ah_to_chg_end);
        fprintf('SOC1: %.3f%%, SOC2: %.3f%%\n', SOC1, SOC2);
        fprintf('dSOC: %.3f%%\n', SOC2 - SOC1);
        fprintf('========================\n\n');
    end
    
    % SOC change used for Qmax (from chg_start to chg_end)
    dSOC_q = SOC2 - SOC1;   % percent

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
    Results(k).befChg = befChg;  % 충전 전 휴지구간 마지막 시점
    Results(k).aftChg = aftChg;  % 충전 후 휴지구간 마지막 시점
    Results(k).V1 = Vcell_avg(chg_start); Results(k).V2 = Vcell_avg(chg_end);
    Results(k).I1 = I_cell(chg_start); Results(k).I2 = I_cell(chg_end);
    Results(k).SOC1 = SOC1; Results(k).SOC2 = SOC2; Results(k).dSOC = abs(dSOC_q);
    Results(k).SOC1_raw = SOC_befChg_raw; Results(k).SOC2_raw = SOC_aftChg_raw; % Raw SOC data
    Results(k).Q_Ah_cell = Q_Ah_cell; Results(k).Qmax_cell_Ah = Qmax_cell_Ah;

    dSOC_row = abs(dSOC_q); % percent, chg_start->chg_end
    rows(end+1, :) = {k, datestr(t(chg_start)), datestr(t(chg_end)), ...
        SOC1, SOC_befChg_raw, Vcell_avg(chg_start), I_cell(chg_start), SOC2, SOC_aftChg_raw, Vcell_avg(chg_end), I_cell(chg_end), Q_Ah_cell, dSOC_row, Qmax_cell_Ah}; %#ok<AGROW>
end

%% Visualization
fig = figure('Name','FieldQmax Rack01 (Newdata)','NumberTitle','off');
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
% SOC overlay per charge segment
axes(ax1); yyaxis right; ylabel('SOC [%]');
for k = 1:numel(Results)
    cstart = Results(k).chg_start; cend = Results(k).chg_end;
    befChg = Results(k).befChg; aftChg = Results(k).aftChg;
    SOC1 = Results(k).SOC1; SOC2 = Results(k).SOC2;
    
    % Recalculate variables for this charge segment
    V_befChg = Vcell_avg(befChg);
    V_aftChg = Vcell_avg(aftChg);
    SOC_befChg = SOC_from_OCV(V_befChg);
    SOC_aftChg = SOC_from_OCV(V_aftChg);
    Ah_total = trapz(t_sec(befChg:aftChg), I_cell(befChg:aftChg)) / 3600;
    
    % compute SOC across rest-charge-rest using new formula
    t_idx_full = befChg:aftChg;
    soc_full = nan(size(t_idx_full));
    
    % SOC at specific time points
    soc_full(1) = SOC_befChg;                    % befChg 시점
    soc_full(chg_start - befChg + 1) = SOC1;     % chg_start 시점  
    soc_full(chg_end - befChg + 1) = SOC2;       % chg_end 시점
    soc_full(end) = SOC_aftChg;                  % aftChg 시점
    
    % 충전 구간에서만 SOC(t) 공식 사용
    for i = (chg_start - befChg + 1):(chg_end - befChg + 1)
        t_idx = befChg + i - 1;
        Ah_to_t = trapz(t_sec(befChg:t_idx), I_cell(befChg:t_idx)) / 3600;
        soc_full(i) = (SOC_aftChg - SOC_befChg) * (Ah_to_t / Ah_total) + SOC_befChg;
    end
    
    plot(t(t_idx_full), soc_full, '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 3);
    % markers at charge start/end (red dots)
    plot(t(cstart), SOC1, 'o', 'Color', [1 0 0], 'MarkerFaceColor',[1 0 0], 'MarkerSize', 8);
    plot(t(cend), SOC2, 'o', 'Color', [1 0 0], 'MarkerFaceColor',[1 0 0], 'MarkerSize', 8);
end
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
yyaxis left;
% unify axis colors to black
ax2.XColor = [0 0 0]; if numel(ax2.YAxis)>=1, ax2.YAxis(1).Color = [0 0 0]; end
if numel(ax2.YAxis)>=2, ax2.YAxis(2).Color = [0 0 0]; end

% Right column: summary table spanning both rows
axTbl = nexttile(tl, 2, [2 1]);
posTbl = axTbl.Position; delete(axTbl);
varNames = {'Start','End','SOC1[%]','Raw_SOC1[%]','OCV1[V]','I1[A]','SOC2[%]','Raw_SOC2[%]','OCV2[V]','I2[A]','∫I dt [Ah]','ΔSOC[%]','Qmax_cell [Ah]','SOH [%]', 'Raw SOH (%)'};
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
        tblData{4,s} = Results(s).SOC1_raw;
        tblData{5,s} = Results(s).V1;
        tblData{6,s} = Results(s).I1;
        tblData{7,s} = Results(s).SOC2;
        tblData{8,s} = Results(s).SOC2_raw;
        tblData{9,s} = Results(s).V2;
        tblData{10,s} = Results(s).I2;
        tblData{11,s} = Results(s).Q_Ah_cell;
        tblData{12,s} = Results(s).dSOC;
        tblData{13,s} = Results(s).Qmax_cell_Ah;
        tblData{14,s} = (Results(s).Qmax_cell_Ah/64)*100;
        tblData{15,s} = SOH_end; % Use actual SOH from data
    end
end
uit = uitable(fig, 'Data', tblData, 'ColumnName', colNames, 'RowName', varNames, 'Units','normalized', 'Position', posTbl);
uit.ColumnWidth = {300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300}; % 고정 너비로 설정
set(findall(fig,'-property','FontSize'),'FontSize',12);

%% Local function: find contiguous true segments
function segs = local_find_segments(mask)
    mask = mask(:)';
    d = diff([false mask false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    segs = [starts(:) ends(:)];
end


