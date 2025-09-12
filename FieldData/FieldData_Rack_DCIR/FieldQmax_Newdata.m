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
min_charge_sec = 300;       % 5 minutes
dt = 1;                     % s

Ns = 17*14; Np = 2;         % config (ref)

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
        % t = datetime(2025,07,11) + D.Date_Time;6
        
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

%% Iterate segments
Results = struct();
for k = 1:size(chgSegs,1)
    cs = chgSegs(k,1); ce = chgSegs(k,2);
    prevIdleIdx = find(idleSegs(:,2) < cs, 1, 'last');
    nextIdleIdx = find(idleSegs(:,1) > ce, 1, 'first');
    if isempty(prevIdleIdx) || isempty(nextIdleIdx), continue; end
    ibe = idleSegs(prevIdleIdx,2); iae = idleSegs(nextIdleIdx,2);
    ibs = idleSegs(prevIdleIdx,1); ias = idleSegs(nextIdleIdx,1);

    V1 = Vcell_avg(ibe); V2 = Vcell_avg(iae);
    I1 = I_cell(ibe);    I2 = I_cell(iae);
    SOC1 = SOC_from_OCV(V1); SOC2 = SOC_from_OCV(V2);

    % Convert time to seconds for trapz
    t_sec = seconds(t - t(1));
    Q_Ah_cell = abs(trapz(t_sec(cs:ce), I_cell(cs:ce)) / 3600);
    Ah_rest2rest = trapz(t_sec(ibe+1:iae), I_cell(ibe+1:iae)) / 3600;
    Ah_to_chgend = trapz(t_sec(ibe+1:ce), I_cell(ibe+1:ce)) / 3600;
    if Ah_rest2rest==0
        dSOC_q = NaN;
    else
        dSOC_q = (SOC2 - SOC1) * (Ah_to_chgend / Ah_rest2rest);
    end
    if ~isnan(dSOC_q) && dSOC_q~=0
        Qmax_cell_Ah = Q_Ah_cell / (abs(dSOC_q)/100);
    else
        Qmax_cell_Ah = NaN;
    end

    Results(k).chg_start = cs; Results(k).chg_end = ce;
    Results(k).idle_before_start = ibs; Results(k).idle_before_end = ibe;
    Results(k).idle_after_start  = ias; Results(k).idle_after_end  = iae;
    Results(k).SOC1 = SOC1; Results(k).SOC2 = SOC2; Results(k).dSOC = abs(dSOC_q);
    Results(k).V1 = V1; Results(k).V2 = V2; Results(k).I1 = I1; Results(k).I2 = I2;
    Results(k).Q_Ah_cell = Q_Ah_cell; Results(k).Qmax_cell_Ah = Qmax_cell_Ah;
end

%% Plot & table (same style)
fig = figure('Name','FieldQmax Newdata (Rack01)','NumberTitle','off');
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

% Current + SOC
ax1 = nexttile(tl,1); hold(ax1,'on'); grid(ax1,'on');
plot(ax1, t, I_cell, '-', 'Color',[0.8 0.3 0.1],'LineWidth',1.2);
title(ax1,'Cell Current I_{cell} [A]'); xlabel(ax1,'Time'); ylabel(ax1,'I_{cell} [A]');
ylim([-200 150]);
for k = 1:size(chgSegs,1), xs=[t(chgSegs(k,1)) t(chgSegs(k,2))]; area(ax1,xs,[max(I_rack) max(I_rack)],'FaceColor',[0.9 0.95 1.0],'EdgeColor','none','ShowBaseLine','off'); end
uistack(findobj(ax1,'Type','line'),'top');
axes(ax1); yyaxis right; ylabel('SOC [%]');
validIdx = find(arrayfun(@(s) isfield(s,'chg_start'), Results));
for ii = 1:numel(validIdx)
    k = validIdx(ii);
    cs=Results(k).chg_start; ce=Results(k).chg_end; ibs=Results(k).idle_before_start; ibe=Results(k).idle_before_end; ias=Results(k).idle_after_start; iae=Results(k).idle_after_end;
    SOC1=Results(k).SOC1; SOC2=Results(k).SOC2; Iseg=I_cell(cs:ce);
    if isempty(Iseg), continue; end
    cumAh=cumtrapz(t_sec(cs:ce), Iseg)/3600; if cumAh(end)==0, continue; end
    soc_ch = SOC1 + (SOC2-SOC1)*(cumAh/cumAh(end));
    idx = ibs:iae; soc_full = nan(size(idx));
    if ibs<cs, soc_full(1:(cs-ibs)) = SOC1; end
    soc_full((cs-ibs+1):(ce-ibs+1)) = soc_ch;
    if ce<iae, soc_full((ce-ibs+2):end) = SOC2; end
    plot(t(idx), soc_full, '-', 'Color',[0.45 0.7 0.2],'LineWidth',1.5);
    plot(t(cs), SOC1, 'o','Color',[0.45 0.7 0.2],'MarkerFaceColor',[0.45 0.7 0.2]);
    plot(t(ce), soc_ch(end), 's','Color',[0.45 0.7 0.2],'MarkerFaceColor',[0.45 0.7 0.2]);
end
ylim([10 70]); yyaxis left; ylim([-300 150]); ax1.XColor=[0 0 0]; if numel(ax1.YAxis)>=1, ax1.YAxis(1).Color=[0 0 0]; end; if numel(ax1.YAxis)>=2, ax1.YAxis(2).Color=[0 0 0]; end

% Voltage (avg cell) + SOC
ax2 = nexttile(tl,3); hold(ax2,'on'); grid(ax2,'on');
plot(ax2, t, Vcell_avg, '-', 'Color',[0.5 0.2 0.6],'LineWidth',1.2);
title(ax2,'Average Cell Voltage V_{cell,avg} [V]'); xlabel(ax2,'Time'); ylabel(ax2,'V_{cell,avg} [V]');
ylim([2.7 4.3]);
for k = 1:size(chgSegs,1), xs=[t(chgSegs(k,1)) t(chgSegs(k,2))]; area(ax2,xs,[max(Vcell_avg) max(Vcell_avg)],'FaceColor',[0.95 0.95 0.9],'EdgeColor','none','ShowBaseLine','off'); end
uistack(findobj(ax2,'Type','line'),'top');
axes(ax2); yyaxis right; ylabel('SOC [%]');
validIdx = find(arrayfun(@(s) isfield(s,'chg_start'), Results));
for ii = 1:numel(validIdx)
    k = validIdx(ii);
    cs=Results(k).chg_start; ce=Results(k).chg_end; ibs=Results(k).idle_before_start; ibe=Results(k).idle_before_end; ias=Results(k).idle_after_start; iae=Results(k).idle_after_end;
    SOC1=Results(k).SOC1; SOC2=Results(k).SOC2; Iseg=I_cell(cs:ce);
    if isempty(Iseg), continue; end
    cumAh=cumtrapz(t_sec(cs:ce), Iseg)/3600; if cumAh(end)==0, continue; end
    soc_ch = SOC1 + (SOC2-SOC1)*(cumAh/cumAh(end));
    idx = ibs:iae; soc_full = nan(size(idx));
    if ibs<cs, soc_full(1:(cs-ibs)) = SOC1; end
    soc_full((cs-ibs+1):(ce-ibs+1)) = soc_ch;
    if ce<iae, soc_full((ce-ibs+2):end) = SOC2; end
    plot(t(idx), soc_full, '-', 'Color',[0.45 0.7 0.2],'LineWidth',1.5);
    plot(t(cs), SOC1, 'o','Color',[0.45 0.7 0.2],'MarkerFaceColor',[0.45 0.7 0.2]);
    plot(t(ce), soc_ch(end), 's','Color',[0.45 0.7 0.2],'MarkerFaceColor',[0.45 0.7 0.2]);
end
ylim([10 70]); yyaxis left; ax2.XColor=[0 0 0]; if numel(ax2.YAxis)>=1, ax2.YAxis(1).Color=[0 0 0]; end; if numel(ax2.YAxis)>=2, ax2.YAxis(2).Color=[0 0 0]; end

% Table right
axTbl = nexttile(tl,2,[2 1]); posTbl = axTbl.Position; delete(axTbl);
vars = {'Start','End','SOC1[%]','OCV1[V]','I1[A]','SOC2[%]','OCV2[V]','I2[A]','∫I dt [Ah]','ΔSOC[%]','Qmax_cell [Ah]', 'SOH [%]', 'Raw SOH [%]'};
numSeg = numel(Results); data = cell(numel(vars), max(numSeg,1)); cols = cell(1, max(numSeg,1));
if numSeg==0
    cols{1}='Rack01';
else
    for s=1:numSeg
        cols{s}=sprintf('Rack%02d',s);
        data{1,s}=datestr(t(Results(s).chg_start)); data{2,s}=datestr(t(Results(s).chg_end));
        % Convert all numeric values to double, handle NaN
        data{3,s}=double(Results(s).SOC1); data{4,s}=double(Results(s).V1); data{5,s}=double(Results(s).I1);
        data{6,s}=double(Results(s).SOC2); data{7,s}=double(Results(s).V2); data{8,s}=double(Results(s).I2);
        data{9,s}=double(Results(s).Q_Ah_cell); data{10,s}=double(Results(s).dSOC); 
        data{11,s} = double(Results(s).Qmax_cell_Ah);
        data{12,s} = double((Results(s).Qmax_cell_Ah/64)*100);
        data{13,s}=94; % 2023
        % data{13,s}=92.5; % 2025
    end
end
uit = uitable(fig,'Data',data,'ColumnName',cols,'RowName',vars,'Units','normalized','Position',posTbl);
uit.ColumnWidth='auto'; uit.FontName='Consolas';

set(findall(fig,'-property','FontSize'),'FontSize',12);

%% helper
function segs = local_find_segments(mask)
    mask = mask(:)'; d = diff([false mask false]);
    segs = [find(d==1)' (find(d==-1)-1)'];
end


