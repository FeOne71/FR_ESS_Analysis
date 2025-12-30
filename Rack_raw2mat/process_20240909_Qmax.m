%% process_20240909_Qmax.m - Process 20240909 CSV and visualize Qmax
% Converts CSV to MAT format and performs Qmax analysis with visualization

clear all; close all; clc;

%% Parameters
Cnom = 128;                         % Rack nominal Capacity (Ah) (참조용)
C_cell_Ah = 64;                     % Cell capacity (Ah)
thr_A = C_cell_Ah * 0.045;          % Idle threshold (A)  C_cell_Ah * 0.05
Np = 2;                             % parallel cells (2P)
dt = 1;                             % s (고정 가정)

%% Paths
csvFile = '20240909_KIMJ_LGE_01_01_01_KENTECH.csv';
ocvFile = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\OCV_integrated.mat';
base_date = datetime(2024, 09, 09);
min_charge_sec = 300;  % 5 minutes minimum charge duration

% Create Qmax/charge folder for results
if ~exist('Qmax', 'dir')
    mkdir('Qmax');
end
if ~exist('Qmax/charge', 'dir')
    mkdir('Qmax/charge');
end

%% Step 1: Convert CSV to MAT format (similar to New_rack2mat.m)
fprintf('=== Step 1: Converting CSV to MAT format ===\n');

% Read CSV file
T = readtable(csvFile, 'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
if isempty(T)
    error('Empty table for file %s', csvFile);
end

header_row = T.Properties.VariableNames;
wanted_vars_raw = { ...
    'Date_Time', 'Rack_Status', 'Charge_Mode', 'Discharge_Mode', 'CB_Feedback', ...
    'Contactor_Feedback', 'CellBalancing', 'SOC_BMS', 'SOH_BMS', 'DCCurrent', ...
    'DCPower', 'DCchgPowerLimit', 'DCdchgPowerLimit', 'CVavg', 'CVmax', 'CVmin', ...
    'CVmaxNoMBMS', 'CVmaxNoCell', 'CVminNoMBMS', 'CVminNoCell', ...
    'MTavg', 'MTmax', 'MTmin', 'MTmaxNoMBMS', 'MTminNoMBMS'};

selected_cols = [];
for k = 1:length(wanted_vars_raw)
    idx = find(strcmp(strtrim(header_row), wanted_vars_raw{k}), 1);
    if ~isempty(idx)
        selected_cols(end+1) = idx;
    end
end

if isempty(selected_cols)
    error('No matching columns in file %s', csvFile);
end

T_sel = T(:, selected_cols);
clean_var_names = cell(size(selected_cols));
for i = 1:length(selected_cols)
    clean_var_names{i} = convertVariableName(header_row{selected_cols(i)}, false);
end
T_sel.Properties.VariableNames = clean_var_names;

for j = 1:length(clean_var_names)
    var_name = clean_var_names{j};
    data_col = T_sel.(var_name);
    if iscell(data_col)
        numeric_data = nan(size(data_col));
        for k = 1:length(data_col)
            if ~isempty(data_col{k}) && ~isnan(str2double(data_col{k}))
                numeric_data(k) = str2double(data_col{k});
            end
        end
        T_sel.(var_name) = numeric_data;
    elseif isa(data_col, 'duration')
        % duration 타입을 그대로 저장
        if strcmp(var_name, 'DateTime') || strcmp(var_name, 'Date_Time')
            T_sel.Date_Time = data_col;
        else
            T_sel.(var_name) = seconds(data_col);
        end
    else
        T_sel.(var_name) = double(T_sel.(var_name));
    end
end

% Raw 구조체 생성 
Raw = struct();
field_names = T_sel.Properties.VariableNames;
for i = 1:length(field_names)
    field = convertVariableName(field_names{i}, false);
    Raw.(field) = T_sel.(field_names{i});
end

% Save MAT file
mat_file = 'Raw_20240909.mat';
save(mat_file, 'Raw', '-v7.3');
fprintf('MAT file saved: %s\n', mat_file);

%% Step 2: Qmax Analysis (similar to FieldQmax_Newdata.m)
fprintf('\n=== Step 2: Qmax Analysis ===\n');

% Load data
D = Raw;

% Time
if isduration(D.Date_Time)
    t = base_date + D.Date_Time;
else
    t = datetime(D.Date_Time);
end
t0 = t(1); tsec = seconds(t - t0);

% Signals
I_rack = D.DCCurrent(:);
Vcell_avg = D.CVavg(:);
P_rack_kW = D.DCPower(:) / 1000;  % W -> kW

% Cell-level signals
I_cell = I_rack / Np;  % A per cell
V_cell = Vcell_avg;    % V per cell

% SOC (raw from BMS)
SOC_BMS = D.SOC_BMS(:);
SOH_raw = D.SOH_BMS(:);           % Raw SOH from BMS (%)
SOH_end = SOH_raw(end);          % Get the last SOH value

% Load OCV data
T_ocv = load(ocvFile);
OCV_data = T_ocv.OCV_data;

% Build inverse OCV->SOC using avg_ocv_rpt0 and soc_grid (data already clean)
ocv = OCV_data.avg_ocv_rpt0(:);
soc = OCV_data.soc_grid(:);
SOC_from_OCV = @(v) interp1(ocv, soc, v, 'linear', 'extrap');

%% 전압 필터링
window_size_samples = 30; % 30초 윈도우 (1초 샘플링 가정)
Vcell_avg_filtered = movmean(Vcell_avg, window_size_samples);

%% 원본 전류로 충전/방전/유휴 구간 계산
isIdle = abs(I_cell) < thr_A / Np;
isChg = I_cell > thr_A / Np;
isDchg = I_cell < -thr_A / Np;
idleSegs = local_find_segments(isIdle);
chgSegs = local_find_segments(isChg);
dchgSegs = local_find_segments(isDchg);

% Filter charge segments by minimum duration
valid_segs = [];
for k = 1:size(chgSegs,1)
    seg_duration = chgSegs(k,2) - chgSegs(k,1) + 1;
    if seg_duration >= min_charge_sec
        valid_segs(end+1) = k;
    end
end
chgSegs = chgSegs(valid_segs,:);

% Filter discharge segments by minimum duration
valid_segs = [];
for k = 1:size(dchgSegs,1)
    seg_duration = dchgSegs(k,2) - dchgSegs(k,1) + 1;
    if seg_duration >= min_charge_sec
        valid_segs(end+1) = k;
    end
end
dchgSegs = dchgSegs(valid_segs,:);

fprintf('Found %d charge segments (>= %d sec)\n', size(chgSegs,1), min_charge_sec);
fprintf('Found %d discharge segments (>= %d sec)\n', size(dchgSegs,1), min_charge_sec);

% Process each charge segment
Results = [];
rows = {};
result_idx = 0;  % Track actual result index
for k = 1:size(chgSegs,1)
    chg_start = chgSegs(k,1);
    chg_end = chgSegs(k,2);
    
    % Find idle segments around charge
    prevIdleIdx = find(idleSegs(:,2) < chg_start, 1, 'last');
    nextIdleIdx = find(idleSegs(:,1) > chg_end, 1, 'first');
    
    if isempty(prevIdleIdx) || isempty(nextIdleIdx)
        if isempty(prevIdleIdx)
            fprintf('Skipping segment %d: no idle segment before charge\n', k);
        end
        if isempty(nextIdleIdx)
            fprintf('Skipping segment %d: no idle segment after charge\n', k);
        end
        continue;
    end

    % 충전 시작 직전 휴지 구간의 마지막 시점
    befChg = idleSegs(prevIdleIdx, 2);
    continuous_idle_start = idleSegs(prevIdleIdx, 1);
    
    % 실제 연속 휴지구간이 충분히 긴지 확인 (최소 5분)
    actual_rest_duration_sec = befChg - continuous_idle_start + 1;
    if actual_rest_duration_sec < 300
        fprintf('Skipping segment %d: rest duration too short (%.1f min)\n', k, actual_rest_duration_sec/60);
        continue;
    end
    
    % 충전 종료 후부터 실제 idle threshold를 만족하는 연속 구간 찾기
    aftChg = chg_end;
    for i = chg_end+1:idleSegs(nextIdleIdx,2)
        if abs(I_cell(i)) >= thr_A / Np
            break;
        end
        aftChg = i;
    end
    
    % SOC2 휴지구간 최소 시간 확인 (3분 이상)
    soc2_rest_duration = aftChg - chg_end;
    if soc2_rest_duration < 180
        fprintf('Skipping segment %d: SOC2 rest duration too short (%.1f min)\n', k, soc2_rest_duration/60);
        continue;
    end
    
    fprintf('\n=== Segment %d ===\n', k);
    fprintf('Charge: idx %d-%d (%.1f min)\n', chg_start, chg_end, (chg_end-chg_start+1)/60);
    fprintf('SOC1 time: %s\n', datestr(t(befChg), 'yyyy-mm-dd HH:MM:SS'));
    fprintf('SOC2 time: %s\n', datestr(t(aftChg), 'yyyy-mm-dd HH:MM:SS'));
    
    % SOC at idle end points
    V_befChg = Vcell_avg(befChg);
    V_aftChg = Vcell_avg(aftChg);
    I_befChg = I_cell(befChg);
    I_aftChg = I_cell(aftChg);

    % SOC at idle end points (OCV-based and BMS)
    SOC1 = SOC_from_OCV(V_befChg);
    SOC2 = SOC_from_OCV(V_aftChg);
    SOC1_raw = SOC_BMS(befChg);
    SOC2_raw = SOC_BMS(aftChg);

    % integrate cell current
    t_sec = seconds(t - t(1));
    Q_Ah_cell = abs(trapz(t_sec(chg_start:chg_end), I_cell(chg_start:chg_end)) / 3600);

    % SOC change used for Qmax
    dSOC_q = SOC2 - SOC1;
    dSOC_q_BMS = SOC2_raw - SOC1_raw;

    % Qmax calculation
    if ~isnan(dSOC_q) && dSOC_q ~= 0
        Qmax_cell_Ah = Q_Ah_cell / (abs(dSOC_q)/100);
    else
        Qmax_cell_Ah = NaN;
    end
    
    if ~isnan(dSOC_q_BMS) && dSOC_q_BMS ~= 0
        Qmax_cell_Ah_BMS = Q_Ah_cell / (abs(dSOC_q_BMS)/100);
    else
        Qmax_cell_Ah_BMS = NaN;
    end
    
    % Calculate rest period durations
    rest1_start_idx = continuous_idle_start;
    rest1_end_idx = befChg;
    rest1_duration_sec = actual_rest_duration_sec;
    rest1_duration = duration(0, 0, rest1_duration_sec);
    
    rest2_start_idx = idleSegs(nextIdleIdx, 1);
    rest2_end_idx = idleSegs(nextIdleIdx, 2);
    rest2_duration_sec = rest2_end_idx - rest2_start_idx + 1;
    rest2_duration = duration(0, 0, rest2_duration_sec);

    % store (use result_idx to avoid empty elements)
    result_idx = result_idx + 1;
    Results(result_idx).idx = k;
    Results(result_idx).chg_start = chg_start;
    Results(result_idx).chg_end = chg_end;
    Results(result_idx).befChg = befChg;
    Results(result_idx).aftChg = aftChg;
    Results(result_idx).V1 = Vcell_avg(befChg);
    Results(result_idx).V2 = Vcell_avg(aftChg);
    Results(result_idx).I1 = I_cell(befChg);
    Results(result_idx).I2 = I_cell(aftChg);
    Results(result_idx).SOC1 = SOC1;
    Results(result_idx).SOC2 = SOC2;
    Results(result_idx).dSOC = abs(dSOC_q);
    Results(result_idx).SOC1_raw = SOC1_raw;
    Results(result_idx).SOC2_raw = SOC2_raw;
    Results(result_idx).Q_Ah_cell = Q_Ah_cell;
    Results(result_idx).Qmax_cell_Ah = Qmax_cell_Ah;
    Results(result_idx).Qmax_cell_Ah_BMS = Qmax_cell_Ah_BMS;
    Results(result_idx).rest1_duration = rest1_duration;
    Results(result_idx).rest2_duration = rest2_duration;
    Results(result_idx).rest1_start_idx = rest1_start_idx;
    Results(result_idx).rest1_end_idx = rest1_end_idx;
    Results(result_idx).rest2_start_idx = rest2_start_idx;
    Results(result_idx).rest2_end_idx = rest2_end_idx;
    Results(result_idx).continuous_idle_start = continuous_idle_start;
    Results(result_idx).actual_rest_duration_sec = actual_rest_duration_sec;

    dSOC_row = abs(dSOC_q);
    dSOC_row_BMS = abs(dSOC_q_BMS);
    Results(result_idx).type = 'charge';  % Mark as charge segment
    
    dSOC_row = abs(dSOC_q);
    dSOC_row_BMS = abs(dSOC_q_BMS);
    rows(end+1, :) = {k, datestr(t(chg_start)), datestr(t(chg_end)), ...
        char(rest1_duration), char(rest2_duration), SOC1, SOC2, dSOC_row, ...
        SOC1_raw, SOC2_raw, dSOC_row_BMS, Qmax_cell_Ah, Qmax_cell_Ah_BMS, ...
        Q_Ah_cell, Vcell_avg(chg_start), Vcell_avg(chg_end), ...
        I_cell(chg_start), I_cell(chg_end), SOH_end, SOH_end};
end

%% Process Discharge Segments
fprintf('\n=== Processing Discharge Segments ===\n');
dchg_result_start_idx = result_idx;  % Track where discharge segments start

for k = 1:size(dchgSegs,1)
    dchg_start = dchgSegs(k,1);
    dchg_end = dchgSegs(k,2);
    
    % Find idle segments around discharge
    prevIdleIdx = find(idleSegs(:,2) < dchg_start, 1, 'last');
    nextIdleIdx = find(idleSegs(:,1) > dchg_end, 1, 'first');
    
    if isempty(prevIdleIdx) || isempty(nextIdleIdx)
        if isempty(prevIdleIdx)
            fprintf('Skipping discharge segment %d: no idle segment before discharge\n', k);
        end
        if isempty(nextIdleIdx)
            fprintf('Skipping discharge segment %d: no idle segment after discharge\n', k);
        end
        continue;
    end

    % 방전 시작 직전 휴지 구간의 마지막 시점
    befDchg = idleSegs(prevIdleIdx, 2);
    continuous_idle_start = idleSegs(prevIdleIdx, 1);
    
    % 실제 연속 휴지구간이 충분히 긴지 확인 (최소 5분)
    actual_rest_duration_sec = befDchg - continuous_idle_start + 1;
    if actual_rest_duration_sec < 300
        fprintf('Skipping discharge segment %d: rest duration too short (%.1f min)\n', k, actual_rest_duration_sec/60);
        continue;
    end
    
    % 방전 종료 후부터 실제 idle threshold를 만족하는 연속 구간 찾기
    aftDchg = dchg_end;
    for i = dchg_end+1:idleSegs(nextIdleIdx,2)
        if abs(I_cell(i)) >= thr_A / Np
            break;
        end
        aftDchg = i;
    end
    
    % SOC2 휴지구간 최소 시간 확인 (3분 이상)
    soc2_rest_duration = aftDchg - dchg_end;
    if soc2_rest_duration < 180
        fprintf('Skipping discharge segment %d: SOC2 rest duration too short (%.1f min)\n', k, soc2_rest_duration/60);
        continue;
    end
    
    fprintf('\n=== Discharge Segment %d ===\n', k);
    fprintf('Discharge: idx %d-%d (%.1f min)\n', dchg_start, dchg_end, (dchg_end-dchg_start+1)/60);
    fprintf('SOC1 time: %s\n', datestr(t(befDchg), 'yyyy-mm-dd HH:MM:SS'));
    fprintf('SOC2 time: %s\n', datestr(t(aftDchg), 'yyyy-mm-dd HH:MM:SS'));
    
    % SOC at idle end points
    V_befDchg = Vcell_avg(befDchg);
    V_aftDchg = Vcell_avg(aftDchg);
    I_befDchg = I_cell(befDchg);
    I_aftDchg = I_cell(aftDchg);

    % SOC at idle end points (OCV-based and BMS)
    SOC1 = SOC_from_OCV(V_befDchg);
    SOC2 = SOC_from_OCV(V_aftDchg);
    SOC1_raw = SOC_BMS(befDchg);
    SOC2_raw = SOC_BMS(aftDchg);

    % integrate cell current (discharge is negative, so use abs)
    t_sec = seconds(t - t(1));
    Q_Ah_cell = abs(trapz(t_sec(dchg_start:dchg_end), I_cell(dchg_start:dchg_end)) / 3600);

    % SOC change used for Qmax (discharge: SOC1 > SOC2, so dSOC is negative)
    dSOC_q = SOC1 - SOC2;  % For discharge, SOC decreases
    dSOC_q_BMS = SOC1_raw - SOC2_raw;

    % Qmax calculation
    if ~isnan(dSOC_q) && dSOC_q ~= 0
        Qmax_cell_Ah = Q_Ah_cell / (abs(dSOC_q)/100);
    else
        Qmax_cell_Ah = NaN;
    end
    
    if ~isnan(dSOC_q_BMS) && dSOC_q_BMS ~= 0
        Qmax_cell_Ah_BMS = Q_Ah_cell / (abs(dSOC_q_BMS)/100);
    else
        Qmax_cell_Ah_BMS = NaN;
    end
    
    % Calculate rest period durations
    rest1_start_idx = continuous_idle_start;
    rest1_end_idx = befDchg;
    rest1_duration_sec = actual_rest_duration_sec;
    rest1_duration = duration(0, 0, rest1_duration_sec);
    
    rest2_start_idx = idleSegs(nextIdleIdx, 1);
    rest2_end_idx = idleSegs(nextIdleIdx, 2);
    rest2_duration_sec = rest2_end_idx - rest2_start_idx + 1;
    rest2_duration = duration(0, 0, rest2_duration_sec);

    % store (use result_idx to avoid empty elements)
    result_idx = result_idx + 1;
    Results(result_idx).idx = 100 + k;  % Use 100+ for discharge segments
    Results(result_idx).type = 'discharge';  % Mark as discharge segment
    Results(result_idx).chg_start = dchg_start;  % Reuse field name for discharge start
    Results(result_idx).chg_end = dchg_end;  % Reuse field name for discharge end
    Results(result_idx).befChg = befDchg;  % Reuse field name
    Results(result_idx).aftChg = aftDchg;  % Reuse field name
    Results(result_idx).V1 = Vcell_avg(befDchg);
    Results(result_idx).V2 = Vcell_avg(aftDchg);
    Results(result_idx).I1 = I_cell(befDchg);
    Results(result_idx).I2 = I_cell(aftDchg);
    Results(result_idx).SOC1 = SOC1;
    Results(result_idx).SOC2 = SOC2;
    Results(result_idx).dSOC = abs(dSOC_q);
    Results(result_idx).SOC1_raw = SOC1_raw;
    Results(result_idx).SOC2_raw = SOC2_raw;
    Results(result_idx).Q_Ah_cell = Q_Ah_cell;
    Results(result_idx).Qmax_cell_Ah = Qmax_cell_Ah;
    Results(result_idx).Qmax_cell_Ah_BMS = Qmax_cell_Ah_BMS;
    Results(result_idx).rest1_duration = rest1_duration;
    Results(result_idx).rest2_duration = rest2_duration;
    Results(result_idx).rest1_start_idx = rest1_start_idx;
    Results(result_idx).rest1_end_idx = rest1_end_idx;
    Results(result_idx).rest2_start_idx = rest2_start_idx;
    Results(result_idx).rest2_end_idx = rest2_end_idx;
    Results(result_idx).continuous_idle_start = continuous_idle_start;
    Results(result_idx).actual_rest_duration_sec = actual_rest_duration_sec;

    dSOC_row = abs(dSOC_q);
    dSOC_row_BMS = abs(dSOC_q_BMS);
    rows(end+1, :) = {100+k, datestr(t(dchg_start)), datestr(t(dchg_end)), ...
        char(rest1_duration), char(rest2_duration), SOC1, SOC2, dSOC_row, ...
        SOC1_raw, SOC2_raw, dSOC_row_BMS, Qmax_cell_Ah, Qmax_cell_Ah_BMS, ...
        Q_Ah_cell, Vcell_avg(dchg_start), Vcell_avg(dchg_end), ...
        I_cell(dchg_start), I_cell(dchg_end), SOH_end, SOH_end};
end

%% Add Manual Segment 04 (10:42:47 ~ 12:53:46)
fprintf('\n=== Adding Manual Segment 04 ===\n');
target_start_time = datetime(2024, 09, 09, 10, 42, 47);
target_end_time = datetime(2024, 09, 09, 12, 53, 46);

% Find indices for the time range
chg_start_idx = find(t >= target_start_time, 1, 'first');
chg_end_idx = find(t <= target_end_time, 1, 'last');

if isempty(chg_start_idx) || isempty(chg_end_idx) || chg_start_idx >= chg_end_idx
    fprintf('Warning: Could not find valid time range for Segment 04\n');
else
    % Find idle segments around this charge period
    prevIdleIdx = find(idleSegs(:,2) < chg_start_idx, 1, 'last');
    nextIdleIdx = find(idleSegs(:,1) > chg_end_idx, 1, 'first');
    
    if isempty(prevIdleIdx) || isempty(nextIdleIdx)
        fprintf('Warning: Could not find idle segments around Segment 04\n');
    else
        % 충전 시작 직전 휴지 구간의 마지막 시점
        befChg = idleSegs(prevIdleIdx, 2);
        continuous_idle_start = idleSegs(prevIdleIdx, 1);
        
        % 실제 연속 휴지구간이 충분히 긴지 확인 (최소 5분)
        actual_rest_duration_sec = befChg - continuous_idle_start + 1;
        if actual_rest_duration_sec < 300
            fprintf('Warning: Rest duration too short for Segment 04 (%.1f min)\n', actual_rest_duration_sec/60);
        else
            % 충전 종료 후부터 실제 idle threshold를 만족하는 연속 구간 찾기
            aftChg = chg_end_idx;
            for i = chg_end_idx+1:idleSegs(nextIdleIdx,2)
                if abs(I_cell(i)) >= thr_A / Np
                    break;
                end
                aftChg = i;
            end
            
            % SOC2 휴지구간 최소 시간 확인 (3분 이상)
            soc2_rest_duration = aftChg - chg_end_idx;
            if soc2_rest_duration < 180
                fprintf('Warning: SOC2 rest duration too short for Segment 04 (%.1f min)\n', soc2_rest_duration/60);
            else
                fprintf('Segment 04: idx %d-%d (%.1f min)\n', chg_start_idx, chg_end_idx, (chg_end_idx-chg_start_idx+1)/60);
                fprintf('SOC1 time: %s\n', datestr(t(befChg), 'yyyy-mm-dd HH:MM:SS'));
                fprintf('SOC2 time: %s\n', datestr(t(aftChg), 'yyyy-mm-dd HH:MM:SS'));
                
                % SOC at idle end points
                V_befChg = Vcell_avg(befChg);
                V_aftChg = Vcell_avg(aftChg);
                I_befChg = I_cell(befChg);
                I_aftChg = I_cell(aftChg);
                
                % SOC at idle end points (OCV-based and BMS)
                SOC1 = SOC_from_OCV(V_befChg);
                SOC2 = SOC_from_OCV(V_aftChg);
                SOC1_raw = SOC_BMS(befChg);
                SOC2_raw = SOC_BMS(aftChg);
                
                % integrate cell current
                t_sec = seconds(t - t(1));
                Q_Ah_cell = abs(trapz(t_sec(chg_start_idx:chg_end_idx), I_cell(chg_start_idx:chg_end_idx)) / 3600);
                
                % SOC change used for Qmax
                dSOC_q = SOC2 - SOC1;
                dSOC_q_BMS = SOC2_raw - SOC1_raw;
                
                % Qmax calculation
                if ~isnan(dSOC_q) && dSOC_q ~= 0
                    Qmax_cell_Ah = Q_Ah_cell / (abs(dSOC_q)/100);
                else
                    Qmax_cell_Ah = NaN;
                end
                
                if ~isnan(dSOC_q_BMS) && dSOC_q_BMS ~= 0
                    Qmax_cell_Ah_BMS = Q_Ah_cell / (abs(dSOC_q_BMS)/100);
                else
                    Qmax_cell_Ah_BMS = NaN;
                end
                
                % Calculate rest period durations
                rest1_start_idx = continuous_idle_start;
                rest1_end_idx = befChg;
                rest1_duration_sec = actual_rest_duration_sec;
                rest1_duration = duration(0, 0, rest1_duration_sec);
                
                rest2_start_idx = idleSegs(nextIdleIdx, 1);
                rest2_end_idx = idleSegs(nextIdleIdx, 2);
                rest2_duration_sec = rest2_end_idx - rest2_start_idx + 1;
                rest2_duration = duration(0, 0, rest2_duration_sec);
                
                % store as Segment 04
                result_idx = result_idx + 1;
                Results(result_idx).idx = 4;  % Manual segment 04
                Results(result_idx).type = 'charge';  % Mark as charge segment
                Results(result_idx).chg_start = chg_start_idx;
                Results(result_idx).chg_end = chg_end_idx;
                Results(result_idx).befChg = befChg;
                Results(result_idx).aftChg = aftChg;
                Results(result_idx).V1 = Vcell_avg(befChg);
                Results(result_idx).V2 = Vcell_avg(aftChg);
                Results(result_idx).I1 = I_cell(befChg);
                Results(result_idx).I2 = I_cell(aftChg);
                Results(result_idx).SOC1 = SOC1;
                Results(result_idx).SOC2 = SOC2;
                Results(result_idx).dSOC = abs(dSOC_q);
                Results(result_idx).SOC1_raw = SOC1_raw;
                Results(result_idx).SOC2_raw = SOC2_raw;
                Results(result_idx).Q_Ah_cell = Q_Ah_cell;
                Results(result_idx).Qmax_cell_Ah = Qmax_cell_Ah;
                Results(result_idx).Qmax_cell_Ah_BMS = Qmax_cell_Ah_BMS;
                Results(result_idx).rest1_duration = rest1_duration;
                Results(result_idx).rest2_duration = rest2_duration;
                Results(result_idx).rest1_start_idx = rest1_start_idx;
                Results(result_idx).rest1_end_idx = rest1_end_idx;
                Results(result_idx).rest2_start_idx = rest2_start_idx;
                Results(result_idx).rest2_end_idx = rest2_end_idx;
                Results(result_idx).continuous_idle_start = continuous_idle_start;
                Results(result_idx).actual_rest_duration_sec = actual_rest_duration_sec;
                
                dSOC_row = abs(dSOC_q);
                dSOC_row_BMS = abs(dSOC_q_BMS);
                rows(end+1, :) = {4, datestr(t(chg_start_idx)), datestr(t(chg_end_idx)), ...
                    char(rest1_duration), char(rest2_duration), SOC1, SOC2, dSOC_row, ...
                    SOC1_raw, SOC2_raw, dSOC_row_BMS, Qmax_cell_Ah, Qmax_cell_Ah_BMS, ...
                    Q_Ah_cell, Vcell_avg(chg_start_idx), Vcell_avg(chg_end_idx), ...
                    I_cell(chg_start_idx), I_cell(chg_end_idx), SOH_end, SOH_end};
                
                fprintf('Segment 04 added successfully\n');
            end
        end
    end
end

% Sort Results by idx to ensure proper ordering (charge first, then discharge)
if ~isempty(Results)
    % Separate charge and discharge
    charge_mask = strcmp({Results.type}, 'charge');
    discharge_mask = strcmp({Results.type}, 'discharge');
    
    charge_results = Results(charge_mask);
    discharge_results = Results(discharge_mask);
    
    % Sort charge by idx
    if ~isempty(charge_results)
        [~, sort_idx] = sort([charge_results.idx]);
        charge_results = charge_results(sort_idx);
    end
    
    % Sort discharge by idx
    if ~isempty(discharge_results)
        [~, sort_idx] = sort([discharge_results.idx]);
        discharge_results = discharge_results(sort_idx);
    end
    
    % Combine: charge first, then discharge
    Results = [charge_results, discharge_results];
end

%% Visualization
fprintf('\n=== Step 3: Visualization ===\n');

year = '2024';
fig = figure('Name', sprintf('FieldQmax Rack01 (%s)', year), 'NumberTitle', 'off');
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (1,1) Current with SOC overlay
ax1 = nexttile(tl, 1); hold(ax1, 'on'); grid(ax1, 'on');
plot(ax1, t, I_cell, '-', 'Color', [0.8 0.3 0.1], 'LineWidth', 2, 'MarkerSize', 2);

if ~isempty(Results)
    min_time = t(min([Results.befChg])) - hours(1);
    max_time = t(max([Results.aftChg])) + hours(1);
    xlim(ax1, [min_time, max_time]);
end

title(ax1, 'Cell Current I_{cell} [A]');
xlabel(ax1, 'Time');
ylabel(ax1, 'I_{cell} [A]');
ylim([-200 150]);

% Add vertical dashed lines for SOC1, SOC2
for k = 1:numel(Results)
    befChg = Results(k).befChg;
    aftChg = Results(k).aftChg;
    xline(ax1, t(befChg), '--', 'Color', [0 0 1], 'LineWidth', 2, 'Alpha', 0.7);
    xline(ax1, t(aftChg), '--', 'Color', [0 0 1], 'LineWidth', 2, 'Alpha', 0.7);
end

% SOC overlay
axes(ax1); yyaxis right; ylabel('SOC [%]');
for k = 1:numel(Results)
    cstart = Results(k).chg_start;
    cend = Results(k).chg_end;
    befChg = Results(k).befChg;
    aftChg = Results(k).aftChg;
    SOC1 = Results(k).SOC1;
    SOC2 = Results(k).SOC2;
    
    plot([t(befChg) t(cstart)], [SOC1 SOC1], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 2);
    plot([t(cstart) t(cend)], [SOC1 SOC2], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 2);
    plot([t(cend) t(aftChg)], [SOC2 SOC2], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 2);
    plot(t(befChg), SOC1, 'o', 'Color', [0 0 1], 'MarkerFaceColor', [0 0 1], 'MarkerSize', 2);
    plot(t(aftChg), SOC2, 'o', 'Color', [0 0 1], 'MarkerFaceColor', [0 0 1], 'MarkerSize', 2);
end
yyaxis left;
ax1.XColor = [0 0 0];
if numel(ax1.YAxis) >= 1, ax1.YAxis(1).Color = [0 0 0]; end
if numel(ax1.YAxis) >= 2, ax1.YAxis(2).Color = [0 0 0]; end

% (2,1) Voltage with SOC overlay
ax2 = nexttile(tl, 3); hold(ax2, 'on'); grid(ax2, 'on');
plot(ax2, t, Vcell_avg, '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 2, 'MarkerSize', 2);

if ~isempty(Results)
    min_time = t(min([Results.befChg])) - hours(1);
    max_time = t(max([Results.aftChg])) + hours(1);
    xlim(ax2, [min_time, max_time]);
end

title(ax2, 'Average Cell Voltage V_{cell,avg} [V]');
xlabel(ax2, 'Time');
ylabel(ax2, 'V_{cell,avg} [V]');
ylim([2.7 4.3]);

% Add vertical dashed lines
for k = 1:numel(Results)
    befChg = Results(k).befChg;
    aftChg = Results(k).aftChg;
    xline(ax2, t(befChg), '--', 'Color', [0 0 1], 'LineWidth', 2, 'Alpha', 0.7);
    xline(ax2, t(aftChg), '--', 'Color', [0 0 1], 'LineWidth', 2, 'Alpha', 0.7);
end

% SOC overlay
axes(ax2); yyaxis right; ylabel('SOC [%]');
for k = 1:numel(Results)
    cstart = Results(k).chg_start;
    cend = Results(k).chg_end;
    befChg = Results(k).befChg;
    aftChg = Results(k).aftChg;
    SOC1 = Results(k).SOC1;
    SOC2 = Results(k).SOC2;
    
    plot([t(befChg) t(cstart)], [SOC1 SOC1], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 2);
    plot([t(cstart) t(cend)], [SOC1 SOC2], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 2);
    plot([t(cend) t(aftChg)], [SOC2 SOC2], '-', 'Color', [0.45 0.7 0.2], 'LineWidth', 2);
    plot(t(befChg), SOC1, 'o', 'Color', [0 0 1], 'MarkerFaceColor', [0 0 1], 'MarkerSize', 2);
    plot(t(aftChg), SOC2, 'o', 'Color', [0 0 1], 'MarkerFaceColor', [0 0 1], 'MarkerSize', 2);
end
yyaxis left;
ax2.XColor = [0 0 0];
if numel(ax2.YAxis) >= 1, ax2.YAxis(1).Color = [0 0 0]; end
if numel(ax2.YAxis) >= 2, ax2.YAxis(2).Color = [0 0 0]; end

% Right column: summary table
axTbl = nexttile(tl, 2, [2 1]);
posTbl = axTbl.Position; delete(axTbl);
varNames = {'Start', 'End', 'Rest1 [HH:MM:SS]', 'Rest2 [HH:MM:SS]', 'SOC1(EST)[%]', ...
    'SOC2(EST)[%]', 'ΔSOC(EST)[%]', 'SOC1(Raw)[%]', 'SOC2(Raw)[%]', 'ΔSOC(Raw)[%]', ...
    'Qmax(EST) [Ah]', 'Qmax(BMS) [Ah]', '∫I dt [Ah]', 'OCV1[V]', 'OCV2[V]', ...
    'I1[A]', 'I2[A]', 'SOH(Raw) [%]', 'SOH(EST) [%]', 'SOH(BMS) [%]'};
numSeg = numel(Results);
tblData = cell(numel(varNames), max(numSeg, 1));
colNames = cell(1, max(numSeg, 1));
if numSeg == 0
    colNames{1} = 'Segment01';
else
    for s = 1:numSeg
        % Use actual idx and type from Results
        if isfield(Results(s), 'type') && strcmp(Results(s).type, 'discharge')
            colNames{s} = sprintf('Dchg%02d', Results(s).idx - 100);
        else
            colNames{s} = sprintf('Chg%02d', Results(s).idx);
        end
        tblData{1,s} = datestr(t(Results(s).befChg));
        tblData{2,s} = datestr(t(Results(s).aftChg));
        tblData{3,s} = char(Results(s).rest1_duration);
        tblData{4,s} = char(Results(s).rest2_duration);
        tblData{5,s} = sprintf('%.4f', Results(s).SOC1);
        tblData{6,s} = sprintf('%.4f', Results(s).SOC2);
        tblData{7,s} = sprintf('%.4f', Results(s).dSOC);
        tblData{8,s} = sprintf('%.4f', Results(s).SOC1_raw);
        tblData{9,s} = sprintf('%.4f', Results(s).SOC2_raw);
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
uit = uitable(fig, 'Data', tblData, 'ColumnName', colNames, 'RowName', varNames, ...
    'Units', 'normalized', 'Position', posTbl);
uit.ColumnWidth = {300};
uit.ColumnFormat = {'char', 'char', 'char', 'char', 'numeric', 'numeric', 'numeric', ...
    'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', ...
    'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'};
set(findall(fig, '-property', 'FontSize'), 'FontSize', 12);
set(findall(fig, '-property', 'FontWeight'), 'FontWeight', 'bold');

% Save figure
saveas(fig, sprintf('Qmax/charge/FieldQmax_Figure_%s.fig', year));

% Save table data to mat file
T = array2table(tblData', 'VariableNames', varNames, 'RowNames', colNames);
save(sprintf('Qmax/charge/FieldQmax_Results_%s.mat', year), 'Results', 'T', 'tblData', 'varNames', 'colNames');

fprintf('Visualization saved: Qmax/charge/FieldQmax_Figure_%s.fig\n', year);
fprintf('Results saved: Qmax/charge/FieldQmax_Results_%s.mat\n', year);

%% Helper Functions

function clean_name = convertVariableName(original_name, is_total)
    clean_name = original_name;
    clean_name = strrep(clean_name, ' ', '');
    clean_name = strrep(clean_name, '.', '');
    clean_name = strrep(clean_name, '(%)', 'Pct');
    clean_name = strrep(clean_name, '(A)', '_A');
    clean_name = strrep(clean_name, '(V)', '_V');
    clean_name = strrep(clean_name, '(kW)', '_kW');
    clean_name = strrep(clean_name, '(Wh)', '_Wh');
    clean_name = strrep(clean_name, '(mOhm)', '_mOhm');
    clean_name = strrep(clean_name, '(mA)', '_mA');
    clean_name = strrep(clean_name, '(oC)', '_degC');
    clean_name = regexprep(clean_name, '[^a-zA-Z0-9_]', '');
    if is_total
        clean_name = ['Total_' clean_name];
    end
    if isempty(clean_name) || isstrprop(clean_name(1), 'digit')
        clean_name = ['Var_' clean_name];
    end
end

function segs = local_find_segments(mask)
    mask = mask(:)';
    d = diff([false mask false]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    segs = [starts(:) ends(:)];
end

