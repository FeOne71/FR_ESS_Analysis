% Phase1_Data_Parsing.m
% Parses raw field data and segments into Charge/Discharge/Idle events
% Applies 30-second minimum duration filter
% Outputs event structures to avoid reprocessing large raw files
% EXCLUDES the 4 RPT dates and runs on the rest of the field data
% Np = 2, Not yet implemented. Still using rack current.

clear; clc; close all;

%% Configuration
fieldProcessDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Field_Process';
if ~exist(fieldProcessDir, 'dir')
    mkdir(fieldProcessDir);
end

% Directory to save parsed events
outputDir = fullfile(fieldProcessDir, 'Parsed_Events');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Base directories for Raw Data
oldDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old';
newDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New';

% Find all Raw_*.mat files
oldFiles = dir(fullfile(oldDir, '**', 'Raw_*.mat'));
newFiles = dir(fullfile(newDir, '**', 'Raw_*.mat'));
allRawFiles = [oldFiles; newFiles];

% RPT Dates to Exclude
% rpt_dates = {'20210603', '20231016', '20240909', '20250711'};
rpt_dates = {'20231016','20240909','20250711'};

dataFiles = {};
dates = {};

for i = 1:length(allRawFiles)
    fname = allRawFiles(i).name;
    % Extract date string from filename, e.g., Raw_20210701.mat -> 20210701
    token = regexp(fname, 'Raw_(\d{8})\.mat', 'tokens', 'once');
    if ~isempty(token)
        d_str = token{1};
        if ~ismember(d_str, rpt_dates)
            dataFiles{end+1} = fullfile(allRawFiles(i).folder, fname);
            dates{end+1} = d_str;
        end
    end
end

disp(['Total non-RPT field data files found: ', num2str(length(dataFiles))]);

% Parameters
I_threshold = 0.5;   % [A] Threshold for charge/discharge
min_duration_sec = 30; % [s] Minimum duration for an event
dt = 1;              % [s] assumed sampling time

%% Process each file
for i = 1:length(dataFiles)
    inFile = dataFiles{i};
    date_str = dates{i};
    
    fprintf('Processing %s ...\n', inFile);
    if ~exist(inFile, 'file')
        fprintf('  File not found! Skipping.\n');
        continue;
    end
    
    S = load(inFile);
    
    % Handle Old vs New data formats
    is_old = false;
    if isfield(S, 'Raw')
        if isstruct(S.Raw) && isfield(S.Raw, 'Rack01')
            D = S.Raw.Rack01;
            is_old = true;
        else
            D = S.Raw;
        end
    else
        fprintf('  Unknown data format. Skipping.\n');
        continue;
    end
    
    % Base date estimation from date_str (e.g. 20210701)
    yyyy = str2double(date_str(1:4));
    mm = str2double(date_str(5:6));
    dd = str2double(date_str(7:8));
    base_date = datetime(yyyy, mm, dd);
    
    % Time extraction
    if isfield(D, 'Time')
        t = datetime(D.Time);
    elseif isfield(D, 'Date_Time')
        if isduration(D.Date_Time)
            t = base_date + D.Date_Time;
        else
            t = datetime(D.Date_Time);
        end
    else
        fprintf('  Time field not found. Skipping.\n');
        continue;
    end
    t_sec = seconds(t - t(1));
    
    % Signal extraction 
    if is_old
        I_rack = D.DCCurrent_A(:);
        V_avg = D.AverageCV_V(:);
        SOC_bms = D.SOCPct(:);
        if isfield(D, 'AverageMT_C')
            T_avg = D.AverageMT_C(:);
        else
            T_avg = zeros(size(I_rack)); % Dummy if not available
        end
        if isfield(D, 'SOHPct')
            SOH_bms = D.SOHPct(:);
        else
            SOH_bms = NaN(size(I_rack));
        end
    else
        I_rack = D.DCCurrent(:);
        V_avg = D.CVavg(:);
        SOC_bms = D.SOC_BMS(:);
        if isfield(D, 'MTavg')
            T_avg = D.MTavg(:);
        else
            T_avg = zeros(size(I_rack));
        end
        if isfield(D, 'SOH_BMS')
            SOH_bms = D.SOH_BMS(:);
        else
            SOH_bms = NaN(size(I_rack));
        end
    end
    
    % --- Event Detection Logic (from EventDetection.m) ---
    idle_thr_rack = 64 * 0.02 * 2; % 2.56A (Assuming Np=2 for rack current)
    min_duration_sec = 30;
    min_idle_duration_sec = 60;
    movmean_window = 5;
    
    % Smooth data
    I_smooth = movmean(I_rack, movmean_window);
    
    % Find events using EventDetection.m logic
    [chgSegs, dchgSegs] = find_events_simplified(I_smooth, idle_thr_rack, min_duration_sec, t_sec, min_idle_duration_sec);
    
    fprintf('  Found %d Charge events (from EventDetection.m logic)\n', size(chgSegs,1));
    fprintf('  Found %d Discharge events (from EventDetection.m logic)\n', size(dchgSegs,1));
    
    % Pack into structures
    events.Charge = pack_segments(chgSegs, t, t_sec, V_avg, I_rack, T_avg, SOC_bms, SOH_bms);
    events.Discharge = pack_segments(dchgSegs, t, t_sec, V_avg, I_rack, T_avg, SOC_bms, SOH_bms);
    events.Idle = struct(); % Idle skipped natively in this logic
    
    % Save parsed events
    outFile = fullfile(outputDir, sprintf('Parsed_Events_%s.mat', date_str));
    save(outFile, 'events', '-v7.3');
    fprintf('  Saved to %s\n', outFile);
end

fprintf('\nPhase 1 Data Parsing Complete.\n');

%% Helper Functions
function [chg_events, dch_events] = find_events_simplified(I_smooth, idle_thr, min_duration, t_sec, min_idle_duration_sec)
    % 이벤트 검출: idle 상태의 마지막 시점부터 이벤트로 저장
    % idle_threshold 이상이면서 동일한 부호로 지속되는 경우만 검출
    % 추가 조건: 이벤트 시작 전 휴지 구간이 최소 min_idle_duration_sec 이상이어야 함
    chg_events = [];
    dch_events = [];
    n = length(I_smooth);
    i = 1;
    
    while i <= n
        % 휴지 상태에서 벗어난 경우 (|I| >= idle_thr)
        if abs(I_smooth(i)) >= idle_thr
            event_start = i;  
            idle_start_idx = [];  
            
            % i 이전에서 idle 상태(|I| < idle_thr)의 마지막 지점 찾기
            for k = i-1:-1:1
                if abs(I_smooth(k)) < idle_thr
                    event_start = k;  
                    idle_start_idx = 1;  
                    for m = event_start-1:-1:1
                        if abs(I_smooth(m)) >= idle_thr
                            idle_start_idx = m + 1;  
                            break;
                        end
                    end
                    break;
                end
            end
            
            % 만약 처음부터 idle_thr 이상이면 event_start = 1이고 휴지 구간 없음
            if event_start == 1 || isempty(idle_start_idx)
                i = i + 1;
                continue;
            end
            
            % 이벤트 시작 전 휴지 구간 지속 시간 확인
            if idle_start_idx >= event_start
                i = i + 1;
                continue;
            end
            
            % 휴지 구간 지속 시간 계산 (초 단위)
            idle_duration_sec = t_sec(event_start) - t_sec(idle_start_idx);
            
            % 휴지 구간이 최소 시간 미만이면 이벤트 제외
            if idle_duration_sec < min_idle_duration_sec
                i = i + 1;
                continue;
            end
            
            event_type = sign(I_smooth(i)); % 1: charge, -1: discharge
            
            % 동일한 부호로 지속되는 동안 이벤트 끝 찾기
            j = i + 1;
            while j <= n
                if (event_type > 0 && I_smooth(j) >= idle_thr) || (event_type < 0 && I_smooth(j) <= -idle_thr)
                    j = j + 1;
                else
                    break;
                end
            end
            event_end = j - 1;
            
            % 최소 지속 시간 확인
            duration_samples = event_end - event_start + 1;
            if duration_samples >= min_duration
                if event_type > 0
                    chg_events = [chg_events; event_start, event_end];
                else
                    dch_events = [dch_events; event_start, event_end];
                end
            end
            i = event_end + 1; % 다음 탐색 위치로 점프
        else
            i = i + 1;
        end
    end
end

function packed = pack_segments(segs, t, t_sec, V, I, T, SOC, SOH)
    packed = struct();
    for k = 1:size(segs, 1)
        s_idx = segs(k, 1);
        e_idx = segs(k, 2);
        
        evt = struct();
        evt.indices = [s_idx, e_idx];
        evt.t = t(s_idx:e_idx);
        evt.t_sec = t_sec(s_idx:e_idx);
        evt.V = V(s_idx:e_idx);
        evt.I = I(s_idx:e_idx);
        evt.T = T(s_idx:e_idx);
        evt.SOC = SOC(s_idx:e_idx);
        evt.SOH = SOH(s_idx:e_idx);
        
        evtName = sprintf('Evt%04d', k);
        packed.(evtName) = evt;
    end
end
