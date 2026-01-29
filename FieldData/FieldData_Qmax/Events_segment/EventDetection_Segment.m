%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EventDetection_Segment.m
% 이벤트 검출 및 세그먼트 추출 (필드 데이터용)
% 
% 목적: 
% - Raw 데이터에서 충/방전 이벤트 구간 검출
% - 스무딩: 윈도우 크기 3 (인덱스 3개)로 원시 전압, 전류 스무딩
% - 이벤트 검출: Simple_EventDetection처럼 idle -> active 전환 감지
% - 필터링: duration 300초 이상, 3초 ~ end-3초 구간에서 표준편차 기반 필터링
% - 이벤트 시작점: idle -> Active에서 idle 마지막 값
%
% 출력:
% - Events_segment_*.mat (각 날짜별 이벤트 세그먼트 데이터)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Parameters
% =========================================================================
C_cell_Ah = 64;                     % Cell capacity (Ah)
idle_thr_cell = C_cell_Ah * 0.05;  % Idle threshold (A) - 3.2A
Np = 2;                             % parallel cells (2P)
dt = 1;                             % s (고정 가정)

% 이벤트 검출 파라미터
min_duration_sec = 100;             % 최소 이벤트 지속 시간 (초)
smoothing_window = 10;                % 스무딩 윈도우 크기 (인덱스 3개)
std_calc_trim_sec = 10;              % 표준편차 계산 구간: 시작+3초 ~ 끝-3초
max_I_std = 5.0;                    % 최대 전류 표준편차 임계값 (A)

% 데이터 디렉토리 (EventDetection.m과 동일)
scriptDir = fileparts(mfilename('fullpath'));
% 프로젝트 루트로 이동 (Events_segment -> FieldData_Qmax -> FieldData -> KEPCO_ESS_Local -> Rack_raw2mat)
% EventDetection.m은 Reference 폴더에 있어서 '..', '..'이지만, Events_segment는 '..', '..', '..'이 필요
projectRoot = fullfile(scriptDir, '..', '..', '..');
dataDir = fullfile(projectRoot, 'Rack_raw2mat');

% RPT 평가 날짜 제외 (이 날짜들은 분석에서 제외)
RPT_dates = {'2021-06-03', '2023-10-16', '2024-09-09', '2025-07-11'};
RPT_dateSet = containers.Map();
for i = 1:length(RPT_dates)
    RPT_dateSet(RPT_dates{i}) = true;
end

% 연도 및 데이터 타입
yearList = {'2021', '2022', '2023', '2024', '2025'};
dataTypes = {'New', 'Old'};
rackNames_all = {'Rack01'};

% 출력 디렉토리 (현재 폴더의 Results)
saveDir = fullfile(pwd, 'Results');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% 월별 계층적 구조로 결과 변수 초기화 (충전/방전 분리)
all_events = struct();
for r_idx = 1:length(rackNames_all)
    rackName = rackNames_all{r_idx};
    all_events.(rackName) = struct();
    for y_idx = 1:length(yearList)
        year_key = ['Y' yearList{y_idx}];
        all_events.(rackName).(year_key) = struct();
    end
end

fprintf('=== Event Detection and Segment Extraction ===\n');
fprintf('Data directory: %s\n', dataDir);
fprintf('Smoothing window: %d (indices)\n', smoothing_window);
fprintf('Idle threshold: %.2f A\n', idle_thr_cell);
fprintf('Min event duration: %d sec\n', min_duration_sec);
fprintf('Max I std (3s ~ end-3s): %.2f A\n', max_I_std);
fprintf('Excluding RPT dates: %s\n', strjoin(RPT_dates, ', '));
fprintf('\n');

%% Process all files (excluding RPT dates)
for type_idx = 1:length(dataTypes)
    type = dataTypes{type_idx};
    typePath = fullfile(dataDir, type);
    if ~exist(typePath, 'dir')
        continue;
    end
    
    fprintf('\n################### Processing Data Type: %s ###################\n', type);
    
    for year_idx = 1:length(yearList)
        year = yearList{year_idx};
        yearPath = fullfile(typePath, year);
        if ~exist(yearPath, 'dir')
            continue;
        end
        
        fprintf('\n================== Processing Year: %s ==================\n', year);
        
        monthDirs = dir(yearPath);
        monthDirs = monthDirs([monthDirs.isdir] & ~ismember({monthDirs.name}, {'.', '..'}));
        
        for m = 1:length(monthDirs)
            month = monthDirs(m).name;
            if isempty(regexp(month, '^\d{6}$', 'once'))
                continue;
            end
            
            monthPath = fullfile(yearPath, month);
            matFiles = dir(fullfile(monthPath, '*.mat'));
            if length(matFiles) > 0
                fprintf('\n--- Processing Month: %s (%d files) ---\n', month, length(matFiles));
            end
            
            for f = 1:length(matFiles)
                matFilePath = fullfile(monthPath, matFiles(f).name);
                
                % 파일명에서 날짜 추출
                [~, filename_only, ~] = fileparts(matFiles(f).name);
                dateMatch = regexp(filename_only, '\d{8}', 'once', 'match');
                if isempty(dateMatch)
                    fprintf('  Skipping %s (no date in filename)\n', matFiles(f).name);
                    continue;
                end
                
                % 날짜 문자열 생성 (YYYY-MM-DD)
                fileDate_str = sprintf('%s-%s-%s', dateMatch(1:4), dateMatch(5:6), dateMatch(7:8));
                
                % RPT 날짜인지 확인
                if isKey(RPT_dateSet, fileDate_str)
                    fprintf('  Skipping %s (RPT date: %s)\n', matFiles(f).name, fileDate_str);
                    continue;
                end
                
                fprintf('  Processing: %s (date: %s)\n', matFiles(f).name, fileDate_str);
                
                %% Load data
                try
                    S = load(matFilePath);
                catch ME
                    fprintf('  Error loading %s: %s\n', matFiles(f).name, ME.message);
                    continue;
                end
                
                % 데이터 타입 결정
                if strcmpi(type, 'Old')
                    dataType = 'old';
                else
                    dataType = 'new';
                end
                
                % base_date는 파일명에서 추출한 날짜 사용
                fileDate = datetime(dateMatch, 'InputFormat', 'yyyyMMdd');
                base_date = fileDate;
                
                if strcmp(dataType, 'old')
                    % Old data format
                    if ~isfield(S, 'Raw')
                        fprintf('  Error: Raw not found in %s\n', matFiles(f).name);
                        continue;
                    end
                    Raw_Rack = S.Raw;
                    if ~isfield(Raw_Rack, 'Rack01')
                        fprintf('  Error: Raw_Rack.Rack01 not found in %s\n', matFiles(f).name);
                        continue;
                    end
                    rackData = Raw_Rack.Rack01;
        
                    % Time 처리 (EventDetection.m 방식)
                    if isdatetime(rackData.Time)
                        t = rackData.Time;
                    elseif iscell(rackData.Time)
                        % 셀 배열인 경우 첫 번째 요소 확인 후 변환
                        if isdatetime(rackData.Time{1})
                            t = [rackData.Time{:}];
                        elseif isnumeric(rackData.Time{1})
                            t = datetime(cell2mat(rackData.Time), 'ConvertFrom', 'datenum');
                        else
                            t = datetime(rackData.Time);
                        end
                    elseif isnumeric(rackData.Time)
                        % 숫자 형식인 경우 datenum으로 변환
                        t = datetime(rackData.Time, 'ConvertFrom', 'datenum');
                    else
                        % 기타 형식은 datetime 생성자로 변환
                        t = datetime(rackData.Time);
                    end
                    
                    % Signals
                    I_rack = rackData.DCCurrent_A(:);
                    Vcell_avg = rackData.AverageCV_V(:);
                    SOC_raw = rackData.SOCPct(:);  % Raw SOC from BMS (%)
                    P_rack_kW = rackData.DCPower_kW(:);
                    T = rackData.AverageMT_degC(:);
                    
                else
                    % New data format
                    rackData = S.Raw;
                    
                    % Time 처리 (EventDetection.m 방식)
                    if isduration(rackData.Date_Time)
                        if ~isempty(fileDate)
                            % duration을 datetime으로 변환 (파일명 날짜 + duration)
                            t = fileDate + rackData.Date_Time;
                        else
                            % 파일명 날짜가 없으면 base_date 사용
                            t = base_date + rackData.Date_Time;
                        end
                    elseif isdatetime(rackData.Date_Time)
                        t = rackData.Date_Time;
                        % 날짜가 없고 시간만 있는 경우 파일명 날짜와 결합
                        if ~isempty(fileDate) && any(t.Year == 0)
                            t = fileDate + timeofday(t);
                        end
                    else
                        % 기타 형식은 원본 그대로 사용 (이미 datetime일 수 있음)
                        t = rackData.Date_Time;
                    end
                    
                    % 최종 확인: t가 datetime인지 확인
                    if ~isdatetime(t)
                        fprintf('  Error: t is not datetime after processing. Type: %s\n', class(t));
                        continue;
                    end
                    
                    % Signals
                    I_rack = rackData.DCCurrent(:);
                    Vcell_avg = rackData.CVavg(:);
                    SOC_raw = rackData.SOC_BMS(:);  % Raw SOC from BMS (%)
                    P_rack_kW = rackData.DCPower(:);
                    T = rackData.MTavg(:);
                end
                
                % 전력 단위 변환 및 셀 단위 변환
                P_rack_W = P_rack_kW * 1000;
                Ns = 14 * 17;
                P_cell_W = P_rack_W / (Ns * Np);
                
                % 데이터 길이 확인 및 정렬
                vars = {I_rack, Vcell_avg, SOC_raw, P_cell_W, T};
                min_len = min([length(t), cellfun(@length, vars)]);
                if min_len < min_duration_sec
                    fprintf('  Skipping %s (data too short: %d samples)\n', matFiles(f).name, min_len);
                    continue;
                end
                
                t = t(1:min_len);
                I_rack = I_rack(1:min_len);
                Vcell_avg = Vcell_avg(1:min_len);
                SOC_raw = SOC_raw(1:min_len);
                P_cell_W = P_cell_W(1:min_len);
                T = T(1:min_len);
                
                % Cell-level signals
                I_cell_raw = I_rack / Np;  % A per cell
                
                % 1. 스무딩: 윈도우 크기 3 (인덱스 3개)로 원시 전압, 전류 스무딩
                if length(I_cell_raw) >= smoothing_window
                    I_cell_smooth = movmean(I_cell_raw, smoothing_window);
                    V_cell_smooth = movmean(Vcell_avg, smoothing_window);
                else
                    I_cell_smooth = I_cell_raw;
                    V_cell_smooth = Vcell_avg;
                end
                
                % 2. 이벤트 검출: Simple_EventDetection처럼 idle -> active 전환 감지
                % 이벤트 검출에는 smooth 데이터 사용
                [charge_events, discharge_events] = detect_events_idle_to_active(...
                    I_cell_smooth, idle_thr_cell, min_duration_sec, t, ...
                    std_calc_trim_sec, max_I_std);
                
                fprintf('    Found %d charge events, %d discharge events\n', ...
                    size(charge_events, 1), size(discharge_events, 1));
                
                %% Save events (충전/방전 분리, 구조체 배열로 저장)
                % year_key와 month_key 생성
                year_key = ['Y' year];
                month_key = ['M' month(end-1:end)];
                
                % file_info 생성
                file_info = struct('year', year, 'dataType', type, 'filename', matFiles(f).name);
                
                % 이벤트 저장: 충전/방전 분리, 구조체 배열로 저장
                rackName = rackNames_all{1};
                
                % 월별 구조체 초기화 (charge, discharge 구조체)
                if ~isfield(all_events.(rackName).(year_key), month_key)
                    all_events.(rackName).(year_key).(month_key) = struct();
                    all_events.(rackName).(year_key).(month_key).charge = struct();
                    all_events.(rackName).(year_key).(month_key).discharge = struct();
                end
                
                % 기존 충전/방전 이벤트 개수 확인
                charge_fields = fieldnames(all_events.(rackName).(year_key).(month_key).charge);
                discharge_fields = fieldnames(all_events.(rackName).(year_key).(month_key).discharge);
                
                % 충전 이벤트 개수 계산
                chg_count = 0;
                for f = 1:length(charge_fields)
                    if startsWith(charge_fields{f}, 'chg_evt')
                        chg_count = chg_count + 1;
                    end
                end
                
                % 방전 이벤트 개수 계산
                dchg_count = 0;
                for f = 1:length(discharge_fields)
                    if startsWith(discharge_fields{f}, 'dchg_evt')
                        dchg_count = dchg_count + 1;
                    end
                end
                
                % 충전 이벤트 저장 (charge 구조체 안에: chg_evt1, chg_evt2, ...)
                for k = 1:size(charge_events, 1)
                    s = charge_events(k, 1);
                    e = charge_events(k, 2);
                    evt = create_event_struct('charge', s, e, t, I_cell_raw, I_cell_smooth, Vcell_avg, V_cell_smooth, SOC_raw, P_cell_W, T, file_info);
                    if ~isempty(evt)
                        chg_count = chg_count + 1;
                        field_name = sprintf('chg_evt%d', chg_count);
                        all_events.(rackName).(year_key).(month_key).charge.(field_name) = evt;
                    end
                end
                
                % 방전 이벤트 저장 (discharge 구조체 안에: dchg_evt1, dchg_evt2, ...)
                for k = 1:size(discharge_events, 1)
                    s = discharge_events(k, 1);
                    e = discharge_events(k, 2);
                    evt = create_event_struct('discharge', s, e, t, I_cell_raw, I_cell_smooth, Vcell_avg, V_cell_smooth, SOC_raw, P_cell_W, T, file_info);
                    if ~isempty(evt)
                        dchg_count = dchg_count + 1;
                        field_name = sprintf('dchg_evt%d', dchg_count);
                        all_events.(rackName).(year_key).(month_key).discharge.(field_name) = evt;
                    end
                end
                
            end % file loop
        end % month loop
    end % year loop
end % dataType loop

%% 파일 저장 (EventDetection.m과 동일)
save(fullfile(saveDir, 'all_events_raw_cell_level.mat'), 'all_events', '-v7.3');

fprintf('\n=== Processing Complete ===\n');
fprintf('Raw cell-level events saved to: %s\n', fullfile(saveDir, 'all_events_raw_cell_level.mat'));

%% Helper Functions
% =========================================================================

function evt = create_event_struct(type, s, e, t_datetime, I_raw, I_smooth, V_raw, V_smooth, soc, P_cell_W, T, file_info)
    % EventDetection.m의 create_event_struct와 동일한 구조
    if (e-s) < 1
        evt = [];
        return;
    end
    
    delta_I_cell = I_smooth(s+1) - I_smooth(s);
    if abs(delta_I_cell) < 1e-3
        evt = [];
        return;
    end
    
    % 최대 전류 확인 (절댓값, smooth 데이터 사용)
    max_current = max(abs(I_smooth(s:e)));
    min_current_threshold = 64 * 0.1;  % 6.4A
    if max_current < min_current_threshold
        evt = [];
        return;
    end
    
    evt = struct();
    evt.type = type;
    evt.file_info = file_info;
    evt.time_seq_datetime = t_datetime(s:e);
    
    % Raw 데이터 저장
    evt.current_seq_cell_A_raw = I_raw(s:e);
    evt.voltage_seq_cell_V_raw = V_raw(s:e);
    
    % Smooth 데이터 저장
    evt.current_seq_cell_A = I_smooth(s:e);
    evt.voltage_seq_cell_V = V_smooth(s:e);
    
    evt.soc_seq_pct = soc(s:e);
    evt.power_seq_cell_W = P_cell_W(s:e);
    evt.temp_seq_C = T(s:e);
    
    % DCIR 계산은 smooth 데이터 사용
    evt.dcir_cell_mOhm = (V_smooth(s+1) - V_smooth(s)) / delta_I_cell * 1000;
end

function [charge_events, discharge_events] = detect_events_idle_to_active(...
    I_cell, idle_thr, min_duration_sec, t_datetime, trim_sec, max_I_std)
    % 이벤트 검출: idle -> active 전환 감지
    % - Simple_EventDetection 방식
    % - 이벤트 시작점: idle -> Active에서 idle 마지막 값
    % - 필터링: duration 300초 이상, 3초 ~ end-3초 구간에서 표준편차 기반 필터링
    
    charge_events = [];
    discharge_events = [];
    n = length(I_cell);
    
    if n < 10
        return;
    end
    
    % 1. Idle -> Active 전환 감지
    is_idle = abs(I_cell) < idle_thr;
    is_active = abs(I_cell) >= idle_thr;
    idle_to_active = find(is_idle(1:end-1) & is_active(2:end));
    
    if isempty(idle_to_active)
        return;
    end
    
    % 2. 각 전환점에서 이벤트 검출
    for k = 1:length(idle_to_active)
        % idle -> active 전환점: idle_to_active(k)는 idle 마지막 인덱스
        % active 시작점: idle_to_active(k) + 1
        idle_last_idx = idle_to_active(k);
        active_start_idx = idle_last_idx + 1;
        
        % 이벤트 타입 판별 (충전 또는 방전)
        event_type = sign(I_cell(active_start_idx));
        if event_type == 0
            continue;  % 전류가 0이면 스킵
        end
        
        % 3. 동일한 부호로 지속되는 동안 이벤트 끝 찾기
        active_end_idx = active_start_idx;
        while active_end_idx < n
            if (event_type > 0 && I_cell(active_end_idx + 1) >= idle_thr) || ...
               (event_type < 0 && I_cell(active_end_idx + 1) <= -idle_thr)
                active_end_idx = active_end_idx + 1;
            else
                break;
            end
        end
        
        % 이벤트 시작점: idle 마지막 값 (idle_last_idx)
        event_start = idle_last_idx;
        event_end = active_end_idx;
        
        % 4. 최소 지속 시간 확인
        if isdatetime(t_datetime) && length(t_datetime) >= event_end
            duration_sec = seconds(t_datetime(event_end) - t_datetime(event_start));
        else
            duration_sec = event_end - event_start;
        end
        
        if duration_sec < min_duration_sec
            continue;  % 최소 지속 시간 미만이면 스킵
        end
        
        % 5. 표준편차 기반 필터링: 3초 ~ end-3초 구간
        t_start_calc = t_datetime(event_start) + seconds(trim_sec);
        t_end_calc = t_datetime(event_end) - seconds(trim_sec);
        
        if t_end_calc <= t_start_calc
            continue;  % 계산 구간이 없으면 스킵
        end
        
        % 계산 구간 인덱스 찾기
        calc_idx = find(t_datetime >= t_start_calc & t_datetime <= t_end_calc);
        calc_idx = calc_idx(calc_idx >= event_start & calc_idx <= event_end);
        
        if length(calc_idx) < 5
            continue;  % 데이터 포인트가 너무 적으면 스킵
        end
        
        % 표준편차 계산
        I_std_val = std(I_cell(calc_idx));
        
        if I_std_val > max_I_std
            continue;  % 표준편차가 임계값 초과하면 스킵
        end
        
        % 6. 이벤트 저장
        if event_type > 0
            charge_events = [charge_events; event_start, event_end];
        else
            discharge_events = [discharge_events; event_start, event_end];
        end
    end
end
