%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 파일명: EventDetection.m (v16.2 - 전체 코드 최종)
% 기능:
% - 이벤트 시작점을 '상태 변화 직전의 마지막 휴지 시점'으로 정의
% - 셀 단위 변환 및 월별 계층 저장 기능
% - 수정된 그룹화 파라미터 적용
% MATLAB R2022a 이상 권장
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% 1. 기본 설정 (고정 및 사용자 지정)
% =========================================================================
% 스크립트 위치를 기준으로 상대 경로 설정
scriptDir = fileparts(mfilename('fullpath'));
% --- 고정 경로 ---
saveDir = fullfile(scriptDir, 'EventsResults');  % 현재 스크립트 위치의 EventsResults 폴더

% --- 사용자 설정 경로 및 파라미터 ---
% 프로젝트 루트로 이동 (FieldData\FieldData_Rpeak -> FieldData -> KEPCO_ESS_Local -> Rack_raw2mat)
projectRoot = fullfile(scriptDir, '..', '..');
dataDir = fullfile(projectRoot, 'Rack_raw2mat');
yearList = {'2021', '2022', '2023', '2024', '2025'};
rackNames_all = {'Rack01'};

% --- 배터리 구성 파라미터 ---
Ns = 14 * 17; Np = 2;
C_nom_cell = 64;

% --- 이벤트 검출 파라미터 ---
idle_thr_cell = 64*0.03;  % 휴지 상태 전류 임계값 (1.92A)
min_duration_sec = 30;  % 이벤트로 인정할 최소 지속 시간 (초)

% =========================================================================

%% 2. 스크립트 실행 부분
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

% 월별 계층적 구조로 결과 변수 초기화
all_events = struct();
for r_idx = 1:length(rackNames_all)
    rackName = rackNames_all{r_idx};
    all_events.(rackName) = struct();
    for y_idx = 1:length(yearList)
        year_key = ['Y' yearList{y_idx}];
        all_events.(rackName).(year_key) = struct();
    end
end

dataTypes = {'New', 'Old'};

for type_idx = 1:length(dataTypes)
    type = dataTypes{type_idx};
    typePath = fullfile(dataDir, type);
    if ~exist(typePath, 'dir'), continue; end
    
    fprintf('\n################### Processing Data Type: %s ###################\n', type);

    for year_idx = 1:length(yearList)
        year = yearList{year_idx};
        year_key = ['Y' year];
        yearPath = fullfile(typePath, year);
        if ~exist(yearPath, 'dir'), continue; end
        
        fprintf('\n================== Processing Year: %s ==================\n', year);

        monthDirs = dir(yearPath);
        monthDirs = monthDirs([monthDirs.isdir] & ~ismember({monthDirs.name}, {'.', '..'}));

        for m = 1:length(monthDirs)
            month = monthDirs(m).name;
            if isempty(regexp(month, '^\d{6}$', 'once')), continue; end
            
            month_key = ['M' month(end-1:end)];
            monthPath = fullfile(yearPath, month);
            matFiles = dir(fullfile(monthPath, '*.mat'));
            fprintf('\n--- Processing Month: %s (%d files) ---\n', month, length(matFiles));

            for f = 1:length(matFiles)
                matFilePath = fullfile(monthPath, matFiles(f).name);
                fprintf('Loading file: %s ... ', matFiles(f).name);
                data = load(matFilePath);
                
                file_info = struct('year', year, 'dataType', type, 'filename', matFiles(f).name);

                if strcmp(type, 'New')
                    all_events = process_rack_data(data, 'Rack01', type, file_info, year_key, month_key, all_events, idle_thr_cell, min_duration_sec);
                elseif strcmp(type, 'Old')
                    for r_idx = 1:length(rackNames_all)
                        rackName = rackNames_all{r_idx};
                        all_events = process_rack_data(data, rackName, type, file_info, year_key, month_key, all_events, idle_thr_cell, min_duration_sec);
                    end
                end

            end % file loop
        end % month loop
    end % year loop
end % dataType loop

%% 3. 파일 저장
save(fullfile(saveDir, 'all_events_raw_cell_level.mat'), 'all_events', '-v7.3');

fprintf('\nProcessing complete!\n');
fprintf('Raw cell-level events saved to: %s\n', fullfile(saveDir, 'all_events_raw_cell_level.mat'));


%% Helper Functions
% =========================================================================
function all_events_out = process_rack_data(data, rackName, type, file_info, year_key, month_key, all_events_in, idle_thr_cell, min_duration_sec)
    all_events_out = all_events_in;
    Ns = 14 * 17; Np = 2;
    
    if strcmp(type, 'New')
            rackData = data.Raw; 
            % New 타입: Date_Time이 시간만 있을 수 있으므로 파일명에서 날짜 추출
            filename = file_info.filename;
            fileDate = [];
            % 파일명에서 날짜 추출 (예: Raw_20230907.mat -> 20230907)
            % 파일명에서 8자리 숫자 찾기 (YYYYMMDD)
            dateMatch = regexp(filename, '\d{8}', 'once', 'match');
            if ~isempty(dateMatch)
                fileDate = datetime(dateMatch, 'InputFormat', 'yyyyMMdd');
            else
                fileDate = [];
            end
            
            % Date_Time 처리 - duration이면 파일명 날짜와 결합하여 datetime으로 변환
            if isduration(rackData.Date_Time)
                if ~isempty(fileDate)
                    % duration을 datetime으로 변환 (파일명 날짜 + duration)
                    t_datetime = fileDate + rackData.Date_Time;
                else
                    % 파일명 날짜가 없으면 에러 (duration만으로는 datetime 생성 불가)
                    error('Date_Time is duration but cannot extract date from filename: %s', filename);
                end
            elseif isdatetime(rackData.Date_Time)
                t_datetime = rackData.Date_Time;
                % 날짜가 없고 시간만 있는 경우 파일명 날짜와 결합
                if ~isempty(fileDate) && any(t_datetime.Year == 0)
                    t_datetime = fileDate + timeofday(t_datetime);
                end
            else
                % 기타 형식은 원본 그대로 사용 (이미 datetime일 수 있음)
                t_datetime = rackData.Date_Time;
            end
            
            % 최종 확인: t_datetime이 datetime인지 확인
            if ~isdatetime(t_datetime)
                error('t_datetime is not datetime after processing. Type: %s', class(t_datetime));
            end
            I_rack = rackData.DCCurrent; V_cell = rackData.CVavg; soc = rackData.SOC_BMS; P_rack_kW = rackData.DCPower; T = rackData.MTavg;
        else
            % Old 타입: Time을 datetime으로 변환
            rackData = data.Raw.(rackName);
            if isdatetime(rackData.Time)
                t_datetime = rackData.Time;
            elseif iscell(rackData.Time)
                % 셀 배열인 경우 첫 번째 요소 확인 후 변환
                if isdatetime(rackData.Time{1})
                    t_datetime = [rackData.Time{:}];
                elseif isnumeric(rackData.Time{1})
                    t_datetime = datetime(cell2mat(rackData.Time), 'ConvertFrom', 'datenum');
                else
                    t_datetime = datetime(rackData.Time);
                end
            elseif isnumeric(rackData.Time)
                % 숫자 형식인 경우 datenum으로 변환
                t_datetime = datetime(rackData.Time, 'ConvertFrom', 'datenum');
            else
                % 기타 형식은 datetime 생성자로 변환
                t_datetime = datetime(rackData.Time);
            end
            I_rack = rackData.DCCurrent_A; V_cell = rackData.AverageCV_V; soc = rackData.SOCPct; P_rack_kW = rackData.DCPower_kW; T = rackData.AverageMT_degC;
        end
        P_rack_W = P_rack_kW * 1000;
        vars = {I_rack, V_cell, soc, P_rack_W, T};
        min_len = min([length(t_datetime), cellfun(@length, vars)]);
        t_datetime = t_datetime(1:min_len); I_rack = I_rack(1:min_len); V_cell = V_cell(1:min_len); soc = soc(1:min_len); P_rack_W = P_rack_W(1:min_len); T = T(1:min_len);
        if min_len < min_duration_sec, return; end
        
        % 원본 데이터 준비 (raw)
        I_raw = I_rack / Np;
        V_raw = V_cell;
        P_cell_W = P_rack_W / (Ns * Np);
        
        % SG 필터 적용 (Savitzky-Golay): 노이즈 제거하되 엣지 유지
        % 윈도우 5 (앞뒤 2초), 차수 1: 저항값 왜곡 최소화하면서 노이즈 제거
        % 데이터가 1초 주기이므로 윈도우 5가 가장 이상적
        if length(I_raw) >= 5
            I_smooth = sgolayfilt(I_raw, 1, 5);  % 전류 필터링
            V_smooth = sgolayfilt(V_raw, 1, 5);  % 전압 필터링 (위상 맞춤)
        else
            I_smooth = I_raw;
            V_smooth = V_raw;
        end
        
        % 이벤트 검출: smooth 데이터 사용
        [chg_indices, dch_indices] = find_events_simplified(I_smooth, idle_thr_cell, min_duration_sec);
        
        % 이벤트 저장: raw와 smooth 둘 다 저장
        for k = 1:size(chg_indices, 1)
            s = chg_indices(k, 1); e = chg_indices(k, 2);
            evt = create_event_struct('charge', s, e, t_datetime, I_raw, I_smooth, V_raw, V_smooth, soc, P_cell_W, T, file_info);
            if ~isempty(evt), if ~isfield(all_events_out.(rackName).(year_key), month_key), all_events_out.(rackName).(year_key).(month_key) = {}; end; all_events_out.(rackName).(year_key).(month_key){end+1} = evt; end
        end
        for k = 1:size(dch_indices, 1)
            s = dch_indices(k, 1); e = dch_indices(k, 2);
            evt = create_event_struct('discharge', s, e, t_datetime, I_raw, I_smooth, V_raw, V_smooth, soc, P_cell_W, T, file_info);
            if ~isempty(evt), if ~isfield(all_events_out.(rackName).(year_key), month_key), all_events_out.(rackName).(year_key).(month_key) = {}; end; all_events_out.(rackName).(year_key).(month_key){end+1} = evt; end
        end
        fprintf('%s: Found %d Charge, %d Discharge events.\n', rackName, size(chg_indices, 1), size(dch_indices, 1));
end

function [chg_events, dch_events] = find_events_simplified(I_cell, idle_thr, min_duration)
    % 이벤트 검출: idle 상태의 마지막 시점부터 이벤트로 저장
    % idle_threshold 이상이면서 동일한 부호로 지속되는 경우만 검출
    chg_events = [];
    dch_events = [];
    n = length(I_cell);
    i = 1;
    
    while i <= n
        % 휴지 상태에서 벗어난 경우 (|I| >= idle_thr)
        if abs(I_cell(i)) >= idle_thr
            % idle 상태의 마지막 시점 찾기 (현재 지점 이전에서)
            event_start = i;  % 기본값: idle_thr 이상인 첫 지점
            % i 이전에서 idle 상태(|I| < idle_thr)의 마지막 지점 찾기
            for k = i-1:-1:1
                if abs(I_cell(k)) < idle_thr
                    event_start = k;  % idle 상태의 마지막 시점부터 저장
                    break;
                end
            end
            % 만약 처음부터 idle_thr 이상이면 event_start = 1
            if event_start > i
                event_start = 1;
            end
            
            event_type = sign(I_cell(i)); % 1: charge, -1: discharge
            
            % 동일한 부호로 지속되는 동안 이벤트 끝 찾기
            j = i + 1;
            while j <= n
                % 같은 방향(부호)이고 idle_threshold 이상이면 계속
                if (event_type > 0 && I_cell(j) >= idle_thr) || (event_type < 0 && I_cell(j) <= -idle_thr)
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


function evt = create_event_struct(type, s, e, t_datetime, I_raw, I_smooth, V_raw, V_smooth, soc, P_cell_W, T, file_info)
    if (e-s) < 1, evt = []; return; end
    delta_I_cell = I_smooth(s+1) - I_smooth(s);
    if abs(delta_I_cell) < 1e-3, evt = []; return; end
    
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

% AnalyzeResistanceFeatures