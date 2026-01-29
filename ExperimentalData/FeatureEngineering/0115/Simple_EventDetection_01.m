%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 01_Simple_EventDetection.m
% 목적: Raw 데이터에서 충/방전 이벤트 구간만 심플하게 잘라내어 저장
% 기준: 
%   1. Idle -> Active 전환
%   2. 지속시간 10초 이상
%   3. 전류 표준편차(Fluctuation) 허용 범위 이내
%
% [수정] 모든 채널에서 동일한 이벤트를 찾아서 동일한 번호 할당
%        기준 채널(ch9)의 이벤트를 기준으로 다른 채널과 매칭
%
% 출력:
%   - Lab_DC_Events_Raw_*cyc.mat (각 사이클별 Raw 이벤트 데이터)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 설정
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
outputDir = fullfile(pwd, 'Results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% 파라미터
Cnom = 64; 
current_threshold = Cnom * 0.01; % 0.64A
min_duration = 30;               % 최소 30초
max_I_std = 1.5;                 % 전류 표준편차 허용 범위 2A
reference_channel = 'ch9';       % 기준 채널 (이 채널의 이벤트를 기준으로 매칭)
time_tolerance = 5.0;            % 이벤트 매칭 시 시간 허용 오차 (초)

% 이벤트 타입 선택: 'Charge', 'Discharge', 'Both'
event_type_selection = 'Charge';  % 'Charge': 충전만, 'Discharge': 방전만, 'Both': 모두

% 이벤트 타입 선택 검증
if ~ismember(event_type_selection, {'Charge', 'Discharge', 'Both'})
    error('event_type_selection must be ''Charge'', ''Discharge'', or ''Both''');
end

fprintf('=== Simple Event Detection (Multi-Channel Synchronized) ===\n');
fprintf('Data directory: %s\n', dataDir);
fprintf('Output directory: %s\n', outputDir);
fprintf('Reference channel: %s\n', reference_channel);
fprintf('Event type selection: %s\n', event_type_selection);

%% 사이클 자동 감지 및 루프
matFiles = dir(fullfile(dataDir, 'parsedDriveCycle_*cyc_filtered.mat'));
if isempty(matFiles)
    error('No data files found in %s', dataDir);
end

fprintf('Found %d cycle files\n', length(matFiles));

% 디버깅용: 이벤트 카운트 저장 구조체
eventCounts = struct();
cycleNameMap = containers.Map();
cycleNameReverseMap = containers.Map();

% 필터링 통계 수집용 구조체
filterStatsCollection = struct();

for i = 1:length(matFiles)
    fileName = matFiles(i).name;
    token = regexp(fileName, 'parsedDriveCycle_(\d+cyc)_filtered', 'tokens');
    if isempty(token), continue; end
    cycleType = token{1}{1};
    
    fprintf('\n--- Processing %s ---\n', cycleType);
    load(fullfile(dataDir, fileName));
    
    % 데이터 변수명 찾기
    varName = sprintf('parsedDriveCycle_%s', cycleType);
    if ~exist(varName, 'var')
        fprintf('  WARNING: Variable %s not found. Skipping.\n', varName);
        continue;
    end
    eval(sprintf('data_var = %s;', varName));
    
    % 결과 저장용 구조체
    rawEvents = struct();
    
    % 채널 목록
    channels = fieldnames(data_var);
    fprintf('  Channels found: %s\n', strjoin(channels, ', '));
    
    % 기준 채널 찾기
    ref_channel_idx = [];
    for ch = 1:length(channels)
        match = regexp(channels{ch}, ['^' reference_channel], 'once');
        if ~isempty(match)
            ref_channel_idx = ch;
            break;
        end
    end
    
    if isempty(ref_channel_idx)
        fprintf('  WARNING: Reference channel %s not found. Skipping.\n', reference_channel);
        continue;
    end
    
    ref_channel_name = channels{ref_channel_idx};
    fprintf('  Reference channel: %s\n', ref_channel_name);
    
    % SOC -> Profile 루프로 변경
    % 먼저 SOC와 Profile 목록을 기준 채널에서 가져옴
    if ~isfield(data_var, ref_channel_name) || ~isstruct(data_var.(ref_channel_name))
        continue;
    end
    
    socs = fieldnames(data_var.(ref_channel_name));
    
    for s = 1:length(socs)
        socName = socs{s};
        if ~isfield(data_var.(ref_channel_name), socName) || ~isstruct(data_var.(ref_channel_name).(socName))
            continue;
        end
        
        profs = fieldnames(data_var.(ref_channel_name).(socName));
        
        for p = 1:length(profs)
            profName = profs{p};
            if ~isfield(data_var.(ref_channel_name).(socName), profName)
                continue;
            end
            
            fprintf('    Processing %s %s...\n', socName, profName);
            
            % 1단계: 기준 채널(ch9)에서 이벤트 검출
            ref_data = data_var.(ref_channel_name).(socName).(profName);
            if ~isfield(ref_data, 'I') || ~isfield(ref_data, 'V') || ~isfield(ref_data, 't')
                continue;
            end
            
            ref_I = ref_data.I;
            ref_V = ref_data.V;
            ref_t = ref_data.t;
            
            if isa(ref_t, 'duration')
                ref_t_seconds = seconds(ref_t);
            else
                ref_t_seconds = ref_t;
            end
            
            if length(ref_t_seconds) > 1
                ref_dt = median(diff(ref_t_seconds));
            else
                ref_dt = 1.0;
            end
            
            if ref_dt <= 0 || isnan(ref_dt)
                continue;
            end
            
            [ref_events, ref_filter_stats] = detect_events_in_channel(ref_I, ref_V, ref_t_seconds, current_threshold, ...
                                                   min_duration, max_I_std, ref_dt);
            
            % event_type_selection에 맞는 통계만 선택
            if strcmp(event_type_selection, 'Charge')
                ref_filter_stats.total_candidates = ref_filter_stats.total_candidates_charge;
                ref_filter_stats.filtered_by_duration = ref_filter_stats.filtered_by_duration_charge;
                ref_filter_stats.filtered_by_std = ref_filter_stats.filtered_by_std_charge;
                ref_filter_stats.filtered_by_direction = ref_filter_stats.filtered_by_direction_charge;
                ref_filter_stats.final_count = ref_filter_stats.final_count_charge;
            elseif strcmp(event_type_selection, 'Discharge')
                ref_filter_stats.total_candidates = ref_filter_stats.total_candidates_discharge;
                ref_filter_stats.filtered_by_duration = ref_filter_stats.filtered_by_duration_discharge;
                ref_filter_stats.filtered_by_std = ref_filter_stats.filtered_by_std_discharge;
                ref_filter_stats.filtered_by_direction = 0; % Discharge는 방향 체크 없음
                ref_filter_stats.final_count = ref_filter_stats.final_count_discharge;
            else % 'Both'
                ref_filter_stats.total_candidates = ref_filter_stats.total_candidates_charge + ref_filter_stats.total_candidates_discharge;
                ref_filter_stats.filtered_by_duration = ref_filter_stats.filtered_by_duration_charge + ref_filter_stats.filtered_by_duration_discharge;
                ref_filter_stats.filtered_by_std = ref_filter_stats.filtered_by_std_charge + ref_filter_stats.filtered_by_std_discharge;
                ref_filter_stats.filtered_by_direction = ref_filter_stats.filtered_by_direction_charge;
                ref_filter_stats.final_count = ref_filter_stats.final_count_charge + ref_filter_stats.final_count_discharge;
            end
            
            % 필터링 통계 저장 (기준 채널)
            statKey_raw = sprintf('%s_%s_%s_%s', cycleType, ref_channel_name, socName, profName);
            % 유효한 필드명으로 변환 (숫자로 시작하면 안 됨)
            statKey = regexprep(statKey_raw, '^(\d)', 'x$1'); % 숫자로 시작하면 앞에 'x' 추가
            statKey = regexprep(statKey, '[^a-zA-Z0-9_]', '_'); % 특수문자를 언더스코어로 변환
            filterStatsCollection.(statKey) = ref_filter_stats;
            filterStatsCollection.(statKey).channel = ref_channel_name;
            filterStatsCollection.(statKey).cycle = cycleType;
            filterStatsCollection.(statKey).soc = socName;
            filterStatsCollection.(statKey).profile = profName;
            
            if isempty(ref_events)
                fprintf('      No events found in reference channel. Skipping.\n');
                continue;
            end
            
            fprintf('      Reference channel: %d events detected\n', length(ref_events));
            
            % 기준 채널 이벤트 저장 (이벤트 번호는 순서대로 1, 2, 3, ...)
            for evt_num = 1:length(ref_events)
                evt = ref_events(evt_num);
                evtField = sprintf('event%d', evt_num);
                structName = sprintf('%s_%s', ref_channel_name, evt.type);
                
                if ~isfield(rawEvents, structName)
                    rawEvents.(structName) = struct();
                end
                if ~isfield(rawEvents.(structName), socName)
                    rawEvents.(structName).(socName) = struct();
                end
                if ~isfield(rawEvents.(structName).(socName), profName)
                    rawEvents.(structName).(socName).(profName) = struct();
                end
                
                rawEvents.(structName).(socName).(profName).(evtField).t = evt.t_seg;
                rawEvents.(structName).(socName).(profName).(evtField).V = evt.V_seg;
                rawEvents.(structName).(socName).(profName).(evtField).I = evt.I_seg;
                rawEvents.(structName).(socName).(profName).(evtField).duration = evt.duration;
            end
            
            % 디버깅용: 기준 채널 이벤트 카운트
            match_ref = regexp(ref_channel_name, '^(ch\d+)', 'tokens', 'once');
            if ~isempty(match_ref)
                chNumOnly_ref = match_ref{1};
            else
                chNumOnly_ref = ref_channel_name;
            end
            
            if ~isKey(cycleNameMap, cycleType)
                if ~isempty(regexp(cycleType, '^\d', 'once'))
                    validCycleName = ['cyc_' cycleType];
                else
                    validCycleName = cycleType;
                end
                cycleNameMap(cycleType) = validCycleName;
                cycleNameReverseMap(validCycleName) = cycleType;
            else
                validCycleName = cycleNameMap(cycleType);
            end
            
            % 선택된 타입에 맞는 이벤트만 카운트
            if strcmp(event_type_selection, 'Both')
                selected_count = length(ref_events);
            else
                selected_count = sum(strcmp({ref_events.type}, event_type_selection));
            end
            
            if selected_count > 0
                if ~isfield(eventCounts, profName)
                    eventCounts.(profName) = struct();
                end
                if ~isfield(eventCounts.(profName), chNumOnly_ref)
                    eventCounts.(profName).(chNumOnly_ref) = struct();
                end
                if ~isfield(eventCounts.(profName).(chNumOnly_ref), validCycleName)
                    eventCounts.(profName).(chNumOnly_ref).(validCycleName) = 0;
                end
                eventCounts.(profName).(chNumOnly_ref).(validCycleName) = ...
                    eventCounts.(profName).(chNumOnly_ref).(validCycleName) + selected_count;
            end
            
            % 2단계: 다른 모든 채널에서 이벤트 검출 및 매칭
            for ch = 1:length(channels)
                chName = channels{ch};
                
                % 기준 채널은 이미 처리했으므로 스킵
                if strcmp(chName, ref_channel_name)
                    continue;
                end
                
                if ~isfield(data_var, chName) || ~isstruct(data_var.(chName))
                    continue;
                end
                if ~isfield(data_var.(chName), socName) || ~isstruct(data_var.(chName).(socName))
                    continue;
                end
                if ~isfield(data_var.(chName).(socName), profName)
                    continue;
                end
                
                ch_data = data_var.(chName).(socName).(profName);
                if ~isfield(ch_data, 'I') || ~isfield(ch_data, 'V') || ~isfield(ch_data, 't')
                    continue;
                end
                
                ch_I = ch_data.I;
                ch_V = ch_data.V;
                ch_t = ch_data.t;
                
                if isa(ch_t, 'duration')
                    ch_t_seconds = seconds(ch_t);
                else
                    ch_t_seconds = ch_t;
                end
                
                if length(ch_t_seconds) > 1
                    ch_dt = median(diff(ch_t_seconds));
                else
                    ch_dt = 1.0;
                end
                
                if ch_dt <= 0 || isnan(ch_dt)
                    continue;
                end
                
                % 채널에서 이벤트 검출
                [ch_events, ch_filter_stats] = detect_events_in_channel(ch_I, ch_V, ch_t_seconds, current_threshold, ...
                                                      min_duration, max_I_std, ch_dt);
                
                % event_type_selection에 맞는 통계만 선택
                if strcmp(event_type_selection, 'Charge')
                    ch_filter_stats.total_candidates = ch_filter_stats.total_candidates_charge;
                    ch_filter_stats.filtered_by_duration = ch_filter_stats.filtered_by_duration_charge;
                    ch_filter_stats.filtered_by_std = ch_filter_stats.filtered_by_std_charge;
                    ch_filter_stats.filtered_by_direction = ch_filter_stats.filtered_by_direction_charge;
                    ch_filter_stats.final_count = ch_filter_stats.final_count_charge;
                elseif strcmp(event_type_selection, 'Discharge')
                    ch_filter_stats.total_candidates = ch_filter_stats.total_candidates_discharge;
                    ch_filter_stats.filtered_by_duration = ch_filter_stats.filtered_by_duration_discharge;
                    ch_filter_stats.filtered_by_std = ch_filter_stats.filtered_by_std_discharge;
                    ch_filter_stats.filtered_by_direction = 0; % Discharge는 방향 체크 없음
                    ch_filter_stats.final_count = ch_filter_stats.final_count_discharge;
                else % 'Both'
                    ch_filter_stats.total_candidates = ch_filter_stats.total_candidates_charge + ch_filter_stats.total_candidates_discharge;
                    ch_filter_stats.filtered_by_duration = ch_filter_stats.filtered_by_duration_charge + ch_filter_stats.filtered_by_duration_discharge;
                    ch_filter_stats.filtered_by_std = ch_filter_stats.filtered_by_std_charge + ch_filter_stats.filtered_by_std_discharge;
                    ch_filter_stats.filtered_by_direction = ch_filter_stats.filtered_by_direction_charge;
                    ch_filter_stats.final_count = ch_filter_stats.final_count_charge + ch_filter_stats.final_count_discharge;
                end
                
                % 필터링 통계 저장 (다른 채널)
                statKey_raw = sprintf('%s_%s_%s_%s', cycleType, chName, socName, profName);
                % 유효한 필드명으로 변환 (숫자로 시작하면 안 됨)
                statKey = regexprep(statKey_raw, '^(\d)', 'x$1'); % 숫자로 시작하면 앞에 'x' 추가
                statKey = regexprep(statKey, '[^a-zA-Z0-9_]', '_'); % 특수문자를 언더스코어로 변환
                filterStatsCollection.(statKey) = ch_filter_stats;
                filterStatsCollection.(statKey).channel = chName;
                filterStatsCollection.(statKey).cycle = cycleType;
                filterStatsCollection.(statKey).soc = socName;
                filterStatsCollection.(statKey).profile = profName;
                
                if isempty(ch_events)
                    continue;
                end
                
                % 기준 채널의 이벤트와 매칭
                % 각 채널 이벤트를 기준 이벤트와 매칭하고, 매칭된 이벤트 번호 사용
                matched_events = containers.Map('KeyType', 'double', 'ValueType', 'any');
                
                for ch_evt_idx = 1:length(ch_events)
                    ch_evt = ch_events(ch_evt_idx);
                    match_idx = match_event_to_reference(ch_evt.start_time, ch_evt.end_time, ...
                                                         ref_events, time_tolerance);
                    
                    if match_idx > 0
                        % 매칭 성공: 기준 이벤트 번호 사용
                        if ~isKey(matched_events, match_idx)
                            matched_events(match_idx) = ch_evt_idx;
                        else
                            % 이미 매칭된 경우, 겹침이 더 큰 것으로 선택
                            existing_ch_evt_idx = matched_events(match_idx);
                            existing_ch_evt = ch_events(existing_ch_evt_idx);
                            ref_evt = ref_events(match_idx);
                            
                            % 새로운 매칭의 겹침 계산
                            overlap_new = max(0, min(ch_evt.end_time, ref_evt.end_time) - ...
                                             max(ch_evt.start_time, ref_evt.start_time));
                            overlap_existing = max(0, min(existing_ch_evt.end_time, ref_evt.end_time) - ...
                                                   max(existing_ch_evt.start_time, ref_evt.start_time));
                            
                            if overlap_new > overlap_existing
                                matched_events(match_idx) = ch_evt_idx;
                            end
                        end
                    end
                end
                
                % 매칭된 이벤트 저장
                match_keys = matched_events.keys;
                for k = 1:length(match_keys)
                    ref_evt_num = match_keys{k};
                    ch_evt_idx = matched_events(ref_evt_num);
                    ch_evt = ch_events(ch_evt_idx);
                    
                    evtField = sprintf('event%d', ref_evt_num);
                    structName = sprintf('%s_%s', chName, ch_evt.type);
                    
                    if ~isfield(rawEvents, structName)
                        rawEvents.(structName) = struct();
                    end
                    if ~isfield(rawEvents.(structName), socName)
                        rawEvents.(structName).(socName) = struct();
                    end
                    if ~isfield(rawEvents.(structName).(socName), profName)
                        rawEvents.(structName).(socName).(profName) = struct();
                    end
                    
                    rawEvents.(structName).(socName).(profName).(evtField).t = ch_evt.t_seg;
                    rawEvents.(structName).(socName).(profName).(evtField).V = ch_evt.V_seg;
                    rawEvents.(structName).(socName).(profName).(evtField).I = ch_evt.I_seg;
                    rawEvents.(structName).(socName).(profName).(evtField).duration = ch_evt.duration;
                end
                
                % 디버깅용: 채널 이벤트 카운트 (매칭된 선택 타입 이벤트만)
                match_ch = regexp(chName, '^(ch\d+)', 'tokens', 'once');
                if ~isempty(match_ch)
                    chNumOnly = match_ch{1};
                else
                    chNumOnly = chName;
                end
                
                matched_selected_count = 0;
                for k = 1:length(match_keys)
                    ref_evt_num = match_keys{k};
                    ch_evt_idx = matched_events(ref_evt_num);
                    if strcmp(event_type_selection, 'Both') || ...
                       strcmp(ch_events(ch_evt_idx).type, event_type_selection)
                        matched_selected_count = matched_selected_count + 1;
                    end
                end
                
                if matched_selected_count > 0
                    if ~isfield(eventCounts, profName)
                        eventCounts.(profName) = struct();
                    end
                    if ~isfield(eventCounts.(profName), chNumOnly)
                        eventCounts.(profName).(chNumOnly) = struct();
                    end
                    if ~isfield(eventCounts.(profName).(chNumOnly), validCycleName)
                        eventCounts.(profName).(chNumOnly).(validCycleName) = 0;
                    end
                    eventCounts.(profName).(chNumOnly).(validCycleName) = ...
                        eventCounts.(profName).(chNumOnly).(validCycleName) + matched_selected_count;
                end
            end
        end
    end
    
    % 저장
    savePath = fullfile(outputDir, sprintf('Lab_DC_Events_Raw_%s.mat', cycleType));
    save(savePath, 'rawEvents');
    fprintf('  Saved events to: %s\n', savePath);
end

%% 디버깅: 필터링 통계 표시 (DC별)
fprintf('\n=== 디버깅: 필터링 통계 (DC별) ===\n');

statsFields = fieldnames(filterStatsCollection);
if ~isempty(statsFields)
    % DC별로 통계 집계 (구조체 사용)
    dcStatsStruct = struct();
    
    for s = 1:length(statsFields)
        statKey = statsFields{s};
        stats = filterStatsCollection.(statKey);
        
        if isfield(stats, 'profile') && isfield(stats, 'total_candidates')
            profName = stats.profile;
            
            % DC별로 그룹화
            if ~isfield(dcStatsStruct, profName)
                dcStatsStruct.(profName) = struct();
                dcStatsStruct.(profName).total_candidates = 0;
                dcStatsStruct.(profName).filtered_by_duration = 0;
                dcStatsStruct.(profName).filtered_by_std = 0;
                dcStatsStruct.(profName).filtered_by_direction = 0;
                dcStatsStruct.(profName).final_count = 0;
            end
            
            % 통계 합산
            dcStatsStruct.(profName).total_candidates = dcStatsStruct.(profName).total_candidates + stats.total_candidates;
            dcStatsStruct.(profName).filtered_by_duration = dcStatsStruct.(profName).filtered_by_duration + stats.filtered_by_duration;
            dcStatsStruct.(profName).filtered_by_std = dcStatsStruct.(profName).filtered_by_std + stats.filtered_by_std;
            dcStatsStruct.(profName).filtered_by_direction = dcStatsStruct.(profName).filtered_by_direction + stats.filtered_by_direction;
            dcStatsStruct.(profName).final_count = dcStatsStruct.(profName).final_count + stats.final_count;
        end
    end
    
    dcNames = fieldnames(dcStatsStruct);
    if ~isempty(dcNames)
        % DC 목록 정렬
        dcNamesSorted = sort(dcNames);
        
        % 표 헤더 (엑셀 복사용 쉼표 구분)
        fprintf('\n필터링 통계 (DC별 집계) - 엑셀 복사용 (쉼표 구분):\n');
        fprintf('DC,총후보,제거:시간,제거:표준편차,제거:방향,최종\n');
        
        % 각 DC별 통계 출력 (쉼표 구분)
        for d = 1:length(dcNamesSorted)
            dcName = dcNamesSorted{d};
            dcStats = dcStatsStruct.(dcName);
            
            fprintf('%s,%d,%d,%d,%d,%d\n', ...
                dcName, ...
                dcStats.total_candidates, ...
                dcStats.filtered_by_duration, ...
                dcStats.filtered_by_std, ...
                dcStats.filtered_by_direction, ...
                dcStats.final_count);
        end
        
        % 전체 합계
        totalAll = struct();
        totalAll.total_candidates = 0;
        totalAll.filtered_by_duration = 0;
        totalAll.filtered_by_std = 0;
        totalAll.filtered_by_direction = 0;
        totalAll.final_count = 0;
        
        for d = 1:length(dcNamesSorted)
            dcName = dcNamesSorted{d};
            dcStats = dcStatsStruct.(dcName);
            totalAll.total_candidates = totalAll.total_candidates + dcStats.total_candidates;
            totalAll.filtered_by_duration = totalAll.filtered_by_duration + dcStats.filtered_by_duration;
            totalAll.filtered_by_std = totalAll.filtered_by_std + dcStats.filtered_by_std;
            totalAll.filtered_by_direction = totalAll.filtered_by_direction + dcStats.filtered_by_direction;
            totalAll.final_count = totalAll.final_count + dcStats.final_count;
        end
        
        fprintf('전체,%d,%d,%d,%d,%d\n', ...
            totalAll.total_candidates, ...
            totalAll.filtered_by_duration, ...
            totalAll.filtered_by_std, ...
            totalAll.filtered_by_direction, ...
            totalAll.final_count);
        
        % 비율 표시
        if totalAll.total_candidates > 0
            fprintf('\n비율 (전체):\n');
            fprintf('%-40s: %8d (%5.1f%%)\n', '총 후보 이벤트 수', totalAll.total_candidates, 100.0);
            fprintf('%-40s: %8d (%5.1f%%)\n', '제거: 지속시간 부족', totalAll.filtered_by_duration, totalAll.filtered_by_duration/totalAll.total_candidates*100);
            fprintf('%-40s: %8d (%5.1f%%)\n', '제거: 전류 표준편차 초과', totalAll.filtered_by_std, totalAll.filtered_by_std/totalAll.total_candidates*100);
            fprintf('%-40s: %8d (%5.1f%%)\n', '제거: 방향 불일치', totalAll.filtered_by_direction, totalAll.filtered_by_direction/totalAll.total_candidates*100);
            fprintf('%-40s: %8d (%5.1f%%)\n', '최종 채택 이벤트 수', totalAll.final_count, totalAll.final_count/totalAll.total_candidates*100);
        end
    end
end

%% 디버깅: 각 DC별 이벤트 검출 현황 표시
fprintf('\n=== 디버깅: DC별 이벤트 검출 현황 ===\n');
fprintf('(저장된 파일에서 실제 이벤트 개수 계산)\n');

% 저장된 파일들에서 실제 이벤트 개수 계산
filteredEventCounts = struct();
for i = 1:length(matFiles)
    fileName = matFiles(i).name;
    token = regexp(fileName, 'parsedDriveCycle_(\d+cyc)_filtered', 'tokens');
    if isempty(token), continue; end
    cycleType = token{1}{1};
    
    % 사이클명을 유효한 필드명으로 변환
    if ~isKey(cycleNameMap, cycleType)
        if ~isempty(regexp(cycleType, '^\d', 'once'))
            validCycleName = ['cyc_' cycleType];
        else
            validCycleName = cycleType;
        end
        cycleNameMap(cycleType) = validCycleName;
        cycleNameReverseMap(validCycleName) = cycleType;
    else
        validCycleName = cycleNameMap(cycleType);
    end
    
    % 저장된 파일 로드
    savePath = fullfile(outputDir, sprintf('Lab_DC_Events_Raw_%s.mat', cycleType));
    if ~exist(savePath, 'file')
        continue;
    end
    
    tempData = load(savePath, 'rawEvents');
    if ~isfield(tempData, 'rawEvents')
        continue;
    end
    rawEvents_loaded = tempData.rawEvents;
    
    % 선택된 타입의 이벤트만 처리
    allStructNames = fieldnames(rawEvents_loaded);
    if strcmp(event_type_selection, 'Charge')
        selectedStructNames = allStructNames(contains(allStructNames, '_Charge'));
    elseif strcmp(event_type_selection, 'Discharge')
        selectedStructNames = allStructNames(contains(allStructNames, '_Discharge'));
    else % 'Both'
        selectedStructNames = allStructNames;
    end
    
    for cs = 1:length(selectedStructNames)
        structName = selectedStructNames{cs};
        % 채널명 추출
        match = regexp(structName, '^(ch\d+)_', 'tokens', 'once');
        if isempty(match), continue; end
        chNumOnly = match{1};
        
        socs = fieldnames(rawEvents_loaded.(structName));
        for s = 1:length(socs)
            socName = socs{s};
            profs = fieldnames(rawEvents_loaded.(structName).(socName));
            for p = 1:length(profs)
                profName = profs{p};
                
                % 이벤트 개수 계산 (각 SOC별로 개별 계산)
                events = fieldnames(rawEvents_loaded.(structName).(socName).(profName));
                eventCount = length(events);
                
                % 구조체에 저장 (모든 SOC의 이벤트를 누적)
                if ~isfield(filteredEventCounts, profName)
                    filteredEventCounts.(profName) = struct();
                end
                if ~isfield(filteredEventCounts.(profName), chNumOnly)
                    filteredEventCounts.(profName).(chNumOnly) = struct();
                end
                if ~isfield(filteredEventCounts.(profName).(chNumOnly), validCycleName)
                    filteredEventCounts.(profName).(chNumOnly).(validCycleName) = 0;
                end
                filteredEventCounts.(profName).(chNumOnly).(validCycleName) = ...
                    filteredEventCounts.(profName).(chNumOnly).(validCycleName) + eventCount;
            end
        end
    end
end

% 테이블 생성 및 표시
dcNames = fieldnames(filteredEventCounts);
for d = 1:length(dcNames)
    dcName = dcNames{d};
    fprintf('\n--- DC: %s ---\n', dcName);
    
    % 채널 목록 수집
    chNames = fieldnames(filteredEventCounts.(dcName));
    if isempty(chNames)
        fprintf('  이벤트 없음\n');
        continue;
    end
    
    % 사이클 목록 수집 (모든 채널에서 사용된 사이클)
    allCyclesValid = {};
    allCyclesDisplay = {};
    for c = 1:length(chNames)
        chName = chNames{c};
        if isfield(filteredEventCounts.(dcName), chName) && isstruct(filteredEventCounts.(dcName).(chName))
            cyclesValid = fieldnames(filteredEventCounts.(dcName).(chName));
            allCyclesValid = unique([allCyclesValid; cyclesValid]);
        end
    end
    
    if isempty(allCyclesValid)
        fprintf('  사이클 데이터 없음\n');
        continue;
    end
    
    % 유효 필드명을 원본 사이클명으로 역변환
    for cy = 1:length(allCyclesValid)
        validName = allCyclesValid{cy};
        if isKey(cycleNameReverseMap, validName)
            allCyclesDisplay{cy} = cycleNameReverseMap(validName);
        else
            allCyclesDisplay{cy} = validName;
        end
    end
    
    % 테이블 데이터 생성
    tableData = zeros(length(chNames), length(allCyclesValid));
    for c = 1:length(chNames)
        chName = chNames{c};
        for cy = 1:length(allCyclesValid)
            cycleNameValid = allCyclesValid{cy};
            if isfield(filteredEventCounts.(dcName), chName) && isstruct(filteredEventCounts.(dcName).(chName))
                if isfield(filteredEventCounts.(dcName).(chName), cycleNameValid)
                    tableData(c, cy) = filteredEventCounts.(dcName).(chName).(cycleNameValid);
                end
            end
        end
    end
    
    % 테이블 생성 및 표시
    T = array2table(tableData, 'RowNames', chNames, 'VariableNames', allCyclesDisplay);
    fprintf('\n');
    disp(T);
end

fprintf('\n=== Event Detection Complete ===\n');
fprintf('All results saved to: %s\n', outputDir);

%% DC별 전류 프로파일 시각화 (SOC70만)
fprintf('\n=== Creating Current Profile Visualization (SOC70) ===\n');

figuresDir = fullfile(outputDir, 'figures', 'CurrentProfiles');
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

targetSOC = {'SOC90','SOC70','SOC50'};  % SOC70만 시각화

% 기준 채널에서 데이터 추출
ref_channel_match = regexp(reference_channel, '^(ch\d+)', 'tokens', 'once');
if isempty(ref_channel_match)
    fprintf('  WARNING: Cannot extract reference channel number. Skipping visualization.\n');
else
    ref_ch_num = ref_channel_match{1};
    
    % 선택된 이벤트 타입에 맞는 구조체 이름 생성
    if strcmp(event_type_selection, 'Charge')
        eventTypes = {'Charge'};
    elseif strcmp(event_type_selection, 'Discharge')
        eventTypes = {'Discharge'};
    else % 'Both'
        eventTypes = {'Charge', 'Discharge'};
    end
    
    for etIdx = 1:length(eventTypes)
        eventType = eventTypes{etIdx};
        
        % 모든 사이클 파일에서 데이터 수집
        allCycleData = struct();
        
        for i = 1:length(matFiles)
            fileName = matFiles(i).name;
            token = regexp(fileName, 'parsedDriveCycle_(\d+cyc)_filtered', 'tokens');
            if isempty(token), continue; end
            cycleType = token{1}{1};
            
            savePath = fullfile(outputDir, sprintf('Lab_DC_Events_Raw_%s.mat', cycleType));
            if ~exist(savePath, 'file')
                continue;
            end
            
            tempData = load(savePath, 'rawEvents');
            if ~isfield(tempData, 'rawEvents')
                continue;
            end
            rawEvents_loaded = tempData.rawEvents;
            
            % 저장된 구조체 이름 찾기 (reference_channel로 시작하고 eventType으로 끝나는 것)
            allStructNames = fieldnames(rawEvents_loaded);
            structName = '';
            for s = 1:length(allStructNames)
                structName_candidate = allStructNames{s};
                % reference_channel로 시작하고 eventType으로 끝나는지 확인
                if startsWith(structName_candidate, reference_channel) && ...
                   endsWith(structName_candidate, ['_' eventType])
                    structName = structName_candidate;
                    break;
                end
            end
            
            if isempty(structName)
                continue;
            end
            
            % 실제 SOC 필드명 확인 (디버깅)
            availableSOCs = fieldnames(rawEvents_loaded.(structName));
            if ~ismember(targetSOC, availableSOCs)
                fprintf('    WARNING: %s not found in %s. Available SOCs: %s\n', ...
                        targetSOC, structName, strjoin(availableSOCs, ', '));
                continue;
            end
            
            profs = fieldnames(rawEvents_loaded.(structName).(targetSOC));
            
            for p = 1:length(profs)
                profName = profs{p};
                
                if ~isfield(rawEvents_loaded.(structName).(targetSOC), profName)
                    continue;
                end
                
                events = fieldnames(rawEvents_loaded.(structName).(targetSOC).(profName));
                
                if isempty(events)
                    continue;
                end
                
                % 첫 번째 이벤트의 데이터 사용 (대표 이벤트)
                if ~isfield(allCycleData, profName)
                    allCycleData.(profName) = struct();
                end
                
                evtData = rawEvents_loaded.(structName).(targetSOC).(profName).(events{1});
                
                if isfield(evtData, 'I') && isfield(evtData, 't')
                    % cycleType을 유효한 필드명으로 변환 (숫자로 시작하면 cyc_ 접두사 추가)
                    if ~isempty(regexp(cycleType, '^\d', 'once'))
                        validCycleFieldName = ['cyc_' cycleType];
                    else
                        validCycleFieldName = cycleType;
                    end
                    allCycleData.(profName).(validCycleFieldName) = evtData;
                end
            end
        end
        
        % DC1-8 고정 서브플롯 생성
        allDCProfiles = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};
        dcProfiles = fieldnames(allCycleData);
        
        fprintf('  Creating visualization for %s events\n', eventType);
        
        % 서브플롯으로 DC1-8 모두 표시 (2행 4열)
        nCols = 4;
        nRows = 2;
        
        fig = figure('Name', sprintf('Current Profiles - %s Events (%s)', eventType, targetSOC), ...
                     'Position', [100, 100, 1600, 900], 'Visible', 'on');
        
        colors = lines(length(matFiles));  % 사이클별 색상
        
        for dcIdx = 1:length(allDCProfiles)
            profName = allDCProfiles{dcIdx};
            
            subplot(nRows, nCols, dcIdx);
            hold on;
            
            legendEntries = {};
            legendHandles = [];
            
            % 데이터가 있는지 확인
            if ismember(profName, dcProfiles)
                % 각 사이클의 데이터 플롯
                cycles = fieldnames(allCycleData.(profName));
                
            if ~isempty(cycles)
                % 사이클 정렬 (cyc_ 접두사 처리)
                cycle_nums = zeros(length(cycles), 1);
                for c = 1:length(cycles)
                    cyc_str = cycles{c};
                    % cyc_ 접두사 제거하여 숫자 추출
                    if startsWith(cyc_str, 'cyc_')
                        cyc_str_clean = cyc_str(5:end);
                    else
                        cyc_str_clean = cyc_str;
                    end
                    num_match = regexp(cyc_str_clean, '(\d+)cyc', 'tokens');
                    if ~isempty(num_match)
                        cycle_nums(c) = str2double(num_match{1}{1});
                    else
                        cycle_nums(c) = 999;
                    end
                end
                [~, cyc_sort_idx] = sort(cycle_nums);
                sortedCycles = cycles(cyc_sort_idx);
                
                for cycIdx = 1:length(sortedCycles)
                    cycleFieldName = sortedCycles{cycIdx};
                    
                    if ~isfield(allCycleData.(profName), cycleFieldName)
                        continue;
                    end
                    
                    evtData = allCycleData.(profName).(cycleFieldName);
                    
                    if isfield(evtData, 'I') && isfield(evtData, 't') && length(evtData.t) > 0
                        t_rel = evtData.t - evtData.t(1);  % 상대 시간
                        
                        % 표시용 사이클 이름 (cyc_ 접두사 제거)
                        if startsWith(cycleFieldName, 'cyc_')
                            cycleDisplayName = cycleFieldName(5:end);
                        else
                            cycleDisplayName = cycleFieldName;
                        end
                        
                        h = plot(t_rel, evtData.I, '-', 'LineWidth', 1.5, ...
                                 'Color', colors(min(cycIdx, size(colors, 1)), :), ...
                                 'DisplayName', cycleDisplayName);
                        legendHandles(end+1) = h;
                        legendEntries{end+1} = cycleDisplayName;
                    end
                end
                end
            end
            
            xlabel('Time (s)', 'FontSize', 10, 'FontWeight', 'bold');
            ylabel('Current (A)', 'FontSize', 10, 'FontWeight', 'bold');
            title(profName, 'FontSize', 11, 'FontWeight', 'bold');
            grid on;
            
            if ~isempty(legendHandles)
                legend(legendHandles, legendEntries, 'Location', 'best', 'FontSize', 8);
            else
                % 데이터가 없으면 메시지 표시
                text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'Units', 'normalized', ...
                     'FontSize', 12, 'Color', [0.5 0.5 0.5]);
            end
            
            hold off;
        end
        
        % 전체 제목
        sgtitle(sprintf('Current Profiles by DC Profile - %s Events (%s, Channel %s)', ...
                        eventType, targetSOC, ref_ch_num), ...
                'FontSize', 14, 'FontWeight', 'bold');
        
        % 저장 (fig만)
        saveFileName = sprintf('CurrentProfiles_%s_%s_%s', eventType, targetSOC, ref_ch_num);
        saveFileName = regexprep(saveFileName, '[^a-zA-Z0-9_]', '_');
        savePath_fig = fullfile(figuresDir, [saveFileName, '.fig']);
        
        saveas(fig, savePath_fig);
        fprintf('  Saved: %s\n', savePath_fig);
        
        close(fig);
    end
end

fprintf('\n=== Current Profile Visualization Complete ===\n');

%% 로컬 함수들

function [events, filter_stats] = detect_events_in_channel(I, V, t_seconds, current_threshold, min_duration, max_I_std, dt)
    % 단일 채널에서 이벤트를 검출하고 반환
    % 반환값: 
    %   events: events 구조체 배열 (각 요소는 start_time, end_time, start_idx, end_idx, type, I_seg, V_seg, t_seg, duration 포함)
    %   filter_stats: 필터링 통계 구조체
    
    events = struct('start_time', {}, 'end_time', {}, 'start_idx', {}, 'end_idx', {}, ...
                    'type', {}, 'I_seg', {}, 'V_seg', {}, 't_seg', {}, 'duration', {});
    
    % 필터링 통계 초기화 (타입별로 분리)
    filter_stats = struct();
    filter_stats.total_candidates_charge = 0;
    filter_stats.total_candidates_discharge = 0;
    filter_stats.filtered_by_duration_charge = 0;
    filter_stats.filtered_by_duration_discharge = 0;
    filter_stats.filtered_by_std_charge = 0;
    filter_stats.filtered_by_std_discharge = 0;
    filter_stats.filtered_by_direction_charge = 0;
    filter_stats.final_count_charge = 0;
    filter_stats.final_count_discharge = 0;
    
    if length(t_seconds) < 10 || length(I) < 10 || length(V) < 10
        return;
    end
    
    % 1. 노이즈 제거 (3초 이동평균)
    target_smoothing_time = 3;
    moving_avg_window = max(1, round(target_smoothing_time / dt));
    if length(I) >= moving_avg_window
        I_filt = movmean(I, moving_avg_window);
        V_filt = movmean(V, moving_avg_window);
    else
        I_filt = movmean(I, length(I));
        V_filt = movmean(V, length(V));
    end
    
    % 2. 이벤트 감지 (Idle -> Active) - 필터링된 데이터 I_filt 사용
    is_idle = abs(I_filt) < current_threshold;
    is_driving = abs(I_filt) >= current_threshold;
    idle_to_driving = find(is_idle(1:end-1) & is_driving(2:end));
    
    if isempty(idle_to_driving)
        filter_stats.total_candidates_charge = 0;
        filter_stats.total_candidates_discharge = 0;
        return;
    end
    
    evt_count = 0;
    evt_count_charge = 0;
    evt_count_discharge = 0;
    for k = 1:length(idle_to_driving)
        idx1 = idle_to_driving(k);
        start_driving_idx = idx1 + 1;
        
        % 이벤트 타입 판별 - 필터링된 데이터 I_filt 사용
        event_type = sign(I_filt(start_driving_idx));
        % 주석: event_type == 0 체크는 불필요 (abs(I_filt) >= threshold이므로 0일 수 없음)
        
        % 타입별 후보 카운트
        if event_type > 0
            filter_stats.total_candidates_charge = filter_stats.total_candidates_charge + 1;
        else
            filter_stats.total_candidates_discharge = filter_stats.total_candidates_discharge + 1;
        end
        
        % driving 구간의 끝 찾기
        driving_end_idx = start_driving_idx;
        while driving_end_idx <= length(I_filt)
            if abs(I_filt(driving_end_idx)) >= current_threshold && sign(I_filt(driving_end_idx)) == event_type
                driving_end_idx = driving_end_idx + 1;
            else
                break;
            end
        end
        driving_end_idx = driving_end_idx - 1;
        
        start_idx = idx1;
        end_idx = driving_end_idx;
        
        % 지속시간 체크
        driving_time = t_seconds(driving_end_idx) - t_seconds(start_driving_idx);
        if driving_time < min_duration
            if event_type > 0
                filter_stats.filtered_by_duration_charge = filter_stats.filtered_by_duration_charge + 1;
            else
                filter_stats.filtered_by_duration_discharge = filter_stats.filtered_by_duration_discharge + 1;
            end
            continue;
        end
        
        % 세그먼트 추출 (필터링된 데이터 사용)
        I_seg = I_filt(start_idx:end_idx);
        V_seg = V_filt(start_idx:end_idx);
        t_seg = t_seconds(start_idx:end_idx);
        
        % 전류 안정성 체크
        I_std_val = std(I_seg);
        if I_std_val > max_I_std
            if event_type > 0
                filter_stats.filtered_by_std_charge = filter_stats.filtered_by_std_charge + 1;
            else
                filter_stats.filtered_by_std_discharge = filter_stats.filtered_by_std_discharge + 1;
            end
            continue;
        end
        
        % 방향 체크 (충전 이벤트는 대부분 양수여야 함)
        if event_type > 0
            positive_ratio = sum(I_seg > 0) / length(I_seg);
            if positive_ratio < 0.5
                filter_stats.filtered_by_direction_charge = filter_stats.filtered_by_direction_charge + 1;
                continue;
            end
        end
        
        % 이벤트 타입 문자열
        if event_type > 0
            type_str = 'Charge';
        else
            type_str = 'Discharge';
        end
        
        % 이벤트 저장
        evt_count = evt_count + 1;
        if event_type > 0
            evt_count_charge = evt_count_charge + 1;
        else
            evt_count_discharge = evt_count_discharge + 1;
        end
        events(evt_count).start_time = t_seconds(start_driving_idx);
        events(evt_count).end_time = t_seconds(driving_end_idx);
        events(evt_count).start_idx = start_idx;
        events(evt_count).end_idx = end_idx;
        events(evt_count).type = type_str;
        events(evt_count).I_seg = I_seg;
        events(evt_count).V_seg = V_seg;
        events(evt_count).t_seg = t_seg;
        events(evt_count).duration = driving_time;
    end
    
    filter_stats.final_count_charge = evt_count_charge;
    filter_stats.final_count_discharge = evt_count_discharge;
end

function match_idx = match_event_to_reference(event_start, event_end, ref_events, time_tolerance)
    % 단일 이벤트를 기준 이벤트들과 매칭
    % 반환값: 매칭된 기준 이벤트의 인덱스 (매칭 실패 시 0)
    
    match_idx = 0;
    min_overlap = 0;
    
    for i = 1:length(ref_events)
        ref_start = ref_events(i).start_time;
        ref_end = ref_events(i).end_time;
        
        % 시간 겹침 계산 (overlap ratio)
        overlap_start = max(event_start, ref_start);
        overlap_end = min(event_end, ref_end);
        overlap_time = max(0, overlap_end - overlap_start);
        event_duration = event_end - event_start;
        ref_duration = ref_end - ref_start;
        avg_duration = (event_duration + ref_duration) / 2;
        
        if overlap_time > 0 && avg_duration > 0
            overlap_ratio = overlap_time / avg_duration;
        else
            % 겹침이 없으면 시작 시간 차이로 판단
            % time_diff = min(abs(event_start - ref_start), abs(event_end - ref_end));
            % if time_diff <= time_tolerance
            %     overlap_ratio = 0.5; % 겹침은 없지만 시간적으로 가까움
            % else
                overlap_ratio = 0;
            % end
        end
        
        % 최대 겹침을 가진 이벤트 선택 (최소 30% 겹침 필요)
        if overlap_ratio > min_overlap && overlap_ratio >= 0.3
            min_overlap = overlap_ratio;
            match_idx = i;
        end
    end
    
    % 겹침이 부족한 경우 시작 시간 차이만으로 매칭 시도
    % if match_idx == 0
    %     min_time_diff = inf;
    %     for i = 1:length(ref_events)
    %         ref_start = ref_events(i).start_time;
    %         time_diff_start = abs(event_start - ref_start);
    %         if time_diff_start < min_time_diff && time_diff_start <= time_tolerance
    %             min_time_diff = time_diff_start;
    %             match_idx = i;
    %         end
    %     end
    % end
end
