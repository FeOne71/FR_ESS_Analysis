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

% 파라미터
Cnom = 64; 
current_threshold = Cnom * 0.01; % 0.64A
% min_duration = 30;               % 기본값 (paramSets에서 덮어씀)
% max_I_std = 2.0;                 % 기본값 (paramSets에서 덮어씀)ㅅ
% time_tolerance = 5.0;            % 이벤트 매칭 시 시간 허용 오차 (초)

% 이벤트 타입 선택: 'Charge', 'Discharge', 'Both'
reference_channel = 'ch11';       % 기준 채널 (이 채널의 이벤트를 기준으로 매칭)
reference_channel = regexprep(reference_channel, '^ch0+(\d+)$', 'ch$1');
event_type_selection = 'Both';  % 'Charge': 충전만, 'Discharge': 방전만, 'Both': 모두

% 파라미터 조합 (Systematic Combination)
min_duration_list = [10, 30, 60];
max_I_std_list = [2.0, 1.0, 0.5];
% max_I_range_list = [3.0, 1.5, 0.5];

paramSets = struct('min_duration', {}, 'max_I_std', {}, 'max_I_range', {});
setIdx = 0;
for md = min_duration_list
    for ms = max_I_std_list
        for mr = 0
            setIdx = setIdx + 1;
            paramSets(setIdx).min_duration = md;
            paramSets(setIdx).max_I_std = ms;
            paramSets(setIdx).max_I_range = mr;
        end
    end
end

% 이벤트 타입 선택 검증
if ~ismember(event_type_selection, {'Charge', 'Discharge', 'Both'})
    error('event_type_selection must be ''Charge'', ''Discharge'', or ''Both''');
end

fprintf('=== Simple Event Detection (Multi-Channel Synchronized) ===\n');
fprintf('Data directory: %s\n', dataDir);
fprintf('Reference channel: %s\n', reference_channel);
fprintf('Event type selection: %s\n', event_type_selection);

%% 사이클 자동 감지 및 루프
matFiles = dir(fullfile(dataDir, 'parsedDriveCycle_*cyc_filtered.mat'));
if isempty(matFiles)
    error('No data files found in %s', dataDir);
end

fprintf('Found %d cycle files\n', length(matFiles));

% 병렬 실행 (Parallel Computing Toolbox 필요)
for p = 1:length(paramSets)
    min_duration = paramSets(p).min_duration;
    max_I_std = paramSets(p).max_I_std;
    max_I_range = paramSets(p).max_I_range;
    paramLabel = sprintf('min%d_std%s_rng%s_%s', ...
        min_duration, num2str(max_I_std), num2str(max_I_range), event_type_selection);
    paramLabel = regexprep(paramLabel, '\.', 'p');
    paramLabel = regexprep(paramLabel, '[^a-zA-Z0-9_]', '_');
    outputDir = fullfile(pwd, 'Results', paramLabel);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    fprintf('\n=== ParamSet: %s (min_duration=%d, max_I_std=%.3f, max_I_range=%.3f) ===\n', ...
            paramLabel, min_duration, max_I_std, max_I_range);
    fprintf('Output directory: %s\n', outputDir);
    
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
    dataStruct = load(fullfile(dataDir, fileName));
    
    % 데이터 변수명 찾기
    varName = sprintf('parsedDriveCycle_%s', cycleType);
    if ~isfield(dataStruct, varName)
        fprintf('  WARNING: Variable %s not found. Skipping.\n', varName);
        continue;
    end
    data_var = dataStruct.(varName);
    
    % 결과 저장용 구조체
    rawEvents = struct();
    rejectedEvents = struct();
    
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
        
        for pIdx = 1:length(profs)
            profName = profs{pIdx};
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
            
            debug_tag = sprintf('%s %s %s %s', cycleType, socName, profName, ref_channel_name);
            [ref_events, ref_filter_stats, ref_rejected] = detect_events_in_channel(ref_I, ref_V, ref_t_seconds, current_threshold, ...
                                                   min_duration, max_I_std, max_I_range, ref_dt, debug_tag, event_type_selection);
            
            % event_type_selection에 맞는 통계만 선택
            if strcmp(event_type_selection, 'Charge')
                ref_filter_stats.total_candidates = ref_filter_stats.total_candidates_charge;
                ref_filter_stats.filtered_by_duration = ref_filter_stats.filtered_by_duration_charge;
                ref_filter_stats.filtered_by_std = ref_filter_stats.filtered_by_std_charge;
                ref_filter_stats.filtered_by_direction = ref_filter_stats.filtered_by_direction_charge;
                ref_filter_stats.filtered_by_range = ref_filter_stats.filtered_by_range_charge;
                ref_filter_stats.final_count = ref_filter_stats.final_count_charge;
                ref_filter_stats.fail_duration = ref_filter_stats.filtered_by_duration_charge;
                ref_filter_stats.fail_std = ref_filter_stats.filtered_by_std_charge;
                ref_filter_stats.fail_direction = ref_filter_stats.filtered_by_direction_charge;
                ref_filter_stats.fail_range = ref_filter_stats.filtered_by_range_charge;
            elseif strcmp(event_type_selection, 'Discharge')
                ref_filter_stats.total_candidates = ref_filter_stats.total_candidates_discharge;
                ref_filter_stats.filtered_by_duration = ref_filter_stats.filtered_by_duration_discharge;
                ref_filter_stats.filtered_by_std = ref_filter_stats.filtered_by_std_discharge;
                ref_filter_stats.filtered_by_direction = 0; % Discharge는 방향 체크 없음
                ref_filter_stats.filtered_by_range = ref_filter_stats.filtered_by_range_discharge;
                ref_filter_stats.final_count = ref_filter_stats.final_count_discharge;
                ref_filter_stats.fail_duration = ref_filter_stats.filtered_by_duration_discharge;
                ref_filter_stats.fail_std = ref_filter_stats.filtered_by_std_discharge;
                ref_filter_stats.fail_direction = 0;
                ref_filter_stats.fail_range = ref_filter_stats.filtered_by_range_discharge;
            else % 'Both'
                ref_filter_stats.total_candidates = ref_filter_stats.total_candidates_charge + ref_filter_stats.total_candidates_discharge;
                ref_filter_stats.filtered_by_duration = ref_filter_stats.filtered_by_duration_charge + ref_filter_stats.filtered_by_duration_discharge;
                ref_filter_stats.filtered_by_std = ref_filter_stats.filtered_by_std_charge + ref_filter_stats.filtered_by_std_discharge;
                ref_filter_stats.filtered_by_direction = ref_filter_stats.filtered_by_direction_charge;
                ref_filter_stats.filtered_by_range = ref_filter_stats.filtered_by_range_charge + ref_filter_stats.filtered_by_range_discharge;
                ref_filter_stats.final_count = ref_filter_stats.final_count_charge + ref_filter_stats.final_count_discharge;
                ref_filter_stats.fail_duration = ref_filter_stats.filtered_by_duration;
                ref_filter_stats.fail_std = ref_filter_stats.filtered_by_std;
                ref_filter_stats.fail_direction = ref_filter_stats.filtered_by_direction;
                ref_filter_stats.fail_range = ref_filter_stats.filtered_by_range;
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
            
            % 기준 채널 rejected 이벤트 저장
            for r = 1:length(ref_rejected)
                rej_evt = ref_rejected(r);
                rejField = sprintf('reject%d', r);
                rejStructName = sprintf('%s_%s_rejected', ref_channel_name, rej_evt.type);
                
                if ~isfield(rejectedEvents, rejStructName)
                    rejectedEvents.(rejStructName) = struct();
                end
                if ~isfield(rejectedEvents.(rejStructName), socName)
                    rejectedEvents.(rejStructName).(socName) = struct();
                end
                if ~isfield(rejectedEvents.(rejStructName).(socName), profName)
                    rejectedEvents.(rejStructName).(socName).(profName) = struct();
                end
                
                rejectedEvents.(rejStructName).(socName).(profName).(rejField).t = rej_evt.t_seg;
                rejectedEvents.(rejStructName).(socName).(profName).(rejField).V = rej_evt.V_seg;
                rejectedEvents.(rejStructName).(socName).(profName).(rejField).I = rej_evt.I_seg;
                rejectedEvents.(rejStructName).(socName).(profName).(rejField).I_raw = rej_evt.I_raw_seg;
                rejectedEvents.(rejStructName).(socName).(profName).(rejField).reason = rej_evt.reason;
                rejectedEvents.(rejStructName).(socName).(profName).(rejField).std_calc_indices = rej_evt.std_calc_indices;
                rejectedEvents.(rejStructName).(socName).(profName).(rejField).duration = rej_evt.duration;
            end
            
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
                debug_tag = sprintf('%s %s %s %s', cycleType, socName, profName, chName);
                [ch_events, ch_filter_stats, ch_rejected] = detect_events_in_channel(ch_I, ch_V, ch_t_seconds, current_threshold, ...
                                                      min_duration, max_I_std, max_I_range, ch_dt, debug_tag, event_type_selection);
                
                % event_type_selection에 맞는 통계만 선택
                if strcmp(event_type_selection, 'Charge')
                    ch_filter_stats.total_candidates = ch_filter_stats.total_candidates_charge;
                    ch_filter_stats.filtered_by_duration = ch_filter_stats.filtered_by_duration_charge;
                    ch_filter_stats.filtered_by_std = ch_filter_stats.filtered_by_std_charge;
                    ch_filter_stats.filtered_by_direction = ch_filter_stats.filtered_by_direction_charge;
                    ch_filter_stats.filtered_by_range = ch_filter_stats.filtered_by_range_charge;
                    ch_filter_stats.final_count = ch_filter_stats.final_count_charge;
                    ch_filter_stats.fail_duration = ch_filter_stats.filtered_by_duration_charge;
                    ch_filter_stats.fail_std = ch_filter_stats.filtered_by_std_charge;
                    ch_filter_stats.fail_direction = ch_filter_stats.filtered_by_direction_charge;
                    ch_filter_stats.fail_range = ch_filter_stats.filtered_by_range_charge;
                elseif strcmp(event_type_selection, 'Discharge')
                    ch_filter_stats.total_candidates = ch_filter_stats.total_candidates_discharge;
                    ch_filter_stats.filtered_by_duration = ch_filter_stats.filtered_by_duration_discharge;
                    ch_filter_stats.filtered_by_std = ch_filter_stats.filtered_by_std_discharge;
                    ch_filter_stats.filtered_by_direction = 0; % Discharge는 방향 체크 없음
                    ch_filter_stats.filtered_by_range = ch_filter_stats.filtered_by_range_discharge;
                    ch_filter_stats.final_count = ch_filter_stats.final_count_discharge;
                    ch_filter_stats.fail_duration = ch_filter_stats.filtered_by_duration_discharge;
                    ch_filter_stats.fail_std = ch_filter_stats.filtered_by_std_discharge;
                    ch_filter_stats.fail_direction = 0;
                    ch_filter_stats.fail_range = ch_filter_stats.filtered_by_range_discharge;
                else % 'Both'
                    ch_filter_stats.total_candidates = ch_filter_stats.total_candidates_charge + ch_filter_stats.total_candidates_discharge;
                    ch_filter_stats.filtered_by_duration = ch_filter_stats.filtered_by_duration_charge + ch_filter_stats.filtered_by_duration_discharge;
                    ch_filter_stats.filtered_by_std = ch_filter_stats.filtered_by_std_charge + ch_filter_stats.filtered_by_std_discharge;
                    ch_filter_stats.filtered_by_direction = ch_filter_stats.filtered_by_direction_charge;
                    ch_filter_stats.filtered_by_range = ch_filter_stats.filtered_by_range_charge + ch_filter_stats.filtered_by_range_discharge;
                    ch_filter_stats.final_count = ch_filter_stats.final_count_charge + ch_filter_stats.final_count_discharge;
                    ch_filter_stats.fail_duration = ch_filter_stats.filtered_by_duration;
                    ch_filter_stats.fail_std = ch_filter_stats.filtered_by_std;
                    ch_filter_stats.fail_direction = ch_filter_stats.filtered_by_direction;
                    ch_filter_stats.fail_range = ch_filter_stats.filtered_by_range;
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
                
                % 다른 채널 rejected 이벤트 저장
                for r = 1:length(ch_rejected)
                    rej_evt = ch_rejected(r);
                    rejField = sprintf('reject%d', r);
                    rejStructName = sprintf('%s_%s_rejected', chName, rej_evt.type);
                    
                    if ~isfield(rejectedEvents, rejStructName)
                        rejectedEvents.(rejStructName) = struct();
                    end
                    if ~isfield(rejectedEvents.(rejStructName), socName)
                        rejectedEvents.(rejStructName).(socName) = struct();
                    end
                    if ~isfield(rejectedEvents.(rejStructName).(socName), profName)
                        rejectedEvents.(rejStructName).(socName).(profName) = struct();
                    end
                    
                    rejectedEvents.(rejStructName).(socName).(profName).(rejField).t = rej_evt.t_seg;
                    rejectedEvents.(rejStructName).(socName).(profName).(rejField).V = rej_evt.V_seg;
                    rejectedEvents.(rejStructName).(socName).(profName).(rejField).I = rej_evt.I_seg;
                    rejectedEvents.(rejStructName).(socName).(profName).(rejField).I_raw = rej_evt.I_raw_seg;
                    rejectedEvents.(rejStructName).(socName).(profName).(rejField).reason = rej_evt.reason;
                    rejectedEvents.(rejStructName).(socName).(profName).(rejField).std_calc_indices = rej_evt.std_calc_indices;
                    rejectedEvents.(rejStructName).(socName).(profName).(rejField).duration = rej_evt.duration;
                end
                
                if isempty(ch_events)
                    continue;
                end
                
                % 모든 채널 이벤트 저장 (매칭 없이 채널 자체 기준)
                for evt_num = 1:length(ch_events)
                    ch_evt = ch_events(evt_num);
                    evtField = sprintf('event%d', evt_num);
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
                
                % 디버깅용: 채널 이벤트 카운트 (선택 타입 이벤트만)
                match_ch = regexp(chName, '^(ch\d+)', 'tokens', 'once');
                if ~isempty(match_ch)
                    chNumOnly = match_ch{1};
                else
                    chNumOnly = chName;
                end
                
                if strcmp(event_type_selection, 'Both')
                    selected_count = length(ch_events);
                else
                    selected_count = sum(strcmp({ch_events.type}, event_type_selection));
                end
                
                if selected_count > 0
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
                        eventCounts.(profName).(chNumOnly).(validCycleName) + selected_count;
                end
            end
        end
    end
    
    % 저장
    savePath = fullfile(outputDir, sprintf('Lab_DC_Events_Raw_%s.mat', cycleType));
    parsave(savePath, rawEvents, rejectedEvents);
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
                dcStatsStruct.(profName).filtered_by_range = 0;
                dcStatsStruct.(profName).final_count = 0;
                dcStatsStruct.(profName).fail_duration = 0;
                dcStatsStruct.(profName).fail_std = 0;
                dcStatsStruct.(profName).fail_direction = 0;
                dcStatsStruct.(profName).fail_range = 0;
            end
            
            % 통계 합산
            dcStatsStruct.(profName).total_candidates = dcStatsStruct.(profName).total_candidates + stats.total_candidates;
            dcStatsStruct.(profName).filtered_by_duration = dcStatsStruct.(profName).filtered_by_duration + stats.filtered_by_duration;
            dcStatsStruct.(profName).filtered_by_std = dcStatsStruct.(profName).filtered_by_std + stats.filtered_by_std;
            dcStatsStruct.(profName).filtered_by_direction = dcStatsStruct.(profName).filtered_by_direction + stats.filtered_by_direction;
            if isfield(stats, 'filtered_by_range')
                dcStatsStruct.(profName).filtered_by_range = dcStatsStruct.(profName).filtered_by_range + stats.filtered_by_range;
            end
            dcStatsStruct.(profName).final_count = dcStatsStruct.(profName).final_count + stats.final_count;
            if isfield(stats, 'fail_duration')
                dcStatsStruct.(profName).fail_duration = dcStatsStruct.(profName).fail_duration + stats.fail_duration;
                dcStatsStruct.(profName).fail_std = dcStatsStruct.(profName).fail_std + stats.fail_std;
                dcStatsStruct.(profName).fail_direction = dcStatsStruct.(profName).fail_direction + stats.fail_direction;
                if isfield(stats, 'fail_range')
                    dcStatsStruct.(profName).fail_range = dcStatsStruct.(profName).fail_range + stats.fail_range;
                end
            end
        end
    end
    
    dcNames = fieldnames(dcStatsStruct);
    if ~isempty(dcNames)
        % DC 목록 정렬
        dcNamesSorted = sort(dcNames);
        
        % 표 헤더 (엑셀 복사용 쉼표 구분)
        fprintf('\n필터링 통계 (DC별 집계) - 엑셀 복사용 (쉼표 구분):\n');
        fprintf('DC,총후보,제거:시간,제거:표준편차,제거:범위,제거:방향,최종\n');
        
        % 각 DC별 통계 출력 (쉼표 구분)
        for d = 1:length(dcNamesSorted)
            dcName = dcNamesSorted{d};
            dcStats = dcStatsStruct.(dcName);
            
            fprintf('%s,%d,%d,%d,%d,%d,%d\n', ...
                dcName, ...
                dcStats.total_candidates, ...
                dcStats.filtered_by_duration, ...
                dcStats.filtered_by_std, ...
                dcStats.filtered_by_range, ...
                dcStats.filtered_by_direction, ...
                dcStats.final_count);
        end
        
        % 전체 합계
        totalAll = struct();
        totalAll.total_candidates = 0;
        totalAll.filtered_by_duration = 0;
        totalAll.filtered_by_std = 0;
        totalAll.filtered_by_direction = 0;
        totalAll.filtered_by_range = 0;
        totalAll.final_count = 0;
        totalAll.fail_duration = 0;
        totalAll.fail_std = 0;
        totalAll.fail_direction = 0;
        totalAll.fail_range = 0;
        
        for d = 1:length(dcNamesSorted)
            dcName = dcNamesSorted{d};
            dcStats = dcStatsStruct.(dcName);
            totalAll.total_candidates = totalAll.total_candidates + dcStats.total_candidates;
            totalAll.filtered_by_duration = totalAll.filtered_by_duration + dcStats.filtered_by_duration;
            totalAll.filtered_by_std = totalAll.filtered_by_std + dcStats.filtered_by_std;
            totalAll.filtered_by_direction = totalAll.filtered_by_direction + dcStats.filtered_by_direction;
            totalAll.filtered_by_range = totalAll.filtered_by_range + dcStats.filtered_by_range;
            totalAll.final_count = totalAll.final_count + dcStats.final_count;
            totalAll.fail_duration = totalAll.fail_duration + dcStats.fail_duration;
            totalAll.fail_std = totalAll.fail_std + dcStats.fail_std;
            totalAll.fail_direction = totalAll.fail_direction + dcStats.fail_direction;
            totalAll.fail_range = totalAll.fail_range + dcStats.fail_range;
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
            fprintf('%-40s: %8d (%5.1f%%)\n', '제거: 지속시간 부족', totalAll.fail_duration, totalAll.fail_duration/totalAll.total_candidates*100);
            fprintf('%-40s: %8d (%5.1f%%)\n', '제거: 전류 표준편차 초과', totalAll.fail_std, totalAll.fail_std/totalAll.total_candidates*100);
            fprintf('%-40s: %8d (%5.1f%%)\n', '제거: 전류 범위 초과', totalAll.fail_range, totalAll.fail_range/totalAll.total_candidates*100);
            fprintf('%-40s: %8d (%5.1f%%)\n', '제거: 방향 불일치', totalAll.fail_direction, totalAll.fail_direction/totalAll.total_candidates*100);
            fprintf('%-40s: %8d (%5.1f%%)\n', '최종 채택 이벤트 수', totalAll.final_count, totalAll.final_count/totalAll.total_candidates*100);
        end

        % 채널별 필터링 통계 테이블 저장 (txt)
        channelSet = {};
        for s = 1:length(statsFields)
            stats = filterStatsCollection.(statsFields{s});
            if isfield(stats, 'channel')
                channelSet{end+1} = stats.channel; %#ok<AGROW>
            end
        end
        channelSet = unique(channelSet);

        for cIdx = 1:length(channelSet)
            channelName = channelSet{cIdx};
            statsDir = fullfile(outputDir, 'tables', channelName);
            if ~exist(statsDir, 'dir')
                mkdir(statsDir);
            end

            if strcmp(event_type_selection, 'Both')
                dcStatsCharge = struct();
                dcStatsDischarge = struct();
                for s = 1:length(statsFields)
                    statKey = statsFields{s};
                    stats = filterStatsCollection.(statKey);
                    if ~isfield(stats, 'profile')
                        continue;
                    end
                    if ~isfield(stats, 'channel') || ~strcmpi(stats.channel, channelName)
                        continue;
                    end
                    profName = stats.profile;
                    if ~isfield(dcStatsCharge, profName)
                        dcStatsCharge.(profName) = struct('total_candidates',0,'filtered_by_duration',0, ...
                            'filtered_by_std',0,'filtered_by_range',0,'filtered_by_direction',0,'final_count',0);
                    end
                    if ~isfield(dcStatsDischarge, profName)
                        dcStatsDischarge.(profName) = struct('total_candidates',0,'filtered_by_duration',0, ...
                            'filtered_by_std',0,'filtered_by_range',0,'filtered_by_direction',0,'final_count',0);
                    end
                    if isfield(stats, 'total_candidates_charge')
                        dcStatsCharge.(profName).total_candidates = dcStatsCharge.(profName).total_candidates + stats.total_candidates_charge;
                        dcStatsCharge.(profName).filtered_by_duration = dcStatsCharge.(profName).filtered_by_duration + stats.filtered_by_duration_charge;
                        dcStatsCharge.(profName).filtered_by_std = dcStatsCharge.(profName).filtered_by_std + stats.filtered_by_std_charge;
                        dcStatsCharge.(profName).filtered_by_direction = dcStatsCharge.(profName).filtered_by_direction + stats.filtered_by_direction_charge;
                        if isfield(stats, 'filtered_by_range_charge')
                            dcStatsCharge.(profName).filtered_by_range = dcStatsCharge.(profName).filtered_by_range + stats.filtered_by_range_charge;
                        end
                        dcStatsCharge.(profName).final_count = dcStatsCharge.(profName).final_count + stats.final_count_charge;
                    end
                    if isfield(stats, 'total_candidates_discharge')
                        dcStatsDischarge.(profName).total_candidates = dcStatsDischarge.(profName).total_candidates + stats.total_candidates_discharge;
                        dcStatsDischarge.(profName).filtered_by_duration = dcStatsDischarge.(profName).filtered_by_duration + stats.filtered_by_duration_discharge;
                        dcStatsDischarge.(profName).filtered_by_std = dcStatsDischarge.(profName).filtered_by_std + stats.filtered_by_std_discharge;
                        if isfield(stats, 'filtered_by_range_discharge')
                            dcStatsDischarge.(profName).filtered_by_range = dcStatsDischarge.(profName).filtered_by_range + stats.filtered_by_range_discharge;
                        end
                        if isfield(stats, 'filtered_by_direction_discharge')
                            dcStatsDischarge.(profName).filtered_by_direction = dcStatsDischarge.(profName).filtered_by_direction + stats.filtered_by_direction_discharge;
                        end
                        dcStatsDischarge.(profName).final_count = dcStatsDischarge.(profName).final_count + stats.final_count_discharge;
                    end
                end

                dcNames = fieldnames(dcStatsCharge);
                if ~isempty(dcNames)
                    dcNamesSorted = sort(dcNames);
                    statsPath = fullfile(statsDir, 'filtering_stats_charge.txt');
                    fid = fopen(statsPath, 'w');
                    fprintf(fid, 'DC,총후보(Idle->Charge),제거:시간,제거:표준편차,제거:범위,제거:방향,최종\n');
                    totalRef = struct('total_candidates',0,'filtered_by_duration',0,'filtered_by_std',0, ...
                                      'filtered_by_range',0,'filtered_by_direction',0,'final_count',0);
                    for d = 1:length(dcNamesSorted)
                        dcName = dcNamesSorted{d};
                        dcStats = dcStatsCharge.(dcName);
                        fprintf(fid, '%s,%d,%d,%d,%d,%d,%d\n', ...
                            dcName, dcStats.total_candidates, dcStats.filtered_by_duration, ...
                            dcStats.filtered_by_std, dcStats.filtered_by_range, ...
                            dcStats.filtered_by_direction, dcStats.final_count);
                        totalRef.total_candidates = totalRef.total_candidates + dcStats.total_candidates;
                        totalRef.filtered_by_duration = totalRef.filtered_by_duration + dcStats.filtered_by_duration;
                        totalRef.filtered_by_std = totalRef.filtered_by_std + dcStats.filtered_by_std;
                        totalRef.filtered_by_range = totalRef.filtered_by_range + dcStats.filtered_by_range;
                        totalRef.filtered_by_direction = totalRef.filtered_by_direction + dcStats.filtered_by_direction;
                        totalRef.final_count = totalRef.final_count + dcStats.final_count;
                    end
                    fprintf(fid, '전체,%d,%d,%d,%d,%d,%d\n', ...
                        totalRef.total_candidates, totalRef.filtered_by_duration, ...
                        totalRef.filtered_by_std, totalRef.filtered_by_range, ...
                        totalRef.filtered_by_direction, totalRef.final_count);
                    fclose(fid);
                end

                dcNames = fieldnames(dcStatsDischarge);
                if ~isempty(dcNames)
                    dcNamesSorted = sort(dcNames);
                    statsPath = fullfile(statsDir, 'filtering_stats_discharge.txt');
                    fid = fopen(statsPath, 'w');
                    fprintf(fid, 'DC,총후보(Idle->Discharge),제거:시간,제거:표준편차,제거:범위,제거:방향,최종\n');
                    totalRef = struct('total_candidates',0,'filtered_by_duration',0,'filtered_by_std',0, ...
                                      'filtered_by_range',0,'filtered_by_direction',0,'final_count',0);
                    for d = 1:length(dcNamesSorted)
                        dcName = dcNamesSorted{d};
                        dcStats = dcStatsDischarge.(dcName);
                        fprintf(fid, '%s,%d,%d,%d,%d,%d,%d\n', ...
                            dcName, dcStats.total_candidates, dcStats.filtered_by_duration, ...
                            dcStats.filtered_by_std, dcStats.filtered_by_range, ...
                            dcStats.filtered_by_direction, dcStats.final_count);
                        totalRef.total_candidates = totalRef.total_candidates + dcStats.total_candidates;
                        totalRef.filtered_by_duration = totalRef.filtered_by_duration + dcStats.filtered_by_duration;
                        totalRef.filtered_by_std = totalRef.filtered_by_std + dcStats.filtered_by_std;
                        totalRef.filtered_by_range = totalRef.filtered_by_range + dcStats.filtered_by_range;
                        totalRef.filtered_by_direction = totalRef.filtered_by_direction + dcStats.filtered_by_direction;
                        totalRef.final_count = totalRef.final_count + dcStats.final_count;
                    end
                    fprintf(fid, '전체,%d,%d,%d,%d,%d,%d\n', ...
                        totalRef.total_candidates, totalRef.filtered_by_duration, ...
                        totalRef.filtered_by_std, totalRef.filtered_by_range, ...
                        totalRef.filtered_by_direction, totalRef.final_count);
                    fclose(fid);
                end
            else
                dcStatsStruct = struct();
                for s = 1:length(statsFields)
                    statKey = statsFields{s};
                    stats = filterStatsCollection.(statKey);
                    if ~isfield(stats, 'profile') || ~isfield(stats, 'total_candidates')
                        continue;
                    end
                    if ~isfield(stats, 'channel') || ~strcmpi(stats.channel, channelName)
                        continue;
                    end
                    profName = stats.profile;
                    if ~isfield(dcStatsStruct, profName)
                        dcStatsStruct.(profName) = struct();
                        dcStatsStruct.(profName).total_candidates = 0;
                        dcStatsStruct.(profName).filtered_by_duration = 0;
                        dcStatsStruct.(profName).filtered_by_std = 0;
                        dcStatsStruct.(profName).filtered_by_direction = 0;
                        dcStatsStruct.(profName).filtered_by_range = 0;
                        dcStatsStruct.(profName).final_count = 0;
                    end
                    dcStatsStruct.(profName).total_candidates = dcStatsStruct.(profName).total_candidates + stats.total_candidates;
                    dcStatsStruct.(profName).filtered_by_duration = dcStatsStruct.(profName).filtered_by_duration + stats.filtered_by_duration;
                    dcStatsStruct.(profName).filtered_by_std = dcStatsStruct.(profName).filtered_by_std + stats.filtered_by_std;
                    dcStatsStruct.(profName).filtered_by_direction = dcStatsStruct.(profName).filtered_by_direction + stats.filtered_by_direction;
                    if isfield(stats, 'filtered_by_range')
                        dcStatsStruct.(profName).filtered_by_range = dcStatsStruct.(profName).filtered_by_range + stats.filtered_by_range;
                    end
                    dcStatsStruct.(profName).final_count = dcStatsStruct.(profName).final_count + stats.final_count;
                end

                dcNames = fieldnames(dcStatsStruct);
                if ~isempty(dcNames)
                    dcNamesSorted = sort(dcNames);
                    statsPath = fullfile(statsDir, sprintf('filtering_stats_%s.txt', lower(event_type_selection)));
                    fid = fopen(statsPath, 'w');
                    if strcmp(event_type_selection, 'Charge')
                        headerLabel = 'Idle->Charge';
                    else
                        headerLabel = 'Idle->Discharge';
                    end
                    fprintf(fid, 'DC,총후보(%s),제거:시간,제거:표준편차,제거:범위,제거:방향,최종\n', headerLabel);
                    totalRef = struct('total_candidates',0,'filtered_by_duration',0,'filtered_by_std',0, ...
                                      'filtered_by_range',0,'filtered_by_direction',0,'final_count',0);
                    for d = 1:length(dcNamesSorted)
                        dcName = dcNamesSorted{d};
                        dcStats = dcStatsStruct.(dcName);
                        fprintf(fid, '%s,%d,%d,%d,%d,%d,%d\n', ...
                            dcName, dcStats.total_candidates, dcStats.filtered_by_duration, ...
                            dcStats.filtered_by_std, dcStats.filtered_by_range, ...
                            dcStats.filtered_by_direction, dcStats.final_count);
                        totalRef.total_candidates = totalRef.total_candidates + dcStats.total_candidates;
                        totalRef.filtered_by_duration = totalRef.filtered_by_duration + dcStats.filtered_by_duration;
                        totalRef.filtered_by_std = totalRef.filtered_by_std + dcStats.filtered_by_std;
                        totalRef.filtered_by_range = totalRef.filtered_by_range + dcStats.filtered_by_range;
                        totalRef.filtered_by_direction = totalRef.filtered_by_direction + dcStats.filtered_by_direction;
                        totalRef.final_count = totalRef.final_count + dcStats.final_count;
                    end
                    fprintf(fid, '전체,%d,%d,%d,%d,%d,%d\n', ...
                        totalRef.total_candidates, totalRef.filtered_by_duration, ...
                        totalRef.filtered_by_std, totalRef.filtered_by_range, ...
                        totalRef.filtered_by_direction, totalRef.final_count);
                    fclose(fid);
                end
            end
        end

        % 전체 채널 통합표 저장
        statsDir = fullfile(outputDir, 'tables');
        if ~exist(statsDir, 'dir')
            mkdir(statsDir);
        end

        if strcmp(event_type_selection, 'Both')
            combinedCharge = struct();
            combinedDischarge = struct();
            for s = 1:length(statsFields)
                statKey = statsFields{s};
                stats = filterStatsCollection.(statKey);
                if ~isfield(stats, 'profile')
                    continue;
                end
                profName = stats.profile;
                if ~isfield(combinedCharge, profName)
                    combinedCharge.(profName) = struct('total_candidates',0,'filtered_by_duration',0, ...
                        'filtered_by_std',0,'filtered_by_range',0,'filtered_by_direction',0,'final_count',0);
                end
                if ~isfield(combinedDischarge, profName)
                    combinedDischarge.(profName) = struct('total_candidates',0,'filtered_by_duration',0, ...
                        'filtered_by_std',0,'filtered_by_range',0,'filtered_by_direction',0,'final_count',0);
                end
                if isfield(stats, 'total_candidates_charge')
                    combinedCharge.(profName).total_candidates = combinedCharge.(profName).total_candidates + stats.total_candidates_charge;
                    combinedCharge.(profName).filtered_by_duration = combinedCharge.(profName).filtered_by_duration + stats.filtered_by_duration_charge;
                    combinedCharge.(profName).filtered_by_std = combinedCharge.(profName).filtered_by_std + stats.filtered_by_std_charge;
                    combinedCharge.(profName).filtered_by_direction = combinedCharge.(profName).filtered_by_direction + stats.filtered_by_direction_charge;
                    if isfield(stats, 'filtered_by_range_charge')
                        combinedCharge.(profName).filtered_by_range = combinedCharge.(profName).filtered_by_range + stats.filtered_by_range_charge;
                    end
                    combinedCharge.(profName).final_count = combinedCharge.(profName).final_count + stats.final_count_charge;
                end
                if isfield(stats, 'total_candidates_discharge')
                    combinedDischarge.(profName).total_candidates = combinedDischarge.(profName).total_candidates + stats.total_candidates_discharge;
                    combinedDischarge.(profName).filtered_by_duration = combinedDischarge.(profName).filtered_by_duration + stats.filtered_by_duration_discharge;
                    combinedDischarge.(profName).filtered_by_std = combinedDischarge.(profName).filtered_by_std + stats.filtered_by_std_discharge;
                    if isfield(stats, 'filtered_by_range_discharge')
                        combinedDischarge.(profName).filtered_by_range = combinedDischarge.(profName).filtered_by_range + stats.filtered_by_range_discharge;
                    end
                    if isfield(stats, 'filtered_by_direction_discharge')
                        combinedDischarge.(profName).filtered_by_direction = combinedDischarge.(profName).filtered_by_direction + stats.filtered_by_direction_discharge;
                    end
                    combinedDischarge.(profName).final_count = combinedDischarge.(profName).final_count + stats.final_count_discharge;
                end
            end

            dcNames = fieldnames(combinedCharge);
            if ~isempty(dcNames)
                dcNamesSorted = sort(dcNames);
                statsPath = fullfile(statsDir, 'filtering_stats_allchannels_charge.txt');
                fid = fopen(statsPath, 'w');
                fprintf(fid, 'DC,총후보(Idle->Charge),제거:시간,제거:표준편차,제거:범위,제거:방향,최종\n');
                totalRef = struct('total_candidates',0,'filtered_by_duration',0,'filtered_by_std',0, ...
                                  'filtered_by_range',0,'filtered_by_direction',0,'final_count',0);
                for d = 1:length(dcNamesSorted)
                    dcName = dcNamesSorted{d};
                    dcStats = combinedCharge.(dcName);
                    fprintf(fid, '%s,%d,%d,%d,%d,%d,%d\n', ...
                        dcName, dcStats.total_candidates, dcStats.filtered_by_duration, ...
                        dcStats.filtered_by_std, dcStats.filtered_by_range, ...
                        dcStats.filtered_by_direction, dcStats.final_count);
                    totalRef.total_candidates = totalRef.total_candidates + dcStats.total_candidates;
                    totalRef.filtered_by_duration = totalRef.filtered_by_duration + dcStats.filtered_by_duration;
                    totalRef.filtered_by_std = totalRef.filtered_by_std + dcStats.filtered_by_std;
                    totalRef.filtered_by_range = totalRef.filtered_by_range + dcStats.filtered_by_range;
                    totalRef.filtered_by_direction = totalRef.filtered_by_direction + dcStats.filtered_by_direction;
                    totalRef.final_count = totalRef.final_count + dcStats.final_count;
                end
                fprintf(fid, '전체,%d,%d,%d,%d,%d,%d\n', ...
                    totalRef.total_candidates, totalRef.filtered_by_duration, ...
                    totalRef.filtered_by_std, totalRef.filtered_by_range, ...
                    totalRef.filtered_by_direction, totalRef.final_count);
                fclose(fid);
            end

            dcNames = fieldnames(combinedDischarge);
            if ~isempty(dcNames)
                dcNamesSorted = sort(dcNames);
                statsPath = fullfile(statsDir, 'filtering_stats_allchannels_discharge.txt');
                fid = fopen(statsPath, 'w');
                fprintf(fid, 'DC,총후보(Idle->Discharge),제거:시간,제거:표준편차,제거:범위,제거:방향,최종\n');
                totalRef = struct('total_candidates',0,'filtered_by_duration',0,'filtered_by_std',0, ...
                                  'filtered_by_range',0,'filtered_by_direction',0,'final_count',0);
                for d = 1:length(dcNamesSorted)
                    dcName = dcNamesSorted{d};
                    dcStats = combinedDischarge.(dcName);
                    fprintf(fid, '%s,%d,%d,%d,%d,%d,%d\n', ...
                        dcName, dcStats.total_candidates, dcStats.filtered_by_duration, ...
                        dcStats.filtered_by_std, dcStats.filtered_by_range, ...
                        dcStats.filtered_by_direction, dcStats.final_count);
                    totalRef.total_candidates = totalRef.total_candidates + dcStats.total_candidates;
                    totalRef.filtered_by_duration = totalRef.filtered_by_duration + dcStats.filtered_by_duration;
                    totalRef.filtered_by_std = totalRef.filtered_by_std + dcStats.filtered_by_std;
                    totalRef.filtered_by_range = totalRef.filtered_by_range + dcStats.filtered_by_range;
                    totalRef.filtered_by_direction = totalRef.filtered_by_direction + dcStats.filtered_by_direction;
                    totalRef.final_count = totalRef.final_count + dcStats.final_count;
                end
                fprintf(fid, '전체,%d,%d,%d,%d,%d,%d\n', ...
                    totalRef.total_candidates, totalRef.filtered_by_duration, ...
                    totalRef.filtered_by_std, totalRef.filtered_by_range, ...
                    totalRef.filtered_by_direction, totalRef.final_count);
                fclose(fid);
            end
        else
            combinedStats = struct();
            for s = 1:length(statsFields)
                statKey = statsFields{s};
                stats = filterStatsCollection.(statKey);
                if ~isfield(stats, 'profile') || ~isfield(stats, 'total_candidates')
                    continue;
                end
                profName = stats.profile;
                if ~isfield(combinedStats, profName)
                    combinedStats.(profName) = struct('total_candidates',0,'filtered_by_duration',0, ...
                        'filtered_by_std',0,'filtered_by_range',0,'filtered_by_direction',0,'final_count',0);
                end
                combinedStats.(profName).total_candidates = combinedStats.(profName).total_candidates + stats.total_candidates;
                combinedStats.(profName).filtered_by_duration = combinedStats.(profName).filtered_by_duration + stats.filtered_by_duration;
                combinedStats.(profName).filtered_by_std = combinedStats.(profName).filtered_by_std + stats.filtered_by_std;
                combinedStats.(profName).filtered_by_direction = combinedStats.(profName).filtered_by_direction + stats.filtered_by_direction;
                if isfield(stats, 'filtered_by_range')
                    combinedStats.(profName).filtered_by_range = combinedStats.(profName).filtered_by_range + stats.filtered_by_range;
                end
                combinedStats.(profName).final_count = combinedStats.(profName).final_count + stats.final_count;
            end

            dcNames = fieldnames(combinedStats);
            if ~isempty(dcNames)
                dcNamesSorted = sort(dcNames);
                statsPath = fullfile(statsDir, sprintf('filtering_stats_allchannels_%s.txt', lower(event_type_selection)));
                fid = fopen(statsPath, 'w');
                if strcmp(event_type_selection, 'Charge')
                    headerLabel = 'Idle->Charge';
                else
                    headerLabel = 'Idle->Discharge';
                end
                fprintf(fid, 'DC,총후보(%s),제거:시간,제거:표준편차,제거:범위,제거:방향,최종\n', headerLabel);
                totalRef = struct('total_candidates',0,'filtered_by_duration',0,'filtered_by_std',0, ...
                                  'filtered_by_range',0,'filtered_by_direction',0,'final_count',0);
                for d = 1:length(dcNamesSorted)
                    dcName = dcNamesSorted{d};
                    dcStats = combinedStats.(dcName);
                    fprintf(fid, '%s,%d,%d,%d,%d,%d,%d\n', ...
                        dcName, dcStats.total_candidates, dcStats.filtered_by_duration, ...
                        dcStats.filtered_by_std, dcStats.filtered_by_range, ...
                        dcStats.filtered_by_direction, dcStats.final_count);
                    totalRef.total_candidates = totalRef.total_candidates + dcStats.total_candidates;
                    totalRef.filtered_by_duration = totalRef.filtered_by_duration + dcStats.filtered_by_duration;
                    totalRef.filtered_by_std = totalRef.filtered_by_std + dcStats.filtered_by_std;
                    totalRef.filtered_by_range = totalRef.filtered_by_range + dcStats.filtered_by_range;
                    totalRef.filtered_by_direction = totalRef.filtered_by_direction + dcStats.filtered_by_direction;
                    totalRef.final_count = totalRef.final_count + dcStats.final_count;
                end
                fprintf(fid, '전체,%d,%d,%d,%d,%d,%d\n', ...
                    totalRef.total_candidates, totalRef.filtered_by_duration, ...
                    totalRef.filtered_by_std, totalRef.filtered_by_range, ...
                    totalRef.filtered_by_direction, totalRef.final_count);
                fclose(fid);
            end
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
            for pIdx = 1:length(profs)
                profName = profs{pIdx};
                
            % 이벤트 개수 계산 (SOC별로 개별 계산)
                events = fieldnames(rawEvents_loaded.(structName).(socName).(profName));
                eventCount = length(events);
                
            % 구조체에 저장 (SOC별로 저장)
                if ~isfield(filteredEventCounts, profName)
                    filteredEventCounts.(profName) = struct();
                end
                if ~isfield(filteredEventCounts.(profName), chNumOnly)
                    filteredEventCounts.(profName).(chNumOnly) = struct();
                end
            if ~isfield(filteredEventCounts.(profName).(chNumOnly), socName)
                filteredEventCounts.(profName).(chNumOnly).(socName) = struct();
            end
            if ~isfield(filteredEventCounts.(profName).(chNumOnly).(socName), validCycleName)
                filteredEventCounts.(profName).(chNumOnly).(socName).(validCycleName) = 0;
            end
            filteredEventCounts.(profName).(chNumOnly).(socName).(validCycleName) = ...
                filteredEventCounts.(profName).(chNumOnly).(socName).(validCycleName) + eventCount;
            end
        end
    end
end

% 테이블 생성 및 저장 (txt, 채널별 요약)
dcNames = fieldnames(filteredEventCounts);
dcOrder = {'DC1','DC2','DC3','DC4','DC5','DC6','DC7','DC8'};
dcNames = dcOrder(ismember(dcOrder, dcNames));

% txt 저장 경로 (시각화가 생략되어도 보장)
if ~exist('figuresDir', 'var') || isempty(figuresDir)
    figuresDir = fullfile(outputDir, 'figures', 'CurrentProfiles');
end
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

% 채널/사이클/SOC 목록 수집
allChannels = {};
allSocs = {};
allCyclesValid = {};
allCyclesDisplay = {};
for d = 1:length(dcNames)
    dcName = dcNames{d};
    chNames = fieldnames(filteredEventCounts.(dcName));
    allChannels = unique([allChannels; chNames]);
    for c = 1:length(chNames)
        chName = chNames{c};
        socs = fieldnames(filteredEventCounts.(dcName).(chName));
        allSocs = unique([allSocs; socs]);
        for s = 1:length(socs)
            cyclesValid = fieldnames(filteredEventCounts.(dcName).(chName).(socs{s}));
            allCyclesValid = unique([allCyclesValid; cyclesValid]);
        end
    end
end

for cy = 1:length(allCyclesValid)
    validName = allCyclesValid{cy};
    if isKey(cycleNameReverseMap, validName)
        allCyclesDisplay{cy} = cycleNameReverseMap(validName);
    else
        allCyclesDisplay{cy} = validName;
    end
end

for c = 1:length(allChannels)
    chName = allChannels{c};
    txtPath = fullfile(figuresDir, sprintf('%s_event_counts.txt', chName));
    fid = fopen(txtPath, 'w');
    fprintf(fid, 'DC,Channel,SOC');
    for cy = 1:length(allCyclesDisplay)
        fprintf(fid, ',%s', allCyclesDisplay{cy});
    end
    fprintf(fid, '\n');
    
    for d = 1:length(dcNames)
        dcName = dcNames{d};
        for s = 1:length(allSocs)
            socName = allSocs{s};
            fprintf(fid, '%s,%s,%s', dcName, chName, socName);
            for cy = 1:length(allCyclesValid)
                cycleNameValid = allCyclesValid{cy};
                countVal = 0;
                if isfield(filteredEventCounts.(dcName), chName) && ...
                   isfield(filteredEventCounts.(dcName).(chName), socName) && ...
                   isfield(filteredEventCounts.(dcName).(chName).(socName), cycleNameValid)
                    countVal = filteredEventCounts.(dcName).(chName).(socName).(cycleNameValid);
                end
                fprintf(fid, ',%d', countVal);
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
end

% 리포트용 전체 조합 테이블 저장
summaryPath = fullfile(outputDir, 'event_counts_table.txt');
fid = fopen(summaryPath, 'w');
fprintf(fid, 'DC,Cycle,SOC,Channel,Count\n');
for d = 1:length(dcNames)
    dcName = dcNames{d};
    for c = 1:length(allChannels)
        chName = allChannels{c};
        for s = 1:length(allSocs)
            socName = allSocs{s};
            for cy = 1:length(allCyclesValid)
                cycleNameValid = allCyclesValid{cy};
                cycleDisplay = allCyclesDisplay{cy};
                countVal = 0;
                if isfield(filteredEventCounts.(dcName), chName) && ...
                   isfield(filteredEventCounts.(dcName).(chName), socName) && ...
                   isfield(filteredEventCounts.(dcName).(chName).(socName), cycleNameValid)
                    countVal = filteredEventCounts.(dcName).(chName).(socName).(cycleNameValid);
                end
                fprintf(fid, '%s,%s,%s,%s,%d\n', dcName, cycleDisplay, socName, chName, countVal);
            end
        end
    end
end
fclose(fid);

% 사이클별 이벤트 개수 일치 여부 및 보정 검사 리포트
reportPath = fullfile(outputDir, 'missing_event_report.txt');
fid = fopen(reportPath, 'w');
fprintf(fid, 'DC,SOC,Channel,BaseCycle,BaseCount,Cycle,EventIdx,Detected,Reason,Duration,Std,Range,DirectionRatio\n');

cyclesWanted = {'0cyc','200cyc','400cyc','600cyc','800cyc'};
rt = readtable(summaryPath, 'Delimiter', ',', 'ReadVariableNames', true, 'VariableNamingRule', 'preserve');

if ~isempty(rt)
    combos = unique(rt(:, {'DC','SOC','Channel'}), 'rows');
    for r = 1:height(combos)
        dcName = combos.('DC'){r};
        socName = combos.('SOC'){r};
        chName = combos.('Channel'){r};
        
        % 기준 사이클(0cyc) 이벤트 개수
        baseRows = rt(strcmp(rt.('DC'), dcName) & strcmp(rt.('SOC'), socName) & ...
                      strcmp(rt.('Channel'), chName) & strcmp(rt.('Cycle'), '0cyc'), :);
        if isempty(baseRows)
            continue;
        end
        baseCount = baseRows.Count(1);
        if baseCount == 0
            continue;
        end
        
        % 기준 이벤트 구간 로드
        basePath = fullfile(outputDir, 'Lab_DC_Events_Raw_0cyc.mat');
        if ~exist(basePath, 'file')
            continue;
        end
        baseData = load(basePath, 'rawEvents');
        if ~isfield(baseData, 'rawEvents')
            continue;
        end
        rawEvents_loaded = baseData.rawEvents;
        
        structName = sprintf('%s_%s', chName, event_type_selection);
        if ~isfield(rawEvents_loaded, structName) || ...
           ~isfield(rawEvents_loaded.(structName), socName) || ...
           ~isfield(rawEvents_loaded.(structName).(socName), dcName)
            continue;
        end
        baseEvents = rawEvents_loaded.(structName).(socName).(dcName);
        baseEventFields = fieldnames(baseEvents);
        
        for cy = 1:length(cyclesWanted)
            cycleName = cyclesWanted{cy};
            rows = rt(strcmp(rt.('DC'), dcName) & strcmp(rt.('SOC'), socName) & ...
                      strcmp(rt.('Channel'), chName) & strcmp(rt.('Cycle'), cycleName), :);
            detectedCount = 0;
            if ~isempty(rows)
                detectedCount = rows.Count(1);
            end
            
            if detectedCount == baseCount
                continue;
            end
            
            % 해당 사이클 원시 데이터 로드
            dataPath = fullfile(dataDir, sprintf('parsedDriveCycle_%s_filtered.mat', cycleName));
            if ~exist(dataPath, 'file')
                continue;
            end
            dataStruct = load(dataPath);
            varName = sprintf('parsedDriveCycle_%s', cycleName);
            if ~isfield(dataStruct, varName)
                continue;
            end
            data_var = dataStruct.(varName);
            if ~isfield(data_var, chName) || ~isfield(data_var.(chName), socName) || ...
               ~isfield(data_var.(chName).(socName), dcName)
                continue;
            end
            
            seg = data_var.(chName).(socName).(dcName);
            if ~isfield(seg, 'I') || ~isfield(seg, 'V') || ~isfield(seg, 't')
                continue;
            end
            I = seg.I;
            V = seg.V;
            t = seg.t;
            if isa(t, 'duration')
                t_seconds = seconds(t);
            else
                t_seconds = t;
            end
            
            for e = 1:length(baseEventFields)
                evt = baseEvents.(baseEventFields{e});
                if ~isfield(evt, 't') || isempty(evt.t)
                    continue;
                end
                t0 = evt.t(1) - 2;
                t1 = evt.t(end) + 2;
                
                idx = find(t_seconds >= t0 & t_seconds <= t1);
                if isempty(idx)
                    fprintf(fid, '%s,%s,%s,0cyc,%d,%s,%d,0,%s,0,0,0,0\n', ...
                        dcName, socName, chName, baseCount, cycleName, e, 'no_window_data');
                    continue;
                end
                
                I_seg = I(idx);
                t_seg = t_seconds(idx);
                
                [pass, reason, duration, stdVal, rangeVal, dirRatio] = evaluate_event_window( ...
                    I_seg, t_seg, event_type_selection, min_duration, max_I_std, max_I_range, ...
                    std_calc_mode, trim_sec, fixed_duration);
                
                fprintf(fid, '%s,%s,%s,0cyc,%d,%s,%d,%d,%s,%.3f,%.3f,%.3f,%.3f\n', ...
                    dcName, socName, chName, baseCount, cycleName, e, pass, reason, ...
                    duration, stdVal, rangeVal, dirRatio);
            end
        end
    end
end
fclose(fid);

fprintf('\n=== Event Detection Complete ===\n');
fprintf('All results saved to: %s\n', outputDir);

%% DC별 전류 프로파일 시각화 (SOC90/70/50)
fprintf('\n=== Creating Current Profile Visualization (SOC90/70/50) ===\n');

% 기준 채널에서 데이터 추출
ref_channel_match = regexp(reference_channel, '^(ch\d+)', 'tokens', 'once');

figuresDir = fullfile(outputDir, 'figures', 'CurrentProfiles');
if ~isempty(ref_channel_match)
    figuresDir = fullfile(figuresDir, reference_channel);
end
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

targetSOCs = {'SOC90', 'SOC70', 'SOC50'};  % 시각화할 SOC 목록

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
    
    for socIdx = 1:length(targetSOCs)
        targetSOC = targetSOCs{socIdx};
        if isstring(targetSOC)
            targetSOC = char(targetSOC);
        elseif iscell(targetSOC)
            targetSOC = targetSOC{1};
        end
        if ~isvarname(targetSOC)
            targetSOC = regexprep(targetSOC, '[^a-zA-Z0-9_]', '_');
            if ~isempty(targetSOC) && ~isletter(targetSOC(1))
                targetSOC = ['SOC_' targetSOC];
            end
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
            
            for pIdx = 1:length(profs)
                profName = profs{pIdx};
                
                if ~isfield(rawEvents_loaded.(structName).(targetSOC), profName)
                    continue;
                end
                
                events = fieldnames(rawEvents_loaded.(structName).(targetSOC).(profName));
                
                if isempty(events)
                    continue;
                end
                
                % 모든 이벤트의 데이터 사용
                if ~isfield(allCycleData, profName)
                    allCycleData.(profName) = struct();
                end
                
                % cycleType을 유효한 필드명으로 변환 (숫자로 시작하면 cyc_ 접두사 추가)
                if ~isempty(regexp(cycleType, '^\d', 'once'))
                    validCycleFieldName = ['cyc_' cycleType];
                else
                    validCycleFieldName = cycleType;
                end
                
                if ~isfield(allCycleData.(profName), validCycleFieldName)
                    allCycleData.(profName).(validCycleFieldName) = {};
                end
                
                for e = 1:length(events)
                    evtData = rawEvents_loaded.(structName).(targetSOC).(profName).(events{e});
                    if isfield(evtData, 'I') && isfield(evtData, 't')
                        allCycleData.(profName).(validCycleFieldName){end+1} = evtData;
                    end
                end
            end
        end
        
        % DC1-8 고정 서브플롯 생성
        allDCProfiles = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};
        dcProfiles = fieldnames(allCycleData);
        
        fprintf('  Creating visualization for %s events (%s)\n', eventType, targetSOC);
        
        % 서브플롯으로 DC1-8 모두 표시 (2행 4열)
        nCols = 4;
        nRows = 2;
        
        fig = figure('Name', sprintf('Current Profiles - %s Events (%s)', eventType, targetSOC), ...
                     'Position', [100, 100, 1600, 900], 'Visible', 'on');
        
        % DC별 이벤트 개수 집계
        dcEventCounts = struct();
        totalEventsAll = 0;
        for dcIdx = 1:length(allDCProfiles)
            profName = allDCProfiles{dcIdx};
            count = 0;
            if isfield(allCycleData, profName)
                cycles = fieldnames(allCycleData.(profName));
                for cIdx = 1:length(cycles)
                    evtDataList = allCycleData.(profName).(cycles{cIdx});
                    if isstruct(evtDataList)
                        evtDataList = {evtDataList};
                    end
                    count = count + length(evtDataList);
                end
            end
            dcEventCounts.(profName) = count;
            totalEventsAll = totalEventsAll + count;
        end
        
        % 색상은 MATLAB 기본 color order에 맡김 (이벤트마다 자동으로 변경)
        
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
                    
                    evtDataList = allCycleData.(profName).(cycleFieldName);
                    if isstruct(evtDataList)
                        evtDataList = {evtDataList};
                    end
                    
                    % 표시용 사이클 이름 (cyc_ 접두사 제거)
                    if startsWith(cycleFieldName, 'cyc_')
                        cycleDisplayName = cycleFieldName(5:end);
                    else
                        cycleDisplayName = cycleFieldName;
                    end
                    
                    legendAdded = false;
                    for e = 1:length(evtDataList)
                        evtData = evtDataList{e};
                        if ~isfield(evtData, 'I') || ~isfield(evtData, 't') || isempty(evtData.t)
                            continue;
                        end
                        t_rel = evtData.t - evtData.t(1);  % 상대 시간
                        
                        if ~legendAdded
                            h = plot(t_rel, evtData.I, '-', 'LineWidth', 1.5, ...
                                     'DisplayName', cycleDisplayName);
                            legendHandles(end+1) = h;
                            legendEntries{end+1} = cycleDisplayName;
                            legendAdded = true;
                        else
                            plot(t_rel, evtData.I, '-', 'LineWidth', 1.0, ...
                                 'HandleVisibility', 'off');
                        end
                    end
                end
                end
            end
            
            xlabel('Time (s)', 'FontSize', 10, 'FontWeight', 'bold');
            ylabel('Current (A)', 'FontSize', 10, 'FontWeight', 'bold');
            
            % 사이클별 이벤트 개수 (0,200,400,600,800)
            cyclesWanted = {'0cyc','200cyc','400cyc','600cyc','800cyc'};
            countsPerCycle = zeros(1, length(cyclesWanted));
            if isfield(allCycleData, profName)
                for cyIdx = 1:length(cyclesWanted)
                    cycName = cyclesWanted{cyIdx};
                    cycField = ['cyc_' cycName];
                    evtList = {};
                    if isfield(allCycleData.(profName), cycField)
                        evtList = allCycleData.(profName).(cycField);
                    elseif isfield(allCycleData.(profName), cycName)
                        evtList = allCycleData.(profName).(cycName);
                    end
                    if isstruct(evtList)
                        evtList = {evtList};
                    end
                    countsPerCycle(cyIdx) = numel(evtList);
                end
            end
            title(sprintf('%s (n=%d,%d,%d,%d,%d)', profName, countsPerCycle), ...
                  'FontSize', 11, 'FontWeight', 'bold');
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
end

fprintf('\n=== Current Profile Visualization Complete ===\n');

% %% 탈락 원인 분석 시각화
% fprintf('\n=== Creating Failure Analysis Visualization ===\n');
% failureDir = fullfile(outputDir, 'figures', 'FailureAnalysis');
% if ~exist(failureDir, 'dir')
%     mkdir(failureDir);
% end
% 
% if ~isempty(ref_channel_match)
%     for socIdx = 1:length(targetSOCs)
%         targetSOC = targetSOCs{socIdx};
%         if isstring(targetSOC)
%             targetSOC = char(targetSOC);
%         end
%         
%         for etIdx = 1:length(eventTypes)
%             eventType = eventTypes{etIdx};
%             
%             passCycleData = struct();
%             rejectCycleData = struct();
%             
%             for i = 1:length(matFiles)
%                 fileName = matFiles(i).name;
%                 token = regexp(fileName, 'parsedDriveCycle_(\d+cyc)_filtered', 'tokens');
%                 if isempty(token), continue; end
%                 cycleType = token{1}{1};
%                 
%                 savePath = fullfile(outputDir, sprintf('Lab_DC_Events_Raw_%s.mat', cycleType));
%                 if ~exist(savePath, 'file')
%                     continue;
%                 end
%                 
%                 tempData = load(savePath, 'rawEvents', 'rejectedEvents');
%                 if ~isfield(tempData, 'rawEvents') || ~isfield(tempData, 'rejectedEvents')
%                     continue;
%                 end
%                 rawEvents_loaded = tempData.rawEvents;
%                 rejectedEvents_loaded = tempData.rejectedEvents;
%                 
%                 % pass/reject 구조체 이름 찾기
%                 allStructNames = fieldnames(rawEvents_loaded);
%                 structNamePass = '';
%                 for s = 1:length(allStructNames)
%                     structName_candidate = allStructNames{s};
%                     if startsWith(structName_candidate, reference_channel) && ...
%                        endsWith(structName_candidate, ['_' eventType])
%                         structNamePass = structName_candidate;
%                         break;
%                     end
%                 end
%                 allRejectNames = fieldnames(rejectedEvents_loaded);
%                 structNameReject = '';
%                 for s = 1:length(allRejectNames)
%                     structName_candidate = allRejectNames{s};
%                     if startsWith(structName_candidate, reference_channel) && ...
%                        endsWith(structName_candidate, ['_' eventType '_rejected'])
%                         structNameReject = structName_candidate;
%                         break;
%                     end
%                 end
%                 
%                 if isempty(structNamePass)
%                     continue;
%                 end
%                 
%                 availableSOCs = fieldnames(rawEvents_loaded.(structNamePass));
%                 if ~ismember(targetSOC, availableSOCs)
%                     continue;
%                 end
%                 
%                 % cycleType 유효 필드명
%                 if ~isempty(regexp(cycleType, '^\d', 'once'))
%                     validCycleFieldName = ['cyc_' cycleType];
%                 else
%                     validCycleFieldName = cycleType;
%                 end
%                 
%                 profs = fieldnames(rawEvents_loaded.(structNamePass).(targetSOC));
%                 for p = 1:length(profs)
%                     profName = profs{p};
%                     
%                     % pass 이벤트 수집
%                     eventsPass = fieldnames(rawEvents_loaded.(structNamePass).(targetSOC).(profName));
%                     if ~isfield(passCycleData, profName)
%                         passCycleData.(profName) = struct();
%                     end
%                     if ~isfield(passCycleData.(profName), validCycleFieldName)
%                         passCycleData.(profName).(validCycleFieldName) = {};
%                     end
%                     for e = 1:length(eventsPass)
%                         evtData = rawEvents_loaded.(structNamePass).(targetSOC).(profName).(eventsPass{e});
%                         if isfield(evtData, 'I') && isfield(evtData, 't')
%                             passCycleData.(profName).(validCycleFieldName){end+1} = evtData;
%                         end
%                     end
%                     
%                     % reject 이벤트 수집
%                     if ~isempty(structNameReject) && isfield(rejectedEvents_loaded.(structNameReject), targetSOC) && ...
%                        isfield(rejectedEvents_loaded.(structNameReject).(targetSOC), profName)
%                         eventsReject = fieldnames(rejectedEvents_loaded.(structNameReject).(targetSOC).(profName));
%                         if ~isfield(rejectCycleData, profName)
%                             rejectCycleData.(profName) = struct();
%                         end
%                         if ~isfield(rejectCycleData.(profName), validCycleFieldName)
%                             rejectCycleData.(profName).(validCycleFieldName) = {};
%                         end
%                         for e = 1:length(eventsReject)
%                             evtData = rejectedEvents_loaded.(structNameReject).(targetSOC).(profName).(eventsReject{e});
%                             if isfield(evtData, 'I') && isfield(evtData, 't')
%                                 rejectCycleData.(profName).(validCycleFieldName){end+1} = evtData;
%                             end
%                         end
%                     end
%                 end
%             end
%             
%             % 시각화
%             allDCProfiles = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};
%             fig = figure('Name', sprintf('Failure Analysis - %s (%s)', eventType, targetSOC), ...
%                          'Position', [100, 100, 1600, 900], 'Visible', 'off');
%             nCols = 4;
%             nRows = 2;
%             
%             for dcIdx = 1:length(allDCProfiles)
%                 profName = allDCProfiles{dcIdx};
%                 subplot(nRows, nCols, dcIdx);
%                 hold on;
%                 
%                 passLegendAdded = false;
%                 rejectLegendAdded = false;
%                 stdLegendAdded = false;
%                 
%                 if isfield(passCycleData, profName)
%                     cycles = fieldnames(passCycleData.(profName));
%                     for cycIdx = 1:length(cycles)
%                         cycleFieldName = cycles{cycIdx};
%                         evtList = passCycleData.(profName).(cycleFieldName);
%                         for e = 1:length(evtList)
%                             evtData = evtList{e};
%                             t_rel = evtData.t - evtData.t(1);
%                             if ~passLegendAdded
%                                 plot(t_rel, evtData.I, 'b-', 'LineWidth', 1.0, 'DisplayName', 'pass');
%                                 passLegendAdded = true;
%                             else
%                                 plot(t_rel, evtData.I, 'b-', 'LineWidth', 1.0, 'HandleVisibility', 'off');
%                             end
%                         end
%                     end
%                 end
%                 
%                 if isfield(rejectCycleData, profName)
%                     cycles = fieldnames(rejectCycleData.(profName));
%                     for cycIdx = 1:length(cycles)
%                         cycleFieldName = cycles{cycIdx};
%                         evtList = rejectCycleData.(profName).(cycleFieldName);
%                         for e = 1:length(evtList)
%                             evtData = evtList{e};
%                             t_rel = evtData.t - evtData.t(1);
%                             if ~rejectLegendAdded
%                                 plot(t_rel, evtData.I, 'r-', 'LineWidth', 1.0, 'DisplayName', 'reject');
%                                 rejectLegendAdded = true;
%                             else
%                                 plot(t_rel, evtData.I, 'r-', 'LineWidth', 1.0, 'HandleVisibility', 'off');
%                             end
%                             
%                             if isfield(evtData, 'reason') && strcmp(evtData.reason, 'std')
%                                 if isfield(evtData, 'I_raw') && ~isempty(evtData.I_raw)
%                                     plot(t_rel, evtData.I_raw, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
%                                 end
%                                 if isfield(evtData, 'std_calc_indices') && all(~isnan(evtData.std_calc_indices))
%                                     idx1 = evtData.std_calc_indices(1);
%                                     idx2 = evtData.std_calc_indices(2);
%                                     if idx1 >= 1 && idx2 <= length(evtData.t)
%                                         x1 = evtData.t(idx1) - evtData.t(1);
%                                         x2 = evtData.t(idx2) - evtData.t(1);
%                                         yL = ylim;
%                                         patch([x1 x2 x2 x1], [yL(1) yL(1) yL(2) yL(2)], ...
%                                               [1 0 0], 'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
%                                               'HandleVisibility', 'off');
%                                         if ~stdLegendAdded
%                                             plot([x1 x2], [yL(2) yL(2)], 'r-', 'LineWidth', 3, 'DisplayName', 'std window');
%                                             stdLegendAdded = true;
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%                 
%                 title(profName);
%                 xlabel('Time (s)');
%                 ylabel('Current (A)');
%                 grid on;
%                 if passLegendAdded || rejectLegendAdded || stdLegendAdded
%                     legend('Location', 'best', 'FontSize', 8);
%                 end
%                 hold off;
%             end
%             
%             sgtitle(sprintf('Failure Analysis - %s (%s, Channel %s)', eventType, targetSOC, ref_ch_num), ...
%                     'FontSize', 14, 'FontWeight', 'bold');
%             saveFileName = sprintf('FailureAnalysis_%s_%s_%s', eventType, targetSOC, ref_ch_num);
%             saveFileName = regexprep(saveFileName, '[^a-zA-Z0-9_]', '_');
%             savePath_fig = fullfile(failureDir, [saveFileName, '.fig']);
%             saveas(fig, savePath_fig);
%             fprintf('  Saved failure analysis: %s\n', savePath_fig);
%             close(fig);
%         end
%     end
% end

fprintf('\n분석 가이드: 특정 DC에서 탈락자가 많은 이유는 무엇인가? (예: 특정 SOC에서의 전류 노이즈, 급격한 부하 변동 등)\n');

end % paramSets loop

%% 로컬 함수들

function [events, filter_stats, rejected_events] = detect_events_in_channel(I, V, t_seconds, current_threshold, min_duration, max_I_std, max_I_range, dt, debug_tag, event_type_selection)
    % 단일 채널에서 이벤트를 검출하고 반환
    % 반환값: 
    %   events: events 구조체 배열 (각 요소는 start_time, end_time, start_idx, end_idx, type, I_seg, V_seg, t_seg, duration 포함)
    %   filter_stats: 필터링 통계 구조체
    
    if nargin < 9
        debug_tag = '';
    end
    if nargin < 10
        event_type_selection = 'Both';
    end
    
    events = struct('start_time', {}, 'end_time', {}, 'start_idx', {}, 'end_idx', {}, ...
                    'type', {}, 'I_seg', {}, 'V_seg', {}, 't_seg', {}, 'duration', {}, ...
                    'reason', {}, 'std_calc_indices', {});
    rejected_events = struct('start_time', {}, 'end_time', {}, 'start_idx', {}, 'end_idx', {}, ...
                             'type', {}, 'I_seg', {}, 'V_seg', {}, 't_seg', {}, 'duration', {}, ...
                             'reason', {}, 'std_calc_indices', {}, 'I_raw_seg', {}, 'I_filt_seg', {});
    
    % 필터링 통계 초기화 (타입별로 분리)
    filter_stats = struct();
    filter_stats.total_candidates_charge = 0;
    filter_stats.total_candidates_discharge = 0;
    filter_stats.filtered_by_duration_charge = 0;
    filter_stats.filtered_by_duration_discharge = 0;
    filter_stats.filtered_by_std_charge = 0;
    filter_stats.filtered_by_std_discharge = 0;
    filter_stats.filtered_by_direction_charge = 0;
    filter_stats.filtered_by_range_charge = 0;
    filter_stats.filtered_by_range_discharge = 0;
    filter_stats.final_count_charge = 0;
    filter_stats.final_count_discharge = 0;
    % 필터링 실패 사유 카운트
    filter_stats.fail_duration = 0;
    filter_stats.fail_std = 0;
    filter_stats.fail_direction = 0;
    filter_stats.fail_range = 0;
    
    % 표준편차 계산 구간 설정
    std_calc_mode = 'trimmed'; % 'trimmed' | 'fixed_start' | 'full'
    trim_sec = 3;
    fixed_duration = 30;
    debug_std_stats = true; % 표준편차 구간 통계 출력
    
    if length(t_seconds) < 10 || length(I) < 10 || length(V) < 10
        return;
    end
    
    % 1. 이벤트 감지용 데이터 (1초 간격으로 Decimation)
    decim_factor = 1;
    if abs(dt - 0.1) < 1e-6
        decim_factor = 10;
    end
    I_detect_raw = I(1:decim_factor:end);
    V_detect_raw = V(1:decim_factor:end);
    t_detect = t_seconds(1:decim_factor:end);
    dt_detect = dt * decim_factor;
    
    % 2. 이동평균 (3초 윈도우) - 1초 간격 데이터에만 적용
    target_window_sec = 3;
    window_size = max(1, round(target_window_sec / dt_detect));
    if length(I_detect_raw) >= window_size
        I_detect = movmean(I_detect_raw, window_size);
        V_detect = movmean(V_detect_raw, window_size);
    else
        I_detect = movmean(I_detect_raw, length(I_detect_raw));
        V_detect = movmean(V_detect_raw, length(V_detect_raw));
    end
    
    % 세그먼트 추출용 원본 데이터 (필터 미적용)
    I_filt = I;
    V_filt = V;
    
    % 3. 이벤트 감지 (Idle -> Active) - Decimation된 데이터 사용
    is_idle = abs(I_detect) < current_threshold;
    is_driving = abs(I_detect) >= current_threshold;
    idle_to_driving = find(is_idle(1:end-1) & is_driving(2:end));
    
    if isempty(idle_to_driving)
        filter_stats.total_candidates_charge = 0;
        filter_stats.total_candidates_discharge = 0;
        return;
    end
    
    evt_count = 0;
    evt_count_charge = 0;
    evt_count_discharge = 0;
    rej_count = 0;
    for k = 1:length(idle_to_driving)
        idx1 = idle_to_driving(k);
        start_driving_idx = idx1 + 1;
        
        % 이벤트 타입 판별 - 필터링된 데이터 I_filt 사용
        event_type = sign(I_detect(start_driving_idx));
        % 주석: event_type == 0 체크는 불필요 (abs(I_filt) >= threshold이므로 0일 수 없음)
        if event_type > 0
            type_label = 'Charge';
        else
            type_label = 'Discharge';
        end
        
        % 타입별 후보 카운트
        if event_type > 0
            filter_stats.total_candidates_charge = filter_stats.total_candidates_charge + 1;
        else
            filter_stats.total_candidates_discharge = filter_stats.total_candidates_discharge + 1;
        end
        
        % driving 구간의 끝 찾기
        driving_end_idx = start_driving_idx;
        while driving_end_idx <= length(I_detect)
            if abs(I_detect(driving_end_idx)) >= current_threshold && sign(I_detect(driving_end_idx)) == event_type
                driving_end_idx = driving_end_idx + 1;
            else
                break;
            end
        end
        driving_end_idx = driving_end_idx - 1;
        
        % Decimation 인덱스를 원본 인덱스로 변환
        start_idx = (idx1 - 1) * decim_factor + 1;
        end_idx = (driving_end_idx - 1) * decim_factor + 1;
        start_driving_idx_full = (start_driving_idx - 1) * decim_factor + 1;
        driving_end_idx_full = end_idx;
        
        % 경계 보정
        start_idx = max(1, min(start_idx, length(I_filt)));
        end_idx = max(1, min(end_idx, length(I_filt)));
        start_driving_idx_full = max(1, min(start_driving_idx_full, length(I_filt)));
        driving_end_idx_full = max(1, min(driving_end_idx_full, length(I_filt)));
        
        % 세그먼트 추출 (필터링된 데이터 사용)
        I_seg = I_filt(start_idx:end_idx);
        V_seg = V_filt(start_idx:end_idx);
        t_seg = t_seconds(start_idx:end_idx);
        I_raw_seg = I(start_idx:end_idx);
        
        % 지속시간 체크
        driving_time = t_detect(driving_end_idx) - t_detect(start_driving_idx);
        if driving_time < min_duration
            if event_type > 0
                filter_stats.filtered_by_duration_charge = filter_stats.filtered_by_duration_charge + 1;
            else
                filter_stats.filtered_by_duration_discharge = filter_stats.filtered_by_duration_discharge + 1;
            end
            filter_stats.fail_duration = filter_stats.fail_duration + 1;
            rej_count = rej_count + 1;
            rejected_events(rej_count).start_time = t_seconds(start_driving_idx_full);
            rejected_events(rej_count).end_time = t_seconds(driving_end_idx_full);
            rejected_events(rej_count).start_idx = start_idx;
            rejected_events(rej_count).end_idx = end_idx;
            rejected_events(rej_count).type = type_label;
            rejected_events(rej_count).I_seg = I_seg;
            rejected_events(rej_count).V_seg = V_seg;
            rejected_events(rej_count).t_seg = t_seg;
            rejected_events(rej_count).duration = driving_time;
            rejected_events(rej_count).reason = 'duration';
            rejected_events(rej_count).std_calc_indices = [NaN, NaN];
            rejected_events(rej_count).I_raw_seg = I_raw_seg;
            rejected_events(rej_count).I_filt_seg = I_seg;
            continue;
        end
        
        % 전류 안정성 체크 (표준편차 계산 구간 선택)
        switch std_calc_mode
            case 'trimmed'
                t_start = t_seg(1) + trim_sec;
                t_end = t_seg(end) - trim_sec;
            case 'fixed_start'
                t_start = t_seg(1) + trim_sec;
                t_end = t_seg(1) + fixed_duration;
            otherwise % 'full'
                t_start = t_seg(1);
                t_end = t_seg(end);
        end
        
        if t_end > t_seg(end)
            t_end = t_seg(end);
        end
        
        calc_idx = find(t_seg >= t_start & t_seg <= t_end);
        if numel(calc_idx) < 5
            I_std_val = inf;
            calc_idx_start = NaN;
            calc_idx_end = NaN;
        else
            calc_idx_start = calc_idx(1);
            calc_idx_end = calc_idx(end);
            I_std_val = std(I_seg(calc_idx_start:calc_idx_end));
            % 디버그 출력은 통과 이벤트만 아래에서 수행
        end
        
        if I_std_val > max_I_std
            if event_type > 0
                filter_stats.filtered_by_std_charge = filter_stats.filtered_by_std_charge + 1;
            else
                filter_stats.filtered_by_std_discharge = filter_stats.filtered_by_std_discharge + 1;
            end
            filter_stats.fail_std = filter_stats.fail_std + 1;
            rej_count = rej_count + 1;
            rejected_events(rej_count).start_time = t_seconds(start_driving_idx_full);
            rejected_events(rej_count).end_time = t_seconds(driving_end_idx_full);
            rejected_events(rej_count).start_idx = start_idx;
            rejected_events(rej_count).end_idx = end_idx;
            rejected_events(rej_count).type = type_label;
            rejected_events(rej_count).I_seg = I_seg;
            rejected_events(rej_count).V_seg = V_seg;
            rejected_events(rej_count).t_seg = t_seg;
            rejected_events(rej_count).duration = driving_time;
            rejected_events(rej_count).reason = 'std';
            rejected_events(rej_count).std_calc_indices = [calc_idx_start, calc_idx_end];
            rejected_events(rej_count).I_raw_seg = I_raw_seg;
            rejected_events(rej_count).I_filt_seg = I_seg;
            continue;
        end
        
        % 전류 범위(최대-최소) 체크 - CC처럼 일정한 전류만 통과
        I_range_val = max(I_seg(calc_idx_start:calc_idx_end)) - min(I_seg(calc_idx_start:calc_idx_end));
        if max_I_range > 0 && I_range_val > max_I_range
            if event_type > 0
                filter_stats.filtered_by_range_charge = filter_stats.filtered_by_range_charge + 1;
            else
                filter_stats.filtered_by_range_discharge = filter_stats.filtered_by_range_discharge + 1;
            end
            filter_stats.fail_range = filter_stats.fail_range + 1;
            rej_count = rej_count + 1;
            rejected_events(rej_count).start_time = t_seconds(start_driving_idx);
            rejected_events(rej_count).end_time = t_seconds(driving_end_idx);
            rejected_events(rej_count).start_idx = start_idx;
            rejected_events(rej_count).end_idx = end_idx;
            rejected_events(rej_count).type = type_label;
            rejected_events(rej_count).I_seg = I_seg;
            rejected_events(rej_count).V_seg = V_seg;
            rejected_events(rej_count).t_seg = t_seg;
            rejected_events(rej_count).duration = driving_time;
            rejected_events(rej_count).reason = 'range';
            rejected_events(rej_count).std_calc_indices = [calc_idx_start, calc_idx_end];
            rejected_events(rej_count).I_raw_seg = I_raw_seg;
            rejected_events(rej_count).I_filt_seg = I_seg;
            continue;
        end
        
        % 방향 체크 (충전 이벤트는 대부분 양수여야 함)
        if event_type > 0
            positive_ratio = sum(I_seg > 0) / length(I_seg);
            if positive_ratio < 0.5
                filter_stats.filtered_by_direction_charge = filter_stats.filtered_by_direction_charge + 1;
                filter_stats.fail_direction = filter_stats.fail_direction + 1;
                rej_count = rej_count + 1;
                rejected_events(rej_count).start_time = t_seconds(start_driving_idx);
                rejected_events(rej_count).end_time = t_seconds(driving_end_idx);
                rejected_events(rej_count).start_idx = start_idx;
                rejected_events(rej_count).end_idx = end_idx;
                rejected_events(rej_count).type = 'Charge';
                rejected_events(rej_count).I_seg = I_seg;
                rejected_events(rej_count).V_seg = V_seg;
                rejected_events(rej_count).t_seg = t_seg;
                rejected_events(rej_count).duration = driving_time;
                rejected_events(rej_count).reason = 'direction';
                rejected_events(rej_count).std_calc_indices = [calc_idx_start, calc_idx_end];
                rejected_events(rej_count).I_raw_seg = I_raw_seg;
                rejected_events(rej_count).I_filt_seg = I_seg;
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
        events(evt_count).start_time = t_seconds(start_driving_idx_full);
        events(evt_count).end_time = t_seconds(driving_end_idx_full);
        events(evt_count).start_idx = start_idx;
        events(evt_count).end_idx = end_idx;
        events(evt_count).type = type_str;
        events(evt_count).I_seg = I_seg;
        events(evt_count).V_seg = V_seg;
        events(evt_count).t_seg = t_seg;
        events(evt_count).duration = driving_time;
        events(evt_count).reason = 'pass';
        events(evt_count).std_calc_indices = [calc_idx_start, calc_idx_end];
        
        % % 통과 이벤트 std 디버그 출력 (식별자 포함)
        % if debug_std_stats && ~isnan(calc_idx_start) && ~isnan(calc_idx_end)
        %     if strcmp(event_type_selection, 'Charge') && ~strcmp(type_str, 'Charge')
        %         continue;
        %     elseif strcmp(event_type_selection, 'Discharge') && ~strcmp(type_str, 'Discharge')
        %         continue;
        %     end
        %     I_slice = I_seg(calc_idx_start:calc_idx_end);
        %     fprintf('STD-DEBUG-PASS [%s] evt=%d %s calc_idx=%d:%d, t=%.2f..%.2f, mean=%.3f, std=%.3f, min=%.3f, max=%.3f\n', ...
        %             type_str, evt_count, debug_tag, calc_idx_start, calc_idx_end, ...
        %             t_seg(calc_idx_start), t_seg(calc_idx_end), ...
        %             mean(I_slice), std(I_slice), min(I_slice), max(I_slice));
        % end
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

function parsave(savePath, rawEvents, rejectedEvents)
    save(savePath, 'rawEvents', 'rejectedEvents');
end

function [pass, reason, duration, stdVal, rangeVal, dirRatio] = evaluate_event_window( ...
    I_seg, t_seg, event_type_selection, min_duration, max_I_std, max_I_range, std_calc_mode, trim_sec, fixed_duration)
    pass = 0;
    reason = 'pass';
    duration = t_seg(end) - t_seg(1);
    stdVal = inf;
    rangeVal = max(I_seg) - min(I_seg);
    dirRatio = 0;
    
    if duration < min_duration
        reason = 'duration';
        return;
    end
    
    % 표준편차 계산 구간
    switch std_calc_mode
        case 'trimmed'
            t0 = t_seg(1) + trim_sec;
            t1 = t_seg(end) - trim_sec;
            calc_idx = find(t_seg >= t0 & t_seg <= t1);
        case 'fixed_start'
            t0 = t_seg(1) + trim_sec;
            t1 = t0 + fixed_duration;
            calc_idx = find(t_seg >= t0 & t_seg <= t1);
        otherwise
            calc_idx = 1:length(t_seg);
    end
    
    if numel(calc_idx) < 5
        reason = 'std';
        return;
    end
    
    stdVal = std(I_seg(calc_idx));
    if stdVal > max_I_std
        reason = 'std';
        return;
    end
    
    if max_I_range > 0 && rangeVal > max_I_range
        reason = 'range';
        return;
    end
    
    if strcmp(event_type_selection, 'Charge')
        dirRatio = sum(I_seg > 0) / length(I_seg);
        if dirRatio < 0.5
            reason = 'direction';
            return;
        end
    end
    
    pass = 1;
end
