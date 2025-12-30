%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 01_FieldData_EventAnalysis.m
% 이벤트 검출 및 Rchg 계산 (필드 데이터용)
% 
% 목적: 
% - all_events_raw_cell_level.mat에서 이벤트 로드
% - 필터링 조건 적용 (min_duration, delta_I, current_variation 등)
% - 각 이벤트에서 Rchg 계산 (1s, 5s, 10s, 30s)
% - 그룹화 (V_I_SOC)
%
% 출력:
% - FieldData_Rack01_charge_Events.mat (또는 discharge)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Configuration - User Settings
% =========================================================================
% 스크립트 위치를 기준으로 상대 경로 설정
scriptDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(scriptDir, 'EventsResults');
outputDir = fullfile(resultsDir, 'FieldData_Analysis');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

rawEventsFile = 'all_events_raw_cell_level.mat';

targetRack = 'Rack01';
targetEventType = 'charge';  % 'charge' or 'discharge'

% 필터링 파라미터
min_event_duration_sec = 60;  % 최소 이벤트 길이 (초)
min_delta_I = 64 * 0.10;     % 최소 전류 변화 (A)
current_variation_threshold = 3.0;  % 전류 변동률 임계값 (%)
min_events_per_year = 1;  % 모든 연도에서 최소 이벤트 수 (시각화/상관분석을 위해 모든 연도에 데이터 필요)

% 그룹화 파라미터 (AnalyzeResistanceFeatures.m과 동일)
voltage_bin_size = 0.01;  % 전압 bin 크기 (V) - AnalyzeResistanceFeatures.m과 동일
current_bin_size = 2;     % 전류 bin 크기 (A) - AnalyzeResistanceFeatures.m과 동일
% soc_bin_size = 1;        % SOC bin 크기 (%) - 사용 안 함 (SOC는 그룹화에서 제외)
% =========================================================================

fprintf('=== Field Data Event Analysis ===\n');
fprintf('Target Rack: %s\n', targetRack);
fprintf('Target Event Type: %s\n', targetEventType);
fprintf('Data directory: %s\n', resultsDir);

%% Load Raw Events Data
fprintf('\n=== Loading Raw Events Data ===\n');
rawEventsPath = fullfile(resultsDir, rawEventsFile);
if ~exist(rawEventsPath, 'file')
    error('Raw events file not found: %s', rawEventsPath);
end

load(rawEventsPath, 'all_events');
fprintf('Data loaded successfully.\n');

if ~isfield(all_events, targetRack)
    error('Target Rack "%s" not found in data.', targetRack);
end

rackData = all_events.(targetRack);
available_year_keys = fieldnames(rackData);
available_year_keys = available_year_keys(startsWith(available_year_keys, 'Y'));
fprintf('Available years: %s\n', strjoin(available_year_keys, ', '));

%% Collect and Filter Events
fprintf('\n=== Collecting and Filtering Events ===\n');

% 모든 이벤트 수집
all_events_list = {};
for y_idx = 1:length(available_year_keys)
    year_key = available_year_keys{y_idx};
    if ~isfield(rackData, year_key)
        continue;
    end
    
    year_data = rackData.(year_key);
    month_keys = fieldnames(year_data);
    
    for m_idx = 1:length(month_keys)
        month_key = month_keys{m_idx};
        if ~isfield(year_data, month_key)
            continue;
        end
        
        month_events = year_data.(month_key);
        if iscell(month_events)
            for e_idx = 1:length(month_events)
                evt = month_events{e_idx};
                if isstruct(evt) && isfield(evt, 'type')
                    if strcmpi(evt.type, targetEventType)
                        all_events_list{end+1} = evt;
                    end
                end
            end
        end
    end
end

fprintf('Total events collected: %d\n', length(all_events_list));

%% Filter Events
fprintf('\n=== Filtering Events ===\n');

filter_stats = struct();
filter_stats.total = length(all_events_list);
filter_stats.no_voltage = 0;
filter_stats.no_current = 0;
filter_stats.length_too_short = 0;
filter_stats.delta_I_too_small = 0;
filter_stats.current_variation_too_high = 0;
filter_stats.wrong_current_sign = 0;  % 충전/방전 전류 부호 불일치
filter_stats.passed = 0;

filtered_events = {};
grouped_events_map = containers.Map('KeyType', 'char', 'ValueType', 'any');

for e_idx = 1:length(all_events_list)
    evt = all_events_list{e_idx};
    
    % 전압 데이터 확인
    if ~isfield(evt, 'voltage_seq_cell_V') || isempty(evt.voltage_seq_cell_V)
        filter_stats.no_voltage = filter_stats.no_voltage + 1;
        continue;
    end
    
    % 전류 데이터 확인
    if ~isfield(evt, 'current_seq_cell_A') || isempty(evt.current_seq_cell_A)
        filter_stats.no_current = filter_stats.no_current + 1;
        continue;
    end
    
    voltage_seq = evt.voltage_seq_cell_V;
    current_seq = evt.current_seq_cell_A;
    
    % delta_I 필터링: 시작점부터 인덱스 10 이내의 전류 변화가 6.4A 이상인지 확인 (AnalyzeResistanceFeatures.m과 동일)
    if length(current_seq) < 6  % 최소 21개 포인트 필요 (시작점 + 20)
        filter_stats.length_too_short = filter_stats.length_too_short + 1;
        continue;
    end
    
    % 시작점과 인덱스 10 사이의 전류 변화 확인: I(20) - I(1) >= min_delta_I (AnalyzeResistanceFeatures.m과 동일)
    idx_5 = min(5, length(current_seq));  % 인덱스 10
    delta_I = current_seq(idx_5) - current_seq(1);
    if delta_I < min_delta_I
        filter_stats.delta_I_too_small = filter_stats.delta_I_too_small + 1;
        continue;
    end
    
    % 평균 전류 (10초부터 종료인덱스-5까지) 및 전류 변동성 필터링
    % 시간 시퀀스에서 상대 시간 계산
    time_seq = evt.time_seq_datetime;
    if isdatetime(time_seq(1))
        time_start = time_seq(1);
        t_rel = seconds(time_seq - time_start);
    elseif isduration(time_seq(1))
        t_rel = seconds(time_seq - time_seq(1));
    else
        % Numeric인 경우 datenum 또는 초 단위로 가정
        if time_seq(1) > 730000
            % datenum 형식
            dt = datetime(time_seq, 'ConvertFrom', 'datenum');
            time_start = dt(1);
            t_rel = seconds(dt - time_start);
        else
            % 초 단위로 가정
            t_rel = time_seq - time_seq(1);
        end
    end
    
    % 이벤트 지속 시간 확인 (초 단위)
    if length(t_rel) < 2
        filter_stats.length_too_short = filter_stats.length_too_short + 1;
        continue;
    end
    event_duration = t_rel(end) - t_rel(1);
    if event_duration < min_event_duration_sec
        filter_stats.length_too_short = filter_stats.length_too_short + 1;
        continue;
    end
    
    % 10초 이상인 첫 번째 인덱스 찾기
    idx_10s = find(t_rel >= 10.0, 1);
    
    % 종료 인덱스-5 계산
    idx_end_minus_5 = length(current_seq) - 5;
    
    if isempty(idx_10s) || idx_end_minus_5 < idx_10s || idx_end_minus_5 < 1
        % 10초 이상인 지점이 없거나 종료인덱스-5가 10초 지점보다 앞이면 스킵
        filter_stats.current_variation_too_high = filter_stats.current_variation_too_high + 1;
        continue;
    end
    
    % 10초부터 종료인덱스-5까지의 전류 데이터
    stable_current = abs(current_seq(idx_10s:idx_end_minus_5));
        avg_current = mean(stable_current);
        
    % 전류 변동성 필터링: 변동률이 임계값 이내인 경우만 통과
        if avg_current > 0
            current_max = max(stable_current);
            current_min = min(stable_current);
            current_variation = (current_max - current_min) / avg_current * 100;
            if current_variation > current_variation_threshold
                filter_stats.current_variation_too_high = filter_stats.current_variation_too_high + 1;
            continue;
        end
    else
        filter_stats.current_variation_too_high = filter_stats.current_variation_too_high + 1;
        continue;
    end
    
    % 충전/방전 이벤트 필터링: 첫 인덱스 이후 전류 부호 확인
    % 충전: 첫 인덱스 이후 전류가 항상 양수
    % 방전: 첫 인덱스 이후 전류가 항상 음수
    if length(current_seq) > 1
        current_after_first = current_seq(2:end);
        
        if strcmpi(targetEventType, 'charge')
            % 충전 이벤트: 첫 인덱스 이후 전류가 모두 양수여야 함
            if any(current_after_first <= 0)
                filter_stats.wrong_current_sign = filter_stats.wrong_current_sign + 1;
                continue;
            end
        elseif strcmpi(targetEventType, 'discharge')
            % 방전 이벤트: 첫 인덱스 이후 전류가 모두 음수여야 함
            if any(current_after_first >= 0)
                filter_stats.wrong_current_sign = filter_stats.wrong_current_sign + 1;
                continue;
            end
        end
    end
    
    % 그룹화: 전압, 전류, SOC bin 계산
    voltage_bin = floor(voltage_seq(1) / voltage_bin_size) * voltage_bin_size;
    current_bin = floor(avg_current / current_bin_size) * current_bin_size;
    
    % SOC bin 계산 (시작 SOC 사용) - AnalyzeResistanceFeatures.m과 동일
    % 그룹 키 생성 (전압, 전류만 사용, SOC 제외)
    group_key = sprintf('V%.2f_I%d', voltage_bin, round(current_bin));
    
    % 이벤트에 평균 전류값 저장 (그룹 필터링에서 사용)
    evt.avg_current_10s_to_end5 = avg_current;
    
    % 이벤트 저장
    if ~isKey(grouped_events_map, group_key)
        grouped_events_map(group_key) = {evt};
    else
        current_events = grouped_events_map(group_key);
        current_events{end+1} = evt;
        grouped_events_map(group_key) = current_events;
    end
    
    filtered_events{end+1} = evt;
    filter_stats.passed = filter_stats.passed + 1;
end

% 필터링 통계 출력
fprintf('\n--- Filtering Statistics ---\n');
fprintf('Total events: %d\n', filter_stats.total);
fprintf('  - No voltage data: %d\n', filter_stats.no_voltage);
fprintf('  - No current data: %d\n', filter_stats.no_current);
fprintf('  - Length < 21 points: %d\n', filter_stats.length_too_short);
fprintf('  - Duration < %d sec: %d\n', min_event_duration_sec, filter_stats.length_too_short);
fprintf('  - delta_I <= %.2f A: %d\n', min_delta_I, filter_stats.delta_I_too_small);
fprintf('  - Current variation > %.1f%%: %d\n', current_variation_threshold, filter_stats.current_variation_too_high);
fprintf('  - Wrong current sign (charge/discharge): %d\n', filter_stats.wrong_current_sign);
fprintf('  - Passed all filters: %d\n', filter_stats.passed);
fprintf('Total groups created: %d\n', grouped_events_map.Count);
fprintf('Grouping parameters: voltage_bin=%.2fV, current_bin=%dA (SOC excluded)\n', ...
    voltage_bin_size, current_bin_size);
fprintf('Current averaging range: 10s to end-5s\n');

%% Debug: Check year distribution of all filtered events
fprintf('\n--- Debug: Year distribution of all filtered events ---\n');
all_filtered_year_counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
for yk = available_year_keys'
    all_filtered_year_counts(yk{1}) = 0;
end

for e_idx = 1:length(filtered_events)
    evt = filtered_events{e_idx};
    
    % 연도 추출
    year_key = 'Unknown';
    if isfield(evt, 'time_seq_datetime') && ~isempty(evt.time_seq_datetime)
        time_point = evt.time_seq_datetime(1);
        if isdatetime(time_point)
            event_year = year(time_point);
            year_key = sprintf('Y%d', event_year);
        elseif isduration(time_point)
            if isfield(evt, 'file_info') && isfield(evt.file_info, 'filename')
                filename = evt.file_info.filename;
                if length(filename) >= 4
                    try
                        event_year = str2double(filename(1:4));
                        if ~isnan(event_year)
                            year_key = sprintf('Y%d', event_year);
                        end
                    catch
                    end
                end
            end
        elseif isnumeric(time_point)
            try
                if time_point > 730000
                    dt = datetime(time_point, 'ConvertFrom', 'datenum');
                    event_year = year(dt);
                    year_key = sprintf('Y%d', event_year);
                end
            catch
            end
        end
    end
    if strcmp(year_key, 'Unknown') && isfield(evt, 'file_info') && isfield(evt.file_info, 'filename')
        filename = evt.file_info.filename;
        if length(filename) >= 4
            try
                event_year = str2double(filename(1:4));
                if ~isnan(event_year) && event_year >= 2000 && event_year <= 2100
                    year_key = sprintf('Y%d', event_year);
                end
            catch
            end
        end
    end
    
    if isKey(all_filtered_year_counts, year_key)
        all_filtered_year_counts(year_key) = all_filtered_year_counts(year_key) + 1;
    end
end

fprintf('Total filtered events: %d\n', length(filtered_events));
for yk = available_year_keys'
    fprintf('  %s: %d events\n', yk{1}, all_filtered_year_counts(yk{1}));
end
fprintf('\n');

%% Filter Groups: Remove events with current deviating > 2A from group mean
fprintf('\n=== Filtering Groups (Remove Events with Current Deviation > 2A) ===\n');
fprintf('Removing events where average current deviates > 2A from group mean...\n');

% 그룹별 평균 전류값 계산 및 필터링
group_keys_all = keys(grouped_events_map);
current_deviation_threshold = 2.0;  % 2A 임계값
filtered_grouped_events_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
removed_by_current_deviation = 0;

for g_idx = 1:length(group_keys_all)
    group_key = group_keys_all{g_idx};
    group_events = grouped_events_map(group_key);
    
    % 그룹 내 모든 이벤트의 평균 전류값 추출
    avg_currents = [];
    for e_idx = 1:length(group_events)
        evt = group_events{e_idx};
        if isfield(evt, 'avg_current_10s_to_end5')
            avg_currents(end+1) = evt.avg_current_10s_to_end5;
        end
    end
    
    if isempty(avg_currents)
        % 평균 전류값이 없는 경우 그룹 전체 제외
        continue;
    end
    
    % 그룹 평균 전류값 계산
    group_mean_current = mean(avg_currents);
    
    % 그룹 평균에서 2A 이상 벗어난 이벤트 제외
    filtered_group_events = {};
    for e_idx = 1:length(group_events)
        evt = group_events{e_idx};
        if isfield(evt, 'avg_current_10s_to_end5')
            current_deviation = abs(evt.avg_current_10s_to_end5 - group_mean_current);
            if current_deviation <= current_deviation_threshold
                filtered_group_events{end+1} = evt;
            else
                removed_by_current_deviation = removed_by_current_deviation + 1;
            end
        else
            % 평균 전류값이 없는 경우 제외
            removed_by_current_deviation = removed_by_current_deviation + 1;
        end
    end
    
    % 필터링된 이벤트가 있는 경우만 저장
    if length(filtered_group_events) > 0
        filtered_grouped_events_map(group_key) = filtered_group_events;
    end
end

fprintf('Removed %d events due to current deviation > 2A from group mean\n', removed_by_current_deviation);
fprintf('Total groups before current deviation filtering: %d\n', length(group_keys_all));
fprintf('Total groups after current deviation filtering: %d\n', filtered_grouped_events_map.Count);

% 필터링된 그룹으로 교체
grouped_events_map = filtered_grouped_events_map;

%% Filter Groups: All years must have events
fprintf('\n=== Filtering Groups (All Years Must Have Events) ===\n');
fprintf('Filtering groups with >= %d events for each available year...\n', min_events_per_year);
fprintf('Available years: %s\n', strjoin(available_year_keys, ', '));

qualified_groups_map = containers.Map('KeyType', 'char', 'ValueType', 'any');

group_keys_all = keys(grouped_events_map);
fprintf('Total groups before filtering: %d\n', length(group_keys_all));

% 디버깅: 첫 번째 그룹의 연도 추출 테스트
if length(group_keys_all) > 0
    test_group_key = group_keys_all{1};
    test_group_events = grouped_events_map(test_group_key);
    fprintf('\n--- Debug: Testing year extraction for first group "%s" (%d events) ---\n', ...
        test_group_key, length(test_group_events));
    for test_idx = 1:min(3, length(test_group_events))
        test_evt = test_group_events{test_idx};
        fprintf('  Event %d:\n', test_idx);
        if isfield(test_evt, 'time_seq_datetime') && ~isempty(test_evt.time_seq_datetime)
            time_point = test_evt.time_seq_datetime(1);
            fprintf('    time_seq_datetime(1) type: %s\n', class(time_point));
            if isdatetime(time_point)
                fprintf('    Extracted year: %d\n', year(time_point));
            elseif isduration(time_point)
                fprintf('    Duration format, checking file_info...\n');
            elseif isnumeric(time_point)
                fprintf('    Numeric format, value: %f\n', time_point);
            end
        else
            fprintf('    No time_seq_datetime field\n');
        end
        if isfield(test_evt, 'file_info') && isfield(test_evt.file_info, 'filename')
            fprintf('    file_info.filename: %s\n', test_evt.file_info.filename);
        else
            fprintf('    No file_info.filename field\n');
        end
    end
    fprintf('\n');
end

% 디버깅: 그룹별 이벤트 수 및 연도 분포 확인
fprintf('\n--- Debug: Group event counts and year distribution ---\n');
group_event_counts = zeros(length(group_keys_all), 1);
for debug_g_idx = 1:min(10, length(group_keys_all))
    debug_group_key = group_keys_all{debug_g_idx};
    debug_group_events = grouped_events_map(debug_group_key);
    group_event_counts(debug_g_idx) = length(debug_group_events);
    
    % 연도별 분포 확인
    debug_year_counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for yk = available_year_keys'
        debug_year_counts(yk{1}) = 0;
    end
    
    for debug_e_idx = 1:length(debug_group_events)
        debug_evt = debug_group_events{debug_e_idx};
        debug_year_key = 'Unknown';
        if isfield(debug_evt, 'time_seq_datetime') && ~isempty(debug_evt.time_seq_datetime)
            debug_time_point = debug_evt.time_seq_datetime(1);
            if isdatetime(debug_time_point)
                debug_event_year = year(debug_time_point);
                debug_year_key = sprintf('Y%d', debug_event_year);
            elseif isduration(debug_time_point)
                if isfield(debug_evt, 'file_info') && isfield(debug_evt.file_info, 'filename')
                    debug_filename = debug_evt.file_info.filename;
                    if length(debug_filename) >= 4
                        try
                            debug_event_year = str2double(debug_filename(1:4));
                            if ~isnan(debug_event_year)
                                debug_year_key = sprintf('Y%d', debug_event_year);
                            end
                        catch
                        end
                    end
                end
            elseif isnumeric(debug_time_point)
                try
                    if debug_time_point > 730000
                        debug_dt = datetime(debug_time_point, 'ConvertFrom', 'datenum');
                        debug_event_year = year(debug_dt);
                        debug_year_key = sprintf('Y%d', debug_event_year);
                    end
                catch
                end
            end
        end
        if strcmp(debug_year_key, 'Unknown') && isfield(debug_evt, 'file_info') && isfield(debug_evt.file_info, 'filename')
            debug_filename = debug_evt.file_info.filename;
            if length(debug_filename) >= 4
                try
                    debug_event_year = str2double(debug_filename(1:4));
                    if ~isnan(debug_event_year) && debug_event_year >= 2000 && debug_event_year <= 2100
                        debug_year_key = sprintf('Y%d', debug_event_year);
                    end
                catch
                end
            end
        end
        if isKey(debug_year_counts, debug_year_key)
            debug_year_counts(debug_year_key) = debug_year_counts(debug_year_key) + 1;
        end
    end
    
    year_dist_str = {};
    for yk = available_year_keys'
        year_dist_str{end+1} = sprintf('%s:%d', yk{1}, debug_year_counts(yk{1}));
    end
    fprintf('  Group %d: "%s" - Total: %d events, Years: %s\n', ...
        debug_g_idx, debug_group_key, length(debug_group_events), strjoin(year_dist_str, ', '));
end
fprintf('\n');

for g_idx = 1:length(group_keys_all)
    group_key = group_keys_all{g_idx};
    group_events = grouped_events_map(group_key);
    
    % 연도별 이벤트 수 계산
    year_counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for yk = available_year_keys'
        year_counts(yk{1}) = 0;
    end
    
    for e_idx = 1:length(group_events)
        evt = group_events{e_idx};
        
        % 연도 추출 (기존 AnalyzeResistanceFeatures.m 로직 참고)
        year_key = 'Unknown';
        if isfield(evt, 'time_seq_datetime') && ~isempty(evt.time_seq_datetime)
            time_point = evt.time_seq_datetime(1);
            if isdatetime(time_point)
                event_year = year(time_point);
                year_key = sprintf('Y%d', event_year);
            elseif isduration(time_point)
                % duration 형식인 경우 file_info에서 연도 추출
                if isfield(evt, 'file_info') && isfield(evt.file_info, 'filename')
                    filename = evt.file_info.filename;
                    if length(filename) >= 4
                        try
                            event_year = str2double(filename(1:4));
                            if ~isnan(event_year)
                                year_key = sprintf('Y%d', event_year);
                            end
                        catch
                        end
                    end
                end
            elseif isnumeric(time_point)
                % 숫자 형식인 경우 (datenum 등)
                try
                    if time_point > 730000  % datenum 형식 (2000년 이후)
                        dt = datetime(time_point, 'ConvertFrom', 'datenum');
                        event_year = year(dt);
                        year_key = sprintf('Y%d', event_year);
                    end
                catch
                end
            end
        end
        
        % file_info에서 연도 추출 시도 (time_seq_datetime이 없는 경우)
        if strcmp(year_key, 'Unknown') && isfield(evt, 'file_info') && isfield(evt.file_info, 'filename')
            filename = evt.file_info.filename;
            if length(filename) >= 4
                try
                    event_year = str2double(filename(1:4));
                    if ~isnan(event_year) && event_year >= 2000 && event_year <= 2100
                        year_key = sprintf('Y%d', event_year);
                    end
                catch
                end
            end
        end
        
        if isKey(year_counts, year_key)
            year_counts(year_key) = year_counts(year_key) + 1;
        end
    end
    
    % 모든 연도에 최소 이벤트 수가 있는지 확인
    is_qualified = true;
    missing_years = {};
    for yk = available_year_keys'
        if year_counts(yk{1}) < min_events_per_year
            is_qualified = false;
            missing_years{end+1} = sprintf('%s (count: %d)', yk{1}, year_counts(yk{1}));
        end
    end
    
    if is_qualified
        qualified_groups_map(group_key) = group_events;
    else
        fprintf('  - Group "%s" disqualified. Missing/insufficient data in: %s\n', ...
            group_key, strjoin(missing_years, ', '));
    end
end

fprintf('Total groups after filtering: %d\n', qualified_groups_map.Count);

if qualified_groups_map.Count == 0
    fprintf('\n=== ERROR: No groups met the criteria ===\n');
    fprintf('Available years: %s\n', strjoin(available_year_keys, ', '));
    fprintf('Required: >= %d events per year\n', min_events_per_year);
    fprintf('\nSuggestions:\n');
    fprintf('  1. Lower min_events_per_year (currently: %d)\n', min_events_per_year);
    fprintf('  2. Check if data exists for all years\n');
    fprintf('  3. Review filtering criteria (delta_I, current_variation, etc.)\n');
    error('No groups met the criteria. Try lowering min_events_per_year.');
end

%% Create Group-wise Current-Time Plots
fprintf('\n=== Creating Group-wise Current-Time Plots ===\n');
groupCheckDir = fullfile(outputDir, 'figures', 'events check', 'by_group');
if ~exist(groupCheckDir, 'dir')
    mkdir(groupCheckDir);
end

group_keys = keys(qualified_groups_map);
fprintf('Creating plots for %d qualified groups...\n', length(group_keys));

for g_idx = 1:length(group_keys)
    group_key = group_keys{g_idx};
    group_events = qualified_groups_map(group_key);
    
    fprintf('  Group %d/%d: %s (%d events)\n', g_idx, length(group_keys), group_key, length(group_events));
    
    % 그룹별로 연도별 이벤트 분류
    group_year_events_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for yk = available_year_keys'
        group_year_events_map(yk{1}) = {};
    end
    
    for e_idx = 1:length(group_events)
        evt = group_events{e_idx};
        
        % 연도 추출
        year_key = 'Unknown';
        if isfield(evt, 'time_seq_datetime') && ~isempty(evt.time_seq_datetime)
            time_point = evt.time_seq_datetime(1);
            if isdatetime(time_point)
                event_year = year(time_point);
                year_key = sprintf('Y%d', event_year);
            elseif isduration(time_point)
                if isfield(evt, 'file_info') && isfield(evt.file_info, 'filename')
                    filename = evt.file_info.filename;
                    if length(filename) >= 4
                        try
                            event_year = str2double(filename(1:4));
                            if ~isnan(event_year) && event_year >= 2000 && event_year <= 2100
                                year_key = sprintf('Y%d', event_year);
                            end
                        catch
                        end
                    end
                end
            elseif isnumeric(time_point)
                try
                    if time_point > 730000
                        dt = datetime(time_point, 'ConvertFrom', 'datenum');
                        event_year = year(dt);
                        year_key = sprintf('Y%d', event_year);
                    end
                catch
                end
            end
        end
        if strcmp(year_key, 'Unknown') && isfield(evt, 'file_info') && isfield(evt.file_info, 'filename')
            filename = evt.file_info.filename;
            if length(filename) >= 4
                try
                    event_year = str2double(filename(1:4));
                    if ~isnan(event_year) && event_year >= 2000 && event_year <= 2100
                        year_key = sprintf('Y%d', event_year);
                    end
                catch
                end
            end
        end
        
        if isKey(group_year_events_map, year_key)
            year_events = group_year_events_map(year_key);
            year_events{end+1} = evt;
            group_year_events_map(year_key) = year_events;
        end
    end
    
    % 그룹별 figure 생성 (연도별 서브플롯)
    num_years = length(available_year_keys);
    fig_group = figure('Name', sprintf('Current-Time Plot - %s', group_key), ...
                       'Position', [100 + g_idx*50, 100 + g_idx*50, 1600, 400*ceil(num_years/2)], ...
                       'Visible', 'on');
    
    for y_idx = 1:num_years
        year_key = available_year_keys{y_idx};
        year_events = group_year_events_map(year_key);
        
        % 서브플롯 생성
        ax = subplot(ceil(num_years/2), 2, y_idx);
        hold(ax, 'on');
        
        num_events = length(year_events);
        if num_events == 0
            text(0.5, 0.5, sprintf('No events for %s', year_key), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 14, 'FontWeight', 'bold', 'Units', 'normalized');
        else
            % 모든 이벤트를 같은 서브플롯에 표시
            num_plotted = 0;
            for e_idx = 1:num_events
                evt = year_events{e_idx};
                
                if isfield(evt, 'current_seq_cell_A') && ~isempty(evt.current_seq_cell_A) && ...
                   isfield(evt, 'time_seq_datetime') && ~isempty(evt.time_seq_datetime)
                    
                    current_seq = evt.current_seq_cell_A;
                    time_seq = evt.time_seq_datetime;
                    
                    % 상대 시간 계산
                    if isdatetime(time_seq(1))
                        time_start = time_seq(1);
                        t_rel = seconds(time_seq - time_start);
                    elseif isduration(time_seq(1))
                        t_rel = seconds(time_seq - time_seq(1));
                    else
                        t_rel = (1:length(time_seq))' - 1;  % 1초 간격 가정
                    end
                    
                    plot(ax, t_rel, current_seq, '-', 'LineWidth', 1.0);
                    num_plotted = num_plotted + 1;
                end
            end
        end
        
        % 연도 레이블
        year_str = strrep(year_key, 'Y', '');
        title(ax, sprintf('Year %s (n=%d events)', year_str, num_events), ...
              'FontSize', 12, 'FontWeight', 'bold');
        xlabel(ax, 'Time (s)', 'FontSize', 10);
        ylabel(ax, 'Current (A)', 'FontSize', 10);
        grid(ax, 'on');
        hold(ax, 'off');
    end
    
    % 전체 제목 (소수점 포함하여 표시)
    group_title_safe = strrep(group_key, '_', ' ');
    sgtitle(sprintf('Current vs Time - %s', group_title_safe), ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    % 저장 (전류-시간 그래프)
    group_filename_safe = strrep(group_key, '.', '_');
    group_filename_safe = strrep(group_filename_safe, '-', '_N');
    group_filename_safe = strrep(group_filename_safe, ' ', '_');
    savePath = fullfile(groupCheckDir, sprintf('CurrentTime_Group_%s.fig', group_filename_safe));
    saveas(fig_group, savePath);
    fprintf('    Saved: %s\n', savePath);
    close(fig_group);
    
    % 전압-시간 그래프 생성
    fig_group_voltage = figure('Name', sprintf('Voltage-Time Plot - %s', group_key), ...
                               'Position', [100 + g_idx*50 + 50, 100 + g_idx*50, 1600, 400*ceil(num_years/2)], ...
                               'Visible', 'on');
    
    for y_idx = 1:num_years
        year_key = available_year_keys{y_idx};
        year_events = group_year_events_map(year_key);
        
        % 서브플롯 생성
        ax = subplot(ceil(num_years/2), 2, y_idx);
        hold(ax, 'on');
        
        num_events = length(year_events);
        if num_events == 0
            text(0.5, 0.5, sprintf('No events for %s', year_key), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 14, 'FontWeight', 'bold', 'Units', 'normalized');
        else
            % 모든 이벤트를 같은 서브플롯에 표시
            num_plotted = 0;
            for e_idx = 1:num_events
                evt = year_events{e_idx};
                
                if isfield(evt, 'voltage_seq_cell_V') && ~isempty(evt.voltage_seq_cell_V) && ...
                   isfield(evt, 'time_seq_datetime') && ~isempty(evt.time_seq_datetime)
                    
                    voltage_seq = evt.voltage_seq_cell_V;
                    time_seq = evt.time_seq_datetime;
                    
                    % 상대 시간 계산
                    if isdatetime(time_seq(1))
                        time_start = time_seq(1);
                        t_rel = seconds(time_seq - time_start);
                    elseif isduration(time_seq(1))
                        t_rel = seconds(time_seq - time_seq(1));
                    else
                        t_rel = (1:length(time_seq))' - 1;  % 1초 간격 가정
                    end
                    
                    plot(ax, t_rel, voltage_seq, '-', 'LineWidth', 1.0);
                    num_plotted = num_plotted + 1;
                end
            end
        end
        
        % 연도 레이블
        year_str = strrep(year_key, 'Y', '');
        title(ax, sprintf('Year %s (n=%d events)', year_str, num_events), ...
              'FontSize', 12, 'FontWeight', 'bold');
        xlabel(ax, 'Time (s)', 'FontSize', 10);
        ylabel(ax, 'Voltage (V)', 'FontSize', 10);
        grid(ax, 'on');
        hold(ax, 'off');
    end
    
    % 전체 제목 (소수점 포함하여 표시)
    sgtitle(sprintf('Voltage vs Time - %s', group_title_safe), ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    % 저장 (전압-시간 그래프)
    savePath_voltage = fullfile(groupCheckDir, sprintf('VoltageTime_Group_%s.fig', group_filename_safe));
    saveas(fig_group_voltage, savePath_voltage);
    fprintf('    Saved: %s\n', savePath_voltage);
    close(fig_group_voltage);
end

fprintf('Group-wise plots saved to: %s\n', groupCheckDir);

%% Calculate Rchg for Each Event (Only Qualified Groups)
fprintf('\n=== Calculating Rchg for Each Event (Qualified Groups Only) ===\n');

resultStructName = sprintf('FieldData_%s_%s', targetRack, targetEventType);
eval(sprintf('%s = struct();', resultStructName));

fprintf('Processing %d qualified groups...\n', length(group_keys));

for g_idx = 1:length(group_keys)
    group_key = group_keys{g_idx};
    group_events = qualified_groups_map(group_key);
    
    fprintf('  Group %d/%d: %s (%d events)\n', g_idx, length(group_keys), group_key, length(group_events));
    
    % 그룹별 구조체 초기화 (안전한 필드명으로 변환)
    group_key_safe = strrep(group_key, '.', '_');
    group_key_safe = strrep(group_key_safe, '-', '_N');  % 마이너스 기호를 _N으로 변환 (예: I-5 -> I_N5)
    group_key_safe = strrep(group_key_safe, ' ', '_');     % 공백을 언더스코어로 변환
    eval(sprintf('if ~isfield(%s, ''%s''); %s.(''%s'') = struct(); end', ...
        resultStructName, group_key_safe, resultStructName, group_key_safe));
    
    % 연도별로 분류
    year_groups = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for e_idx = 1:length(group_events)
        evt = group_events{e_idx};
        
        % 연도 추출 (기존 AnalyzeResistanceFeatures.m 로직 참고)
        year_key = 'Unknown';
        if isfield(evt, 'time_seq_datetime') && ~isempty(evt.time_seq_datetime)
            time_point = evt.time_seq_datetime(1);
            if isdatetime(time_point)
                event_year = year(time_point);
                year_key = sprintf('Y%d', event_year);
            elseif isduration(time_point)
                % duration 형식인 경우 file_info에서 연도 추출
                if isfield(evt, 'file_info') && isfield(evt.file_info, 'filename')
                    filename = evt.file_info.filename;
                    if length(filename) >= 4
                        try
                            event_year = str2double(filename(1:4));
                            if ~isnan(event_year)
                                year_key = sprintf('Y%d', event_year);
                            end
                        catch
                        end
                    end
                end
            elseif isnumeric(time_point)
                % 숫자 형식인 경우 (datenum 등)
                try
                    if time_point > 730000  % datenum 형식 (2000년 이후)
                        dt = datetime(time_point, 'ConvertFrom', 'datenum');
                        event_year = year(dt);
                        year_key = sprintf('Y%d', event_year);
                    end
                catch
                end
            end
        end
        
        % file_info에서 연도 추출 시도 (time_seq_datetime이 없는 경우)
        if strcmp(year_key, 'Unknown') && isfield(evt, 'file_info') && isfield(evt.file_info, 'filename')
            filename = evt.file_info.filename;
            if length(filename) >= 4
                try
                    event_year = str2double(filename(1:4));
                    if ~isnan(event_year) && event_year >= 2000 && event_year <= 2100
                        year_key = sprintf('Y%d', event_year);
                    end
                catch
                end
            end
        end
        
        if ~isKey(year_groups, year_key)
            year_groups(year_key) = {};
        end
        year_events = year_groups(year_key);
        year_events{end+1} = evt;
        year_groups(year_key) = year_events;
    end
    
    % 연도별로 이벤트 저장 및 Rchg 계산
    year_keys = keys(year_groups);
    for y_idx = 1:length(year_keys)
        year_key = year_keys{y_idx};
        year_events = year_groups(year_key);
        
        year_key_safe = year_key;
        if ~startsWith(year_key_safe, 'Y')
            year_key_safe = ['Y' year_key_safe];
        end
        
        eval(sprintf('if ~isfield(%s.(''%s''), ''%s''); %s.(''%s'').(''%s'') = struct(); end', ...
            resultStructName, group_key_safe, year_key_safe, resultStructName, group_key_safe, year_key_safe));
        
        for e_idx = 1:length(year_events)
            evt = year_events{e_idx};
            
            % Rchg 계산 (1s, 3s, 5s, 10s, 30s, 60s)
            [R_chg_1s, R_chg_3s, R_chg_5s, R_chg_10s, R_chg_30s, R_chg_60s] = calculate_dynamic_resistance_v4(evt);
            
            % 이벤트 저장
            evtName = sprintf('event%d', e_idx);
            eval(sprintf('%s.(''%s'').(''%s'').(''%s'').R_chg_1s = R_chg_1s;', ...
                resultStructName, group_key_safe, year_key_safe, evtName));
            eval(sprintf('%s.(''%s'').(''%s'').(''%s'').R_chg_3s = R_chg_3s;', ...
                resultStructName, group_key_safe, year_key_safe, evtName));
            eval(sprintf('%s.(''%s'').(''%s'').(''%s'').R_chg_5s = R_chg_5s;', ...
                resultStructName, group_key_safe, year_key_safe, evtName));
            eval(sprintf('%s.(''%s'').(''%s'').(''%s'').R_chg_10s = R_chg_10s;', ...
                resultStructName, group_key_safe, year_key_safe, evtName));
            eval(sprintf('%s.(''%s'').(''%s'').(''%s'').R_chg_30s = R_chg_30s;', ...
                resultStructName, group_key_safe, year_key_safe, evtName));
            eval(sprintf('%s.(''%s'').(''%s'').(''%s'').R_chg_60s = R_chg_60s;', ...
                resultStructName, group_key_safe, year_key_safe, evtName));
            
            % 원본 데이터 저장
            eval(sprintf('%s.(''%s'').(''%s'').(''%s'').voltage_seq = evt.voltage_seq_cell_V;', ...
                resultStructName, group_key_safe, year_key_safe, evtName));
            eval(sprintf('%s.(''%s'').(''%s'').(''%s'').current_seq = evt.current_seq_cell_A;', ...
                resultStructName, group_key_safe, year_key_safe, evtName));
            eval(sprintf('%s.(''%s'').(''%s'').(''%s'').soc_seq = evt.soc_seq_pct;', ...
                resultStructName, group_key_safe, year_key_safe, evtName));
            eval(sprintf('%s.(''%s'').(''%s'').(''%s'').time_seq = evt.time_seq_datetime;', ...
                resultStructName, group_key_safe, year_key_safe, evtName));
            
            % 타임스탬프 저장
            if isfield(evt, 'time_seq_datetime') && ~isempty(evt.time_seq_datetime)
                time_point = evt.time_seq_datetime(1);
                if isdatetime(time_point)
                    eval(sprintf('%s.(''%s'').(''%s'').(''%s'').timestamp = time_point;', ...
                        resultStructName, group_key_safe, year_key_safe, evtName));
                end
            end
        end
    end
end

%% Save Results
fprintf('\n=== Saving Results ===\n');
savePath = fullfile(outputDir, sprintf('%s_Events.mat', resultStructName));
eval(sprintf('save(''%s'', ''%s'');', savePath, resultStructName));
fprintf('Results saved to: %s\n', savePath);

fprintf('\n=== Event Analysis Complete ===\n');

%% Helper Function: Calculate Dynamic Resistance
function [R_chg_1s, R_chg_3s, R_chg_5s, R_chg_10s, R_chg_30s, R_chg_60s] = calculate_dynamic_resistance_v4(event_struct)
    % EventDetection.m에서 이미 SG 필터(Savitzky-Golay, 윈도우 5, 차수 1)로 필터링된 데이터를 사용
    R_chg_1s = NaN; R_chg_3s = NaN; R_chg_5s = NaN; R_chg_10s = NaN; R_chg_30s = NaN; R_chg_60s = NaN;
    V = event_struct.voltage_seq_cell_V; I = event_struct.current_seq_cell_A;
    if length(V) < 2, return; end
    
    % 시간 시퀀스에서 상대 시간 계산
    time_seq = event_struct.time_seq_datetime;
    if isdatetime(time_seq(1))
        time_start = time_seq(1);
        t_rel = seconds(time_seq - time_start);
    elseif isduration(time_seq(1))
        t_rel = seconds(time_seq - time_seq(1));
    else
        t_rel = time_seq - time_seq(1);
    end
    
    % 1초 후 저항 계산
    idx_1s = find(t_rel >= 1.0, 1);
    if ~isempty(idx_1s) && idx_1s > 1 && idx_1s <= length(V) && idx_1s <= length(I)
        delta_V_1s = V(idx_1s) - V(1);
        delta_I_1s = I(idx_1s) - I(1);
        if abs(delta_I_1s) > 1e-6
            R_chg_1s = (delta_V_1s / delta_I_1s) * 1000;
        end
    end
    
    % 3초 후 저항 계산
    idx_3s = find(t_rel >= 3.0, 1);
    if ~isempty(idx_3s) && idx_3s > 1 && idx_3s <= length(V) && idx_3s <= length(I)
        delta_V_3s = V(idx_3s) - V(1);
        delta_I_3s = I(idx_3s) - I(1);
        if abs(delta_I_3s) > 1e-6
            R_chg_3s = (delta_V_3s / delta_I_3s) * 1000;
        end
    end
    
    % 5초 후 저항 계산
    idx_5s = find(t_rel >= 5.0, 1);
    if ~isempty(idx_5s) && idx_5s > 1 && idx_5s <= length(V) && idx_5s <= length(I)
        delta_V_5s = V(idx_5s) - V(1);
        delta_I_5s = I(idx_5s) - I(1);
        if abs(delta_I_5s) > 1e-6
            R_chg_5s = (delta_V_5s / delta_I_5s) * 1000;
        end
    end
    
    % 10초 후 저항 계산
    idx_10s = find(t_rel >= 10.0, 1);
    if ~isempty(idx_10s) && idx_10s > 1 && idx_10s <= length(V) && idx_10s <= length(I)
        delta_V_10s = V(idx_10s) - V(1);
        delta_I_10s = I(idx_10s) - I(1);
        if abs(delta_I_10s) > 1e-6
            R_chg_10s = (delta_V_10s / delta_I_10s) * 1000;
        end
    end
    
    % 30초 후 저항 계산
    idx_30s = find(t_rel >= 30.0, 1);
    if ~isempty(idx_30s) && idx_30s > 1 && idx_30s <= length(V) && idx_30s <= length(I)
        delta_V_30s = V(idx_30s) - V(1);
        delta_I_30s = I(idx_30s) - I(1);
        if abs(delta_I_30s) > 1e-6
            R_chg_30s = (delta_V_30s / delta_I_30s) * 1000;
        end
    end
    
    % 60초 후 저항 계산
    idx_60s = find(t_rel >= 60.0, 1);
    if ~isempty(idx_60s) && idx_60s > 1 && idx_60s <= length(V) && idx_60s <= length(I)
        delta_V_60s = V(idx_60s) - V(1);
        delta_I_60s = I(idx_60s) - I(1);
        if abs(delta_I_60s) > 1e-6
            R_chg_60s = (delta_V_60s / delta_I_60s) * 1000;
        end
    end
    
    % 음수 저항 제거
    if ~isnan(R_chg_1s) && R_chg_1s < 0, R_chg_1s = NaN; end
    if ~isnan(R_chg_3s) && R_chg_3s < 0, R_chg_3s = NaN; end
    if ~isnan(R_chg_5s) && R_chg_5s < 0, R_chg_5s = NaN; end
    if ~isnan(R_chg_10s) && R_chg_10s < 0, R_chg_10s = NaN; end
    if ~isnan(R_chg_30s) && R_chg_30s < 0, R_chg_30s = NaN; end
    if ~isnan(R_chg_60s) && R_chg_60s < 0, R_chg_60s = NaN; end
end

