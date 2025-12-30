%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 파일명: AnalyzeResistanceFeatures.m (v4.3 - 최종 안정화 및 시각화 개선)
% 기능:
% - [핵심 수정] boxplot의 'GroupOrder' 타입 오류 해결
% - [유지] 정확한 그룹 선정, 디버깅 출력, 저항 계산 로직
% MATLAB R2022b 이상 권장
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% 1. 분석 설정
% =========================================================================
% 스크립트 위치를 기준으로 상대 경로 설정
scriptDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(scriptDir, 'EventsResults');  % 현재 디렉토리의 EventsResults 폴더
saveDir = fullfile(resultsDir, 'Rchg_ver01');  % EventsResults/Rchg_ver01 폴더에 저장
rawEventsFile = 'all_events_raw_cell_level.mat';

targetRack = 'Rack01';
targetEventType = 'charge';
min_events_per_year = 1;  % 모든 연도에서 최소 이벤트 수 (낮춤: 3 -> 1)
min_event_duration_sec = 30;  % 최소 이벤트 길이 (초)
min_delta_I = 64*0.10;

% 그룹화 파라미터
voltage_bin_size = 0.01;  % 전압 bin 크기 (V) - 0.01V 간격
current_bin_size = 5;      % 전류 bin 크기 (A) - 10A 간격
soc_bin_size = 5;          % SOC bin 크기 (%) - 5% 간격

savePlots = true;
% =========================================================================

%% 2. 데이터 로드 및 그룹화
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

% EventType별 폴더 생성 (charge 또는 discharge)
eventTypeDir = fullfile(saveDir, targetEventType);
if ~exist(eventTypeDir, 'dir')
    mkdir(eventTypeDir);
end

fprintf('Loading raw events data...\n');
try, load(fullfile(resultsDir, rawEventsFile), 'all_events'); catch ME, error('Failed to load %s. Error: %s', rawEventsFile, ME.message); end
fprintf('Data loaded.\n');

if ~isfield(all_events, targetRack), error('Target Rack "%s" not found.', targetRack); end
rackData = all_events.(targetRack);
available_year_keys = fieldnames(rackData);
available_year_keys = available_year_keys(startsWith(available_year_keys, 'Y'));

if isempty(available_year_keys), error('No year data found for Rack "%s".', targetRack); end
fprintf('Found data for years: %s\n', strjoin(available_year_keys, ', '));

% 모든 이벤트 수집 및 그룹화
fprintf('Collecting and grouping events...\n');
all_collected_events = {};
for y_idx = 1:length(available_year_keys)
    year_key = available_year_keys{y_idx};
    if ~isstruct(rackData.(year_key)), continue; end
    month_keys = fieldnames(rackData.(year_key));
    for m_idx = 1:length(month_keys)
        month_key = month_keys{m_idx};
        if isempty(rackData.(year_key).(month_key)), continue; end
        month_events = rackData.(year_key).(month_key);
        for e_idx = 1:length(month_events)
            evt = month_events{e_idx};
            if isempty(evt) || ~isstruct(evt), continue; end
            if ~isfield(evt, 'type') || ~strcmp(evt.type, targetEventType), continue; end
            
            % 이벤트 길이 확인
            if isfield(evt, 'time_seq_datetime') && ~isempty(evt.time_seq_datetime)
                time_seq = evt.time_seq_datetime;
                if isdatetime(time_seq(1))
                    event_duration = seconds(time_seq(end) - time_seq(1));
                elseif isduration(time_seq(1))
                    event_duration = seconds(time_seq(end) - time_seq(1));
                else
                    try
                        if isnumeric(time_seq(1))
                            event_duration = time_seq(end) - time_seq(1);
                        else
                            continue;
                        end
                    catch
                        continue;
                    end
                end
                if event_duration < min_event_duration_sec, continue; end
            else
                continue;
            end
            
            all_collected_events{end+1} = evt;
        end
    end
end
fprintf('Collected %d events (after filtering by duration >= %d sec).\n', length(all_collected_events), min_event_duration_sec);

% 그룹화: 전압(시작), 평균 전류, 평균 power
fprintf('Grouping events by Voltage, Average Current, and Average Power...\n');
grouped_events_map = containers.Map('KeyType', 'char', 'ValueType', 'any');

% 디버깅: 필터링 단계별 카운터
filter_stats = struct();
filter_stats.total = length(all_collected_events);
filter_stats.no_voltage = 0;
filter_stats.no_current = 0;
filter_stats.length_too_short_21 = 0;
filter_stats.delta_I_too_small = 0;
filter_stats.length_too_short_12 = 0;
filter_stats.current_variation_too_high = 0;
filter_stats.passed = 0;

for i = 1:length(all_collected_events)
    evt = all_collected_events{i};
    
    % 시작 전압
    if ~isfield(evt, 'voltage_seq_cell_V') || isempty(evt.voltage_seq_cell_V)
        filter_stats.no_voltage = filter_stats.no_voltage + 1;
        continue;
    end
    start_voltage = evt.voltage_seq_cell_V(1);
    voltage_bin = floor(start_voltage / voltage_bin_size) * voltage_bin_size;
    
    % delta_I 필터링: 시작점부터 인덱스 20 이내의 전류 변화가 6.4A 이상인지 확인
    if ~isfield(evt, 'current_seq_cell_A') || isempty(evt.current_seq_cell_A)
        filter_stats.no_current = filter_stats.no_current + 1;
        continue;
    end
    current_seq = evt.current_seq_cell_A;
    if length(current_seq) < 21  % 최소 21개 포인트 필요 (시작점 + 20)
        filter_stats.length_too_short_21 = filter_stats.length_too_short_21 + 1;
        continue;
    end
    
    % 시작점과 인덱스 20 사이의 전류 변화 확인: I(20) - I(1) > 6.4A
    idx_20 = min(20, length(current_seq));  % 인덱스 20
    delta_I = current_seq(idx_20) - current_seq(1);
    if delta_I < min_delta_I  
        filter_stats.delta_I_too_small = filter_stats.delta_I_too_small + 1;
        continue;
    end
    
    % 평균 전류 (10번째 인덱스부터 end-10까지) 및 전류 변동성 필터링
    if length(current_seq) < 20  % 최소 20개 포인트 필요 (10부터 end-10까지)
        filter_stats.length_too_short_12 = filter_stats.length_too_short_12 + 1;
        continue;
    end
    stable_current = abs(current_seq(10:end-10));
    avg_current = mean(stable_current);
    
    % 전류 변동성 필터링: 변동률이 10% 이내인 경우만 통과 (완화: 10% -> 30%)
    if avg_current > 0
        current_max = max(stable_current);
        current_min = min(stable_current);
        current_variation = (current_max - current_min) / avg_current * 100;
        if current_variation > 5.0  % 30% 초과면 제외 (완화: 10% -> 30%)
            filter_stats.current_variation_too_high = filter_stats.current_variation_too_high + 1;
            continue;
        end
    end
    
    current_bin = floor(avg_current / current_bin_size) * current_bin_size;
    
    % SOC bin 계산 (시작 SOC 사용)
    soc_bin = NaN;
    if isfield(evt, 'soc_seq_pct') && ~isempty(evt.soc_seq_pct)
        start_soc = evt.soc_seq_pct(1);
        if ~isnan(start_soc)
            soc_bin = floor(start_soc / soc_bin_size) * soc_bin_size;
        end
    end
    
    % 그룹 키 생성 (전압, 전류, SOC 사용)
    % 전압은 0.01V 간격이므로 소수점 2자리까지 표시
    if ~isnan(soc_bin)
        group_key = sprintf('V%.2f_I%d_SOC%d', voltage_bin, round(current_bin), round(soc_bin));
    else
        % SOC가 없으면 기존 방식 사용
        group_key = sprintf('V%.2f_I%d', voltage_bin, round(current_bin));
    end
    
    if ~isKey(grouped_events_map, group_key)
        grouped_events_map(group_key) = {evt};
    else
        current_events = grouped_events_map(group_key);
        current_events{end+1} = evt;
        grouped_events_map(group_key) = current_events;
    end
    filter_stats.passed = filter_stats.passed + 1;
end

% 필터링 통계 출력
fprintf('\n--- Filtering Statistics ---\n');
fprintf('Total events: %d\n', filter_stats.total);
fprintf('  - No voltage data: %d\n', filter_stats.no_voltage);
fprintf('  - No current data: %d\n', filter_stats.no_current);
fprintf('  - Length < 21 points: %d\n', filter_stats.length_too_short_21);
fprintf('  - delta_I <= %d: %d\n', min_delta_I, filter_stats.delta_I_too_small);
fprintf('  - Length < 12 points: %d\n', filter_stats.length_too_short_12);
fprintf('  - Current variation > 5%%: %d\n', filter_stats.current_variation_too_high);
fprintf('  - Passed all filters: %d\n', filter_stats.passed);
fprintf('Total groups created: %d\n', length(grouped_events_map));
fprintf('\n');

% 연도별 통계 계산
fprintf('Calculating statistics per group and year...\n');
group_list = {};
group_keys_all = keys(grouped_events_map);
for g_idx = 1:length(group_keys_all)
    group_key = group_keys_all{g_idx};
    events_in_group = grouped_events_map(group_key);
    
    % 연도별 카운트
    stats = struct('total_count', length(events_in_group));
    for yk_init = available_year_keys', stats.(yk_init{1}) = 0; end
    
    for e_idx = 1:length(events_in_group)
        evt = events_in_group{e_idx};
        % 연도 추출
        if isfield(evt, 'time_seq_datetime') && ~isempty(evt.time_seq_datetime)
            time_point = evt.time_seq_datetime(1);
            if isdatetime(time_point)
                event_year = year(time_point);
            elseif isduration(time_point)
                if isfield(evt, 'file_info') && isfield(evt.file_info, 'filename')
                    filename = evt.file_info.filename;
                    if length(filename) >= 4
                        try
                            event_year = str2double(filename(1:4));
                        catch
                            continue;
                        end
                    else
                        continue;
                    end
                else
                    continue;
                end
            else
                continue;
            end
            year_key = sprintf('Y%d', event_year);
            if isfield(stats, year_key)
                stats.(year_key) = stats.(year_key) + 1;
            end
        end
    end
    
    group_list{end+1, 1} = group_key;
    group_list{end, 2} = stats;
end

% 연도별 최소 이벤트 수 조건 확인
fprintf('Filtering groups with >= %d events for each available year...\n', min_events_per_year);
fprintf('Total groups before filtering: %d\n', size(group_list, 1));

% 디버깅: 전류값별 그룹 분포 확인 (전체 그룹)
fprintf('\n--- Current bin distribution (all groups - before filtering) ---\n');
current_bins_all = containers.Map('KeyType', 'char', 'ValueType', 'double');
for g_idx = 1:size(group_list, 1)
    group_key = group_list{g_idx, 1};
    % 그룹 키에서 전류값 추출 (예: "V3.79_I10_SOC60" -> "I10_SOC60")
    if contains(group_key, '_I')
        parts = strsplit(group_key, '_I');
        if length(parts) >= 2
            current_bin_str = ['I' parts{2}];
            if isKey(current_bins_all, current_bin_str)
                current_bins_all(current_bin_str) = current_bins_all(current_bin_str) + 1;
            else
                current_bins_all(current_bin_str) = 1;
            end
        end
    end
end
current_bin_keys = keys(current_bins_all);
current_bin_counts = zeros(length(current_bin_keys), 1);
for k = 1:length(current_bin_keys)
    current_bin_counts(k) = current_bins_all(current_bin_keys{k});
end
[~, sort_idx_current] = sort(current_bin_counts, 'descend');
fprintf('Current Bin | Number of Groups (All)\n');
fprintf('-------------------------------\n');
for k = 1:length(current_bin_keys)
    idx = sort_idx_current(k);
    fprintf('  %-8s | %d groups\n', current_bin_keys{idx}, current_bin_counts(idx));
end
fprintf('\n');

% 디버깅: 상위 그룹들의 연도별 통계 출력
if size(group_list, 1) > 0
    fprintf('\n--- Top 10 groups by total count (for debugging) ---\n');
    % total_count 추출
    total_counts = zeros(size(group_list, 1), 1);
    for debug_idx = 1:size(group_list, 1)
        stats_debug = group_list{debug_idx, 2};
        total_counts(debug_idx) = stats_debug.total_count;
    end
    [~, sort_idx_debug] = sort(total_counts, 'descend');
    for debug_idx = 1:min(10, size(group_list, 1))
        idx = sort_idx_debug(debug_idx);
        group_key_debug = group_list{idx, 1};
        stats_debug = group_list{idx, 2};
        fprintf('  Group "%s": Total=%d, ', group_key_debug, stats_debug.total_count);
        year_counts_str = {};
        for yk = available_year_keys'
            year_counts_str{end+1} = sprintf('%s=%d', yk{1}, stats_debug.(yk{1}));
        end
        fprintf('%s\n', strjoin(year_counts_str, ', '));
    end
    fprintf('\n');
else
    fprintf('No groups found after grouping. Check filtering criteria.\n');
end

qualified_groups = {};
for i = 1:size(group_list, 1)
    group_key = group_list{i, 1};
    stats = group_list{i, 2};
    is_qualified = true;
    for y_idx = 1:length(available_year_keys)
        year_key = available_year_keys{y_idx};
        if stats.(year_key) < min_events_per_year
            is_qualified = false;
            fprintf('  - Group "%s" disqualified. Reason: Failed at %s (count: %d)\n', group_key, year_key, stats.(year_key));
            break;
        end
    end
    if is_qualified
        qualified_groups{end+1, 1} = group_key;
        qualified_groups{end, 2} = stats.total_count;
    end
end

if isempty(qualified_groups)
    fprintf('\n=== No groups met the criteria ===\n');
    fprintf('Available years: %s\n', strjoin(available_year_keys, ', '));
    fprintf('Required: >= %d events per year\n', min_events_per_year);
    fprintf('\nSuggestions:\n');
    fprintf('  1. Lower min_events_per_year (currently: %d)\n', min_events_per_year);
    fprintf('  2. Check if data exists for all years\n');
    fprintf('  3. Review filtering criteria (delta_I, current_variation, etc.)\n');
    error('No groups met the criteria. Try lowering min_events_per_year.');
end

[~, sorted_idx] = sort(cell2mat(qualified_groups(:,2)), 'descend');
top_groups_info = qualified_groups(sorted_idx, :);
% 모든 qualified 그룹 처리 (numTopGroups 제한 제거)
top_group_keys = top_groups_info(:,1);

% 디버깅: Qualified 그룹의 전류값별 분포 확인
fprintf('\n--- Current bin distribution (qualified groups - after filtering) ---\n');
current_bins_qualified = containers.Map('KeyType', 'char', 'ValueType', 'double');
for g_idx = 1:length(top_group_keys)
    group_key = top_group_keys{g_idx};
    % 그룹 키에서 전류값 추출 (예: "V3.79_I10_SOC60" -> "I10_SOC60")
    if contains(group_key, '_I')
        parts = strsplit(group_key, '_I');
        if length(parts) >= 2
            current_bin_str = ['I' parts{2}];
            if isKey(current_bins_qualified, current_bin_str)
                current_bins_qualified(current_bin_str) = current_bins_qualified(current_bin_str) + 1;
            else
                current_bins_qualified(current_bin_str) = 1;
            end
        end
    end
end
current_bin_keys_qual = keys(current_bins_qualified);
current_bin_counts_qual = zeros(length(current_bin_keys_qual), 1);
for k = 1:length(current_bin_keys_qual)
    current_bin_counts_qual(k) = current_bins_qualified(current_bin_keys_qual{k});
end
[~, sort_idx_current_qual] = sort(current_bin_counts_qual, 'descend');
fprintf('Current Bin | Number of Groups (Qualified)\n');
fprintf('-------------------------------\n');
for k = 1:length(current_bin_keys_qual)
    idx = sort_idx_current_qual(k);
    fprintf('  %-8s | %d groups\n', current_bin_keys_qual{idx}, current_bin_counts_qual(idx));
end
fprintf('\n');

fprintf('\nAll %d qualified groups selected for analysis:\n', length(top_group_keys));
for i = 1:length(top_group_keys), fprintf('  %d. %s (Total: %d events)\n', i, top_group_keys{i}, top_groups_info{i,2}); end
fprintf('\n=== Starting analysis for ALL %d groups ===\n', length(top_group_keys));


%% 3. 결과 저장용 구조체 초기화
results_data = struct();

%% 4. 각 Top 그룹에 대해 분석 및 시각화 수행
for g = 1:length(top_group_keys)
    currentTargetGroupKey = top_group_keys{g};
    fprintf('\n========== Processing Group %d/%d: %s ==========\n', g, length(top_group_keys), currentTargetGroupKey);
    
    % 그룹별 결과 구조체 초기화
    group_key_safe = strrep(currentTargetGroupKey, '.', '_');  % 구조체 필드명에 . 사용 불가
    results_data.(group_key_safe) = struct();
    
    % Current bin 추출 (예: "V3.79_I10_SOC60" -> "I10")
    current_bin_str = '';
    if contains(currentTargetGroupKey, '_I')
        parts = strsplit(currentTargetGroupKey, '_I');
        if length(parts) >= 2
            % "10_SOC60" -> "I10"
            soc_part = parts{2};
            if contains(soc_part, '_SOC')
                current_num = strsplit(soc_part, '_SOC');
                current_bin_str = ['I' current_num{1}];
            else
                current_bin_str = ['I' soc_part];
            end
        end
    end
    
    % Current bin 폴더 생성 (EventType 폴더 안에)
    if ~isempty(current_bin_str)
        current_bin_folder = fullfile(eventTypeDir, current_bin_str);
        if ~exist(current_bin_folder, 'dir')
            mkdir(current_bin_folder);
        end
        % 그룹별 폴더 생성 (Current bin 폴더 안에)
        group_folder = fullfile(current_bin_folder, group_key_safe);
    else
        % Current bin을 추출할 수 없으면 기존 방식 사용
        group_folder = fullfile(eventTypeDir, group_key_safe);
    end
    if ~exist(group_folder, 'dir')
        mkdir(group_folder);
    end
    
    % 데이터 추출 (그룹화된 이벤트에서)
    all_target_events = grouped_events_map(currentTargetGroupKey);
    
    % 연도별 카운트
    year_event_counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for yk = available_year_keys', year_event_counts(yk{1}) = 0; end
    
    for e_idx = 1:length(all_target_events)
        evt = all_target_events{e_idx};
        if isfield(evt, 'time_seq_datetime') && ~isempty(evt.time_seq_datetime)
            time_point = evt.time_seq_datetime(1);
            if isdatetime(time_point)
                event_year = year(time_point);
                year_key = sprintf('Y%d', event_year);
                if isKey(year_event_counts, year_key)
                    year_event_counts(year_key) = year_event_counts(year_key) + 1;
                end
            elseif isduration(time_point)
                if isfield(evt, 'file_info') && isfield(evt.file_info, 'filename')
                    filename = evt.file_info.filename;
                    if length(filename) >= 4
                        try
                            event_year = str2double(filename(1:4));
                            year_key = sprintf('Y%d', event_year);
                            if isKey(year_event_counts, year_key)
                                year_event_counts(year_key) = year_event_counts(year_key) + 1;
                            end
                        catch
                        end
                    end
                end
            end
        end
    end
    
    fprintf('Events extracted by year:\n');
    for yk = available_year_keys'
        fprintf('  %s: %d events\n', yk{1}, year_event_counts(yk{1}));
    end

    % 피처 계산 (타임스탬프, R_chg_1s, R_chg_5s, R_chg_10s, R_chg_30s, 이벤트 시퀀스 데이터 추출)
    num_events = length(all_target_events);
    timestamps = NaT(num_events, 1);
    R_chg_1s_vals = NaN(num_events, 1);
    R_chg_5s_vals = NaN(num_events, 1);
    R_chg_10s_vals = NaN(num_events, 1);
    R_chg_30s_vals = NaN(num_events, 1);
    
    % 이벤트 시퀀스 데이터 저장 (시간-전류 시각화용)
    event_sequences = cell(num_events, 1);  % 각 이벤트의 {time_relative_sec, current, year}

    for i = 1:num_events
        evt = all_target_events{i};
        
        % 이벤트 길이 확인 (30초 이상 필터링)
        if isfield(evt, 'time_seq_datetime') && ~isempty(evt.time_seq_datetime)
            time_seq = evt.time_seq_datetime;
            if isdatetime(time_seq(1))
                event_duration = seconds(time_seq(end) - time_seq(1));
            elseif isduration(time_seq(1))
                event_duration = seconds(time_seq(end) - time_seq(1));
            else
                % 숫자 형식인 경우
                try
                    if isnumeric(time_seq(1))
                        event_duration = time_seq(end) - time_seq(1);
                    else
                        continue;  % 길이 확인 불가능하면 스킵
                    end
                catch
                    continue;
                end
            end
            
            % 최소 길이 미만 이벤트는 제외
            if event_duration < min_event_duration_sec
                continue;
            end
        else
            continue;  % 시간 시퀀스가 없으면 스킵
        end
        
        % 이벤트 발생 시간 추출 (R_chg_1s 계산 및 연도 분류용)
        if isfield(evt, 'time_seq_datetime') && ~isempty(evt.time_seq_datetime)
            time_point = evt.time_seq_datetime(1);
            if isdatetime(time_point)
                timestamps(i) = time_point;
            elseif isnumeric(time_point)
                try
                    timestamps(i) = datetime(time_point, 'ConvertFrom', 'datenum');
                catch
                    try
                        timestamps(i) = datetime(time_point, 'ConvertFrom', 'excel');
                    catch
                    end
                end
            elseif isduration(time_point)
                try
                    if isfield(evt, 'file_info') && isfield(evt.file_info, 'filename')
                        filename = evt.file_info.filename;
                        if length(filename) >= 8
                            fileDate = datetime(filename(1:8), 'InputFormat', 'yyyyMMdd');
                            timestamps(i) = fileDate + time_point;
                        end
                    end
                catch
                end
            else
                try
                    if isfield(evt, 'file_info') && isfield(evt.file_info, 'filename')
                        filename = evt.file_info.filename;
                        if length(filename) >= 8
                            fileDate = datetime(filename(1:8), 'InputFormat', 'yyyyMMdd');
                            timestamps(i) = fileDate;
                        end
                    end
                catch
                end
            end
        end
        
        % R_chg 계산 (1s, 5s, 10s, 30s)
        [R_chg_1s, R_chg_5s, R_chg_10s, R_chg_30s] = calculate_dynamic_resistance_v4(evt);
        R_chg_1s_vals(i) = R_chg_1s;
        R_chg_5s_vals(i) = R_chg_5s;
        R_chg_10s_vals(i) = R_chg_10s;
        R_chg_30s_vals(i) = R_chg_30s;
        
        % 디버깅: 처음 몇 개 이벤트만 상세 정보 출력
        if i <= 3 && isnan(R_chg_1s)
            if isfield(evt, 'time_seq_datetime') && ~isempty(evt.time_seq_datetime)
                time_seq = evt.time_seq_datetime;
                if isdatetime(time_seq(1))
                    t_rel = seconds(time_seq - time_seq(1));
                elseif isduration(time_seq(1))
                    t_rel = seconds(time_seq - time_seq(1));
                else
                    t_rel = time_seq - time_seq(1);
                end
                idx_1s = find(t_rel >= 1.0, 1);
                if isfield(evt, 'current_seq_cell_A') && ~isempty(evt.current_seq_cell_A)
                    I = evt.current_seq_cell_A;
                    if ~isempty(idx_1s) && idx_1s > 1 && idx_1s <= length(I)
                        delta_I = abs(I(idx_1s) - I(1));
                        fprintf('  [Debug Event %d] idx_1s=%d, delta_I=%.4f, length(I)=%d\n', i, idx_1s, delta_I, length(I));
                    else
                        fprintf('  [Debug Event %d] idx_1s empty or invalid, length(I)=%d\n', i, length(I));
                    end
                end
            end
        end
        
        % 이벤트 시퀀스 데이터 추출 (시간-전류, 시간-전압 시각화용)
        if isfield(evt, 'time_seq_datetime') && isfield(evt, 'current_seq_cell_A') && ...
           isfield(evt, 'voltage_seq_cell_V') && ~isempty(evt.time_seq_datetime) && ...
           ~isempty(evt.current_seq_cell_A) && ~isempty(evt.voltage_seq_cell_V)
            current_seq = evt.current_seq_cell_A;
            voltage_seq = evt.voltage_seq_cell_V;
            
            % 이벤트 시작점을 0초로 맞춤 (상대 시간 계산)
            if isdatetime(time_seq(1))
                time_start = time_seq(1);
                time_relative_sec = seconds(time_seq - time_start);
            elseif isduration(time_seq(1))
                time_relative_sec = seconds(time_seq - time_seq(1));
            else
                % 숫자나 다른 형식인 경우
                try
                    if isnumeric(time_seq(1))
                        time_relative_sec = time_seq - time_seq(1);
                    else
                        continue;
                    end
                catch
                    continue;
                end
            end
            
            % 연도 추출
            event_year = NaN;
            if ~isnat(timestamps(i))
                event_year = year(timestamps(i));
            elseif isfield(evt, 'file_info') && isfield(evt.file_info, 'filename')
                filename = evt.file_info.filename;
                if length(filename) >= 4
                    try
                        event_year = str2double(filename(1:4));
                    catch
                    end
                end
            end
            
            % SOC 추출 (시작 SOC)
            start_soc = NaN;
            if isfield(evt, 'soc_seq_pct') && ~isempty(evt.soc_seq_pct)
                start_soc = evt.soc_seq_pct(1);
            end
            
            event_sequences{i} = struct('time_sec', time_relative_sec, 'current', current_seq, ...
                'voltage', voltage_seq, 'year', event_year, 'start_soc', start_soc);
        end
    end
    
    
    valid_idx = ~isnat(timestamps) & ~isnan(R_chg_1s_vals);
    
    % 디버깅 정보 출력
    fprintf('Total events: %d\n', num_events);
    fprintf('Valid timestamps: %d\n', sum(~isnat(timestamps)));
    fprintf('Valid R_chg_1s: %d\n', sum(~isnan(R_chg_1s_vals)));
    fprintf('Valid events (both timestamp and R_chg_1s): %d\n', sum(valid_idx));
    
    % 원본 인덱스 저장 (나중에 저장할 때 사용)
    original_indices = find(valid_idx);
    
    feature_table = table(timestamps(valid_idx), R_chg_1s_vals(valid_idx), R_chg_5s_vals(valid_idx), ...
        R_chg_10s_vals(valid_idx), R_chg_30s_vals(valid_idx), ...
        'VariableNames', {'timestamp', 'R_chg_1s', 'R_chg_5s', 'R_chg_10s', 'R_chg_30s'});
    feature_table = sortrows(feature_table, 'timestamp');
    
    % 정렬 후 원본 인덱스도 함께 정렬
    [~, sort_order] = sortrows(feature_table, 'timestamp');
    original_indices = original_indices(sort_order);
    
    % 연도 컬럼 추가
    if height(feature_table) > 0
        feature_table.Year = categorical(year(feature_table.timestamp));
    else
        fprintf('Skipping visualization: feature_table is empty.\n'); continue;
    end
    
    % 유효한 이벤트 시퀀스만 필터링
    valid_sequences = event_sequences(valid_idx);
    
    % 연도별 유효 데이터 확인
    if height(feature_table) > 0
        year_counts = countcats(feature_table.Year);
        unique_years = categories(feature_table.Year);
        fprintf('Valid events by year in feature_table:\n');
        for y_idx = 1:length(unique_years)
            fprintf('  %s: %d events\n', unique_years{y_idx}, year_counts(y_idx));
        end
    end
    
    if height(feature_table) < 20, fprintf('Skipping visualization due to insufficient data (need >= 20, got %d).\n', height(feature_table)); continue; end

    % 연도별 이상치 제거 (IQR 방법) - 모든 R_chg 값에 대해 적용
    fprintf('Removing outliers by year using IQR method (for all R_chg values)...\n');
    outlier_idx = false(height(feature_table), 1);
    unique_years = categories(feature_table.Year);
    r_chg_vars_outlier = {'R_chg_1s', 'R_chg_5s', 'R_chg_10s', 'R_chg_30s'};
    
    for r_var_idx = 1:length(r_chg_vars_outlier)
        var_name = r_chg_vars_outlier{r_var_idx};
        if ~ismember(var_name, feature_table.Properties.VariableNames)
            continue;
        end
        
        fprintf('  Processing outliers for %s...\n', var_name);
        for y_idx = 1:length(unique_years)
            year_val = unique_years{y_idx};
            year_mask = feature_table.Year == year_val;
            year_data = feature_table.(var_name)(year_mask);
            year_data = year_data(~isnan(year_data));
            
            if sum(~isnan(year_data)) > 10  % 최소 10개 이상의 데이터가 있을 때만 이상치 제거
                Q1 = prctile(year_data, 25);
                Q3 = prctile(year_data, 75);
                IQR = Q3 - Q1;
                lower_bound = Q1 - 1.5 * IQR;
                upper_bound = Q3 + 1.5 * IQR;
                
                year_outliers = year_mask & (feature_table.(var_name) < lower_bound | feature_table.(var_name) > upper_bound);
                outlier_idx = outlier_idx | year_outliers;
                if sum(year_outliers) > 0
                    fprintf('    %s: Removed %d outliers (bounds: [%.4f, %.4f])\n', ...
                        year_val, sum(year_outliers), lower_bound, upper_bound);
                end
            end
        end
    end
    
    feature_table = feature_table(~outlier_idx, :);
    fprintf('After outlier removal: %d events remaining.\n', height(feature_table));
    
    if height(feature_table) < 10, fprintf('Skipping visualization due to insufficient data after outlier removal.\n'); continue; end

    % 디버깅 정보 출력
    fprintf('%d valid events processed for this group (after outlier removal).\n', height(feature_table));
    fprintf('--- Analysis Summary for Group: %s ---\n', currentTargetGroupKey);
    fprintf('  Feature    |    Mean    |    Std     |    Min     |    Max     | Valid Pts\n');
    fprintf('---------------------------------------------------------------------------\n');
    stats_1s = [mean(feature_table.R_chg_1s, 'omitnan'), std(feature_table.R_chg_1s, 'omitnan'), min(feature_table.R_chg_1s), max(feature_table.R_chg_1s), sum(~isnan(feature_table.R_chg_1s))];
    stats_5s = [mean(feature_table.R_chg_5s, 'omitnan'), std(feature_table.R_chg_5s, 'omitnan'), min(feature_table.R_chg_5s), max(feature_table.R_chg_5s), sum(~isnan(feature_table.R_chg_5s))];
    stats_10s = [mean(feature_table.R_chg_10s, 'omitnan'), std(feature_table.R_chg_10s, 'omitnan'), min(feature_table.R_chg_10s), max(feature_table.R_chg_10s), sum(~isnan(feature_table.R_chg_10s))];
    stats_30s = [mean(feature_table.R_chg_30s, 'omitnan'), std(feature_table.R_chg_30s, 'omitnan'), min(feature_table.R_chg_30s), max(feature_table.R_chg_30s), sum(~isnan(feature_table.R_chg_30s))];
    fprintf('  R_chg_1s  | %10.4f | %10.4f | %10.4f | %10.4f | %d\n', stats_1s);
    fprintf('  R_chg_5s  | %10.4f | %10.4f | %10.4f | %10.4f | %d\n', stats_5s);
    fprintf('  R_chg_10s | %10.4f | %10.4f | %10.4f | %10.4f | %d\n', stats_10s);
    fprintf('  R_chg_30s | %10.4f | %10.4f | %10.4f | %10.4f | %d\n', stats_30s);
    fprintf('---------------------------------------------------------------------------\n');
    
    %% 연도별 통계 분석 및 p-value 계산
    fprintf('\n========== Yearly Statistical Analysis for Group: %s ==========\n', currentTargetGroupKey);
    unique_years = categories(feature_table.Year);
    num_years = length(unique_years);
    
    % 연도별 기본 통계 (모든 R_chg 값에 대해)
    r_chg_vars_stats = {'R_chg_1s', 'R_chg_5s', 'R_chg_10s', 'R_chg_30s'};
    r_chg_labels_stats = {'R_chg_1s', 'R_chg_5s', 'R_chg_10s', 'R_chg_30s'};
    
    yearly_data = cell(num_years, 1);
    
    for r_stat_idx = 1:length(r_chg_vars_stats)
        var_name = r_chg_vars_stats{r_stat_idx};
        if ~ismember(var_name, feature_table.Properties.VariableNames)
            continue;
        end
        
        fprintf('\n--- Yearly Statistics for %s ---\n', r_chg_labels_stats{r_stat_idx});
        fprintf('Year    |    Mean    |    Std     |  Median   |     N     |   Min    |   Max\n');
        fprintf('---------------------------------------------------------------------------\n');
        
        yearly_means = zeros(num_years, 1);
        yearly_stds = zeros(num_years, 1);
        yearly_medians = zeros(num_years, 1);
        yearly_counts = zeros(num_years, 1);
        
        for y_idx = 1:num_years
            year_val = unique_years{y_idx};
            year_mask = feature_table.Year == year_val;
            year_data = feature_table.(var_name)(year_mask);
            year_data = year_data(~isnan(year_data));
            
            if r_stat_idx == 1  % R_chg_1s만 yearly_data에 저장 (정규성 검정 등에 사용)
                yearly_data{y_idx} = year_data;
            end
            
            yearly_means(y_idx) = mean(year_data);
            yearly_stds(y_idx) = std(year_data);
            yearly_medians(y_idx) = median(year_data);
            yearly_counts(y_idx) = length(year_data);
            
            fprintf('%s  | %10.4f | %10.4f | %10.4f | %8d | %8.4f | %8.4f\n', ...
                char(year_val), yearly_means(y_idx), yearly_stds(y_idx), yearly_medians(y_idx), ...
                yearly_counts(y_idx), min(year_data), max(year_data));
        end
        fprintf('---------------------------------------------------------------------------\n');
    end
    
    % 정규성 검정 (Shapiro-Wilk 또는 Kolmogorov-Smirnov)
    fprintf('\n--- Normality Test (Kolmogorov-Smirnov) ---\n');
    normality_pvalues = zeros(num_years, 1);
    all_normal = true;
    for y_idx = 1:num_years
        year_val = unique_years{y_idx};
        year_data = yearly_data{y_idx};
        if length(year_data) >= 3
            [~, p_norm] = kstest((year_data - mean(year_data)) / std(year_data));
            normality_pvalues(y_idx) = p_norm;
            is_normal = p_norm > 0.05;
            if ~is_normal, all_normal = false; end
            norm_str = '(Normal)';
            if ~is_normal, norm_str = '(Non-normal)'; end
            fprintf('  %s: p = %.4f %s\n', char(year_val), p_norm, norm_str);
        end
    end
    
    % 연도별 차이 검정
    fprintf('\n--- Statistical Test for Yearly Differences ---\n');
    fprintf('두 가지 방법으로 비교합니다:\n');
    fprintf('  1. 전체 비교: 모든 연도를 한꺼번에 비교 (ANOVA/Kruskal-Wallis)\n');
    fprintf('  2. 쌍별 비교: 각 연도 쌍별로 비교 (t-test/Mann-Whitney)\n\n');
    
    if num_years >= 2
        % ========== 1. 전체 비교 (모든 연도를 한꺼번에) ==========
        fprintf('--- 1. Overall Comparison (All Years Together) ---\n');
        fprintf('질문: "모든 연도들 간에 통계적으로 유의한 차이가 있는가?"\n');
        
        % 모든 연도 데이터 준비
        all_values = [];
        all_groups = [];
        for y_idx = 1:num_years
            year_data = yearly_data{y_idx};
            all_values = [all_values; year_data];
            all_groups = [all_groups; repmat(y_idx, length(year_data), 1)];
        end
        
        % ANOVA (분산분석) 사용 - 모든 경우에 적용
        fprintf('Using One-way ANOVA (Analysis of Variance)\n');
        fprintf('  ANOVA는 모든 그룹의 평균이 같은지 검정합니다.\n');
        fprintf('  귀무가설(H0): 모든 연도의 평균이 같다\n');
        fprintf('  대립가설(H1): 적어도 하나의 연도는 다른 연도와 평균이 다르다\n');
        if num_years == 2
            fprintf('  (2개 그룹인 경우에도 ANOVA 사용 - t-test와 동일한 결과)\n');
        end
        [p_overall, ~, anova_stats] = anova1(all_values, all_groups, 'off');
        test_name_overall = 'One-way ANOVA';
        
        fprintf('  Test: %s\n', test_name_overall);
        fprintf('  p-value: %.6f\n', p_overall);
        
        if p_overall < 0.001
            fprintf('  Interpretation: *** Highly significant difference (p < 0.001)\n');
            fprintf('  Conclusion: 모든 연도들 간에 통계적으로 매우 유의한 차이가 있습니다.\n');
            fprintf('              적어도 하나의 연도는 다른 연도들과 다릅니다.\n');
        elseif p_overall < 0.01
            fprintf('  Interpretation: ** Very significant difference (p < 0.01)\n');
            fprintf('  Conclusion: 모든 연도들 간에 통계적으로 매우 유의한 차이가 있습니다.\n');
        elseif p_overall < 0.05
            fprintf('  Interpretation: * Significant difference (p < 0.05)\n');
            fprintf('  Conclusion: 모든 연도들 간에 통계적으로 유의한 차이가 있습니다.\n');
            fprintf('              적어도 하나의 연도는 다른 연도들과 다릅니다.\n');
        else
            fprintf('  Interpretation: Not significant (p >= 0.05)\n');
            fprintf('  Conclusion: 모든 연도들 간에 통계적으로 유의한 차이가 없습니다.\n');
            fprintf('              모든 연도의 평균(또는 중앙값)이 같다고 볼 수 있습니다.\n');
        end
        fprintf('\n');
        
        % ========== 2. 쌍별 비교 (각 연도 쌍별로) ==========
        fprintf('--- 2. Pairwise Comparisons (Each Year Pair) ---\n');
        fprintf('질문: "어떤 연도 쌍들 간에 통계적으로 유의한 차이가 있는가?"\n');
        fprintf('p-value 해석:\n');
        fprintf('  - p < 0.05: 두 연도 간에 통계적으로 유의한 차이가 있습니다 (우연이 아닙니다)\n');
        fprintf('  - p >= 0.05: 두 연도 간에 통계적으로 유의한 차이가 없습니다 (우연일 수 있습니다)\n');
        fprintf('  - p < 0.001: 매우 강한 유의성 (***)\n');
        fprintf('  - p < 0.01: 강한 유의성 (**)\n');
        fprintf('  - p < 0.05: 유의성 (*)\n\n');
        % 연도 쌍별 비교 수행
        num_comparisons = 0;
        comparison_results = {};
        
        for i = 1:num_years-1
            for j = i+1:num_years
                year1_val = unique_years{i};
                year2_val = unique_years{j};
                data1 = yearly_data{i};
                data2 = yearly_data{j};
                
                % 각 쌍에 대해 적절한 검정 선택
                if all_normal
                    % 정규 분포: t-test
                    [~, p_pair] = ttest2(data1, data2);
                    test_name = 'Two-sample t-test';
                    mean_diff = mean(data1) - mean(data2);
                else
                    % 비정규 분포: Mann-Whitney U test (Wilcoxon rank-sum test)
                    [p_pair, ~] = ranksum(data1, data2);
                    test_name = 'Mann-Whitney U test';
                    mean_diff = median(data1) - median(data2);
                end
                
                num_comparisons = num_comparisons + 1;
                comparison_results{num_comparisons, 1} = char(year1_val);
                comparison_results{num_comparisons, 2} = char(year2_val);
                comparison_results{num_comparisons, 3} = p_pair;
                comparison_results{num_comparisons, 4} = mean_diff;
                comparison_results{num_comparisons, 5} = test_name;
            end
        end
        
        % 결과 출력
        fprintf('--- Pairwise Comparison Results ---\n');
        fprintf('Year1 vs Year2 |    p-value    |   Mean Diff   | Significance | Test\n');
        fprintf('---------------------------------------------------------------------------\n');
        
        for comp_idx = 1:num_comparisons
            year1 = comparison_results{comp_idx, 1};
            year2 = comparison_results{comp_idx, 2};
            p_val = comparison_results{comp_idx, 3};
            mean_diff = comparison_results{comp_idx, 4};
            test_name = comparison_results{comp_idx, 5};
            
            % 유의성 표시
            if p_val < 0.001
                sig_str = '*** (p<0.001)';
            elseif p_val < 0.01
                sig_str = '** (p<0.01)';
            elseif p_val < 0.05
                sig_str = '* (p<0.05)';
            else
                sig_str = 'ns (p>=0.05)';
            end
            
            fprintf('%s vs %s | %13.6f | %13.4f | %-13s | %s\n', ...
                year1, year2, p_val, mean_diff, sig_str, test_name);
        end
        fprintf('---------------------------------------------------------------------------\n');
        
        % Bonferroni 보정 (다중 비교 보정)
        if num_comparisons > 1
            fprintf('\n--- Bonferroni Correction for Multiple Comparisons ---\n');
            fprintf('여러 쌍을 비교할 때는 우연히 유의한 결과가 나올 확률이 증가합니다.\n');
            fprintf('  예: 3개 연도를 비교하면 3쌍을 비교하게 되고,\n');
            fprintf('      각 비교에서 5%% 확률로 우연히 유의한 결과가 나올 수 있습니다.\n');
            fprintf('      따라서 전체적으로는 5%%보다 높은 확률로 오류가 발생할 수 있습니다.\n\n');
            fprintf('Bonferroni 보정: 유의수준을 비교 횟수로 나눕니다.\n');
            fprintf('  원래 유의수준: 0.05\n');
            fprintf('  비교 횟수: %d\n', num_comparisons);
            bonferroni_alpha = 0.05 / num_comparisons;
            fprintf('  보정된 유의수준: %.6f\n', bonferroni_alpha);
            fprintf('  (보정 후에도 유의하려면 p < %.6f 이어야 합니다)\n\n', bonferroni_alpha);
            
            fprintf('--- Bonferroni-Corrected Results ---\n');
            fprintf('Year1 vs Year2 |    p-value    | Bonferroni Corrected\n');
            fprintf('---------------------------------------------------------------------------\n');
            for comp_idx = 1:num_comparisons
                year1 = comparison_results{comp_idx, 1};
                year2 = comparison_results{comp_idx, 2};
                p_val = comparison_results{comp_idx, 3};
                
                if p_val < bonferroni_alpha
                    sig_str = sprintf('Significant (p<%.6f)', bonferroni_alpha);
                else
                    sig_str = 'Not significant';
                end
                
                fprintf('%s vs %s | %13.6f | %s\n', year1, year2, p_val, sig_str);
            end
            fprintf('---------------------------------------------------------------------------\n');
        end
        
        % 전체 비교와 쌍별 비교 결과 요약
        fprintf('\n--- Summary ---\n');
        fprintf('전체 비교 결과: ');
        if p_overall < 0.05
            fprintf('유의한 차이 있음 (p=%.6f) → 어떤 연도들 간에는 차이가 있습니다.\n', p_overall);
            fprintf('쌍별 비교 결과: 위의 표를 참고하여 어떤 연도 쌍들 간에 차이가 있는지 확인하세요.\n');
        else
            fprintf('유의한 차이 없음 (p=%.6f) → 모든 연도들이 비슷합니다.\n', p_overall);
            fprintf('쌍별 비교는 전체 비교가 유의할 때만 의미가 있습니다.\n');
        end
        
    else
        fprintf('  Not enough years for statistical comparison (need >= 2 years).\n');
    end
    
    fprintf('\n===========================================================================\n');

    % Figure 1: 그룹화된 이벤트별로 subplot 형태로 연도별로 시간-전류 시각화
    % (이벤트 내부 시간, 시작점을 0초로 맞춤, 연도별로 그룹화)
    
    % 그룹 키에서 전압과 전류 값 추출
    voltage_val = NaN;
    current_val = NaN;
    if contains(currentTargetGroupKey, 'V') && contains(currentTargetGroupKey, 'I')
        parts = strsplit(currentTargetGroupKey, '_');
        for p = 1:length(parts)
            if startsWith(parts{p}, 'V')
                voltage_val = str2double(parts{p}(2:end));
            elseif startsWith(parts{p}, 'I')
                current_val = str2double(parts{p}(2:end));
            end
        end
    end
    
    % 연도별로 이벤트 그룹화
    year_bins = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % SOC 범위 계산을 위한 변수
    all_soc_values = [];
    
    for seq_idx = 1:length(valid_sequences)
        if ~isempty(valid_sequences{seq_idx}) && isfield(valid_sequences{seq_idx}, 'year')
            seq = valid_sequences{seq_idx};
            if ~isnan(seq.year)
                year_key = sprintf('Y%d', seq.year);
                
                if ~isKey(year_bins, year_key)
                    year_bins(year_key) = [];
                end
                year_bins(year_key) = [year_bins(year_key), seq_idx];
                
                % SOC 값 수집
                if isfield(seq, 'start_soc') && ~isnan(seq.start_soc)
                    all_soc_values = [all_soc_values, seq.start_soc];
                end
            end
        end
    end
    
    % SOC 범위 문자열 생성 (소수점 3째자리까지)
    soc_range_str = '';
    if ~isempty(all_soc_values)
        min_soc = min(all_soc_values);
        max_soc = max(all_soc_values);
        soc_range_str = sprintf('SOC %.3f-%.3f%%', min_soc, max_soc);
    end
    
    % 연도 키를 정렬
    year_bin_keys = keys(year_bins);
    year_nums = zeros(length(year_bin_keys), 1);
    for k = 1:length(year_bin_keys)
        key = year_bin_keys{k};
        % "Y" 다음 숫자 추출
        num_str = regexp(key, 'Y(\d+)', 'tokens');
        if ~isempty(num_str)
            year_nums(k) = str2double(num_str{1}{1});
        end
    end
    [~, sort_idx] = sort(year_nums);
    sorted_year_bin_keys = year_bin_keys(sort_idx);
    num_year_bins = length(sorted_year_bin_keys);
    
    fig1 = figure('Name', ['Time-Current: ' currentTargetGroupKey], 'Position', [100, 100, 1600, 400*ceil(num_year_bins/2)]);
    plotTitle1 = sprintf('Time-Current Plot - %s - Group: %s', upper(targetEventType), strrep(currentTargetGroupKey, '_', ' '));
    sgtitle(plotTitle1, 'FontSize', 16, 'FontWeight', 'bold');
    
    % 색상 맵 생성 (HSV 색상 공간 사용)
    num_colors = 20;  % 최대 20개 색상
    color_map = hsv(num_colors);
    
    for year_idx = 1:num_year_bins
        year_bin_key = sorted_year_bin_keys{year_idx};
        seq_indices = year_bins(year_bin_key);
        
        if ~isempty(seq_indices)
            ax = subplot(ceil(num_year_bins/2), 2, year_idx);
            hold(ax, 'on');
            
            % 해당 연도의 모든 이벤트 시퀀스 플롯
            num_plotted = 0;
            for idx = 1:length(seq_indices)
                seq_idx = seq_indices(idx);
                if seq_idx <= length(valid_sequences) && ~isempty(valid_sequences{seq_idx})
                    seq = valid_sequences{seq_idx};
                    if isfield(seq, 'current')
                        % 전체 시간 표시 (필터링 제거)
                        time_filtered = seq.time_sec;
                        current_filtered = seq.current;
                        
                        % 이동평균 제거 - 원본 데이터 직접 사용
                        
                        % 자동 색상 할당 (이벤트마다 다른 색상)
                        color_idx = mod(num_plotted, num_colors) + 1;
                        plot(ax, time_filtered, current_filtered, '-', 'LineWidth', 0.5, 'Color', color_map(color_idx, :));
                        num_plotted = num_plotted + 1;
                    end
                end
            end
            
            % 연도 레이블 생성
            year_str = strrep(year_bin_key, 'Y', '');
            title(ax, sprintf('Year %s (n=%d events)', year_str, num_plotted));
            ylabel(ax, 'Current (A)');
            xlabel(ax, 'Time from Event Start (sec)');
            grid(ax, 'on');
            hold(ax, 'off');
        end
    end
    
    % Figure 1b: 그룹화된 이벤트별로 subplot 형태로 연도별로 시간-전압 시각화
    fig1b = figure('Name', ['Time-Voltage: ' currentTargetGroupKey], 'Position', [150, 150, 1600, 400*ceil(num_year_bins/2)]);
    plotTitle1b = sprintf('Time-Voltage Plot - %s - Group: %s', upper(targetEventType), strrep(currentTargetGroupKey, '_', ' '));
    sgtitle(plotTitle1b, 'FontSize', 16, 'FontWeight', 'bold');
    
    for year_idx = 1:num_year_bins
        year_bin_key = sorted_year_bin_keys{year_idx};
        seq_indices = year_bins(year_bin_key);
        
        if ~isempty(seq_indices)
            ax = subplot(ceil(num_year_bins/2), 2, year_idx);
            hold(ax, 'on');
            
            % 해당 연도의 모든 이벤트 시퀀스 플롯
            num_plotted = 0;
            for idx = 1:length(seq_indices)
                seq_idx = seq_indices(idx);
                if seq_idx <= length(valid_sequences) && ~isempty(valid_sequences{seq_idx})
                    seq = valid_sequences{seq_idx};
                    if isfield(seq, 'voltage')
                        % 전체 시간 표시 (필터링 제거)
                        time_filtered = seq.time_sec;
                        voltage_filtered = seq.voltage;
                        
                        % 이동평균 제거 - 원본 데이터 직접 사용
                        
                        % 자동 색상 할당 (이벤트마다 다른 색상)
                        color_idx = mod(num_plotted, num_colors) + 1;
                        plot(ax, time_filtered, voltage_filtered, '-', 'LineWidth', 0.5, 'Color', color_map(color_idx, :));
                        num_plotted = num_plotted + 1;
                    end
                end
            end
            
            % 연도 레이블 생성
            year_str = strrep(year_bin_key, 'Y', '');
            title(ax, sprintf('Year %s (n=%d events)', year_str, num_plotted));
            ylabel(ax, 'Voltage (V)');
            xlabel(ax, 'Time from Event Start (sec)');
            grid(ax, 'on');
            hold(ax, 'off');
        end
    end
    
    % Figure 2: 그룹화된 이벤트를 subplot으로 여러 Rchg 박스플롯 생성
    fig2 = figure('Name', ['Boxplot: ' currentTargetGroupKey], 'Position', [200, 200, 1600, 1200]);
    plotTitle2 = sprintf('Rchg - %s - Group: %s', upper(targetEventType), strrep(currentTargetGroupKey, '_', ' '));
    sgtitle(plotTitle2, 'FontSize', 16, 'FontWeight', 'bold');
    
    unique_years = categories(feature_table.Year);
    sorted_year_labels = sort(unique_years);
    num_years = length(unique_years);
    r_chg_vars = {'R_chg_1s', 'R_chg_5s', 'R_chg_10s', 'R_chg_30s'};
    r_chg_labels = {'R_{chg, 1s}', 'R_{chg, 5s}', 'R_{chg, 10s}', 'R_{chg, 30s}'};
    
    for r_idx = 1:length(r_chg_vars)
        ax = subplot(2, 2, r_idx);
        var_name = r_chg_vars{r_idx};
        if ismember(var_name, feature_table.Properties.VariableNames)
            % 박스플롯 생성
            boxplot(feature_table.(var_name), feature_table.Year, 'GroupOrder', cellstr(sorted_year_labels));
            
            % 연도별 평균과 표준편차 계산
            yearly_means = zeros(num_years, 1);
            yearly_stds = zeros(num_years, 1);
            yearly_data = cell(num_years, 1);
            
            for y_idx = 1:num_years
                year_val = sorted_year_labels{y_idx};
                year_mask = feature_table.Year == year_val;
                year_data = feature_table.(var_name)(year_mask);
                year_data = year_data(~isnan(year_data));
                
                yearly_data{y_idx} = year_data;
                yearly_means(y_idx) = mean(year_data);
                yearly_stds(y_idx) = std(year_data);
            end
            
            % 제목 설정 (p-value 제거)
            title(r_chg_labels{r_idx});
            
            ylabel('Resistance (m\Omega)');
            xlabel('Year');
            grid on;
            
            % 평균 ± std를 박스 옆에 텍스트로 표시
            hold(ax, 'on');
            x_positions = 1:num_years;
            for y_idx = 1:num_years
                x_pos = x_positions(y_idx);
                mean_val = yearly_means(y_idx);
                std_val = yearly_stds(y_idx);
                event_count = length(yearly_data{y_idx});
                
                % 평균 ± std 텍스트 위치: 박스플롯의 오른쪽에 표시 (소수점 4째자리까지)
                text(x_pos + 0.35, mean_val, sprintf('%.4f±%.4f', mean_val, std_val), ...
                    'FontSize', 8, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
                
                % 각 연도별 이벤트 개수를 박스 위에 크게 표시
                y_limits = ylim(ax);
                text_y_pos = y_limits(2) * 0.95;  % 상단에서 5% 아래 위치
                text(x_pos, text_y_pos, sprintf('n=%d', event_count), ...
                    'FontSize', 14, 'FontWeight', 'bold', ...
                    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
            end
            hold(ax, 'off');
        end
    end
    
    % 그래프 저장
    if savePlots
        plotFileName1 = sprintf('%s_%s_%s_TimeCurrent.fig', targetRack, targetEventType, currentTargetGroupKey);
        plotFileName1b = sprintf('%s_%s_%s_TimeVoltage.fig', targetRack, targetEventType, currentTargetGroupKey);
        plotFileName2 = sprintf('%s_%s_%s_Boxplot.fig', targetRack, targetEventType, currentTargetGroupKey);
        savefig(fig1, fullfile(group_folder, plotFileName1));
        savefig(fig1b, fullfile(group_folder, plotFileName1b));
        savefig(fig2, fullfile(group_folder, plotFileName2));
        fprintf('[FIG SAVED] Group %d/%d: %s -> %s\n', g, length(top_group_keys), currentTargetGroupKey, group_folder);
        fprintf('  - %s\n  - %s\n  - %s\n', plotFileName1, plotFileName1b, plotFileName2);
        
        % 모든 figure 닫기
        close(fig1);
        close(fig1b);
        close(fig2);
    else
        fprintf('[FIG NOT SAVED] Group %d/%d: %s (savePlots=false)\n', g, length(top_group_keys), currentTargetGroupKey);
    end
    
    % 그룹별 데이터 저장 (연도별, 이벤트별)
    fprintf('\nSaving group data to results structure...\n');
    unique_years_table = categories(feature_table.Year);
    
    for y_idx = 1:length(unique_years_table)
        year_val = unique_years_table{y_idx};
        year_key = char(year_val);
        % Y2021 형식으로 변환 (숫자만 있으면 Y 추가)
        if ~startsWith(year_key, 'Y')
            year_key_safe = ['Y' year_key];
        else
            year_key_safe = year_key;
        end
        
        % 연도별 구조체 초기화
        if ~isfield(results_data.(group_key_safe), year_key_safe)
            results_data.(group_key_safe).(year_key_safe) = struct();
        end
        
        % 해당 연도의 이벤트 인덱스 찾기
        year_mask = feature_table.Year == year_val;
        year_table_indices = find(year_mask);
        
        event_counter = 0;
        for table_idx = 1:length(year_table_indices)
            table_row = year_table_indices(table_idx);
            seq_idx = table_row;  % valid_sequences와 feature_table은 같은 순서
            
            if seq_idx <= length(valid_sequences) && ~isempty(valid_sequences{seq_idx})
                seq = valid_sequences{seq_idx};
                if isfield(seq, 'voltage') && isfield(seq, 'current')
                    event_counter = event_counter + 1;
                    event_key = sprintf('Event%d', event_counter);
                    
                    % 원본 이벤트에서 SOC 시계열 가져오기
                    orig_evt_idx = original_indices(table_row);
                    if orig_evt_idx <= length(all_target_events)
                        orig_evt = all_target_events{orig_evt_idx};
                        if isfield(orig_evt, 'soc_seq_pct') && ~isempty(orig_evt.soc_seq_pct)
                            soc_seq = orig_evt.soc_seq_pct;
                        else
                            soc_seq = [];
                        end
                    else
                        soc_seq = [];
                    end
                    
                    % 전력 시계열 계산 (V * I)
                    power_seq = seq.voltage .* seq.current;
                    
                    % 표준편차 계산
                    max_p_std = NaN;
                    max_I_std = NaN;
                    if ~isempty(power_seq) && length(power_seq) > 1
                        max_p_std = std(power_seq);
                    end
                    if isfield(seq, 'current') && length(seq.current) > 1
                        max_I_std = std(seq.current);
                    end
                    
                    % 이벤트 데이터 저장
                    results_data.(group_key_safe).(year_key_safe).(event_key) = struct();
                    results_data.(group_key_safe).(year_key_safe).(event_key).voltage_seq = seq.voltage;
                    results_data.(group_key_safe).(year_key_safe).(event_key).current_seq = seq.current;
                    results_data.(group_key_safe).(year_key_safe).(event_key).soc_seq = soc_seq;
                    results_data.(group_key_safe).(year_key_safe).(event_key).time_seq = seq.time_sec;
                    results_data.(group_key_safe).(year_key_safe).(event_key).timestamp = feature_table.timestamp(table_row);
                    results_data.(group_key_safe).(year_key_safe).(event_key).R_chg_1s = feature_table.R_chg_1s(table_row);
                    results_data.(group_key_safe).(year_key_safe).(event_key).R_chg_5s = feature_table.R_chg_5s(table_row);
                    results_data.(group_key_safe).(year_key_safe).(event_key).R_chg_10s = feature_table.R_chg_10s(table_row);
                    results_data.(group_key_safe).(year_key_safe).(event_key).R_chg_30s = feature_table.R_chg_30s(table_row);
                    results_data.(group_key_safe).(year_key_safe).(event_key).max_p_std = max_p_std;
                    results_data.(group_key_safe).(year_key_safe).(event_key).max_I_std = max_I_std;
                end
            end
        end
        fprintf('  %s: %d events saved\n', year_key_safe, event_counter);
    end
    
    fprintf('Completed processing group %d/%d: %s\n', g, length(top_group_keys), currentTargetGroupKey);
end

fprintf('\n=== Finished processing ALL %d groups ===\n', length(top_group_keys));
fprintf('Summary: %d groups processed, %d groups should have fig files saved\n', length(top_group_keys), length(top_group_keys));

% 전체 결과를 mat 파일로 저장 (EventType 폴더 안에)
fprintf('\nSaving all results to mat file...\n');
resultsFileName = sprintf('%s_%s_resistance_results.mat', targetRack, targetEventType);
resultsFilePath = fullfile(eventTypeDir, resultsFileName);
save(resultsFilePath, 'results_data', '-v7.3');
fprintf('Results saved to: %s\n', resultsFilePath);

fprintf('\nAll processing finished.\n');

%% Helper Function (v4.3)
function [R_chg_1s, R_chg_5s, R_chg_10s, R_chg_30s] = calculate_dynamic_resistance_v4(event_struct)
    % EventDetection.m에서 이미 SG 필터(Savitzky-Golay, 윈도우 5, 차수 1)로 필터링된 데이터를 사용
    % SG 필터는 엣지(급격한 변화)를 유지하면서 노이즈를 제거하므로,
    % t=0 시점의 정확한 타이밍이 보존되어 R_chg_1s 계산이 정확해짐
    R_chg_1s = NaN; R_chg_5s = NaN; R_chg_10s = NaN; R_chg_30s = NaN;
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
        % 숫자 형식인 경우 (초 단위로 가정)
        t_rel = time_seq - time_seq(1);
    end
    
    % 1초 후 저항 계산: 시작점(0초)과 1초 후 지점 비교
    idx_1s = find(t_rel >= 1.0, 1);
    if ~isempty(idx_1s) && idx_1s > 1 && idx_1s <= length(V) && idx_1s <= length(I)
        delta_V_1s = V(idx_1s) - V(1);
        delta_I_1s = I(idx_1s) - I(1);
        if abs(delta_I_1s) > 1e-6  % 전류 변화가 0이 아닌 경우만 계산
            R_chg_1s = (delta_V_1s / delta_I_1s) * 1000;
        end
    end
    
    % 5초 후 저항 계산: 시작점(0초)과 5초 후 지점 비교
    idx_5s = find(t_rel >= 5.0, 1);
    if ~isempty(idx_5s) && idx_5s > 1 && idx_5s <= length(V) && idx_5s <= length(I)
        delta_V_5s = V(idx_5s) - V(1);
        delta_I_5s = I(idx_5s) - I(1);
        if abs(delta_I_5s) > 1e-6  % 전류 변화가 0이 아닌 경우만 계산
            R_chg_5s = (delta_V_5s / delta_I_5s) * 1000;
        end
    end
    
    % 10초 후 저항 계산: 시작점(0초)과 10초 후 지점 비교
    idx_10s = find(t_rel >= 10.0, 1);
    if ~isempty(idx_10s) && idx_10s > 1 && idx_10s <= length(V) && idx_10s <= length(I)
        delta_V_10s = V(idx_10s) - V(1);
        delta_I_10s = I(idx_10s) - I(1);
        if abs(delta_I_10s) > 1e-6  % 전류 변화가 0이 아닌 경우만 계산
            R_chg_10s = (delta_V_10s / delta_I_10s) * 1000;
        end
    end
    
    % 30초 후 저항 계산: 시작점(0초)과 30초 후 지점 비교
    idx_30s = find(t_rel >= 30.0, 1);
    if ~isempty(idx_30s) && idx_30s > 1 && idx_30s <= length(V) && idx_30s <= length(I)
        delta_V_30s = V(idx_30s) - V(1);
        delta_I_30s = I(idx_30s) - I(1);
        if abs(delta_I_30s) > 1e-6  % 전류 변화가 0이 아닌 경우만 계산
            R_chg_30s = (delta_V_30s / delta_I_30s) * 1000;
        end
    end
    
    if ~isnan(R_chg_1s) && R_chg_1s < 0, R_chg_1s = NaN; end
    if ~isnan(R_chg_5s) && R_chg_5s < 0, R_chg_5s = NaN; end
    if ~isnan(R_chg_10s) && R_chg_10s < 0, R_chg_10s = NaN; end
    if ~isnan(R_chg_30s) && R_chg_30s < 0, R_chg_30s = NaN; end
end