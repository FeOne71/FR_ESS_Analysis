%% 다년도 데이터 분석 (2021, 2022, 2023년 6월)
% 목적: 연도별 저항 변화율 계산 → 열화 라벨 생성

clear; clc; close all;

%% 경로 및 파라미터 설정
baseDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old';
savePath = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_0113';

% 분석 연도
years = [2021, 2022, 2023];

% 파라미터 (v3 조건)
Cnom_cell = 64;
Cnom_rack = 128;
idle_threshold = 0.01 * Cnom_cell;
min_duration = 60;
max_I_std = 1.5;
std_check_window = [3, 60];
R_timepoints = [1, 3, 5, 10, 30, 60];

fprintf('========================================\n');
fprintf('  다년도 데이터 분석 시작\n');
fprintf('========================================\n');

fprintf('분석 연도: %s\n', mat2str(years));
fprintf('조건: Duration≥60s, Std(3-60s)≤1.5A\n');

%% 각 연도별 분석
all_results = struct();

for y_idx = 1:length(years)
    year = years(y_idx);
    fprintf('\n=== %d년 6월 분석 ===\n', year);
    
    dataPath = fullfile(baseDataPath, sprintf('%d', year), sprintf('%d06', year));
    
    % 폴더 존재 확인
    if ~exist(dataPath, 'dir')
        fprintf('  경고: %s 폴더가 없습니다!\n', dataPath);
        continue;
    end
    
    % 데이터 로드
    fprintf('  데이터 로딩...\n');
    matFiles = dir(fullfile(dataPath, '*.mat'));
    
    if isempty(matFiles)
        fprintf('  경고: mat 파일이 없습니다!\n');
        continue;
    end
    
    allData = struct();
    allData.voltage = [];
    allData.current = [];
    allData.soc = [];
    allData.temperature = [];
    
    for day = 1:length(matFiles)
        filePath = fullfile(dataPath, matFiles(day).name);
        data = load(filePath);
        rack01 = data.Raw.Rack01;
        
        allData.voltage = [allData.voltage; rack01.AverageCV_V];
        allData.current = [allData.current; rack01.DCCurrent_A];
        allData.soc = [allData.soc; rack01.SOCPct];
        allData.temperature = [allData.temperature; rack01.AverageMT_degC];
    end
    
    fprintf('    데이터: %d개 (%.1f일)\n', length(allData.current), length(allData.current)/86400);
    
    % 이벤트 검출
    fprintf('  이벤트 검출...\n');
    
    is_idle = abs(allData.current) < idle_threshold;
    idle_to_active = find(is_idle(1:end-1) & ~is_idle(2:end));
    
    events = struct();
    event_count = 0;
    filtered = struct();
    filtered.too_short = 0;
    filtered.high_std = 0;
    filtered.insufficient = 0;
    
    for k = 1:length(idle_to_active)
        idx_idle = idle_to_active(k);
        idx_start = idx_idle + 1;
        
        if idx_start > length(allData.current), continue; end
        
        event_type = sign(allData.current(idx_start));
        if event_type == 0, continue; end
        
        % 종료점
        idx_end = idx_start;
        for idx = idx_start+1:length(allData.current)
            I_current = allData.current(idx);
            if abs(I_current) < idle_threshold || ...
               sign(allData.current(idx_start)) * allData.current(idx) < 0
                idx_end = idx - 1;
                break;
            end
            idx_end = idx;
        end
        
        duration_sec = idx_end - idx_start + 1;
        if duration_sec < min_duration
            filtered.too_short = filtered.too_short + 1;
            continue;
        end
        
        % 전류 std 체크
        check_idx_start = idx_start + std_check_window(1);
        check_idx_end = idx_start + std_check_window(2);
        if check_idx_end > idx_end
            filtered.insufficient = filtered.insufficient + 1;
            continue;
        end
        
        I_check = allData.current(check_idx_start:check_idx_end);
        I_std_3_60 = std(I_check);
        if I_std_3_60 > max_I_std
            filtered.high_std = filtered.high_std + 1;
            continue;
        end
        
        I_avg_3_60 = mean(abs(I_check));
        
        % 이벤트 저장
        event_count = event_count + 1;
        events(event_count).idx_idle = idx_idle;
        events(event_count).idx_start = idx_start;
        events(event_count).idx_end = idx_end;
        events(event_count).duration = duration_sec;
        events(event_count).type = event_type;
        events(event_count).V_idle = allData.voltage(idx_idle);
        events(event_count).I_idle = allData.current(idx_idle);
        
        % 저항 계산
        for t_idx = 1:length(R_timepoints)
            t_sec = R_timepoints(t_idx);
            idx_t = idx_start + t_sec;
            if idx_t <= idx_end
                V_t = allData.voltage(idx_t);
                I_t = allData.current(idx_t);
                dV = V_t - events(event_count).V_idle;
                dI = I_t - events(event_count).I_idle;
                if abs(dI) > 0.1
                    R_t = abs(dV / dI) * 1000;
                    events(event_count).(sprintf('R_%ds', t_sec)) = R_t;
                else
                    events(event_count).(sprintf('R_%ds', t_sec)) = NaN;
                end
            else
                events(event_count).(sprintf('R_%ds', t_sec)) = NaN;
            end
        end
        
        % 특성 저장
        events(event_count).I_avg_3_60s = I_avg_3_60;
        events(event_count).I_std_3_60s = I_std_3_60;
        events(event_count).SOC_mean = mean(allData.soc(idx_start:idx_end));
        events(event_count).T_mean = mean(allData.temperature(idx_start:idx_end));
        events(event_count).C_rate = (I_avg_3_60 / 2) / Cnom_cell;
    end
    
    fprintf('    이벤트: %d개 (필터: %d+%d+%d)\n', ...
            event_count, filtered.too_short, filtered.insufficient, filtered.high_std);
    
    % 충전/방전 분리
    charge_events = events([events.type] > 0);
    discharge_events = events([events.type] < 0);
    
    fprintf('    충전: %d개, 방전: %d개\n', length(charge_events), length(discharge_events));
    
    % 통계 출력
    if ~isempty(charge_events)
        fprintf('    충전 R_10s: %.3f ± %.3f mOhm\n', ...
                mean([charge_events.R_10s], 'omitnan'), ...
                std([charge_events.R_10s], 'omitnan'));
    end
    
    % 결과 저장
    year_result = struct();
    year_result.year = year;
    year_result.events = events;
    year_result.charge_events = charge_events;
    year_result.discharge_events = discharge_events;
    year_result.filtered_stats = filtered;
    
    all_results.(sprintf('year%d', year)) = year_result;
end

%% 연도별 비교
fprintf('\n========================================\n');
fprintf('         연도별 비교\n');
fprintf('========================================\n');

% 메인 그룹 (SOC 55-65%, 25-30A) 추출
soc_min = 55; soc_max = 65;
curr_min = 25; curr_max = 30;

comparison = struct();

for y_idx = 1:length(years)
    year = years(y_idx);
    year_field = sprintf('year%d', year);
    
    if ~isfield(all_results, year_field)
        fprintf('%d년: 데이터 없음\n', year);
        continue;
    end
    
    charge_events = all_results.(year_field).charge_events;
    
    % 메인 그룹 필터링
    main_group = [];
    for i = 1:length(charge_events)
        evt = charge_events(i);
        if evt.SOC_mean >= soc_min && evt.SOC_mean <= soc_max && ...
           evt.I_avg_3_60s >= curr_min && evt.I_avg_3_60s < curr_max
            main_group = [main_group; evt];
        end
    end
    
    if isempty(main_group)
        fprintf('%d년: 메인 그룹 없음\n', year);
        continue;
    end
    
    fprintf('%d년 (N=%d):\n', year, length(main_group));
    
    % 시간별 저항 통계
    for t_idx = 1:length(R_timepoints)
        t_sec = R_timepoints(t_idx);
        R_field = sprintf('R_%ds', t_sec);
        R_values = [main_group.(R_field)];
        R_mean = mean(R_values, 'omitnan');
        R_std = std(R_values, 'omitnan');
        
        comparison.(year_field).(R_field).mean = R_mean;
        comparison.(year_field).(R_field).std = R_std;
        comparison.(year_field).(R_field).n = sum(~isnan(R_values));
        
        fprintf('  R_%2ds: %.3f ± %.3f mOhm\n', t_sec, R_mean, R_std);
    end
    fprintf('\n');
end

%% 열화율 계산 (2021 대비)
fprintf('========================================\n');
fprintf('         열화율 (2021년 대비)\n');
fprintf('========================================\n');

if isfield(comparison, 'year2021')
    base_year = 'year2021';
    
    for y_idx = 2:length(years)
        year = years(y_idx);
        year_field = sprintf('year%d', year);
        
        if ~isfield(comparison, year_field)
            continue;
        end
        
        fprintf('%d년:\n', year);
        
        for t_idx = 1:length(R_timepoints)
            t_sec = R_timepoints(t_idx);
            R_field = sprintf('R_%ds', t_sec);
            
            R_base = comparison.(base_year).(R_field).mean;
            R_current = comparison.(year_field).(R_field).mean;
            
            degradation_rate = (R_current - R_base) / R_base * 100;
            
            fprintf('  R_%2ds: %.3f → %.3f mOhm (+%.2f%%)\n', ...
                    t_sec, R_base, R_current, degradation_rate);
        end
        fprintf('\n');
    end
end

%% 저장
results_multi_year = struct();
results_multi_year.all_results = all_results;
results_multi_year.comparison = comparison;
results_multi_year.parameters.years = years;
results_multi_year.parameters.main_group_criteria.soc = [soc_min, soc_max];
results_multi_year.parameters.main_group_criteria.current = [curr_min, curr_max];

save(fullfile(savePath, 'MultiYear_Analysis_2021_2023.mat'), 'results_multi_year', '-v7.3');

fprintf('========================================\n');
fprintf('저장: MultiYear_Analysis_2021_2023.mat\n');
fprintf('========================================\n');
