%% 필드 데이터 이벤트 추출 및 시간별 저항 계산
% 목적: 열화 모드 정량화를 위한 물리 기반 지표 추출
% 최종 목표: 물리 모델 + AI 하이브리드 모델 개발

clear; clc; close all;

%% 경로 설정
dataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old\2021\202106';
savePath = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_0113';

fprintf('========================================\\n');
fprintf('  이벤트 기반 저항 추출 (수정 버전)\\n');
fprintf('========================================\\n\\n');

%% 파라미터 설정
Cnom_cell = 64;        % Cell 용량 [Ah]
Cnom_rack = 128;       % Rack 용량 [Ah]

% 이벤트 검출 조건
idle_threshold = 0.01 * Cnom_rack;  % 1.28A
min_duration = 30;                   % 최소 30초
max_I_std = 1.5;                     % 3-30초 구간 전류 std

% 저항 계산 시점 [초]
R_timepoints = [1, 3, 5, 10, 30, 60];

% 전류 그룹 정의 [A] (Rack 기준)
current_groups = struct();
current_groups.ranges = [20 30; 30 50; 50 100; 100 200; 200 350];
current_groups.names = {'20-30A', '30-50A', '50-100A', '100-200A', '200-350A'};

fprintf('파라미터:\\n');
fprintf('  Idle threshold: %.2f A\\n', idle_threshold);
fprintf('  최소 지속시간: %d 초\\n', min_duration);
fprintf('  전류 std 허용: %.2f A\\n', max_I_std);
fprintf('  저항 계산 시점: %s\\n', mat2str(R_timepoints));
fprintf('  전류 그룹: %d개\\n\\n', size(current_groups.ranges, 1));

%% 데이터 로드 (30일 전체)
fprintf('[1단계] 데이터 로딩 중...\\n');

matFiles = dir(fullfile(dataPath, '*.mat'));
numDays = length(matFiles);
fprintf('  총 파일 수: %d개\\n', numDays);

% 전체 데이터 통합
allData = struct();
allData.time = [];
allData.voltage = [];
allData.current = [];
allData.soc = [];
allData.temperature = [];
allData.dayIndex = [];

for day = 1:numDays
    if mod(day, 5) == 0
        fprintf('  진행: %d/%d\\n', day, numDays);
    end
    
    filePath = fullfile(dataPath, matFiles(day).name);
    data = load(filePath);
    rack01 = data.Raw.Rack01;
    
    allData.time = [allData.time; rack01.Time];
    allData.voltage = [allData.voltage; rack01.AverageCV_V];
    allData.current = [allData.current; rack01.DCCurrent_A];
    allData.soc = [allData.soc; rack01.SOCPct];
    allData.temperature = [allData.temperature; rack01.AverageMT_degC];
    allData.dayIndex = [allData.dayIndex; ones(length(rack01.SOCPct), 1) * day];
end

fprintf('\\n데이터 로드 완료!\\n');
fprintf('  총 데이터 포인트: %d개\\n', length(allData.current));
fprintf('  기간: %.2f일\\n\\n', length(allData.current) / 86400);

%% 이벤트 검출 (올바른 정의)
fprintf('[2단계] 이벤트 검출 중...\\n');

% Idle 판별
is_idle = abs(allData.current) < idle_threshold;

% Idle -> Active 전환점 찾기
idle_to_active = find(is_idle(1:end-1) & ~is_idle(2:end));
fprintf('  Idle->Active 전환: %d개\\n', length(idle_to_active));

% 이벤트 구조체
events = struct();
event_count = 0;

filtered_stats = struct();
filtered_stats.too_short = 0;
filtered_stats.high_std = 0;
filtered_stats.insufficient_length = 0;
filtered_stats.sign_change = 0;

for k = 1:length(idle_to_active)
    % Idle 마지막 인덱스
    idx_idle = idle_to_active(k);
    idx_start = idx_idle + 1;
    
    if idx_start > length(allData.current)
        continue;
    end
    
    % 이벤트 타입 판별
    event_type = sign(allData.current(idx_start));
    if event_type == 0
        continue;
    end
    
    % 종료 조건 찾기
    idx_end = idx_start;
    for idx = idx_start+1:length(allData.current)
        I_current = allData.current(idx);
        
        % 조건 1: Idle 범위 내로 돌아옴
        if abs(I_current) < idle_threshold
            idx_end = idx - 1;
            break;
        end
        
        % 조건 2: 전류 방향 반전 (충전->방전 or 방전->충전)
        if sign(allData.current(idx_start)) * allData.current(idx) < 0
            idx_end = idx - 1;
            break;
        end
        
        idx_end = idx;
    end
    
    % 지속시간 체크
    duration_sec = idx_end - idx_start + 1;
    if duration_sec < min_duration
        filtered_stats.too_short = filtered_stats.too_short + 1;
        continue;
    end
    
    % 60초 이상 데이터 필요 (저항 계산용)
    if duration_sec < 60
        filtered_stats.insufficient_length = filtered_stats.insufficient_length + 1;
        continue;
    end
    
    % 3-30초 구간 전류 표준편차 체크
    check_idx_start = idx_start + 3;
    check_idx_end = idx_start + 30;
    
    if check_idx_end > idx_end
        filtered_stats.insufficient_length = filtered_stats.insufficient_length + 1;
        continue;
    end
    
    I_check = allData.current(check_idx_start:check_idx_end);
    I_std_3_30 = std(I_check);
    
    if I_std_3_30 > max_I_std
        filtered_stats.high_std = filtered_stats.high_std + 1;
        continue;
    end
    
    % 3-30초 평균 전류 (전류 그룹 분류용)
    I_avg_3_30 = mean(abs(I_check));
    
    % 이벤트 저장
    event_count = event_count + 1;
    
    events(event_count).idx_idle = idx_idle;
    events(event_count).idx_start = idx_start;
    events(event_count).idx_end = idx_end;
    events(event_count).duration = duration_sec;
    events(event_count).type = event_type;
    
    % Idle 상태 (기준점)
    events(event_count).V_idle = allData.voltage(idx_idle);
    events(event_count).I_idle = allData.current(idx_idle);
    
    % 저항 계산 (1, 3, 5, 10, 30, 60초)
    for t_idx = 1:length(R_timepoints)
        t_sec = R_timepoints(t_idx);
        idx_t = idx_start + t_sec;
        
        if idx_t <= idx_end
            V_t = allData.voltage(idx_t);
            I_t = allData.current(idx_t);
            
            dV = V_t - events(event_count).V_idle;
            dI = I_t - events(event_count).I_idle;
            
            if abs(dI) > 0.1
                R_t = abs(dV / dI) * 1000;  % mOhm
                events(event_count).(sprintf('R_%ds', t_sec)) = R_t;
            else
                events(event_count).(sprintf('R_%ds', t_sec)) = NaN;
            end
        else
            events(event_count).(sprintf('R_%ds', t_sec)) = NaN;
        end
    end
    
    % 전류 특성
    events(event_count).I_avg_3_30s = I_avg_3_30;
    events(event_count).I_std_3_30s = I_std_3_30;
    events(event_count).I_max = max(abs(allData.current(idx_start:idx_end)));
    
    % 전압 특성
    events(event_count).V_mean = mean(allData.voltage(idx_start:idx_end));
    events(event_count).V_start = allData.voltage(idx_start);
    events(event_count).V_end = allData.voltage(idx_end);
    events(event_count).V_drop = events(event_count).V_start - events(event_count).V_end;
    
    % SOC 특성
    events(event_count).SOC_start = allData.soc(idx_start);
    events(event_count).SOC_mean = mean(allData.soc(idx_start:idx_end));
    events(event_count).SOC_end = allData.soc(idx_end);
    
    % 온도 특성
    events(event_count).T_start = allData.temperature(idx_start);
    events(event_count).T_end = allData.temperature(idx_end);
    events(event_count).T_rise = events(event_count).T_end - events(event_count).T_start;
    events(event_count).T_rise_rate = events(event_count).T_rise / (duration_sec / 60);  % °C/min
    
    % 전류 그룹 분류
    current_group_idx = 0;
    for g = 1:size(current_groups.ranges, 1)
        if I_avg_3_30 >= current_groups.ranges(g, 1) && ...
           I_avg_3_30 < current_groups.ranges(g, 2)
            current_group_idx = g;
            break;
        end
    end
    events(event_count).current_group = current_group_idx;
    
    % Cell 레벨 C-rate (2P 구조)
    I_cell = I_avg_3_30 / 2;
    events(event_count).C_rate = I_cell / Cnom_cell;
end

fprintf('\\n이벤트 검출 완료!\\n');
fprintf('  총 이벤트: %d개\\n', event_count);
fprintf('  필터링 (30초 미만): %d개\\n', filtered_stats.too_short);
fprintf('  필터링 (60초 미만): %d개\\n', filtered_stats.insufficient_length);
fprintf('  필터링 (전류 std 초과): %d개\\n\\n', filtered_stats.high_std);

%% 이벤트 분석
fprintf('[3단계] 이벤트 분석\\n\\n');

if event_count == 0
    fprintf('이벤트가 없습니다!\\n');
    return;
end

% 충전/방전 분리
charge_events = events([events.type] > 0);
discharge_events = events([events.type] < 0);

fprintf('=== 이벤트 타입별 통계 ===\\n');
fprintf('  충전 이벤트: %d개\\n', length(charge_events));
fprintf('  방전 이벤트: %d개\\n\\n', length(discharge_events));

% 전류 그룹별 통계
fprintf('=== 전류 그룹별 분포 ===\\n');
for g = 1:length(current_groups.names)
    n_events = sum([events.current_group] == g);
    fprintf('  %s: %d개\\n', current_groups.names{g}, n_events);
end
fprintf('  기타: %d개\\n\\n', sum([events.current_group] == 0));

% 충전 이벤트 통계
if ~isempty(charge_events)
    fprintf('=== 충전 이벤트 통계 ===\\n');
    fprintf('  평균 지속시간: %.1f초\\n', mean([charge_events.duration]));
    fprintf('  평균 전류(3-30s): %.2f A\\n', mean([charge_events.I_avg_3_30s]));
    fprintf('  평균 SOC: %.2f%%\\n', mean([charge_events.SOC_mean]));
    fprintf('  평균 C-rate: %.4f C\\n', mean([charge_events.C_rate]));
    fprintf('\\n  시간별 평균 저항:\\n');
    for t_idx = 1:length(R_timepoints)
        t_sec = R_timepoints(t_idx);
        R_field = sprintf('R_%ds', t_sec);
        R_values = [charge_events.(R_field)];
        fprintf('    R_%ds: %.2f ± %.2f mOhm (N=%d)\\n', ...
                t_sec, mean(R_values, 'omitnan'), std(R_values, 'omitnan'), ...
                sum(~isnan(R_values)));
    end
    fprintf('\\n');
end

% 방전 이벤트 통계
if ~isempty(discharge_events)
    fprintf('=== 방전 이벤트 통계 ===\\n');
    fprintf('  평균 지속시간: %.1f초\\n', mean([discharge_events.duration]));
    fprintf('  평균 전류(3-30s): %.2f A\\n', mean([discharge_events.I_avg_3_30s]));
    fprintf('  평균 SOC: %.2f%%\\n', mean([discharge_events.SOC_mean]));
    fprintf('  평균 C-rate: %.4f C\\n', mean([discharge_events.C_rate]));
    fprintf('\\n  시간별 평균 저항:\\n');
    for t_idx = 1:length(R_timepoints)
        t_sec = R_timepoints(t_idx);
        R_field = sprintf('R_%ds', t_sec);
        R_values = [discharge_events.(R_field)];
        fprintf('    R_%ds: %.2f ± %.2f mOhm (N=%d)\\n', ...
                t_sec, mean(R_values, 'omitnan'), std(R_values, 'omitnan'), ...
                sum(~isnan(R_values)));
    end
    fprintf('\\n');
end

%% 시각화
fprintf('[4단계] 시각화 생성 중...\\n');

figDir = fullfile(savePath, 'figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

% 1. 시간별 저항 변화
fig1 = figure('Position', [100, 100, 1200, 600]);

subplot(1, 2, 1);
if ~isempty(charge_events)
    R_matrix_charge = zeros(length(charge_events), length(R_timepoints));
    for i = 1:length(charge_events)
        for t_idx = 1:length(R_timepoints)
            R_field = sprintf('R_%ds', R_timepoints(t_idx));
            R_matrix_charge(i, t_idx) = charge_events(i).(R_field);
        end
    end
    
    plot(R_timepoints, R_matrix_charge, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
    hold on;
    R_mean = mean(R_matrix_charge, 1, 'omitnan');
    plot(R_timepoints, R_mean, 'k-', 'LineWidth', 3);
    hold off;
end
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Resistance [mOhm]', 'FontSize', 12, 'FontWeight', 'bold');
title('Charge Events', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend([repmat({''}, 1, size(R_matrix_charge, 1)), {'Mean'}], 'Location', 'best');

subplot(1, 2, 2);
if ~isempty(discharge_events)
    R_matrix_discharge = zeros(length(discharge_events), length(R_timepoints));
    for i = 1:length(discharge_events)
        for t_idx = 1:length(R_timepoints)
            R_field = sprintf('R_%ds', R_timepoints(t_idx));
            R_matrix_discharge(i, t_idx) = discharge_events(i).(R_field);
        end
    end
    
    plot(R_timepoints, R_matrix_discharge, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
    hold on;
    R_mean = mean(R_matrix_discharge, 1, 'omitnan');
    plot(R_timepoints, R_mean, 'k-', 'LineWidth', 3);
    hold off;
end
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Resistance [mOhm]', 'FontSize', 12, 'FontWeight', 'bold');
title('Discharge Events', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend([repmat({''}, 1, size(R_matrix_discharge, 1)), {'Mean'}], 'Location', 'best');

sgtitle('Time-dependent Resistance (2021 June)', 'FontSize', 16, 'FontWeight', 'bold');
saveas(fig1, fullfile(figDir, 'Resistance_vs_Time_2021Jun.fig'));
fprintf('  저장: Resistance_vs_Time_2021Jun.fig\\n');
close(fig1);

% 2. 전류 그룹별 저항 분포 (10초 저항 기준)
fig2 = figure('Position', [100, 100, 1200, 600]);

subplot(1, 2, 1);
if ~isempty(charge_events)
    groups_c = [charge_events.current_group];
    R10_c = [charge_events.R_10s];
    
    for g = 1:length(current_groups.names)
        idx_g = (groups_c == g) & ~isnan(R10_c);
        if sum(idx_g) > 0
            boxplot(R10_c(idx_g), 'Positions', g, 'Widths', 0.5);
            hold on;
        end
    end
    hold off;
    
    set(gca, 'XTick', 1:length(current_groups.names), ...
             'XTickLabel', current_groups.names, 'XTickLabelRotation', 45);
end
ylabel('R_{10s} [mOhm]', 'FontSize', 12, 'FontWeight', 'bold');
title('Charge Events', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

subplot(1, 2, 2);
if ~isempty(discharge_events)
    groups_d = [discharge_events.current_group];
    R10_d = [discharge_events.R_10s];
    
    for g = 1:length(current_groups.names)
        idx_g = (groups_d == g) & ~isnan(R10_d);
        if sum(idx_g) > 0
            boxplot(R10_d(idx_g), 'Positions', g, 'Widths', 0.5);
            hold on;
        end
    end
    hold off;
    
    set(gca, 'XTick', 1:length(current_groups.names), ...
             'XTickLabel', current_groups.names, 'XTickLabelRotation', 45);
end
ylabel('R_{10s} [mOhm]', 'FontSize', 12, 'FontWeight', 'bold');
title('Discharge Events', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

sgtitle('Resistance by Current Group (2021 June)', 'FontSize', 16, 'FontWeight', 'bold');
saveas(fig2, fullfile(figDir, 'Resistance_by_CurrentGroup_2021Jun.fig'));
fprintf('  저장: Resistance_by_CurrentGroup_2021Jun.fig\\n');
close(fig2);

%% 결과 저장
fprintf('\\n[5단계] 결과 저장 중...\\n');

results = struct();
results.events = events;
results.charge_events = charge_events;
results.discharge_events = discharge_events;
results.parameters.idle_threshold = idle_threshold;
results.parameters.min_duration = min_duration;
results.parameters.max_I_std = max_I_std;
results.parameters.R_timepoints = R_timepoints;
results.parameters.current_groups = current_groups;
results.filtered_stats = filtered_stats;

save(fullfile(savePath, 'Events_Resistance_2021Jun_v2.mat'), 'results', '-v7.3');
fprintf('  저장: Events_Resistance_2021Jun_v2.mat\\n');

fprintf('\\n========================================\\n');
fprintf('         분석 완료!\\n');
fprintf('========================================\\n');
