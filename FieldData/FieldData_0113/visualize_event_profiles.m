%% 이벤트 전류-전압 프로파일 시각화
% 목적: 이벤트 품질 확인 및 이상 데이터 탐지

clear; clc; close all;

%% 데이터 로드
savePath = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_0113';
dataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old\2021\202106';

load(fullfile(savePath, 'Events_Resistance_2021Jun_v2.mat'), 'results');

fprintf('========================================\n');
fprintf('  이벤트 프로파일 시각화\n');
fprintf('========================================\n\n');

events = results.events;
charge_events = results.charge_events;

% 전체 데이터 로드 (인덱스로 접근하기 위해)
fprintf('전체 데이터 로딩 중...\n');
matFiles = dir(fullfile(dataPath, '*.mat'));

allData = struct();
allData.voltage = [];
allData.current = [];
allData.time_seconds = [];

for day = 1:length(matFiles)
    if mod(day, 5) == 0
        fprintf('  %d/%d\n', day, length(matFiles));
    end
    filePath = fullfile(dataPath, matFiles(day).name);
    data = load(filePath);
    rack01 = data.Raw.Rack01;
    
    allData.voltage = [allData.voltage; rack01.AverageCV_V];
    allData.current = [allData.current; rack01.DCCurrent_A];
    allData.time_seconds = [allData.time_seconds; (1:length(rack01.AverageCV_V))'];
end

fprintf('데이터 로드 완료\n\n');

%% 메인 그룹 이벤트 시각화 (SOC 55-65%, 15-30A)
fprintf('[메인 그룹 분석]\n');

% 메인 그룹 선택
soc_min = 50; soc_max = 70;
curr_min = 15; curr_max = 30;

valid_idx = [];
for i = 1:length(charge_events)
    evt = charge_events(i);
    if evt.SOC_mean >= soc_min && evt.SOC_mean <= soc_max && ...
       evt.I_avg_3_30s >= curr_min && evt.I_avg_3_30s < curr_max
        valid_idx = [valid_idx; i];
    end
end

main_events = charge_events(valid_idx);
fprintf('메인 그룹 이벤트: %d개\n\n', length(main_events));

%% 저항값 이상치 찾기
fprintf('[1단계] 저항값 이상치 탐지\n');

% R_10s 기준
R10_values = [main_events.R_10s];
R10_mean = mean(R10_values, 'omitnan');
R10_std = std(R10_values, 'omitnan');

fprintf('R_10s 통계:\n');
fprintf('  평균: %.3f mOhm\n', R10_mean);
fprintf('  표준편차: %.3f mOhm\n', R10_std);
fprintf('  범위: [%.3f ~ %.3f]\n\n', min(R10_values), max(R10_values));

% 이상치 기준: mean ± 2*std
outlier_threshold_low = R10_mean - 2 * R10_std;
outlier_threshold_high = R10_mean + 2 * R10_std;

outlier_idx = find(R10_values < outlier_threshold_low | R10_values > outlier_threshold_high);
normal_idx = find(R10_values >= outlier_threshold_low & R10_values <= outlier_threshold_high);

fprintf('이상치 탐지 (mean ± 2σ):\n');
fprintf('  정상: %d개\n', length(normal_idx));
fprintf('  이상치: %d개\n', length(outlier_idx));

if ~isempty(outlier_idx)
    fprintf('\\n이상치 목록:\n');
    for i = 1:length(outlier_idx)
        evt_idx = outlier_idx(i);
        fprintf('  이벤트 %d: R_10s = %.3f mOhm\n', ...
                evt_idx, main_events(evt_idx).R_10s);
    end
end
fprintf('\n');

%% 샘플 이벤트 시각화 (정상 5개 + 이상치 전부)
fprintf('[2단계] 이벤트 프로파일 시각화\n');

figDir = fullfile(savePath, 'figures');

% 정상 이벤트 5개 랜덤 선택
if length(normal_idx) >= 5
    sample_normal = randsample(normal_idx, 5);
else
    sample_normal = normal_idx;
end

% 이상치는 모두 선택 (최대 10개)
sample_outlier = outlier_idx(1:min(10, length(outlier_idx)));

sample_events = [sample_normal; sample_outlier];
fprintf('시각화 대상: 정상 %d개 + 이상치 %d개\n\n', ...
        length(sample_normal), length(sample_outlier));

% Figure 생성
n_events = length(sample_events);
n_cols = 3;
n_rows = ceil(n_events / n_cols);

fig = figure('Position', [50, 50, 1600, 300*n_rows]);

for i = 1:n_events
    evt_idx = sample_events(i);
    evt = main_events(evt_idx);
    
    % 데이터 추출
    idx_start = evt.idx_start;
    idx_end = evt.idx_end;
    idx_idle = evt.idx_idle;
    
    I_event = allData.current(idx_start:idx_end);
    V_event = allData.voltage(idx_start:idx_end);
    t_event = (0:length(I_event)-1)';
    
    % 서브플롯
    subplot(n_rows, n_cols, i);
    
    yyaxis left
    plot(t_event, I_event, 'b-', 'LineWidth', 1.5);
    ylabel('Current [A]', 'FontSize', 10);
    ylim([min(I_event)-5, max(I_event)+5]);
    
    yyaxis right
    plot(t_event, V_event, 'r-', 'LineWidth', 1.5);
    ylabel('Voltage [V]', 'FontSize', 10);
    
    xlabel('Time [s]', 'FontSize', 10);
    
    % 제목 (정상 vs 이상치 표시)
    if ismember(evt_idx, outlier_idx)
        title_str = sprintf('Event %d [OUTLIER]\\nR_{10s}=%.2f mOhm', ...
                            evt_idx, evt.R_10s);
        title(title_str, 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');
    else
        title_str = sprintf('Event %d [NORMAL]\\nR_{10s}=%.2f mOhm', ...
                            evt_idx, evt.R_10s);
        title(title_str, 'FontSize', 10, 'Color', 'k');
    end
    
    grid on;
end

sgtitle('Event Profiles: Main Group (SOC 55-65%, 15-30A)', ...
        'FontSize', 16, 'FontWeight', 'bold');

saveas(fig, fullfile(figDir, 'Event_Profiles_MainGroup_Sample.fig'));
fprintf('저장: Event_Profiles_MainGroup_Sample.fig\n');

close(fig);

%% 전류 그룹별 전체 이벤트 시각화
fprintf('\n[3단계] 전류 그룹별 전체 이벤트 오버레이\n');

% 전류 그룹 정의
current_groups = {[15, 20], [20, 25], [25, 30]};
group_names = {'15-20A', '20-25A', '25-30A'};

fig2 = figure('Position', [100, 100, 1400, 900]);

for g = 1:length(current_groups)
    curr_range = current_groups{g};
    
    % 해당 그룹 이벤트 찾기
    group_events_idx = [];
    for i = 1:length(main_events)
        if main_events(i).I_avg_3_30s >= curr_range(1) && ...
           main_events(i).I_avg_3_30s < curr_range(2)
            group_events_idx = [group_events_idx; i];
        end
    end
    
    fprintf('  그룹 %s: %d개 이벤트\n', group_names{g}, length(group_events_idx));
    
    if isempty(group_events_idx)
        continue;
    end
    
    % 전류 서브플롯
    subplot(3, 2, 2*g-1);
    hold on;
    for i = 1:length(group_events_idx)
        evt_idx = group_events_idx(i);
        evt = main_events(evt_idx);
        
        I_event = allData.current(evt.idx_start:evt.idx_end);
        t_event = (0:length(I_event)-1)';
        
        % 이상치는 빨간색
        if ismember(evt_idx, outlier_idx)
            plot(t_event, I_event, 'r-', 'LineWidth', 1.5);
        else
            plot(t_event, I_event, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
        end
    end
    hold off;
    xlabel('Time [s]', 'FontSize', 11);
    ylabel('Current [A]', 'FontSize', 11);
    title(sprintf('Current: %s (N=%d)', group_names{g}, length(group_events_idx)), ...
          'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % 전압 서브플롯
    subplot(3, 2, 2*g);
    hold on;
    for i = 1:length(group_events_idx)
        evt_idx = group_events_idx(i);
        evt = main_events(evt_idx);
        
        V_event = allData.voltage(evt.idx_start:evt.idx_end);
        t_event = (0:length(V_event)-1)';
        
        if ismember(evt_idx, outlier_idx)
            plot(t_event, V_event, 'r-', 'LineWidth', 1.5);
        else
            plot(t_event, V_event, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
        end
    end
    hold off;
    xlabel('Time [s]', 'FontSize', 11);
    ylabel('Voltage [V]', 'FontSize', 11);
    title(sprintf('Voltage: %s (N=%d)', group_names{g}, length(group_events_idx)), ...
          'FontSize', 12, 'FontWeight', 'bold');
    grid on;
end

sgtitle('Event Profiles by Current Group (Red = Outlier)', ...
        'FontSize', 16, 'FontWeight', 'bold');

saveas(fig2, fullfile(figDir, 'Event_Profiles_by_CurrentGroup.fig'));
fprintf('\\n저장: Event_Profiles_by_CurrentGroup.fig\n');

close(fig2);

fprintf('\n========================================\n');
fprintf('         완료!\n');
fprintf('========================================\n');
