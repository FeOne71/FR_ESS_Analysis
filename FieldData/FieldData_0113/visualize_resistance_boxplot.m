%% 시간별 저항 박스플롯 시각화 (그룹별)
% 목적: 비슷한 조건(SOC, 전류)끼리 그룹화하여 시간별 저항 변화 확인

clear; clc; close all;

%% 데이터 로드
savePath = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_0113';
load(fullfile(savePath, 'Events_Resistance_2021Jun_v2.mat'), 'results');

fprintf('========================================\\n');
fprintf('  시간별 저항 박스플롯 시각화\\n');
fprintf('========================================\\n\\n');

events = results.events;
charge_events = results.charge_events;
discharge_events = results.discharge_events;
R_timepoints = results.parameters.R_timepoints;

fprintf('로드 완료:\\n');
fprintf('  총 이벤트: %d개\\n', length(events));
fprintf('  충전: %d개, 방전: %d개\\n\\n', length(charge_events), length(discharge_events));

%% 그룹화 기준 설정
% SOC 그룹
soc_groups = struct();
soc_groups.centers = [40, 50, 60, 70];
soc_groups.tolerance = 5;  % ±5%
soc_groups.names = {'SOC 35-45%', 'SOC 45-55%', 'SOC 55-65%', 'SOC 65-75%'};

% 전류 그룹 (Rack 기준)
current_groups = struct();
current_groups.ranges = [15 30; 30 50; 50 100];
current_groups.names = {'15-30A', '30-50A', '50-100A'};

fprintf('그룹화 기준:\\n');
fprintf('  SOC: %d개 그룹\\n', length(soc_groups.centers));
fprintf('  전류: %d개 그룹\\n\\n', size(current_groups.ranges, 1));

%% 충전 이벤트 그룹화
fprintf('[1단계] 충전 이벤트 그룹화\\n');

charge_grouped = struct();
group_count = 0;

for soc_idx = 1:length(soc_groups.centers)
    soc_center = soc_groups.centers(soc_idx);
    soc_min = soc_center - soc_groups.tolerance;
    soc_max = soc_center + soc_groups.tolerance;
    
    for curr_idx = 1:size(current_groups.ranges, 1)
        curr_min = current_groups.ranges(curr_idx, 1);
        curr_max = current_groups.ranges(curr_idx, 2);
        
        % 조건에 맞는 이벤트 찾기
        valid_idx = [];
        for i = 1:length(charge_events)
            evt = charge_events(i);
            if evt.SOC_mean >= soc_min && evt.SOC_mean <= soc_max && ...
               evt.I_avg_3_30s >= curr_min && evt.I_avg_3_30s < curr_max
                valid_idx = [valid_idx; i];
            end
        end
        
        if ~isempty(valid_idx)
            group_count = group_count + 1;
            charge_grouped(group_count).soc_idx = soc_idx;
            charge_grouped(group_count).curr_idx = curr_idx;
            charge_grouped(group_count).soc_name = soc_groups.names{soc_idx};
            charge_grouped(group_count).curr_name = current_groups.names{curr_idx};
            charge_grouped(group_count).events = charge_events(valid_idx);
            charge_grouped(group_count).n_events = length(valid_idx);
            
            fprintf('  그룹 %d: %s × %s → %d개\\n', ...
                    group_count, soc_groups.names{soc_idx}, ...
                    current_groups.names{curr_idx}, length(valid_idx));
        end
    end
end

fprintf('\\n충전 이벤트: 총 %d개 그룹\\n\\n', group_count);

%% 방전 이벤트 그룹화
fprintf('[2단계] 방전 이벤트 그룹화\\n');

discharge_grouped = struct();
group_count_d = 0;

for soc_idx = 1:length(soc_groups.centers)
    soc_center = soc_groups.centers(soc_idx);
    soc_min = soc_center - soc_groups.tolerance;
    soc_max = soc_center + soc_groups.tolerance;
    
    for curr_idx = 1:size(current_groups.ranges, 1)
        curr_min = current_groups.ranges(curr_idx, 1);
        curr_max = current_groups.ranges(curr_idx, 2);
        
        valid_idx = [];
        for i = 1:length(discharge_events)
            evt = discharge_events(i);
            if evt.SOC_mean >= soc_min && evt.SOC_mean <= soc_max && ...
               evt.I_avg_3_30s >= curr_min && evt.I_avg_3_30s < curr_max
                valid_idx = [valid_idx; i];
            end
        end
        
        if ~isempty(valid_idx)
            group_count_d = group_count_d + 1;
            discharge_grouped(group_count_d).soc_idx = soc_idx;
            discharge_grouped(group_count_d).curr_idx = curr_idx;
            discharge_grouped(group_count_d).soc_name = soc_groups.names{soc_idx};
            discharge_grouped(group_count_d).curr_name = current_groups.names{curr_idx};
            discharge_grouped(group_count_d).events = discharge_events(valid_idx);
            discharge_grouped(group_count_d).n_events = length(valid_idx);
            
            fprintf('  그룹 %d: %s × %s → %d개\\n', ...
                    group_count_d, soc_groups.names{soc_idx}, ...
                    current_groups.names{curr_idx}, length(valid_idx));
        end
    end
end

fprintf('\\n방전 이벤트: 총 %d개 그룹\\n\\n', group_count_d);

%% 충전 이벤트 시각화
fprintf('[3단계] 충전 이벤트 박스플롯 생성\\n');

figDir = fullfile(savePath, 'figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

% 그룹별로 figure 생성
for g = 1:length(charge_grouped)
    grp = charge_grouped(g);
    
    if grp.n_events < 3
        fprintf('  그룹 %d: 이벤트 부족 (N=%d), 스킵\\n', g, grp.n_events);
        continue;
    end
    
    fprintf('  그룹 %d 시각화: %s × %s (N=%d)\\n', ...
            g, grp.soc_name, grp.curr_name, grp.n_events);
    
    % 시간별 저항 데이터 수집
    R_data = zeros(grp.n_events, length(R_timepoints));
    for i = 1:grp.n_events
        evt = grp.events(i);
        for t_idx = 1:length(R_timepoints)
            R_field = sprintf('R_%ds', R_timepoints(t_idx));
            R_data(i, t_idx) = evt.(R_field);
        end
    end
    
    % Figure 생성
    fig = figure('Position', [100, 100, 1000, 600]);
    
    % 박스플롯
    boxplot(R_data, 'Labels', {'1s', '3s', '5s', '10s', '30s', '60s'});
    
    xlabel('Time', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Resistance [mOhm]', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Charge Events: %s, %s (N=%d)', ...
                  grp.soc_name, grp.curr_name, grp.n_events), ...
          'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % 통계 표시 (텍스트 박스)
    stats_text = {};
    for t_idx = 1:length(R_timepoints)
        R_vals = R_data(:, t_idx);
        R_vals = R_vals(~isnan(R_vals));
        if ~isempty(R_vals)
            stats_text{end+1} = sprintf('R_%ds: %.2f±%.2f', ...
                                         R_timepoints(t_idx), ...
                                         mean(R_vals), std(R_vals));
        end
    end
    
    annotation('textbox', [0.15, 0.75, 0.25, 0.15], ...
               'String', stats_text, 'FontSize', 9, ...
               'BackgroundColor', 'w', 'EdgeColor', 'k');
    
    % 저장
    filename = sprintf('Charge_Resistance_Boxplot_SOC%d_I%d.fig', ...
                       grp.soc_idx, grp.curr_idx);
    saveas(fig, fullfile(figDir, filename));
    fprintf('    저장: %s\\n', filename);
    
    close(fig);
end

fprintf('\\n충전 이벤트 시각화 완료\\n\\n');
