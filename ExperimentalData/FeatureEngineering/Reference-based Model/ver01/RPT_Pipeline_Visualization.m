%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Pipeline Visualization
% Phase별 전처리 및 피처 추출 과정 시각화
% 데이터소스, Master Ruler, 피처 추출, 라벨, 정규화
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
warning on;

%% ========================================================================
% Configuration
% ========================================================================
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_static_mat = fullfile(baseDir, 'Capacity_Trend_Figures', 'Capacity_Data_Static.mat');
path_rpt_vq_mat = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');
path_feature_mat = fullfile(baseDir, 'FeatureEngineering', 'Lab_RPT_Analysis', 'Dataset', 'Feature_Matrix_Final.mat');
saveDir = fullfile(baseDir, 'FeatureEngineering', 'Lab_RPT_Analysis', 'Pipeline_Visualizations');

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

fprintf('Loading Data...\n');
load(path_static_mat, 'allChannelsCapacity');
load(path_rpt_vq_mat, 'RPT_VQ_grid');
load(path_feature_mat, 'FeatureTable');
fprintf('Data Loaded.\n');

channels = fieldnames(allChannelsCapacity);
cyc_fields = fieldnames(RPT_VQ_grid);

%% ========================================================================
% Phase 0: 데이터소스 시각화
% ========================================================================
fprintf('\n=== Phase 0: 데이터소스 시각화 ===\n');

% --- 1. 원본 Raw vs 0.001V 보간 비교 ---
fig1 = figure('Position', [50, 50, 1600, 800], 'Name', 'Phase 0-1: Raw vs Interpolated Data');

% 샘플 데이터: Ch09, cyc0, 0.5C charge
sample_data = RPT_VQ_grid.cyc0.Ch09.c05_charge;

subplot(2,3,1);
plot(sample_data.Q_raw, sample_data.V_raw, 'o-', 'MarkerSize', 3, 'LineWidth', 1);
ylabel('Voltage (V)', 'FontWeight', 'bold');
xlabel('Capacity (Ah)', 'FontWeight', 'bold');
title('원본 Raw 데이터 (불규칙 샘플링)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend(sprintf('Raw 포인트: %d개', length(sample_data.V_raw)), 'Location', 'best');

subplot(2,3,2);
plot(sample_data.Q,sample_data.V_grid, '-', 'LineWidth', 1.5);
ylabel('Voltage (V)', 'FontWeight', 'bold');
xlabel('Capacity (Ah)', 'FontWeight', 'bold');
title('0.001V 균일 그리드 보간', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend(sprintf('보간 포인트: %d개', length(sample_data.V_grid)), 'Location', 'best');

subplot(2,3,3);
hold on;
plot(sample_data.Q_raw,sample_data.V_raw, 'o', 'MarkerSize', 4, 'DisplayName', 'Raw');
plot( sample_data.Q,sample_data.V_grid, '-', 'LineWidth', 1.5, 'DisplayName', '0.001V 보간');
ylabel('Voltage (V)', 'FontWeight', 'bold');
xlabel('Capacity (Ah)', 'FontWeight', 'bold');
title('오버레이 비교', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;

% 샘플링 간격 분포
subplot(2,3,4);
raw_intervals = diff(sample_data.V_raw);
histogram(raw_intervals, 30, 'FaceColor', [0.8 0.3 0.3]);
xlabel('전압 간격 (V)', 'FontWeight', 'bold');
ylabel('빈도', 'FontWeight', 'bold');
title(sprintf('Raw 샘플링 간격 분포\n평균: %.4fV, 표준편차: %.4fV', mean(raw_intervals), std(raw_intervals)), ...
    'FontSize', 11, 'FontWeight', 'bold');
grid on;

subplot(2,3,5);
grid_intervals = diff(sample_data.V_grid);
histogram(grid_intervals, 30, 'FaceColor', [0.3 0.8 0.3]);
xlabel('전압 간격 (V)', 'FontWeight', 'bold');
ylabel('빈도', 'FontWeight', 'bold');
title(sprintf('보간 후 간격 분포\n평균: %.6fV, 표준편차: %.8fV', mean(grid_intervals), std(grid_intervals)), ...
    'FontSize', 11, 'FontWeight', 'bold');
grid on;

% 보간 품질 (오차)
subplot(2,3,6);
Q_check = interp1(sample_data.V_raw, sample_data.Q_raw, sample_data.V_grid, 'linear', 'extrap');
interp_error = abs(sample_data.Q - Q_check);
plot(sample_data.V_grid, interp_error * 1000, 'LineWidth', 1.5);  % mAh 단위
xlabel('Voltage (V)', 'FontWeight', 'bold');
ylabel('보간 오차 (mAh)', 'FontWeight', 'bold');
title(sprintf('보간 품질 검증\n최대 오차: %.3f mAh', max(interp_error)*1000), ...
    'FontSize', 11, 'FontWeight', 'bold');
grid on;

sgtitle('Phase 0-1: 원본 Raw vs 0.001V 보간 데이터 비교 (Ch09, cyc0, 0.5C Charge)', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig1, fullfile(saveDir, 'Phase0_1_Raw_vs_Interpolated.fig'));
close(fig1);

% --- 2. 전체 채널 V-Q 곡선 오버뷰 ---
fig2 = figure('Position', [50, 50, 1600, 900], 'Name', 'Phase 0-2: All Channels VQ Overview');

% Charge 곡선
subplot(1,2,1);
hold on;
colors = lines(length(cyc_fields));
for c = 1:length(cyc_fields)
    cyc = cyc_fields{c};
    for i = 1:length(channels)
        ch = channels{i};
        if isfield(RPT_VQ_grid.(cyc), ch) && isfield(RPT_VQ_grid.(cyc).(ch), 'c05_charge')
            data = RPT_VQ_grid.(cyc).(ch).c05_charge;
            if i == 1
                plot(data.Q, data.V_grid, 'Color', colors(c,:), 'LineWidth', 1, 'DisplayName', cyc);
            else
                plot(data.Q, data.V_grid, 'Color', colors(c,:), 'LineWidth', 1, 'HandleVisibility', 'off');
            end
        end
    end
end
xlabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
title('Charge V-Q (0.5C, All Channels)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

% Discharge 곡선
subplot(1,2,2);
hold on;
for c = 1:length(cyc_fields)
    cyc = cyc_fields{c};
    for i = 1:length(channels)
        ch = channels{i};
        if isfield(RPT_VQ_grid.(cyc), ch) && isfield(RPT_VQ_grid.(cyc).(ch), 'c05_discharge')
            data = RPT_VQ_grid.(cyc).(ch).c05_discharge;
            if i == 1
                plot(data.Q, data.V_grid, 'Color', colors(c,:), 'LineWidth', 1, 'DisplayName', cyc);
            else
                plot(data.Q, data.V_grid, 'Color', colors(c,:), 'LineWidth', 1, 'HandleVisibility', 'off');
            end
        end
    end
end
ylabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
title('Discharge V-Q (0.5C, All Channels)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

sgtitle('Phase 0-2: 전체 V-Q 곡선 오버뷰 (8 Channels × 6 Cycles)', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig2, fullfile(saveDir, 'Phase0_2_All_VQ_Overview.fig'));
close(fig2);

% --- 3. 분석 전압 윈도우 표시 ---
fig3 = figure('Position', [50, 50, 1600, 600], 'Name', 'Phase 0-3: Voltage Windows');

win_chg_min = 3.7;
win_chg_max = 3.95;
win_dch_min = 3.75;
win_dch_max = 3.88;

% Charge 윈도우
subplot(1,2,1);
hold on;
for c = 1:length(cyc_fields)
    cyc = cyc_fields{c};
    if isfield(RPT_VQ_grid.(cyc), 'Ch09') && isfield(RPT_VQ_grid.(cyc).Ch09, 'c05_charge')
        data = RPT_VQ_grid.(cyc).Ch09.c05_charge;
        plot(data.Q, data.V_grid, 'Color', colors(c,:), 'LineWidth', 1.5, 'DisplayName', cyc);
    end
end
% 윈도우 영역 강조
xl = xlim;
fill([xl(1) xl(2) xl(2) xl(1)], [win_chg_min win_chg_min win_chg_max win_chg_max], ...
    'g', 'FaceAlpha', 0.1, 'EdgeColor', 'g', 'LineWidth', 2, 'LineStyle', '--', ...
    'DisplayName', sprintf('분석 윈도우\n[%.2f-%.2fV]', win_chg_min, win_chg_max));
xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
title('Charge 분석 윈도우', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

% Discharge 윈도우
subplot(1,2,2);
hold on;
for c = 1:length(cyc_fields)
    cyc = cyc_fields{c};
    if isfield(RPT_VQ_grid.(cyc), 'Ch09') && isfield(RPT_VQ_grid.(cyc).Ch09, 'c05_discharge')
        data = RPT_VQ_grid.(cyc).Ch09.c05_discharge;
        plot(data.Q, data.V_grid, 'Color', colors(c,:), 'LineWidth', 1.5, 'DisplayName', cyc);
    end
end
xl = xlim;
fill([xl(1) xl(2) xl(2) xl(1)], [win_dch_min win_dch_min win_dch_max win_dch_max], ...
    'r', 'FaceAlpha', 0.1, 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--', ...
    'DisplayName', sprintf('분석 윈도우\n[%.2f-%.2fV]', win_dch_min, win_dch_max));
xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
title('Discharge 분석 윈도우', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

sgtitle('Phase 0-3: dQ/dV 분석 전압 윈도우 (Ch09, 0.5C)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig3, fullfile(saveDir, 'Phase0_3_Voltage_Windows.fig'));
close(fig3);

% --- 4. 사이클별 용량 트렌드 ---
fig4 = figure('Position', [50, 50, 1400, 600], 'Name', 'Phase 0-4: Capacity Trends');

subplot(1,2,1);
hold on;
for i = 1:length(channels)
    ch = channels{i};
    cycles = allChannelsCapacity.(ch).cycles;
    capacities = cellfun(@max, allChannelsCapacity.(ch).Q);
    plot(cycles, capacities, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', ch);
end
xlabel('Cycle Number', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Maximum Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
title('사이클별 용량 감소 (SOH)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

subplot(1,2,2);
hold on;
for i = 1:length(channels)
    ch = channels{i};
    cycles = allChannelsCapacity.(ch).cycles;
    capacities = cellfun(@max, allChannelsCapacity.(ch).Q);
    if length(cycles) > 1
        % 정규화 (cyc0 = 100%)
        cap_normalized = (capacities / capacities(1)) * 100;
        plot(cycles, cap_normalized, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', ch);
    end
end
xlabel('Cycle Number', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Capacity Retention (%)', 'FontWeight', 'bold', 'FontSize', 12);
title('용량 보존율', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
ylim([85 105]);
hold off;

sgtitle('Phase 0-4: 사이클별 용량 열화 트렌드 (All Channels)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig4, fullfile(saveDir, 'Phase0_4_Capacity_Trends.fig'));
close(fig4);

fprintf('Phase 0 완료.\n');

%% ========================================================================
% Phase 1: Master Ruler 시각화
% ========================================================================
fprintf('\n=== Phase 1: Master Ruler 생성 시각화 ===\n');

% Master Ruler 재계산 (시각화용)
num_segments = 5;
win_chg_min = 3.7; win_chg_max = 3.95;
win_dch_min = 3.75; win_dch_max = 3.88;

%% 마스터 룰러 (Static Capacity 방전 데이터 기반, 충전/방전 동일 소스)
standard_V_grid = RPT_VQ_grid.cyc0.(channels{1}).Static.V_grid;

% 충전 윈도우
mask_chg = standard_V_grid >= win_chg_min & standard_V_grid <= win_chg_max;
V_standard_chg = standard_V_grid(mask_chg);

% 방전 윈도우
mask_dch = standard_V_grid >= win_dch_min & standard_V_grid <= win_dch_max;
V_standard_dch = standard_V_grid(mask_dch);

all_T_grids_chg = zeros(length(channels), length(V_standard_chg));
all_T_grids_dch = zeros(length(channels), length(V_standard_dch));
valid_ch = 0;

for i = 1:length(channels)
    ch = channels{i};
    if ~isfield(RPT_VQ_grid.cyc0, ch) || ~isfield(RPT_VQ_grid.cyc0.(ch), 'Static')
        continue;
    end
    valid_ch = valid_ch + 1;
    data_s = RPT_VQ_grid.cyc0.(ch).Static;
    V_s = data_s.V_grid;
    t_s = data_s.t;
    
    % 충전 윈도우
    mc = V_s >= win_chg_min & V_s <= win_chg_max;
    all_T_grids_chg(valid_ch, :) = interp1(V_s(mc), t_s(mc), V_standard_chg, 'linear')';
    
    % 방전 윈도우
    md = V_s >= win_dch_min & V_s <= win_dch_max;
    all_T_grids_dch(valid_ch, :) = interp1(V_s(md), t_s(md), V_standard_dch, 'linear')';
end

all_T_grids_chg = all_T_grids_chg(1:valid_ch, :);
all_T_grids_dch = all_T_grids_dch(1:valid_ch, :);

% 충전 Master Ruler
avg_T_chg = mean(all_T_grids_chg, 1, 'omitnan');
avg_T_chg = fillmissing(avg_T_chg, 'linear', 'EndValues', 'nearest');
[avg_T_chg, uid_c] = unique(avg_T_chg, 'stable');
V_chg_u = V_standard_chg(uid_c);
T_start_chg = min(avg_T_chg); T_end_chg = max(avg_T_chg);
target_Ts_chg = linspace(T_start_chg, T_end_chg, num_segments + 1);
Global_V_bounds_chg = interp1(avg_T_chg, V_chg_u, target_Ts_chg, 'linear');
Global_V_bounds_chg = sort(Global_V_bounds_chg, 'ascend');  % 충전: 오름차순

% 방전 Master Ruler
avg_T_dch = mean(all_T_grids_dch, 1, 'omitnan');
avg_T_dch = fillmissing(avg_T_dch, 'linear', 'EndValues', 'nearest');
[avg_T_dch, uid_d] = unique(avg_T_dch, 'stable');
V_dch_u = V_standard_dch(uid_d);
target_Ts_dch = linspace(min(avg_T_dch), max(avg_T_dch), num_segments + 1);
Global_V_bounds_dch = interp1(avg_T_dch, V_dch_u, target_Ts_dch, 'linear');

fprintf('  Master Ruler (Static): Chg %.4f~%.4fV, Dch %.4f~%.4fV\n', ...
    Global_V_bounds_chg(1), Global_V_bounds_chg(end), ...
    Global_V_bounds_dch(1), Global_V_bounds_dch(end));

% --- 1. Step 1: 각 채널의 원본 V-t 데이터 추출 ---
fig5 = figure('Position', [50, 50, 1600, 900], 'Name', 'Phase 1-1: Master Ruler Step 1');

% 상단: 8개 채널의 V-t 곡선 (원본)
% 상단: 8개 채널의 V-t 곡선 (원본 - 모든 사이클)
subplot(2,2,[1,2]);
hold on;
ch_colors = lines(8);
cyc_fields = fieldnames(RPT_VQ_grid);

for i = 1:length(channels)
    ch = channels{i};
    % Loop through all cycles for this channel
    for c = 1:length(cyc_fields)
        cyc = cyc_fields{c};
        if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && ...
           isfield(RPT_VQ_grid.(cyc).(ch), 'c05_charge')
            
            data = RPT_VQ_grid.(cyc).(ch).c05_charge;
            mask = data.V_grid >= win_chg_min & data.V_grid <= win_chg_max;
            V_ch = data.V_grid(mask);
            t_ch = data.t(mask);
            
            % Use transparency for cycles
            base_color = ch_colors(i,:);
            if strcmp(cyc, 'cyc0')
                plot(V_ch, t_ch, 'Color', base_color, 'LineWidth', 1.5, 'DisplayName', ch);
            else
                plot(V_ch, t_ch, 'Color', [base_color 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
            end
        end
    end
end
xlabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 12);
title('Step 1: 각 채널의 원본 V → t 관계 (모든 사이클 포함)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9, 'NumColumns', 2);
grid on;
hold off;

% 하단 좌: 특정 전압에서의 도달 시간 비교
subplot(2,2,3);
V_samples = [3.70, 3.75, 3.80, 3.85, 3.90, 3.95];
t_samples = zeros(length(channels), length(V_samples));
for i = 1:length(channels)
    for j = 1:length(V_samples)
        idx = find(V_standard_chg >= V_samples(j), 1);
        if ~isempty(idx)
            t_samples(i,j) = all_T_grids_chg(i, idx);
        end
    end
end
bar(t_samples');
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.2fV', x), V_samples, 'UniformOutput', false));
xlabel('전압', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('도달 시간 (s)', 'FontWeight', 'bold', 'FontSize', 11);
title('각 전압에 도달하는 시간 (채널별 차이, Cyc0)', 'FontSize', 12, 'FontWeight', 'bold');
legend(channels, 'Location', 'best', 'FontSize', 8, 'NumColumns', 2);
grid on;

% 하단 우: 소요 시간 분포
subplot(2,2,4);
durations = zeros(length(channels), 1);
for i = 1:length(channels)
    durations(i) = abs(all_T_grids_chg(i, end) - all_T_grids_chg(i, 1));
end
bar(durations);
set(gca, 'XTickLabel', channels);
xlabel('채널', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('총 소요 시간 (s)', 'FontWeight', 'bold', 'FontSize', 11);
title(sprintf('전압 범위(%.2f-%.2fV) 통과 시간 (Cyc0)', win_chg_min, win_chg_max), 'FontSize', 12, 'FontWeight', 'bold');
if range(durations) > 0
    ylim([min(durations)*0.95, max(durations)*1.05]);
end
grid on;
for i = 1:length(channels)
    text(i, durations(i), sprintf('%.0fs', durations(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
end

sgtitle('Phase 1-1: Master Ruler 생성 - Step 1 (원본 데이터)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig5, fullfile(saveDir, 'Phase1_1_MasterRuler_Step1.fig'));
close(fig5);

% --- 2. Step 2: 시간 그리드 평균화 과정 ---
fig6 = figure('Position', [50, 50, 1600, 900], 'Name', 'Phase 1-2: Master Ruler Step 2');

% 상단: 개별 채널+사이클 + 평균 V-t 곡선
subplot(2,2,[1,2]);
hold on;
for i = 1:size(all_T_grids_chg, 1)
    plot(V_standard_chg, all_T_grids_chg(i,:), 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
plot(V_standard_chg, avg_T_chg, 'r-', 'LineWidth', 3, 'DisplayName', '평균 (대표 배터리)');
xlabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 12);
title('Step 2: 시간 그리드 평균화 (8개 채널 → 대표값)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
grid on;
hold off;

% 하단 좌: 특정 전압에서의 평균 계산 시각화
subplot(2,2,3);
V_demo = [3.70, 3.80, 3.90];
V_demo_idx = zeros(size(V_demo));
for i = 1:length(V_demo)
    V_demo_idx(i) = find(V_standard_chg >= V_demo(i), 1);
end

t_demo = zeros(length(channels), length(V_demo));
for i = 1:length(channels)
    for j = 1:length(V_demo)
        t_demo(i,j) = all_T_grids_chg(i, V_demo_idx(j));
    end
end

hold on;
for j = 1:length(V_demo)
    scatter(ones(8,1)*j, t_demo(:,j), 80, 'filled', 'MarkerFaceAlpha', 0.6);
    t_mean = mean(t_demo(:,j));
    plot([j-0.3, j+0.3], [t_mean, t_mean], 'r-', 'LineWidth', 3);
    text(j+0.35, t_mean, sprintf('평균=%.0fs', t_mean), 'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold');
end
set(gca, 'XTick', 1:length(V_demo), 'XTickLabel', arrayfun(@(x) sprintf('%.2fV', x), V_demo, 'UniformOutput', false));
xlabel('전압', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('시간 (s)', 'FontWeight', 'bold', 'FontSize', 11);
title('평균 계산 예시 (3개 전압 포인트)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
hold off;

% 하단 우: 표준편차 분석
subplot(2,2,4);
std_T = std(all_T_grids_chg, 0, 1);
plot(V_standard_chg, std_T, 'b-', 'LineWidth', 2);
xlabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('표준편차 (s)', 'FontWeight', 'bold', 'FontSize', 11);
title(sprintf('전체 데이터 시간 변동성\n평균 표준편차: %.fs', mean(std_T)), 'FontSize', 12, 'FontWeight', 'bold');
grid on;

sgtitle('Phase 1-2: Master Ruler 생성 - Step 2 (평균화)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig6, fullfile(saveDir, 'Phase1_2_MasterRuler_Step2.fig'));
close(fig6);

% --- 3. Step 3-4: 시간 균등 분할 ---
fig7 = figure('Position', [50, 50, 1600, 900], 'Name', 'Phase 1-3: Master Ruler Step 3-4');

% 상단: 평균 V-t 곡선에 시간 경계 표시
subplot(2,1,1);
hold on;
plot(V_standard_chg, avg_T_chg, 'b-', 'LineWidth', 2.5, 'DisplayName', '평균 V-t');
yline(T_start_chg, 'g--', 'LineWidth', 2, 'Label', sprintf('시작: %.0fs', T_start_chg));
for i = 2:length(target_Ts_chg)-1
    yline(target_Ts_chg(i), 'r--', 'LineWidth', 1.5, ...
        'Label', sprintf('경계%d: %.0fs', i, target_Ts_chg(i)));
end
yline(T_end_chg, 'g--', 'LineWidth', 2, 'Label', sprintf('종료: %.0fs', T_end_chg));
xlabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 12);
title(sprintf('Step 3: 시간 균등 분할 (%.0fs ÷ 5 = %.0fs/세그먼트)', ...
    T_end_chg-T_start_chg, (T_end_chg-T_start_chg)/num_segments), 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

% 하단: 세그먼트별 시간 & 전압 정보
subplot(2,1,2);
seg_info = zeros(num_segments, 3);  % [시간길이, 전압시작, 전압끝]
for i = 1:num_segments
    seg_info(i,1) = target_Ts_chg(i+1) - target_Ts_chg(i);
    seg_info(i,2) = Global_V_bounds_chg(i);
    seg_info(i,3) = Global_V_bounds_chg(i+1);
end

hold on;
% 시간 길이 (bar)
yyaxis left;
bar(seg_info(:,1), 'FaceColor', [0.3 0.6 0.9], 'FaceAlpha', 0.7);
ylabel('시간 길이 (s)', 'FontWeight', 'bold', 'FontSize', 11);
ylim([0, max(seg_info(:,1))*1.2]);

% 전압 범위 (scatter)
yyaxis right;
V_ranges = seg_info(:,3) - seg_info(:,2);
scatter(1:num_segments, V_ranges*1000, 120, 'r', 'filled', 'MarkerFaceAlpha', 0.7);
ylabel('전압 간격 (mV)', 'FontWeight', 'bold', 'FontSize', 11);

set(gca, 'XTick', 1:num_segments, 'XTickLabel', arrayfun(@(x) sprintf('Seg%d', x), 1:num_segments, 'UniformOutput', false));
xlabel('세그먼트', 'FontWeight', 'bold', 'FontSize', 11);
title('Step 4: 각 세그먼트의 시간(균등) vs 전압(불균등) 비교', 'FontSize', 12, 'FontWeight', 'bold');
legend({'시간 길이 (균등!)', '전압 간격 (불균등)'}, 'Location', 'best', 'FontSize', 10);
grid on;
hold off;

sgtitle('Phase 1-3: Master Ruler 생성 - Step 3-4 (시간 분할 & 전압 변환)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig7, fullfile(saveDir, 'Phase1_3_MasterRuler_Step3_4.fig'));
close(fig7);

% --- 4. 최종 결과 (충전/방전 세그먼트 - 전체 범위, 전체 채널) ---
fig8 = figure('Position', [50, 50, 1800, 900], 'Name', 'Phase 1-4: Master Ruler Final Result');

% 세그먼트 색상 정의
seg_colors = [
    0.6 0.8 1.0;  % 파란색 계열 - Seg1
    1.0 0.7 0.7;  % 빨간색 계열 - Seg2
    0.7 0.7 1.0;  % 파란색 계열 - Seg3
    1.0 0.8 0.6;  % 주황색 계열 - Seg4
    0.8 0.6 0.9;  % 보라색 계열 - Seg5
];

ch_colors = lines(length(channels));

%% === 상단 좌: 충전 V-Q 전체 (3.0-4.2V) - 전체 채널/사이클 ===
subplot(2,2,1);
hold on;
% 전체 채널/사이클 곡선 그리기
for ch_idx = 1:length(channels)
    ch = channels{ch_idx};
    for cyc_idx = 1:length(cyc_fields)
        cyc = cyc_fields{cyc_idx};
        if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && ...
           isfield(RPT_VQ_grid.(cyc).(ch), 'c05_charge')
            data = RPT_VQ_grid.(cyc).(ch).c05_charge;
            if strcmp(cyc, 'cyc0')
                plot(data.Q, data.V_grid, 'Color', ch_colors(ch_idx,:), 'LineWidth', 1.5);
            else
                plot(data.Q, data.V_grid, 'Color', [ch_colors(ch_idx,:) 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
            end
        end
    end
end
xl = xlim;

% 세그먼트 색상 칠하기
for i = 1:num_segments
    V_start = Global_V_bounds_chg(i);
    V_end = Global_V_bounds_chg(i+1);
    patch([xl(1) xl(2) xl(2) xl(1)], [V_start V_start V_end V_end], ...
        seg_colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% 세그먼트 경계선
for i = 1:length(Global_V_bounds_chg)
    line(xl, [Global_V_bounds_chg(i) Global_V_bounds_chg(i)], ...
        'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
end

xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
title('Charge Voltage Segment (All Cycles)', 'FontSize', 13, 'FontWeight', 'bold');
ylim([3.0 4.2]);
grid on;
hold off;

%% === 상단 우: 방전 V-Q 전체 (3.0-4.2V) - 전체 채널/사이클 ===
subplot(2,2,2);
hold on;
% 전체 채널/사이클 곡선 그리기
for ch_idx = 1:length(channels)
    ch = channels{ch_idx};
    for cyc_idx = 1:length(cyc_fields)
        cyc = cyc_fields{cyc_idx};
        if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && ...
           isfield(RPT_VQ_grid.(cyc).(ch), 'c05_discharge')
            data = RPT_VQ_grid.(cyc).(ch).c05_discharge;
            if strcmp(cyc, 'cyc0')
                plot(data.Q, data.V_grid, 'Color', ch_colors(ch_idx,:), 'LineWidth', 1.5);
            else
                plot(data.Q, data.V_grid, 'Color', [ch_colors(ch_idx,:) 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
            end
        end
    end
end
xl = xlim;

% 세그먼트 색상 칠하기
for i = 1:num_segments
    V_start = Global_V_bounds_dch(i);
    V_end = Global_V_bounds_dch(i+1);
    patch([xl(1) xl(2) xl(2) xl(1)], [V_start V_start V_end V_end], ...
        seg_colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% 세그먼트 경계선
for i = 1:length(Global_V_bounds_dch)
    line(xl, [Global_V_bounds_dch(i) Global_V_bounds_dch(i)], ...
        'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
end

xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
title('Discharge Voltage Segment (All Cycles)', 'FontSize', 13, 'FontWeight', 'bold');
ylim([3.0 4.2]);
grid on;
hold off;

%% === 하단 좌: 충전 확대 (3.70-3.95V) - 전체 채널/사이클 ===
subplot(2,2,3);
hold on;
for ch_idx = 1:length(channels)
    ch = channels{ch_idx};
    for cyc_idx = 1:length(cyc_fields)
        cyc = cyc_fields{cyc_idx};
        if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && ...
           isfield(RPT_VQ_grid.(cyc).(ch), 'c05_charge')
            data = RPT_VQ_grid.(cyc).(ch).c05_charge;
            mask = data.V_grid >= 3.70 & data.V_grid <= 3.95;
            if strcmp(cyc, 'cyc0')
                plot(data.Q(mask), data.V_grid(mask), 'Color', ch_colors(ch_idx,:), 'LineWidth', 1.5);
            else
                plot(data.Q(mask), data.V_grid(mask), 'Color', [ch_colors(ch_idx,:) 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
            end
        end
    end
end
xl = xlim;

% 세그먼트 색상 및 경계
for i = 1:num_segments
    V_start = Global_V_bounds_chg(i);
    V_end = Global_V_bounds_chg(i+1);
    if V_end >= 3.70 && V_start <= 3.95
        V_plot_start = max(V_start, 3.70);
        V_plot_end = min(V_end, 3.95);
        patch([xl(1) xl(2) xl(2) xl(1)], ...
            [V_plot_start V_plot_start V_plot_end V_plot_end], ...
            seg_colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
end

for i = 1:length(Global_V_bounds_chg)
    if Global_V_bounds_chg(i) >= 3.70 && Global_V_bounds_chg(i) <= 3.95
        line(xl, [Global_V_bounds_chg(i) Global_V_bounds_chg(i)], ...
            'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
    end
end

xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 11);
title(sprintf('Charge Segment Zoom (All Cycles)'), 'FontSize', 12, 'FontWeight', 'bold');
ylim([3.70 3.95]);
grid on;
hold off;

%% === 하단 우: 방전 확대 (3.75-3.88V) - 전체 채널/사이클 ===
subplot(2,2,4);
hold on;
for ch_idx = 1:length(channels)
    ch = channels{ch_idx};
    for cyc_idx = 1:length(cyc_fields)
        cyc = cyc_fields{cyc_idx};
        if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && ...
           isfield(RPT_VQ_grid.(cyc).(ch), 'c05_discharge')
            data = RPT_VQ_grid.(cyc).(ch).c05_discharge;
            mask = data.V_grid >= win_dch_min & data.V_grid <= win_dch_max;
            if strcmp(cyc, 'cyc0')
                plot(data.Q(mask), data.V_grid(mask), 'Color', ch_colors(ch_idx,:), 'LineWidth', 1.5);
            else
                plot(data.Q(mask), data.V_grid(mask), 'Color', [ch_colors(ch_idx,:) 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
            end
        end
    end
end
xl = xlim;

% 세그먼트 색상 및 경계
for i = 1:num_segments
    V_start = Global_V_bounds_dch(i);
    V_end = Global_V_bounds_dch(i+1);
    if V_end >= win_dch_min && V_start <= win_dch_max
        V_plot_start = max(V_start, win_dch_min);
        V_plot_end = min(V_end, win_dch_max);
        patch([xl(1) xl(2) xl(2) xl(1)], ...
            [V_plot_start V_plot_start V_plot_end V_plot_end], ...
            seg_colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
end

for i = 1:length(Global_V_bounds_dch)
    if Global_V_bounds_dch(i) >= win_dch_min && Global_V_bounds_dch(i) <= win_dch_max
        line(xl, [Global_V_bounds_dch(i) Global_V_bounds_dch(i)], ...
            'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
    end
end

xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 11);
title(sprintf('Discharge Segment Zoom (All Cycles)'), 'FontSize', 12, 'FontWeight', 'bold');
ylim([win_dch_min win_dch_max]);
grid on;
hold off;

sgtitle('Phase 1-4: Master Ruler Final Result (All Channels)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig8, fullfile(saveDir, 'Phase1_4_MasterRuler_Final.fig'));
close(fig8);

fprintf('Phase 1 완료.\n');

%% ========================================================================
% Phase 2: 피처 추출 시각화 (All Channels)
%% ========================================================================

% fprintf('\n=== Phase 2: 피처 추출 시각화 (All Channels & All Cycles) ===\n');
% 
% for ch_i = 1:length(channels)
%     ch = channels{ch_i};
%     for cyc_i = 1:length(cyc_fields)
%         cyc = cyc_fields{cyc_i};
%         fprintf('  Processing %s - %s...\n', ch, cyc);
% 
%         % --- 1. dQ/dV 필터링 효과 ---
%         fig7 = figure('Position', [50, 50, 1600, 900], 'Name', sprintf('Phase 2-1: Filtering Effect - %s %s', ch, cyc));
% 
%         if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && isfield(RPT_VQ_grid.(cyc).(ch), 'c05_charge')
%             sample_data = RPT_VQ_grid.(cyc).(ch).c05_charge;
%             V_grid = sample_data.V_grid;
%             Q_val = sample_data.Q;
%             [V_u, uid] = unique(V_grid);
%             Q_u = Q_val(uid);
% 
%             % 수치 미분
%             dV = gradient(V_u);
%             dQ = gradient(Q_u);
%             dV(dV==0) = NaN;
%             dQdV_raw = dQ ./ dV;
% 
%             % Smoothing 필터 (통일된 Moving Average)
%             dQdV_filt21 = movmean(dQdV_raw, 21);
%             dQdV_filt51 = movmean(dQdV_raw, 51);
% 
%             % 윈도우 영역
%             mask_win = V_u >= win_chg_min & V_u <= win_chg_max;
% 
%             subplot(2,2,1);
%             plot(V_u, dQdV_raw, 'LineWidth', 1);
%             xlabel('Voltage (V)', 'FontWeight', 'bold');
%             ylabel('dQ/dV (Ah/V)', 'FontWeight', 'bold');
%             title('원본 dQ/dV (노이즈 多)', 'FontSize', 12, 'FontWeight', 'bold');
%             grid on;
%             xlim([3.0 4.2]);
% 
%             subplot(2,2,2);
%             plot(V_u(mask_win), dQdV_raw(mask_win), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
%             xlabel('Voltage (V)', 'FontWeight', 'bold');
%             ylabel('dQ/dV (Ah/V)', 'FontWeight', 'bold');
%             title('분석 영역 (Raw)', 'FontSize', 12, 'FontWeight', 'bold');
%             grid on;
%             xlim([win_chg_min win_chg_max]);
% 
%             subplot(2,2,3);
%             plot(V_u(mask_win), dQdV_filt51(mask_win), 'LineWidth', 1.5);
%             xlabel('Voltage (V)', 'FontWeight', 'bold');
%             ylabel('dQ/dV (Ah/V)', 'FontWeight', 'bold');
%             title('MovMean (Window=51, 과도 스무딩)', 'FontSize', 12, 'FontWeight', 'bold');
%             grid on;
%             xlim([win_chg_min win_chg_max]);
% 
%             subplot(2,2,4);
%             hold on;
%             plot(V_u(mask_win), dQdV_raw(mask_win), 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'DisplayName', 'Raw');
%             plot(V_u(mask_win), dQdV_filt21(mask_win), 'b-', 'LineWidth', 2, 'DisplayName', 'MovMean-21 (사용)');
%             xlabel('Voltage (V)', 'FontWeight', 'bold');
%             ylabel('dQ/dV (Ah/V)', 'FontWeight', 'bold');
%             title('필터링 효과 비교 (분석 윈도우)', 'FontSize', 12, 'FontWeight', 'bold');
%             legend('Location', 'best', 'FontSize', 10);
%             grid on;
%             hold off;
% 
%             sgtitle(sprintf('Phase 2-1: 필터링 효과 (%s, %s, 0.5C Charge)', ch, cyc), ...
%                 'FontSize', 14, 'FontWeight', 'bold');
%             saveas(fig7, fullfile(saveDir, sprintf('Phase2_1_Filtering_Effect_%s_%s.fig', ch, cyc)));
%             close(fig7);
%         else
%             fprintf('  Skipping %s %s (Phase 2-1: No data)\n', ch, cyc);
%             close(fig7);
%         end
% 
%         % --- 2. 피크 검출 결과 ---
%         fig8 = figure('Position', [50, 50, 1400, 600], 'Name', sprintf('Phase 2-2: Peak Detection - %s %s', ch, cyc));
% 
%         if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && isfield(RPT_VQ_grid.(cyc).(ch), 'c05_charge')
%             % Charge 피크
%             pk_height_chg = max(dQdV_filt21(mask_win));
%             pk_area_chg = trapz(V_u(mask_win), dQdV_filt21(mask_win));
%             [~, pk_idx] = max(dQdV_filt21(mask_win));
%             V_win = V_u(mask_win);
%             pk_V_chg = V_win(pk_idx);
% 
%             subplot(1,2,1);
%             hold on;
%             plot(V_u(mask_win), dQdV_filt21(mask_win), 'b-', 'LineWidth', 2);
%             plot(pk_V_chg, pk_height_chg, 'ro', 'MarkerSize', 12, 'LineWidth', 2);
%             plot([V_win(1) V_win(end)], [0 0], 'k--', 'LineWidth', 1);
%             fill([V_win; flipud(V_win)], [dQdV_filt21(mask_win); zeros(sum(mask_win),1)], ...
%                 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%             text(pk_V_chg, pk_height_chg*1.1, sprintf('피크\n%.3fV\n%.1f Ah/V', pk_V_chg, pk_height_chg), ...
%                 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
%             xlabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
%             ylabel('dQ/dV (Ah/V)', 'FontWeight', 'bold', 'FontSize', 12);
%             title(sprintf('Charge 피크 검출 (%s, %s)\n높이: %.1f Ah/V, 면적: %.2f Ah', ch, cyc, pk_height_chg, pk_area_chg), ...
%                 'FontSize', 12, 'FontWeight', 'bold');
%             grid on;
%             hold off;
% 
%             % Discharge도 동일하게
%             if isfield(RPT_VQ_grid.(cyc).(ch), 'c05_discharge')
%                 sample_data_dch = RPT_VQ_grid.(cyc).(ch).c05_discharge;
%                 V_grid_dch = sample_data_dch.V_grid;
%                 Q_val_dch = sample_data_dch.Q;
%                 [V_u_dch, uid_dch] = unique(V_grid_dch);
%                 Q_u_dch = Q_val_dch(uid_dch);
%                 dV_dch = gradient(V_u_dch);
%                 dQ_dch = gradient(Q_u_dch);
%                 dV_dch(dV_dch==0) = NaN;
%                 dQdV_raw_dch = dQ_dch ./ dV_dch;
%                 dQdV_filt_dch = movmean(dQdV_raw_dch, 21);
%                 mask_win_dch = V_u_dch <= win_dch_max & V_u_dch >= win_dch_min;
% 
%                 if sum(mask_win_dch) > 0
%                     pk_height_dch = max(dQdV_filt_dch(mask_win_dch));
%                     pk_area_dch = trapz(V_u_dch(mask_win_dch), dQdV_filt_dch(mask_win_dch));
% 
%                     subplot(1,2,2);
%                     hold on;
%                     plot(V_u_dch(mask_win_dch), dQdV_filt_dch(mask_win_dch), 'r-', 'LineWidth', 2);
%                     xlabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
%                     ylabel('dQ/dV (Ah/V)', 'FontWeight', 'bold', 'FontSize', 12);
%                     title(sprintf('Discharge 피크 (%s, %s)\n높이: %.1f Ah/V, 면적: %.2f Ah', ch, cyc, pk_height_dch, abs(pk_area_dch)), ...
%                         'FontSize', 12, 'FontWeight', 'bold');
%                     grid on;
%                     hold off;
%                 end
%             end
% 
%             sgtitle(sprintf('Phase 2-2: dQ/dV 피크 특성 추출 (%s, %s, 0.5C)', ch, cyc), 'FontSize', 14, 'FontWeight', 'bold');
%             saveas(fig8, fullfile(saveDir, sprintf('Phase2_2_Peak_Detection_%s_%s.fig', ch, cyc)));
%             close(fig8);
%         else
%             fprintf('  Skipping Peak Detection for %s %s (No data)\n', ch, cyc);
%             close(fig8);
%         end
%     end
% end
% fprintf('Phase 2 완료.\n');

%% ========================================================================
% Phase 3: 라벨 생성 시각화 (All Channels Overlay)
%% ========================================================================
fprintf('\n=== Phase 3: 라벨 생성 시각화 (All Channels Overlay) ===\n');

% FeatureTable에서 라벨 추출
fig9 = figure('Position', [50, 50, 1600, 500], 'Name', 'Phase 3: Label Trends (All Channels)');

% Subplot 설정
subplot(1,3,1); hold on;
title('SOH 열화 트렌드', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Cycle Number', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('SOH - Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
grid on;

subplot(1,3,2); hold on;
title('LLI (리튬 손실) 트렌드', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Cycle Number', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('LLI (%)', 'FontWeight', 'bold', 'FontSize', 12);
grid on;

subplot(1,3,3); hold on;
title('LAM (활물질 손실) 트렌드', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Cycle Number', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('LAM (%)', 'FontWeight', 'bold', 'FontSize', 12);
grid on;

% 각 채널별로 Plot
ch_colors = lines(length(channels));

for ch_i = 1:length(channels)
    ch = channels{ch_i};
    
    % 해당 채널의 0.5C 데이터 추출
    ch_idx = strcmp(FeatureTable.CellID, ch) & FeatureTable.CrateNum == 0.5;
    if sum(ch_idx) == 0, continue; end
    
    ch_data = FeatureTable(ch_idx, :);
    cycles = ch_data.Cycle;
    SOH = ch_data.Y_Labels(:,1);
    LLI = ch_data.Y_Labels(:,2);
    LAM = ch_data.Y_Labels(:,3);
    
    subplot(1,3,1);
    plot(cycles, SOH, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', ch_colors(ch_i,:), 'DisplayName', ch);
    
    subplot(1,3,2);
    plot(cycles, LLI, 's-', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', ch_colors(ch_i,:), 'DisplayName', ch);
    
    subplot(1,3,3);
    plot(cycles, LAM, '^-', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', ch_colors(ch_i,:), 'DisplayName', ch);
end

% 범례 추가
subplot(1,3,1); legend('Location', 'best', 'FontSize', 8, 'NumColumns', 2);
subplot(1,3,2); legend('Location', 'best', 'FontSize', 8, 'NumColumns', 2);
subplot(1,3,3); legend('Location', 'best', 'FontSize', 8, 'NumColumns', 2);

sgtitle('Phase 3: 라벨 (SOH, LLI, LAM) 사이클별 트렌드 (All Channels)', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig9, fullfile(saveDir, 'Phase3_Label_Trends_All_Channels.fig'));
close(fig9);

fprintf('Phase 3 완료.\n');

%% ========================================================================
% Phase 4: 정규화 시각화
% ========================================================================
fprintf('\n=== Phase 4: 정규화 시각화 ===\n');

fig10 = figure('Position', [50, 50, 1600, 900], 'Name', 'Phase 4: Normalization');

% 0.5C 데이터만
crate_idx = FeatureTable.CrateNum == 0.5;
X_raw = FeatureTable.X_Features(crate_idx, :);
X_norm = FeatureTable.X_Normalized(crate_idx, :);

% 정규화 전
subplot(2,1,1);
boxplot(X_raw, 'Labels', {'Seg1','Seg2','Seg3','Seg4','Seg5','Seg1','Seg2','Seg3','Seg4','Seg5','PkH','PkA','PkH','PkA'});
ylabel('Feature Value (Raw)', 'FontWeight', 'bold', 'FontSize', 12);
title('정규화 전 피처 분포 (0.5C)', 'FontSize', 13, 'FontWeight', 'bold');
grid on;

% 정규화 후
subplot(2,1,2);
boxplot(X_norm, 'Labels', {'Seg1','Seg2','Seg3','Seg4','Seg5','Seg1','Seg2','Seg3','Seg4','Seg5','PkH','PkA','PkH','PkA'});
ylabel('Feature Value (Normalized)', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Features (Chg dQ 1-5, Dch dQ 1-5, Chg dQdV 2, Dch dQdV 2)', 'FontWeight', 'bold', 'FontSize', 11);
title('정규화 후 피처 분포 (0.5C, Z-score)', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
ylim([-4 4]);

sgtitle('Phase 4: C-rate별 독립 정규화 (0.5C 예시)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig10, fullfile(saveDir, 'Phase4_Normalization.fig'));
close(fig10);

fprintf('Phase 4 완료.\n');

%% ========================================================================
% Box 3: 방전 V-t 곡선 (Static, 8채널) — 시간기반 균등 분할
% ========================================================================
fprintf('\n=== Box 3: 방전 V-t 곡선 — Static Capacity (시간 균등 분할) ===\n');

% 8채널 색상
ch_colors = lines(length(channels));

% 8채널 Static 데이터 수집 (cyc0)
all_V_full = cell(length(channels), 1);
all_t_full = cell(length(channels), 1);
all_V_win = cell(length(channels), 1);
all_t_win = cell(length(channels), 1);

for i = 1:length(channels)
    ch = channels{i};
    if ~isfield(RPT_VQ_grid.cyc0.(ch), 'Static'), continue; end
    data_s = RPT_VQ_grid.cyc0.(ch).Static;
    V_s = data_s.V_grid;
    t_s = data_s.t;
    
    % 시간 오름차순 정렬 (Static V_grid는 전압 오름차순, 시간은 역순)
    [t_s, si] = sort(t_s, 'ascend');
    V_s = V_s(si);
    
    all_V_full{i} = V_s;
    all_t_full{i} = t_s;
    
    % 윈도우 내 추출
    mask_w = V_s <= win_dch_max & V_s >= win_dch_min;
    all_V_win{i} = V_s(mask_w);
    all_t_win{i} = t_s(mask_w);
end

% 평균 V-t 관계 (윈도우 내, 공통 전압 그리드에 시간 보간 후 평균)
V_common = RPT_VQ_grid.cyc0.(channels{1}).Static.V_grid;
mask_common = V_common <= win_dch_max & V_common >= win_dch_min;
V_common_win = V_common(mask_common);

T_interp_all = zeros(length(channels), length(V_common_win));
valid_cnt = 0;
for i = 1:length(channels)
    if isempty(all_V_win{i}), continue; end
    valid_cnt = valid_cnt + 1;
    T_interp_all(valid_cnt, :) = interp1(all_V_win{i}, all_t_win{i}, V_common_win, 'linear', NaN)';
end
T_interp_all = T_interp_all(1:valid_cnt, :);
avg_T_win = mean(T_interp_all, 1, 'omitnan');

% 시간 오름차순 정렬 (V_common_win은 오름차순이므로 t는 내림차순)
[avg_T_win_sorted, si_avg] = sort(avg_T_win, 'ascend');
V_avg_sorted = V_common_win(si_avg);

% 시간 5등분 경계
T_total = avg_T_win_sorted(end) - avg_T_win_sorted(1);
t_divisions = linspace(avg_T_win_sorted(1), avg_T_win_sorted(end), num_segments + 1);

% 각 분할점에서의 전압값
V_at_divisions = interp1(avg_T_win_sorted, V_avg_sorted, t_divisions, 'linear');

% 세그먼트 색상
seg_colors = [0.90 0.30 0.30;   % Seg1: 빨강
              0.90 0.60 0.15;   % Seg2: 주황
              0.20 0.70 0.30;   % Seg3: 초록
              0.20 0.50 0.85;   % Seg4: 파랑
              0.60 0.30 0.80];  % Seg5: 보라

fig_box3 = figure('Position', [50, 50, 1300, 700], 'Name', 'Box 3: Static Discharge V-t');
hold on;

% 세그먼트별 색상 밴드
for seg = 1:num_segments
    t1 = t_divisions(seg);
    t2 = t_divisions(seg+1);
    fill([t1, t2, t2, t1], [win_dch_min, win_dch_min, win_dch_max, win_dch_max], ...
        seg_colors(seg, :), 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end

% 8채널 개별 V-t 곡선 (연한 색, 가는 선)
for i = 1:length(channels)
    if isempty(all_t_win{i}), continue; end
    plot(all_t_win{i}, all_V_win{i}, '-', 'Color', [ch_colors(i,:), 0.35], ...
        'LineWidth', 1.0, 'DisplayName', channels{i});
end

% 평균 V-t 곡선 (굵은 검정)
plot(avg_T_win_sorted, V_avg_sorted, 'k-', 'LineWidth', 3, ...
    'DisplayName', 'Average (8 ch)');

% 5등분 수직 점선
for i = 1:length(t_divisions)
    xline(t_divisions(i), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    if i <= num_segments
        t_mid = (t_divisions(i) + t_divisions(i+1)) / 2;
        text(t_mid, win_dch_max + 0.005, sprintf('Seg %d', i), ...
            'HorizontalAlignment', 'center', 'FontSize', 13, ...
            'FontWeight', 'bold', 'Color', seg_colors(i, :));
    end
end

% T_total 표시
text(mean(avg_T_win_sorted), win_dch_min - 0.012, ...
    sprintf('T_{window} = %.1f s → 5등분', T_total), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, ...
    'FontWeight', 'bold', 'BackgroundColor', [1 1 0.85]);

xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Voltage (V)', 'FontSize', 14, 'FontWeight', 'bold');
title({'Box 3: Static Capacity 방전 V-t — 시간 기반 균등 분할'; ...
       '(8채널 cyc0 평균, 분석 윈도우 내)'}, ...
    'FontSize', 15, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 9, 'NumColumns', 3);
set(gca, 'FontSize', 12);
ylim([win_dch_min - 0.02, win_dch_max + 0.015]);
grid on;
hold off;

saveas(fig_box3, fullfile(saveDir, 'Box3_Vt_Time_Division.fig'));
close(fig_box3);
fprintf('  Box 3 저장 완료.\n');

%% ========================================================================
% Box 4: 전압 매핑 과정 (Time → Curve → Voltage)
% ========================================================================
fprintf('\n=== Box 4: 전압 매핑 과정 (Mapping Concept) ===\n');

fig_box4 = figure('Position', [50, 50, 1300, 750], 'Name', 'Box 4: Voltage Mapping');
hold on;

% 세그먼트별 색상 밴드
for seg = 1:num_segments
    t1 = t_divisions(seg);
    t2 = t_divisions(seg+1);
    fill([t1, t2, t2, t1], [win_dch_min, win_dch_min, win_dch_max, win_dch_max], ...
        seg_colors(seg, :), 'FaceAlpha', 0.10, 'EdgeColor', 'none');
end

% 평균 V-t 곡선 (굵은 검정)
plot(avg_T_win_sorted, V_avg_sorted, 'k-', 'LineWidth', 3);

% 각 분할점에서 ㄱ자 화살표: t축 → 곡선 위 → V축
t_left = avg_T_win_sorted(1) - T_total*0.02;
for i = 1:length(t_divisions)
    t_pt = t_divisions(i);
    V_pt = V_at_divisions(i);
    
    % 마커 색상: 양 끝=검정, 내부=빨강
    if i == 1 || i == length(t_divisions)
        mk_color = [0.1 0.1 0.1];
    else
        mk_color = [0.85 0.15 0.15];
    end
    
    % 수직 점선: t축 → 곡선 위 점
    plot([t_pt, t_pt], [win_dch_min, V_pt], ':', 'Color', mk_color, 'LineWidth', 1.8);
    
    % 수평 점선: 곡선 위 점 → V축 방향
    plot([t_left, t_pt], [V_pt, V_pt], ':', 'Color', mk_color, 'LineWidth', 1.8);
    
    % 곡선 위 교차점 마커
    plot(t_pt, V_pt, 'o', 'MarkerSize', 10, 'MarkerFaceColor', mk_color, ...
        'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
    
    % V축 옆에 전압값 표시
    text(t_left - T_total*0.01, V_pt, ...
        sprintf('%.3f V', V_pt), ...
        'HorizontalAlignment', 'right', 'FontSize', 10, ...
        'FontWeight', 'bold', 'Color', mk_color, ...
        'BackgroundColor', [1 1 1]);
    
    % t축 아래에 시간 라벨
    text(t_pt, win_dch_min - 0.008, sprintf('t_%d', i-1), ...
        'HorizontalAlignment', 'center', 'FontSize', 11, ...
        'FontWeight', 'bold', 'Color', mk_color);
end

% 알고리즘 설명 박스
text(avg_T_win_sorted(end) - T_total*0.05, win_dch_min + 0.015, ...
    {'알고리즘:', ...
     '  1. 윈도우 내 시간을 5등분', ...
     '  2. 각 t 지점에서 V(t) 읽기', ...
     '  3. 읽은 V가 구간 경계 V_{grid}'}, ...
    'HorizontalAlignment', 'right', 'FontSize', 10, ...
    'FontWeight', 'bold', 'BackgroundColor', [1 1 0.9], ...
    'EdgeColor', [0.5 0.5 0.5], 'Margin', 6);

xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Voltage (V)', 'FontSize', 14, 'FontWeight', 'bold');
title({'Box 4: 전압 경계 매핑 — Time 균등분할 → V_{grid} 결정'; ...
       '(Static Capacity 방전, 8채널 평균)'}, ...
    'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12);
ylim([win_dch_min - 0.025, win_dch_max + 0.015]);
xlim([t_left - T_total*0.08, avg_T_win_sorted(end) + T_total*0.05]);
grid on;
hold off;

saveas(fig_box4, fullfile(saveDir, 'Box4_Voltage_Mapping.fig'));
close(fig_box4);
fprintf('  Box 4 저장 완료.\n');

%% ========================================================================
% 완료
% ========================================================================
fprintf('\n=== 모든 Phase 시각화 완료 ===\n');
fprintf('저장 위치: %s\n', saveDir);
fprintf('생성된 파일:\n');
fprintf('  - Phase 0: 4개 (데이터소스)\n');
fprintf('  - Phase 1: 2개 (Master Ruler)\n');
fprintf('  - Phase 2: 2개 (피처 추출)\n');
fprintf('  - Phase 3: 1개 (라벨)\n');
fprintf('  - Phase 4: 1개 (정규화)\n');
fprintf('  - Box 3:  1개 (V-t 시간 균등 분할)\n');
fprintf('  - Box 4:  1개 (전압 매핑 과정)\n');
