%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 파일명: Final_Correlation_Analysis_v19.m
% 기능:
%   - Part 1: RPT 데이터로 열화 지표(LLI, LAM, Capacity) 정량화
%   - Part 2: 주행부하 데이터에서 '초기 전류가 안정적인' 충전 이벤트만 필터링하여 동적 지표 추출
%   - Part 2.5: [핵심 수정] 안정적인 방식으로 통계량을 계산하여 디버깅 정보 출력
%   - Part 3, 4, 5: 상관관계, Heatmap, DC별 동향 분석
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
warning off;

%% ========================================================================
%  1. 기본 설정 (사용자 수정 영역)
% =========================================================================

% --- 저장 경로 ---
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Model_Data\Final_Correlation';

% --- 데이터 경로 ---
rptDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated';
driveCycleDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';

% --- 동적 지표 추출 파라미터 ---
idle_thr_A_event = 5;       
min_duration_sec_event = 35; 
min_delta_I_A_event = (64 * 0.1);

% 초기 전류 안정성 필터 파라미터
stable_current_duration_s = 30; 
stable_current_std_threshold_A = 0.5; 

% 이상치 제거를 위한 유효성 검사 범위
valid_range.R_1s = [0, 10];       % mOhm
valid_range.R_5s = [0, 10];       % mOhm
valid_range.R_30s = [0, 10];      % mOhm
valid_range.dVdt_1s = [0, 100];  % mV/s
valid_range.dVdt_5s = [0, 50];   % mV/s
valid_range.dVdt_20s = [0, 20];  % mV/s
% =========================================================================

if ~exist(saveDir, 'dir'); mkdir(saveDir); end
fprintf('모든 결과는 다음 폴더에 저장됩니다:\n%s\n\n', saveDir);


%% ========================================================================
%  Part 1: RPT 기반 열화 지표 정량화 (LLI, LAM, Capacity)
% =========================================================================
fprintf('################### Part 1: 열화 지표(LLI, LAM, Capacity) 정량화 시작 ###################\n');

rptMatFile = fullfile(rptDataPath, 'OCV_integrated.mat');
load(rptMatFile, 'OCV_data');
all_fields = fieldnames(OCV_data);
q_grid_fields = all_fields(startsWith(all_fields, 'q_grid_rpt'));
cycle_keys_str = cellfun(@(s) s(11:end), q_grid_fields, 'UniformOutput', false);
[~, sort_idx] = sort(cellfun(@str2double, cycle_keys_str));
cycle_keys = cycle_keys_str(sort_idx);
fprintf('분석 대상 RPT 사이클: %s\n', strjoin(cycle_keys, ', '));

V_PEAK_SEARCH_MIN = 3.3; V_PEAK_SEARCH_MAX = 3.8; MIN_PEAK_PROMINENCE = 0.001;
degradation_modes_table = table(); 

base_key = cycle_keys{1};
Q_base = OCV_data.(['q_grid_rpt' base_key]); V_base = OCV_data.(['avg_ocv_rpt' base_key]);
[dQdV_base, V_mid_base] = calculate_dQdV(Q_base, V_base);
[V_peak_base, dQdV_peak_base, ~] = find_main_peak(V_mid_base, dQdV_base, 0, V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);
Cap_base = OCV_data.(['mean_capacity_rpt' base_key]);

new_row = table(str2double(base_key), 0, 0, Cap_base, 'VariableNames', {'Cycle', 'LLI_mV', 'LAM_pct', 'Capacity_Ah'});
degradation_modes_table = [degradation_modes_table; new_row];

for i = 2:length(cycle_keys)
    curr_key = cycle_keys{i};
    Q_curr = OCV_data.(['q_grid_rpt' curr_key]); V_curr = OCV_data.(['avg_ocv_rpt' curr_key]);
    [dQdV_curr, V_mid_curr] = calculate_dQdV(Q_curr, V_curr);
    [V_peak_curr, dQdV_peak_curr, ~] = find_main_peak(V_mid_curr, dQdV_curr, V_peak_base, V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);
    Cap_curr = OCV_data.(['mean_capacity_rpt' curr_key]);
    
    if ~isnan(V_peak_curr)
        [lli_V, lam_rate] = quantify_lli_lam(V_peak_base, dQdV_peak_base, V_peak_curr, dQdV_peak_curr);
        new_row = table(str2double(curr_key), lli_V * 1000, lam_rate * 100, Cap_curr, 'VariableNames', {'Cycle', 'LLI_mV', 'LAM_pct', 'Capacity_Ah'});
        degradation_modes_table = [degradation_modes_table; new_row];
    else
        fprintf('  > 경고: %s 사이클에서 dQ/dV 피크를 찾지 못해 LLI/LAM 분석에서 제외합니다.\n', curr_key);
    end
end
fprintf('LLI, LAM, Capacity 정량화 완료.\n');


%% ========================================================================
%  Part 2: 주행부하 동적 지표 추출 (개별 이벤트)
% =========================================================================
fprintf('\n\n################### Part 2: 동적 지표 추출 시작 ###################\n');

all_events_table = table();
extracted_events_by_dc = struct(); 

for i = 1:length(cycle_keys)
    cycle_key = cycle_keys{i};
    fprintf('\n--- 처리 중인 사이클: %s ---\n', cycle_key);
    
    matFileName = fullfile(driveCycleDataDir, sprintf('parsedDriveCycle_%scyc_filtered.mat', cycle_key));
    if ~exist(matFileName, 'file'), fprintf('  > 경고: %s 파일을 찾을 수 없어 건너뜁니다.\n', matFileName); continue; end
    
    data = load(matFileName); data_struct_name = fieldnames(data); drive_data = data.(data_struct_name{1});
    
    channel_names = fieldnames(drive_data);
    for ch_idx = 1:length(channel_names)
        channel_data = drive_data.(channel_names{ch_idx});
        soc_levels = {'SOC90', 'SOC70', 'SOC50'};
        
        for soc_idx = 1:length(soc_levels)
            soc_level = soc_levels{soc_idx};
            if ~isfield(channel_data, soc_level), continue; end
            
            profile_names = fieldnames(channel_data.(soc_level));
            for p_idx = 1:length(profile_names)
                profile_name = profile_names{p_idx};
                profile_data = channel_data.(soc_level).(profile_name);
                
                I_profile = profile_data.I; V_profile = profile_data.V; t_profile = profile_data.t;
                if isduration(t_profile), t_profile = seconds(t_profile); end

                chg_indices = find_pure_events(I_profile, t_profile, idle_thr_A_event, min_duration_sec_event, min_delta_I_A_event, 1);
                
                for k = 1:size(chg_indices, 1)
                    s = chg_indices(k, 1); e = chg_indices(k, 2);
                    evt = struct('time_seq_s', t_profile(s:e), 'current_seq_A', I_profile(s:e), 'voltage_seq_V', V_profile(s:e));
                    
                    if ~is_current_stable(evt, stable_current_duration_s, stable_current_std_threshold_A), continue; end

                    [features, is_valid] = extract_dynamic_features(evt, valid_range);
                    if is_valid
                        new_row = table(str2double(cycle_key), string(profile_name), string(soc_level), ...
                            features.R_1s, features.R_5s, features.R_30s, ...
                            features.dVdt_1s, features.dVdt_5s, features.dVdt_20s, ...
                            'VariableNames', {'Cycle', 'DC_Profile', 'SOC_Level', 'R_1s', 'R_5s', 'R_30s', 'dVdt_1s', 'dVdt_5s', 'dVdt_20s'});
                        all_events_table = [all_events_table; new_row];
                        
                        cycle_field = sprintf('C%s', cycle_key);
                        if ~isfield(extracted_events_by_dc, profile_name), extracted_events_by_dc.(profile_name) = struct(); end
                        if ~isfield(extracted_events_by_dc.(profile_name), cycle_field), extracted_events_by_dc.(profile_name).(cycle_field) = {}; end
                        extracted_events_by_dc.(profile_name).(cycle_field){end+1} = evt;
                    end
                end
            end
        end
    end
end
fprintf('\n총 %d개의 유효 충전 이벤트에서 동적 지표 추출 완료.\n', height(all_events_table));
save(fullfile(saveDir, 'extracted_events.mat'), 'extracted_events_by_dc');
fprintf('추출된 이벤트 데이터가 extracted_events.mat 파일로 저장되었습니다.\n');

%% ========================================================================
%  Part 2.5: 통계 분석 결과 디버깅
% =========================================================================
fprintf('\n\n################### Part 2.5: 통계 분석 디버깅 시작 ###################\n');

if isempty(all_events_table)
    fprintf('분석할 유효 이벤트가 없어 통계 분석을 건너뜁니다.\n');
else
    % ★★★★★ [수정] groupsummary를 안정적인 방식으로 분리 실행 ★★★★★
    fprintf('\n--- [디버깅] 사이클별 동적 지표 요약 통계 ---\n');
    summary_vars = {'R_1s', 'R_5s', 'R_30s', 'dVdt_1s', 'dVdt_5s', 'dVdt_20s'};
    
    % 평균, 표준편차, 개수를 각각 계산
    mean_stats = groupsummary(all_events_table, 'Cycle', 'mean', summary_vars);
    std_stats = groupsummary(all_events_table, 'Cycle', 'std', summary_vars);
    count_stats = groupsummary(all_events_table, 'Cycle', 'numel', summary_vars{1}); % numel은 변수 하나만 지정해도 됨
    count_stats.Properties.VariableNames{end} = 'EventCount'; % 변수 이름 변경
    
    % 계산된 통계 테이블들을 'Cycle'을 기준으로 결합
    summary_stats = join(mean_stats, std_stats);
    summary_stats = join(summary_stats, count_stats(:, {'Cycle', 'EventCount'}));
    
    disp(summary_stats);

    fprintf('\n--- [디버깅] DC 프로파일 및 사이클별 R_1s 평균값 요약 ---\n');
    summary_by_dc_cycle = groupsummary(all_events_table, {'DC_Profile', 'Cycle'}, 'mean', 'R_1s');
    disp(head(summary_by_dc_cycle, 16));
end

%% ========================================================================
%  Part 3: 열화 모드와 동적 지표 상관관계 분석 (산점도)
% =========================================================================
fprintf('\n\n################### Part 3: 상관관계 산점도 분석 시작 ###################\n');

final_table = outerjoin(all_events_table, degradation_modes_table, 'Keys', 'Cycle', 'MergeKeys', true);
final_table = rmmissing(final_table);
if height(final_table) < 2, fprintf('상관분석을 위한 데이터가 부족합니다.\n'); return; end

fig_corr_R = figure('Name', 'Resistance Correlation', 'Position', [100 100, 1600, 800]);
t1 = tiledlayout(3,3, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t1, 'Correlation: Degradation Indicators vs. Dynamic Resistances', 'FontSize', 16, 'FontWeight', 'bold');
x_vars = {'LLI_mV', 'LAM_pct', 'Capacity_Ah'};
x_labels = {'LLI (Voltage Shift) [mV]', 'LAM (Peak Height Loss) [%]', 'Capacity [Ah]'};
y_vars_R = {'R_1s', 'R_5s', 'R_30s'};
y_labels_R = {'R_{1s} [m\Omega]', 'R_{5s} [m\Omega]', 'R_{30s} [m\Omega]'};
for i = 1:3, for j = 1:3, nexttile; plot_correlation(final_table.(x_vars{j}), final_table.(y_vars_R{i}), x_labels{j}, y_labels_R{i}); end, end
savefig(fig_corr_R, fullfile(saveDir, 'Correlation_Scatter_Resistance.fig'));
fprintf('저항 상관관계 산점도 그래프 저장 완료.\n');
close(fig_corr_R);

fig_corr_dVdt = figure('Name', 'dVdt Correlation', 'Position', [150 150, 1600, 800]);
t2 = tiledlayout(3,3, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t2, 'Correlation: Degradation Indicators vs. dV/dt', 'FontSize', 16, 'FontWeight', 'bold');
y_vars_dVdt = {'dVdt_1s', 'dVdt_5s', 'dVdt_20s'};
y_labels_dVdt = {'dV/dt_{1s} [mV/s]', 'dV/dt_{5s} [mV/s]', 'dV/dt_{20s} [mV/s]'};
for i = 1:3, for j = 1:3, nexttile; plot_correlation(final_table.(x_vars{j}), final_table.(y_vars_dVdt{i}), x_labels{j}, y_labels_dVdt{i}); end, end
savefig(fig_corr_dVdt, fullfile(saveDir, 'Correlation_Scatter_dVdt.fig'));
fprintf('dV/dt 상관관계 산점도 그래프 저장 완료.\n');
close(fig_corr_dVdt);


%% ========================================================================
%  Part 4: Pearson 상관계수 행렬 시각화 (Heatmap)
% =========================================================================
fprintf('\n\n################### Part 4: Pearson 상관계수 행렬 분석 시작 ###################\n');
fig_heatmap = figure('Name', 'Pearson Correlation Matrix', 'Position', [200 200, 800, 700]);
corr_table = final_table(:, {'LLI_mV', 'LAM_pct', 'Capacity_Ah', 'R_1s', 'R_5s', 'R_30s', 'dVdt_1s', 'dVdt_5s', 'dVdt_20s'});
corr_matrix = corrcoef(table2array(corr_table));
h = heatmap(corr_table.Properties.VariableNames, corr_table.Properties.VariableNames, corr_matrix, 'Colormap', jet);
h.Title = 'Pearson Correlation Matrix of All Indicators'; h.FontSize = 12;
fprintf('Pearson 상관계수 행렬 그래프 저장 완료.\n');
savefig(fig_heatmap, fullfile(saveDir, 'Correlation_Heatmap.fig'));


%% ========================================================================
%  Part 5: DC 프로파일별/사이클별 동향 분석
% =========================================================================
fprintf('\n\n################### Part 5: DC 프로파일별 동향 분석 시작 ###################\n');
fig_dc_trends = figure('Name', 'Trends by DC Profile', 'Position', [250 250, 1600, 900]);
t3 = tiledlayout(2,3, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t3, 'Dynamic Feature Trends by DC Profile over Cycles', 'FontSize', 16, 'FontWeight', 'bold');

dc_profiles = unique(all_events_table.DC_Profile);
dynamic_features = {'R_1s', 'R_5s', 'R_30s', 'dVdt_1s', 'dVdt_5s', 'dVdt_20s'};
feature_labels = {'R_{1s} [m\Omega]', 'R_{5s} [m\Omega]', 'R_{30s} [m\Omega]', 'dV/dt_{1s} [mV/s]', 'dV/dt_{5s} [mV/s]', 'dV/dt_{20s} [mV/s]'};
colors = lines(length(dc_profiles));

for i = 1:length(dynamic_features)
    nexttile; hold on; grid on;
    feature = dynamic_features{i};
    for j = 1:length(dc_profiles)
        dc_profile = dc_profiles(j);
        subset = all_events_table(all_events_table.DC_Profile == dc_profile, :);
        summary = groupsummary(subset, 'Cycle', 'mean', feature);
        plot(summary.Cycle, summary.(['mean_' feature]), 'o-', 'LineWidth', 1.5, 'Color', colors(j,:), 'DisplayName', char(dc_profile));
    end
    xlabel('Cycle'); ylabel(feature_labels{i}); title(strrep(feature, '_', ' '));
    if i == 1, legend('show', 'Location', 'best'); end
end
fprintf('DC 프로파일별 동향 분석 그래프 저장 완료.\n');
savefig(fig_dc_trends, fullfile(saveDir, 'DC_Profile_Trends.fig'));


fprintf('\n\n################### 모든 분석 완료 ###################\n');


%% ========================================================================
%                         🌟 서브 함수 (Helper Functions) 🌟
% =========================================================================

function [dQdV_AhV, V_mid] = calculate_dQdV(Q_grid_Ah, V_ocv)
    dQ = diff(Q_grid_Ah); dV = diff(V_ocv); valid_indices = abs(dV) > 1e-6; 
    dQdV_AhV = NaN(size(dV)); dQdV_AhV(valid_indices) = dQ(valid_indices) ./ dV(valid_indices);
    V_mid = V_ocv(1:end-1) + dV/2;
end

function [V_peak_main, dQdV_peak_main, Peak_Table] = find_main_peak(V_mid, dQdV, V_peak_ref, V_MIN, V_MAX, MIN_PROM)
    valid_mask = ~isnan(dQdV) & (V_mid >= V_MIN) & (V_mid <= V_MAX);
    V_peak_main=NaN; dQdV_peak_main=NaN; Peak_Table=table();
    if isempty(V_mid(valid_mask)), return; end
    [pks, locs, ~, prom] = findpeaks(dQdV(valid_mask), V_mid(valid_mask), 'MinPeakProminence', MIN_PROM);
    if isempty(pks), return; end
    Peak_Table = table(locs', pks', prom', 'VariableNames', {'V_Peak', 'dQdV_Peak', 'Prominence'});
    Peak_Table = sortrows(Peak_Table, 'Prominence', 'descend');
    if V_peak_ref > 0, [~, idx] = min(abs(Peak_Table.V_Peak - V_peak_ref)); else, idx = 1; end
    V_peak_main = Peak_Table.V_Peak(idx);
    dQdV_peak_main = Peak_Table.dQdV_Peak(idx);
end

function [LLI_shift_V, LAM_loss_rate] = quantify_lli_lam(V_peak_base, dQdV_peak_base, V_peak_deg, dQdV_peak_deg)
    LLI_shift_V = V_peak_base - V_peak_deg;
    if dQdV_peak_base > 1e-6, LAM_loss_rate = (dQdV_peak_base - dQdV_peak_deg) / dQdV_peak_base; else, LAM_loss_rate = NaN; end
end

function events = find_pure_events(I, t, idle_thr, min_dur_s, min_delta_I, target_state)
    events = [];
    states = zeros(size(I));
    states(I >= idle_thr) = 1; states(I <= -idle_thr) = -1;
    change_points = [1; find(diff(states) ~= 0) + 1; length(I)+1];

    for i = 1:length(change_points)-1
        s = change_points(i);
        e = change_points(i+1) - 1;
        if s > e, continue; end
        current_state = states(s);
        if current_state == target_state
            duration = t(e) - t(s);
            if duration >= min_dur_s && max(abs(I(s:e))) > min_delta_I
                start_idx = max(1, s-1);
                events = [events; start_idx, e];
            end
        end
    end
end

function is_stable = is_current_stable(evt, duration_s, std_threshold)
    is_stable = false;
    t_rel = evt.time_seq_s - evt.time_seq_s(1);
    start_idx = find(t_rel > 0, 1);
    if isempty(start_idx), return; end
    end_idx = find(t_rel >= duration_s, 1);
    if isempty(end_idx), return; end
    current_segment = evt.current_seq_A(start_idx:end_idx);
    if std(current_segment) < std_threshold
        is_stable = true;
    end
end

function [features, is_valid] = extract_dynamic_features(evt, valid_range)
    features = struct(); is_valid = false;
    t_rel = evt.time_seq_s - evt.time_seq_s(1); V_seq = evt.voltage_seq_V; I_seq = evt.current_seq_A;

    idx_1s = find(t_rel >= 1, 1); idx_5s = find(t_rel >= 5, 1);
    idx_20s = find(t_rel >= 20, 1); idx_30s = find(t_rel >= 30, 1);
    if isempty(idx_1s) || isempty(idx_5s) || isempty(idx_20s) || isempty(idx_30s), return; end
    
    delta_I_1s = I_seq(idx_1s) - I_seq(1); delta_I_5s = I_seq(idx_5s) - I_seq(1); delta_I_30s = I_seq(idx_30s) - I_seq(1);
    if abs(delta_I_1s) < 1 || abs(delta_I_5s) < 1 || abs(delta_I_30s) < 1, return; end

    features.R_1s = (V_seq(idx_1s) - V_seq(1)) / delta_I_1s * 1000;
    features.R_5s = (V_seq(idx_5s) - V_seq(1)) / delta_I_5s * 1000;
    features.R_30s = (V_seq(idx_30s) - V_seq(1)) / delta_I_30s * 1000;

    features.dVdt_1s = (V_seq(idx_1s) - V_seq(1)) / t_rel(idx_1s) * 1000;
    features.dVdt_5s = (V_seq(idx_5s) - V_seq(1)) / t_rel(idx_5s) * 1000;
    features.dVdt_20s = (V_seq(idx_20s) - V_seq(1)) / t_rel(idx_20s) * 1000;
    
    is_valid = (features.R_1s > valid_range.R_1s(1) && features.R_1s < valid_range.R_1s(2)) && ...
               (features.R_5s > valid_range.R_5s(1) && features.R_5s < valid_range.R_5s(2)) && ...
               (features.R_30s > valid_range.R_30s(1) && features.R_30s < valid_range.R_30s(2)) && ...
               (features.dVdt_1s > valid_range.dVdt_1s(1) && features.dVdt_1s < valid_range.dVdt_1s(2)) && ...
               (features.dVdt_5s > valid_range.dVdt_5s(1) && features.dVdt_5s < valid_range.dVdt_5s(2)) && ...
               (features.dVdt_20s > valid_range.dVdt_20s(1) && features.dVdt_20s < valid_range.dVdt_20s(2));
end

function plot_correlation(x_data, y_data, x_label, y_label)
    scatter(x_data, y_data, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
    hold on; grid on;
    p = polyfit(x_data, y_data, 1);
    x_fit = linspace(min(x_data), max(x_data), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    R = corrcoef(x_data, y_data);
    R_val = R(1,2);
    text(0.1, 0.9, sprintf('R = %.3f', R_val), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel(x_label);
    ylabel(y_label);
end