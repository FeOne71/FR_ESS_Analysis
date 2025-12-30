%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 파일명: Integrated_Degradation_Analysis_v5.m
% 기능:
%   - Part 1: RPT OCV 데이터를 이용한 dQ/dV 열화 메커니즘 분석 (LLI, LAM)
%   - Part 2: 주행부하 데이터의 충전 이벤트를 이용한 동적 열화 지표 분석 (R, dV/dt)
%   - 수정사항 1: SOC 계산 함수를 V -> SOC (SOC = f(OCV))로 올바르게 재정의
%   - 수정사항 2: 후기 휴지기 탐색 로직을 '마지막 유효 휴지기'로 수정
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
warning off;

%% ========================================================================
%  1. 기본 설정 (사용자 수정 영역)
% =========================================================================

% --- 저장 경로 ---
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Model_Data';

% --- Part 1: dQ/dV 분석용 데이터 경로 ---
rptDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated';

% --- Part 2: 주행부하 분석용 데이터 경로 ---
driveCycleDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';

% --- Part 2: SOC 계산 파라미터 ---
rest_current_threshold_A = 0.1; % 휴지기로 판단할 전류 임계값 [A]
min_rest_duration_sec = 60;   % 유효한 휴지기로 인정할 최소 지속 시간 [초]

% --- Part 2: 이벤트 검출 파라미터 ---
idle_thr_A_event = 5;       % 이벤트 구분을 위한 휴지 상태 전류 임계값 [A]
min_duration_sec_event = 10;% 이벤트로 인정할 최소 지속 시간 [초]
min_delta_I_A_event = (64 * 0.1); % 이벤트로 인정할 최소 전류 변화량 (셀 단위)

% --- Part 2: 열화 지표 추출 필터링 조건 ---
soc_range_filter = [0, 100];    % 분석할 SOC 범위 [%]
min_event_duration_for_analysis = 25; % 저항 분석을 위한 최소 이벤트 길이 [초]

% --- Part 2: 이상치(Outlier) 제거를 위한 유효성 검사 범위 ---
valid_range.R_CL = [0, 5];       % 단위: mOhm
valid_range.R_LLI = [-2, 2];     % 단위: mOhm
valid_range.R_LAM = [-2, 2];     % 단위: mOhm
valid_range.dVdt_0_1s = [0, 100];  % 단위: mV/s
valid_range.dVdt_1_5s = [0, 20];   % 단위: mV/s
valid_range.dVdt_5_20s = [0, 10];  % 단위: mV/s
% =========================================================================

if ~exist(saveDir, 'dir'); mkdir(saveDir); end
fprintf('모든 결과는 다음 폴더에 저장됩니다:\n%s\n\n', saveDir);


%% ========================================================================
%  Part 1: dQ/dV 기반 열화 메커니즘 분석
% =========================================================================
fprintf('################### Part 1: dQ/dV 분석 시작 ###################\n');

% --- 1a. RPT 데이터 로딩 ---
rptMatFile = fullfile(rptDataPath, 'OCV_integrated.mat');
if ~exist(rptMatFile, 'file'), error('RPT 데이터 파일이 없습니다: %s', rptMatFile); end
load(rptMatFile, 'OCV_data');
fprintf('OCV_integrated.mat 로드 완료.\n');

% --- 1b. RPT 데이터 및 ★★★ SOC = f(OCV) 함수 생성 ★★★ ---
all_fields = fieldnames(OCV_data);
q_grid_fields = all_fields(startsWith(all_fields, 'q_grid_rpt'));
cycle_keys_str = cellfun(@(s) s(11:end), q_grid_fields, 'UniformOutput', false);
[~, sort_idx] = sort(cellfun(@str2double, cycle_keys_str));
cycle_keys = cycle_keys_str(sort_idx);
fprintf('분석 대상 RPT 사이클: %s\n', strjoin(cycle_keys, ', '));

Q_data = struct(); V_data = struct(); Cap_data = struct(); soc_functions = struct();
for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    Q_data.(['c' key]) = OCV_data.(['q_grid_rpt' key]);
    V_data.(['c' key]) = OCV_data.(['avg_ocv_rpt' key]);
    Cap_data.(['c' key]) = OCV_data.(['mean_capacity_rpt' key]);
    
    % [수정] SOC = f(OCV) 함수를 올바르게 정의
    soc_grid = OCV_data.soc_grid;
    avg_ocv = OCV_data.(['avg_ocv_rpt' key]);
    
    % interp1을 위해 OCV 벡터가 단조 증가하도록 만듦
    [unique_ocv, idx] = unique(avg_ocv);
    unique_soc = soc_grid(idx);
    
    soc_functions.(['c' key]) = @(v_query) interp1(unique_ocv, unique_soc, v_query, 'linear', 'extrap');
end

% --- 1c. dQ/dV 계산 및 피크 정량화 ---
V_PEAK_SEARCH_MIN = 3.3; V_PEAK_SEARCH_MAX = 3.8; MIN_PEAK_PROMINENCE = 0.001;
dQdV_results = struct(); V_mid_points = struct();
for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    [dQdV_results.(['c' key]), V_mid_points.(['c' key])] = calculate_dQdV(Q_data.(['c' key]), V_data.(['c' key]));
end

V_peaks = struct(); dQdV_peaks = struct(); Peak_Tables = struct(); LLI_LAM_results = struct();
base_key = cycle_keys{1};
[V_peaks.(['c' base_key]), dQdV_peaks.(['c' base_key]), Peak_Tables.(['c' base_key])] = ...
    find_main_peak(V_mid_points.(['c' base_key]), dQdV_results.(['c' base_key]), 0, V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);

for i = 2:length(cycle_keys)
    prev_key = cycle_keys{i-1}; curr_key = cycle_keys{i};
    [V_peaks.(['c' curr_key]), dQdV_peaks.(['c' curr_key]), Peak_Tables.(['c' curr_key])] = ...
        find_main_peak(V_mid_points.(['c' curr_key]), dQdV_results.(['c' curr_key]), V_peaks.(['c' prev_key]), V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);
    [lli_V, lam_rate] = quantify_lli_lam(V_peaks.(['c' base_key]), dQdV_peaks.(['c' base_key]), V_peaks.(['c' curr_key]), dQdV_peaks.(['c' curr_key]));
    LLI_LAM_results.(['c' curr_key]) = struct('LLI_V', lli_V, 'LAM_rate', lam_rate);
end
fprintf('dQ/dV 계산 및 피크 정량화 완료.\n');

% --- 1d. 디버깅: 피크 테이블 출력 ---
fprintf('\n--- [디버깅] dQ/dV 피크 식별 결과 ---\n');
for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    fprintf('\n--- %s cyc Primary Peak ---\n', key);
    if ~isempty(Peak_Tables.(['c' key])); disp(Peak_Tables.(['c' key])); else; fprintf('피크를 찾지 못했습니다.\n'); end
end

% --- 1e. 시각화 및 저장 ---
fig1 = figure('Name', 'dQdV Degradation Analysis', 'Position', [100 100 1200 700]);
hold on; grid on; colors = lines(length(cycle_keys));
for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    plot(V_mid_points.(['c' key]), dQdV_results.(['c' key]) * 1000, 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', sprintf('%s cyc', key));
end
for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    plot_peak(V_peaks.(['c' key]), dQdV_peaks.(['c' key]), 'o', colors(i,:), sprintf('Peak @ %.3fV', V_peaks.(['c' key])));
end
title('dQ/dV Curve Comparison: Degradation Progression'); xlabel('Voltage [V]'); ylabel('dQ/dV [mAh/V]');
legend('Location', 'northeast'); xlim([V_PEAK_SEARCH_MIN-0.1 V_PEAK_SEARCH_MAX+0.1]);
summary_str = {'\bfDegradation Analysis (vs. 0 cyc)'}; summary_str{end+1} = ' ';
for i = 2:length(cycle_keys)
    curr_key = cycle_keys{i}; res = LLI_LAM_results.(['c' curr_key]);
    summary_str{end+1} = sprintf('\\bf\\color[rgb]{0,0.4,0.8}--- [%s cyc] ---', curr_key);
    summary_str{end+1} = sprintf('LLI (Voltage Shift): %.4f V', res.LLI_V);
    summary_str{end+1} = sprintf('LAM (Peak-Height Loss): %.2f %%', res.LAM_rate * 100);
end
annotation('textbox', [0.15, 0.6, 0.3, 0.3], 'String', summary_str, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k');
hold off;
savefig(fig1, fullfile(saveDir, 'dQdV_Degradation_Analysis.fig'));
fprintf('\nPart 1: dQ/dV 분석 그래프 저장 완료.\n');
close(fig1);

%% ========================================================================
%  Part 2: 동적 지표 기반 열화 추세 분석
% =========================================================================
fprintf('\n\n################### Part 2: 동적 지표 분석 시작 ###################\n');

driveCycleFiles = dir(fullfile(driveCycleDataDir, 'parsedDriveCycle_*_filtered.mat'));
fprintf('%d개의 사이클 단위 주행부하 파일을 찾았습니다.\n', length(driveCycleFiles));

feature_data = table();
for f = 1:length(driveCycleFiles)
    matFilePath = fullfile(driveCycleDataDir, driveCycleFiles(f).name);
    fprintf('\n--- 처리 중인 파일: %s ---\n', driveCycleFiles(f).name);
    cycle_str_token = regexp(driveCycleFiles(f).name, '(\d+cyc)', 'tokens');
    if isempty(cycle_str_token), continue; end
    cycle_key = strrep(cycle_str_token{1}{1}, 'cyc', '');
    if ~isfield(soc_functions, ['c' cycle_key])
        fprintf('경고: %s 사이클에 해당하는 SOC 함수를 찾을 수 없습니다. 건너뜁니다.\n', cycle_key); continue;
    end
    current_soc_func = soc_functions.(['c' cycle_key]);
    data = load(matFilePath); data_struct_name = fieldnames(data); drive_data = data.(data_struct_name{1});
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
                fprintf('  > 처리 중: %s - %s - %s\n', channel_names{ch_idx}, soc_level, profile_name);
                V_profile = profile_data.V; I_profile = profile_data.I; t_profile = profile_data.t;
                [SOC_profile, soc_calc_success] = calculate_soc_two_point(...
                    V_profile, I_profile, t_profile, current_soc_func, rest_current_threshold_A, min_rest_duration_sec);
                if ~soc_calc_success, continue; end
                [chg_indices, ~] = find_events_by_state_transition(I_profile, idle_thr_A_event, min_duration_sec_event, min_delta_I_A_event);
                for k = 1:size(chg_indices, 1)
                    s = chg_indices(k, 1); e = chg_indices(k, 2);
                    evt = create_event_struct('charge', s, e, t_profile, I_profile, V_profile, SOC_profile, str2double(cycle_key));
                    if isempty(evt), continue; end
                    avg_soc = mean(evt.soc_seq_pct); duration = evt.time_seq_s(end) - evt.time_seq_s(1);
                    if ~(avg_soc >= soc_range_filter(1) && avg_soc <= soc_range_filter(2) && duration >= min_event_duration_for_analysis), continue; end
                    [features, is_valid] = extract_dynamic_features(evt, valid_range);
                    if is_valid
                        temp_table = table(evt.cycle, avg_soc, features.R_CL, features.R_LLI, features.R_LAM, ...
                            features.dVdt_0_1s, features.dVdt_1_5s, features.dVdt_5_20s, ...
                            'VariableNames', {'Cycle', 'SOC', 'R_CL', 'R_LLI', 'R_LAM', 'dVdt_0_1s', 'dVdt_1_5s', 'dVdt_5_20s'});
                        feature_data = [feature_data; temp_table];
                    end
                end
            end
        end
    end
end
fprintf('\n총 %d개의 유효 이벤트에서 동적 열화 지표를 추출했습니다.\n', height(feature_data));

fprintf('\n--- [디버깅] 사이클별 열화 지표 평균 ---\n');
if ~isempty(feature_data), disp(groupsummary(feature_data, 'Cycle', 'mean')); else, fprintf('분석할 유효 이벤트가 없습니다.\n'); return; end

ma_window = min(50, height(feature_data));
fig2 = figure('Name', 'Resistance Feature Trends', 'Position', [100 100 1200 800]); tiledlayout(3,1, 'Padding', 'compact');
resistance_features = {'R_CL', 'R_LLI', 'R_LAM'}; y_labels_R = {'R_{CL} (m\Omega)', 'R_{LLI} (m\Omega)', 'R_{LAM} (m\Omega)'};
for i = 1:3
    nexttile; feature_name = resistance_features{i};
    scatter(feature_data.Cycle, feature_data.(feature_name), 30, feature_data.SOC, 'filled', 'MarkerFaceAlpha', 0.5);
    hold on; grid on; box on; if ma_window > 1, plot(feature_data.Cycle, movmean(feature_data.(feature_name), ma_window), 'r-', 'LineWidth', 2.5); end
    title(strrep(feature_name, '_', '-')); ylabel(y_labels_R{i}); xlabel('Cycle'); h = colorbar; ylabel(h, 'SOC (%)');
end
sgtitle('Dynamic Resistance Trends over Cycles', 'FontSize', 16);
savefig(fig2, fullfile(saveDir, 'Resistance_Trends.fig')); fprintf('\nPart 2: 저항 트렌드 그래프 저장 완료.\n'); close(fig2);

fig3 = figure('Name', 'dVdt Feature Trends', 'Position', [150 150 1200 800]); tiledlayout(3,1, 'Padding', 'compact');
dvdt_features = {'dVdt_0_1s', 'dVdt_1_5s', 'dVdt_5_20s'}; y_labels_dVdt = {'dV/dt_{0-1s} (mV/s)', 'dV/dt_{1-5s} (mV/s)', 'dV/dt_{5-20s} (mV/s)'};
for i = 1:3
    nexttile; feature_name = dvdt_features{i};
    scatter(feature_data.Cycle, feature_data.(feature_name), 30, feature_data.SOC, 'filled', 'MarkerFaceAlpha', 0.5);
    hold on; grid on; box on; if ma_window > 1, plot(feature_data.Cycle, movmean(feature_data.(feature_name), ma_window), 'b-', 'LineWidth', 2.5); end
    title(strrep(feature_name, '_', '-')); ylabel(y_labels_dVdt{i}); xlabel('Cycle'); h = colorbar; ylabel(h, 'SOC (%)');
end
sgtitle('Voltage Response (dV/dt) Trends over Cycles', 'FontSize', 16);
savefig(fig3, fullfile(saveDir, 'dVdt_Trends.fig')); fprintf('Part 2: dV/dt 트렌드 그래프 저장 완료.\n'); close(fig3);

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
    if isempty(V_mid(valid_mask)), V_peak_main=NaN; dQdV_peak_main=NaN; Peak_Table=table(); return; end
    [pks, locs, ~, prom] = findpeaks(dQdV(valid_mask), V_mid(valid_mask), 'MinPeakProminence', MIN_PROM);
    if isempty(pks), V_peak_main=NaN; dQdV_peak_main=NaN; Peak_Table=table(); return; end
    Peak_Table = table(locs', pks', prom', 'VariableNames', {'V_Peak', 'dQdV_Peak', 'Prominence'});
    Peak_Table = sortrows(Peak_Table, 'Prominence', 'descend');
    if V_peak_ref > 0, [~, idx] = min(abs(Peak_Table.V_Peak - V_peak_ref)); else, idx = 1; end
    V_peak_main = Peak_Table.V_Peak(idx); dQdV_peak_main = Peak_Table.dQdV_Peak(idx);
end

function [LLI_shift_V, LAM_loss_rate] = quantify_lli_lam(V_peak_base, dQdV_peak_base, V_peak_deg, dQdV_peak_deg)
    LLI_shift_V = V_peak_base - V_peak_deg;
    if dQdV_peak_base > 1e-6, LAM_loss_rate = (dQdV_peak_base - dQdV_peak_deg) / dQdV_peak_base; else, LAM_loss_rate = NaN; end
end

function plot_peak(V_peak, dQdV_peak, marker, color, label_text)
    if isfinite(V_peak) && isfinite(dQdV_peak)
        plot(V_peak, dQdV_peak * 1000, marker, 'MarkerSize', 10, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k');
        text(V_peak, dQdV_peak * 1000 * 1.05, label_text, 'FontSize', 9, 'Color', 'k', 'HorizontalAlignment', 'center');
    end
end

function [chg_events, dch_events] = find_events_by_state_transition(I_profile, idle_thr, min_duration, min_delta_I)
    chg_events = []; dch_events = []; states = zeros(size(I_profile));
    states(I_profile >= idle_thr) = 1; states(I_profile <= -idle_thr) = -1;
    start_points = find(diff(states) ~= 0 & states(1:end-1) == 0) + 1;
    end_points = find(diff(states) ~= 0 & states(2:end) == 0);
    current_event_start = 1;
    for i = 1:length(start_points)
        s_actual = start_points(i); if s_actual < current_event_start, continue; end
        possible_ends = end_points(end_points > s_actual); if isempty(possible_ends), continue; end
        e_actual = possible_ends(1); s_defined = s_actual - 1; e_defined = e_actual;
        if (e_defined - s_defined + 1) >= min_duration
            max_abs_I = max(abs(I_profile(s_defined:e_defined)));
            if states(s_actual) == 1 && max_abs_I > min_delta_I
                chg_events = [chg_events; s_defined, e_defined]; current_event_start = e_defined + 1;
            elseif states(s_actual) == -1 && max_abs_I > min_delta_I
                dch_events = [dch_events; s_defined, e_defined]; current_event_start = e_defined + 1;
            end
        end
    end
end

function [soc_profile, success] = calculate_soc_two_point(V, I, t, soc_func, rest_thr, min_rest_dur)
    soc_profile = []; success = false;
    if isduration(t), t = seconds(t); end
    t = t(:); t_relative = t - t(1);

    is_rest = abs(I) < rest_thr;
    rest_blocks = bwconncomp(is_rest);
    
    initial_rest_indices = [];
    final_rest_indices = [];

    % 1. 초기 휴지기 찾기
    if ~isempty(rest_blocks.PixelIdxList) && rest_blocks.PixelIdxList{1}(1) <= 5
        indices = rest_blocks.PixelIdxList{1};
        duration = t_relative(indices(end)) - t_relative(indices(1));
        if duration >= min_rest_dur, initial_rest_indices = indices; end
    end
    if isempty(initial_rest_indices), fprintf('      > Debug: 유효한 초기 휴지기 탐색 실패.\n'); return; end
    
    % 2. ★★★ 후기 휴지기 찾기 (데이터의 '마지막' 유효 휴지기) ★★★
    if rest_blocks.NumObjects > 1
        % 마지막 블록부터 거꾸로 탐색
        for i = rest_blocks.NumObjects:-1:1
            indices = rest_blocks.PixelIdxList{i};
            % 초기 휴지기 블록이 아니어야 함
            if ~isequal(indices, initial_rest_indices)
                duration = t_relative(indices(end)) - t_relative(indices(1));
                if duration >= min_rest_dur
                    final_rest_indices = indices;
                    break; % 마지막 유효 휴지기를 찾으면 탐색 종료
                end
            end
        end
    end
    if isempty(final_rest_indices), fprintf('      > Debug: 유효한 후기 휴지기 탐색 실패.\n'); return; end
    
    % 3. 기준점 t1, SOC1 / t2, SOC2 설정
    idx1 = initial_rest_indices(end);
    v1_start_idx = find(t_relative >= (t_relative(idx1) - 10), 1);
    V1 = mean(V(v1_start_idx:idx1));
    SOC1 = soc_func(V1);
    
    idx2 = final_rest_indices(end);
    v2_start_idx = find(t_relative >= (t_relative(idx2) - 10), 1);
    V2 = mean(V(v2_start_idx:idx2));
    SOC2 = soc_func(V2);

    if isnan(SOC1) || isnan(SOC2)
        fprintf('      > Debug: OCV-SOC 변환 실패. V1=%.3f, V2=%.3f\n', V1, V2); return;
    end
    
    % 4. 전하량 계산 및 SOC 보간
    delta_Q_total = trapz(t(idx1:idx2), I(idx1:idx2)) / 3600;
    if abs(delta_Q_total) < 1e-3, delta_Q_total = 1e-3; end

    soc_profile = NaN(size(t)); soc_profile(1:idx1) = SOC1; soc_profile(idx2:end) = SOC2;
    for i = (idx1+1):(idx2-1)
        delta_Q_partial = trapz(t(idx1:i), I(idx1:i)) / 3600;
        soc_profile(i) = SOC1 + (SOC2 - SOC1) * (delta_Q_partial / delta_Q_total);
    end
    soc_profile = max(0, min(100, soc_profile));
    success = true;
    fprintf('      > SOC 계산 성공 (%.1f%% -> %.1f%%)\n', SOC1, SOC2);
end

function evt = create_event_struct(type, s, e, t_s, I_A, V_V, soc_pct, cycle_num)
    if (e-s) < 1, evt = []; return; end
    evt = struct('type', type, 'cycle', cycle_num, ...
        'time_seq_s', t_s(s:e), 'current_seq_cell_A', I_A(s:e), ...
        'voltage_seq_cell_V', V_V(s:e), 'soc_seq_pct', soc_pct(s:e));
end

function [features, is_valid] = extract_dynamic_features(evt, valid_range)
    features = struct(); is_valid = false;
    t_rel = evt.time_seq_s - evt.time_seq_s(1); V_seq = evt.voltage_seq_cell_V; I_seq = evt.current_seq_cell_A;
    idx_1s = find(t_rel >= 1, 1); idx_5s = find(t_rel >= 5, 1); idx_20s = find(t_rel >= 20, 1);
    if isempty(idx_1s) || isempty(idx_5s) || isempty(idx_20s), return; end
    if abs(I_seq(1)) > 0.1 || abs(I_seq(idx_1s)) < 1 || abs(I_seq(idx_5s)) < 1 || abs(I_seq(idx_20s)) < 1, return; end
    delta_I_1s = I_seq(idx_1s) - I_seq(1); delta_I_5s = I_seq(idx_5s) - I_seq(1); delta_I_20s = I_seq(idx_20s) - I_seq(1);
    if abs(delta_I_1s) < 1 || abs(delta_I_5s) < 1 || abs(delta_I_20s) < 1, return; end
    R_1s = (V_seq(idx_1s) - V_seq(1)) / delta_I_1s * 1000;
    R_5s = (V_seq(idx_5s) - V_seq(1)) / delta_I_5s * 1000;
    R_20s = (V_seq(idx_20s) - V_seq(1)) / delta_I_20s * 1000;
    features.R_CL = R_1s; features.R_LLI = R_5s - R_1s; features.R_LAM = R_20s - R_5s;
    features.dVdt_0_1s = (V_seq(idx_1s) - V_seq(1)) / t_rel(idx_1s) * 1000;
    features.dVdt_1_5s = (V_seq(idx_5s) - V_seq(idx_1s)) / (t_rel(idx_5s) - t_rel(idx_1s)) * 1000;
    features.dVdt_5_20s = (V_seq(idx_20s) - V_seq(idx_5s)) / (t_rel(idx_20s) - t_rel(idx_5s)) * 1000;
    is_valid = (features.R_CL > valid_range.R_CL(1) && features.R_CL < valid_range.R_CL(2)) && ...
               (features.R_LLI > valid_range.R_LLI(1) && features.R_LLI < valid_range.R_LLI(2)) && ...
               (features.R_LAM > valid_range.R_LAM(1) && features.R_LAM < valid_range.R_LAM(2)) && ...
               (features.dVdt_0_1s > valid_range.dVdt_0_1s(1) && features.dVdt_0_1s < valid_range.dVdt_0_1s(2)) && ...
               (features.dVdt_1_5s > valid_range.dVdt_1_5s(1) && features.dVdt_1_5s < valid_range.dVdt_1_5s(2)) && ...
               (features.dVdt_5_20s > valid_range.dVdt_5_20s(1) && features.dVdt_5_20s < valid_range.dVdt_5_20s(2));
end