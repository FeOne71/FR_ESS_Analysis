%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FR ESS 데이터 통합 분석 스크립트 (v11 - 최종 오류 수정)
%
% 주의사항:
% 이 스크립트에서 생성한 저항은 0.1초 저항값. Not 1초 저항값
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. 기본 설정 및 경로
% =========================================================================
ocv_data_folder = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated';
ocv_data_path = fullfile(ocv_data_folder, 'OCV_integrated.mat');
drive_cycle_folder = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
save_dir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureExtraction\Results';
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

TARGET_CHANNEL_PREFIX = 'ch9_Drive_'; 
TARGET_SOC_LEVELS = {'SOC90', 'SOC70', 'SOC50'};
V_PEAK_SEARCH_MIN = 3.3; V_PEAK_SEARCH_MAX = 3.8; MIN_PEAK_PROMINENCE = 0.001;
IDLE_CURRENT_THRESHOLD = 64 * 0.05; 
MIN_EVENT_DURATION_S = 5;

%% 2. 데이터 로딩
% =========================================================================
fprintf('=== 1. 기본 데이터 로딩 중... ===\n');
if ~exist(ocv_data_path, 'file'), error('OCV 데이터 파일 경로를 확인하세요: %s', ocv_data_path); end
load(ocv_data_path); 
if ~exist('OCV_data', 'var'), error("'OCV_data' 변수를 찾을 수 없습니다."); end
drive_cycle_files = dir(fullfile(drive_cycle_folder, 'parsedDriveCycle_*_filtered.mat'));
if isempty(drive_cycle_files), error('지정된 폴더에 파싱된 주행 부하 파일이 없습니다.'); end
fprintf('%d개의 주행 부하 사이클 파일을 찾았습니다.\n데이터 로딩 완료.\n\n', length(drive_cycle_files));

%% 3. dQ/dV 분석 및 열화 모드 정량화
% =========================================================================
fprintf('=== 2. dQ/dV 분석 및 열화 모드 정량화 시작 ===\n');
all_fields = fieldnames(OCV_data);
q_grid_fields = all_fields(startsWith(all_fields, 'q_grid_rpt'));
cycle_keys_str = cellfun(@(s) s(11:end), q_grid_fields, 'UniformOutput', false);
[~, sort_idx] = sort(cellfun(@str2double, cycle_keys_str));
cycle_keys = cycle_keys_str(sort_idx);
Q_data = struct(); V_data = struct();
for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    Q_data.(['c' key]) = OCV_data.(['q_grid_rpt' key]);
    V_data.(['c' key]) = OCV_data.(['avg_ocv_rpt' key]);
end
dQdV_results = struct(); V_mid_points = struct();
for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    [dQdV_results.(['c' key]), V_mid_points.(['c' key])] = calculate_dQdV(Q_data.(['c' key]), V_data.(['c' key]));
end
V_peaks = struct(); dQdV_peaks = struct(); LLI_LAM_results = struct('LLI_V', [], 'LAM_rate', []);
base_key = cycle_keys{1};
[V_peaks.(['c' base_key]), dQdV_peaks.(['c' base_key]), ~] = find_main_peak(V_mid_points.(['c' base_key]), dQdV_results.(['c' base_key]), 0, V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);
for i = 2:length(cycle_keys)
    prev_key = cycle_keys{i-1}; curr_key = cycle_keys{i};
    [V_peaks.(['c' curr_key]), dQdV_peaks.(['c' curr_key]), ~] = find_main_peak(V_mid_points.(['c' curr_key]), dQdV_results.(['c' curr_key]), V_peaks.(['c' prev_key]), V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);
    [lli_V, lam_rate] = quantify_lli_lam(V_peaks.(['c' prev_key]), dQdV_peaks.(['c' prev_key]), V_peaks.(['c' curr_key]), dQdV_peaks.(['c' curr_key]));
    LLI_LAM_results.LLI_V(i-1) = lli_V; LLI_LAM_results.LAM_rate(i-1) = lam_rate;
end
degradation_quantified = struct('Cycles', {cycle_keys}, 'V_peaks', V_peaks, 'dQdV_peaks', dQdV_peaks, 'LLI_LAM', LLI_LAM_results);
fprintf('dQ/dV 분석 완료.\n');

%% 4. 모든 주행 부하 파일에 대해 반복 처리
% =========================================================================
fprintf('\n=== 3. 모든 주행 부하 프로파일 분석 시작 ===\n');
drive_cycle_analysis_all_cycles = struct();

for file_idx = 1:length(drive_cycle_files)
    current_filename = drive_cycle_files(file_idx).name;
    current_filepath = fullfile(drive_cycle_folder, current_filename);
    match = regexp(current_filename, 'parsedDriveCycle_(\w+)_filtered.mat', 'tokens');
    cycle_str = match{1}{1}; cycle_num_str = regexp(cycle_str, '\d+', 'match', 'once');
    fprintf('\n<<<<< %s 주행 부하 데이터 처리 중... >>>>>\n', cycle_str);
    
    loaded_data = load(current_filepath);
    expected_var_name = ['parsedDriveCycle_', cycle_str];
    if isfield(loaded_data, expected_var_name), drive_cycle_data_struct = loaded_data.(expected_var_name);
    else, fprintf('  경고: %s 파일 내에 예상 변수 ''%s''가 없습니다.\n', current_filename, expected_var_name); continue; end
    
    target_channel_name = [TARGET_CHANNEL_PREFIX, cycle_str];
    if ~isfield(drive_cycle_data_struct, target_channel_name), fprintf('  경고: %s 채널 없음.\n', target_channel_name); continue; end
    dc_data_all = drive_cycle_data_struct.(target_channel_name);
    
    ocv_field_name = ['OCV_integrated_', cycle_num_str];
    if isfield(OCV_data, ocv_field_name), ocv_ref_data = OCV_data.(ocv_field_name); fprintf('  %s OCV 데이터를 사용합니다.\n', ocv_field_name);
    else, ocv_ref_data = OCV_data.OCV_integrated_0; fprintf('  경고: %s OCV 데이터 없음. 0cyc OCV로 대체합니다.\n', ocv_field_name); end
    
    drive_cycle_events_this_cycle = struct();
    for soc_idx = 1:length(TARGET_SOC_LEVELS)
        soc_level = TARGET_SOC_LEVELS{soc_idx};
        fprintf('-- 처리 중: %s --\n', soc_level);
        drive_cycle_events_this_cycle.(soc_level) = struct();
        
        for dc_idx = 1:8
            dc_name = sprintf('DC%d', dc_idx);
            if ~isfield(dc_data_all.(soc_level), dc_name) || isempty(dc_data_all.(soc_level).(dc_name).V), continue; end
            
            dc_data = dc_data_all.(soc_level).(dc_name);
            V_measured = dc_data.V; I_measured = dc_data.I; t_measured = dc_data.t;
            
            [SOC_profile, ~] = calculate_soc_profile(V_measured, I_measured, 0.1, ocv_ref_data);
            [charge_event_indices, discharge_event_indices] = find_events_by_transition(I_measured, IDLE_CURRENT_THRESHOLD, MIN_EVENT_DURATION_S / 0.1);
            fprintf('  %s: 충전 %d개, 방전 %d개 이벤트 검출.\n', dc_name, size(charge_event_indices, 1), size(discharge_event_indices, 1));
            
            drive_cycle_events_this_cycle.(soc_level).(dc_name).charge_events = process_events(charge_event_indices, V_measured, I_measured, t_measured, SOC_profile);
            drive_cycle_events_this_cycle.(soc_level).(dc_name).discharge_events = process_events(discharge_event_indices, V_measured, I_measured, t_measured, SOC_profile);
            drive_cycle_events_this_cycle.(soc_level).(dc_name).SOC_profile_full = SOC_profile;
        end
    end
    drive_cycle_analysis_all_cycles.(['cyc' cycle_num_str]) = drive_cycle_events_this_cycle;
end
fprintf('\n모든 주행 부하 프로파일 분석 완료.\n\n');

%% 5. 최종 결과 저장
% =========================================================================
final_results = struct('degradation_quantified', degradation_quantified, 'drive_cycle_analysis', drive_cycle_analysis_all_cycles);
save_path = fullfile(save_dir, 'FR_ESS_Analysis_Results_All_Cycles.mat');
save(save_path, 'final_results', '-v7.3');
fprintf('=== 모든 분석 완료 ===\n최종 결과가 저장되었습니다: %s\n', save_path);

%% Helper Functions
% =========================================================================

function events_list = process_events(event_indices, V_measured, I_measured, t_measured, SOC_profile)
    events_list = {};
    for k = 1:size(event_indices, 1)
        s_idx = event_indices(k, 1);
        e_idx = event_indices(k, 2);
        
        if (s_idx + 1) <= e_idx
            delta_V = V_measured(s_idx + 1) - V_measured(s_idx);
            delta_I = I_measured(s_idx + 1) - I_measured(s_idx);
            
            if abs(delta_I) > 0.1
                R_transient_mOhm = (delta_V / delta_I) * 1000;
            else
                R_transient_mOhm = NaN;
            end
            
            event_data = struct(...
                'start_index', s_idx, 'end_index', e_idx, 'duration_s', (e_idx - s_idx) * 0.1, ...
                'R_transient_mOhm', R_transient_mOhm, 'SOC_start', SOC_profile(s_idx), ...
                'time_seq_s', t_measured(s_idx:e_idx), 'voltage_seq_V', V_measured(s_idx:e_idx), ...
                'current_seq_A', I_measured(s_idx:e_idx), 'soc_seq_pct', SOC_profile(s_idx:e_idx));
            events_list{end+1} = event_data;
        end
    end
end

function [dQdV_AhV, V_mid] = calculate_dQdV(Q_grid_Ah, V_ocv)
    dQ = diff(Q_grid_Ah); dV = diff(V_ocv); valid_indices = abs(dV) > 1e-6;
    dQdV_AhV = NaN(size(dV)); dQdV_AhV(valid_indices) = dQ(valid_indices) ./ dV(valid_indices);
    V_mid = V_ocv(1:end-1) + dV/2;
end

% [수정됨] 함수 정의부에서 '~'를 'Peak_Table'로 최종 수정
function [V_peak_main, dQdV_peak_main, Peak_Table] = find_main_peak(V_mid, dQdV, V_peak_ref, V_SEARCH_MIN, V_SEARCH_MAX, MIN_PROMINENCE)
    valid_mask = ~isnan(dQdV) & (V_mid >= V_SEARCH_MIN) & (V_mid <= V_SEARCH_MAX);
    V_search = V_mid(valid_mask); dQdV_search = dQdV(valid_mask); Peak_Table = table();
    if isempty(V_search), V_peak_main=NaN; dQdV_peak_main=NaN; return; end
    [pks, locs, ~, prom] = findpeaks(dQdV_search, V_search, 'MinPeakProminence', MIN_PROMINENCE);
    if isempty(pks), V_peak_main=NaN; dQdV_peak_main=NaN; return; end
    Peak_Table = table(locs', pks', prom', 'VariableNames', {'V_Peak', 'dQdV_Peak', 'Prominence'});
    if V_peak_ref > 0, [~, nearest_idx] = min(abs(Peak_Table.V_Peak - V_peak_ref));
    else, [~, nearest_idx] = max(Peak_Table.Prominence); end
    V_peak_main = Peak_Table.V_Peak(nearest_idx); dQdV_peak_main = Peak_Table.dQdV_Peak(nearest_idx);
end

function [LLI_shift_V, LAM_loss_rate] = quantify_lli_lam(V_peak_start, dQdV_peak_start, V_peak_end, dQdV_peak_end)
    LLI_shift_V = V_peak_start - V_peak_end;
    if dQdV_peak_start > 1e-6, LAM_loss_rate = (dQdV_peak_start - dQdV_peak_end) / dQdV_peak_start;
    else, LAM_loss_rate = NaN; end
end

% [수정됨] 함수 전체를 더 안정적이고 정확한 로직으로 변경
function [SOC, battery_capacity] = calculate_soc_profile(V, I, dt, OCV_data)
    % OCV 데이터로부터 SOC를 역으로 찾는 함수를 정의합니다.
    soc_grid = OCV_data.SOC_grid; 
    ocv_values = OCV_data.V_avg_SOC;
    [ocv_values_sorted, uniqueIdx] = unique(ocv_values, 'stable');
    soc_grid_sorted = soc_grid(uniqueIdx);
    inverse_OCV_func = @(voltage) interp1(ocv_values_sorted, soc_grid_sorted, voltage, 'linear', 'extrap');

    % [수정됨] 초기 휴지 구간을 찾는 로직 개선
    rest_mask = abs(I) <= 2.0;
    % 첫 번째 휴지 지점을 찾습니다.
    start_of_initial_rest = find(rest_mask, 1, 'first');
    
    % 만약 휴지 구간이 전혀 없다면, 첫 포인트를 기준으로 삼습니다.
    if isempty(start_of_initial_rest)
        initial_rest_end = 1;
        V_ocv_initial = V(1);
    else
        % 첫 휴지 구간이 끝나는 지점을 찾습니다.
        first_non_rest = find(~rest_mask(start_of_initial_rest:end), 1, 'first');
        if isempty(first_non_rest)
            % 파일 끝까지 휴지 상태인 경우
            initial_rest_end = length(I);
        else
            % 첫 비-휴지 지점 바로 앞이 휴지 구간의 끝입니다.
            initial_rest_end = start_of_initial_rest + first_non_rest - 2;
        end
        
        % [수정됨] OCV 계산을 위해 "정확한 휴지 구간"의 전압 평균을 사용합니다.
        V_ocv_initial = mean(V(start_of_initial_rest:initial_rest_end));
    end
    
    % 초기 SOC를 계산합니다.
    SOC_initial = inverse_OCV_func(V_ocv_initial);
    
    % 배터리 용량 및 SOC 벡터 초기화
    battery_capacity = OCV_data.mean_capacity;
    N = length(V); 
    SOC = zeros(N, 1);
    
    % [수정됨] 초기 휴지 구간 전체에 초기 SOC 값을 할당합니다.
    SOC(1:initial_rest_end) = SOC_initial;
    
    % 전류 적산을 통해 나머지 SOC를 계산합니다.
    for i = (initial_rest_end + 1):N
        SOC_change = (I(i) * dt) / (battery_capacity * 3600) * 100;
        SOC(i) = SOC(i-1) + SOC_change;
    end
    
    % SOC 값을 0~100 사이로 제한합니다.
    SOC = max(0, min(100, SOC));
end

function [chg_events, dch_events] = find_events_by_transition(I, idle_thr, min_duration_points)
    chg_events = []; dch_events = [];
    states = zeros(size(I));
    states(I > idle_thr) = 1; states(I < -idle_thr) = -1;
    
    charge_start_points = find(states(1:end-1) == 0 & states(2:end) == 1) + 1;
    for s = charge_start_points'
        end_candidate = find(states(s:end) ~= 1, 1, 'first');
        if isempty(end_candidate), e_B = length(states);
        else, e_B = s + end_candidate - 2; end
        s_A = s - 1;
        if (e_B > s_A) && ((e_B - s_A + 1) >= min_duration_points) && (isempty(chg_events) || s_A > chg_events(end, 2))
             chg_events = [chg_events; s_A, e_B];
        end
    end

    discharge_start_points = find(states(1:end-1) == 0 & states(2:end) == -1) + 1;
    for s = discharge_start_points'
        end_candidate = find(states(s:end) ~= -1, 1, 'first');
        if isempty(end_candidate), e_B = length(states);
        else, e_B = s + end_candidate - 2; end
        s_A = s - 1;
        if (e_B > s_A) && ((e_B - s_A + 1) >= min_duration_points) && (isempty(dch_events) || s_A > dch_events(end, 2))
             dch_events = [dch_events; s_A, e_B];
        end
    end
end