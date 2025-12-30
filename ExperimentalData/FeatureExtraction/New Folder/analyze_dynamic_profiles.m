%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 동적 주행 부하 기반 '유사 EIS' 분석 스크립트 (v6 - 상세 디버깅 최종 복원)
%
% 변경사항:
% 1. [디버깅 강화] 각 분석 그룹(예: SOC90-DC1)별로 상관 분석에 사용된
%    실제 데이터 테이블(사이클별 임피던스, LLI, LAM)을 모두 출력
% 2. 분석 구조를 사용자의 의도에 맞게 (SOC레벨 -> DC프로파일) 명확히 유지
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. 설정 및 데이터 로딩
% =========================================================================
ocv_data_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\OCV_integrated.mat';
load(ocv_data_path);
if ~exist('OCV_data', 'var'), error("'OCV_data' 변수를 찾을 수 없습니다."); end

drive_cycle_folder = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';

% --- 분석 파라미터 ---
TARGET_CHANNEL_PREFIX = 'ch9_Drive_';
TARGET_SOC_LEVELS = {'SOC90', 'SOC70', 'SOC50'};
TARGET_FREQUENCIES = [1.0, 0.5, 0.1, 0.05];
V_PEAK_SEARCH_MIN = 3.3; V_PEAK_SEARCH_MAX = 3.8; MIN_PEAK_PROMINENCE = 0.001;

%% 2. 통합 OCV 데이터 기반 LLI/LAM 정량화
% =========================================================================
fprintf('=== 1. 통합 OCV 데이터로부터 LLI/LAM 정량화 중... ===\n');
% (이전과 동일)
all_fields = fieldnames(OCV_data);
q_grid_fields = all_fields(startsWith(all_fields, 'q_grid_rpt'));
cycle_keys_str = cellfun(@(s) s(11:end), q_grid_fields, 'UniformOutput', false);
[~, sort_idx] = sort(cellfun(@str2double, cycle_keys_str));
cycle_keys = cycle_keys_str(sort_idx);
Q_data = struct(); V_data = struct();
for i = 1:length(cycle_keys), key = cycle_keys{i}; Q_data.(['c' key]) = OCV_data.(['q_grid_rpt' key]); V_data.(['c' key]) = OCV_data.(['avg_ocv_rpt' key]); end
dQdV_results = struct(); V_mid_points = struct(); V_peaks = struct(); dQdV_peaks = struct();
for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    [dQdV_results.(['c' key]), V_mid_points.(['c' key])] = calculate_dQdV(Q_data.(['c' key]), V_data.(['c' key]));
    if i == 1, ref_peak_V = 0; else, ref_peak_V = V_peaks.(['c' cycle_keys{i-1}]); end
    [V_peaks.(['c' key]), dQdV_peaks.(['c' key]), ~] = find_main_peak(V_mid_points.(['c' key]), dQdV_results.(['c' key]), ref_peak_V, V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);
end
lli_values = zeros(length(cycle_keys), 1); lam_values = zeros(length(cycle_keys), 1);
base_peak_V = V_peaks.c0; base_peak_dQdV = dQdV_peaks.c0;
for i = 2:length(cycle_keys)
    key = cycle_keys{i};
    lli_values(i) = base_peak_V - V_peaks.(['c' key]);
    lam_values(i) = (base_peak_dQdV - dQdV_peaks.(['c' key])) / base_peak_dQdV;
end
all_cycles = cellfun(@str2double, cycle_keys);
degradation_modes = table(all_cycles, lli_values, lam_values, 'VariableNames', {'Cycle', 'LLI', 'LAM'});

%% 3. SOC/DC 프로파일별 '유사 임피던스' 특징 추출
% =========================================================================
fprintf('\n=== 2. 주행 부하 프로파일로부터 유사 임피던스 특징 추출 중... ===\n');
% (이전과 동일)
feature_list = [];
drive_cycle_files = dir(fullfile(drive_cycle_folder, 'parsedDriveCycle_*_filtered.mat'));
for file_idx = 1:length(drive_cycle_files)
    current_filename = drive_cycle_files(file_idx).name;
    current_filepath = fullfile(drive_cycle_folder, current_filename);
    match = regexp(current_filename, 'parsedDriveCycle_(\w+)_filtered.mat', 'tokens');
    cycle_str = match{1}{1}; cycle_num = str2double(regexp(cycle_str, '\d+', 'match'));
    loaded_data = load(current_filepath);
    expected_var_name = ['parsedDriveCycle_', cycle_str];
    if ~isfield(loaded_data, expected_var_name), continue; end
    drive_cycle_data_struct = loaded_data.(expected_var_name);
    target_channel_name = [TARGET_CHANNEL_PREFIX, cycle_str];
    if ~isfield(drive_cycle_data_struct, target_channel_name), continue; end
    dc_data_all = drive_cycle_data_struct.(target_channel_name);
    for soc_name_cell = TARGET_SOC_LEVELS
        soc_name = soc_name_cell{1};
        for dc_idx = 1:8
            dc_name = sprintf('DC%d', dc_idx);
            if ~isfield(dc_data_all.(soc_name), dc_name) || isempty(dc_data_all.(soc_name).(dc_name).V), continue; end
            dc_data = dc_data_all.(soc_name).(dc_name);
            I_ac = dc_data.I - mean(dc_data.I);
            V_ac = dc_data.V - mean(dc_data.V);
            Fs = 10;
            [Z_est, f_axis] = tfestimate(I_ac, V_ac, hanning(length(I_ac)), [], [], Fs);
            impedance_features = zeros(1, length(TARGET_FREQUENCIES));
            for f_idx = 1:length(TARGET_FREQUENCIES)
                [~, idx] = min(abs(f_axis - TARGET_FREQUENCIES(f_idx)));
                impedance_features(f_idx) = abs(Z_est(idx));
            end
            new_row = [{cycle_num, soc_name, dc_name}, num2cell(impedance_features)];
            feature_list = [feature_list; new_row];
        end
    end
end
feature_table = cell2table(feature_list, 'VariableNames', {'Cycle', 'SOC_Level', 'DC_Profile', 'Z_1Hz', 'Z_0_5Hz', 'Z_0_1Hz', 'Z_0_05Hz'});

%% 4. [핵심 변경] 상세 디버깅이 포함된 상관 분석
% =========================================================================
fprintf('\n=== 3. 상관 분석 수행 중... ===\n');
final_summary_table = table();
resistance_types = {'Z_1Hz', 'Z_0_5Hz', 'Z_0_1Hz', 'Z_0_05Hz'};

% SOC/DC별로 분리하여 분석 (원래 구조 유지)
for soc_idx = 1:length(TARGET_SOC_LEVELS)
    soc_name = TARGET_SOC_LEVELS{soc_idx};
    fprintf('\n--------------------------------------------------\n');
    fprintf('           분석 시작: %s 그룹 (%d/%d)\n', soc_name, soc_idx, length(TARGET_SOC_LEVELS));
    fprintf('--------------------------------------------------\n');
    
    for dc_idx = 1:8
        dc_name = sprintf('DC%d', dc_idx);
        group_data = feature_table(strcmp(feature_table.SOC_Level, soc_name) & strcmp(feature_table.DC_Profile, dc_name), :);
        if height(group_data) < 2, continue; end
        
        % [디버깅 추가] 분석 대상 그룹 정보 및 실제 데이터 테이블 출력
        fprintf('\n--- 분석 중: [%s, %s] 그룹 ---\n', soc_name, dc_name);
        analysis_table_raw = innerjoin(group_data, degradation_modes, 'Keys', 'Cycle');
        disp('  >> 분석용 데이터 테이블:');
        disp(analysis_table_raw);
        
        for r_type_cell = resistance_types
            r_type = r_type_cell{1};
            analysis_table = analysis_table_raw(:, {'Cycle', r_type, 'LLI', 'LAM'});
            
            mdl_lli = fitlm(analysis_table, [r_type ' ~ LLI']);
            mdl_lam = fitlm(analysis_table, [r_type ' ~ LAM']);
            
            % 각 사이클별로 결과 저장
            for cycle_idx = 1:height(analysis_table)
                cycle_num = analysis_table.Cycle(cycle_idx);
                summary_row = {cycle_num, soc_name, dc_name, r_type, mdl_lli.Rsquared.Ordinary, mdl_lam.Rsquared.Ordinary, mdl_lli.Coefficients.pValue(2), mdl_lam.Coefficients.pValue(2)};
                final_summary_table = [final_summary_table; summary_row];
            end
        end
    end
end

%% 5. 최종 결과 요약
% =========================================================================
fprintf('\n\n======================================================================\n');
fprintf('                최종 분석 결과 요약 (사이클별 분석)\n');
fprintf('======================================================================\n');
if isempty(final_summary_table)
    fprintf('분석에 유효한 그룹이 없습니다.\n');
else
    final_summary_table.Properties.VariableNames = {'Cycle', 'SOC_Level', 'DC_Profile', 'Impedance_Feature', 'R2_vs_LLI', 'R2_vs_LAM', 'PValue_vs_LLI', 'PValue_vs_LAM'};
    
    % 사이클별로 정렬 (오름차순), 그 다음 SOC, DC, Impedance 순으로 정렬
    soc_numeric = cellfun(@(x) str2double(x(4:end)), final_summary_table.SOC_Level);
    dc_numeric = cellfun(@(x) str2double(x(3:end)), final_summary_table.DC_Profile);
    imp_order = {'Z_1Hz', 'Z_0_5Hz', 'Z_0_1Hz', 'Z_0_05Hz'};
    [~, imp_numeric] = ismember(final_summary_table.Impedance_Feature, imp_order);
    
    % 복합 정렬: Cycle (오름차순), SOC (내림차순), DC (오름차순), Impedance (오름차순)
    [~, sort_idx] = sortrows([final_summary_table.Cycle, -soc_numeric, dc_numeric, imp_numeric]);
    final_summary_table = final_summary_table(sort_idx, :);
    
    disp(final_summary_table);
end

%% Helper Functions
% (이전과 동일)
% ...
function [dQdV_AhV, V_mid] = calculate_dQdV(Q_grid_Ah, V_ocv)
    dQ = diff(Q_grid_Ah); dV = diff(V_ocv); valid_indices = abs(dV) > 1e-6;
    dQdV_AhV = NaN(size(dV)); dQdV_AhV(valid_indices) = dQ(valid_indices) ./ dV(valid_indices);
    V_mid = V_ocv(1:end-1) + dV/2;
end

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