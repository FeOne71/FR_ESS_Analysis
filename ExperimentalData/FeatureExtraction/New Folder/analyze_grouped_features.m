%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 그룹화 기반 열화 모드 상관 분석 스크립트 (v13 - 이벤트 DB 저장)
%
% 변경사항:
% 1. [DB 저장 기능 추가] 분석에 사용된 '안정적인 이벤트'의 전체 데이터를
%    (SOC레벨/사이클) 계층 구조에 맞춰 'Filtered_Event_Database.mat'로 저장
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. 설정 및 데이터 로딩
% =========================================================================
analysis_data_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureExtraction\Results\FR_ESS_Analysis_Results_All_Cycles.mat';
load(analysis_data_path);
ocv_data_folder = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated';
ocv_data_path = fullfile(ocv_data_folder, 'OCV_integrated.mat');
load(ocv_data_path);
saveDir_Features = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureExtraction\Results';
if ~exist(saveDir_Features, 'dir'); mkdir(saveDir_Features); end

% --- 필터링 및 그룹화 파라미터 ---
MIN_DURATION_FILTER = 10; 
STABILITY_THRESHOLD = 0.05;
CURRENT_BINS = 5:5:60;
FINE_SOC_BINS = 0:5:100;

%% 2. [핵심 변경] 이벤트 필터링, 특성 추출 및 DB 저장
% =========================================================================
fprintf('=== 1. SOC 레벨별로 안정적인 이벤트 필터링 및 DB 저장 중... ===\n');
soc_level_tables = struct();
% [추가됨] 최종 이벤트 데이터베이스를 위한 구조체 초기화
Filtered_Events_DB = struct();

SOC_LEVEL_NAMES = {'SOC90', 'SOC70', 'SOC50'};
cycles = fieldnames(final_results.drive_cycle_analysis);

for soc_name_cell = SOC_LEVEL_NAMES
    soc_name = soc_name_cell{1};
    event_feature_list = [];
    Filtered_Events_DB.(soc_name) = struct(); % SOC 레벨별 구조체 생성
    
    for c_idx = 1:length(cycles)
        cycle_key = cycles{c_idx};
        cycle_num = str2double(regexp(cycle_key, '\d+', 'match'));
        
        cycle_db_name = ['Cycle_' num2str(cycle_num)];
        Filtered_Events_DB.(soc_name).(cycle_db_name) = {}; % 사이클별 셀 배열 초기화
        
        if isfield(final_results.drive_cycle_analysis.(cycle_key), soc_name)
            current_soc_level_data = final_results.drive_cycle_analysis.(cycle_key).(soc_name);
            dc_profiles = fieldnames(current_soc_level_data);
            
            for d_idx = 1:length(dc_profiles)
                if isfield(current_soc_level_data.(dc_profiles{d_idx}), 'charge_events')
                    events = current_soc_level_data.(dc_profiles{d_idx}).charge_events;
                    for e_idx = 1:length(events)
                        event = events{e_idx};
                        if event.duration_s >= MIN_DURATION_FILTER && length(event.current_seq_A) >= 101
                            current_window = event.current_seq_A(31:101);
                            if (std(current_window) / mean(current_window)) < STABILITY_THRESHOLD
                                
                                % ----------------- 필터링 통과한 이벤트 -----------------
                                
                                % 1. 특성 추출 (상관 분석용)
                                R = @(t_idx) (event.voltage_seq_V(t_idx) - event.voltage_seq_V(1)) / (event.current_seq_A(t_idx) - event.current_seq_A(1)) * 1000;
                                R_1s = R(11); R_5s = R(51); R_10s = R(101);
                                new_row = {cycle_num, event.SOC_start, mean(current_window), R_1s, R_5s, R_10s};
                                event_feature_list = [event_feature_list; new_row];
                                
                                % 2. [추가됨] DB 저장을 위한 데이터 구조체 생성
                                event_to_save = struct();
                                event_to_save.SOC = event.soc_seq_pct;
                                event_to_save.Voltage = event.voltage_seq_V;
                                event_to_save.Current = event.current_seq_A;
                                event_to_save.Power = event.voltage_seq_V .* event.current_seq_A; % Power 계산
                                event_to_save.Duration_s = event.duration_s;
                                event_to_save.Resistances_mOhm = struct('R_1s', R_1s, 'R_5s', R_5s, 'R_10s', R_10s);
                                
                                % 최종 DB에 추가
                                Filtered_Events_DB.(soc_name).(cycle_db_name){end+1} = event_to_save;
                            end
                        end
                    end
                end
            end
        end
    end
    
    if ~isempty(event_feature_list)
        stable_events_table = cell2table(event_feature_list, 'VariableNames', {'Cycle', 'SOC_start', 'Stable_Current_A', 'R_1s', 'R_5s', 'R_10s'});
        soc_level_tables.(soc_name) = stable_events_table;
    end
end

% [추가됨] 최종 이벤트 데이터베이스 파일로 저장
db_save_path = fullfile(saveDir_Features, 'Filtered_Event_Database.mat');
save(db_save_path, 'Filtered_Events_DB', '-v7.3');
fprintf('\n>>> 필터링된 전체 이벤트 데이터가 다음 파일에 저장되었습니다:\n    %s\n', db_save_path);

%% 3. 상세 디버깅 및 자동 분석
% =========================================================================
% (이후 분석 로직은 이전과 동일)
fprintf('\n=== 2. 모든 유효 그룹에 대한 자동 상관 분석 수행 ===\n');
% ... (이하 동일)
final_summary_table = table();
all_cycles = [0, 200, 400, 600];
lli_full = [0; cumsum(final_results.degradation_quantified.LLI_LAM.LLI_V)'];
lam_full = [0; final_results.degradation_quantified.LLI_LAM.LAM_rate'];
degradation_modes = table(all_cycles', lli_full, lam_full, 'VariableNames', {'Cycle', 'LLI_Cumulative_V', 'LAM_Rate'});
resistance_types = {'R_1s', 'R_5s', 'R_10s'};

for soc_name_cell = SOC_LEVEL_NAMES
    soc_name = soc_name_cell{1};
    if ~isfield(soc_level_tables, soc_name), continue; end
    
    fprintf('\n--------------------------------------------------\n');
    fprintf('           분석 시작: %s 그룹\n', soc_name);
    fprintf('--------------------------------------------------\n');
    
    event_feature_table = soc_level_tables.(soc_name);
    
    event_feature_table.Current_Group_Idx = discretize(event_feature_table.Stable_Current_A, CURRENT_BINS);
    event_feature_table.Fine_SOC_Group_Idx = discretize(event_feature_table.SOC_start, FINE_SOC_BINS);
    
    summary_counts = groupsummary(event_feature_table, {'Cycle', 'Current_Group_Idx', 'Fine_SOC_Group_Idx'});
    valid_counts = summary_counts(~isnan(summary_counts.Current_Group_Idx) & ~isnan(summary_counts.Fine_SOC_Group_Idx), :);
    
    if ~isempty(valid_counts)
        valid_counts.Current_Group = CURRENT_BINS(valid_counts.Current_Group_Idx)';
        valid_counts.Fine_SOC_Group = FINE_SOC_BINS(valid_counts.Fine_SOC_Group_Idx)';
        
        fprintf('--- [%s] 그룹별 이벤트 수 요약 ---\n', soc_name);
        unique_current_groups_in_data = unique(valid_counts.Current_Group);
        for uc = unique_current_groups_in_data'
            fprintf('>> Current Group: %d-%d A\n', uc, uc+5);
            pivot_data = valid_counts(valid_counts.Current_Group == uc, :);
            pivot_table = unstack(pivot_data, 'GroupCount', 'Cycle', 'GroupingVariables', 'Fine_SOC_Group');
            pivot_table.Properties.VariableNames = regexprep(pivot_table.Properties.VariableNames, 'x_(\d+)', 'Cycle_$1');
            disp(pivot_table);
        end
    else
        fprintf('[%s]에서 유효한 그룹이 없습니다.\n', soc_name);
        continue;
    end

    for r_type_cell = resistance_types
        r_type = r_type_cell{1};
        grouped_stats = groupsummary(event_feature_table, {'Cycle', 'Current_Group_Idx', 'Fine_SOC_Group_Idx'}, 'mean', r_type);
        valid_rows = ~isnan(grouped_stats.Current_Group_Idx) & ~isnan(grouped_stats.Fine_SOC_Group_Idx);
        grouped_stats = grouped_stats(valid_rows, :);
        if isempty(grouped_stats), continue; end
        
        grouped_stats.Current_Group = CURRENT_BINS(grouped_stats.Current_Group_Idx)';
        grouped_stats.Fine_SOC_Group = FINE_SOC_BINS(grouped_stats.Fine_SOC_Group_Idx)';
        
        unique_groups = unique(grouped_stats(:, {'Current_Group', 'Fine_SOC_Group'}), 'rows');
        
        for i = 1:height(unique_groups)
            current_current_group = unique_groups.Current_Group(i);
            current_fine_soc_group = unique_groups.Fine_SOC_Group(i);
            
            target_group_data = grouped_stats(grouped_stats.Current_Group == current_current_group & grouped_stats.Fine_SOC_Group == current_fine_soc_group, :);
            if height(target_group_data) < 2, continue; end
            
            analysis_table = innerjoin(target_group_data, degradation_modes, 'Keys', 'Cycle');
            analysis_table.Properties.VariableNames{['mean_' r_type]} = 'Mean_Resistance';
            
            mdl_lli = fitlm(analysis_table, 'Mean_Resistance ~ LLI_Cumulative_V');
            mdl_lam = fitlm(analysis_table, 'Mean_Resistance ~ LAM_Rate');
            mdl_multi = fitlm(analysis_table, 'Mean_Resistance ~ LLI_Cumulative_V + LAM_Rate');
            
            summary_row = {soc_name, current_current_group, current_fine_soc_group, r_type, height(target_group_data), ...
                           mdl_lli.Rsquared.Ordinary, mdl_lam.Rsquared.Ordinary, mdl_multi.Rsquared.Ordinary};
            final_summary_table = [final_summary_table; summary_row];
        end
    end
end

%% 4. 최종 분석 결과 요약
% =========================================================================
fprintf('\n\n==========================================================================================\n');
fprintf('                                 최종 분석 결과 요약\n');
fprintf('==========================================================================================\n');
if isempty(final_summary_table)
    fprintf('분석에 유효한 그룹(2개 이상 사이클 데이터 보유)이 없습니다.\n');
else
    final_summary_table.Properties.VariableNames = {'SOC_Level', 'Current_A', 'Fine_SOC', 'Resistance', 'Data_Points', 'R2_vs_LLI', 'R2_vs_LAM', 'R2_vs_LLI_LAM'};
    final_summary_table = sortrows(final_summary_table, 'R2_vs_LLI_LAM', 'descend');
    disp(final_summary_table);
end