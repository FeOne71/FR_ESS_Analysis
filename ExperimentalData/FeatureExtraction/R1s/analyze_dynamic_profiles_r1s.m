%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 동적 주행 부하 기반 R1s 분석 스크립트 (v1 - R1s 기반 분석)
%
% 변경사항:
% 1. 주파수 분석 대신 R1s 계산으로 변경
% 2. R1s 계산 공식: R1s = dV / dI (온도 보정 없음)
% 3. SOC 범위 ±3%로 필터링
% 4. LLI/LAM과의 상관 분석 구조 유지
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
TARGET_SOC_LEVELS = {'SOC50'};
SOC_TOLERANCE = 3; % SOC 범위 ±3%
V_PEAK_SEARCH_MIN = 3.3; V_PEAK_SEARCH_MAX = 3.8; MIN_PEAK_PROMINENCE = 0.001;

%% 2. 통합 OCV 데이터 기반 LLI/LAM 정량화
% =========================================================================
fprintf('=== 1. 통합 OCV 데이터로부터 LLI/LAM 정량화 중... ===\n');
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

%% 3. SOC/DC 프로파일별 R1s 특징 추출
% =========================================================================
fprintf('\n=== 2. 주행 부하 프로파일로부터 R1s 특징 추출 중... ===\n');
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
        target_soc = str2double(soc_name(4:end)); % SOC90 -> 90
        
        for dc_idx = 1:8
            dc_name = sprintf('DC%d', dc_idx);
            if ~isfield(dc_data_all.(soc_name), dc_name) || isempty(dc_data_all.(soc_name).(dc_name).V), continue; end
            dc_data = dc_data_all.(soc_name).(dc_name);
            
            % R1s 계산을 위한 데이터 준비
            V = dc_data.V;
            I = dc_data.I;
            
            % SOC 계산 (OCV 데이터 기반) - 0사이클 OCV 데이터 사용
            ocv_ref_data = OCV_data.OCV_integrated_0;
            [SOC, ~] = calculate_soc_profile(V, I, 0.1, ocv_ref_data);
            
            % SOC 범위 ±3% 필터링
            soc_mask = (SOC >= (target_soc - SOC_TOLERANCE)) & (SOC <= (target_soc + SOC_TOLERANCE));
            if sum(soc_mask) < 10, continue; end % 최소 10개 포인트 필요
            
            V_filtered = V(soc_mask);
            I_filtered = I(soc_mask);
            
            % R1s 계산: dV/dI (1초 간격 = 10개 점 간격)
            dt_1s = 10; % 0.1초 * 10 = 1초
            dV = V_filtered(1+dt_1s:end) - V_filtered(1:end-dt_1s);
            dI = I_filtered(1+dt_1s:end) - I_filtered(1:end-dt_1s);
            
            % 유효한 데이터만 선택 (dI ≠ 0, NaN이 아닌 경우)
            valid_idx = (dI ~= 0) & ~isnan(dV) & ~isnan(dI);
            if sum(valid_idx) < 5, continue; end % 최소 5개 포인트 필요
            
            dV_valid = dV(valid_idx);
            dI_valid = dI(valid_idx);
            
            % R1s 계산 (mΩ 단위)
            R1s_values = (dV_valid ./ dI_valid) * 1000;
            
            % 양수 R1s만 선택
            positive_idx = R1s_values > 0;
            if sum(positive_idx) < 3, continue; end % 최소 3개 양수 값 필요
            
            R1s_positive = R1s_values(positive_idx);
            
            % R1s 통계량 계산
            R1s_mean = mean(R1s_positive);
            R1s_std = std(R1s_positive);
            R1s_median = median(R1s_positive);
            R1s_count = length(R1s_positive);
            
            new_row = {cycle_num, soc_name, dc_name, R1s_mean, R1s_std, R1s_median, R1s_count};
            feature_list = [feature_list; new_row];
        end
    end
end

feature_table = cell2table(feature_list, 'VariableNames', {'Cycle', 'SOC_Level', 'DC_Profile', 'R1s_Mean', 'R1s_Std', 'R1s_Median', 'R1s_Count'});

%% 4. 상세 디버깅이 포함된 상관 분석
% =========================================================================
fprintf('\n=== 3. R1s와 LLI/LAM 상관 분석 수행 중... ===\n');
final_summary_table = table();
r1s_features = {'R1s_Mean', 'R1s_Std', 'R1s_Median'};

for soc_name_cell = TARGET_SOC_LEVELS
    soc_name = soc_name_cell{1};
    fprintf('\n--------------------------------------------------\n');
    fprintf('           분석 시작: %s 그룹\n', soc_name);
    fprintf('--------------------------------------------------\n');
    
    for dc_idx = 1:8
        dc_name = sprintf('DC%d', dc_idx);
        group_data = feature_table(strcmp(feature_table.SOC_Level, soc_name) & strcmp(feature_table.DC_Profile, dc_name), :);
        if height(group_data) < 2, continue; end
        
        % [디버깅 간소화] 분석 대상 그룹 정보만 출력
        fprintf('\n--- 분석 중: [%s, %s] 그룹 (%d개 사이클) ---\n', soc_name, dc_name, height(group_data));
        analysis_table_raw = innerjoin(group_data, degradation_modes, 'Keys', 'Cycle');
        
        for r1s_type_cell = r1s_features
            r1s_type = r1s_type_cell{1};
            analysis_table = analysis_table_raw(:, {'Cycle', r1s_type, 'LLI', 'LAM'});
            
            % NaN 값 제거
            valid_rows = ~any(isnan(table2array(analysis_table)), 2);
            analysis_table = analysis_table(valid_rows, :);
            
            if height(analysis_table) < 2, continue; end
            
            mdl_lli = fitlm(analysis_table, [r1s_type ' ~ LLI']);
            mdl_lam = fitlm(analysis_table, [r1s_type ' ~ LAM']);
            
            summary_row = {soc_name, dc_name, r1s_type, height(analysis_table), mdl_lli.Rsquared.Ordinary, mdl_lam.Rsquared.Ordinary};
            final_summary_table = [final_summary_table; summary_row];
        end
    end
end

%% 5. 최종 결과 요약
% =========================================================================
fprintf('\n\n======================================================================\n');
fprintf('                최종 R1s 분석 결과 요약 (SOC/DC 분리 분석)\n');
fprintf('======================================================================\n');
if isempty(final_summary_table)
    fprintf('분석에 유효한 그룹이 없습니다.\n');
else
    final_summary_table.Properties.VariableNames = {'SOC_Level', 'DC_Profile', 'R1s_Feature', 'Data_Points', 'R2_vs_LLI', 'R2_vs_LAM'};
    final_summary_table = sortrows(final_summary_table, 'R2_vs_LAM', 'descend');
    disp(final_summary_table);
    
    % R1s 특징별 최고 상관관계 출력
    fprintf('\n=== R1s 특징별 최고 상관관계 (LAM 기준) ===\n');
    for r1s_type_cell = r1s_features
        r1s_type = r1s_type_cell{1};
        type_data = final_summary_table(strcmp(final_summary_table.R1s_Feature, r1s_type), :);
        if ~isempty(type_data)
            [~, best_idx] = max(type_data.R2_vs_LAM);
            best_row = type_data(best_idx, :);
            fprintf('%s: %s-%s (R²=%.4f)\n', r1s_type, best_row.SOC_Level{1}, best_row.DC_Profile{1}, best_row.R2_vs_LAM);
        end
    end
end

%% 6. SOC별 개별 R1s 값 Scatter Plot 시각화 (Subplot 사용)
% =========================================================================
fprintf('\n=== 4. SOC별 개별 R1s 값 Scatter Plot 시각화 중... ===\n');

% R1s 개별 값들을 저장할 구조체
individual_r1s_data = struct();

% 각 SOC별로 figure 생성
for soc_name_cell = TARGET_SOC_LEVELS
    soc_name = soc_name_cell{1};
    target_soc = str2double(soc_name(4:end)); % SOC90 -> 90
    
    % 해당 SOC의 모든 DC/Cycle 조합 데이터 찾기
    soc_data = feature_table(strcmp(feature_table.SOC_Level, soc_name), :);
    
    if height(soc_data) == 0
        fprintf('❌ %s: feature_table에 데이터가 없습니다.\n', soc_name);
        continue; 
    end
    
    fprintf('✅ %s: %d개 DC/Cycle 조합 발견\n', soc_name, height(soc_data));
    
    % 새로운 figure 생성 (8x4 subplot으로 확장)
    figure('Position', [100, 100, 2400, 2000]);
    sgtitle(sprintf('%s (%.0f±%d%%): 개별 R1s 값 분포 (DC별/Cycle별)', soc_name, target_soc, SOC_TOLERANCE), 'FontSize', 16, 'FontWeight', 'bold');
    
    subplot_count = 0;
    
    % DC1~DC8에 대해 subplot 생성
    for dc_idx = 1:8
        dc_name = sprintf('DC%d', dc_idx);
        
        % 해당 SOC/DC 조합의 모든 사이클 데이터 찾기
        soc_dc_data = soc_data(strcmp(soc_data.DC_Profile, dc_name), :);
        
        if height(soc_dc_data) == 0
            fprintf('  ❌ %s: %s 데이터 없음\n', soc_name, dc_name);
            continue; 
        end
        
        fprintf('  ✅ %s: %d개 사이클 발견\n', dc_name, height(soc_dc_data));
        
        % 각 사이클별로 개별 R1s 값들 추출
        for row_idx = 1:height(soc_dc_data)
            cycle_num = soc_dc_data.Cycle(row_idx);
            subplot_count = subplot_count + 1;
            
            if subplot_count > 32, break; end % 최대 32개 subplot (8x4)
            
            subplot(8, 4, subplot_count);
            
            % 원본 주행 사이클 데이터에서 해당 SOC/DC/Cycle 조합 찾기
            R1s_positive = [];
            for file_idx = 1:length(drive_cycle_files)
                current_filename = drive_cycle_files(file_idx).name;
                current_filepath = fullfile(drive_cycle_folder, current_filename);
                match = regexp(current_filename, 'parsedDriveCycle_(\w+)_filtered.mat', 'tokens');
                cycle_str = match{1}{1}; 
                file_cycle_num = str2double(regexp(cycle_str, '\d+', 'match'));
                
                if file_cycle_num ~= cycle_num, continue; end
                
                loaded_data = load(current_filepath);
                expected_var_name = ['parsedDriveCycle_', cycle_str];
                if ~isfield(loaded_data, expected_var_name), continue; end
                drive_cycle_data_struct = loaded_data.(expected_var_name);
                target_channel_name = [TARGET_CHANNEL_PREFIX, cycle_str];
                if ~isfield(drive_cycle_data_struct, target_channel_name), continue; end
                dc_data_all = drive_cycle_data_struct.(target_channel_name);
                
                if ~isfield(dc_data_all.(soc_name), dc_name) || isempty(dc_data_all.(soc_name).(dc_name).V), continue; end
                dc_data = dc_data_all.(soc_name).(dc_name);
                
                % R1s 계산을 위한 데이터 준비
                V = dc_data.V;
                I = dc_data.I;
                
                % SOC 계산 (OCV 데이터 기반) - 0사이클 OCV 데이터 사용
                ocv_ref_data = OCV_data.OCV_integrated_0;
                [SOC, ~] = calculate_soc_profile(V, I, 0.1, ocv_ref_data);
                
                % SOC 범위 ±3% 필터링
                soc_mask = (SOC >= (target_soc - SOC_TOLERANCE)) & (SOC <= (target_soc + SOC_TOLERANCE));
                if sum(soc_mask) < 10, continue; end % 최소 10개 포인트 필요
                
                V_filtered = V(soc_mask);
                I_filtered = I(soc_mask);
                
                % R1s 계산: dV/dI (1초 간격 = 10개 점 간격)
                dt_1s = 10; % 0.1초 * 10 = 1초
                dV = V_filtered(1+dt_1s:end) - V_filtered(1:end-dt_1s);
                dI = I_filtered(1+dt_1s:end) - I_filtered(1:end-dt_1s);
                
                % 유효한 데이터만 선택 (dI ≠ 0, NaN이 아닌 경우)
                valid_idx = (dI ~= 0) & ~isnan(dV) & ~isnan(dI);
                if sum(valid_idx) < 5, continue; end % 최소 5개 포인트 필요
                
                dV_valid = dV(valid_idx);
                dI_valid = dI(valid_idx);
                
                % R1s 계산 (mΩ 단위)
                R1s_values = (dV_valid ./ dI_valid) * 1000;
                
                % [디버깅] R1s 계산 검증
                if ~isempty(R1s_values)
                    fprintf('  %s-%s-Cycle%d: dV 범위 %.4f~%.4f V, dI 범위 %.4f~%.4f A, R1s 범위 %.1f~%.1f mΩ\n', ...
                            soc_name, dc_name, cycle_num, min(dV_valid), max(dV_valid), min(dI_valid), max(dI_valid), min(R1s_values), max(R1s_values));
                end
                
                % 양수 R1s만 선택
                positive_idx = R1s_values > 0;
                if sum(positive_idx) < 3, continue; end % 최소 3개 양수 값 필요
                
                R1s_positive = R1s_values(positive_idx);
                
                % 개별 R1s 값들을 저장
                data_key = sprintf('%s_%s_Cycle%d', soc_name, dc_name, cycle_num);
                individual_r1s_data.(data_key) = R1s_positive;
                
                break; % 해당 사이클을 찾았으므로 루프 종료
            end
            
            % Subplot에 dV vs dI scatter plot 표시
            if ~isempty(R1s_positive)
                % dV vs dI scatter plot (R1s = dV/dI)
                scatter(dI_valid(positive_idx), dV_valid(positive_idx), 30, 'filled', 'MarkerFaceColor', [0.2, 0.6, 0.8]);
                
                % 통계 정보 계산
                mean_r1s = mean(R1s_positive);
                std_r1s = std(R1s_positive);
                median_r1s = median(R1s_positive);
                
                % 직선 피팅 (y = ax 형태, 여기서 a = R1s)
                x_data = dI_valid(positive_idx);  % dI
                y_data = dV_valid(positive_idx);  % dV
                
                % 선형 회귀 (y = ax + b, 여기서 b ≈ 0이어야 함)
                p = polyfit(x_data, y_data, 1);
                slope = p(1);  % 기울기 (R1s)
                intercept = p(2);  % y절편 (b)
                
                % 피팅된 직선 계산
                x_fit = linspace(min(x_data), max(x_data), 100);
                y_fit = polyval(p, x_fit);
                
                % 직선 표시
                hold on;
                plot(x_fit, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('R1s=%.4f mΩ', slope*1000));
                hold off;
                
                % 제목 및 라벨
                title(sprintf('%s-Cycle%d\nR1s: %.4f mΩ', dc_name, cycle_num, slope*1000), 'FontSize', 10, 'FontWeight', 'bold');
                xlabel('dI (A)', 'FontSize', 8);
                ylabel('dV (V)', 'FontSize', 8);
                grid on;
                
                % X, Y축 범위 설정 (0,0을 중심으로)
                x_range = max(abs(x_data)) * 1.1;
                y_range = max(abs(y_data)) * 1.1;
                
                xlim([-x_range, x_range]);
                ylim([-y_range, y_range]);
                
                % 원점 (0,0) 표시
                hold on;
                plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'black');
                hold off;
                
                % 텍스트 정보 추가
                text(0.05, 0.95, sprintf('N=%d\nR1s=%.4f', length(R1s_positive), slope*1000), ...
                     'Units', 'normalized', 'VerticalAlignment', 'top', ...
                     'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 8);
            else
                % 데이터가 없는 경우
                title(sprintf('%s-Cycle%d\nNo Data', dc_name, cycle_num), 'FontSize', 10);
                xlabel('Data Points', 'FontSize', 8);
                ylabel('R1s (mΩ)', 'FontSize', 8);
                text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'red');
            end
        end
        
        if subplot_count >= 32, break; end % 최대 32개 subplot
    end
    
    % figure 저장
    save_dir = 'Results/R1s_ScatterPlots';
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    save_filename = sprintf('R1s_Scatter_%s.fig', soc_name);
    save_path = fullfile(save_dir, save_filename);
    saveas(gcf, save_path);
    fprintf('저장됨: %s\n', save_filename);
    
    % 메모리 절약을 위해 figure 닫기
    close(gcf);
end

fprintf('SOC별 개별 R1s 값 Scatter Plot 시각화 완료!\n');

%% 6-2. SOC별 전체 dV/dI Scatter Plot 시각화 (Tolerance 없음)
% =========================================================================
fprintf('\n=== 6-2. SOC별 전체 dV/dI Scatter Plot 시각화 (Tolerance 없음) ===\n');

% 각 SOC별로 figure 생성
for soc_name_cell = TARGET_SOC_LEVELS
    soc_name = soc_name_cell{1};
    target_soc = str2double(soc_name(4:end)); % SOC90 -> 90
    
    % 새로운 figure 생성 (8x4 subplot: 행=DC, 열=Cycle)
    figure('Position', [100, 100, 2400, 2000]);
    sgtitle(sprintf('%s (전체 SOC): dV/dI Scatter Plot (행=DC, 열=Cycle)', soc_name), 'FontSize', 16, 'FontWeight', 'bold');
    
    % 고정된 사이클 순서
    cycle_order = [0, 200, 400, 600];
    
    % DC1~DC8 (행), Cycle 0,200,400,600 (열)로 subplot 생성
    for dc_idx = 1:8
        dc_name = sprintf('DC%d', dc_idx);
        
        for cycle_col = 1:4
            cycle_num = cycle_order(cycle_col);
            
            % subplot 위치 계산 (행=DC, 열=Cycle)
            subplot(8, 4, (dc_idx-1)*4 + cycle_col);
            
            % 원본 주행 사이클 데이터에서 해당 SOC/DC/Cycle 조합 찾기
            R1s_positive = [];
            dV_valid = [];
            dI_valid = [];
            
            for file_idx = 1:length(drive_cycle_files)
                current_filename = drive_cycle_files(file_idx).name;
                current_filepath = fullfile(drive_cycle_folder, current_filename);
                match = regexp(current_filename, 'parsedDriveCycle_(\w+)_filtered.mat', 'tokens');
                cycle_str = match{1}{1}; 
                file_cycle_num = str2double(regexp(cycle_str, '\d+', 'match'));
                
                if file_cycle_num ~= cycle_num, continue; end
                
                loaded_data = load(current_filepath);
                expected_var_name = ['parsedDriveCycle_', cycle_str];
                if ~isfield(loaded_data, expected_var_name), continue; end
                drive_cycle_data_struct = loaded_data.(expected_var_name);
                target_channel_name = [TARGET_CHANNEL_PREFIX, cycle_str];
                if ~isfield(drive_cycle_data_struct, target_channel_name), continue; end
                dc_data_all = drive_cycle_data_struct.(target_channel_name);
                
                if ~isfield(dc_data_all.(soc_name), dc_name) || isempty(dc_data_all.(soc_name).(dc_name).V), continue; end
                dc_data = dc_data_all.(soc_name).(dc_name);
                
                % R1s 계산을 위한 데이터 준비
                V = dc_data.V;
                I = dc_data.I;
                
                % SOC 계산 (OCV 데이터 기반) - 0사이클 OCV 데이터 사용
                ocv_ref_data = OCV_data.OCV_integrated_0;
                [SOC, ~] = calculate_soc_profile(V, I, 0.1, ocv_ref_data);
                
                % SOC tolerance 없이 모든 데이터 사용
                V_filtered = V;
                I_filtered = I;
                
                % R1s 계산: dV/dI (1초 간격 = 10개 점 간격)
                dt_1s = 10; % 0.1초 * 10 = 1초
                dV = V_filtered(1+dt_1s:end) - V_filtered(1:end-dt_1s);
                dI = I_filtered(1+dt_1s:end) - I_filtered(1:end-dt_1s);
                
                % 유효한 데이터만 선택 (dI ≠ 0, NaN이 아닌 경우)
                valid_idx = (dI ~= 0) & ~isnan(dV) & ~isnan(dI);
                if sum(valid_idx) < 5, continue; end % 최소 5개 포인트 필요
                
                dV_valid = dV(valid_idx);
                dI_valid = dI(valid_idx);
                
                % R1s 계산 (mΩ 단위)
                R1s_values = (dV_valid ./ dI_valid) * 1000;
                
                % 양수 R1s만 선택
                positive_idx = R1s_values > 0;
                if sum(positive_idx) < 3, continue; end % 최소 3개 양수 값 필요
                
                R1s_positive = R1s_values(positive_idx);
                
                break; % 해당 사이클을 찾았으므로 루프 종료
            end
            
            % Subplot에 dV vs dI scatter plot 표시
            if ~isempty(R1s_positive) && ~isempty(dV_valid) && ~isempty(dI_valid)
                % dV vs dI scatter plot (R1s = dV/dI)
                scatter(dI_valid(positive_idx), dV_valid(positive_idx), 30, 'filled', 'MarkerFaceColor', [0.2, 0.6, 0.8]);
                
                % 통계 정보 계산
                mean_r1s = mean(R1s_positive);
                std_r1s = std(R1s_positive);
                median_r1s = median(R1s_positive);
                
                % 직선 피팅 (y = ax 형태, 여기서 a = R1s)
                x_data = dI_valid(positive_idx);  % dI
                y_data = dV_valid(positive_idx);  % dV
                
                % 선형 회귀 (y = ax + b, 여기서 b ≈ 0이어야 함)
                p = polyfit(x_data, y_data, 1);
                slope = p(1);  % 기울기 (R1s)
                intercept = p(2);  % y절편 (b)
                
                % 피팅된 직선 계산
                x_fit = linspace(min(x_data), max(x_data), 100);
                y_fit = polyval(p, x_fit);
                
                % 직선 표시
                hold on;
                plot(x_fit, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('R1s=%.4f mΩ', slope*1000));
                hold off;
                
                % 제목 및 라벨
                title(sprintf('%s-Cycle%d\nR1s: %.4f mΩ', dc_name, cycle_num, slope*1000), 'FontSize', 10, 'FontWeight', 'bold');
                xlabel('dI (A)', 'FontSize', 8);
                ylabel('dV (V)', 'FontSize', 8);
                grid on;
                
                % X, Y축 범위 설정 (0,0을 중심으로)
                x_range = max(abs(x_data)) * 1.1;
                y_range = max(abs(y_data)) * 1.1;
                
                xlim([-x_range, x_range]);
                ylim([-y_range, y_range]);
                
                % 원점 (0,0) 표시
                hold on;
                plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'black');
                hold off;
                
                % 텍스트 정보 추가
                text(0.05, 0.95, sprintf('N=%d\nR1s=%.4f', length(R1s_positive), slope*1000), ...
                     'Units', 'normalized', 'VerticalAlignment', 'top', ...
                     'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 8);
            else
                % 데이터가 없는 경우
                title(sprintf('%s-Cycle%d\nNo Data', dc_name, cycle_num), 'FontSize', 10);
                xlabel('dI (A)', 'FontSize', 8);
                ylabel('dV (V)', 'FontSize', 8);
                text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'red');
            end
        end
    end
    
    % figure 저장
    save_dir = 'Results/R1s_ScatterPlots';
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    save_filename = sprintf('R1s_Scatter_%s_AllSOC.fig', soc_name);
    save_path = fullfile(save_dir, save_filename);
    saveas(gcf, save_path);
    fprintf('저장됨: %s\n', save_filename);
    
    % 메모리 절약을 위해 figure 닫기
    close(gcf);
end

fprintf('SOC별 전체 dV/dI Scatter Plot 시각화 완료!\n');

%% 6-3. 전체 SOC 데이터 분석 디버깅
% =========================================================================
fprintf('\n=== 6-3. 전체 SOC 데이터 분석 디버깅 ===\n');

% 전체 SOC 데이터를 저장할 구조체
all_soc_r1s_data = struct();

% 각 SOC별로 전체 데이터 분석
for soc_name_cell = TARGET_SOC_LEVELS
    soc_name = soc_name_cell{1};
    target_soc = str2double(soc_name(4:end)); % SOC90 -> 90
    
    fprintf('\n--- %s (전체 SOC) 분석 ---\n', soc_name);
    
    % 해당 SOC의 모든 DC/Cycle 조합 데이터 찾기
    soc_data = feature_table(strcmp(feature_table.SOC_Level, soc_name), :);
    
    if height(soc_data) == 0
        fprintf('❌ %s: feature_table에 데이터가 없습니다.\n', soc_name);
        continue; 
    end
    
    fprintf('✅ %s: %d개 DC/Cycle 조합 발견\n', soc_name, height(soc_data));
    
    % 고정된 사이클 순서
    cycle_order = [0, 200, 400, 600];
    
    % DC별, Cycle별 통계 저장
    dc_cycle_stats = [];
    
    % DC1~DC8, Cycle 0,200,400,600 분석
    for dc_idx = 1:8
        dc_name = sprintf('DC%d', dc_idx);
        
        for cycle_col = 1:4
            cycle_num = cycle_order(cycle_col);
            
            % 원본 주행 사이클 데이터에서 해당 SOC/DC/Cycle 조합 찾기
            R1s_positive = [];
            dV_valid = [];
            dI_valid = [];
            soc_range = [];
            
            for file_idx = 1:length(drive_cycle_files)
                current_filename = drive_cycle_files(file_idx).name;
                current_filepath = fullfile(drive_cycle_folder, current_filename);
                match = regexp(current_filename, 'parsedDriveCycle_(\w+)_filtered.mat', 'tokens');
                cycle_str = match{1}{1}; 
                file_cycle_num = str2double(regexp(cycle_str, '\d+', 'match'));
                
                if file_cycle_num ~= cycle_num, continue; end
                
                loaded_data = load(current_filepath);
                expected_var_name = ['parsedDriveCycle_', cycle_str];
                if ~isfield(loaded_data, expected_var_name), continue; end
                drive_cycle_data_struct = loaded_data.(expected_var_name);
                target_channel_name = [TARGET_CHANNEL_PREFIX, cycle_str];
                if ~isfield(drive_cycle_data_struct, target_channel_name), continue; end
                dc_data_all = drive_cycle_data_struct.(target_channel_name);
                
                if ~isfield(dc_data_all.(soc_name), dc_name) || isempty(dc_data_all.(soc_name).(dc_name).V), continue; end
                dc_data = dc_data_all.(soc_name).(dc_name);
                
                % R1s 계산을 위한 데이터 준비
                V = dc_data.V;
                I = dc_data.I;
                
                % SOC 계산 (OCV 데이터 기반) - 0사이클 OCV 데이터 사용
                ocv_ref_data = OCV_data.OCV_integrated_0;
                [SOC, ~] = calculate_soc_profile(V, I, 0.1, ocv_ref_data);
                
                % SOC tolerance 없이 모든 데이터 사용
                V_filtered = V;
                I_filtered = I;
                
                % R1s 계산: dV/dI (1초 간격 = 10개 점 간격)
                dt_1s = 10; % 0.1초 * 10 = 1초
                dV = V_filtered(1+dt_1s:end) - V_filtered(1:end-dt_1s);
                dI = I_filtered(1+dt_1s:end) - I_filtered(1:end-dt_1s);
                
                % 유효한 데이터만 선택 (dI ≠ 0, NaN이 아닌 경우)
                valid_idx = (dI ~= 0) & ~isnan(dV) & ~isnan(dI);
                if sum(valid_idx) < 5, continue; end % 최소 5개 포인트 필요
                
                dV_valid = dV(valid_idx);
                dI_valid = dI(valid_idx);
                
                % R1s 계산 (mΩ 단위)
                R1s_values = (dV_valid ./ dI_valid) * 1000;
                
                % 양수 R1s만 선택
                positive_idx = R1s_values > 0;
                if sum(positive_idx) < 3, continue; end % 최소 3개 양수 값 필요
                
                R1s_positive = R1s_values(positive_idx);
                soc_range = [min(SOC), max(SOC)];
                
                break; % 해당 사이클을 찾았으므로 루프 종료
            end
            
            % 통계 계산 및 저장
            if ~isempty(R1s_positive)
                % 직선 피팅
                x_data = dI_valid(positive_idx);
                y_data = dV_valid(positive_idx);
                p = polyfit(x_data, y_data, 1);
                slope = p(1);
                
                % 통계 정보
                r1s_mean = mean(R1s_positive);
                r1s_std = std(R1s_positive);
                r1s_median = median(R1s_positive);
                data_count = length(R1s_positive);
                
                % SOC 범위 정보
                soc_min = soc_range(1);
                soc_max = soc_range(2);
                soc_span = soc_max - soc_min;
                
                % 통계 저장
                dc_cycle_stats = [dc_cycle_stats; {dc_name, cycle_num, data_count, r1s_mean, r1s_std, r1s_median, slope*1000, soc_min, soc_max, soc_span}];
                
                fprintf('  %s-Cycle%d: N=%d, R1s=%.4f±%.4f mΩ, SOC=%.1f~%.1f%% (span=%.1f%%)\n', ...
                        dc_name, cycle_num, data_count, r1s_mean, r1s_std, soc_min, soc_max, soc_span);
            else
                fprintf('  %s-Cycle%d: 데이터 없음\n', dc_name, cycle_num);
            end
        end
    end
    
    % 전체 SOC 데이터 저장
    all_soc_r1s_data.(soc_name) = dc_cycle_stats;
    
    % SOC별 요약 통계
    if ~isempty(dc_cycle_stats)
        all_r1s_values = cell2mat(dc_cycle_stats(:, 4)); % R1s_mean
        all_data_counts = cell2mat(dc_cycle_stats(:, 3)); % data_count
        all_soc_spans = cell2mat(dc_cycle_stats(:, 10)); % soc_span
        
        fprintf('\n%s 전체 요약:\n', soc_name);
        fprintf('  총 데이터 포인트: %d개\n', sum(all_data_counts));
        fprintf('  R1s 평균: %.4f ± %.4f mΩ\n', mean(all_r1s_values), std(all_r1s_values));
        fprintf('  R1s 범위: %.4f ~ %.4f mΩ\n', min(all_r1s_values), max(all_r1s_values));
        fprintf('  SOC 범위 평균: %.1f ± %.1f%%\n', mean(all_soc_spans), std(all_soc_spans));
        fprintf('  SOC 범위 최대: %.1f%%\n', max(all_soc_spans));
    end
end

% 전체 SOC vs Tolerance 적용 비교
fprintf('\n=== 전체 SOC vs Tolerance 적용 비교 ===\n');
for soc_name_cell = TARGET_SOC_LEVELS
    soc_name = soc_name_cell{1};
    
    % Tolerance 적용 데이터 (기존 feature_table)
    tolerance_data = feature_table(strcmp(feature_table.SOC_Level, soc_name), :);
    
    % 전체 SOC 데이터
    all_soc_data = all_soc_r1s_data.(soc_name);
    
    if ~isempty(tolerance_data) && ~isempty(all_soc_data)
        tolerance_r1s = tolerance_data.R1s_Mean * 1000; % mΩ 단위
        all_soc_r1s = cell2mat(all_soc_data(:, 4)); % R1s_mean
        
        fprintf('\n%s:\n', soc_name);
        fprintf('  Tolerance 적용: %d개 조합, R1s=%.4f±%.4f mΩ\n', height(tolerance_data), mean(tolerance_r1s), std(tolerance_r1s));
        fprintf('  전체 SOC: %d개 조합, R1s=%.4f±%.4f mΩ\n', height(all_soc_data), mean(all_soc_r1s), std(all_soc_r1s));
        fprintf('  차이: %.4f mΩ (%.1f%%)\n', mean(all_soc_r1s) - mean(tolerance_r1s), (mean(all_soc_r1s) - mean(tolerance_r1s))/mean(tolerance_r1s)*100);
    end
end

fprintf('\n전체 SOC 데이터 분석 디버깅 완료!\n');

%% 6-4. 전체 SOC R1s 결과 테이블 생성
% =========================================================================
fprintf('\n=== 6-4. 전체 SOC R1s 결과 테이블 생성 중... ===\n');

% 전체 SOC 사이클별 R1s 결과를 저장할 테이블 생성
all_soc_cycle_r1s_table = table();

% 각 SOC별로 테이블 생성
for soc_name_cell = TARGET_SOC_LEVELS
    soc_name = soc_name_cell{1};
    target_soc = str2double(soc_name(4:end)); % SOC90 -> 90
    
    fprintf('\n%s (전체 SOC) 테이블 생성 중...\n', soc_name);
    
    % 전체 SOC 데이터 가져오기
    all_soc_data = all_soc_r1s_data.(soc_name);
    
    if isempty(all_soc_data), continue; end
    
    % 사이클별로 그룹화
    cycles = unique(cell2mat(all_soc_data(:, 2))); % Cycle 컬럼
    
    for cycle_num = cycles'
        % 해당 사이클의 모든 DC 데이터
        cycle_mask = cell2mat(all_soc_data(:, 2)) == cycle_num;
        cycle_data = all_soc_data(cycle_mask, :);
        
        if isempty(cycle_data), continue; end
        
        % R1s 값들 추출 (mΩ 단위)
        r1s_values = cell2mat(cycle_data(:, 4)); % R1s_mean
        r1s_std_values = cell2mat(cycle_data(:, 5)); % R1s_std
        data_counts = cell2mat(cycle_data(:, 3)); % data_count
        soc_spans = cell2mat(cycle_data(:, 10)); % soc_span
        
        % 통계 계산
        r1s_mean = mean(r1s_values);
        r1s_std = std(r1s_values);
        total_n = sum(data_counts);
        avg_soc_span = mean(soc_spans);
        
        % 0사이클일 때의 R1s 값 찾기 (기준값)
        if cycle_num == 0
            r1s_reference = r1s_mean; % 0사이클의 평균 R1s를 기준으로 사용
        else
            % 0사이클 데이터가 없으면 현재 값 사용
            r1s_reference = r1s_mean;
        end
        
        % 테이블에 행 추가
        new_row = table(target_soc, r1s_reference, total_n, r1s_mean, r1s_std, avg_soc_span, ...
                       'VariableNames', {'SOC_target', 'R1s', 'N', 'R1s_mean', 'R1s_std', 'SOC_span'});
        all_soc_cycle_r1s_table = [all_soc_cycle_r1s_table; new_row];
        
        fprintf('  Cycle %d: R1s=%.4f mΩ, N=%d, Mean=%.4f±%.4f mΩ, SOC_span=%.1f%%\n', ...
                cycle_num, r1s_reference, total_n, r1s_mean, r1s_std, avg_soc_span);
    end
end

% 테이블 정렬 (SOC_target, R1s 기준)
all_soc_cycle_r1s_table = sortrows(all_soc_cycle_r1s_table, {'SOC_target', 'R1s'}, {'ascend', 'descend'});

% 결과 출력
fprintf('\n=== 전체 SOC 사이클별 R1s 결과 테이블 ===\n');
disp(all_soc_cycle_r1s_table);

%% 7. 사이클별 R1s 결과 테이블 생성
% =========================================================================
fprintf('\n=== 7. 사이클별 R1s 결과 테이블 생성 중... ===\n');

% 사이클별 R1s 결과를 저장할 테이블 생성
cycle_r1s_table = table();

% 각 SOC별로 테이블 생성
for soc_name_cell = TARGET_SOC_LEVELS
    soc_name = soc_name_cell{1};
    target_soc = str2double(soc_name(4:end)); % SOC90 -> 90
    
    fprintf('\n%s (%.0f±%d%%) 처리 중...\n', soc_name, target_soc, SOC_TOLERANCE);
    
    % 해당 SOC의 모든 사이클 데이터 찾기
    soc_data = feature_table(strcmp(feature_table.SOC_Level, soc_name), :);
    
    if height(soc_data) == 0, continue; end
    
    % 사이클별로 그룹화
    cycles = unique(soc_data.Cycle);
    
    for cycle_num = cycles'
        % 해당 사이클의 모든 DC 데이터
        cycle_data = soc_data(soc_data.Cycle == cycle_num, :);
        
        if height(cycle_data) == 0, continue; end
        
        % R1s 값들 추출 (mΩ 단위)
        r1s_values = cycle_data.R1s_Mean * 1000; % mΩ 단위
        r1s_std_values = cycle_data.R1s_Std * 1000; % mΩ 단위
        data_counts = cycle_data.R1s_Count;
        
        % 통계 계산
        r1s_mean = mean(r1s_values);
        r1s_std = std(r1s_values);
        total_n = sum(data_counts);
        
        % 0사이클일 때의 R1s 값 찾기 (기준값)
        if cycle_num == 0
            r1s_reference = r1s_mean; % 0사이클의 평균 R1s를 기준으로 사용
        else
            % 0사이클 데이터가 없으면 현재 값 사용
            r1s_reference = r1s_mean;
        end
        
        % 테이블에 행 추가
        new_row = table(target_soc, r1s_reference, total_n, r1s_mean, r1s_std, ...
                       'VariableNames', {'SOC_target', 'R1s', 'N', 'R1s_mean', 'R1s_std'});
        cycle_r1s_table = [cycle_r1s_table; new_row];
        
        fprintf('  Cycle %d: R1s=%.4f mΩ, N=%d, Mean=%.4f±%.4f mΩ\n', ...
                cycle_num, r1s_reference, total_n, r1s_mean, r1s_std);
    end
end

% 테이블 정렬 (SOC_target, R1s 기준)
cycle_r1s_table = sortrows(cycle_r1s_table, {'SOC_target', 'R1s'}, {'ascend', 'descend'});

% 결과 출력
fprintf('\n=== 사이클별 R1s 결과 테이블 ===\n');
disp(cycle_r1s_table);

%% 8. 결과 저장
% =========================================================================
save_path = 'Results/R1s_Analysis_Results.mat';
if ~exist('Results', 'dir'), mkdir('Results'); end
save(save_path, 'feature_table', 'degradation_modes', 'final_summary_table', 'cycle_r1s_table', 'all_soc_cycle_r1s_table', 'all_soc_r1s_data', '-v7.3');
fprintf('\n분석 결과가 저장되었습니다: %s\n', save_path);

% CSV 파일로도 저장
csv_path = 'Results/R1s_Cycle_Results.csv';
writetable(cycle_r1s_table, csv_path);
fprintf('사이클별 R1s 결과가 CSV로 저장되었습니다: %s\n', csv_path);

% 전체 SOC 결과도 CSV로 저장
all_soc_csv_path = 'Results/R1s_AllSOC_Cycle_Results.csv';
writetable(all_soc_cycle_r1s_table, all_soc_csv_path);
fprintf('전체 SOC 사이클별 R1s 결과가 CSV로 저장되었습니다: %s\n', all_soc_csv_path);

%% Helper Functions
% =========================================================================

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

function [SOC, battery_capacity] = calculate_soc_profile(V, I, dt, ocv_ref_data)
    % OCV 데이터로부터 SOC를 역으로 찾는 함수를 정의합니다.
    soc_grid = ocv_ref_data.SOC_grid; 
    ocv_values = ocv_ref_data.V_avg_SOC;
    [ocv_values_sorted, uniqueIdx] = unique(ocv_values, 'stable');
    soc_grid_sorted = soc_grid(uniqueIdx);
    inverse_OCV_func = @(voltage) interp1(ocv_values_sorted, soc_grid_sorted, voltage, 'linear', 'extrap');

    % 초기 휴지 구간을 찾는 로직 개선
    rest_mask = abs(I) <= 64*0.05;
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
        
        % OCV 계산을 위해 "정확한 휴지 구간"의 전압 평균을 사용합니다.
        V_ocv_initial = mean(V(start_of_initial_rest:initial_rest_end));
    end
    
    % 초기 SOC를 계산합니다.
    SOC_initial = inverse_OCV_func(V_ocv_initial);
    
    % 배터리 용량 및 SOC 벡터 초기화
    battery_capacity = ocv_ref_data.mean_capacity;
    N = length(V); 
    SOC = zeros(N, 1);
    
    % 초기 휴지 구간 전체에 초기 SOC 값을 할당합니다.
    SOC(1:initial_rest_end) = SOC_initial;
    
    % 전류 적산을 통해 나머지 SOC를 계산합니다.
    for i = (initial_rest_end + 1):N
        SOC_change = (I(i) * dt) / (battery_capacity * 3600) * 100;
        SOC(i) = SOC(i-1) + SOC_change;
    end
    
    % SOC 값을 0~100 사이로 제한합니다.
    SOC = max(0, min(100, SOC));
end

%% 5. 통계 분석 디버깅
% =========================================================================
fprintf('\n=== 5. R1s 통계 분석 디버깅 ===\n');

% 기본 통계 정보
fprintf('\n--- 기본 통계 정보 ---\n');
fprintf('총 데이터 개수: %d\n', height(feature_table));
fprintf('SOC 레벨: %s\n', strjoin(unique(feature_table.SOC_Level), ', '));
fprintf('DC 프로파일: %s\n', strjoin(unique(feature_table.DC_Profile), ', '));
fprintf('사이클: %s\n', strjoin(string(unique(feature_table.Cycle)), ', '));

% SOC별 통계 요약
fprintf('\n--- SOC별 통계 요약 ---\n');
soc_levels = unique(feature_table.SOC_Level);
for i = 1:length(soc_levels)
    soc_data = feature_table(strcmp(feature_table.SOC_Level, soc_levels{i}), :);
    r1s_values = soc_data.R1s_Mean * 1000; % mΩ 단위
    fprintf('%s: %d개, R1s = %.4f ± %.4f mΩ (%.4f ~ %.4f)\n', ...
        soc_levels{i}, height(soc_data), mean(r1s_values), std(r1s_values), ...
        min(r1s_values), max(r1s_values));
end

% DC별 통계 요약
fprintf('\n--- DC별 통계 요약 ---\n');
dc_profiles = unique(feature_table.DC_Profile);
for i = 1:length(dc_profiles)
    dc_data = feature_table(strcmp(feature_table.DC_Profile, dc_profiles{i}), :);
    r1s_values = dc_data.R1s_Mean * 1000; % mΩ 단위
    fprintf('%s: %d개, R1s = %.4f ± %.4f mΩ (%.4f ~ %.4f)\n', ...
        dc_profiles{i}, height(dc_data), mean(r1s_values), std(r1s_values), ...
        min(r1s_values), max(r1s_values));
end

% 이상치 탐지
fprintf('\n--- 이상치 탐지 ---\n');
all_r1s = feature_table.R1s_Mean * 1000; % mΩ 단위
Q1 = quantile(all_r1s, 0.25);
Q3 = quantile(all_r1s, 0.75);
IQR = Q3 - Q1;
lower_bound = Q1 - 1.5 * IQR;
upper_bound = Q3 + 1.5 * IQR;

outliers = all_r1s < lower_bound | all_r1s > upper_bound;
fprintf('전체 R1s 데이터: %d개\n', length(all_r1s));
fprintf('이상치 개수: %d개 (%.1f%%)\n', sum(outliers), sum(outliers)/length(all_r1s)*100);
fprintf('이상치 범위: < %.4f mΩ 또는 > %.4f mΩ\n', lower_bound, upper_bound);

% 데이터 품질 체크
fprintf('\n--- 데이터 품질 체크 ---\n');
fprintf('R1s 값이 0인 데이터: %d개\n', sum(feature_table.R1s_Mean == 0));
fprintf('R1s 값이 음수인 데이터: %d개\n', sum(feature_table.R1s_Mean < 0));
fprintf('R1s 값이 NaN인 데이터: %d개\n', sum(isnan(feature_table.R1s_Mean)));

% 요약
fprintf('\n--- 요약 ---\n');
fprintf(' 총 %d개의 SOC-DC-Cycle 조합 분석 완료\n', height(feature_table));
fprintf(' R1s 평균: %.4f ± %.4f mΩ\n', mean(all_r1s), std(all_r1s));
fprintf(' R1s 범위: %.4f ~ %.4f mΩ\n', min(all_r1s), max(all_r1s));
fprintf(' 이상치 비율: %.1f%%\n', sum(outliers)/length(all_r1s)*100);
