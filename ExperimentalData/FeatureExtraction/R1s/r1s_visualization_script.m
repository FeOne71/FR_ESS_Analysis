%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R1s Visualization Script
%
% Purpose: Load R1s data from ECM fitting results and visualize
% Function: Generate dV vs dI scatter plots for each cycle
% SOC Filtering: SOC50 ± 1%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. Settings and Data Loading
% =========================================================================
fprintf('=== R1s Visualization Script Start ===\n');

% Data path settings
results_path = 'Results/ECM_Fitting/ECM_2RC_Fitting_Results.mat';

% Load result data
if ~exist(results_path, 'file')
    error('ECM fitting result file not found: %s', results_path);
end

load(results_path);
fprintf('ECM fitting results loaded: %s\n', results_path);

%% 2. R1s 데이터 추출 및 정리
% =========================================================================
fprintf('\n=== R1s 데이터 추출 ===\n');

% 결과 키 목록 생성
result_keys = fieldnames(ecm_results);
valid_results = {};
r1s_data = [];

for i = 1:length(result_keys)
    key = result_keys{i};
    result = ecm_results.(key);
    
    % R1s 계산에 필요한 데이터가 있는지 확인
    if isfield(result, 'V_measured') && isfield(result, 'I_measured') && isfield(result, 'SOC') && ~isempty(result.V_measured)
        valid_results{end+1} = key;
        
        % R1s 계산 (전체 구간, SOC 필터링 없음)
        fprintf('  Cycle %d: 전체 데이터 %d개 사용\n', result.Cycle, length(result.V_measured));
        
        if length(result.V_measured) >= 20
            V_filtered = result.V_measured;
            I_filtered = result.I_measured;
            
            % R1s 계산 (1초 간격)
            dt_1s = 10;
            dV = V_filtered(1+dt_1s:end) - V_filtered(1:end-dt_1s);
            dI = I_filtered(1+dt_1s:end) - I_filtered(1:end-dt_1s);
            
            % 유효한 데이터 필터링
            valid_idx = (dI ~= 0) & ~isnan(dV) & ~isnan(dI);
            dV_valid = dV(valid_idx);
            dI_valid = dI(valid_idx);
            
            if length(dV_valid) >= 5
                % 선형 회귀 피팅: dV = a * dI (y = ax)
                p = polyfit(dI_valid, dV_valid, 1);
                slope = p(1);
                
                % R1s (mΩ 단위)
                R1s_val = slope * 1000;
                
                % R1s 데이터 수집
                r1s_data = [r1s_data; result.SOC_target, result.Cycle, R1s_val, R1s_val, 0, R1s_val, length(dV_valid)];
            end
        end
    end
end

if isempty(valid_results)
    error('시각화할 R1s 데이터가 없습니다.');
end

fprintf('유효한 R1s 데이터: %d개\n', length(valid_results));

%% 3. 사이클별 그룹화
% =========================================================================
fprintf('\n=== 사이클별 그룹화 ===\n');

% 사이클 번호 추출
cycles = unique([ecm_results.(valid_results{1}).Cycle]);
for i = 2:length(valid_results)
    cycles = [cycles, ecm_results.(valid_results{i}).Cycle];
end
cycles = unique(cycles);

fprintf('사이클 번호: %s\n', mat2str(cycles));

% 사이클별로 그룹화
cycle_data = struct();
for i = 1:length(valid_results)
    key = valid_results{i};
    cycle_num = ecm_results.(key).Cycle;
    
    if ~isfield(cycle_data, sprintf('cycle_%d', cycle_num))
        cycle_data.(sprintf('cycle_%d', cycle_num)) = {};
    end
    cycle_data.(sprintf('cycle_%d', cycle_num)){end+1} = key;
end

%% 4. R1s 시각화 생성
% =========================================================================
fprintf('\n=== R1s 시각화 생성 ===\n');

% 시각화 생성
figure('Position', [100, 100, 1600, 1000]);

% subplot 설정
num_cycles = length(cycles);
if num_cycles <= 4
    subplot_rows = 2;
    subplot_cols = 2;
else
    subplot_rows = 2;
    subplot_cols = ceil(num_cycles / 2);
end

for i = 1:length(cycles)
    cycle_num = cycles(i);
    cycle_key = sprintf('cycle_%d', cycle_num);
    
    subplot(subplot_rows, subplot_cols, i);
    
    if isfield(cycle_data, cycle_key) && ~isempty(cycle_data.(cycle_key))
        % 해당 사이클의 데이터가 있는 경우
        key = cycle_data.(cycle_key){1}; % 첫 번째 결과 사용
        result = ecm_results.(key);
        
        % R1s 데이터 재계산 (시각화용, 전체 구간)
        if length(result.V_measured) >= 20
            V_filtered = result.V_measured;
            I_filtered = result.I_measured;
            
            % R1s 계산 (1초 간격)
            dt_1s = 10;
            dV = V_filtered(1+dt_1s:end) - V_filtered(1:end-dt_1s);
            dI = I_filtered(1+dt_1s:end) - I_filtered(1:end-dt_1s);
            
            % 유효한 데이터 필터링
            valid_idx = (dI ~= 0) & ~isnan(dV) & ~isnan(dI);
            dV_valid = dV(valid_idx);
            dI_valid = dI(valid_idx);
            
            if length(dV_valid) >= 5
                % 선형 회귀 피팅: dV = a * dI (y = ax)
                p = polyfit(dI_valid, dV_valid, 1);
                slope = p(1);
                % intercept는 사용하지 않음 (y = ax 형태)
                
                % scatter plot
                scatter(dI_valid, dV_valid, 20, 'b', 'filled');
                hold on;
                
                % 피팅선 그리기 (y = ax)
                x_fit = linspace(min(dI_valid), max(dI_valid), 100);
                y_fit = slope * x_fit;
                plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
                
                % 축 설정
                xlabel('dI (A)', 'FontWeight', 'bold');
                ylabel('dV (V)', 'FontWeight', 'bold');
                title(sprintf('Cycle %d\nR1s = %.4f mΩ', cycle_num, slope*1000), 'FontWeight', 'bold');
                grid on;
                
                % 범례
                legend(sprintf('N=%d', length(dV_valid)), sprintf('R1s=%.4f mΩ', slope*1000), 'Location', 'best');
                
                fprintf('  Cycle %d: R1s = %.4f mΩ (N=%d)\n', cycle_num, slope*1000, length(dV_valid));
            else
                % 데이터 부족
                text(0.5, 0.5, 'Insufficient Data', 'HorizontalAlignment', 'center', 'FontSize', 12);
                title(sprintf('Cycle %d\nNo Data', cycle_num), 'FontWeight', 'bold');
                fprintf('  Cycle %d: 데이터 부족 (%d개)\n', cycle_num, length(dV_valid));
            end
        else
            % SOC 필터링 데이터 부족
            text(0.5, 0.5, 'Insufficient SOC Data', 'HorizontalAlignment', 'center', 'FontSize', 12);
            title(sprintf('Cycle %d\nNo SOC Data', cycle_num), 'FontWeight', 'bold');
            fprintf('  Cycle %d: SOC 필터링 데이터 부족 (%d개)\n', cycle_num, sum(soc_mask));
        end
    else
        % 해당 사이클의 데이터가 없는 경우
        text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', 'FontSize', 12);
        title(sprintf('Cycle %d\nNo Data', cycle_num), 'FontWeight', 'bold');
        fprintf('  Cycle %d: 데이터 없음\n', cycle_num);
    end
end

% 전체 제목
sgtitle('R1s Analysis: dV vs dI Scatter Plots (전체 구간)', 'FontSize', 16, 'FontWeight', 'bold');

%% 5. 시각화 저장
% =========================================================================
save_dir = 'Results/ECM_Fitting/Visualizations';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

fig_filename = 'R1s_Analysis_Scatter_Plots.fig';
fig_path = fullfile(save_dir, fig_filename);
savefig(fig_path);

fprintf('\nR1s 시각화 저장 완료: %s\n', fig_filename);

%% 6. R1s 통계 요약 (표 형태)
% =========================================================================
fprintf('\n=== R1s 통계 요약 ===\n');

if ~isempty(r1s_data)
    % 사이클별 통계 요약 (표 형태)
    fprintf('\n=== 사이클별 R1s 통계 요약 ===\n');
    fprintf('총 처리된 조합: %d개\n', size(r1s_data, 1));
    
    % 사이클별 분석
    cycles = unique(r1s_data(:, 2));
    
    % 표 헤더
    fprintf('\nCycle | R1s (mΩ)    | N_points | SOC_range\n');
    fprintf('------|-------------|----------|----------\n');
    
    for i = 1:length(cycles)
        cycle_num = cycles(i);
        cycle_mask = r1s_data(:, 2) == cycle_num;
        cycle_data = r1s_data(cycle_mask, :);
        
        % 값 계산
        r1s_val = cycle_data(1, 3);  % R1s 값
        n_points = cycle_data(1, 7); % 데이터 포인트 수
        soc_range = sprintf('%.1f~%.1f', min(result.SOC), max(result.SOC)); % SOC 범위
        
        fprintf('%5d | %11.4f | %8d | %s\n', ...
            cycle_num, r1s_val, n_points, soc_range);
    end
    
    % 전체 통계
    r1s_values = r1s_data(:, 3); % R1s 값들
    fprintf('\n=== 전체 R1s 통계 ===\n');
    fprintf('평균 R1s: %.4f ± %.4f mΩ\n', mean(r1s_values), std(r1s_values));
    fprintf('최소 R1s: %.4f mΩ\n', min(r1s_values));
    fprintf('최대 R1s: %.4f mΩ\n', max(r1s_values));
end

fprintf('\n=== R1s 시각화 완료 ===\n');
