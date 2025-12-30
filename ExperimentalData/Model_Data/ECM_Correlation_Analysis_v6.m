%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 파일명: ECM_Correlation_Analysis_v12.m
% 기능:
%   - Part 1: RPT 데이터로 열화 모드(LLI, LAM) 정량화
%   - Part 2: '채널 10'의 SOC50 주행부하 데이터에 2RC ECM을 피팅하여 파라미터 추출
%   - Part 3: LLI/LAM과 ECM 파라미터(R0, R1, R2) 간의 직접적인 상관관계 분석
% 수정사항:
%   - [핵심] 신뢰도 높은 전역 최적해를 찾기 위해 Part 2의 fmincon 최적화를
%     'MultiStart' 알고리즘으로 전면 교체
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
warning off;

%% ========================================================================
%  1. 기본 설정 (사용자 수정 영역)
% =========================================================================

% --- 저장 경로 ---
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Model_Data\ECM_Correlation_Ch10';

% --- 데이터 경로 ---
rptDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated';
driveCycleDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';

% --- ECM 피팅 파라미터 ---
params_initial = [0.001, 0.003, 0.01, 30, 600];
lb = [0.00005, 0.0001, 0.0005, 0.1, 1];
ub = [0.05, 0.05, 0.1, 500, 5000];

% fmincon 로컬 솔버 옵션
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxIterations', 500, 'MaxFunctionEvaluations', 2000, ...
    'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-8);

% MultiStart 글로벌 솔버 설정
num_start_points = 20; % MultiStart 시도 횟수
% =========================================================================

if ~exist(saveDir, 'dir'); mkdir(saveDir); end
fprintf('모든 결과는 다음 폴더에 저장됩니다:\n%s\n\n', saveDir);


%% ========================================================================
%  Part 1: dQ/dV 기반 열화 모드 정량화
% =========================================================================
fprintf('################### Part 1: 열화 모드(LLI, LAM) 정량화 시작 ###################\n');

% --- 1a. RPT 데이터 로딩 및 양방향 OCV-SOC 함수 생성 ---
rptMatFile = fullfile(rptDataPath, 'OCV_integrated.mat');
load(rptMatFile, 'OCV_data');
all_fields = fieldnames(OCV_data);
q_grid_fields = all_fields(startsWith(all_fields, 'q_grid_rpt'));
cycle_keys_str = cellfun(@(s) s(11:end), q_grid_fields, 'UniformOutput', false);
[~, sort_idx] = sort(cellfun(@str2double, cycle_keys_str));
cycle_keys = cycle_keys_str(sort_idx);
fprintf('분석 대상 RPT 사이클: %s\n', strjoin(cycle_keys, ', '));

Q_data = struct(); V_data = struct(); 
soc_functions = struct(); ocv_functions = struct();

for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    Q_data.(['c' key]) = OCV_data.(['q_grid_rpt' key]);
    V_data.(['c' key]) = OCV_data.(['avg_ocv_rpt' key]);
    soc_grid = OCV_data.soc_grid; avg_ocv = OCV_data.(['avg_ocv_rpt' key]);
    
    [unique_ocv, idx_v] = unique(avg_ocv);
    unique_soc_for_v = soc_grid(idx_v);
    soc_functions.(['c' key]) = @(v_query) interp1(unique_ocv, unique_soc_for_v, v_query, 'linear', 'extrap');
    
    [unique_soc, idx_s] = unique(soc_grid);
    unique_ocv_for_s = avg_ocv(idx_s);
    ocv_functions.(['c' key]) = @(soc_query) interp1(unique_soc, unique_ocv_for_s, soc_query, 'linear', 'extrap');
end

% --- 1b. LLI, LAM 계산 및 테이블 생성 ---
V_PEAK_SEARCH_MIN = 3.3; V_PEAK_SEARCH_MAX = 3.8; MIN_PEAK_PROMINENCE = 0.001;
degradation_modes_table = table();

base_key = cycle_keys{1};
[dQdV_base, V_mid_base] = calculate_dQdV(Q_data.(['c' base_key]), V_data.(['c' base_key]));
[V_peak_base, dQdV_peak_base, ~] = find_main_peak(V_mid_base, dQdV_base, 0, V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);

new_row = table(str2double(base_key), 0, 0, 'VariableNames', {'Cycle', 'LLI_mV', 'LAM_pct'});
degradation_modes_table = [degradation_modes_table; new_row];

for i = 2:length(cycle_keys)
    curr_key = cycle_keys{i};
    [dQdV_curr, V_mid_curr] = calculate_dQdV(Q_data.(['c' curr_key]), V_data.(['c' curr_key]));
    [V_peak_curr, dQdV_peak_curr, ~] = find_main_peak(V_mid_curr, dQdV_curr, V_peak_base, V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);
    
    if ~isnan(V_peak_curr)
        [lli_V, lam_rate] = quantify_lli_lam(V_peak_base, dQdV_peak_base, V_peak_curr, dQdV_peak_curr);
        new_row = table(str2double(curr_key), lli_V * 1000, lam_rate * 100, 'VariableNames', {'Cycle', 'LLI_mV', 'LAM_pct'});
        degradation_modes_table = [degradation_modes_table; new_row];
    else
        fprintf('  > 경고: %s 사이클에서 dQ/dV 피크를 찾지 못해 LLI/LAM 분석에서 제외합니다.\n', curr_key);
    end
end
fprintf('LLI 및 LAM 정량화 완료.\n');


%% ========================================================================
%  Part 2: SOC50 주행부하 ECM 파라미터 추출 (채널 10 한정)
% =========================================================================
fprintf('\n\n################### Part 2: ECM 파라미터 추출 시작 (Channel 10) ###################\n');

ECM_results = table();
for i = 1:length(cycle_keys)
    cycle_key = cycle_keys{i};
    fprintf('\n--- 처리 중인 사이클: %s ---\n', cycle_key);
    
    matFileName = fullfile(driveCycleDataDir, sprintf('parsedDriveCycle_%scyc_filtered.mat', cycle_key));
    if ~exist(matFileName, 'file'), fprintf('  > 경고: %s 파일을 찾을 수 없어 건너뜁니다.\n', matFileName); continue; end
    
    current_soc_func = soc_functions.(['c' cycle_key]); 
    current_ocv_func = ocv_functions.(['c' cycle_key]);
    data = load(matFileName); data_struct_name = fieldnames(data); drive_data = data.(data_struct_name{1});
    
    temp_params = [];
    channel_names = fieldnames(drive_data);
    for ch_idx = 1:length(channel_names)
        if ~contains(channel_names{ch_idx}, 'ch10'), continue; end

        if ~isfield(drive_data.(channel_names{ch_idx}), 'SOC50'), continue; end
        soc50_data = drive_data.(channel_names{ch_idx}).SOC50;
        profile_names = fieldnames(soc50_data);
        
        for p_idx = 1:length(profile_names)
            profile_name = profile_names{p_idx};
            profile_data = soc50_data.(profile_name);
            fprintf('  > 피팅 중: %s - SOC50 - %s\n', channel_names{ch_idx}, profile_name);
            
            V_meas = profile_data.V(:); I_meas = profile_data.I(:); t_meas = profile_data.t(:);
            if isduration(t_meas), t_meas = seconds(t_meas); end
            
            SOC_est = calculate_soc_profile(I_meas, t_meas, V_meas(1), current_soc_func);
            if isempty(SOC_est)
                fprintf('      > 경고: 초기 SOC 추정 실패. 이 프로파일을 건너뜁니다.\n');
                continue;
            end
            
            % ★★★★★ [수정] MultiStart를 이용한 최적화 수행 ★★★★★
            cost_func = @(params) cost_function_2RC(params, I_meas, V_meas, t_meas, SOC_est, current_ocv_func);
            problem = createOptimProblem('fmincon', 'objective', cost_func, 'x0', params_initial, 'lb', lb, 'ub', ub, 'options', options);
            
            ms = MultiStart('Display', 'off');
            [optimal_params, rmse_val] = run(ms, problem, num_start_points);
            
            fprintf('      > 피팅 완료. RMSE: %.4f V\n', rmse_val);
            temp_params = [temp_params; optimal_params];
        end
    end
    
    if ~isempty(temp_params)
        avg_params = mean(temp_params, 1);
        fprintf('  > %s 사이클 평균 ECM 파라미터: R0=%.2fmOhm, R1=%.2fmOhm, R2=%.2fmOhm, t1=%.1fs, t2=%.1fs\n', ...
            cycle_key, avg_params(1)*1000, avg_params(2)*1000, avg_params(3)*1000, avg_params(4), avg_params(5));
        
        new_row = table(str2double(cycle_key), avg_params(1), avg_params(2), avg_params(3), avg_params(4), avg_params(5), ...
            'VariableNames', {'Cycle', 'R0', 'R1', 'R2', 'tau1', 'tau2'});
        ECM_results = [ECM_results; new_row];
    end
end
fprintf('\n모든 사이클의 ECM 파라미터 추출 완료.\n');


%% ========================================================================
%  Part 3: 열화 모드와 ECM 파라미터 상관관계 분석
% =========================================================================
fprintf('\n\n################### Part 3: 상관관계 분석 시작 ###################\n');

% --- 3a. 최종 분석 테이블 생성 ---
final_table = outerjoin(degradation_modes_table, ECM_results, 'Keys', 'Cycle', 'MergeKeys', true);
final_table = rmmissing(final_table);

% mOhm 단위로 변환
final_table.R0_mOhm = final_table.R0 * 1000;
final_table.R1_mOhm = final_table.R1 * 1000;
final_table.R2_mOhm = final_table.R2 * 1000;

fprintf('\n--- 최종 분석 데이터 테이블 ---\n');
disp(final_table);

% --- 3b. 상관관계 시각화 ---
if height(final_table) < 2
    fprintf('상관분석을 위한 데이터가 부족합니다 (최소 2개 사이클 필요).\n');
    return;
end

fig_corr = figure('Name', 'Correlation Analysis (Channel 10)', 'Position', [100 100, 1000, 800]);
t = tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t, 'Correlation between Degradation Modes and ECM Parameters (Channel 10, SOC 50)', 'FontSize', 16, 'FontWeight', 'bold');

x_vars = {'LLI_mV', 'LAM_pct'};
x_labels = {'LLI (Voltage Shift) [mV]', 'LAM (Peak Height Loss) [%]'};
y_vars = {'R0_mOhm', 'R1_mOhm', 'R2_mOhm'};
y_labels = {'R0 [m\Omega]', 'R1 [m\Omega]', 'R2 [m\Omega]'};

nexttile;
plot_correlation(final_table.LLI_mV, final_table.R0_mOhm, x_labels{1}, y_labels{1}, 'R0');
hold on;
plot_correlation(final_table.LLI_mV, final_table.R1_mOhm, x_labels{1}, y_labels{2}, 'R1');
legend('R0 data', 'R0 fit', 'R1 data', 'R1 fit', 'Location', 'best');
title('Impact of LLI on Resistances');

nexttile;
plot_correlation(final_table.LAM_pct, final_table.R0_mOhm, x_labels{2}, y_labels{1}, 'R0');
hold on;
plot_correlation(final_table.LAM_pct, final_table.R2_mOhm, x_labels{2}, y_labels{3}, 'R2');
legend('R0 data', 'R0 fit', 'R2 data', 'R2 fit', 'Location', 'best');
title('Impact of LAM on Resistances');

nexttile([1 2]);
plot(final_table.Cycle, final_table.LLI_mV, 'o-', 'LineWidth', 2, 'DisplayName', 'LLI (mV)');
hold on;
plot(final_table.Cycle, final_table.LAM_pct, 's-', 'LineWidth', 2, 'DisplayName', 'LAM (%)');
ylabel('LLI / LAM');
yyaxis right;
plot(final_table.Cycle, final_table.R0_mOhm, '^-', 'LineWidth', 2, 'DisplayName', 'R0 (mOhm)');
plot(final_table.Cycle, final_table.R1_mOhm, 'd-', 'LineWidth', 2, 'DisplayName', 'R1 (mOhm)');
plot(final_table.Cycle, final_table.R2_mOhm, '*-', 'LineWidth', 2, 'DisplayName', 'R2 (mOhm)');
ylabel('ECM Resistances (mOhm)');
grid on;
title('Trend of All Parameters over Cycles');
xlabel('Cycle');
legend('show', 'Location', 'best');

savefig(fig_corr, fullfile(saveDir, 'ECM_Parameter_Correlation_Ch10.fig'));
fprintf('\nPart 3: 상관관계 분석 그래프 저장 완료.\n');
close(fig_corr);

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

function SOC_est = calculate_soc_profile(I, t, V_initial, soc_func)
    N = length(t); SOC_est = zeros(N, 1);
    initial_soc = soc_func(V_initial);
    if isnan(initial_soc), SOC_est = []; return; end
    SOC_est(1) = initial_soc;
    Q_nominal_As = 64 * 3600;
    for k = 1:N-1
        dt = t(k+1) - t(k);
        SOC_est(k+1) = SOC_est(k) - (I(k) * dt) / Q_nominal_As * 100;
        SOC_est(k+1) = max(0, min(100, SOC_est(k+1)));
    end
end

function rmse = cost_function_2RC(params, I, V_meas, t, SOC_est, ocv_func)
    R0=params(1); R1=params(2); R2=params(3); tau1=params(4); tau2=params(5);
    dt = mean(diff(t)); N = length(t);
    V_model = zeros(N, 1); V_rc1 = 0; V_rc2 = 0;
    V_ocv = ocv_func(SOC_est);
    exp_tau1 = exp(-dt / tau1); exp_tau2 = exp(-dt / tau2);
    for k = 1:N-1
        V_rc1 = V_rc1 * exp_tau1 + R1 * (1 - exp_tau1) * I(k);
        V_rc2 = V_rc2 * exp_tau2 + R2 * (1 - exp_tau2) * I(k);
        V_model(k+1) = V_ocv(k+1) - I(k+1) * R0 - V_rc1 - V_rc2;
    end
    V_model(1) = V_meas(1);
    rmse = sqrt(mean((V_meas(2:end) - V_model(2:end)).^2));
    if any(isnan(V_model)) || rmse > 1, rmse = 1e6; end
end

function plot_correlation(x_data, y_data, x_label, y_label, series_name)
    colors = containers.Map({'R0', 'R1', 'R2'}, {'b', 'g', 'm'});
    markers = containers.Map({'R0', 'R1', 'R2'}, {'o', 's', 'd'});
    color = colors(series_name);
    marker = markers(series_name);
    scatter(x_data, y_data, 100, color, marker, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    p = polyfit(x_data, y_data, 1);
    x_fit = linspace(min(x_data), max(x_data), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, '-', 'Color', color, 'LineWidth', 2);
    R = corrcoef(x_data, y_data);
    R_val = R(1,2);
    text_y_pos = 0.9 - (find(strcmp(series_name, {'R0','R1','R2'}))-1)*0.15;
    text(0.1, text_y_pos, sprintf('%s, R=%.3f', series_name, R_val), 'Units', 'normalized', ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', color);
    xlabel(x_label);
    ylabel(y_label);
    grid on;
end