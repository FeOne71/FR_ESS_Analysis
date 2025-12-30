%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power Spectral Density (PSD) Analysis with Debugging
%
% Purpose: 초보자를 위해 PSD 분석의 각 단계를 설명하고, 중간 결과를
%          출력하여 이해를 돕는 스크립트.
%
% Usage:
%   1. ecm_fitting_drive_cycle_1029.m을 먼저 실행하여 결과 .mat 파일 생성.
%   2. 아래 설정(Configuration) 섹션을 수정.
%   3. 이 스크립트를 실행하고 Command Window에 출력되는 메시지를 확인.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. Configuration
% =========================================================================
fprintf('=== PSD 분석 스크립트 (디버깅 모드) 시작 ===\n');

% --- 사용자 설정 ---
channel_name_for_results = 'ch9';
% 분석할 특정 결과 하나만 선택 (이해를 돕기 위해 하나씩 보는 것을 추천)
% 예: TARGET_RESULTS = {'SOC50_DC1_Cycle0'};
TARGET_RESULTS = {'SOC50_DC1_Cycle0'};
% --- 설정 끝 ---

results_dir = fullfile('Results', 'ECM_Fitting', channel_name_for_results);
results_file = sprintf('ECM_2RC_Fitting_Results_%s.mat', channel_name_for_results);
results_path = fullfile(results_dir, results_file);
save_dir = fullfile('Results', 'PSD_Analysis_Debug', channel_name_for_results);
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% 2. Load Data
% =========================================================================
fprintf('\n[단계 1: 데이터 로드]\n');
fprintf('피팅 결과 파일 로딩: %s\n', results_path);
if ~exist(results_path, 'file'), error('결과 파일이 없습니다!'); end
load(results_path, 'ecm_results');
fprintf('>> 성공: ecm_results 구조체를 로드했습니다.\n');

%% 3. Data Extraction and Pre-analysis
% =========================================================================
if isempty(TARGET_RESULTS)
    result_key = fieldnames(ecm_results);
    result_key = result_key{1}; % 비어있으면 첫 번째 결과 사용
else
    result_key = TARGET_RESULTS{1};
end

fprintf('\n[단계 2: 분석 대상 데이터 준비]\n');
fprintf('분석 대상: %s\n', result_key);

result_data = ecm_results.(result_key);

% 데이터 추출
I = result_data.I_measured;
t = result_data.t_measured;
params_optimal = [result_data.R0, result_data.R1, result_data.R2, result_data.tau1, result_data.tau2];

% --- 디버깅 코드 ---
fprintf('\n--- 디버깅: 데이터 기본 정보 확인 ---\n');
fprintf('  - 전류 데이터(I) 길이: %d 개\n', length(I));
fprintf('  - 시간 데이터(t) 길이: %d 개\n', length(t));
fprintf('  - 시간 범위: %.2f 초 ~ %.2f 초 (총 %.2f 초)\n', t(1), t(end), t(end)-t(1));
fprintf('  - 최적 파라미터 tau1: %.2f 초, tau2: %.2f 초\n', params_optimal(4), params_optimal(5));
% --- 디버깅 끝 ---

if isempty(I) || length(t) < 2
    error('데이터가 너무 짧아 분석할 수 없습니다.');
end

%% 4. Sampling Frequency (fs) Calculation
% =========================================================================
fprintf('\n[단계 3: 샘플링 주파수(fs) 계산]\n');
fprintf('  - PSD 분석을 위해서는 데이터가 얼마나 자주 측정되었는지 알아야 합니다.\n');
fprintf('  - 이것을 샘플링 주파수(fs)라고 하며, 단위는 Hz (1초당 샘플 개수) 입니다.\n');

dt_vector = diff(t); % 시간 데이터 사이의 간격 계산
mean_dt = mean(dt_vector); % 평균 시간 간격
fs = 1 / mean_dt; % 샘플링 주파수는 시간 간격의 역수

% --- 디버깅 코드 ---
fprintf('\n--- 디버깅: fs 계산 과정 확인 ---\n');
fprintf('  - 첫 5개 시간 간격(dt): [%.3f, %.3f, %.3f, %.3f, %.3f] 초\n', dt_vector(1:5));
fprintf('  - 평균 시간 간격(mean_dt): %.4f 초\n', mean_dt);
fprintf('  - 계산된 샘플링 주파수(fs = 1/mean_dt): %.2f Hz\n', fs);
fprintf('>> 성공: 샘플링 주파수가 %.2f Hz로 확정되었습니다.\n', fs);
% --- 디버깅 끝 ---

%% 5. PSD Analysis and Visualization
% =========================================================================
fprintf('\n[단계 4: PSD 분석 및 시각화]\n');
fprintf('  - 이제 준비된 데이터와 fs를 사용하여 주파수 분석을 수행합니다.\n');

fig = figure('Position', [100, 100, 1600, 900], 'Name', sprintf('PSD Analysis - %s', result_key));

% --- 분석 1: Periodogram ---
fprintf('\n  [분석 1: Periodogram - 신호의 평균적인 주파수 특성]\n');
fprintf('    - 이 그래프는 신호 전체에 걸쳐 어떤 주파수 성분이 강한지를 보여줍니다.\n');
fprintf('    - 가로축(Frequency)은 진동 속도, 세로축(Power)은 해당 속도의 에너지 크기를 의미합니다.\n');

ax1 = subplot(2, 1, 1);
[pxx, f] = periodogram(I, [], [], fs, 'power');
plot(f, pxx);
set(gca, 'XScale', 'log', 'YScale', 'log'); % 로그 스케일로 넓은 범위를 관찰
grid on;
title('1. Power Spectral Density (Periodogram)');
xlabel('Frequency (Hz) [로그 스케일]');
ylabel('Power / Frequency');

% --- 디버깅 코드 ---
fprintf('\n--- 디버깅: Periodogram 결과 확인 ---\n');
[max_power, max_idx] = max(pxx);
dominant_freq = f(max_idx);
fprintf('  - 가장 강한 에너지를 가진 주파수: %.4f Hz\n', dominant_freq);
fprintf('    (이 주파수에 해당하는 시간 주기: 1/%.4f = %.2f 초)\n', dominant_freq, 1/dominant_freq);
% --- 디버깅 끝 ---

fprintf('\n    - 그래프에 피팅된 tau 값에 해당하는 주파수(f=1/tau)를 수직선으로 표시합니다.\n');
fprintf('    - 만약 전류 신호의 에너지(그래프의 봉우리)가 이 수직선 근처에 있다면, 모델이 데이터의 특성을 잘 반영한 것입니다.\n');
hold(ax1, 'on');
tau1 = params_optimal(4);
tau2 = params_optimal(5);
f1 = 1 / tau1;
f2 = 1 / tau2;

xline(ax1, f1, '--r', sprintf('f_1=1/\\tau_1 = %.3f Hz', f1), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 2, 'FontSize', 10);
xline(ax1, f2, '--g', sprintf('f_2=1/\\tau_2 = %.4f Hz', f2), 'LabelVerticalAlignment', 'top', 'LineWidth', 2, 'FontSize', 10);
hold(ax1, 'off');

% --- 디버깅 코드 ---
fprintf('\n--- 디버깅: Tau와 주파수 관계 확인 ---\n');
fprintf('  - tau1 = %.2f 초 -> f1 = 1/%.2f = %.3f Hz\n', tau1, tau1, f1);
fprintf('  - tau2 = %.2f 초 -> f2 = 1/%.2f = %.4f Hz\n', tau2, tau2, f2);
fprintf('    >> Periodogram 그래프에서 빨간선(f1)과 초록선(f2) 주변에 에너지 봉우리가 있는지 확인해보세요.\n');
% --- 디버깅 끝 ---

% --- 분석 2: Spectrogram ---
fprintf('\n  [분석 2: Spectrogram - 시간의 흐름에 따른 주파수 특성 변화]\n');
fprintf('    - 이 그래프는 시간에 따라(가로축) 어떤 주파수 성분(세로축)이 강했는지(색상)를 보여줍니다.\n');
fprintf('    - 예를 들어, 주행 초반에만 노란색이 위쪽(고주파)에 집중된다면, 초반에 빠른 변화가 많았다는 의미입니다.\n');

subplot(2, 1, 2);
pspectrum(I, t, 'spectrogram', 'TimeResolution', 60, 'OverlapPercent', 50);
title('2. Time-Frequency Analysis (Spectrogram)');

sgtitle(sprintf('PSD 분석 결과: %s', strrep(result_key, '_', '\_')), 'FontSize', 16, 'FontWeight', 'bold');

%% 6. Save Plot
% =========================================================================
fprintf('\n[단계 5: 결과 저장]\n');
fig_filename = sprintf('%s_PSD_Analysis_Debug.fig', result_key);
fig_path = fullfile(save_dir, fig_filename);
savefig(fig_path);
fprintf('>> 성공: PSD 분석 그래프를 저장했습니다: %s\n', fig_path);

% close(fig);

fprintf('\n=== PSD 분석 완료 ===\n');