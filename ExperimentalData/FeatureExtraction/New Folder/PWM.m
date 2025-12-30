% 1. 저장된 이벤트 데이터베이스 로드
load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureExtraction\Results\Filtered_Event_Database.mat');

% 2. 예시: SOC90, 0 사이클의 첫 번째 안정적인 이벤트 선택
event_data = Filtered_Events_DB.SOC90.Cycle_0{1};
I_signal = event_data.Current;
V_signal = event_data.Voltage;
Fs = 10; % 샘플링 주파수 (10 Hz)

% --- [핵심 개선] DC 성분 제거 ---
% 각 신호에서 평균값을 빼주어 AC 변동 성분만 남깁니다.
I_ac = I_signal - mean(I_signal);
V_ac = V_signal - mean(V_signal);

% 3. 전달 함수(임피던스) 추정 (tfestimate 사용)
% 'Window' 옵션은 Leakage 현상을 줄여주는 중요한 역할을 합니다.
[Z_est, f_axis] = tfestimate(I_ac, V_ac, hanning(length(I_ac)), [], [], Fs);

% 4. 특징 추출: 특정 주파수에서의 임피던스 크기(|Z|) 추출
[~, idx_1Hz] = min(abs(f_axis - 1.0));
[~, idx_0_5Hz] = min(abs(f_axis - 0.5));
[~, idx_0_1Hz] = min(abs(f_axis - 0.1));

Z_1Hz = abs(Z_est(idx_1Hz));
Z_0_5Hz = abs(Z_est(idx_0_5Hz));
Z_0_1Hz = abs(Z_est(idx_0_1Hz));

fprintf('개선된 유사 임피던스 특징 (tfestimate 사용):\n');
fprintf('  |Z(1.0 Hz)| = %.4f Ohm  (%.1f mOhm)\n', Z_1Hz, Z_1Hz*1000);
fprintf('  |Z(0.5 Hz)| = %.4f Ohm  (%.1f mOhm)\n', Z_0_5Hz, Z_0_5Hz*1000);
fprintf('  |Z(0.1 Hz)| = %.4f Ohm  (%.1f mOhm)\n', Z_0_1Hz, Z_0_1Hz*1000);

% 시각화
figure;
plot(f_axis, abs(Z_est));
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)');
ylabel('|Impedance| (Ohm)');
title('전달 함수 추정 기반 유사 임피던스 스펙트럼');
grid on;
xlim([0.01 Fs/2]);