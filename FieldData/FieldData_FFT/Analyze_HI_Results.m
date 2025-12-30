% Analyze_Total_Results.m

% 1. 데이터 로드
load('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_FFT\Operational_Impedance_Spectrum_OldData\Advanced_HI_Results_Total_2023_2024_2025.mat');

% 2. fillmissing으로 NaN 처리 (이전 스크립트에서 이미 처리했다면 생략 가능)
total_table = fillmissing(total_results_table, 'previous', 'DataVariables', @isnumeric);

% 3. 온도와 각 HI의 관계 모델링 (선형 회귀)
mdl_cdrai = fitlm(total_table.Avg_Temp, total_table.CDRAI);
mdl_por = fitlm(total_table.Avg_Temp, total_table.POR);

disp(mdl_cdrai); % 모델 요약 확인 (p-value, R-squared 등)
disp(mdl_por);

% 4. 온도로 인한 영향 예측 및 제거
temp_effect_cdrai = predict(mdl_cdrai, total_table.Avg_Temp);
total_table.CDRAI_Compensated = total_table.CDRAI - temp_effect_cdrai;

temp_effect_por = predict(mdl_por, total_table.Avg_Temp);
total_table.POR_Compensated = total_table.POR - temp_effect_por;

% 5. 온도 보상 전/후 비교 시각화
figure('Position', [100, 100, 1200, 600]);
subplot(2,1,1);
plot(total_table.Date, total_table.CDRAI, 'Color', [0.8 0.8 1]);
hold on;
plot(total_table.Date, movmean(total_table.CDRAI, 30), 'b', 'LineWidth', 2);
title('Original CDRAI (with Seasonality)');
legend('Daily', '30-day Trend'); grid on;

subplot(2,1,2);
plot(total_table.Date, total_table.CDRAI_Compensated, 'Color', [1 0.8 0.8]);
hold on;
plot(total_table.Date, movmean(total_table.CDRAI_Compensated, 30), 'r', 'LineWidth', 2);
title('Temperature-Compensated CDRAI (Pure Degradation Trend)');
legend('Daily Compensated', '30-day Trend'); grid on;

% Analyze_Total_Results.m 스크립트에 추가

% --- 4. 다중 선형 회귀 분석: 온도와 시간의 영향 분리 ---

% 1. '시간' 변수 생성 (첫 날로부터 경과일)
total_table.Time_days = days(total_table.Date - total_table.Date(1));

fprintf('\n--- Multiple Linear Regression Results ---\n');

% 2. CDRAI에 대한 다중 회귀 모델링 (CDRAI ~ 1 + Avg_Temp + Time_days)
mdl_cdrai_multi = fitlm(total_table, 'CDRAI ~ Avg_Temp + Time_days');
disp('--- Model: CDRAI vs. Temp + Time ---');
disp(mdl_cdrai_multi);

% 3. POR에 대한 다중 회귀 모델링 (POR ~ 1 + Avg_Temp + Time_days)
mdl_por_multi = fitlm(total_table, 'POR ~ Avg_Temp + Time_days');
disp('--- Model: POR vs. Temp + Time ---');
disp(mdl_por_multi);