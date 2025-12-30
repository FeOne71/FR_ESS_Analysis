%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build_SOH_Model_v4_Trend.m
%
% [Version 4: 추세선 기반 모델링]
% - 문제점: 일별 R1s 데이터의 노이즈로 인해 예측된 SOH가 비현실적으로 변동함.
% - 해결책: 이동 평균(Moving Average)을 이용해 R1s 데이터의 부드러운 '추세선'을
%           계산하고, 이 추세선을 기반으로 SOH 모델을 구축 및 예측함.
% - 기대효과: 물리적으로 타당한, 부드럽게 감소하는 SOH 예측 곡선 획득.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('====================================================\n');
fprintf('SOH 예측 모델 구축 스크립트 (v4, 추세선 기반) 시작\n');
fprintf('====================================================\n\n');

%% 1. 데이터 로드
% =========================================================================
fprintf('--- [1] 입력 데이터 로드 ---\n');
base_path = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_R1s';
results_file = fullfile(base_path, 'Arrhenius_Analysis_Results', 'Arrhenius_Analysis_Results.mat');
try
    load(results_file);
    fprintf('성공: %s 파일을 로드했습니다.\n', results_file);
catch ME
    error('Arrhenius_Analysis_Results.mat 파일을 찾을 수 없습니다.');
end

all_dates = arrhenius_results.allDates;
all_r1s_corrected = arrhenius_results.r1s_corrected_exp;


%% 2. R1s 데이터의 '추세선' 계산 (핵심 수정사항)
% =========================================================================
fprintf('\n--- [2] R1s 데이터의 추세선 계산 ---\n');

% 이동 평균 윈도우 크기 (일). 값이 클수록 더 부드러운 곡선이 됨.
window_days = 30;
fprintf('이동 평균 윈도우: %d일\n', window_days);

% movmean 함수를 사용하여 추세선 계산
r1s_trend = movmean(all_r1s_corrected, window_days, 'omitnan');

% R1s 원본 데이터와 추세선을 테이블로 묶기
r1s_table = table(all_dates, all_r1s_corrected, r1s_trend, 'VariableNames', {'Date', 'R1s_Daily', 'R1s_Trend'});


%% 3. RPT "Ground Truth" 데이터 정의 및 추세선과 매칭
% =========================================================================
fprintf('\n--- [3] RPT 데이터를 R1s ''추세선''과 매칭 ---\n');
dates_rpt = [datetime('2021-06-07'); datetime('2023-10-16'); datetime('2025-07-11')];
capacity_rpt = [59.8544; 51.8008; 50.3143];
soh_rpt = (capacity_rpt / capacity_rpt(1)) * 100;
ground_truth_table = table(dates_rpt, soh_rpt, 'VariableNames', {'Date', 'SOH'});

model_data = table();
for i = 1:height(ground_truth_table)
    rpt_date = ground_truth_table.Date(i);
    [~, idx] = min(abs(r1s_table.Date - rpt_date));
    
    % 중요: 일별 R1s가 아닌, 해당 날짜의 '추세선 R1s' 값을 사용
    matched_r1s_trend = r1s_table.R1s_Trend(idx);
    
    fprintf('매칭 %d: %s의 R1s 추세선 값 = %.6f Ohm (SOH: %.2f %%)\n', i, datestr(rpt_date), matched_r1s_trend, soh_rpt(i));
    
    temp_row = ground_truth_table(i,:);
    temp_row.R1s_Trend = matched_r1s_trend;
    model_data = [model_data; temp_row];
end


%% 4. SOH 예측 모델 구축 (추세선 기반)
% =========================================================================
fprintf('\n--- [4] 추세선 기반 SOH 예측 모델 구축 ---\n');
mdl_soh = fitlm(model_data.R1s_Trend, model_data.SOH);
disp(mdl_soh);


%% 5. 전체 기간 SOH 예측 및 시각화
% =========================================================================
fprintf('\n--- [5] 전체 기간 SOH 예측 및 결과 시각화 ---\n');

% 중요: 모델 입력으로 전체 기간의 '추세선 R1s' 데이터를 사용
all_soh_predicted = predict(mdl_soh, r1s_table.R1s_Trend);

fprintf('예측된 SOH 최소값: %.2f %%\n', min(all_soh_predicted));
fprintf('예측된 SOH 최대값: %.2f %%\n\n', max(all_soh_predicted));

% --- 최종 결과 시각화 ---
figure('Name', 'SOH Prediction using R1s Trend', 'Position', [100, 100, 1400, 800]);

% 1. 상단 그래프: R1s 데이터와 추세선 비교
subplot(2,1,1);
plot(r1s_table.Date, r1s_table.R1s_Daily * 1000, '.', 'Color', [0.7 0.7 0.7], 'DisplayName', '일별 R1s (노이즈 포함)');
hold on;
plot(r1s_table.Date, r1s_table.R1s_Trend * 1000, 'r-', 'LineWidth', 2.5, 'DisplayName', sprintf('%d일 이동평균 추세선', window_days));
scatter(model_data.Date, model_data.R1s_Trend * 1000, 150, 'b', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', '모델 학습에 사용된 R1s 지점');
title('온도 보정된 R1s의 일별 데이터와 추세선', 'FontSize', 14);
ylabel('R1s_{corrected} (mΩ)');
legend('show');
grid on;

% 2. 하단 그래프: 최종 SOH 예측 결과
subplot(2,1,2);
plot(r1s_table.Date, all_soh_predicted, 'b-', 'LineWidth', 2.5, 'DisplayName', '추세선 기반 SOH 예측');
hold on;
scatter(model_data.Date, model_data.SOH, 150, 'r', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', '실제 RPT 측정값 (Ground Truth)');

title('최종 SOH 예측 결과 (추세선 기반)', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('날짜');
ylabel('SOH (%)');
legend('Location', 'southwest');
grid on;
ylim([min(all_soh_predicted)-2, 101]);

fprintf('====================================================\n');
fprintf('분석이 완료되었습니다. 2개의 그래프를 확인하세요.\n');
fprintf('====================================================\n');