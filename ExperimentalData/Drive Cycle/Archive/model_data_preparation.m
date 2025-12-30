% model_data_preparation.m 실행 환경을 가정합니다.

%% 1. X_Train 데이터 표준화 (Z-Score Normalization)
% X_Train은 이미 double 배열이므로, 추가 변환 없이 바로 zscore 적용

X_Train_Normalized = zscore(X_Train);

fprintf('X_Train 데이터 표준화 완료. (Z-Score)\n');

%% 2. K-Means 클러스터링을 통한 운전 패턴 유형화
% 이 단계는 모델 학습이 실패한 원인(Y 상수 문제)과는 별개로, 
% 현장 DOE 격차 진단을 위해 필요한 과정입니다.

max_k = 10;
sum_D = zeros(max_k, 1);
for k = 1:max_k
    % 10번 반복하여 안정적인 결과 사용 (X_Train_Normalized 사용)
    [~, ~, sumd] = kmeans(X_Train_Normalized, k, 'Replicates', 10, 'MaxIter', 500); 
    sum_D(k) = sum(sumd);
end

% Elbow Method 결과 시각화
figure('Name', 'Elbow Method for Optimal K (FR Drive Cycle Clustering)');
plot(1:max_k, sum_D, 'bo-', 'LineWidth', 2);
xlabel('Number of Clusters (K)');
ylabel('Sum of Squared Distances (Within-Cluster)');
title('DOE Gap Diagnosis: Identify Distinct FR Operation Types');
grid on;

fprintf('Elbow Method 클러스터링 및 시각화 완료. 그래프를 확인하십시오.\n');