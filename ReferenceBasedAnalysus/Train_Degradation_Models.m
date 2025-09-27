% Train_Degradation_Models.m 스크립트의 D 섹션을 아래 코드로 대체합니다.
% X_Train과 Y_Train은 double 배열로 작업 공간에 존재한다고 가정합니다.

%% D. LLI 및 LAM Random Forest 모델 학습

% LLI/LAM Y 라벨 추출 (숫자 인덱스 사용으로 오류 해결)
Y_LLI = Y_Train(:, 1); % LLI_Loss_Rate_Label (첫 번째 열)
Y_LAM = Y_Train(:, 2); % LAM_Loss_Rate_Label (두 번째 열)
Y_CL = Y_Train(:, 3);  % CL_DCIR_Rate_Label (세 번째 열)

fprintf('\n=== LLI 및 LAM Random Forest 모델 학습 시작 ===\n');

% --- 1. LLI (리튬 손실) 예측 모델 학습 ---
rng(3); % 재현성 확보
RF_LLI_Model = fitrensemble(X_Train_Normalized, Y_LLI, ...
                           'Method', 'Bag', ...
                           'NumLearningCycles', 50, ...
                           'Learners', templateTree('MinLeafSize', 5)); 

CV_RF_LLI_Model = crossval(RF_LLI_Model, 'KFold', 5);
LLI_Predicted = kfoldPredict(CV_RF_LLI_Model);

LLI_RMSE = sqrt(mean((LLI_Predicted - Y_LLI).^2));
LLI_R_squared = 1 - (sum((LLI_Predicted - Y_LLI).^2) / sum((Y_LLI - mean(Y_LLI)).^2));

fprintf('1. LLI (Loss of Lithium Inventory) Model 학습 완료.\n');
fprintf('   RMSE: %.6f, R-Squared: %.4f\n', LLI_RMSE, LLI_R_squared);


% --- 2. LAM (활물질 손실) 예측 모델 학습 ---
rng(4); % 재현성 확보
RF_LAM_Model = fitrensemble(X_Train_Normalized, Y_LAM, ...
                           'Method', 'Bag', ...
                           'NumLearningCycles', 50, ...
                           'Learners', templateTree('MinLeafSize', 5));

CV_RF_LAM_Model = crossval(RF_LAM_Model, 'KFold', 5);
LAM_Predicted = kfoldPredict(CV_RF_LAM_Model);

LAM_RMSE = sqrt(mean((LAM_Predicted - Y_LAM).^2));
LAM_R_squared = 1 - (sum((LAM_Predicted - Y_LAM).^2) / sum((Y_LAM - mean(Y_LAM)).^2));

fprintf('2. LAM (Loss of Active Material) Model 학습 완료.\n');
fprintf('   RMSE: %.6f, R-Squared: %.4f\n', LAM_RMSE, LAM_R_squared);


% --- 3. CL (저항 열화) 예측 모델 학습 (참고용) ---
rng(5); 
RF_CL_Model = fitrensemble(X_Train_Normalized, Y_CL, ...
                           'Method', 'Bag', ...
                           'NumLearningCycles', 50, ...
                           'Learners', templateTree('MinLeafSize', 5));

CV_RF_CL_Model = crossval(RF_CL_Model, 'KFold', 5);
CL_Predicted = kfoldPredict(CV_RF_CL_Model);

CL_R_squared = 1 - (sum((CL_Predicted - Y_CL).^2) / sum((Y_CL - mean(Y_CL)).^2));

fprintf('3. CL (Conduction Loss) Model 학습 완료 (참고용).\n');
fprintf('   R-Squared: %.4f\n', CL_R_squared);


%% E. 성능 출력 및 분석 (가장 중요)
fprintf('\n=========================================\n');
fprintf('모델 학습 최종 성능 평가:\n');
fprintf('=========================================\n');
fprintf('1. LLI Model R-Squared: %.4f\n', LLI_R_squared);
fprintf('2. LAM Model R-Squared: %.4f\n', LAM_R_squared);
fprintf('3. CL Model R-Squared (참고): %.4f\n', CL_R_squared);
fprintf('=========================================\n');