%% ML_DeepLearning_v6.m
% =====================================================================
% Phase 4D: 1D-CNN 하이브리드 모델 튜닝 (SOH 추정)
% - 변경점 1 (Target Scaling): 회귀 목표인 SOH를 0~1 사이로 정규화 (SOH/100)
% - 변경점 2 (Network Simplification): 랩 데이터의 작은 샘플 수에 맞춰
%   Conv 레이어를 1단으로 줄이고 파라미터 수를 대폭 축소 (과적합 방지)
% =====================================================================
clear; clc; close all;
rng(42); 

%% 1. 데이터 로드 및 전처리
d_pcc = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4\PCC_v5b_data.mat');

X_chg = d_pcc.dQ_chg_all;  % [N x 11]
c_eff = d_pcc.data_CrateNum;
y_soh = d_pcc.data_SOH;
cell_ids = d_pcc.data_CellID;

valid_chg = ~any(isnan(X_chg), 2) & ~isnan(c_eff) & ~isnan(y_soh);
X_c = X_chg(valid_chg, :); 
c_v = c_eff(valid_chg); 
y_v = y_soh(valid_chg);
cid_v = cell_ids(valid_chg);
N = size(X_c, 1);

% Feature Standardization
mu_X = mean(X_c, 1); sig_X = std(X_c, 0, 1);
mu_C = mean(c_v);    sig_C = std(c_v);

X_norm = (X_c - mu_X) ./ sig_X;
c_norm = (c_v - mu_C) ./ sig_C;

% *** TARGET SCALING (SOH 0~1) ***
y_norm = y_v / 100.0;

% 딥러닝 차원 변환
X_dl = reshape(X_norm', [1, 11, 1, N]); 
C_dl = c_norm; 
y_dl = y_norm;

%% 2. 1D-CNN + C_eff 하이브리드 아키텍처 재설계 (단순화)
lgraph = layerGraph();

% 분기 1: 단순화된 1D-CNN 
layers_CNN = [
    imageInputLayer([1 11 1], 'Name', 'input_dQ', 'Normalization', 'none')
    convolution2dLayer([1 3], 8, 'Padding', 'same', 'Name', 'conv1') % 필터 수 축소(16->8)
    batchNormalizationLayer('Name', 'bn1')
    reluLayer('Name', 'relu1')
    % globalAveragePooling 대신 바로 flatten 하여 공간 정보 유지
    flattenLayer('Name', 'flatten')  
    fullyConnectedLayer(8, 'Name', 'fc_cnn_out')
    reluLayer('Name', 'relu_cnn')
];
lgraph = addLayers(lgraph, layers_CNN);

% 분기 2: C_eff 스칼라 입력
layers_C = [
    featureInputLayer(1, 'Name', 'input_Ceff', 'Normalization', 'none')
    fullyConnectedLayer(4, 'Name', 'fc_c1')
    reluLayer('Name', 'relu_c1')
];
lgraph = addLayers(lgraph, layers_C);

% 병합 및 최종 출력
layers_Concat = [
    concatenationLayer(1, 2, 'Name', 'concat')
    fullyConnectedLayer(8, 'Name', 'fc_merged')
    reluLayer('Name', 'relu_merged')
    fullyConnectedLayer(1, 'Name', 'output')
    regressionLayer('Name', 'regression')
];
lgraph = addLayers(lgraph, layers_Concat);

lgraph = connectLayers(lgraph, 'relu_cnn', 'concat/in1');
lgraph = connectLayers(lgraph, 'relu_c1', 'concat/in2');


%% 3. Group K-Fold 훈련 및 평가
u_cells = unique(cid_v); K = 5;
cv_idx = zeros(N, 1); 
groups = ceil(linspace(0.01, K, length(u_cells)));
for k=1:K, tc=u_cells(groups==k); for ti=1:length(tc), cv_idx(strcmp(cid_v, tc{ti})) = k; end; end

y_pred_dl_norm = nan(N, 1);

options = trainingOptions('adam', ...
    'MiniBatchSize', 32, ...
    'MaxEpochs', 300, ...
    'InitialLearnRate', 0.01, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.2, ...
    'LearnRateDropPeriod', 150, ...
    'Shuffle', 'every-epoch', ...
    'Verbose', false, ...
    'Plots', 'none');

fprintf('=== 1D-CNN Tuned (Phase 4D) Group 5-Fold Training ===\n');
for k = 1:K
    trn = (cv_idx~=k); tst = (cv_idx==k);
    
    ds_X_trn = arrayDatastore(X_dl(:,:,:,trn), 'IterationDimension', 4);
    ds_C_trn = arrayDatastore(C_dl(trn), 'IterationDimension', 1);
    ds_Y_trn = arrayDatastore(y_dl(trn), 'IterationDimension', 1);
    ds_train = combine(ds_X_trn, ds_C_trn, ds_Y_trn);
    
    ds_X_tst = arrayDatastore(X_dl(:,:,:,tst), 'IterationDimension', 4);
    ds_C_tst = arrayDatastore(C_dl(tst), 'IterationDimension', 1);
    ds_test  = combine(ds_X_tst, ds_C_tst);
    
    net_cv = trainNetwork(ds_train, lgraph, options);
    
    y_pred_dl_norm(tst) = predict(net_cv, ds_test);
end

% 역정규화 (0~1 -> 원래 SOH)
y_pred_dl_real = y_pred_dl_norm * 100.0;

% 전체 성능
rmse_total = sqrt(mean((y_v - y_pred_dl_real).^2));
r2_total = 1 - (sum((y_v - y_pred_dl_real).^2) / sum((y_v - mean(y_v)).^2));
fprintf('\n>> Final Tuned 1D-CNN CV Result: R2: %.3f | RMSE: %.3f\n', r2_total, rmse_total);

%% 4. 최종 모델 훈련 및 Parity Plot
fprintf('\nTraining final full model for field application...\n');
ds_X_full = arrayDatastore(X_dl, 'IterationDimension', 4);
ds_C_full = arrayDatastore(C_dl, 'IterationDimension', 1);
ds_Y_full = arrayDatastore(y_dl, 'IterationDimension', 1);
ds_train_full = combine(ds_X_full, ds_C_full, ds_Y_full);

final_net = trainNetwork(ds_train_full, lgraph, options);

% 최종 저장
DL_Model = struct();
DL_Model.net = final_net;
DL_Model.mu_X = mu_X; DL_Model.sig_X = sig_X;
DL_Model.mu_C = mu_C; DL_Model.sig_C = sig_C;

save(fullfile(fileparts(mfilename('fullpath')), 'ML_DL_Model_v6.mat'), 'DL_Model', 'y_pred_dl_real', 'y_v');
fprintf('>> Saved: ML_DL_Model_v6.mat\n');

% Parity Plot
fig_dl = figure('Position', [100,100,600,600], 'Name', '1D-CNN Tuned Parity Plot');
hold on; grid on; axis square;
scatter(y_v, y_pred_dl_real, 40, c_v, 'filled', 'MarkerEdgeColor', 'k');
colormap(turbo); cb=colorbar; ylabel(cb,'C-rate'); clim([0.1 3]);
plot([75 102], [75 102], 'k--', 'LineWidth', 2);
xlabel('True SOH (%)', 'FontWeight', 'bold'); 
ylabel('Predicted DL SOH (%)', 'FontWeight', 'bold');
title(sprintf('1D-CNN Hybrid Model (Tuned)\nR^2 = %.3f, RMSE = %.3f', r2_total, rmse_total), 'FontSize', 14, 'FontWeight','bold');
xlim([75 102]); ylim([75 102]);

saveas(fig_dl, fullfile(fileparts(mfilename('fullpath')), 'CNN_Parity_Plot.fig'));
fprintf('>> Saved: CNN_Parity_Plot.fig\n');
