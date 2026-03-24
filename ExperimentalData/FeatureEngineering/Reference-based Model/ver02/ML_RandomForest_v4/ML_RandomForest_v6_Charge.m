%% ML_RandomForest_v6_Charge.m
% =====================================================================
% Phase 5: NCM 단편화 충전 용량 (Fragmented Charge Capacity) ML 모델 훈련
% - v5_Charge (Equi-capacity 11 Segments) 피처 사용
% - 방전 모델 X, 오직 충전 모델만
% - 전체 상관성이 분산되어 있으므로 (PCC 최대 +0.38), 특정 세그먼트 소수 선택(Track A/B)이 
%   아닌 전체 11개 세그먼트 중 주요 조합 탐색 (PCC Top 4, Top 7, All 11)
% - RF, SVM, GPR 3가지 모델 생성 및 CV 성능 평가
% =====================================================================
clear; clc; close all;
rng(42);

%% 1. 데이터 로드
base_dir  = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
                     'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02', 'ML_RandomForest_v4');
feat_path = fullfile(base_dir, 'RPT_FeatureLabelExtractor_v5_Charge', 'Feature_Matrix_v5_Charge.mat');
label_path = fullfile(base_dir, '..', 'RPT_FeatureLabelExtractor', 'Label_Matrix_ver02.mat');

fprintf('Loading v5 Charge Features...\n');
feat_data = load(feat_path);
label_data = load(label_path);

FeatureTable = feat_data.FeatureTable_v5_Charge;
LabelTable   = label_data.LabelTable_ver02;

% 매칭
[~, ia, ib] = intersect( ...
    strcat(FeatureTable.CellID, num2str(FeatureTable.Cycle), FeatureTable.CrateLabel), ...
    strcat(LabelTable.CellID,   num2str(LabelTable.Cycle),   LabelTable.CrateLabel));

% 11 Segments + C_eff 추출
feature_cols = 5:15; % dQ_chg_S1 ~ dQ_chg_S11
X_seg = table2array(FeatureTable(ia, feature_cols));
C_eff = FeatureTable.C_eff_chg(ia); % col 19 in FeatureTable usually, or by name

y_soh = LabelTable.SOH(ib);
cell_ids = FeatureTable.CellID(ia);

% NaN 제거
valid_idx = ~isnan(y_soh) & ~any(isnan(X_seg), 2) & ~isnan(C_eff);
X_seg    = X_seg(valid_idx, :);
C_eff    = C_eff(valid_idx);
y_soh    = y_soh(valid_idx);
cell_ids = cell_ids(valid_idx);

%% 2. 모델 구성 설정
% 이전 PCC 결과 기반 (Seg 10, 9, 11, 6 이 가장 높았음)
% Top 4: [6, 9, 10, 11]
% Top 7: [4, 5, 6, 7, 8, 9, 10]
% All: 1:11
conf.Top4.segs = [6, 9, 10, 11];
conf.Top7.segs = [4, 5, 6, 7, 8, 9, 10, 11]; % 8장으로 수정하여 Top 8
conf.All.segs  = 1:11;

tracks = {'Top4', 'Top7', 'All'};
model_types = {'RF', 'SVM', 'GPR'};

Final_Models = struct();
Results_CV   = struct();

% CV 그룹 (Group K-Fold by Cell)
u_cells = unique(cell_ids); 
K = 5;
cv_indices = zeros(size(y_soh)); 
groups = ceil(linspace(0.01, K, length(u_cells)));
rng(42); groups = groups(randperm(length(groups))); % Shuffle
for k=1:K
    tc=u_cells(groups==k); 
    for ti=1:length(tc)
        cv_indices(strcmp(cell_ids, tc{ti})) = k; 
    end
end

%% 3. 최종 모델 훈련 및 CV 평가 루프
fprintf('\n=== Phase 5: Model Training (Fragmented Charge) ===\n');

for t_idx = 1:length(tracks)
    trk = tracks{t_idx};
    segs = conf.(trk).segs;
    
    Xc = [X_seg(:, segs), C_eff]; % [선택된 세그먼트들, C_eff]
    
    fprintf('\n>> [%s] Features: Segs [%s] + C_eff (N=%d)\n', ...
        trk, num2str(segs), size(Xc,1));
    
    for mt_idx = 1:length(model_types)
        mtype = model_types{mt_idx};
        y_pred = nan(size(y_soh));
        
        % CV 평가
        for k = 1:K
            trn = (cv_indices ~= k); 
            tst = (cv_indices == k);
            
            switch mtype
                case 'RF'
                    model = TreeBagger(50, Xc(trn,:), y_soh(trn), 'Method','regression','MinLeafSize',5);
                    y_pred(tst) = predict(model, Xc(tst,:));
                case 'SVM'
                    model = fitrsvm(Xc(trn,:), y_soh(trn), 'KernelFunction','gaussian','Standardize',true);
                    y_pred(tst) = predict(model, Xc(tst,:));
                case 'GPR'
                    model = fitrgp(Xc(trn,:), y_soh(trn), 'KernelFunction','squaredexponential','Standardize',true);
                    y_pred(tst) = predict(model, Xc(tst,:));
            end
        end
        
        % CV Metric 계산
        rmse = sqrt(mean((y_soh - y_pred).^2, 'omitnan'));
        SST = sum((y_soh - mean(y_soh,'omitnan')).^2,'omitnan');
        SSE = sum((y_soh - y_pred).^2,'omitnan');
        r2 = 1 - (SSE / SST);
        
        Results_CV.(trk).(mtype).R2 = r2;
        Results_CV.(trk).(mtype).RMSE = rmse;
        Results_CV.(trk).(mtype).y_true = y_soh;
        Results_CV.(trk).(mtype).y_pred = y_pred;
        
        fprintf('  SOH | %3s | R2=%5.3f, RMSE=%5.3f\n', mtype, r2, rmse);
        
        % === 최종 전체 데이터로 모델 재훈련 및 저장 ===
        switch mtype
            case 'RF'
                final_mdl = TreeBagger(50, Xc, y_soh, 'Method','regression');
            case 'SVM'
                final_mdl = fitrsvm(Xc, y_soh, 'KernelFunction','gaussian','Standardize',true);
            case 'GPR'
                final_mdl = fitrgp(Xc, y_soh, 'KernelFunction','squaredexponential','Standardize',true);
        end
        
        Final_Models.(trk).(mtype).model = final_mdl;
        Final_Models.(trk).(mtype).features = segs;
        Final_Models.(trk).(mtype).feature_names = [arrayfun(@(i) sprintf('Seg%d',i), segs, 'UniformOutput',false), {'C_eff'}];
    end
end

%% 4. 결과 저장
saveDir = fullfile(base_dir, 'ML_Models_v6_Charge.mat');
save(saveDir, 'Final_Models', 'Results_CV');
fprintf('\n>> Saved: ML_Models_v6_Charge.mat\n');

%% [선택] 최고 성능 모델 Parity Plot (All GPR/RF)
fig = figure('Position',[100,100,1200,400],'Name','Parity Plot: v6 Fragmented Charge');
y_true_A = Results_CV.All.RF.y_true;
y_pred_A = Results_CV.All.RF.y_pred;
y_true_B = Results_CV.All.GPR.y_true;
y_pred_B = Results_CV.All.GPR.y_pred;

subplot(1,2,1); hold on; grid on;
scatter(y_true_A, y_pred_A, 50, C_eff, 'filled', 'MarkerEdgeColor','k');
colormap(gca, turbo); cb=colorbar; ylabel(cb,'C-rate (C_{eff})'); clim([0.1 3]);
plot([75 102], [75 102], 'r--', 'LineWidth',2);
xlabel('True SOH (%)'); ylabel('Predicted SOH (%)');
r2_a = Results_CV.All.RF.R2;
title(sprintf('RF (All 11 Segs + C_{eff})\nSOH R^2 = %.3f', r2_a));
axis square; xlim([75 102]); ylim([75 102]);

subplot(1,2,2); hold on; grid on;
scatter(y_true_B, y_pred_B, 50, C_eff, 'filled', 'MarkerEdgeColor','k');
colormap(gca, turbo); cb=colorbar; ylabel(cb,'C-rate (C_{eff})'); clim([0.1 3]);
plot([75 102], [75 102], 'r--', 'LineWidth',2);
xlabel('True SOH (%)'); ylabel('Predicted SOH (%)');
r2_b = Results_CV.All.GPR.R2;
title(sprintf('GPR (All 11 Segs + C_{eff})\nSOH R^2 = %.3f', r2_b));
axis square; xlim([75 102]); ylim([75 102]);

saveas(fig, fullfile(base_dir, 'Parity_Plot_v6_Charge.fig'));
saveas(fig, fullfile(base_dir, 'Parity_Plot_v6_Charge.png'));
fprintf('>> Saved Figures: Parity_Plot_v6_Charge.fig / .png\n');
