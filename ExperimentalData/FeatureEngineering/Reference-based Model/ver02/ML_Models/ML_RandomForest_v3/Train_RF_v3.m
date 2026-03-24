% Train_RF_v3.m  — [3.4~4.0]V 11세그먼트 + 부분 세그먼트 지원 (Surrogate Splits)
%
% [v3] 핵심 차이점:
%   1. Feature_Matrix_v3.mat 사용 (Extract_Features_v3.m으로 생성)
%   2. 11 dQ 세그먼트 (chg 11 + dch 11 = 22 features) + C_eff 정규화용
%   3. Surrogate splits 활성화: NaN 피처 허용 → 부분 세그먼트로도 예측 가능
%   4. Label: SOH만
%   5. fold-local normalization + split chg/dch 동일 유지
%   6. NaN 행 제거 안 함 (NaN이 있는 dQ 세그먼트는 학습에서 surrogate로 처리)
%
% [실행 전] Extract_Features_v3.m을 먼저 실행하여 Feature_Matrix_v3.mat 생성 필요

clear; clc; close all;

base_dir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
                    'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
save_dir = fullfile(base_dir, 'ML_RandomForest_v3');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% 1. 데이터 로드 (v3 Feature Matrix)
fprintf('------------------------------------------------------------\n');
fprintf(' [1] Loading Feature_Matrix_v3 + Label_Matrix_ver02\n');
fprintf('------------------------------------------------------------\n');

feat_data  = load(fullfile(save_dir, 'Feature_Matrix_v3.mat'));
label_data = load(fullfile(base_dir, 'RPT_FeatureLabelExtractor', 'Label_Matrix_ver02.mat'));

FT = feat_data.FeatureTable_v3;
LT = label_data.LabelTable_ver02;
V_bounds = feat_data.V_bounds;
all_feature_names = feat_data.feature_names;  % 24개: dQ_chg_1~11, dQ_dch_1~11, C_eff_chg, C_eff_dch
num_segments = length(V_bounds) - 1;

% 라벨: SOH만
label_names = {'SOH'};

% Row alignment (CellID + Cycle + CrateLabel)
[~, ia, ib] = intersect( ...
    strcat(FT.CellID, num2str(FT.Cycle), FT.CrateLabel), ...
    strcat(LT.CellID, num2str(LT.Cycle), LT.CrateLabel));
FT = FT(ia, :);
LT = LT(ib, :);

% dQ 피처만 선택 (C_eff는 정규화용으로만, 피처에서 제외)
dq_feat_names = all_feature_names(1:num_segments*2);  % dQ_chg_1~11 + dQ_dch_1~11
feat_idx = 1:num_segments*2;
X = table2array(FT(:, dq_feat_names));
Y = LT.SOH;
cellIDs = FT.CellID;
crateLabels = FT.CrateLabel;

n_feats = size(X, 2);
N = size(X, 1);

% C_eff 컬럼 (정규화 그룹 결정용)
C_eff_chg_all = FT.C_eff_chg;
C_eff_dch_all = FT.C_eff_dch;

% NaN 행 제거하지 않음! Surrogate splits가 NaN을 처리
% 단, Y(SOH)에 NaN이 있는 행만 제거
nan_y = isnan(Y);
X = X(~nan_y, :); Y = Y(~nan_y); cellIDs = cellIDs(~nan_y);
crateLabels = crateLabels(~nan_y);
N = size(X, 1);

fprintf(' [Data] %d samples | %d features (dQ only) | %d segments\n', N, n_feats, num_segments);
fprintf(' [NaN] X에 NaN 비율: %.1f%% (surrogate splits로 처리)\n', 100*sum(isnan(X(:)))/numel(X));

% 인덱스 정의
chg_idx = 1:num_segments;            % dQ_chg_1~11
dch_idx = (num_segments+1):(2*num_segments);  % dQ_dch_1~11

% C-rate 그룹
crate_groups   = {'c01','c05','c1','c2','c3'};
crate_vals_num = [0.1, 0.5, 1.0, 2.0, 3.0];

%% 2. Group K-Fold CV 인덱스
cells_unique = unique(cellIDs);
num_cells = length(cells_unique);
K_group = min(5, num_cells);
cv_idx = zeros(N, 1);
rng(42);
shuffled = cells_unique(randperm(num_cells));
fold_assign = mod(0:num_cells-1, K_group) + 1;
for i = 1:num_cells
    cv_idx(strcmp(cellIDs, shuffled{i})) = fold_assign(i);
end
K = max(cv_idx);
fprintf(' [CV] %d folds (Group K-Fold, %d cells)\n', K, num_cells);

%% 3. HP 탐색 + CV 평가 + 최종 모델
num_trees_list = [50, 100, 200, 500];
min_leaf_list  = [1, 3, 5];

Results_RF = struct();
Results_RF.feature_names = dq_feat_names(:);
Results_RF.label_names   = label_names;
Results_RF.V_bounds      = V_bounds;

y_pred_cv = zeros(N, 1);

fprintf('\n============================================================\n');
fprintf(' [Target: SOH] (%d dQ segments, surrogate splits)\n', num_segments*2);
fprintf('============================================================\n');

% STEP A: HP 탐색
fprintf(' [A] HP Tuning via %d-Fold CV...\n', K);
best_rmse = inf; best_mt = 100; best_ml = 1;
hp_history = [];

for t = num_trees_list
    for l = min_leaf_list
        rmse_k = zeros(K, 1);
        for k = 1:K
            tr = (cv_idx ~= k); te = ~tr;
            NormS = compute_norm_stats(X(tr,:), crateLabels(tr), chg_idx, dch_idx, crate_groups);
            Xtr = apply_norm_split(X(tr,:), crateLabels(tr), NormS, chg_idx, dch_idx);
            Xte = apply_norm_split(X(te,:), crateLabels(te), NormS, chg_idx, dch_idx);
            rng(42);
            % Surrogate='on' → NaN 피처에서 대리 분기 사용
            rf = fitrensemble(Xtr, Y(tr), 'Method', 'Bag', ...
                'NumLearningCycles', t, ...
                'Learners', templateTree('MinLeafSize', l, 'Surrogate', 'on'));
            rmse_k(k) = sqrt(mean((Y(te) - predict(rf, Xte)).^2));
        end
        avg_r = mean(rmse_k);
        hp_history = [hp_history; t, l, avg_r]; %#ok<AGROW>
        if avg_r < best_rmse, best_rmse = avg_r; best_mt = t; best_ml = l; end
        fprintf('   n_trees=%3d  min_leaf=%d  CV_RMSE=%.4f\n', t, l, avg_r);
    end
end
fprintf(' → Best: n_trees=%d  min_leaf=%d  CV_RMSE=%.4f\n', best_mt, best_ml, best_rmse);

% STEP B: CV 평가 (fold-local)
fprintf(' [B] Group %d-Fold CV Evaluation...\n', K);
imp_folds = zeros(K, n_feats);
for k = 1:K
    tr = (cv_idx ~= k); te = ~tr;
    NormS = compute_norm_stats(X(tr,:), crateLabels(tr), chg_idx, dch_idx, crate_groups);
    Xtr = apply_norm_split(X(tr,:), crateLabels(tr), NormS, chg_idx, dch_idx);
    Xte = apply_norm_split(X(te,:), crateLabels(te), NormS, chg_idx, dch_idx);
    rng(42);
    rf_k = fitrensemble(Xtr, Y(tr), 'Method', 'Bag', ...
        'NumLearningCycles', best_mt, ...
        'Learners', templateTree('MinLeafSize', best_ml, 'Surrogate', 'on'));
    y_pred_cv(te) = predict(rf_k, Xte);
    imp_folds(k,:) = oobPermutedPredictorImportance(rf_k);
end

yt = Y; yp = y_pred_cv;
rmse_cv = sqrt(mean((yt-yp).^2)); mae_cv = mean(abs(yt-yp));
r2_cv = 1 - sum((yt-yp).^2)/sum((yt-mean(yt)).^2);
fprintf(' → CV: R²=%.4f | RMSE=%.4f | MAE=%.4f\n', r2_cv, rmse_cv, mae_cv);

% STEP C: 최종 모델
fprintf(' [C] Training Final Model (100%% data, surrogate splits)...\n');
NormS_all = compute_norm_stats(X, crateLabels, chg_idx, dch_idx, crate_groups);
X_all_norm = apply_norm_split(X, crateLabels, NormS_all, chg_idx, dch_idx);
rng(42);
final_rf = fitrensemble(X_all_norm, Y, 'Method', 'Bag', ...
    'NumLearningCycles', best_mt, ...
    'Learners', templateTree('MinLeafSize', best_ml, 'Surrogate', 'on'));

Results_RF.SOH.Y_true     = yt;
Results_RF.SOH.Y_pred     = yp;
Results_RF.SOH.R2         = r2_cv;
Results_RF.SOH.RMSE       = rmse_cv;
Results_RF.SOH.MAE        = mae_cv;
Results_RF.SOH.Importance = mean(imp_folds, 1);
Results_RF.SOH.hp_history = hp_history;
Results_RF.SOH.Model      = final_rf;
Results_RF.SOH.BestParams = struct('NumTrees', best_mt, 'MinLeaf', best_ml);

%% 4. NormStats 저장
NormStats = NormS_all;
NormStats.crate_groups   = crate_groups;
NormStats.crate_vals_num = crate_vals_num;
NormStats.feature_names  = dq_feat_names;
NormStats.chg_idx = chg_idx;
NormStats.dch_idx = dch_idx;
Results_RF.NormStats = NormStats;

%% 5. 저장
save_path = fullfile(save_dir, 'Result_RF_v3.mat');
save(save_path, 'Results_RF', 'X', 'Y', 'cellIDs', 'crateLabels', 'V_bounds');
fprintf('\n >> Saved: %s\n', save_path);

Vis_RF_v3

%% 로컬 함수
function NormS = compute_norm_stats(X_ref, crateLabels_ref, chg_idx, dch_idx, crate_groups)
    NormS.crate_groups = crate_groups;
    for g = 1:length(crate_groups)
        cr = crate_groups{g}; mask = strcmp(crateLabels_ref, cr);
        if sum(mask) > 1
            mu_c = mean(X_ref(mask,chg_idx),1,'omitnan'); sig_c = std(X_ref(mask,chg_idx),0,1,'omitnan');
            mu_d = mean(X_ref(mask,dch_idx),1,'omitnan'); sig_d = std(X_ref(mask,dch_idx),0,1,'omitnan');
        else
            mu_c=zeros(1,length(chg_idx)); sig_c=ones(1,length(chg_idx));
            mu_d=zeros(1,length(dch_idx)); sig_d=ones(1,length(dch_idx));
        end
        sig_c(sig_c<eps)=1; sig_d(sig_d<eps)=1;
        NormS.(cr).mu_chg=mu_c; NormS.(cr).sigma_chg=sig_c;
        NormS.(cr).mu_dch=mu_d; NormS.(cr).sigma_dch=sig_d;
        NormS.(cr).n=sum(mask);
    end
end

function X_out = apply_norm_split(X_target, crateLabels_target, NormS, chg_idx, dch_idx)
    X_out = X_target;
    for g = 1:length(NormS.crate_groups)
        cr = NormS.crate_groups{g}; mask = strcmp(crateLabels_target, cr);
        if ~any(mask), continue; end
        % NaN은 건드리지 않음 (NaN - mu = NaN → surrogate split에서 처리)
        X_out(mask,chg_idx) = (X_target(mask,chg_idx) - NormS.(cr).mu_chg) ./ (NormS.(cr).sigma_chg + eps);
        X_out(mask,dch_idx) = (X_target(mask,dch_idx) - NormS.(cr).mu_dch) ./ (NormS.(cr).sigma_dch + eps);
    end
end
