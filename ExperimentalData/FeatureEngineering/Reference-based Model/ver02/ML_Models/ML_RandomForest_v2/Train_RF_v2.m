% Train_RF_v2.m  —  SOH only + dQ segments only (10 features)
% [v2 변경사항]
%   - 라벨: SOH만 사용 (LLI, LAM 제거)
%   - 피처: dQ_chg_S1~5 + dQ_dch_S1~5 = 10개 (Peak, Energy, C_eff 제거)
%   - C_eff_chg / C_eff_dch: 정규화 그룹 결정에만 사용 (피처로 입력하지 않음)
%   - fold-local normalization, split chg/dch, K-fold CV HP 탐색 동일 유지
clear; clc; close all;

base_dir    = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
                       'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
loader_path = fullfile(base_dir, 'ML_Models');
addpath(loader_path);

save_dir = fullfile(base_dir, 'ML_RandomForest_v2');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% 1. 데이터 로드
fprintf('------------------------------------------------------------\n');
fprintf(' [1] Loading Data via RPT_ML_DataLoader\n');
fprintf('------------------------------------------------------------\n');
[X_full, Y_full, cv_idx, feature_names_full, label_names_full, cellIDs, crateLabels] = RPT_ML_DataLoader(base_dir);

% --- 피처 선택: dQ 세그먼트만 (10개) ---
dq_feat_names = {'dQ_chg_S1','dQ_chg_S2','dQ_chg_S3','dQ_chg_S4','dQ_chg_S5', ...
                  'dQ_dch_S1','dQ_dch_S2','dQ_dch_S3','dQ_dch_S4','dQ_dch_S5'};
feat_idx = find(ismember(feature_names_full, dq_feat_names));
X = X_full(:, feat_idx);
feature_names = feature_names_full(feat_idx);

% --- 라벨 선택: SOH만 ---
soh_idx = find(strcmp(label_names_full, 'SOH'));
Y = Y_full(:, soh_idx);
label_names = {'SOH'};

% --- C_eff 컬럼 (정규화용, 피처 아님) ---
ceff_chg_col = find(strcmp(feature_names_full, 'C_eff_chg'));
ceff_dch_col = find(strcmp(feature_names_full, 'C_eff_dch'));

n_feats  = size(X, 2);
N        = size(X, 1);
K        = max(cv_idx);

% 피처 인덱스 정의 (10개 피처 기준)
chg_idx  = 1:5;        % dQ_chg_S1~5
dch_idx  = 6:10;       % dQ_dch_S1~5

% C-rate 그룹 정의
crate_groups   = {'c01', 'c05', 'c1', 'c2', 'c3'};
crate_vals_num = [0.1,   0.5,   1.0,  2.0,  3.0];

fprintf(' [Data] %d samples | %d features (dQ only) | 1 label (SOH) | %d folds\n', N, n_feats, K);

%% 2. HP Grid 설정
num_trees_list = [50, 100, 200, 500];
min_leaf_list  = [1, 3, 5];

%% 3. HP 탐색 → CV 평가 → 최종 모델 (SOH only)
Results_RF = struct();
Results_RF.feature_names = feature_names(:);
Results_RF.label_names   = label_names;

y_pred_cv = zeros(N, 1);

fprintf('\n============================================================\n');
fprintf(' [Target: SOH] (dQ features only)\n');
fprintf('============================================================\n');

% STEP A: HP 탐색 (K-fold CV RMSE)
fprintf(' [A] HP Tuning via %d-Fold CV...\n', K);
best_rmse = inf; best_mt = 100; best_ml = 1;
hp_history = [];

for t = num_trees_list
    for l = min_leaf_list
        rmse_k = zeros(K, 1);
        for k = 1:K
            tr = (cv_idx ~= k); te = ~tr;
            NormS_k = compute_norm_stats(X(tr,:), crateLabels(tr), chg_idx, dch_idx, crate_groups);
            Xtr = apply_norm_split(X(tr,:), crateLabels(tr), NormS_k, chg_idx, dch_idx);
            Xte = apply_norm_split(X(te,:), crateLabels(te), NormS_k, chg_idx, dch_idx);
            rng(42);
            rf_hp = fitrensemble(Xtr, Y(tr), 'Method', 'Bag', ...
                'NumLearningCycles', t, 'Learners', templateTree('MinLeafSize', l));
            rmse_k(k) = sqrt(mean((Y(te) - predict(rf_hp, Xte)).^2));
        end
        avg_r = mean(rmse_k);
        hp_history = [hp_history; t, l, avg_r]; %#ok<AGROW>
        if avg_r < best_rmse, best_rmse = avg_r; best_mt = t; best_ml = l; end
        fprintf('   n_trees=%3d  min_leaf=%d  CV_RMSE=%.4f\n', t, l, avg_r);
    end
end
fprintf(' → Best: n_trees=%d  min_leaf=%d  CV_RMSE=%.4f\n', best_mt, best_ml, best_rmse);

% STEP B: CV 최종 평가 (fold-local normalization)
fprintf(' [B] Group %d-Fold CV Evaluation...\n', K);
imp_folds = zeros(K, n_feats);
for k = 1:K
    tr = (cv_idx ~= k); te = ~tr;
    NormS_k = compute_norm_stats(X(tr,:), crateLabels(tr), chg_idx, dch_idx, crate_groups);
    Xtr = apply_norm_split(X(tr,:), crateLabels(tr), NormS_k, chg_idx, dch_idx);
    Xte = apply_norm_split(X(te,:), crateLabels(te), NormS_k, chg_idx, dch_idx);
    rng(42);
    rf_k = fitrensemble(Xtr, Y(tr), 'Method', 'Bag', ...
        'NumLearningCycles', best_mt, 'Learners', templateTree('MinLeafSize', best_ml));
    y_pred_cv(te) = predict(rf_k, Xte);
    imp_folds(k,:) = oobPermutedPredictorImportance(rf_k);
end

yt = Y; yp = y_pred_cv;
rmse_cv = sqrt(mean((yt-yp).^2)); mae_cv = mean(abs(yt-yp));
r2_cv = 1 - sum((yt-yp).^2)/sum((yt-mean(yt)).^2);
fprintf(' → CV: R²=%.4f | RMSE=%.4f | MAE=%.4f\n', r2_cv, rmse_cv, mae_cv);

% STEP C: 최종 모델 (100% 학습)
fprintf(' [C] Training Final Model (100%% lab data)...\n');
NormS_all  = compute_norm_stats(X, crateLabels, chg_idx, dch_idx, crate_groups);
X_all_norm = apply_norm_split(X, crateLabels, NormS_all, chg_idx, dch_idx);
rng(42);
final_rf = fitrensemble(X_all_norm, Y, 'Method', 'Bag', ...
    'NumLearningCycles', best_mt, 'Learners', templateTree('MinLeafSize', best_ml));

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
NormStats = compute_norm_stats(X, crateLabels, chg_idx, dch_idx, crate_groups);
NormStats.crate_groups   = crate_groups;
NormStats.crate_vals_num = crate_vals_num;
NormStats.feature_names  = feature_names;
NormStats.chg_idx        = chg_idx;
NormStats.dch_idx        = dch_idx;
Results_RF.NormStats     = NormStats;

%% 5. 저장 및 시각화
save_path = fullfile(save_dir, 'Result_RF_v2.mat');
save(save_path, 'Results_RF', 'X', 'Y', 'cellIDs', 'crateLabels');
fprintf('\n >> Saved: %s\n', save_path);

Vis_RF_v2

%% 로컬 함수
function NormS = compute_norm_stats(X_ref, crateLabels_ref, chg_idx, dch_idx, crate_groups)
    NormS.crate_groups = crate_groups;
    for g = 1:length(crate_groups)
        cr = crate_groups{g}; mask = strcmp(crateLabels_ref, cr);
        if sum(mask) > 1
            mu_c = mean(X_ref(mask,chg_idx),1,'omitnan'); sig_c = std(X_ref(mask,chg_idx),0,1,'omitnan');
            mu_d = mean(X_ref(mask,dch_idx),1,'omitnan'); sig_d = std(X_ref(mask,dch_idx),0,1,'omitnan');
        else
            mu_c = zeros(1,length(chg_idx)); sig_c = ones(1,length(chg_idx));
            mu_d = zeros(1,length(dch_idx)); sig_d = ones(1,length(dch_idx));
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
        X_out(mask,chg_idx) = (X_target(mask,chg_idx) - NormS.(cr).mu_chg) ./ (NormS.(cr).sigma_chg + eps);
        X_out(mask,dch_idx) = (X_target(mask,dch_idx) - NormS.(cr).mu_dch) ./ (NormS.(cr).sigma_dch + eps);
    end
end
