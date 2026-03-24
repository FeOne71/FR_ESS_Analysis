% Train_RF.m  v3
% [변경] v2 → v3
%   [1] Fold-local normalization: CV 각 fold의 학습 셋에서만 mu/sigma 계산 → data leakage 제거 (논문 방식)
%   [2] Split normalization: chg 피처 → C_eff_chg 그룹 기준, dch 피처 → C_eff_dch 그룹 기준 (분리 정규화)
%   [3] HP 탐색: oobPredict 대신 K-fold CV RMSE 기반으로 변경 (논문 방식)
%
% [최종 모델]: 전체 랩 데이터 100% 사용 (final_rf)
% [NormStats]:  전체 랩 데이터 기반 통계 → 필드 추정 스크립트에서 재사용
%
% 피처 인덱스 (19개 기준, T_avg 제외):
%   chg_idx  = [1:5, 11, 13, 15]     : dQ_chg×5, PkH_chg, PkA_chg, PkPos_chg  (8개)
%   dch_idx  = [6:10, 12, 14, 16, 17]: dQ_dch×5, PkH_dch, PkA_dch, PkPos_dch, Energy_dch (9개)
%   ceff_idx = [18, 19]               : C_eff_chg, C_eff_dch (정규화 없음)

clear; clc; close all;

base_dir    = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
                       'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
loader_path = fullfile(base_dir, 'ML_Models');
addpath(loader_path);

save_dir = fullfile(base_dir, 'ML_RandomForest');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% ============================================================
%  1. 데이터 로드
%% ============================================================
fprintf('------------------------------------------------------------\n');
fprintf(' [1] Loading Data via RPT_ML_DataLoader\n');
fprintf('------------------------------------------------------------\n');
[X, Y, cv_idx, feature_names, label_names, cellIDs, crateLabels] = RPT_ML_DataLoader(base_dir);

n_labels = size(Y, 2);
n_feats  = size(X, 2);
N        = size(X, 1);
K        = max(cv_idx);

% 피처 인덱스 정의 (RPT_FeatureExtractor.m 출력 순서와 일치해야 함)
chg_idx  = [1:5, 11, 13, 15];        % 충전 피처 (8개)
dch_idx  = [6:10, 12, 14, 16, 17];   % 방전 피처 (9개)
ceff_idx = [n_feats-1, n_feats];      % C_eff_chg, C_eff_dch (정규화 제외)

% C-rate 그룹 정의
crate_groups   = {'c01', 'c05', 'c1', 'c2', 'c3'};
crate_vals_num = [0.1,   0.5,   1.0,  2.0,  3.0];

fprintf(' [Data] %d samples | %d features | %d labels | %d folds\n', N, n_feats, n_labels, K);
fprintf(' [chg_idx] %s\n', mat2str(chg_idx));
fprintf(' [dch_idx] %s\n', mat2str(dch_idx));

%% ============================================================
%  2. HP Grid 설정
%% ============================================================
num_trees_list = [50, 100, 200, 500];
min_leaf_list  = [1, 3, 5];

%% ============================================================
%  3. 모델별 학습: HP 탐색 → CV 평가 → 최종 모델
%% ============================================================
Results_RF            = struct();
Results_RF.label_names = label_names;
Results_RF.Merged.feature_names = feature_names(:);

y_pred_cv  = zeros(N, n_labels);

for lbl = 1:n_labels
    fprintf('\n============================================================\n');
    fprintf(' [Target: %s]\n', label_names{lbl});
    fprintf('============================================================\n');

    % --------------------------------------------------------
    % STEP A: HP 탐색 - K-fold CV RMSE 기반 (논문 방식)
    % --------------------------------------------------------
    fprintf(' [A] HP Tuning via %d-Fold CV (%d combinations)...\n', K, ...
        length(num_trees_list) * length(min_leaf_list));

    best_rmse_hp = inf;
    best_mt      = 100;
    best_ml      = 1;
    hp_history   = [];

    for t = num_trees_list
        for l = min_leaf_list
            rmse_k = zeros(K, 1);
            for k = 1:K
                tr = (cv_idx ~= k);
                te = ~tr;
                % fold-local 정규화: 학습 셋에서만 통계 계산
                NormS_k = compute_norm_stats(X(tr,:), crateLabels(tr), chg_idx, dch_idx, crate_groups);
                Xtr = apply_norm_split(X(tr,:), crateLabels(tr), NormS_k, chg_idx, dch_idx);
                Xte = apply_norm_split(X(te,:), crateLabels(te), NormS_k, chg_idx, dch_idx);
                rng(42);
                rf_hp = fitrensemble(Xtr, Y(tr,lbl), 'Method', 'Bag', ...
                    'NumLearningCycles', t, 'Learners', templateTree('MinLeafSize', l));
                rmse_k(k) = sqrt(mean((Y(te,lbl) - predict(rf_hp, Xte)).^2));
            end
            avg_r    = mean(rmse_k);
            hp_history = [hp_history; t, l, avg_r]; %#ok<AGROW>
            if avg_r < best_rmse_hp
                best_rmse_hp = avg_r; best_mt = t; best_ml = l;
            end
            fprintf('   n_trees=%3d  min_leaf=%d  CV_RMSE=%.4f\n', t, l, avg_r);
        end
    end
    fprintf(' → Best: n_trees=%d  min_leaf=%d  CV_RMSE=%.4f\n', best_mt, best_ml, best_rmse_hp);

    % --------------------------------------------------------
    % STEP B: 최적 HP로 Group K-Fold 최종 평가
    %         (fold-local normalization 적용)
    % --------------------------------------------------------
    fprintf(' [B] Group %d-Fold CV Evaluation (fold-local normalization)...\n', K);
    imp_folds = zeros(K, n_feats);

    for k = 1:K
        tr = (cv_idx ~= k);
        te = ~tr;
        NormS_k = compute_norm_stats(X(tr,:), crateLabels(tr), chg_idx, dch_idx, crate_groups);
        Xtr = apply_norm_split(X(tr,:), crateLabels(tr), NormS_k, chg_idx, dch_idx);
        Xte = apply_norm_split(X(te,:), crateLabels(te), NormS_k, chg_idx, dch_idx);
        rng(42);
        rf_k = fitrensemble(Xtr, Y(tr,lbl), 'Method', 'Bag', ...
            'NumLearningCycles', best_mt, 'Learners', templateTree('MinLeafSize', best_ml));
        y_pred_cv(te, lbl) = predict(rf_k, Xte);
        imp_folds(k, :)    = oobPermutedPredictorImportance(rf_k);
    end

    yt      = Y(:, lbl);
    yp      = y_pred_cv(:, lbl);
    rmse_cv = sqrt(mean((yt - yp).^2));
    mae_cv  = mean(abs(yt - yp));
    r2_cv   = 1 - sum((yt - yp).^2) / sum((yt - mean(yt)).^2);
    fprintf(' → CV: R²=%.4f | RMSE=%.4f | MAE=%.4f\n', r2_cv, rmse_cv, mae_cv);

    % --------------------------------------------------------
    % STEP C: 최종 모델 - 전체 랩 데이터 100% 사용
    % --------------------------------------------------------
    fprintf(' [C] Training Final Model (100%% lab data)...\n');
    NormS_all  = compute_norm_stats(X, crateLabels, chg_idx, dch_idx, crate_groups);
    X_all_norm = apply_norm_split(X, crateLabels, NormS_all, chg_idx, dch_idx);
    rng(42);
    final_rf = fitrensemble(X_all_norm, Y(:, lbl), 'Method', 'Bag', ...
        'NumLearningCycles', best_mt, 'Learners', templateTree('MinLeafSize', best_ml));

    Results_RF.Merged.(label_names{lbl}).Y_true     = yt;
    Results_RF.Merged.(label_names{lbl}).Y_pred     = yp;
    Results_RF.Merged.(label_names{lbl}).R2         = r2_cv;
    Results_RF.Merged.(label_names{lbl}).RMSE       = rmse_cv;
    Results_RF.Merged.(label_names{lbl}).MAE        = mae_cv;
    Results_RF.Merged.(label_names{lbl}).Importance = mean(imp_folds, 1);
    Results_RF.Merged.(label_names{lbl}).hp_history = hp_history;
    Results_RF.Merged.(label_names{lbl}).Model      = final_rf;
    Results_RF.Merged.(label_names{lbl}).BestParams = struct('NumTrees', best_mt, 'MinLeaf', best_ml);
end

%% ============================================================
%  4. NormStats 저장 (필드 추정 스크립트에서 재사용)
%     - 전체 랩 데이터 기반 통계
%     - 구조: NormStats.(cr).mu_chg / sigma_chg / mu_dch / sigma_dch
%% ============================================================
NormStats = compute_norm_stats(X, crateLabels, chg_idx, dch_idx, crate_groups);
NormStats.crate_groups   = crate_groups;
NormStats.crate_vals_num = crate_vals_num;
NormStats.feature_names  = feature_names;
NormStats.chg_idx        = chg_idx;
NormStats.dch_idx        = dch_idx;
Results_RF.NormStats     = NormStats;

%% ============================================================
%  5. 저장 및 시각화
%% ============================================================
fprintf('\n------------------------------------------------------------\n');
fprintf(' [4] Exporting All Results\n');
fprintf('------------------------------------------------------------\n');
save_path = fullfile(save_dir, 'Result_RandomForest.mat');
save(save_path, 'Results_RF', 'X', 'Y', 'cellIDs', 'crateLabels');
fprintf(' >> Saved: %s\n', save_path);

Vis_RF


%% ============================================================
%  로컬 함수 (MATLAB R2016b 이후 스크립트 내 local function 지원)
%% ============================================================

function NormS = compute_norm_stats(X_ref, crateLabels_ref, chg_idx, dch_idx, crate_groups)
% compute_norm_stats: 학습 데이터에서 분리 정규화 통계 계산
%   X_ref          : 통계 계산 기준 데이터 (학습 셋 또는 전체)
%   crateLabels_ref: X_ref의 C-rate 그룹 레이블
%   chg_idx        : 충전 피처 인덱스
%   dch_idx        : 방전 피처 인덱스
%   crate_groups   : {'c01','c05','c1','c2','c3'}
%
%   반환 NormS.(cr).mu_chg / sigma_chg / mu_dch / sigma_dch (각 C-rate 그룹별)

    NormS.crate_groups = crate_groups;

    for g = 1:length(crate_groups)
        cr   = crate_groups{g};
        mask = strcmp(crateLabels_ref, cr);

        if sum(mask) > 1
            mu_c  = mean(X_ref(mask, chg_idx), 1, 'omitnan');
            sig_c = std( X_ref(mask, chg_idx), 0, 1, 'omitnan');
            mu_d  = mean(X_ref(mask, dch_idx), 1, 'omitnan');
            sig_d = std( X_ref(mask, dch_idx), 0, 1, 'omitnan');
        else
            % 해당 그룹 샘플 없음 → 정규화 무효 (identity)
            mu_c  = zeros(1, length(chg_idx));
            sig_c = ones( 1, length(chg_idx));
            mu_d  = zeros(1, length(dch_idx));
            sig_d = ones( 1, length(dch_idx));
        end
        sig_c(sig_c < eps) = 1;
        sig_d(sig_d < eps) = 1;

        NormS.(cr).mu_chg    = mu_c;
        NormS.(cr).sigma_chg = sig_c;
        NormS.(cr).mu_dch    = mu_d;
        NormS.(cr).sigma_dch = sig_d;
        NormS.(cr).n         = sum(mask);
    end
end


function X_out = apply_norm_split(X_target, crateLabels_target, NormS, chg_idx, dch_idx)
% apply_norm_split: NormS 통계를 사용해 X_target에 분리 정규화 적용
%   chg 피처  → C_eff_chg 그룹의 mu_chg / sigma_chg
%   dch 피처  → C_eff_dch 그룹의 mu_dch / sigma_dch
%   ceff 피처 → 정규화 없음 (raw 유지)
%
%   [Note] 랩 데이터는 crateLabel이 chg/dch 동일 → 같은 그룹 통계 적용
%           필드 추정에서는 RPT_Field_Estimation_RF.m이 chg/dch 별도 보간 처리

    X_out = X_target;

    for g = 1:length(NormS.crate_groups)
        cr   = NormS.crate_groups{g};
        mask = strcmp(crateLabels_target, cr);
        if ~any(mask), continue; end

        mu_c  = NormS.(cr).mu_chg;
        sig_c = NormS.(cr).sigma_chg;
        mu_d  = NormS.(cr).mu_dch;
        sig_d = NormS.(cr).sigma_dch;

        X_out(mask, chg_idx) = (X_target(mask, chg_idx) - mu_c) ./ (sig_c + eps);
        X_out(mask, dch_idx) = (X_target(mask, dch_idx) - mu_d) ./ (sig_d + eps);
        % ceff_idx: 변경 없음 (raw 유지)
    end
end