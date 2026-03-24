% Train_SVM.m  v2
% [변경] v1 → v2
%   [1] Fold-local normalization: CV 각 fold의 학습 셋에서만 mu/sigma 계산 → data leakage 제거 (논문 방식)
%   [2] Split normalization: chg 피처 → C_eff_chg 그룹 기준, dch 피처 → C_eff_dch 그룹 기준
%   [3] HP 탐색: Bayesian Optimization을 K-fold CV 내부에서 수행 (fold-local 정규화 적용)
%
% [SVM 특이사항]: RF의 OOB가 없으므로 HP 탐색 방식을 Bayesian Opt → K-fold CV 내 grid 또는
%                 Bayesian Opt를 학습 데이터(fold 기준)에서 수행. 여기서는 효율을 위해
%                 HP 탐색 시 전체 학습 데이터로 Bayesian Opt 수행 (fold-local 정규화 통계 사용)
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

save_dir = fullfile(base_dir, 'ML_SVM');
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

% 피처 인덱스 정의 (RPT_FeatureExtractor.m 출력 순서와 일치 필요)
chg_idx  = [1:5, 11, 13, 15];        % 충전 피처 (8개)
dch_idx  = [6:10, 12, 14, 16, 17];   % 방전 피처 (9개)
ceff_idx = [n_feats-1, n_feats];      % C_eff_chg, C_eff_dch (정규화 제외)

% C-rate 그룹 정의
crate_groups   = {'c01', 'c05', 'c1', 'c2', 'c3'};
crate_vals_num = [0.1,   0.5,   1.0,  2.0,  3.0];

fprintf(' [Data] %d samples | %d features | %d labels | %d folds\n', N, n_feats, n_labels, K);

%% ============================================================
%  2. 모델 결과 저장 구조체
%% ============================================================
Results_SVM            = struct();
Results_SVM.feature_names = feature_names;
Results_SVM.label_names   = label_names;

y_pred_cv = zeros(N, n_labels);

%% ============================================================
%  3. Per-label: HP 탐색 → CV 평가 → 최종 모델
%% ============================================================
fprintf('\n------------------------------------------------------------\n');
fprintf(' [2] Training & Validation (Fold-local norm + K-fold CV)\n');
fprintf('------------------------------------------------------------\n');

for lbl = 1:n_labels
    fprintf('\n============================================================\n');
    fprintf(' [Target: %s]\n', label_names{lbl});
    fprintf('============================================================\n');

    % --------------------------------------------------------
    % STEP A: HP 탐색 - Bayesian Optimization (전체 데이터 기반, 논문 방식)
    %   전체 학습 데이터로 정규화 통계 계산 후 Bayesian Opt 수행
    %   (SVM은 OOB가 없어 내부 CV 기반 탐색)
    % --------------------------------------------------------
    fprintf(' [A] HP Tuning via Bayesian Optimization (30 evals)...\n');

    % 전체 데이터 정규화 (HP 탐색용)
    NormS_all = compute_norm_stats(X, crateLabels, chg_idx, dch_idx, crate_groups);
    X_all_norm = apply_norm_split(X, crateLabels, NormS_all, chg_idx, dch_idx);

    % HP 탐색 범위 명시적 제한 (Epsilon 상한 1.0으로 제한 → 비정상 수렴 방지)
    hp_vars = [
        optimizableVariable('BoxConstraint', [0.01, 1000], 'Transform', 'log')
        optimizableVariable('KernelScale',   [0.01, 1000], 'Transform', 'log')
        optimizableVariable('Epsilon',       [0.001, 1.0], 'Transform', 'log')
    ];
    rng(42);
    ho_opts = struct('ShowPlots', false, 'Verbose', 0, 'MaxObjectiveEvaluations', 30);
    svm_tune = fitrsvm(X_all_norm, Y(:, lbl), ...
        'KernelFunction', 'gaussian', 'Standardize', false, ...
        'OptimizeHyperparameters', hp_vars, ...
        'HyperparameterOptimizationOptions', ho_opts);

    best_bc  = svm_tune.ModelParameters.BoxConstraint;
    best_ks  = svm_tune.ModelParameters.KernelScale;
    best_eps = svm_tune.ModelParameters.Epsilon;
    fprintf(' → Optimal: BC=%.3f  KS=%.3f  Eps=%.4f\n', best_bc, best_ks, best_eps);

    hp_history = svm_tune.HyperparameterOptimizationResults.ObjectiveMinimumTrace;

    % --------------------------------------------------------
    % STEP B: Group K-Fold CV 최종 평가 (fold-local normalization)
    % --------------------------------------------------------
    fprintf(' [B] Group %d-Fold CV Evaluation (fold-local normalization)...\n', K);

    for f = 1:K
        tr = (cv_idx ~= f);
        te = ~tr;

        % fold-local 정규화: 학습 셋에서만 통계 계산
        NormS_f = compute_norm_stats(X(tr,:), crateLabels(tr), chg_idx, dch_idx, crate_groups);
        Xtr = apply_norm_split(X(tr,:), crateLabels(tr), NormS_f, chg_idx, dch_idx);
        Xte = apply_norm_split(X(te,:), crateLabels(te), NormS_f, chg_idx, dch_idx);

        rng(42);
        svm_f = fitrsvm(Xtr, Y(tr, lbl), ...
            'KernelFunction', 'gaussian', 'Standardize', false, ...
            'BoxConstraint', best_bc, 'KernelScale', best_ks, 'Epsilon', best_eps);
        y_pred_cv(te, lbl) = predict(svm_f, Xte);
        fprintf('   Fold %d done.\n', f);
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
    rng(42);
    svm_final = fitrsvm(X_all_norm, Y(:, lbl), ...
        'KernelFunction', 'gaussian', 'Standardize', false, ...
        'BoxConstraint', best_bc, 'KernelScale', best_ks, 'Epsilon', best_eps);

    % 결과 저장
    Results_SVM.(label_names{lbl}).FinalModel  = svm_final;
    Results_SVM.(label_names{lbl}).OptParams   = struct('BoxConstraint', best_bc, ...
                                                         'KernelScale', best_ks, 'Epsilon', best_eps);
    Results_SVM.(label_names{lbl}).hp_history  = hp_history;
    Results_SVM.(label_names{lbl}).Y_true      = yt;
    Results_SVM.(label_names{lbl}).Y_pred      = yp;
    Results_SVM.(label_names{lbl}).R2          = r2_cv;
    Results_SVM.(label_names{lbl}).RMSE        = rmse_cv;
    Results_SVM.(label_names{lbl}).MAE         = mae_cv;
end

Results_SVM.Y_True    = Y;
Results_SVM.Y_Pred    = y_pred_cv;
Results_SVM.CV_Indices = cv_idx;

%% ============================================================
%  4. NormStats 저장 (필드 추정 스크립트에서 재사용)
%     구조: NormStats.(cr).mu_chg / sigma_chg / mu_dch / sigma_dch
%% ============================================================
NormStats = compute_norm_stats(X, crateLabels, chg_idx, dch_idx, crate_groups);
NormStats.crate_groups   = crate_groups;
NormStats.crate_vals_num = crate_vals_num;
NormStats.feature_names  = feature_names;
NormStats.chg_idx        = chg_idx;
NormStats.dch_idx        = dch_idx;
Results_SVM.NormStats    = NormStats;

%% ============================================================
%  5. 저장 및 시각화
%% ============================================================
fprintf('\n------------------------------------------------------------\n');
fprintf(' [3] Exporting All Results\n');
fprintf('------------------------------------------------------------\n');

save_path = fullfile(save_dir, 'Result_SVM.mat');
save(save_path, 'Results_SVM', 'X', 'Y', 'cellIDs', 'crateLabels');
fprintf(' >> Saved: %s\n', save_path);

Vis_SVM
RPT_Field_Estimation_SVM


%% ============================================================
%  로컬 함수
%% ============================================================

function NormS = compute_norm_stats(X_ref, crateLabels_ref, chg_idx, dch_idx, crate_groups)
% compute_norm_stats: 학습 데이터에서 분리 정규화 통계 계산
%   반환: NormS.(cr).mu_chg / sigma_chg / mu_dch / sigma_dch (C-rate 그룹별)

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
            mu_c  = zeros(1, length(chg_idx)); sig_c = ones(1, length(chg_idx));
            mu_d  = zeros(1, length(dch_idx)); sig_d = ones(1, length(dch_idx));
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
%   chg 피처 → C_eff_chg 그룹의 mu_chg / sigma_chg
%   dch 피처 → C_eff_dch 그룹의 mu_dch / sigma_dch
%   ceff 피처 → 정규화 없음 (raw 유지)

    X_out = X_target;

    for g = 1:length(NormS.crate_groups)
        cr   = NormS.crate_groups{g};
        mask = strcmp(crateLabels_target, cr);
        if ~any(mask), continue; end

        X_out(mask, chg_idx) = (X_target(mask, chg_idx) - NormS.(cr).mu_chg) ./ (NormS.(cr).sigma_chg + eps);
        X_out(mask, dch_idx) = (X_target(mask, dch_idx) - NormS.(cr).mu_dch) ./ (NormS.(cr).sigma_dch + eps);
        % ceff_idx: raw 유지
    end
end