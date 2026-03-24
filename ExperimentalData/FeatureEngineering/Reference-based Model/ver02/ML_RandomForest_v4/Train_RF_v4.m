% Train_RF_v4.m  — Split Charge/Discharge RF Models
% =====================================================================
% [설계] 충전/방전 분리 모델로 SOH/LLI/LAM 독립 추정
%   - 충전 모델 (RF_chg): 9 피처 → SOH, LLI, LAM
%   - 방전 모델 (RF_dch): 10 피처 → SOH, LLI, LAM
%   - 총 6개 RF 모델 생성
%
% [학습 방법]
%   - Group K-Fold CV (셀 단위 분할)
%   - Fold-local normalization (data leakage 제거)
%   - C-rate 그룹별 정규화 (C_eff 자체는 raw 유지)
%   - HP 탐색: K-fold CV RMSE 기반
%
% [피처 인덱스] (19개 기준, T_avg 제외 — RPT_ML_DataLoader 출력):
%   충전 모델: [1:5, 11, 13, 15, 18]
%     = dQ_chg_S1~5, PkH_chg, PkA_chg, PkPos_chg, C_eff_chg
%   방전 모델: [6:10, 12, 14, 16, 17, 19]
%     = dQ_dch_S1~5, PkH_dch, PkA_dch, PkPos_dch, Energy_dch, C_eff_dch
% =====================================================================

clear; clc; close all;

base_dir    = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
                       'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
loader_path = fullfile(base_dir, 'ML_Models');
addpath(loader_path);

save_dir = fullfile(base_dir, 'ML_RandomForest_v4');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% ============================================================
%  1. 데이터 로드 (ECM 보정된 v4 피처 사용)
%% ============================================================
fprintf('============================================================\n');
fprintf(' [1] Loading Data (ECM-corrected v4 features)\n');
fprintf('============================================================\n');

feat_path  = fullfile(save_dir, 'RPT_FeatureLabelExtractor_v4', 'Feature_Matrix_v4.mat');
label_path = fullfile(base_dir, 'RPT_FeatureLabelExtractor', 'Label_Matrix_ver02.mat');

feat_data  = load(feat_path);   FeatureTable = feat_data.FeatureTable_v4;
label_data = load(label_path);  LabelTable   = label_data.LabelTable_ver02;
feature_names = feat_data.feature_names;

% T_avg 제외
exclude_feats = {'T_chg_avg', 'T_dch_avg'};
keep_mask = ~cellfun(@(f) ismember(f, exclude_feats), feature_names);
feature_names = feature_names(keep_mask);
label_names   = {'SOH', 'LLI', 'LAM'};

% CellID + Cycle + CrateLabel로 행 매칭
[~, ia, ib] = intersect( ...
    strcat(FeatureTable.CellID, num2str(FeatureTable.Cycle), FeatureTable.CrateLabel), ...
    strcat(LabelTable.CellID,   num2str(LabelTable.Cycle),   LabelTable.CrateLabel));
X_all = table2array(FeatureTable(ia, feature_names));
Y_all = table2array(LabelTable(ib, label_names));
cellIDs_all  = FeatureTable.CellID(ia);
crateLabels  = FeatureTable.CrateLabel(ia);

% NaN 제거
nan_rows = any(isnan([X_all, Y_all]), 2);
X = X_all(~nan_rows, :); Y = Y_all(~nan_rows, :);
cellIDs    = cellIDs_all(~nan_rows);
crateLabels = crateLabels(~nan_rows);
fprintf('  -> %d valid samples, %d features\n', size(X,1), size(X,2));

% Group K-Fold (셀 단위)
cells_unique = unique(cellIDs);
K = min(5, length(cells_unique));
rng(42); shuf = cells_unique(randperm(length(cells_unique)));
fold_asgn = mod(0:length(cells_unique)-1, K) + 1;
cv_idx = zeros(size(X,1), 1);
for i = 1:length(cells_unique)
    cv_idx(strcmp(cellIDs, shuf{i})) = fold_asgn(i);
end

n_labels = size(Y, 2);   % 3 (SOH, LLI, LAM)
N        = size(X, 1);
K        = max(cv_idx);

% C-rate 그룹 정의
crate_groups   = {'c01', 'c05', 'c1', 'c2', 'c3'};
crate_vals_num = [0.1,   0.5,   1.0,  2.0,  3.0];

%% ============================================================
%  2. 피처 인덱스 분리 정의
%% ============================================================
% 충전 모델 피처 (19개 전체 피처 중에서 선택)
chg_feat_idx  = [1:5, 11, 13, 15, 18];   % 9 features
chg_norm_idx  = 1:8;                      % 정규화 대상 (9개 중 1~8, C_eff 제외)
chg_ceff_pos  = 9;                        % C_eff_chg 위치 (subspace 내)

% 방전 모델 피처
dch_feat_idx  = [6:10, 12, 14, 16, 17, 19];  % 10 features
dch_norm_idx  = 1:9;                           % 정규화 대상 (10개 중 1~9, C_eff 제외)
dch_ceff_pos  = 10;                            % C_eff_dch 위치 (subspace 내)

chg_feat_names = feature_names(chg_feat_idx);
dch_feat_names = feature_names(dch_feat_idx);

fprintf(' [Charge Model] %d features: %s\n', length(chg_feat_idx), strjoin(chg_feat_names, ', '));
fprintf(' [Discharge Model] %d features: %s\n', length(dch_feat_idx), strjoin(dch_feat_names, ', '));

% 서브스페이스 데이터 추출
X_chg = X(:, chg_feat_idx);  % N × 9
X_dch = X(:, dch_feat_idx);  % N × 10

%% ============================================================
%  3. HP Grid 설정
%% ============================================================
num_trees_list = [50, 100, 200, 500];
min_leaf_list  = [1, 3, 5];

%% ============================================================
%  4. 모델 학습 루프: 충전/방전 × SOH/LLI/LAM
%% ============================================================
model_types = {'chg', 'dch'};
model_data  = {X_chg, X_dch};
model_norm  = {chg_norm_idx, dch_norm_idx};
model_ceff  = {chg_ceff_pos, dch_ceff_pos};
model_fnames = {chg_feat_names, dch_feat_names};

Results_v4 = struct();
Results_v4.label_names = label_names;

for m = 1:2
    mtype      = model_types{m};
    X_sub      = model_data{m};
    norm_idx   = model_norm{m};
    ceff_pos   = model_ceff{m};
    feat_names_sub = model_fnames{m};
    n_feats_sub = size(X_sub, 2);

    Results_v4.(mtype).feature_names = feat_names_sub(:);

    fprintf('\n############################################################\n');
    fprintf(' MODEL TYPE: %s (%d features)\n', upper(mtype), n_feats_sub);
    fprintf('############################################################\n');

    y_pred_cv = zeros(N, n_labels);

    for lbl = 1:n_labels
        fprintf('\n============================================================\n');
        fprintf(' [%s → %s]\n', upper(mtype), label_names{lbl});
        fprintf('============================================================\n');

        % --------------------------------------------------------
        % STEP A: HP 탐색 - K-fold CV RMSE 기반
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
                    % fold-local 정규화
                    NS_k = compute_norm_stats_v4(X_sub(tr,:), crateLabels(tr), norm_idx, crate_groups);
                    Xtr  = apply_norm_v4(X_sub(tr,:), crateLabels(tr), NS_k, norm_idx);
                    Xte  = apply_norm_v4(X_sub(te,:), crateLabels(te), NS_k, norm_idx);
                    rng(42);
                    rf_hp = fitrensemble(Xtr, Y(tr,lbl), 'Method', 'Bag', ...
                        'NumLearningCycles', t, 'Learners', templateTree('MinLeafSize', l));
                    rmse_k(k) = sqrt(mean((Y(te,lbl) - predict(rf_hp, Xte)).^2));
                end
                avg_r = mean(rmse_k);
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
        % --------------------------------------------------------
        fprintf(' [B] Group %d-Fold CV Evaluation...\n', K);
        imp_folds = zeros(K, n_feats_sub);

        for k = 1:K
            tr = (cv_idx ~= k);
            te = ~tr;
            NS_k = compute_norm_stats_v4(X_sub(tr,:), crateLabels(tr), norm_idx, crate_groups);
            Xtr  = apply_norm_v4(X_sub(tr,:), crateLabels(tr), NS_k, norm_idx);
            Xte  = apply_norm_v4(X_sub(te,:), crateLabels(te), NS_k, norm_idx);
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
        NS_all    = compute_norm_stats_v4(X_sub, crateLabels, norm_idx, crate_groups);
        X_all_nm  = apply_norm_v4(X_sub, crateLabels, NS_all, norm_idx);
        rng(42);
        final_rf = fitrensemble(X_all_nm, Y(:, lbl), 'Method', 'Bag', ...
            'NumLearningCycles', best_mt, 'Learners', templateTree('MinLeafSize', best_ml));

        Results_v4.(mtype).(label_names{lbl}).Y_true     = yt;
        Results_v4.(mtype).(label_names{lbl}).Y_pred     = yp;
        Results_v4.(mtype).(label_names{lbl}).R2         = r2_cv;
        Results_v4.(mtype).(label_names{lbl}).RMSE       = rmse_cv;
        Results_v4.(mtype).(label_names{lbl}).MAE        = mae_cv;
        Results_v4.(mtype).(label_names{lbl}).Importance = mean(imp_folds, 1);
        Results_v4.(mtype).(label_names{lbl}).hp_history = hp_history;
        Results_v4.(mtype).(label_names{lbl}).Model      = final_rf;
        Results_v4.(mtype).(label_names{lbl}).BestParams = struct('NumTrees', best_mt, 'MinLeaf', best_ml);
    end

    % NormStats 저장 (전체 랩 데이터 기반)
    NormStats = compute_norm_stats_v4(X_sub, crateLabels, norm_idx, crate_groups);
    NormStats.crate_groups   = crate_groups;
    NormStats.crate_vals_num = crate_vals_num;
    NormStats.feature_names  = feat_names_sub;
    NormStats.norm_idx       = norm_idx;
    NormStats.ceff_pos       = ceff_pos;
    Results_v4.(mtype).NormStats = NormStats;
end

%% ============================================================
%  5. 저장
%% ============================================================
fprintf('\n============================================================\n');
fprintf(' [5] Saving Results\n');
fprintf('============================================================\n');
save_path = fullfile(save_dir, 'Result_RF_v4.mat');
save(save_path, 'Results_v4', 'X', 'Y', 'cellIDs', 'crateLabels');
fprintf(' >> Saved: %s\n', save_path);

%% ============================================================
%  6. Summary 출력
%% ============================================================
fprintf('\n============================================================\n');
fprintf(' SUMMARY: Split Model CV Performance\n');
fprintf('============================================================\n');
fprintf(' %-6s | %-5s | %-8s | %-8s | %-8s\n', 'Model', 'Label', 'R²', 'RMSE', 'MAE');
fprintf(' %s\n', repmat('-', 1, 48));
for m = 1:2
    mtype = model_types{m};
    for lbl = 1:n_labels
        lname = label_names{lbl};
        R = Results_v4.(mtype).(lname);
        fprintf(' %-6s | %-5s | %8.4f | %8.4f | %8.4f\n', ...
            upper(mtype), lname, R.R2, R.RMSE, R.MAE);
    end
end

%% ============================================================
%  7. 시각화
%% ============================================================
Vis_RF_v4


%% ============================================================
%  Local Functions
%% ============================================================

function NS = compute_norm_stats_v4(X_ref, crateLabels_ref, norm_idx, crate_groups)
% 정규화 통계 계산 (C-rate 그룹별, norm_idx 피처만 대상)
    NS.crate_groups = crate_groups;
    for g = 1:length(crate_groups)
        cr   = crate_groups{g};
        mask = strcmp(crateLabels_ref, cr);
        if sum(mask) > 1
            mu  = mean(X_ref(mask, norm_idx), 1, 'omitnan');
            sig = std( X_ref(mask, norm_idx), 0, 1, 'omitnan');
        else
            mu  = zeros(1, length(norm_idx));
            sig = ones( 1, length(norm_idx));
        end
        sig(sig < eps) = 1;
        NS.(cr).mu    = mu;
        NS.(cr).sigma = sig;
        NS.(cr).n     = sum(mask);
    end
end

function X_out = apply_norm_v4(X_target, crateLabels_target, NS, norm_idx)
% C-rate 그룹별 정규화 적용 (norm_idx 피처만 정규화, 나머지 raw)
    X_out = X_target;
    for g = 1:length(NS.crate_groups)
        cr   = NS.crate_groups{g};
        mask = strcmp(crateLabels_target, cr);
        if ~any(mask), continue; end
        mu  = NS.(cr).mu;
        sig = NS.(cr).sigma;
        X_out(mask, norm_idx) = (X_target(mask, norm_idx) - mu) ./ (sig + eps);
    end
end
