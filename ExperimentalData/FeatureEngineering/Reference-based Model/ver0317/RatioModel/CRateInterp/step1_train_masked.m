% step1_train_masked.m  (CRateInterp)
% C-rate별 마스킹 증강 학습 + LOCO-CV + 모델 저장
% 각 C-rate의 43행 → 마스킹 증강 → 1,333행
% 피처: dQ/sum(valid) ratio (부분 세그먼트 대응)
% 검증: LOCO-CV (test = 원본, train = 증강 데이터)

clear; clc; close all;
rng(42);

%% Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
ciDir   = fullfile(verDir, 'RatioModel', 'CRateInterp');
visDir  = fullfile(ciDir, 'Visualization');
if ~exist(visDir,'dir'), mkdir(visDir); end

%% Load Lab Data
d  = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat')); FM = d.FM;
Q_nom = 64;
dQ_c_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i), 3:12, 'UniformOutput', false);
dQ_d_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false);
dQ_c_raw = FM{:, dQ_c_cols};  % 215×10
dQ_d_raw = FM{:, dQ_d_cols};  % 215×11
C_eff_c = FM.C_eff_chg; C_eff_d = FM.C_eff_dch;
y_lab   = FM.Static_Capacity / Q_nom * 100;
cellIDs = FM.CellID;
conditions = FM.Condition;
uConds = unique(conditions);
crate_vals = [0.1, 0.5, 1.0, 2.0, 3.0];
N_mask = 30;  % 마스킹 버전 수

fprintf('Lab: %d rows, %d C-rates, Mask=%d\n', size(dQ_c_raw,1), length(uConds), N_mask);

%% Masking Augmentation 함수
% dQ_raw: (n_segs,) 행벡터, NaN 연속 구간 마스킹 후 ratio 반환
function r = apply_mask_ratio(dq_raw)
    n = length(dq_raw);
    valid = find(~isnan(dq_raw));
    if length(valid) < 2, r = nan(1,n); return; end
    % 랜덤 윈도우: 1~floor(n*0.7) 크기
    win = randi([1, max(1, floor(n*0.7))]);
    st  = randi([1, n-win+1]);
    masked = dq_raw;
    masked(st:st+win-1) = NaN;
    s = sum(masked, 'omitnan');
    if s < 1e-10, r = nan(1,n); return; end
    r = masked / s;
end

%% C-rate별 LOCO-CV + 모델 학습
CRateModels = struct();
fprintf('\n=== C-rate별 마스킹 증강 LOCO-CV ===\n');
fprintf('%-6s %-7s %-7s %-7s\n','Crate','N_aug','RMSE%','R2');

for ci = 1:length(uConds)
    cond = uConds(ci); cr = crate_vals(ci);
    cM   = strcmp(conditions, cond);
    X_c_raw_c = dQ_c_raw(cM,:); X_c_raw_d = dQ_d_raw(cM,:);
    Cc_c = C_eff_c(cM); Cc_d = C_eff_d(cM);
    y_c  = y_lab(cM); ids_c = cellIDs(cM);
    n_c  = size(X_c_raw_c, 1);  % 43
    uCells = unique(ids_c);

    %-- Build augmented dataset for this C-rate --
    X_aug_c = zeros(n_c*(N_mask+1), 10);
    X_aug_d = zeros(n_c*(N_mask+1), 11);
    Cc_aug_c = zeros(n_c*(N_mask+1),1);
    Cc_aug_d = zeros(n_c*(N_mask+1),1);
    y_aug   = zeros(n_c*(N_mask+1),1);
    cell_aug = strings(n_c*(N_mask+1),1);
    row_is_orig = false(n_c*(N_mask+1),1);

    idx = 1;
    for i = 1:n_c
        % 원본 ratio
        rc = X_c_raw_c(i,:) / sum(X_c_raw_c(i,:),'omitnan');
        rd = X_c_raw_d(i,:) / sum(X_c_raw_d(i,:),'omitnan');
        X_aug_c(idx,:) = rc; X_aug_d(idx,:) = rd;
        Cc_aug_c(idx) = Cc_c(i); Cc_aug_d(idx) = Cc_d(i);
        y_aug(idx) = y_c(i); cell_aug(idx) = ids_c(i);
        row_is_orig(idx) = true; idx = idx+1;
        % 마스킹 N_mask개
        for m = 1:N_mask
            X_aug_c(idx,:) = apply_mask_ratio(X_c_raw_c(i,:));
            X_aug_d(idx,:) = apply_mask_ratio(X_c_raw_d(i,:));
            Cc_aug_c(idx) = Cc_c(i); Cc_aug_d(idx) = Cc_d(i);
            y_aug(idx) = y_c(i); cell_aug(idx) = ids_c(i);
            idx = idx+1;
        end
    end
    X_aug = [X_aug_c, X_aug_d, Cc_aug_c, Cc_aug_d];  % n_aug×23

    % Standardize (전체 기준, NaN 제외)
    mu_c  = mean(X_aug,1,'omitnan');
    sig_c = std(X_aug,0,1,'omitnan'); sig_c(sig_c==0)=1;
    X_aug_s = (X_aug - mu_c) ./ sig_c;

    %-- LOCO-CV --
    y_cv_orig = nan(n_c, 1);
    orig_idx  = find(row_is_orig);

    for fold = 1:length(uCells)
        tc   = uCells(fold);
        teM  = strcmp(cell_aug, tc) & row_is_orig;    % test: 해당 셀 원본만
        trM  = ~strcmp(cell_aug, tc);                  % train: 다른 셀 전체 (aug 포함)

        mdl = fitrensemble(X_aug_s(trM,:), y_aug(trM), 'Method','LSBoost',...
            'NumLearningCycles',200,'LearnRate',0.05,...
            'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',3,'Surrogate','on'));

        te_local = find(teM);
        y_pred_te = predict(mdl, X_aug_s(te_local,:));
        % map back to orig index
        orig_in_fold = find(strcmp(ids_c,tc));
        y_cv_orig(orig_in_fold) = y_pred_te;
    end

    err   = y_c - y_cv_orig;
    rmse  = sqrt(mean(err.^2));
    r2    = 1-sum(err.^2)/sum((y_c-mean(y_c)).^2);
    fprintf('%-6.1fC %-7d %-7.3f %-7.4f\n', cr, n_c*(N_mask+1), rmse, r2);

    %-- Train final model on ALL augmented data for this C-rate --
    mdl_final = fitrensemble(X_aug_s, y_aug, 'Method','LSBoost',...
        'NumLearningCycles',200,'LearnRate',0.05,...
        'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',3,'Surrogate','on'));

    cname = matlab.lang.makeValidName(char(cond));
    CRateModels.(cname) = struct('cond',cond,'crate',cr,...
        'mu',mu_c,'sig',sig_c,'mdl',mdl_final,...
        'RMSE',rmse,'R2',r2,'y_orig',y_c,'y_cv',y_cv_orig);
end

%% Save
save(fullfile(ciDir,'crate_models_masked.mat'),...
    'CRateModels','crate_vals','uConds');
fprintf('\nSaved: crate_models_masked.mat\n=== step1_train_masked complete! ===\n');
