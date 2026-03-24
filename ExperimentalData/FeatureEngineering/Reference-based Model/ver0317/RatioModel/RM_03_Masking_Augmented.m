% RM_03_Masking_Augmented.m
% 마스킹 증강 학습: Lab dQ에서 랜덤 연속 세그먼트만 남기고 ratio 계산
% → 모델이 "부분 세그먼트 ratio" 패턴을 학습
% → Field 적용 시 동일 조건

clear; clc; close all;

%% 1. Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
rmDir   = fullfile(verDir, 'RatioModel');
visDir  = fullfile(rmDir, 'Visualization');
rawDir  = fullfile(projDir, 'Rack_raw2mat');
if ~exist(visDir,'dir'), mkdir(visDir); end

pyenv('Version', 'C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

%% 2. Load Lab Data
d = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat'));
FM = d.FM;
mr = load(fullfile(verDir, 'MasterRulers_PseudoOCV.mat'));
V_chg = mr.MasterRuler_ver0317.V_bounds_chg;
V_dch = mr.MasterRuler_ver0317.V_bounds_dch;
Q_nom = 64;

%% 3. Raw dQ (Lab)
dQ_c_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i), 3:12, 'UniformOutput', false);  % 10개
dQ_d_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false);  % 11개
dQ_c_raw = FM{:, dQ_c_cols};
dQ_d_raw = FM{:, dQ_d_cols};
C_eff_c = FM.C_eff_chg;
C_eff_d = FM.C_eff_dch;
y_lab = FM.Static_Capacity / Q_nom * 100;

nRows = size(dQ_c_raw, 1);
nC = size(dQ_c_raw, 2);  % 10
nD = size(dQ_d_raw, 2);  % 11

%% 4. Generate augmented data with masking
rng(42);  % reproducibility
N_aug = 30;  % augmented versions per row

% Pre-generate all contiguous window patterns
% Charge: windows of size 1~10 on 10 segments
chg_patterns = {};
for w = 1:nC
    for s = 1:(nC-w+1)
        chg_patterns{end+1} = s:(s+w-1);
    end
end
% Discharge: windows of size 1~11 on 11 segments
dch_patterns = {};
for w = 1:nD
    for s = 1:(nD-w+1)
        dch_patterns{end+1} = s:(s+w-1);
    end
end

fprintf('Contiguous patterns: charge=%d, discharge=%d\n', length(chg_patterns), length(dch_patterns));

% Build augmented dataset
X_aug = []; y_aug = [];

for i = 1:nRows
    % Original (no masking) → ratio with all segments
    dQ_c_ratio_orig = dQ_c_raw(i,:) / sum(dQ_c_raw(i,:), 'omitnan');
    dQ_d_ratio_orig = dQ_d_raw(i,:) / sum(dQ_d_raw(i,:), 'omitnan');
    X_aug(end+1,:) = [dQ_c_ratio_orig, dQ_d_ratio_orig, C_eff_c(i), C_eff_d(i)];
    y_aug(end+1,1) = y_lab(i);
    
    % N_aug random masking versions
    for a = 1:N_aug
        % Random charge pattern
        cp = chg_patterns{randi(length(chg_patterns))};
        dQ_c_masked = nan(1, nC);
        dQ_c_masked(cp) = dQ_c_raw(i, cp);
        dQ_c_ratio = dQ_c_masked / sum(dQ_c_masked, 'omitnan');
        
        % Random discharge pattern
        dp = dch_patterns{randi(length(dch_patterns))};
        dQ_d_masked = nan(1, nD);
        dQ_d_masked(dp) = dQ_d_raw(i, dp);
        dQ_d_ratio = dQ_d_masked / sum(dQ_d_masked, 'omitnan');
        
        X_aug(end+1,:) = [dQ_c_ratio, dQ_d_ratio, C_eff_c(i), C_eff_d(i)];
        y_aug(end+1,1) = y_lab(i);  % SOH는 동일!
    end
end

fprintf('Augmented dataset: %d rows (original %d × %d augmented)\n', size(X_aug,1), nRows, N_aug+1);
fprintf('NaN ratio in augmented data: %.1f%%\n', sum(isnan(X_aug(:)))/numel(X_aug)*100);

%% 5. Standardize
mu_aug = mean(X_aug, 1, 'omitnan');
sig_aug = std(X_aug, 0, 1, 'omitnan');
sig_aug(sig_aug==0) = 1;
X_aug_s = (X_aug - mu_aug) ./ sig_aug;

% NaN → 0 for tree models that don't handle NaN
X_aug_imp = X_aug_s;
X_aug_imp(isnan(X_aug_imp)) = 0;

%% 6. Train models (XGBoost/LightGBM handle NaN natively)
fprintf('\n=== Training on augmented data ===\n');

% GBM (MATLAB, NaN→0)
Models.GBM = fitrensemble(X_aug_imp, y_aug, 'Method', 'LSBoost', ...
    'NumLearningCycles', 200, 'LearnRate', 0.05, ...
    'Learners', templateTree('MaxNumSplits', 15, 'MinLeafSize', 5, 'Surrogate', 'on'));
fprintf('  GBM trained.\n');

% RF (MATLAB, NaN→0)
Models.RF = fitrensemble(X_aug_imp, y_aug, 'Method', 'Bag', ...
    'NumLearningCycles', 200, ...
    'Learners', templateTree('MaxNumSplits', 20, 'MinLeafSize', 5, 'Surrogate', 'on'));
fprintf('  RF trained.\n');

% XGBoost (NaN passthrough)
Models.XGB = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
    'learning_rate',0.05,'min_child_weight',int32(3),'subsample',0.8,'colsample_bytree',0.8,'random_state',int32(42)));
Models.XGB.fit(py.numpy.array(X_aug_s), py.numpy.array(y_aug));
fprintf('  XGBoost trained.\n');

% LightGBM (NaN passthrough)
Models.LGB = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
    'learning_rate',0.05,'min_child_samples',int32(5),'subsample',0.8,'colsample_bytree',0.8,...
    'random_state',int32(42),'verbose',int32(-1)));
Models.LGB.fit(py.numpy.array(X_aug_s), py.numpy.array(y_aug));
fprintf('  LightGBM trained.\n');

%% 7. Also train baseline (no masking, raw dQ) for comparison
X_base = [dQ_c_raw, dQ_d_raw, C_eff_c, C_eff_d];
mu_base = mean(X_base, 1, 'omitnan');
sig_base = std(X_base, 0, 1, 'omitnan'); sig_base(sig_base==0)=1;
X_base_s = (X_base - mu_base) ./ sig_base;

Models.GBM_base = fitrensemble(X_base_s, y_lab, 'Method', 'LSBoost', ...
    'NumLearningCycles', 200, 'LearnRate', 0.05, ...
    'Learners', templateTree('MaxNumSplits', 15, 'MinLeafSize', 5, 'Surrogate', 'on'));
fprintf('  GBM_baseline trained.\n');

%% 8. Field estimation
dataFiles = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021', 'old', datetime(2021,6,3);
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023', 'new', datetime(2023,10,16);
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024', 'new', datetime(2024,9,9);
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'), 'Y2025', 'new', datetime(2025,7,11);
};

Np = 2; thr_A = Q_nom * 0.05 / Np;
min_chg_sec = [600, 300, 300, 300];
min_dch_sec = [300, 150, 300, 150];
FieldResults = struct();
modelNames = {'GBM','RF','XGBoost','LightGBM'};

for k = 1:size(dataFiles,1)
    fpath=dataFiles{k,1}; yr=dataFiles{k,2}; dtype=dataFiles{k,3}; base_date=dataFiles{k,4};
    if ~exist(fpath,'file'), continue; end
    fprintf('\n=== %s ===\n', yr);

    S = load(fpath);
    if strcmp(dtype,'old')
        D=S.Raw.Rack01;
        if isfield(D,'Time'), t=datetime(D.Time);
        elseif isduration(D.Date_Time), t=base_date+D.Date_Time;
        else, t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent_A(:); V_avg=D.AverageCV_V(:);
        if isfield(D,'SOHPct'), raw_soh=D.SOHPct(:); else, raw_soh=nan(size(I_rack)); end
    else
        D=S.Raw;
        if isduration(D.Date_Time), t=base_date+D.Date_Time; else, t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent(:); V_avg=D.CVavg(:);
        if isfield(D,'SOH_BMS'), raw_soh=D.SOH_BMS(:); else, raw_soh=nan(size(I_rack)); end
    end

    tsec=seconds(t-t(1)); I_cell=I_rack/Np;
    valid_soh=raw_soh; valid_soh(valid_soh<=0)=NaN;
    soh_bms = median(valid_soh, 'omitnan');

    chgSegs=local_find_segments(I_cell>thr_A);
    dchSegs=local_find_segments(I_cell<-thr_A);
    if ~isempty(chgSegs), dur=chgSegs(:,2)-chgSegs(:,1)+1; chgSegs=chgSegs(dur>=min_chg_sec(k),:); end
    if ~isempty(dchSegs), dur=dchSegs(:,2)-dchSegs(:,1)+1; dchSegs=dchSegs(dur>=min_dch_sec(k),:); end
    if isempty(chgSegs), continue; end

    % Charge V-Q (no smoothing)
    [~,bi]=max(chgSegs(:,2)-chgSegs(:,1));
    cs=chgSegs(bi,1); ce=chgSegs(bi,2);
    v_c=V_avg(cs:ce); i_c=I_cell(cs:ce); t_c=tsec(cs:ce);
    q_c=cumtrapz(t_c,abs(i_c))/3600;
    v_mono=v_c(1); q_mono=q_c(1);
    for ii=2:length(v_c), if v_c(ii)>v_mono(end), v_mono(end+1)=v_c(ii); q_mono(end+1)=q_c(ii); end; end
    dQ_chg=nan(1,12);
    if length(v_mono)>1, QR=interp1(v_mono,q_mono,V_chg,'linear',NaN); dQ_chg=abs(diff(QR)); end

    % Discharge V-Q (no smoothing)
    dQ_dch=nan(1,12); C_eff_d_val=NaN;
    if ~isempty(dchSegs)
        [~,bi_d]=max(dchSegs(:,2)-dchSegs(:,1));
        ds=dchSegs(bi_d,1); de=dchSegs(bi_d,2);
        v_d=V_avg(ds:de); i_d=I_cell(ds:de); t_d=tsec(ds:de);
        q_d=cumtrapz(t_d,abs(i_d))/3600;
        v_md=v_d(1); q_md=q_d(1);
        for ii=2:length(v_d), if v_d(ii)<v_md(end), v_md(end+1)=v_d(ii); q_md(end+1)=q_d(ii); end; end
        v_asc=flip(v_md); q_asc=flip(q_md);
        if length(v_asc)>1, QR_d=interp1(v_asc,q_asc,V_dch,'linear',NaN); dQ_dch=abs(diff(QR_d)); end
        C_eff_d_val=mean(abs(i_d))/Q_nom;
    end
    C_eff_c_val=mean(abs(i_c))/Q_nom;

    dQ_c_f = dQ_chg(3:12); dQ_d_f = dQ_dch(1:11);

    % Masking-augmented model: use dQ_ratio
    dQ_c_ratio_f = dQ_c_f / sum(dQ_c_f, 'omitnan');
    dQ_d_ratio_f = dQ_d_f / sum(dQ_d_f, 'omitnan');
    X_field_ratio = [dQ_c_ratio_f, dQ_d_ratio_f, C_eff_c_val, C_eff_d_val];
    X_fr_s = (X_field_ratio - mu_aug) ./ sig_aug;

    % Baseline: raw dQ
    X_field_raw = [dQ_c_f, dQ_d_f, C_eff_c_val, C_eff_d_val];
    X_fb_s = (X_field_raw - mu_base) ./ sig_base;

    fprintf('  BMS=%.1f%%, Chg segs: %s, Dch segs: %s\n', soh_bms, ...
        mat2str(find(~isnan(dQ_c_f))), mat2str(find(~isnan(dQ_d_f))));

    % Predict: masking-augmented models
    X_imp = X_fr_s; X_imp(isnan(X_imp))=0;
    soh_gbm = predict(Models.GBM, X_imp);
    soh_rf = predict(Models.RF, X_imp);
    X_py = py.numpy.array(X_fr_s).reshape(int32(1),int32(-1));
    soh_xgb = double(Models.XGB.predict(X_py));
    soh_lgb = double(Models.LGB.predict(X_py));

    % Predict: baseline (raw dQ, no masking)
    X_base_imp = X_fb_s; X_base_imp(isnan(X_base_imp))=0;
    soh_base = predict(Models.GBM_base, X_base_imp);

    FieldResults.(yr).BMS = soh_bms;
    FieldResults.(yr).GBM_base = soh_base;
    FieldResults.(yr).GBM = soh_gbm;
    FieldResults.(yr).RF = soh_rf;
    FieldResults.(yr).XGBoost = soh_xgb;
    FieldResults.(yr).LightGBM = soh_lgb;

    fprintf('  Baseline GBM(raw):     %.1f%%\n', soh_base);
    fprintf('  Masked   GBM(ratio):   %.1f%%\n', soh_gbm);
    fprintf('  Masked   RF(ratio):    %.1f%%\n', soh_rf);
    fprintf('  Masked   XGB(ratio):   %.1f%%\n', soh_xgb);
    fprintf('  Masked   LGB(ratio):   %.1f%%\n', soh_lgb);
end

%% 9. Summary
fprintf('\n\n========== MASKING AUGMENTED vs BASELINE (SOH%%) ==========\n');
fprintf('%-6s  %-5s  %-10s  %-8s %-8s %-8s %-8s\n', 'Year','BMS','Base(raw)','GBM_M','RF_M','XGB_M','LGB_M');
yrs = fieldnames(FieldResults);
for k = 1:length(yrs)
    r = FieldResults.(yrs{k});
    fprintf('%-6s  %-5.1f  %-10.1f  %-8.1f %-8.1f %-8.1f %-8.1f\n', ...
        yrs{k}, r.BMS, r.GBM_base, r.GBM, r.RF, r.XGBoost, r.LightGBM);
end

%% 10. Visualization
fig = figure('Position', [50 100 1100 500]);
n = length(yrs); x = 1:n;
bms = arrayfun(@(i) FieldResults.(yrs{i}).BMS, 1:n);
base_arr = arrayfun(@(i) FieldResults.(yrs{i}).GBM_base, 1:n);
gbm_arr = arrayfun(@(i) FieldResults.(yrs{i}).GBM, 1:n);
rf_arr = arrayfun(@(i) FieldResults.(yrs{i}).RF, 1:n);
xgb_arr = arrayfun(@(i) FieldResults.(yrs{i}).XGBoost, 1:n);
lgb_arr = arrayfun(@(i) FieldResults.(yrs{i}).LightGBM, 1:n);

hold on; grid on; box on;
plot(x, bms, 'k--s', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'BMS SOH');
plot(x, base_arr, '-o', 'Color', [0.7 0.2 0.2], 'LineWidth', 1.5, 'MarkerSize', 9, 'DisplayName', 'Baseline (raw dQ)');
plot(x, gbm_arr, '-d', 'Color', [0.9 0.5 0.1], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.9 0.5 0.1], 'DisplayName', 'Masked GBM');
plot(x, rf_arr, '-^', 'Color', [0.1 0.7 0.3], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.1 0.7 0.3], 'DisplayName', 'Masked RF');
plot(x, xgb_arr, '-v', 'Color', [0.8 0.2 0.2], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.2 0.2], 'DisplayName', 'Masked XGBoost');
plot(x, lgb_arr, '-p', 'Color', [0.6 0.2 0.8], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.6 0.2 0.8], 'DisplayName', 'Masked LightGBM');

for i=1:n
    text(x(i)+0.05, bms(i)+0.3, sprintf('%.1f', bms(i)), 'FontSize', 8, 'FontWeight', 'bold');
    text(x(i)+0.05, base_arr(i)-0.8, sprintf('%.1f', base_arr(i)), 'FontSize', 7, 'Color', [0.7 0.2 0.2]);
    text(x(i)+0.05, gbm_arr(i)+0.3, sprintf('%.1f', gbm_arr(i)), 'FontSize', 7, 'Color', [0.9 0.5 0.1]);
end

set(gca, 'XTick', x, 'XTickLabel', strrep(yrs, 'Y', ''), 'FontSize', 12);
xlabel('Year', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('SOH (%)', 'FontSize', 13, 'FontWeight', 'bold');
title('Masking-Augmented dQ Ratio vs Baseline (Raw dQ)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 10);
ylim([80 105]);
saveas(fig, fullfile(visDir, 'RM_Masking_Augmented_Comparison.png'));
fprintf('\nSaved: RM_Masking_Augmented_Comparison.png\n');

%% 11. Save
save(fullfile(rmDir, 'RM_Masking_Results.mat'), 'FieldResults', 'Models', 'mu_aug', 'sig_aug');
fprintf('=== Masking Augmented Training Complete! ===\n');

%% Helper
function segs = local_find_segments(mask)
    segs=[]; n=length(mask); i=1;
    while i<=n
        if mask(i), j=i;
            while j<n && mask(j+1), j=j+1; end
            segs=[segs;i,j]; i=j+1;
        else, i=i+1; end
    end
end
