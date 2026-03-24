% ML_02_Field_Estimation.m (v2)
% 수정사항:
%   1. Full 23개 피처로 훈련, 필드에서 NaN은 그대로 전달
%   2. 타겟을 SOH% (= Static_Capacity / 64 * 100)으로 변환
%   3. XGBoost/LightGBM은 NaN 직접 처리, MATLAB 모델은 train mean으로 imputation

clear; clc; close all;

%% 1. Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
rawDir  = fullfile(projDir, 'Rack_raw2mat');
visDir  = fullfile(verDir, 'Visualization');
if ~exist(visDir,'dir'), mkdir(visDir); end

pyenv('Version', 'C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

%% 2. Load Lab Data & Master Ruler
d = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat'));
FM = d.FM;
mr = load(fullfile(verDir, 'MasterRulers_PseudoOCV.mat'));
V_chg = mr.MasterRuler_ver0317.V_bounds_chg;
V_dch = mr.MasterRuler_ver0317.V_bounds_dch;
num_segs = 12;

%% 3. Full feature set (23 features) + SOH% target
feat_names = [arrayfun(@(i) sprintf('dQ_c_%02d',i), 3:12, 'UniformOutput', false), ...
              arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false), ...
              {'C_eff_chg', 'C_eff_dch'}];
Q_nom = 64;  % Nominal capacity (Ah)

X_lab = FM{:, feat_names};
y_lab_pct = FM.Static_Capacity / Q_nom * 100;  % SOH%

% Standardize
mu_lab = mean(X_lab, 1);
sigma_lab = std(X_lab, 0, 1);
sigma_lab(sigma_lab == 0) = 1;
X_lab_s = (X_lab - mu_lab) ./ sigma_lab;

%% 4. Train 5 models on full Lab data (target = SOH%)
fprintf('=== Training models (target=SOH%%, %d rows, %d features) ===\n', size(X_lab,1), size(X_lab,2));

% LASSO
[B, FitInfo] = lasso(X_lab_s, y_lab_pct, 'CV', 4);
idxL = FitInfo.IndexMinMSE;
Models.LASSO.coef = B(:, idxL);
Models.LASSO.intercept = FitInfo.Intercept(idxL);
fprintf('  LASSO trained.\n');

% RF (with Surrogate for NaN handling)
Models.RF = fitrensemble(X_lab_s, y_lab_pct, 'Method', 'Bag', ...
    'NumLearningCycles', 100, ...
    'Learners', templateTree('MaxNumSplits', 20, 'MinLeafSize', 5, 'Surrogate', 'on'));
fprintf('  RF trained.\n');

% GBM
Models.GBM = fitrensemble(X_lab_s, y_lab_pct, 'Method', 'LSBoost', ...
    'NumLearningCycles', 100, 'LearnRate', 0.1, ...
    'Learners', templateTree('MaxNumSplits', 10, 'MinLeafSize', 5, 'Surrogate', 'on'));
fprintf('  GBM trained.\n');

% XGBoost (handles NaN natively)
Models.XGB = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
    'learning_rate',0.1,'min_child_weight',int32(3),'subsample',0.8,'colsample_bytree',0.8,'random_state',int32(42)));
Models.XGB.fit(py.numpy.array(X_lab_s), py.numpy.array(y_lab_pct));
fprintf('  XGBoost trained.\n');

% LightGBM (handles NaN natively)
Models.LGB = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
    'learning_rate',0.1,'min_child_samples',int32(5),'subsample',0.8,'colsample_bytree',0.8,...
    'random_state',int32(42),'verbose',int32(-1)));
Models.LGB.fit(py.numpy.array(X_lab_s), py.numpy.array(y_lab_pct));
fprintf('  LightGBM trained.\n');

%% 5. Field Data
dataFiles = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021', 'old', datetime(2021,6,3);
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023', 'new', datetime(2023,10,16);
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024', 'new', datetime(2024,9,9);
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'), 'Y2025', 'new', datetime(2025,7,11);
};

Np = 2; thr_A = Q_nom * 0.05 / Np;
min_chg_sec = [600, 300, 300, 300];
min_dch_sec = [300, 150, 300, 150];
ma_window = 60;

%% 6. Process each field year
FieldResults = struct();
modelNames = {'LASSO','RF','GBM','XGBoost','LightGBM'};

for k = 1:size(dataFiles, 1)
    fpath = dataFiles{k,1}; yr = dataFiles{k,2};
    dtype = dataFiles{k,3}; base_date = dataFiles{k,4};
    if ~exist(fpath,'file'), fprintf('[SKIP] %s\n', yr); continue; end
    fprintf('\n=== Processing %s ===\n', yr);

    S = load(fpath);
    if strcmp(dtype,'old')
        D = S.Raw.Rack01;
        if isfield(D,'Time'), t=datetime(D.Time);
        elseif isduration(D.Date_Time), t=base_date+D.Date_Time;
        else, t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent_A(:); V_avg=D.AverageCV_V(:);
        if isfield(D,'SOHPct'), raw_soh=D.SOHPct(:); else, raw_soh=nan(size(I_rack)); end
    else
        D = S.Raw;
        if isduration(D.Date_Time), t=base_date+D.Date_Time; else, t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent(:); V_avg=D.CVavg(:);
        if isfield(D,'SOH_BMS'), raw_soh=D.SOH_BMS(:); else, raw_soh=nan(size(I_rack)); end
    end

    tsec=seconds(t-t(1)); I_cell=I_rack/Np;
    valid_soh=raw_soh; valid_soh(valid_soh<=0)=NaN;
    soh_bms = median(valid_soh, 'omitnan');

    % Find charge/discharge segments
    chgSegs=local_find_segments(I_cell>thr_A);
    dchSegs=local_find_segments(I_cell<-thr_A);
    if ~isempty(chgSegs), dur=chgSegs(:,2)-chgSegs(:,1)+1; chgSegs=chgSegs(dur>=min_chg_sec(k),:); end
    if ~isempty(dchSegs), dur=dchSegs(:,2)-dchSegs(:,1)+1; dchSegs=dchSegs(dur>=min_dch_sec(k),:); end
    if isempty(chgSegs), fprintf('  No charge segs!\n'); continue; end

    % Longest charge segment
    [~,bi]=max(chgSegs(:,2)-chgSegs(:,1));
    chg_s=chgSegs(bi,1); chg_e=chgSegs(bi,2);
    v_c=V_avg(chg_s:chg_e); i_c=I_cell(chg_s:chg_e); t_c=tsec(chg_s:chg_e);

    % --- Charge V-Q (no smoothing, Lab-like pipeline) ---
    q_c = cumtrapz(t_c, abs(i_c)) / 3600;
    % Monotonic increasing V filter
    v_mono = v_c(1); q_mono = q_c(1);
    for ii = 2:length(v_c)
        if v_c(ii) > v_mono(end)
            v_mono(end+1) = v_c(ii); q_mono(end+1) = q_c(ii);
        end
    end

    % dQ charge features: direct interp1 at boundary voltages
    dQ_chg = nan(1, num_segs);
    V_min_c = min(v_c); V_max_c = max(v_c);
    if length(v_mono) > 1
        QR = interp1(v_mono, q_mono, V_chg, 'linear', NaN);
        dQ_chg = abs(diff(QR));
    end

    % --- Discharge V-Q (no smoothing, Lab-like pipeline) ---
    dQ_dch = nan(1, num_segs);
    if ~isempty(dchSegs)
        [~,bi_d]=max(dchSegs(:,2)-dchSegs(:,1));
        dch_s=dchSegs(bi_d,1); dch_e=dchSegs(bi_d,2);
        v_d=V_avg(dch_s:dch_e); i_d=I_cell(dch_s:dch_e); t_d=tsec(dch_s:dch_e);

        q_d = cumtrapz(t_d, abs(i_d)) / 3600;
        % Monotonic decreasing V filter → flip to ascending for interp1
        v_mono_d = v_d(1); q_mono_d = q_d(1);
        for ii = 2:length(v_d)
            if v_d(ii) < v_mono_d(end)
                v_mono_d(end+1) = v_d(ii); q_mono_d(end+1) = q_d(ii);
            end
        end
        v_asc = flip(v_mono_d); q_asc = flip(q_mono_d);

        if length(v_asc) > 1
            QR_d = interp1(v_asc, q_asc, V_dch, 'linear', NaN);
            dQ_dch = abs(diff(QR_d));
        end
    end

    % C_eff
    C_eff_chg = mean(abs(i_c)) / Q_nom;
    if ~isempty(dchSegs), C_eff_dch = mean(abs(I_cell(dch_s:dch_e))) / Q_nom;
    else, C_eff_dch = NaN; end

    % Build FULL feature vector (23 features, NaN where unavailable)
    X_field = [dQ_chg(3:12), dQ_dch(1:11), C_eff_chg, C_eff_dch];
    nan_mask = isnan(X_field);
    
    fprintf('  V range: [%.3f, %.3f]V, C_eff=%.3f\n', V_min_c, V_max_c, C_eff_chg);
    fprintf('  Valid segs (chg): %s\n', mat2str(find(~isnan(dQ_chg))));
    fprintf('  NaN features: %d / %d\n', sum(nan_mask), length(X_field));
    fprintf('  BMS SOH: %.1f%%\n', soh_bms);

    % Standardize (NaN positions stay NaN)
    X_field_s = (X_field - mu_lab) ./ sigma_lab;  % NaN stays NaN

    % --- LASSO/RF/GBM: impute NaN with 0 (standardized space, 0 = mean) ---
    X_imputed = X_field_s;
    X_imputed(isnan(X_imputed)) = 0;  % 0 in standardized space = training mean
    
    soh_lasso = X_imputed * Models.LASSO.coef + Models.LASSO.intercept;
    soh_rf    = predict(Models.RF, X_imputed);
    soh_gbm   = predict(Models.GBM, X_imputed);

    % --- XGBoost/LightGBM: pass NaN directly (native handling) ---
    X_py_nan = py.numpy.array(X_field_s).reshape(int32(1), int32(-1));
    soh_xgb = double(Models.XGB.predict(X_py_nan));
    soh_lgb = double(Models.LGB.predict(X_py_nan));

    FieldResults.(yr).BMS_SOH  = soh_bms;
    FieldResults.(yr).LASSO    = soh_lasso;
    FieldResults.(yr).RF       = soh_rf;
    FieldResults.(yr).GBM      = soh_gbm;
    FieldResults.(yr).XGBoost  = soh_xgb;
    FieldResults.(yr).LightGBM = soh_lgb;
    FieldResults.(yr).C_eff    = C_eff_chg;
    FieldResults.(yr).NaN_cnt  = sum(nan_mask);

    fprintf('  SOH%%: LASSO=%.1f, RF=%.1f, GBM=%.1f, XGB=%.1f, LGB=%.1f\n', ...
        soh_lasso, soh_rf, soh_gbm, soh_xgb, soh_lgb);
end

%% 7. Summary
fprintf('\n\n========== FIELD SOH ESTIMATION (SOH %%) ==========\n');
fprintf('%-6s  %-6s  %-7s %-7s %-7s %-7s %-7s  %-5s %-3s\n', ...
    'Year','BMS','LASSO','RF','GBM','XGB','LGB','Ceff','NaN');
yrs = fieldnames(FieldResults);
for k = 1:length(yrs)
    r = FieldResults.(yrs{k});
    fprintf('%-6s  %-6.1f  %-7.1f %-7.1f %-7.1f %-7.1f %-7.1f  %-5.3f %-3d\n', ...
        yrs{k}, r.BMS_SOH, r.LASSO, r.RF, r.GBM, r.XGBoost, r.LightGBM, r.C_eff, r.NaN_cnt);
end

%% 8. Save
save(fullfile(verDir, 'ML_02_FieldResults.mat'), 'FieldResults', 'Models', 'mu_lab', 'sigma_lab', 'feat_names');

%% 9. Visualization
fig = figure('Position', [50 100 1200 500]);
n = length(yrs); x = 1:n;
soh_bms_arr = arrayfun(@(i) FieldResults.(yrs{i}).BMS_SOH, 1:n);
soh_pct = zeros(n,5);
for i = 1:n
    r = FieldResults.(yrs{i});
    soh_pct(i,:) = [r.LASSO, r.RF, r.GBM, r.XGBoost, r.LightGBM];
end

hold on; grid on; box on;
plot(x, soh_bms_arr, 'k--s', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'BMS SOH');
colors = [0.2 0.4 0.85; 0.1 0.7 0.3; 0.9 0.5 0.1; 0.8 0.2 0.2; 0.6 0.2 0.8];
markers = {'o','d','^','v','p'};
for m = 1:5
    plot(x, soh_pct(:,m), ['-' markers{m}], 'Color', colors(m,:), 'LineWidth', 1.8, ...
        'MarkerSize', 9, 'MarkerFaceColor', colors(m,:), 'DisplayName', modelNames{m});
end
for m = 1:5
    for i = 1:n
        text(x(i)+0.03, soh_pct(i,m)+0.3, sprintf('%.1f%%', soh_pct(i,m)), ...
            'FontSize', 7, 'Color', colors(m,:));
    end
end
for i = 1:n
    text(x(i)+0.03, soh_bms_arr(i)+0.3, sprintf('%.0f%%', soh_bms_arr(i)), ...
        'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold');
end

set(gca, 'XTick', x, 'XTickLabel', strrep(yrs, 'Y', ''), 'FontSize', 12);
xlabel('Measurement Year', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('SOH (%)', 'FontSize', 13, 'FontWeight', 'bold');
title('Field SOH Estimation (Full Features + NaN Passthrough) vs BMS', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 10);
ylim([40 105]);

saveas(fig, fullfile(visDir, 'ML_Field_SOH_Comparison.png'));
fprintf('\nSaved: ML_Field_SOH_Comparison.png\n');
fprintf('=== Complete! ===\n');

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
