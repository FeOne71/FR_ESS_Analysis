% DM_02_Field_Estimation.m
% Discharge-Only Field SOH Estimation (2021~2025)
% Full 12 features (dQ_d_01~11 + C_eff_dch), NaN passthrough

clear; clc; close all;

%% 1. Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
dmDir   = fullfile(verDir, 'DischargeModel');
visDir  = fullfile(dmDir, 'Visualization');
rawDir  = fullfile(projDir, 'Rack_raw2mat');
if ~exist(visDir,'dir'), mkdir(visDir); end

pyenv('Version', 'C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

%% 2. Load Lab Data & Master Ruler
d = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat'));
FM = d.FM;
mr = load(fullfile(verDir, 'MasterRulers_PseudoOCV.mat'));
V_dch = mr.MasterRuler_ver0317.V_bounds_dch;
num_segs = 12;

%% 3. Discharge features + SOH%
Q_nom = 64;
feat_names = [arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false), {'C_eff_dch'}];
X_lab = FM{:, feat_names};
y_lab = FM.Static_Capacity / Q_nom * 100;

mu_lab = mean(X_lab, 1);
sigma_lab = std(X_lab, 0, 1);
sigma_lab(sigma_lab == 0) = 1;
X_lab_s = (X_lab - mu_lab) ./ sigma_lab;

%% 4. Train 5 models on full Lab data
fprintf('=== Training Discharge models (%d rows, %d features) ===\n', size(X_lab,1), length(feat_names));

[B, FI] = lasso(X_lab_s, y_lab, 'CV', 4);
idxL = FI.IndexMinMSE;
Models.LASSO.coef = B(:,idxL); Models.LASSO.intercept = FI.Intercept(idxL);
fprintf('  LASSO trained.\n');

Models.RF = fitrensemble(X_lab_s, y_lab, 'Method', 'Bag', ...
    'NumLearningCycles', 100, 'Learners', templateTree('MaxNumSplits', 20, 'MinLeafSize', 5, 'Surrogate', 'on'));
fprintf('  RF trained.\n');

Models.GBM = fitrensemble(X_lab_s, y_lab, 'Method', 'LSBoost', ...
    'NumLearningCycles', 100, 'LearnRate', 0.1, ...
    'Learners', templateTree('MaxNumSplits', 10, 'MinLeafSize', 5, 'Surrogate', 'on'));
fprintf('  GBM trained.\n');

Models.XGB = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
    'learning_rate',0.1,'min_child_weight',int32(3),'subsample',0.8,'colsample_bytree',0.8,'random_state',int32(42)));
Models.XGB.fit(py.numpy.array(X_lab_s), py.numpy.array(y_lab));
fprintf('  XGBoost trained.\n');

Models.LGB = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
    'learning_rate',0.1,'min_child_samples',int32(5),'subsample',0.8,'colsample_bytree',0.8,...
    'random_state',int32(42),'verbose',int32(-1)));
Models.LGB.fit(py.numpy.array(X_lab_s), py.numpy.array(y_lab));
fprintf('  LightGBM trained.\n');

%% 5. Field Data
dataFiles = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021', 'old', datetime(2021,6,3);
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023', 'new', datetime(2023,10,16);
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024', 'new', datetime(2024,9,9);
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'), 'Y2025', 'new', datetime(2025,7,11);
};

Np = 2; thr_A = Q_nom * 0.05 / Np;
min_dch_sec = [300, 150, 300, 150];
ma_window = 30;

%% 6. Process
FieldResults = struct();
modelNames = {'LASSO','RF','GBM','XGBoost','LightGBM'};

for k = 1:size(dataFiles, 1)
    fpath=dataFiles{k,1}; yr=dataFiles{k,2}; dtype=dataFiles{k,3}; base_date=dataFiles{k,4};
    if ~exist(fpath,'file'), fprintf('[SKIP] %s\n', yr); continue; end
    fprintf('\n=== Processing %s ===\n', yr);

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

    % Find discharge segments
    dchSegs = local_find_segments(I_cell < -thr_A);
    if ~isempty(dchSegs)
        dur = dchSegs(:,2)-dchSegs(:,1)+1;
        dchSegs = dchSegs(dur >= min_dch_sec(k), :);
    end

    if isempty(dchSegs), fprintf('  No discharge segments!\n'); continue; end

    % Longest discharge segment
    [~,bi]=max(dchSegs(:,2)-dchSegs(:,1));
    dch_s=dchSegs(bi,1); dch_e=dchSegs(bi,2);
    v_d=V_avg(dch_s:dch_e); i_d=I_cell(dch_s:dch_e); t_d=tsec(dch_s:dch_e);

    V_min_d=min(v_d); V_max_d=max(v_d);

    % --- Discharge V-Q (no smoothing, Lab-like pipeline) ---
    q_d = cumtrapz(t_d, abs(i_d)) / 3600;
    % Monotonic decreasing V filter → flip to ascending for interp1
    v_mono_d = v_d(1); q_mono_d = q_d(1);
    for ii = 2:length(v_d)
        if v_d(ii) < v_mono_d(end)
            v_mono_d(end+1) = v_d(ii); q_mono_d(end+1) = q_d(ii);
        end
    end
    v_asc = flip(v_mono_d); q_asc = flip(q_mono_d);

    % dQ discharge features: direct interp1 at boundary voltages
    dQ_dch = nan(1, num_segs);
    if length(v_asc) > 1
        QR_d = interp1(v_asc, q_asc, V_dch, 'linear', NaN);
        dQ_dch = abs(diff(QR_d));
    end

    C_eff_dch = mean(abs(i_d)) / Q_nom;

    % Build feature vector: dQ_d_01~11, C_eff_dch (12 features)
    X_field = [dQ_dch(1:11), C_eff_dch];
    nan_mask = isnan(X_field);

    fprintf('  Dch V range: [%.3f, %.3f]V, C_eff_dch=%.3f\n', V_min_d, V_max_d, C_eff_dch);
    fprintf('  Valid dch segs: %s\n', mat2str(find(~isnan(dQ_dch))));
    fprintf('  NaN features: %d / %d\n', sum(nan_mask), length(X_field));
    fprintf('  BMS SOH: %.1f%%\n', soh_bms);

    % Standardize (NaN stays NaN)
    X_field_s = (X_field - mu_lab) ./ sigma_lab;

    % LASSO/RF/GBM: impute NaN with 0 (= training mean in z-space)
    X_imp = X_field_s; X_imp(isnan(X_imp)) = 0;
    soh_lasso = X_imp * Models.LASSO.coef + Models.LASSO.intercept;
    soh_rf    = predict(Models.RF, X_imp);
    soh_gbm   = predict(Models.GBM, X_imp);

    % XGBoost/LightGBM: NaN passthrough
    X_py = py.numpy.array(X_field_s).reshape(int32(1),int32(-1));
    soh_xgb = double(Models.XGB.predict(X_py));
    soh_lgb = double(Models.LGB.predict(X_py));

    FieldResults.(yr).BMS_SOH  = soh_bms;
    FieldResults.(yr).LASSO    = soh_lasso;
    FieldResults.(yr).RF       = soh_rf;
    FieldResults.(yr).GBM      = soh_gbm;
    FieldResults.(yr).XGBoost  = soh_xgb;
    FieldResults.(yr).LightGBM = soh_lgb;
    FieldResults.(yr).C_eff    = C_eff_dch;
    FieldResults.(yr).NaN_cnt  = sum(nan_mask);

    fprintf('  SOH%%: LASSO=%.1f, RF=%.1f, GBM=%.1f, XGB=%.1f, LGB=%.1f\n', ...
        soh_lasso, soh_rf, soh_gbm, soh_xgb, soh_lgb);
end

%% 7. Summary
fprintf('\n\n========== DISCHARGE MODEL - FIELD SOH (SOH %%) ==========\n');
fprintf('%-6s  %-6s  %-7s %-7s %-7s %-7s %-7s  %-6s %-3s\n', ...
    'Year','BMS','LASSO','RF','GBM','XGB','LGB','Ceff_d','NaN');
yrs = fieldnames(FieldResults);
for k = 1:length(yrs)
    r = FieldResults.(yrs{k});
    fprintf('%-6s  %-6.1f  %-7.1f %-7.1f %-7.1f %-7.1f %-7.1f  %-6.3f %-3d\n', ...
        yrs{k}, r.BMS_SOH, r.LASSO, r.RF, r.GBM, r.XGBoost, r.LightGBM, r.C_eff, r.NaN_cnt);
end

%% 8. Save
save(fullfile(dmDir, 'DM_FieldResults.mat'), 'FieldResults', 'Models', 'mu_lab', 'sigma_lab', 'feat_names');

%% 9. Visualization
if ~isempty(fieldnames(FieldResults))
    fig = figure('Position', [50 100 1200 500]);
    n = length(yrs); x = 1:n;
    soh_bms_arr = arrayfun(@(i) FieldResults.(yrs{i}).BMS_SOH, 1:n);
    soh_pct = zeros(n,5);
    for i=1:n
        r=FieldResults.(yrs{i});
        soh_pct(i,:)=[r.LASSO, r.RF, r.GBM, r.XGBoost, r.LightGBM];
    end

    hold on; grid on; box on;
    plot(x, soh_bms_arr, 'k--s', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'BMS SOH');
    colors = [0.2 0.4 0.85; 0.1 0.7 0.3; 0.9 0.5 0.1; 0.8 0.2 0.2; 0.6 0.2 0.8];
    markers = {'o','d','^','v','p'};
    for m=1:5
        plot(x, soh_pct(:,m), ['-' markers{m}], 'Color', colors(m,:), 'LineWidth', 1.8, ...
            'MarkerSize', 9, 'MarkerFaceColor', colors(m,:), 'DisplayName', modelNames{m});
    end
    for m=1:5
        for i=1:n
            text(x(i)+0.03, soh_pct(i,m)+0.3, sprintf('%.1f%%', soh_pct(i,m)), 'FontSize', 7, 'Color', colors(m,:));
        end
    end
    for i=1:n
        text(x(i)+0.03, soh_bms_arr(i)+0.3, sprintf('%.1f%%', soh_bms_arr(i)), 'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold');
    end
    set(gca, 'XTick', x, 'XTickLabel', strrep(yrs, 'Y', ''), 'FontSize', 12);
    xlabel('Measurement Year', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('SOH (%)', 'FontSize', 13, 'FontWeight', 'bold');
    title('Discharge-Only Model: Field SOH vs BMS (2021-2025)', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'southwest', 'FontSize', 10);
    ylim([40 105]);
    saveas(fig, fullfile(visDir, 'DM_Field_SOH_Comparison.png'));
    fprintf('\nSaved: DM_Field_SOH_Comparison.png\n');
end
fprintf('=== Discharge Field Estimation Complete! ===\n');

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
