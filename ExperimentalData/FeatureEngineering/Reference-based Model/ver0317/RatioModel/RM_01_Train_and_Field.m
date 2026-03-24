% RM_01_Train_and_Field.m
% dQ Ratio Model: 프로토콜 불변 피처로 Lab 학습 → Field 적용
% 핵심: dQ_ratio(s) = dQ(s) / sum(valid dQ) → 절대 스케일 제거, 형태만 학습
% 
% 충전 피처: dQ_c_ratio_03~12 (10개) + C_eff_chg
% 방전 피처: dQ_d_ratio_01~11 (11개) + C_eff_dch
% 총 23개 피처

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
num_segs = 12;
Q_nom = 64;

%% 3. Convert Lab dQ → dQ_ratio
% Raw dQ columns
dQ_c_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i), 3:12, 'UniformOutput', false);
dQ_d_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false);

dQ_c_raw = FM{:, dQ_c_cols};  % 215 x 10
dQ_d_raw = FM{:, dQ_d_cols};  % 215 x 11
C_eff_chg = FM.C_eff_chg;
C_eff_dch = FM.C_eff_dch;

% dQ_ratio: each dQ / sum of valid dQ in same row
dQ_c_ratio = dQ_c_raw ./ sum(dQ_c_raw, 2, 'omitnan');
dQ_d_ratio = dQ_d_raw ./ sum(dQ_d_raw, 2, 'omitnan');

X_lab = [dQ_c_ratio, dQ_d_ratio, C_eff_chg, C_eff_dch];
y_lab = FM.Static_Capacity / Q_nom * 100;  % SOH%

feat_names = [arrayfun(@(i) sprintf('dQr_c_%02d',i), 3:12, 'UniformOutput', false), ...
              arrayfun(@(i) sprintf('dQr_d_%02d',i), 1:11, 'UniformOutput', false), ...
              {'C_eff_chg', 'C_eff_dch'}];

fprintf('Lab dQ_ratio stats (charge):\n');
fprintf('  Mean ratio per seg: %s\n', mat2str(mean(dQ_c_ratio,1,'omitnan'), 3));
fprintf('  Std:                %s\n', mat2str(std(dQ_c_ratio,0,1,'omitnan'), 3));

% Standardize
mu_lab = mean(X_lab, 1, 'omitnan');
sigma_lab = std(X_lab, 0, 1, 'omitnan');
sigma_lab(sigma_lab == 0) = 1;
X_lab_s = (X_lab - mu_lab) ./ sigma_lab;

%% 4. LOCO-CV
cellIDs = unique(FM.CellID);
nCells = length(cellIDs);
modelNames = {'LASSO','RF','GBM','XGBoost','LightGBM'};
nModels = 5;

fprintf('\n=== dQ Ratio Model LOCO-CV (23 ratio features) ===\n');
Results = table();

for m = 1:nModels
    mName = modelNames{m};
    fprintf('  %s ... ', mName);
    y_true_all = []; y_pred_all = [];
    
    for fold = 1:nCells
        testIdx = FM.CellID == cellIDs(fold);
        trainIdx = ~testIdx;
        X_tr = X_lab_s(trainIdx,:); y_tr = y_lab(trainIdx);
        X_te = X_lab_s(testIdx,:);  y_te = y_lab(testIdx);
        
        switch mName
            case 'LASSO'
                [B,FI] = lasso(X_tr, y_tr, 'CV', 4);
                y_p = X_te * B(:,FI.IndexMinMSE) + FI.Intercept(FI.IndexMinMSE);
            case 'RF'
                mdl = fitrensemble(X_tr, y_tr, 'Method', 'Bag', 'NumLearningCycles', 100, ...
                    'Learners', templateTree('MaxNumSplits', 20, 'MinLeafSize', 5, 'Surrogate', 'on'));
                y_p = predict(mdl, X_te);
            case 'GBM'
                mdl = fitrensemble(X_tr, y_tr, 'Method', 'LSBoost', 'NumLearningCycles', 100, ...
                    'LearnRate', 0.1, 'Learners', templateTree('MaxNumSplits', 10, 'MinLeafSize', 5, 'Surrogate', 'on'));
                y_p = predict(mdl, X_te);
            case 'XGBoost'
                xgb = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
                    'learning_rate',0.1,'min_child_weight',int32(3),'subsample',0.8,'colsample_bytree',0.8,'random_state',int32(42)));
                xgb.fit(py.numpy.array(X_tr), py.numpy.array(y_tr));
                y_p = double(xgb.predict(py.numpy.array(X_te)));
            case 'LightGBM'
                lgb = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
                    'learning_rate',0.1,'min_child_samples',int32(5),'subsample',0.8,'colsample_bytree',0.8,...
                    'random_state',int32(42),'verbose',int32(-1)));
                lgb.fit(py.numpy.array(X_tr), py.numpy.array(y_tr));
                y_p = double(lgb.predict(py.numpy.array(X_te)));
        end
        y_true_all = [y_true_all; y_te]; y_pred_all = [y_pred_all; y_p(:)];
    end
    
    err = y_true_all - y_pred_all;
    rmse = sqrt(mean(err.^2)); mape = mean(abs(err./y_true_all))*100;
    r2 = 1 - sum(err.^2)/sum((y_true_all-mean(y_true_all)).^2);
    fprintf('RMSE=%.3f%%, R2=%.4f\n', rmse, r2);
    Results(m,:) = table({mName}, rmse, mape, r2, 'VariableNames', {'Model','RMSE','MAPE','R2'});
end

disp(Results);
writetable(Results, fullfile(rmDir, 'RM_LOCO_CV_Results.csv'));

%% 5. Train final models on full Lab data
fprintf('\n=== Training final Ratio models ===\n');
[B,FI] = lasso(X_lab_s, y_lab, 'CV', 4);
Models.LASSO.coef = B(:,FI.IndexMinMSE); Models.LASSO.intercept = FI.Intercept(FI.IndexMinMSE);
Models.RF = fitrensemble(X_lab_s, y_lab, 'Method', 'Bag', 'NumLearningCycles', 100, ...
    'Learners', templateTree('MaxNumSplits', 20, 'MinLeafSize', 5, 'Surrogate', 'on'));
Models.GBM = fitrensemble(X_lab_s, y_lab, 'Method', 'LSBoost', 'NumLearningCycles', 100, ...
    'LearnRate', 0.1, 'Learners', templateTree('MaxNumSplits', 10, 'MinLeafSize', 5, 'Surrogate', 'on'));
Models.XGB = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
    'learning_rate',0.1,'min_child_weight',int32(3),'subsample',0.8,'colsample_bytree',0.8,'random_state',int32(42)));
Models.XGB.fit(py.numpy.array(X_lab_s), py.numpy.array(y_lab));
Models.LGB = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
    'learning_rate',0.1,'min_child_samples',int32(5),'subsample',0.8,'colsample_bytree',0.8,...
    'random_state',int32(42),'verbose',int32(-1)));
Models.LGB.fit(py.numpy.array(X_lab_s), py.numpy.array(y_lab));
fprintf('  All models trained.\n');

%% 6. Field Data
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

for k = 1:size(dataFiles,1)
    fpath=dataFiles{k,1}; yr=dataFiles{k,2}; dtype=dataFiles{k,3}; base_date=dataFiles{k,4};
    if ~exist(fpath,'file'), continue; end
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

    chgSegs=local_find_segments(I_cell>thr_A);
    dchSegs=local_find_segments(I_cell<-thr_A);
    if ~isempty(chgSegs), dur=chgSegs(:,2)-chgSegs(:,1)+1; chgSegs=chgSegs(dur>=min_chg_sec(k),:); end
    if ~isempty(dchSegs), dur=dchSegs(:,2)-dchSegs(:,1)+1; dchSegs=dchSegs(dur>=min_dch_sec(k),:); end
    if isempty(chgSegs), fprintf('  No charge segs!\n'); continue; end

    % --- Charge V-Q (raw V, dQ from raw current) ---
    [~,bi]=max(chgSegs(:,2)-chgSegs(:,1));
    cs=chgSegs(bi,1); ce=chgSegs(bi,2);
    v_c=V_avg(cs:ce); i_c=I_cell(cs:ce); t_c=tsec(cs:ce);
    v_c = smoothdata(v_c, 'movmean', 30);
    q_c = cumtrapz(t_c, abs(i_c)) / 3600;
    v_mono=v_c(1); q_mono=q_c(1);
    for ii=2:length(v_c)
        if v_c(ii)>v_mono(end), v_mono(end+1)=v_c(ii); q_mono(end+1)=q_c(ii); end
    end
    dQ_chg = nan(1, num_segs);
    if length(v_mono)>1
        QR = interp1(v_mono, q_mono, V_chg, 'linear', NaN);
        dQ_chg = abs(diff(QR));
    end

    % --- Discharge V-Q (raw V) ---
    dQ_dch = nan(1, num_segs);
    C_eff_dch_val = NaN;
    if ~isempty(dchSegs)
        [~,bi_d]=max(dchSegs(:,2)-dchSegs(:,1));
        ds=dchSegs(bi_d,1); de=dchSegs(bi_d,2);
        v_d=V_avg(ds:de); i_d=I_cell(ds:de); t_d=tsec(ds:de);
        v_d = smoothdata(v_d, 'movmean', 30);
        q_d = cumtrapz(t_d, abs(i_d)) / 3600;
        v_mono_d=v_d(1); q_mono_d=q_d(1);
        for ii=2:length(v_d)
            if v_d(ii)<v_mono_d(end), v_mono_d(end+1)=v_d(ii); q_mono_d(end+1)=q_d(ii); end
        end
        v_asc=flip(v_mono_d); q_asc=flip(q_mono_d);
        if length(v_asc)>1
            QR_d = interp1(v_asc, q_asc, V_dch, 'linear', NaN);
            dQ_dch = abs(diff(QR_d));
        end
        C_eff_dch_val = mean(abs(i_d)) / Q_nom;
    end

    C_eff_chg_val = mean(abs(i_c)) / Q_nom;

    % --- Convert to dQ_ratio (same as Lab) ---
    dQ_c_field = dQ_chg(3:12);  % 10 features
    dQ_d_field = dQ_dch(1:11);  % 11 features
    
    dQ_c_ratio_f = dQ_c_field / sum(dQ_c_field, 'omitnan');
    dQ_d_ratio_f = dQ_d_field / sum(dQ_d_field, 'omitnan');

    X_field = [dQ_c_ratio_f, dQ_d_ratio_f, C_eff_chg_val, C_eff_dch_val];
    nan_mask = isnan(X_field);

    fprintf('  Chg valid segs: %s\n', mat2str(find(~isnan(dQ_chg))));
    fprintf('  Dch valid segs: %s\n', mat2str(find(~isnan(dQ_dch))));
    fprintf('  NaN features: %d / %d\n', sum(nan_mask), length(X_field));
    fprintf('  BMS SOH: %.1f%%\n', soh_bms);
    fprintf('  dQ_c_ratio: %s\n', mat2str(dQ_c_ratio_f, 3));

    % Standardize
    X_field_s = (X_field - mu_lab) ./ sigma_lab;

    % LASSO/RF/GBM: impute NaN → 0 in z-space
    X_imp = X_field_s; X_imp(isnan(X_imp)) = 0;
    soh_lasso = X_imp * Models.LASSO.coef + Models.LASSO.intercept;
    soh_rf = predict(Models.RF, X_imp);
    soh_gbm = predict(Models.GBM, X_imp);

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
    FieldResults.(yr).C_eff_chg = C_eff_chg_val;
    FieldResults.(yr).NaN_cnt  = sum(nan_mask);

    fprintf('  SOH%%: LASSO=%.1f, RF=%.1f, GBM=%.1f, XGB=%.1f, LGB=%.1f\n', ...
        soh_lasso, soh_rf, soh_gbm, soh_xgb, soh_lgb);
end

%% 7. Summary
fprintf('\n\n========== dQ RATIO MODEL - FIELD SOH (SOH %%) ==========\n');
fprintf('%-6s  %-6s  %-7s %-7s %-7s %-7s %-7s  %-3s\n', ...
    'Year','BMS','LASSO','RF','GBM','XGB','LGB','NaN');
yrs = fieldnames(FieldResults);
for k = 1:length(yrs)
    r = FieldResults.(yrs{k});
    fprintf('%-6s  %-6.1f  %-7.1f %-7.1f %-7.1f %-7.1f %-7.1f  %-3d\n', ...
        yrs{k}, r.BMS_SOH, r.LASSO, r.RF, r.GBM, r.XGBoost, r.LightGBM, r.NaN_cnt);
end

%% 8. Save
save(fullfile(rmDir, 'RM_Results.mat'), 'FieldResults', 'Results', 'Models', 'mu_lab', 'sigma_lab', 'feat_names');

%% 9. Visualization
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
title('dQ Ratio Model: Field SOH vs BMS (2021-2025)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 10);
ylim([40 105]);
saveas(fig, fullfile(visDir, 'RM_Field_SOH_Comparison.png'));
fprintf('\nSaved: RM_Field_SOH_Comparison.png\n');
fprintf('=== dQ Ratio Model Complete! ===\n');

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
