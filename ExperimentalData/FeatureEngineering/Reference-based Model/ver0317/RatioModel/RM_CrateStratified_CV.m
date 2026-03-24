% RM_CrateStratified_CV.m
% dQ ratio 피처 + C-rate별 LOCO-CV 검증
% 핵심: 같은 C-rate 내에서만 CV → C-rate 혼재 왜곡 제거

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

dQ_c_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i), 3:12, 'UniformOutput', false);
dQ_d_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false);
dQ_c_raw = FM{:, dQ_c_cols};  % 215 x 10
dQ_d_raw = FM{:, dQ_d_cols};  % 215 x 11
C_eff_c = FM.C_eff_chg;
C_eff_d = FM.C_eff_dch;
y_lab = FM.Static_Capacity / Q_nom * 100;
cellIDs = FM.CellID;
conditions = FM.Condition;  % C-rate 정보

% dQ ratio 변환 (원본, 마스킹 없음)
dQr_c = dQ_c_raw ./ sum(dQ_c_raw, 2, 'omitnan');  % 215 x 10
dQr_d = dQ_d_raw ./ sum(dQ_d_raw, 2, 'omitnan');  % 215 x 11
X_ratio = [dQr_c, dQr_d, C_eff_c, C_eff_d];  % 215 x 23

% 고유 C-rate 목록
uConds = unique(conditions);
fprintf('C-rate conditions: ');
disp(uConds);
fprintf('Total: %d rows, NaN: %d\n', size(X_ratio,1), sum(isnan(X_ratio(:))));

%% 3. C-rate별 LOCO-CV
fprintf('\n=== C-rate Stratified LOCO-CV (dQ Ratio) ===\n');
fprintf('%-15s  %-5s  %-6s  %-6s  %-6s\n', 'Condition', 'N', 'RMSE%', 'MAPE%', 'R2');

all_results = struct();
for ci = 1:length(uConds)
    cond = uConds(ci);
    condMask = strcmp(conditions, cond);
    X_c = X_ratio(condMask,:);
    y_c = y_lab(condMask);
    ids_c = cellIDs(condMask);
    
    uCells = unique(ids_c);
    nCells = length(uCells);
    
    % Standardize on full condition data
    mu_c = mean(X_c, 1, 'omitnan');
    sig_c = std(X_c, 0, 1, 'omitnan'); sig_c(sig_c==0)=1;
    X_c_s = (X_c - mu_c) ./ sig_c;
    
    y_pred_all = nan(size(y_c));
    
    for fold = 1:nCells
        tc = uCells(fold);
        trM = ~strcmp(ids_c, tc);
        teM = strcmp(ids_c, tc);
        X_tr = X_c_s(trM,:); y_tr = y_c(trM);
        X_te = X_c_s(teM,:);
        
        % GBM (handles NaN via surrogate split)
        mdl = fitrensemble(X_tr, y_tr, 'Method', 'LSBoost', ...
            'NumLearningCycles', 200, 'LearnRate', 0.05, ...
            'Learners', templateTree('MaxNumSplits', 10, 'MinLeafSize', 3, 'Surrogate', 'on'));
        y_pred_all(teM) = predict(mdl, X_te);
    end
    
    err = y_c - y_pred_all;
    rmse = sqrt(mean(err.^2, 'omitnan'));
    mape = mean(abs(err./y_c), 'omitnan')*100;
    r2 = 1 - sum(err.^2,'omitnan')/sum((y_c-mean(y_c)).^2,'omitnan');
    fprintf('%-15s  %-5d  %-6.3f  %-6.3f  %-6.4f\n', char(cond), sum(condMask), rmse, mape, r2);
    
    cname = matlab.lang.makeValidName(char(cond));
    all_results.(cname) = struct('Condition', cond, 'N', sum(condMask), 'RMSE', rmse, 'MAPE', mape, 'R2', r2);
end

%% 4. Mixed C-rate LOCO-CV (비교용: 기존 방식)
fprintf('\n--- Mixed C-rate LOCO-CV (기존 비교) ---\n');
uCells_all = unique(cellIDs);
mu_all = mean(X_ratio, 1, 'omitnan'); sig_all = std(X_ratio, 0, 1, 'omitnan'); sig_all(sig_all==0)=1;
X_all_s = (X_ratio - mu_all) ./ sig_all;
y_pred_mix = nan(size(y_lab));

for fold = 1:length(uCells_all)
    tc = uCells_all(fold);
    trM = ~strcmp(cellIDs, tc); teM = strcmp(cellIDs, tc);
    mdl = fitrensemble(X_all_s(trM,:), y_lab(trM), 'Method', 'LSBoost', ...
        'NumLearningCycles', 200, 'LearnRate', 0.05, ...
        'Learners', templateTree('MaxNumSplits', 10, 'MinLeafSize', 3, 'Surrogate', 'on'));
    y_pred_mix(teM) = predict(mdl, X_all_s(teM,:));
end
err_mix = y_lab - y_pred_mix;
rmse_mix = sqrt(mean(err_mix.^2)); mape_mix = mean(abs(err_mix./y_lab))*100;
r2_mix = 1-sum(err_mix.^2)/sum((y_lab-mean(y_lab)).^2);
fprintf('Mixed LOCO-CV: RMSE=%.3f%%, MAPE=%.3f%%, R2=%.4f\n', rmse_mix, mape_mix, r2_mix);

%% 5. Visualization: C-rate CV summary
fig1 = figure('Position', [50 50 1000 400]);
conds = fieldnames(all_results);
rmse_vals = arrayfun(@(i) all_results.(conds{i}).RMSE, 1:length(conds));
r2_vals   = arrayfun(@(i) all_results.(conds{i}).R2,   1:length(conds));
cond_labels = arrayfun(@(i) char(all_results.(conds{i}).Condition), 1:length(conds), 'UniformOutput', false);

subplot(1,2,1);
bar([rmse_vals, rmse_mix], 'FaceColor', [0.3 0.6 0.9]);
set(gca,'XTick',1:length(conds)+1,'XTickLabel',[cond_labels,{'Mixed'}],'FontSize',10,'XTickLabelRotation',30);
ylabel('RMSE (%)'); title('LOCO-CV RMSE by C-rate'); grid on;
for i=1:length(conds)+1
    vals = [rmse_vals, rmse_mix];
    text(i, vals(i)+0.02, sprintf('%.3f',vals(i)),'HorizontalAlignment','center','FontSize',8);
end

subplot(1,2,2);
bar([r2_vals, r2_mix], 'FaceColor', [0.8 0.4 0.2]);
set(gca,'XTick',1:length(conds)+1,'XTickLabel',[cond_labels,{'Mixed'}],'FontSize',10,'XTickLabelRotation',30);
ylabel('R²'); title('LOCO-CV R² by C-rate'); grid on; ylim([-0.2 1.1]);
for i=1:length(conds)+1
    vals = [r2_vals, r2_mix];
    text(i, vals(i)+0.02, sprintf('%.3f',vals(i)),'HorizontalAlignment','center','FontSize',8);
end
sgtitle('dQ Ratio: C-rate Stratified vs Mixed LOCO-CV', 'FontSize', 13, 'FontWeight', 'bold');
saveas(fig1, fullfile(visDir, 'RM_Crate_Stratified_CV.png'));
fprintf('Saved: RM_Crate_Stratified_CV.png\n');

%% 6. Train full model & Field Demo
fprintf('\n=== Training full model (all C-rates, GBM) ===\n');
mdl_full = fitrensemble(X_all_s, y_lab, 'Method', 'LSBoost', ...
    'NumLearningCycles', 200, 'LearnRate', 0.05, ...
    'Learners', templateTree('MaxNumSplits', 10, 'MinLeafSize', 3, 'Surrogate', 'on'));

% XGBoost (NaN native)
mdl_xgb = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
    'learning_rate',0.05,'min_child_weight',int32(3),'subsample',0.8,...
    'colsample_bytree',0.8,'random_state',int32(42)));
mdl_xgb.fit(py.numpy.array(X_all_s), py.numpy.array(y_lab));
fprintf('  GBM + XGBoost trained.\n');

%% 7. Field
dataFiles = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021', 'old', datetime(2021,6,3);
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023', 'new', datetime(2023,10,16);
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024', 'new', datetime(2024,9,9);
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'), 'Y2025', 'new', datetime(2025,7,11);
};
Np=2; thr_A=Q_nom*0.05/Np;
min_chg_sec=[600,300,300,300]; min_dch_sec=[300,150,300,150];
FR=struct();

for k=1:size(dataFiles,1)
    fpath=dataFiles{k,1}; yr=dataFiles{k,2}; dtype=dataFiles{k,3}; bd=dataFiles{k,4};
    if ~exist(fpath,'file'), continue; end
    S=load(fpath);
    if strcmp(dtype,'old')
        D=S.Raw.Rack01;
        if isfield(D,'Time'),t=datetime(D.Time);elseif isduration(D.Date_Time),t=bd+D.Date_Time;else,t=datetime(D.Date_Time);end
        I_r=D.DCCurrent_A(:);V=D.AverageCV_V(:);
        if isfield(D,'SOHPct'),rs=D.SOHPct(:);else,rs=nan(size(I_r));end
    else
        D=S.Raw;
        if isduration(D.Date_Time),t=bd+D.Date_Time;else,t=datetime(D.Date_Time);end
        I_r=D.DCCurrent(:);V=D.CVavg(:);
        if isfield(D,'SOH_BMS'),rs=D.SOH_BMS(:);else,rs=nan(size(I_r));end
    end
    ts=seconds(t-t(1));Ic=I_r/Np;vs=rs;vs(vs<=0)=NaN;bms=median(vs,'omitnan');
    cS=local_find_segments(Ic>thr_A);dS=local_find_segments(Ic<-thr_A);
    if ~isempty(cS),dur=cS(:,2)-cS(:,1)+1;cS=cS(dur>=min_chg_sec(k),:);end
    if ~isempty(dS),dur=dS(:,2)-dS(:,1)+1;dS=dS(dur>=min_dch_sec(k),:);end

    % Charge
    dQ_c_f=nan(1,10); Ceff_c=NaN;
    if ~isempty(cS)
        [~,bi]=max(cS(:,2)-cS(:,1));cs=cS(bi,1);ce=cS(bi,2);
        vc=V(cs:ce);ic=Ic(cs:ce);tc=ts(cs:ce);
        qc=cumtrapz(tc,abs(ic))/3600;
        vm=vc(1);qm=qc(1);for ii=2:length(vc),if vc(ii)>vm(end),vm(end+1)=vc(ii);qm(end+1)=qc(ii);end;end
        dQc=nan(1,12);if length(vm)>1,QR=interp1(vm,qm,V_chg,'linear',NaN);dQc=abs(diff(QR));end
        dQ_c_f=dQc(3:12); Ceff_c=mean(abs(ic))/Q_nom;
    end
    % Discharge
    dQ_d_f=nan(1,11); Ceff_d=NaN;
    if ~isempty(dS)
        [~,bi]=max(dS(:,2)-dS(:,1));ds=dS(bi,1);de=dS(bi,2);
        vd=V(ds:de);id=Ic(ds:de);td=ts(ds:de);
        qd=cumtrapz(td,abs(id))/3600;
        vmd=vd(1);qmd=qd(1);for ii=2:length(vd),if vd(ii)<vmd(end),vmd(end+1)=vd(ii);qmd(end+1)=qd(ii);end;end
        va=flip(vmd);qa=flip(qmd);
        dQd=nan(1,12);if length(va)>1,QRd=interp1(va,qa,V_dch,'linear',NaN);dQd=abs(diff(QRd));end
        dQ_d_f=dQd(1:11); Ceff_d=mean(abs(id))/Q_nom;
    end

    % ratio
    rc = dQ_c_f / sum(dQ_c_f,'omitnan');
    rd = dQ_d_f / sum(dQ_d_f,'omitnan');
    Xf = [rc, rd, Ceff_c, Ceff_d];
    Xfs = (Xf - mu_all) ./ sig_all;
    
    soh_gbm = predict(mdl_full, Xfs);
    Xpy = py.numpy.array(Xfs).reshape(int32(1),int32(-1));
    soh_xgb = double(mdl_xgb.predict(Xpy));
    
    nan_cnt = sum(isnan(Xf));
    fprintf('%s: BMS=%.1f%%, GBM=%.1f%%, XGB=%.1f%%, NaN=%d/23, Ceff=[%.2f,%.2f]\n',...
        yr, bms, soh_gbm, soh_xgb, nan_cnt, Ceff_c, Ceff_d);
    FR.(yr) = struct('BMS',bms,'GBM',soh_gbm,'XGB',soh_xgb,'NaN',nan_cnt,'Ceff_c',Ceff_c,'Ceff_d',Ceff_d);
end

fprintf('\n========== SUMMARY ==========\n');
fprintf('%-6s %-5s %-7s %-7s %-4s %-10s\n','Year','BMS','GBM','XGB','NaN','Ceff(c/d)');
yrs=fieldnames(FR);
for k=1:length(yrs)
    r=FR.(yrs{k});
    fprintf('%-6s %-5.1f %-7.1f %-7.1f %-4d %.2f/%.2f\n',...
        yrs{k},r.BMS,r.GBM,r.XGB,r.NaN,r.Ceff_c,r.Ceff_d);
end
fprintf('=== Done! ===\n');

function segs=local_find_segments(mask)
    segs=[];n=length(mask);i=1;
    while i<=n,if mask(i),j=i;while j<n&&mask(j+1),j=j+1;end;segs=[segs;i,j];i=j+1;else,i=i+1;end;end
end
