% RM_Final_Pipeline.m
% 최종 통합 모델: 충전+방전 23피처, dQ/sum ratio, 마스킹 증강
% XGBoost/LightGBM only (NaN 네이티브, 0 임퓨테이션 금지)
% LOCO-CV + Feature Importance + Field Demo

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
C_eff_c = FM.C_eff_chg; C_eff_d = FM.C_eff_dch;
y_lab = FM.Static_Capacity / Q_nom * 100;
cellIDs = FM.CellID;
nRows = size(FM,1); nC = 10; nD = 11;

feat_names = [arrayfun(@(i) sprintf('dQr_c%02d',i), 3:12, 'UniformOutput', false), ...
              arrayfun(@(i) sprintf('dQr_d%02d',i), 1:11, 'UniformOutput', false), ...
              {'Ceff_c', 'Ceff_d'}];

%% 3. Contiguous masking patterns
chg_pat = {}; for w=1:nC, for s=1:(nC-w+1), chg_pat{end+1}=s:(s+w-1); end; end
dch_pat = {}; for w=1:nD, for s=1:(nD-w+1), dch_pat{end+1}=s:(s+w-1); end; end
fprintf('Patterns: charge=%d, discharge=%d\n', length(chg_pat), length(dch_pat));

%% 4. Build Combined Augmented Data
rng(42); N_aug = 30;
X_aug = []; y_aug = []; cid_aug = string([]);

for i = 1:nRows
    dQc = dQ_c_raw(i,:); dQd = dQ_d_raw(i,:);
    
    % Original (no masking): ratio with all segments
    rc = dQc / sum(dQc, 'omitnan');
    rd = dQd / sum(dQd, 'omitnan');
    X_aug(end+1,:) = [rc, rd, C_eff_c(i), C_eff_d(i)];
    y_aug(end+1,1) = y_lab(i);
    cid_aug(end+1,1) = string(cellIDs(i));
    
    % N_aug masked versions (charge & discharge masked independently)
    for a = 1:N_aug
        % Charge masking
        cp = chg_pat{randi(length(chg_pat))};
        mc = nan(1, nC); mc(cp) = dQc(cp);
        rc_m = mc / sum(mc, 'omitnan');
        
        % Discharge masking
        dp = dch_pat{randi(length(dch_pat))};
        md = nan(1, nD); md(dp) = dQd(dp);
        rd_m = md / sum(md, 'omitnan');
        
        X_aug(end+1,:) = [rc_m, rd_m, C_eff_c(i), C_eff_d(i)];
        y_aug(end+1,1) = y_lab(i);
        cid_aug(end+1,1) = string(cellIDs(i));
    end
end

fprintf('Augmented: %d rows, %d features, %.1f%% NaN\n', ...
    size(X_aug,1), size(X_aug,2), sum(isnan(X_aug(:)))/numel(X_aug)*100);

%% 5. Standardize (NaN stays NaN)
mu_lab = mean(X_aug, 1, 'omitnan');
sig_lab = std(X_aug, 0, 1, 'omitnan'); sig_lab(sig_lab==0)=1;
X_aug_s = (X_aug - mu_lab) ./ sig_lab;

%% 6. LOCO-CV
uCells = unique(cellIDs); nCells = length(uCells);
fprintf('\n=== LOCO-CV (combined 23 features, masking augmented) ===\n');

for m = 1:2
    if m==1, mName='XGBoost'; else, mName='LightGBM'; end
    y_true=[]; y_pred=[];
    
    for fold = 1:nCells
        tc = uCells(fold);
        trM = ~strcmp(cid_aug, tc);
        X_tr = X_aug_s(trM,:); y_tr = y_aug(trM);
        
        % Test on ORIGINAL (unmasked) data for this cell
        teM = strcmp(cellIDs, tc);
        dQc_te = dQ_c_raw(teM,:); dQd_te = dQ_d_raw(teM,:);
        rc_te = dQc_te ./ sum(dQc_te, 2, 'omitnan');
        rd_te = dQd_te ./ sum(dQd_te, 2, 'omitnan');
        X_te = [rc_te, rd_te, C_eff_c(teM), C_eff_d(teM)];
        X_te_s = (X_te - mu_lab) ./ sig_lab;
        y_te = y_lab(teM);
        
        if m==1
            mdl = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
                'learning_rate',0.05,'min_child_weight',int32(3),'subsample',0.8,...
                'colsample_bytree',0.8,'random_state',int32(42)));
        else
            mdl = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
                'learning_rate',0.05,'min_child_samples',int32(5),'subsample',0.8,...
                'colsample_bytree',0.8,'random_state',int32(42),'verbose',int32(-1)));
        end
        mdl.fit(py.numpy.array(X_tr), py.numpy.array(y_tr));
        yp = double(mdl.predict(py.numpy.array(X_te_s)));
        y_true=[y_true;y_te]; y_pred=[y_pred;yp(:)];
    end
    
    err=y_true-y_pred; rmse=sqrt(mean(err.^2));
    mape=mean(abs(err./y_true))*100;
    r2=1-sum(err.^2)/sum((y_true-mean(y_true)).^2);
    fprintf('  %s: RMSE=%.3f%%, MAPE=%.3f%%, R2=%.4f\n', mName, rmse, mape, r2);
end

%% 7. Train final models on ALL augmented data
fprintf('\n=== Training final models ===\n');
XGB = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
    'learning_rate',0.05,'min_child_weight',int32(3),'subsample',0.8,...
    'colsample_bytree',0.8,'random_state',int32(42)));
XGB.fit(py.numpy.array(X_aug_s), py.numpy.array(y_aug));

LGB = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
    'learning_rate',0.05,'min_child_samples',int32(5),'subsample',0.8,...
    'colsample_bytree',0.8,'random_state',int32(42),'verbose',int32(-1)));
LGB.fit(py.numpy.array(X_aug_s), py.numpy.array(y_aug));
fprintf('  XGBoost & LightGBM trained.\n');

%% 8. Feature Importance
fprintf('\n=== Feature Importance ===\n');
fi_xgb = double(XGB.feature_importances_);
fi_lgb = double(LGB.feature_importances_);
fi_lgb_norm = fi_lgb / sum(fi_lgb);  % normalize LGB

fig_fi = figure('Position', [50 50 1200 400]);
subplot(1,2,1);
bar(fi_xgb, 'FaceColor', [0.8 0.3 0.1]); hold on;
set(gca, 'XTick', 1:23, 'XTickLabel', feat_names, 'XTickLabelRotation', 45, 'FontSize', 8);
title('XGBoost Feature Importance', 'FontSize', 12);
ylabel('Importance'); grid on;

subplot(1,2,2);
bar(fi_lgb_norm, 'FaceColor', [0.3 0.6 0.9]); hold on;
set(gca, 'XTick', 1:23, 'XTickLabel', feat_names, 'XTickLabelRotation', 45, 'FontSize', 8);
title('LightGBM Feature Importance', 'FontSize', 12);
ylabel('Importance (normalized)'); grid on;

sgtitle('Feature Importance: dQ Ratio + Masking Augmented', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig_fi, fullfile(visDir, 'RM_Feature_Importance.png'));
fprintf('  Saved: RM_Feature_Importance.png\n');

% Print top features
[~, si] = sort(fi_xgb, 'descend');
fprintf('  XGB Top 5: '); for j=1:5, fprintf('%s(%.3f) ', feat_names{si(j)}, fi_xgb(si(j))); end; fprintf('\n');
[~, si] = sort(fi_lgb_norm, 'descend');
fprintf('  LGB Top 5: '); for j=1:5, fprintf('%s(%.3f) ', feat_names{si(j)}, fi_lgb_norm(si(j))); end; fprintf('\n');

%% 9. Field Demo
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
    fprintf('\n=== %s ===\n', yr);
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

    % Ratio
    rc = dQ_c_f / sum(dQ_c_f, 'omitnan');
    rd = dQ_d_f / sum(dQ_d_f, 'omitnan');
    Xf = [rc, rd, Ceff_c, Ceff_d];
    Xfs = (Xf - mu_lab) ./ sig_lab;

    Xp = py.numpy.array(Xfs).reshape(int32(1),int32(-1));
    soh_xgb = double(XGB.predict(Xp));
    soh_lgb = double(LGB.predict(Xp));

    chg_segs = find(~isnan(dQ_c_f)); dch_segs = find(~isnan(dQ_d_f));
    nan_cnt = sum(isnan(Xf));
    fprintf('  CHG segs: %s, DCH segs: %s, NaN: %d/23\n', mat2str(chg_segs), mat2str(dch_segs), nan_cnt);
    fprintf('  BMS=%.1f%%, XGB=%.1f%%, LGB=%.1f%%\n', bms, soh_xgb, soh_lgb);

    FR.(yr) = struct('BMS',bms,'XGB',soh_xgb,'LGB',soh_lgb,'NaN',nan_cnt);
end

%% 10. Summary
fprintf('\n\n========== FINAL MODEL: FIELD SOH (SOH%%) ==========\n');
fprintf('%-6s  %-5s  %-8s %-8s %-4s\n','Year','BMS','XGBoost','LightGBM','NaN');
yrs=fieldnames(FR);
for k=1:length(yrs)
    r=FR.(yrs{k});
    fprintf('%-6s  %-5.1f  %-8.1f %-8.1f %-4d\n',yrs{k},r.BMS,r.XGB,r.LGB,r.NaN);
end

%% 11. Visualization
fig=figure('Position',[50 100 1000 500]);
n=length(yrs);x=1:n;
bms_a=arrayfun(@(i) FR.(yrs{i}).BMS,1:n);
xgb_a=arrayfun(@(i) FR.(yrs{i}).XGB,1:n);
lgb_a=arrayfun(@(i) FR.(yrs{i}).LGB,1:n);
hold on;grid on;box on;
plot(x,bms_a,'k--s','LineWidth',2.5,'MarkerSize',12,'MarkerFaceColor','k','DisplayName','BMS SOH');
plot(x,xgb_a,'-o','Color',[0.8 0.2 0.2],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor',[0.8 0.2 0.2],'DisplayName','XGBoost');
plot(x,lgb_a,'-d','Color',[0.2 0.5 0.9],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor',[0.2 0.5 0.9],'DisplayName','LightGBM');
for i=1:n
    text(x(i)+0.05,bms_a(i)+0.4,sprintf('%.1f%%',bms_a(i)),'FontSize',9,'FontWeight','bold');
    text(x(i)+0.05,xgb_a(i)-0.7,sprintf('%.1f%%',xgb_a(i)),'FontSize',8,'Color',[0.8 0.2 0.2]);
    text(x(i)+0.05,lgb_a(i)+0.4,sprintf('%.1f%%',lgb_a(i)),'FontSize',8,'Color',[0.2 0.5 0.9]);
end
set(gca,'XTick',x,'XTickLabel',strrep(yrs,'Y',''),'FontSize',12);
xlabel('Year','FontSize',13,'FontWeight','bold');ylabel('SOH (%)','FontSize',13,'FontWeight','bold');
title('Final Model: Masking-Augmented dQ Ratio (Combined 23 Features)','FontSize',14,'FontWeight','bold');
legend('Location','southwest','FontSize',11);ylim([80 105]);
saveas(fig,fullfile(visDir,'RM_Final_Field_SOH.png'));
fprintf('\nSaved: RM_Final_Field_SOH.png\n');
fprintf('=== Final Pipeline Complete! ===\n');

%% Helper
function segs=local_find_segments(mask)
    segs=[];n=length(mask);i=1;
    while i<=n,if mask(i),j=i;while j<n&&mask(j+1),j=j+1;end;segs=[segs;i,j];i=j+1;else,i=i+1;end;end
end
