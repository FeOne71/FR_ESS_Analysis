% RM_FieldAvailable.m
% Feature Set B: Field-Available Segments Only
% 피처: dQ_c_08, dQ_c_09, dQ_d_08, dQ_d_09, C_eff_chg, C_eff_dch (6개, NaN 없음)
% 모델: LASSO, RF, GBM, XGBoost, LightGBM (5종)
% 검증: LOCO-CV (8-Fold), 필드 4개년 시연

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

% Feature Set B: Seg8, Seg9 (charge) + Seg8, Seg9 (discharge) + C_eff x2
X_lab = [FM.dQ_c_08, FM.dQ_c_09, FM.dQ_d_08, FM.dQ_d_09, FM.C_eff_chg, FM.C_eff_dch];
y_lab = FM.Static_Capacity / Q_nom * 100;
cellIDs = FM.CellID;
feat_names = {'dQ_c08','dQ_c09','dQ_d08','dQ_d09','Ceff_c','Ceff_d'};

fprintf('Feature Set B: %d samples, %d features\n', size(X_lab,1), size(X_lab,2));
fprintf('NaN count: %d (should be 0)\n', sum(isnan(X_lab(:))));
fprintf('SOH range: %.1f ~ %.1f%%\n', min(y_lab), max(y_lab));

%% 3. LOCO-CV
uCells = unique(cellIDs);
nCells = length(uCells);
modelNames = {'LASSO','RF','GBM','XGBoost','LightGBM'};
nModels = length(modelNames);

% Store predictions
CV_pred = nan(size(X_lab,1), nModels);
CV_true = y_lab;

fprintf('\n=== LOCO-CV (Feature Set B, 5 models) ===\n');
for fold = 1:nCells
    tc = uCells(fold);
    trM = ~strcmp(cellIDs, tc);
    teM = strcmp(cellIDs, tc);
    X_tr = X_lab(trM,:); y_tr = y_lab(trM);
    X_te = X_lab(teM,:); y_te = y_lab(teM);
    
    % Standardize (Train stats only)
    mu_tr = mean(X_tr); sig_tr = std(X_tr); sig_tr(sig_tr==0)=1;
    X_tr_s = (X_tr - mu_tr) ./ sig_tr;
    X_te_s = (X_te - mu_tr) ./ sig_tr;
    
    % 1. LASSO
    [B, FitInfo] = lasso(X_tr_s, y_tr, 'Alpha', 1, 'NumLambda', 50, 'CV', 5);
    idxMin = FitInfo.IndexMinMSE;
    b0 = FitInfo.Intercept(idxMin);
    CV_pred(teM, 1) = X_te_s * B(:,idxMin) + b0;
    
    % 2. RF
    mdl_rf = fitrensemble(X_tr_s, y_tr, 'Method', 'Bag', 'NumLearningCycles', 200, ...
        'Learners', templateTree('MaxNumSplits', 20, 'MinLeafSize', 3));
    CV_pred(teM, 2) = predict(mdl_rf, X_te_s);
    
    % 3. GBM
    mdl_gbm = fitrensemble(X_tr_s, y_tr, 'Method', 'LSBoost', 'NumLearningCycles', 200, ...
        'LearnRate', 0.05, 'Learners', templateTree('MaxNumSplits', 10, 'MinLeafSize', 3));
    CV_pred(teM, 3) = predict(mdl_gbm, X_te_s);
    
    % 4. XGBoost
    mdl_xgb = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(4),...
        'learning_rate',0.05,'min_child_weight',int32(3),'subsample',0.8,...
        'colsample_bytree',0.8,'random_state',int32(42)));
    mdl_xgb.fit(py.numpy.array(X_tr_s), py.numpy.array(y_tr));
    CV_pred(teM, 4) = double(mdl_xgb.predict(py.numpy.array(X_te_s)));
    
    % 5. LightGBM
    mdl_lgb = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(4),...
        'learning_rate',0.05,'min_child_samples',int32(3),'subsample',0.8,...
        'colsample_bytree',0.8,'random_state',int32(42),'verbose',int32(-1)));
    mdl_lgb.fit(py.numpy.array(X_tr_s), py.numpy.array(y_tr));
    CV_pred(teM, 5) = double(mdl_lgb.predict(py.numpy.array(X_te_s)));
    
    fprintf('  Fold %d (%s): done\n', fold, tc);
end

%% 4. LOCO-CV Results
fprintf('\n=== LOCO-CV Results ===\n');
fprintf('%-10s  %-6s  %-6s  %-6s\n', 'Model', 'RMSE%', 'MAPE%', 'R2');
CV_Results = struct();
for m = 1:nModels
    err = CV_true - CV_pred(:,m);
    rmse = sqrt(mean(err.^2));
    mape = mean(abs(err./CV_true))*100;
    r2 = 1 - sum(err.^2)/sum((CV_true-mean(CV_true)).^2);
    fprintf('%-10s  %-6.3f  %-6.3f  %-6.4f\n', modelNames{m}, rmse, mape, r2);
    CV_Results.(modelNames{m}) = struct('RMSE',rmse,'MAPE',mape,'R2',r2);
end

%% 5. LOCO-CV Visualization
fig_cv = figure('Position',[50 50 900 400]);
subplot(1,2,1);
rmse_vals = arrayfun(@(m) CV_Results.(modelNames{m}).RMSE, 1:nModels);
bar(rmse_vals, 'FaceColor', [0.3 0.6 0.9]);
set(gca,'XTick',1:nModels,'XTickLabel',modelNames,'FontSize',10,'XTickLabelRotation',20);
ylabel('RMSE (%)'); title('LOCO-CV RMSE'); grid on; ylim([0 max(rmse_vals)*1.3]);
for i=1:nModels, text(i, rmse_vals(i)+0.01, sprintf('%.2f%%',rmse_vals(i)),'HorizontalAlignment','center','FontSize',9); end

subplot(1,2,2);
r2_vals = arrayfun(@(m) CV_Results.(modelNames{m}).R2, 1:nModels);
bar(r2_vals, 'FaceColor', [0.8 0.4 0.2]);
set(gca,'XTick',1:nModels,'XTickLabel',modelNames,'FontSize',10,'XTickLabelRotation',20);
ylabel('R²'); title('LOCO-CV R²'); grid on; ylim([min(0,min(r2_vals))-0.05 1.05]);
for i=1:nModels, text(i, r2_vals(i)+0.01, sprintf('%.3f',r2_vals(i)),'HorizontalAlignment','center','FontSize',9); end
sgtitle('Feature Set B: Field-Available (Seg8,9 + C\_eff)', 'FontSize', 13, 'FontWeight','bold');
saveas(fig_cv, fullfile(visDir, 'RM_FieldAvailable_LOCO_CV.png'));

%% 6. Train Final Models on Full Lab Data
fprintf('\n=== Training final models (full 215 rows) ===\n');
mu_all = mean(X_lab); sig_all = std(X_lab); sig_all(sig_all==0)=1;
X_lab_s = (X_lab - mu_all) ./ sig_all;

[B_f, FI_f] = lasso(X_lab_s, y_lab, 'Alpha',1,'NumLambda',50,'CV',5);
b0_f = FI_f.Intercept(FI_f.IndexMinMSE); Bf = B_f(:,FI_f.IndexMinMSE);
mdl_RF_f = fitrensemble(X_lab_s, y_lab,'Method','Bag','NumLearningCycles',200,...
    'Learners',templateTree('MaxNumSplits',20,'MinLeafSize',3));
mdl_GBM_f = fitrensemble(X_lab_s, y_lab,'Method','LSBoost','NumLearningCycles',200,...
    'LearnRate',0.05,'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',3));
mdl_XGB_f = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(4),...
    'learning_rate',0.05,'min_child_weight',int32(3),'subsample',0.8,'colsample_bytree',0.8,'random_state',int32(42)));
mdl_XGB_f.fit(py.numpy.array(X_lab_s), py.numpy.array(y_lab));
mdl_LGB_f = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(4),...
    'learning_rate',0.05,'min_child_samples',int32(3),'subsample',0.8,'colsample_bytree',0.8,...
    'random_state',int32(42),'verbose',int32(-1)));
mdl_LGB_f.fit(py.numpy.array(X_lab_s), py.numpy.array(y_lab));
fprintf('  All 5 models trained.\n');

%% 7. Field Demo
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
    dQ_c_all=nan(1,12); Ceff_c=NaN;
    if ~isempty(cS)
        [~,bi]=max(cS(:,2)-cS(:,1));cs=cS(bi,1);ce=cS(bi,2);
        vc=V(cs:ce);ic=Ic(cs:ce);tc=ts(cs:ce);
        qc=cumtrapz(tc,abs(ic))/3600;
        vm=vc(1);qm=qc(1);for ii=2:length(vc),if vc(ii)>vm(end),vm(end+1)=vc(ii);qm(end+1)=qc(ii);end;end
        if length(vm)>1,QR=interp1(vm,qm,V_chg,'linear',NaN);dQ_c_all=abs(diff(QR));end
        Ceff_c=mean(abs(ic))/Q_nom;
    end
    % Discharge
    dQ_d_all=nan(1,12); Ceff_d=NaN;
    if ~isempty(dS)
        [~,bi]=max(dS(:,2)-dS(:,1));ds=dS(bi,1);de=dS(bi,2);
        vd=V(ds:de);id=Ic(ds:de);td=ts(ds:de);
        qd=cumtrapz(td,abs(id))/3600;
        vmd=vd(1);qmd=qd(1);for ii=2:length(vd),if vd(ii)<vmd(end),vmd(end+1)=vd(ii);qmd(end+1)=qd(ii);end;end
        va=flip(vmd);qa=flip(qmd);
        if length(va)>1,QRd=interp1(va,qa,V_dch,'linear',NaN);dQ_d_all=abs(diff(QRd));end
        Ceff_d=mean(abs(id))/Q_nom;
    end

    % Feature Set B: Seg8=idx8, Seg9=idx9
    Xf = [dQ_c_all(8), dQ_c_all(9), dQ_d_all(8), dQ_d_all(9), Ceff_c, Ceff_d];
    fprintf('  Xf = [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]  (NaN=%d)\n', ...
        Xf(1),Xf(2),Xf(3),Xf(4),Xf(5),Xf(6), sum(isnan(Xf)));
    
    Xfs = (Xf - mu_all) ./ sig_all;
    Xpy = py.numpy.array(Xfs).reshape(int32(1),int32(-1));
    
    soh_lasso = Xfs * Bf + b0_f;
    soh_rf    = predict(mdl_RF_f, Xfs);
    soh_gbm   = predict(mdl_GBM_f, Xfs);
    soh_xgb   = double(mdl_XGB_f.predict(Xpy));
    soh_lgb   = double(mdl_LGB_f.predict(Xpy));
    
    fprintf('  BMS=%.1f | LASSO=%.1f RF=%.1f GBM=%.1f XGB=%.1f LGB=%.1f\n', ...
        bms, soh_lasso, soh_rf, soh_gbm, soh_xgb, soh_lgb);
    FR.(yr)=struct('BMS',bms,'LASSO',soh_lasso,'RF',soh_rf,'GBM',soh_gbm,'XGB',soh_xgb,'LGB',soh_lgb);
end

%% 8. Summary & Plot
fprintf('\n\n========== FIELD SOH (Feature Set B) ==========\n');
fprintf('%-6s  %-5s  %-7s %-7s %-7s %-7s %-7s\n','Year','BMS','LASSO','RF','GBM','XGB','LGB');
yrs=fieldnames(FR);
for k=1:length(yrs)
    r=FR.(yrs{k});
    fprintf('%-6s  %-5.1f  %-7.1f %-7.1f %-7.1f %-7.1f %-7.1f\n',...
        yrs{k},r.BMS,r.LASSO,r.RF,r.GBM,r.XGB,r.LGB);
end

fig2=figure('Position',[50 100 1100 500]);
n=length(yrs);x=1:n;
bms_a=arrayfun(@(i) FR.(yrs{i}).BMS,  1:n);
las_a=arrayfun(@(i) FR.(yrs{i}).LASSO,1:n);
rf_a =arrayfun(@(i) FR.(yrs{i}).RF,   1:n);
gbm_a=arrayfun(@(i) FR.(yrs{i}).GBM,  1:n);
xgb_a=arrayfun(@(i) FR.(yrs{i}).XGB,  1:n);
lgb_a=arrayfun(@(i) FR.(yrs{i}).LGB,  1:n);
hold on;grid on;box on;
plot(x,bms_a,'k--s','LineWidth',2.5,'MarkerSize',12,'MarkerFaceColor','k','DisplayName','BMS');
plot(x,las_a,'-+','Color',[0.6 0.6 0.6],'LineWidth',1.5,'MarkerSize',10,'DisplayName','LASSO');
plot(x,rf_a, '-o','Color',[0.2 0.7 0.3],'LineWidth',2,'MarkerSize',9,'MarkerFaceColor',[0.2 0.7 0.3],'DisplayName','RF');
plot(x,gbm_a,'-d','Color',[0.9 0.5 0.1],'LineWidth',2,'MarkerSize',9,'MarkerFaceColor',[0.9 0.5 0.1],'DisplayName','GBM');
plot(x,xgb_a,'-^','Color',[0.8 0.2 0.2],'LineWidth',2,'MarkerSize',9,'MarkerFaceColor',[0.8 0.2 0.2],'DisplayName','XGBoost');
plot(x,lgb_a,'-v','Color',[0.3 0.3 0.9],'LineWidth',2,'MarkerSize',9,'MarkerFaceColor',[0.3 0.3 0.9],'DisplayName','LightGBM');
for i=1:n
    text(x(i)+0.05,bms_a(i)+0.4,sprintf('%.1f',bms_a(i)),'FontSize',9,'FontWeight','bold');
end
set(gca,'XTick',x,'XTickLabel',strrep(yrs,'Y',''),'FontSize',12);
xlabel('Year','FontSize',13,'FontWeight','bold');ylabel('SOH (%)','FontSize',13,'FontWeight','bold');
title('Feature Set B (Seg8,9 + C\_eff): 5 Models vs BMS','FontSize',14,'FontWeight','bold');
legend('Location','southwest','FontSize',10);ylim([80 105]);
saveas(fig2,fullfile(visDir,'RM_FieldAvailable_Field_SOH.png'));
fprintf('\nSaved. Done!\n');

function segs=local_find_segments(mask)
    segs=[];n=length(mask);i=1;
    while i<=n,if mask(i),j=i;while j<n&&mask(j+1),j=j+1;end;segs=[segs;i,j];i=j+1;else,i=i+1;end;end
end
