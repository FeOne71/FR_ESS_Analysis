% RM_RelativeDQ.m
% Relative dQ: 분모를 항상 존재하는 Seg8+Seg9 합으로 고정
% Lab과 Field에서 동일한 피처 스케일 보장
% 모델: GBM (Surrogate Split, NaN 강건) + XGBoost (비교용)
% 검증: C-rate별 LOCO-CV + 필드 4개년

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
d  = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat'));  FM = d.FM;
mr = load(fullfile(verDir, 'MasterRulers_PseudoOCV.mat'));
V_chg = mr.MasterRuler_ver0317.V_bounds_chg;
V_dch = mr.MasterRuler_ver0317.V_bounds_dch;
Q_nom = 64;

% Seg 인덱스 (dQ_c_raw col 1=Seg03, col 8=Seg10 / dQ_d_raw col 1=Seg01, col 8=Seg08)
% dQ_c_cols: Seg03~12 (10개), dQ_d_cols: Seg01~11 (11개)
dQ_c_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i), 3:12, 'UniformOutput', false);
dQ_d_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false);
dQ_c_raw = FM{:, dQ_c_cols};  % 215×10: col6=Seg08, col7=Seg09
dQ_d_raw = FM{:, dQ_d_cols};  % 215×11: col8=Seg08, col9=Seg09

C_eff_c = FM.C_eff_chg; C_eff_d = FM.C_eff_dch;
y_lab   = FM.Static_Capacity / Q_nom * 100;
cellIDs = FM.CellID;
conditions = FM.Condition;

%% 3. Compute Relative dQ (fixed denominator = Seg8+Seg9)
% Charge: col6=Seg08, col7=Seg09
denom_c = dQ_c_raw(:,6) + dQ_c_raw(:,7);   % always available in Lab
rdQ_c   = dQ_c_raw ./ denom_c;              % 215×10, each col / (Seg08+Seg09)

% Discharge: col8=Seg08, col9=Seg09
denom_d = dQ_d_raw(:,8) + dQ_d_raw(:,9);
rdQ_d   = dQ_d_raw ./ denom_d;              % 215×11

X_lab = [rdQ_c, rdQ_d, C_eff_c, C_eff_d];  % 215×23
feat_names_c = arrayfun(@(i) sprintf('rdQ_c%02d',i), 3:12, 'UniformOutput',false);
feat_names_d = arrayfun(@(i) sprintf('rdQ_d%02d',i), 1:11, 'UniformOutput',false);
feat_names   = [feat_names_c, feat_names_d, {'Ceff_c','Ceff_d'}];

fprintf('RelativeDQ features: %d total, NaN=%d (Lab Seg03-12,01-11)\n', ...
    size(X_lab,2), sum(isnan(X_lab(:))));
fprintf('SOH range: %.1f~%.1f%%\n', min(y_lab), max(y_lab));

%% 4. C-rate Stratified LOCO-CV
uConds = unique(conditions);
fprintf('\n=== C-rate Stratified LOCO-CV (Relative dQ, GBM Surrogate) ===\n');
fprintf('%-8s  %-5s  %-7s %-7s\n','Cond','N','RMSE%','R2');

mu_all  = mean(X_lab,1,'omitnan');
sig_all = std(X_lab,0,1,'omitnan'); sig_all(sig_all==0)=1;
X_lab_s = (X_lab - mu_all) ./ sig_all;

CV_pred_gbm = nan(size(y_lab));
CV_pred_xgb = nan(size(y_lab));

for ci = 1:length(uConds)
    cond = uConds(ci);
    cM = strcmp(conditions, cond);
    X_c = X_lab_s(cM,:); y_c = y_lab(cM); ids_c = cellIDs(cM);
    uCells = unique(ids_c);

    y_gbm_c = nan(size(y_c));
    y_xgb_c = nan(size(y_c));

    for fold = 1:length(uCells)
        tc = uCells(fold);
        trM = ~strcmp(ids_c,tc); teM = strcmp(ids_c,tc);
        X_tr=X_c(trM,:); y_tr=y_c(trM); X_te=X_c(teM,:);

        % GBM Surrogate
        mdl_g = fitrensemble(X_tr, y_tr, 'Method','LSBoost', ...
            'NumLearningCycles',200,'LearnRate',0.05, ...
            'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',3,'Surrogate','on'));
        y_gbm_c(teM) = predict(mdl_g, X_te);

        % XGBoost (NaN native)
        mdl_x = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
            'learning_rate',0.05,'min_child_weight',int32(3),'subsample',0.8,...
            'colsample_bytree',0.8,'random_state',int32(42)));
        mdl_x.fit(py.numpy.array(X_tr), py.numpy.array(y_tr));
        y_xgb_c(teM) = double(mdl_x.predict(py.numpy.array(X_te)));
    end

    r2g  = 1-sum((y_c-y_gbm_c).^2)/sum((y_c-mean(y_c)).^2);
    rmseg= sqrt(mean((y_c-y_gbm_c).^2));
    r2x  = 1-sum((y_c-y_xgb_c).^2)/sum((y_c-mean(y_c)).^2);
    rmsex= sqrt(mean((y_c-y_xgb_c).^2));
    fprintf('%-8s  %-5d  GBM: %-6.3f %-6.4f  |  XGB: %-6.3f %-6.4f\n', ...
        char(cond), sum(cM), rmseg,r2g, rmsex,r2x);

    CV_pred_gbm(cM) = y_gbm_c;
    CV_pred_xgb(cM) = y_xgb_c;
end

% Overall
e_g=y_lab-CV_pred_gbm; e_x=y_lab-CV_pred_xgb;
fprintf('\nOverall GBM: RMSE=%.3f%%, R2=%.4f\n', sqrt(mean(e_g.^2)), ...
    1-sum(e_g.^2)/sum((y_lab-mean(y_lab)).^2));
fprintf('Overall XGB: RMSE=%.3f%%, R2=%.4f\n', sqrt(mean(e_x.^2)), ...
    1-sum(e_x.^2)/sum((y_lab-mean(y_lab)).^2));

%% 5. Train final models (all data)
fprintf('\n=== Training final models ===\n');
mdl_GBM = fitrensemble(X_lab_s, y_lab, 'Method','LSBoost', ...
    'NumLearningCycles',200,'LearnRate',0.05, ...
    'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',3,'Surrogate','on'));
mdl_XGB = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
    'learning_rate',0.05,'min_child_weight',int32(3),'subsample',0.8,...
    'colsample_bytree',0.8,'random_state',int32(42)));
mdl_XGB.fit(py.numpy.array(X_lab_s), py.numpy.array(y_lab));
fprintf('  GBM + XGBoost trained.\n');

%% 6. Field Demo
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
    dQ_c12=nan(1,12); Ceff_c=NaN;
    if ~isempty(cS)
        [~,bi]=max(cS(:,2)-cS(:,1));cs=cS(bi,1);ce=cS(bi,2);
        vc=V(cs:ce);ic=Ic(cs:ce);tc=ts(cs:ce);
        qc=cumtrapz(tc,abs(ic))/3600;
        vm=vc(1);qm=qc(1);for ii=2:length(vc),if vc(ii)>vm(end),vm(end+1)=vc(ii);qm(end+1)=qc(ii);end;end
        if length(vm)>1,QR=interp1(vm,qm,V_chg,'linear',NaN);dQ_c12=abs(diff(QR));end
        Ceff_c=mean(abs(ic))/Q_nom;
    end
    % Discharge
    dQ_d12=nan(1,12); Ceff_d=NaN;
    if ~isempty(dS)
        [~,bi]=max(dS(:,2)-dS(:,1));ds=dS(bi,1);de=dS(bi,2);
        vd=V(ds:de);id=Ic(ds:de);td=ts(ds:de);
        qd=cumtrapz(td,abs(id))/3600;
        vmd=vd(1);qmd=qd(1);for ii=2:length(vd),if vd(ii)<vmd(end),vmd(end+1)=vd(ii);qmd(end+1)=qd(ii);end;end
        va=flip(vmd);qa=flip(qmd);
        if length(va)>1,QRd=interp1(va,qa,V_dch,'linear',NaN);dQ_d12=abs(diff(QRd));end
        Ceff_d=mean(abs(id))/Q_nom;
    end

    % Field relative dQ (same formula as Lab)
    % Charge: Seg03~12 = dQ_c12(3:12), Seg08=col8, Seg09=col9
    dQ_c_f = dQ_c12(3:12);   % 10개 (Seg03~12)
    denom_cf = dQ_c_f(6) + dQ_c_f(7);   % Seg08+Seg09
    rdQ_c_f  = dQ_c_f / denom_cf;

    % Discharge: Seg01~11 = dQ_d12(1:11), Seg08=col8, Seg09=col9
    dQ_d_f = dQ_d12(1:11);   % 11개 (Seg01~11)
    denom_df = dQ_d_f(8) + dQ_d_f(9);   % Seg08+Seg09
    rdQ_d_f  = dQ_d_f / denom_df;

    Xf  = [rdQ_c_f, rdQ_d_f, Ceff_c, Ceff_d];
    Xfs = (Xf - mu_all) ./ sig_all;

    soh_gbm = predict(mdl_GBM, Xfs);
    Xpy = py.numpy.array(Xfs).reshape(int32(1),int32(-1));
    soh_xgb = double(mdl_XGB.predict(Xpy));

    nan_c = sum(isnan(rdQ_c_f)); nan_d = sum(isnan(rdQ_d_f));
    fprintf('%s: BMS=%.1f%%  GBM=%.1f%%  XGB=%.1f%%  NaN(c/d)=%d/%d  Ceff=[%.2f,%.2f]\n',...
        yr,bms,soh_gbm,soh_xgb,nan_c,nan_d,Ceff_c,Ceff_d);
    FR.(yr)=struct('BMS',bms,'GBM',soh_gbm,'XGB',soh_xgb,'NaNc',nan_c,'NaNd',nan_d,'Ceff_c',Ceff_c,'Ceff_d',Ceff_d);
end

%% 7. Summary + Plot
fprintf('\n========== SUMMARY ==========\n');
fprintf('%-6s %-5s %-7s %-7s\n','Year','BMS','GBM','XGB');
yrs=fieldnames(FR);
for k=1:length(yrs)
    r=FR.(yrs{k});
    fprintf('%-6s %-5.1f %-7.1f %-7.1f\n',yrs{k},r.BMS,r.GBM,r.XGB);
end

fig=figure('Position',[50 100 900 450]);
n=length(yrs);x=1:n;
bv=arrayfun(@(i)FR.(yrs{i}).BMS,1:n); gv=arrayfun(@(i)FR.(yrs{i}).GBM,1:n);
xv=arrayfun(@(i)FR.(yrs{i}).XGB,1:n);
hold on;grid on;box on;
plot(x,bv,'k--s','LineWidth',2.5,'MarkerSize',12,'MarkerFaceColor','k','DisplayName','BMS SOH');
plot(x,gv,'-o','Color',[0.2 0.7 0.3],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor',[0.2 0.7 0.3],'DisplayName','GBM (Surrogate)');
plot(x,xv,'-^','Color',[0.8 0.2 0.2],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor',[0.8 0.2 0.2],'DisplayName','XGBoost');
for i=1:n
    text(x(i)+0.05,bv(i)+0.4,sprintf('%.1f%%',bv(i)),'FontSize',9,'FontWeight','bold');
    text(x(i)+0.05,gv(i)-0.7,sprintf('%.1f%%',gv(i)),'FontSize',8,'Color',[0.2 0.7 0.3]);
    text(x(i)+0.05,xv(i)+0.4,sprintf('%.1f%%',xv(i)),'FontSize',8,'Color',[0.8 0.2 0.2]);
end
set(gca,'XTick',x,'XTickLabel',strrep(yrs,'Y',''),'FontSize',12);
xlabel('Year','FontSize',13,'FontWeight','bold');ylabel('SOH (%)','FontSize',13,'FontWeight','bold');
title('Relative dQ (Fixed Denom: Seg8+Seg9): GBM vs XGBoost vs BMS','FontSize',13,'FontWeight','bold');
legend('Location','southwest','FontSize',11);ylim([82 103]);
saveas(fig,fullfile(visDir,'RM_RelativeDQ_Field_SOH.png'));
fprintf('Saved: RM_RelativeDQ_Field_SOH.png\nDone!\n');

function segs=local_find_segments(mask)
    segs=[];n=length(mask);i=1;
    while i<=n,if mask(i),j=i;while j<n&&mask(j+1),j=j+1;end;segs=[segs;i,j];i=j+1;else,i=i+1;end;end
end
