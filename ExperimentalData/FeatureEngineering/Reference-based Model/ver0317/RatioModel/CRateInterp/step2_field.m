% step2_field.m  (CRateInterp)
% 필드 예측: C_eff 기반 C-rate 모델 보간
% crate_models.mat 로드 → C_eff에 해당하는 두 C-rate 모델 선택 → 선형 보간 → SOH

clear; clc; close all;

%% Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
ciDir   = fullfile(verDir, 'RatioModel', 'CRateInterp');
visDir  = fullfile(ciDir, 'Visualization');
rawDir  = fullfile(projDir, 'Rack_raw2mat');
if ~exist(visDir,'dir'), mkdir(visDir); end

pyenv('Version', 'C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

%% Load Models
m = load(fullfile(ciDir,'crate_models.mat'));
CRateModels = m.CRateModels;
crate_vals  = m.crate_vals;   % [0.1, 0.5, 1.0, 2.0, 3.0]
uConds      = m.uConds;
modelNames  = m.modelNames;

%% Retrain Python models if best is XGB/LGB
d  = load(fullfile(verDir,'FeatureMatrix_ver0317.mat')); FM = d.FM;
mr = load(fullfile(verDir,'MasterRulers_PseudoOCV.mat'));
V_chg = mr.MasterRuler_ver0317.V_bounds_chg;
V_dch = mr.MasterRuler_ver0317.V_bounds_dch;
Q_nom = 64;
dQ_c_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i),3:12,'UniformOutput',false);
dQ_d_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i),1:11,'UniformOutput',false);
X_lab_full = [FM{:,dQ_c_cols}, FM{:,dQ_d_cols}, FM.C_eff_chg, FM.C_eff_dch];
y_lab = FM.Static_Capacity / Q_nom * 100;
conditions = FM.Condition;

%% Retrain GBM per C-rate for FIELD prediction (LASSO can't handle NaN)
% Note: LOCO-CV best model is reported for reference, field always uses GBM Surrogate
FieldGBM = struct();
cnames = fieldnames(CRateModels);
fprintf('=== Retraining GBM (Surrogate) per C-rate for field ===\n');
for ci = 1:length(cnames)
    cn   = cnames{ci};
    info = CRateModels.(cn);
    cM   = strcmp(conditions, info.cond);
    X_c  = X_lab_full(cM,:);
    y_c  = y_lab(cM);
    X_c_s = (X_c - info.mu) ./ info.sig;
    mdl_gbm = fitrensemble(X_c_s, y_c, 'Method','LSBoost',...
        'NumLearningCycles',200,'LearnRate',0.05,...
        'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',3,'Surrogate','on'));
    FieldGBM.(cn) = mdl_gbm;
    fprintf('  %s (C-rate=%.1f, best-in-CV: %s) → GBM trained\n', ...
        char(info.cond), info.crate, info.best_model);
end
fprintf('\n');


%% Helper: predict with GBM Surrogate for one C-rate (field safe)
function soh = predict_field_crate(FieldGBM, cn, info, Xf)
    Xfs = (Xf - info.mu) ./ info.sig;
    soh = predict(FieldGBM.(cn), Xfs);
end

%% Helper: find C-rate bracket and interpolate
function soh = interp_crate(crate_vals, cnames, CRateModels, FieldGBM, Xf, ceff)
    ceff = max(crate_vals(1), min(crate_vals(end), ceff));
    idx = find(crate_vals <= ceff, 1,'last');
    if idx >= length(crate_vals), idx = length(crate_vals)-1; end
    idx_lo = idx; idx_hi = idx+1;
    w = (ceff - crate_vals(idx_lo)) / (crate_vals(idx_hi) - crate_vals(idx_lo));

    soh_lo = predict_field_crate(FieldGBM, cnames{idx_lo}, CRateModels.(cnames{idx_lo}), Xf);
    soh_hi = predict_field_crate(FieldGBM, cnames{idx_hi}, CRateModels.(cnames{idx_hi}), Xf);
    soh = (1-w)*soh_lo + w*soh_hi;

    fprintf('    C_eff=%.2f → [%s:%.1f%% × %.2f] + [%s:%.1f%% × %.2f] = %.1f%%\n',...
        ceff, char(CRateModels.(cnames{idx_lo}).cond), soh_lo, 1-w,...
        char(CRateModels.(cnames{idx_hi}).cond), soh_hi, w, soh);
end

%% Field Demo
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
        rs=nan(size(I_r)); if isfield(D,'SOHPct'),rs=D.SOHPct(:);end
    else
        D=S.Raw;
        if isduration(D.Date_Time),t=bd+D.Date_Time;else,t=datetime(D.Date_Time);end
        I_r=D.DCCurrent(:);V=D.CVavg(:);
        rs=nan(size(I_r)); if isfield(D,'SOH_BMS'),rs=D.SOH_BMS(:);end
    end
    ts=seconds(t-t(1));Ic=I_r/Np;vs=rs;vs(vs<=0)=NaN;bms=median(vs,'omitnan');
    cS=local_find_segments(Ic>thr_A);dS=local_find_segments(Ic<-thr_A);
    if ~isempty(cS),dur=cS(:,2)-cS(:,1)+1;cS=cS(dur>=min_chg_sec(k),:);end
    if ~isempty(dS),dur=dS(:,2)-dS(:,1)+1;dS=dS(dur>=min_dch_sec(k),:);end

    % Charge
    dQ_c12=nan(1,12); Ceff_c=NaN;
    if ~isempty(cS)
        [~,bi]=max(cS(:,2)-cS(:,1));cs=cS(bi,1);ce=cS(bi,2);
        vc=V(cs:ce);ic=Ic(cs:ce);tc_=ts(cs:ce);
        qc=cumtrapz(tc_,abs(ic))/3600;
        vm=vc(1);qm=qc(1);for ii=2:length(vc),if vc(ii)>vm(end),vm(end+1)=vc(ii);qm(end+1)=qc(ii);end;end
        if length(vm)>1,QR=interp1(vm,qm,V_chg,'linear',NaN);dQ_c12=abs(diff(QR));end
        Ceff_c=mean(abs(ic))/Q_nom;
    end
    % Discharge
    dQ_d12=nan(1,12); Ceff_d=NaN;
    if ~isempty(dS)
        [~,bi]=max(dS(:,2)-dS(:,1));ds=dS(bi,1);de=dS(bi,2);
        vd=V(ds:de);id=Ic(ds:de);td_=ts(ds:de);
        qd=cumtrapz(td_,abs(id))/3600;
        vmd=vd(1);qmd=qd(1);for ii=2:length(vd),if vd(ii)<vmd(end),vmd(end+1)=vd(ii);qmd(end+1)=qd(ii);end;end
        va=flip(vmd);qa=flip(qmd);
        if length(va)>1,QRd=interp1(va,qa,V_dch,'linear',NaN);dQ_d12=abs(diff(QRd));end
        Ceff_d=mean(abs(id))/Q_nom;
    end

    % Feature vector: dQ ratio (each seg / sum valid) + C_eff
    dQ_c_f = dQ_c12(3:12); dQ_d_f = dQ_d12(1:11);
    rdQ_c_f = dQ_c_f / sum(dQ_c_f,'omitnan');
    rdQ_d_f = dQ_d_f / sum(dQ_d_f,'omitnan');
    Xf = [rdQ_c_f, rdQ_d_f, Ceff_c, Ceff_d];
    nan_cnt = sum(isnan(Xf));
    fprintf('  Ceff_c=%.2f, Ceff_d=%.2f, NaN=%d/23\n', Ceff_c, Ceff_d, nan_cnt);

    % Use mean C_eff for interpolation
    ceff_mean = mean([Ceff_c, Ceff_d], 'omitnan');

    fprintf('  Interpolating (mean Ceff=%.2f):\n', ceff_mean);
    soh_pred = interp_crate(crate_vals, cnames, CRateModels, FieldGBM, Xf, ceff_mean);

    fprintf('  BMS=%.1f%%  Pred=%.1f%%  Err=%.1f%%\n', bms, soh_pred, soh_pred-bms);
    FR.(yr)=struct('BMS',bms,'Pred',soh_pred,'Err',soh_pred-bms,'Ceff_c',Ceff_c,'Ceff_d',Ceff_d,'NaN',nan_cnt);
end

%% Summary
fprintf('\n========== SUMMARY ==========\n');
fprintf('%-6s %-6s %-7s %-7s %-10s\n','Year','BMS','Pred','Err','Ceff(c/d)');
yrs=fieldnames(FR);
for k=1:length(yrs)
    r=FR.(yrs{k});
    fprintf('%-6s %-6.1f %-7.1f %-7.1f %.2f/%.2f\n',yrs{k},r.BMS,r.Pred,r.Err,r.Ceff_c,r.Ceff_d);
end

%% Plot
fig=figure('Position',[50 100 900 450]);
n=length(yrs); x=1:n;
bv=arrayfun(@(i)FR.(yrs{i}).BMS, 1:n);
pv=arrayfun(@(i)FR.(yrs{i}).Pred,1:n);
hold on;grid on;box on;
plot(x,bv,'k--s','LineWidth',2.5,'MarkerSize',12,'MarkerFaceColor','k','DisplayName','BMS SOH');
plot(x,pv,'-o','Color',[0.2 0.6 0.9],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor',[0.2 0.6 0.9],'DisplayName','C-rate Interp Model');
for i=1:n
    text(x(i)+0.05,bv(i)+0.5,sprintf('%.1f%%',bv(i)),'FontSize',9,'FontWeight','bold');
    text(x(i)+0.05,pv(i)-0.8,sprintf('%.1f%%',pv(i)),'FontSize',8,'Color',[0.1 0.5 0.8]);
end
set(gca,'XTick',x,'XTickLabel',strrep(yrs,'Y',''),'FontSize',12);
xlabel('Year','FontSize',13,'FontWeight','bold'); ylabel('SOH (%)','FontSize',13,'FontWeight','bold');
title('C-rate Interpolation Model: Field Demo','FontSize',13,'FontWeight','bold');
legend('Location','southwest','FontSize',11); ylim([82 103]);
saveas(fig,fullfile(visDir,'CRateInterp_Field_SOH.png'));
fprintf('Saved: CRateInterp_Field_SOH.png\nDone!\n');

function segs=local_find_segments(mask)
    segs=[];n=length(mask);i=1;
    while i<=n,if mask(i),j=i;while j<n&&mask(j+1),j=j+1;end;segs=[segs;i,j];i=j+1;else,i=i+1;end;end
end
