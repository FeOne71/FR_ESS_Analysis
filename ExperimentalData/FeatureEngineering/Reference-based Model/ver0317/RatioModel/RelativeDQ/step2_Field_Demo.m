% step2_Field_Demo.m
% RelativeDQ 모델 필드 시연
% rdQ_model.mat 로드 → 필드 4개년 relative dQ 계산 → GBM + XGBoost 예측 → 시각화

clear; clc; close all;

%% Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
rdqDir  = fullfile(verDir, 'RatioModel', 'RelativeDQ');
visDir  = fullfile(rdqDir, 'Visualization');
rawDir  = fullfile(projDir, 'Rack_raw2mat');
if ~exist(visDir,'dir'), mkdir(visDir); end

pyenv('Version', 'C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

%% Load Model (from step1)
m = load(fullfile(rdqDir,'rdQ_model.mat'));
mdl_GBM = m.mdl_GBM_final;
mu_lab  = m.mu_lab;
sig_lab = m.sig_lab;
fprintf('Loaded rdQ_model.mat\n');

%% Load Lab data for XGBoost retraining (Python 모델은 .mat 저장 불가)
d  = load(fullfile(verDir,'FeatureMatrix_ver0317.mat')); FM=d.FM;
mr = load(fullfile(verDir,'MasterRulers_PseudoOCV.mat'));
V_chg = mr.MasterRuler_ver0317.V_bounds_chg;
V_dch = mr.MasterRuler_ver0317.V_bounds_dch;
Q_nom = 64;

dQ_c_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i),3:12,'UniformOutput',false);
dQ_d_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i),1:11,'UniformOutput',false);
dQ_c_raw = FM{:,dQ_c_cols}; dQ_d_raw = FM{:,dQ_d_cols};
denom_c  = dQ_c_raw(:,6)+dQ_c_raw(:,7); rdQ_c = dQ_c_raw./denom_c;
denom_d  = dQ_d_raw(:,8)+dQ_d_raw(:,9); rdQ_d = dQ_d_raw./denom_d;
X_lab    = [rdQ_c,rdQ_d,FM.C_eff_chg,FM.C_eff_dch];
y_lab    = FM.Static_Capacity/Q_nom*100;
X_lab_s  = (X_lab-mu_lab)./sig_lab;

fprintf('Retraining XGBoost...\n');
mdl_XGB = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(5),...
    'learning_rate',0.05,'min_child_weight',int32(3),'subsample',0.8,...
    'colsample_bytree',0.8,'random_state',int32(42)));
mdl_XGB.fit(py.numpy.array(X_lab_s), py.numpy.array(y_lab));
fprintf('  Done.\n');

%% Field Prediction
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
    if ~exist(fpath,'file'), fprintf('%s not found, skip\n',yr); continue; end
    S=load(fpath);
    if strcmp(dtype,'old')
        D=S.Raw.Rack01;
        if isfield(D,'Time'),t=datetime(D.Time);elseif isduration(D.Date_Time),t=bd+D.Date_Time;else,t=datetime(D.Date_Time);end
        I_r=D.DCCurrent_A(:); V=D.AverageCV_V(:);
        rs=nan(size(I_r)); if isfield(D,'SOHPct'),rs=D.SOHPct(:);end
    else
        D=S.Raw;
        if isduration(D.Date_Time),t=bd+D.Date_Time;else,t=datetime(D.Date_Time);end
        I_r=D.DCCurrent(:); V=D.CVavg(:);
        rs=nan(size(I_r)); if isfield(D,'SOH_BMS'),rs=D.SOH_BMS(:);end
    end
    ts=seconds(t-t(1)); Ic=I_r/Np; vs=rs; vs(vs<=0)=NaN; bms=median(vs,'omitnan');
    cS=local_find_segments(Ic>thr_A); dS=local_find_segments(Ic<-thr_A);
    if ~isempty(cS),dur=cS(:,2)-cS(:,1)+1;cS=cS(dur>=min_chg_sec(k),:);end
    if ~isempty(dS),dur=dS(:,2)-dS(:,1)+1;dS=dS(dur>=min_dch_sec(k),:);end

    % Charge V-Q
    dQ_c12=nan(1,12); Ceff_c=NaN;
    if ~isempty(cS)
        [~,bi]=max(cS(:,2)-cS(:,1)); cs=cS(bi,1); ce=cS(bi,2);
        vc=V(cs:ce); ic=Ic(cs:ce); tc_=ts(cs:ce);
        qc=cumtrapz(tc_,abs(ic))/3600;
        vm=vc(1);qm=qc(1);for ii=2:length(vc),if vc(ii)>vm(end),vm(end+1)=vc(ii);qm(end+1)=qc(ii);end;end
        if length(vm)>1, QR=interp1(vm,qm,V_chg,'linear',NaN); dQ_c12=abs(diff(QR)); end
        Ceff_c=mean(abs(ic))/Q_nom;
    end
    % Discharge V-Q
    dQ_d12=nan(1,12); Ceff_d=NaN;
    if ~isempty(dS)
        [~,bi]=max(dS(:,2)-dS(:,1)); ds=dS(bi,1); de=dS(bi,2);
        vd=V(ds:de); id=Ic(ds:de); td_=ts(ds:de);
        qd=cumtrapz(td_,abs(id))/3600;
        vmd=vd(1);qmd=qd(1);for ii=2:length(vd),if vd(ii)<vmd(end),vmd(end+1)=vd(ii);qmd(end+1)=qd(ii);end;end
        va=flip(vmd);qa=flip(qmd);
        if length(va)>1, QRd=interp1(va,qa,V_dch,'linear',NaN); dQ_d12=abs(diff(QRd)); end
        Ceff_d=mean(abs(id))/Q_nom;
    end

    % Relative dQ (same formula as Lab)
    dQ_cf=dQ_c12(3:12); denom_cf=dQ_cf(6)+dQ_cf(7); rdQ_cf=dQ_cf/denom_cf;
    dQ_df=dQ_d12(1:11); denom_df=dQ_df(8)+dQ_df(9); rdQ_df=dQ_df/denom_df;
    Xf=[rdQ_cf,rdQ_df,Ceff_c,Ceff_d]; Xfs=(Xf-mu_lab)./sig_lab;

    soh_gbm=predict(mdl_GBM,Xfs);
    Xpy=py.numpy.array(Xfs).reshape(int32(1),int32(-1));
    soh_xgb=double(mdl_XGB.predict(Xpy));

    nan_c=sum(isnan(rdQ_cf)); nan_d=sum(isnan(rdQ_df));
    fprintf('%s: BMS=%.1f%%  GBM=%.1f%%  XGB=%.1f%%  NaN(c=%d,d=%d)  Ceff=[%.2f,%.2f]\n',...
        yr,bms,soh_gbm,soh_xgb,nan_c,nan_d,Ceff_c,Ceff_d);
    FR.(yr)=struct('BMS',bms,'GBM',soh_gbm,'XGB',soh_xgb,'NaNc',nan_c,'NaNd',nan_d,'Ceff_c',Ceff_c,'Ceff_d',Ceff_d);
end

%% Summary
fprintf('\n========== SUMMARY ==========\n');
fprintf('%-6s %-6s %-7s %-7s %-7s %-7s\n','Year','BMS','GBM','GBM_err','XGB','XGB_err');
yrs=fieldnames(FR);
for k=1:length(yrs)
    r=FR.(yrs{k});
    fprintf('%-6s %-6.1f %-7.1f %-7.1f %-7.1f %-7.1f\n',...
        yrs{k},r.BMS,r.GBM,r.GBM-r.BMS,r.XGB,r.XGB-r.BMS);
end

%% Plot
fig=figure('Position',[50 100 950 480]);
n=length(yrs); x=1:n;
bv=arrayfun(@(i)FR.(yrs{i}).BMS,1:n);
gv=arrayfun(@(i)FR.(yrs{i}).GBM,1:n);
xv=arrayfun(@(i)FR.(yrs{i}).XGB,1:n);
hold on; grid on; box on;
plot(x,bv,'k--s','LineWidth',2.5,'MarkerSize',12,'MarkerFaceColor','k','DisplayName','BMS SOH');
plot(x,gv,'-o','Color',[0.2 0.7 0.3],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor',[0.2 0.7 0.3],'DisplayName','GBM (Surrogate)');
plot(x,xv,'-^','Color',[0.8 0.2 0.2],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor',[0.8 0.2 0.2],'DisplayName','XGBoost');
for i=1:n
    text(x(i)+0.05,bv(i)+0.5,sprintf('%.1f%%',bv(i)),'FontSize',9,'FontWeight','bold');
    text(x(i)+0.05,gv(i)-0.8,sprintf('%.1f%%',gv(i)),'FontSize',8,'Color',[0.2 0.7 0.3]);
    text(x(i)+0.05,xv(i)+0.5,sprintf('%.1f%%',xv(i)),'FontSize',8,'Color',[0.8 0.2 0.2]);
end
set(gca,'XTick',x,'XTickLabel',strrep(yrs,'Y',''),'FontSize',12);
xlabel('Year','FontSize',13,'FontWeight','bold');
ylabel('SOH (%)','FontSize',13,'FontWeight','bold');
title('Relative dQ (Fixed Denom: Seg8+Seg9): Field Demo','FontSize',13,'FontWeight','bold');
legend('Location','southwest','FontSize',11); ylim([82 103]);
saveas(fig,fullfile(visDir,'rdQ_Field_SOH.png'));
fprintf('Saved: rdQ_Field_SOH.png\nDone!\n');

function segs=local_find_segments(mask)
    segs=[];n=length(mask);i=1;
    while i<=n,if mask(i),j=i;while j<n&&mask(j+1),j=j+1;end;segs=[segs;i,j];i=j+1;else,i=i+1;end;end
end
