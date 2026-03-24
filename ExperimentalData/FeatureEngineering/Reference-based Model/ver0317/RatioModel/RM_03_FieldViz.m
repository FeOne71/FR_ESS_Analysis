% RM_03_FieldViz.m
% RM_01보다 깔끔한 필드 시연 시각화 (마스킹 증강 없음)
% 5개 모델 × Y2021/2023/2024/2025 SOH 예측 vs BMS SOH

clear; clc; close all;

%% Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData','FeatureEngineering','Lab_RPT_Analysis','ver0317');
rmDir   = fullfile(verDir, 'RatioModel');
visDir  = fullfile(rmDir, 'Visualization');
rawDir  = fullfile(projDir, 'Rack_raw2mat');
Q_nom   = 64; Np = 2; num_segs = 12;

pyenv('Version','C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

%% Lab Training
d = load(fullfile(verDir,'FeatureMatrix_ver0317.mat')); FM = d.FM;
mr = load(fullfile(verDir,'MasterRulers_PseudoOCV.mat'));
V_chg = mr.MasterRuler_ver0317.V_bounds_chg;
V_dch = mr.MasterRuler_ver0317.V_bounds_dch;

chg_cols = arrayfun(@(i)sprintf('dQ_c_%02d',i),3:12,'UniformOutput',false);
dch_cols = arrayfun(@(i)sprintf('dQ_d_%02d',i),1:11,'UniformOutput',false);
dQ_c_ratio = FM{:,chg_cols} ./ sum(FM{:,chg_cols},2,'omitnan');
dQ_d_ratio = FM{:,dch_cols} ./ sum(FM{:,dch_cols},2,'omitnan');
X_lab = [dQ_c_ratio, dQ_d_ratio, FM.C_eff_chg, FM.C_eff_dch];
y_lab = FM.Static_Capacity / Q_nom * 100;
mu_lab = mean(X_lab,1,'omitnan'); sigma_lab = std(X_lab,0,1,'omitnan');
sigma_lab(sigma_lab==0)=1;
X_lab_s = (X_lab - mu_lab) ./ sigma_lab; X_lab_s(isnan(X_lab_s))=0;

fprintf('Training 5 models on full Lab data...\n');
[B,FI] = lasso(X_lab_s, y_lab, 'CV',4);
M.LASSO.coef = B(:,FI.IndexMinMSE); M.LASSO.ic = FI.Intercept(FI.IndexMinMSE);
M.RF  = fitrensemble(X_lab_s,y_lab,'Method','Bag','NumLearningCycles',100,...
    'Learners',templateTree('MaxNumSplits',20,'MinLeafSize',5,'Surrogate','on'));
M.GBM = fitrensemble(X_lab_s,y_lab,'Method','LSBoost','NumLearningCycles',100,...
    'LearnRate',0.1,'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',5,'Surrogate','on'));
M.XGB = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
    'learning_rate',0.1,'min_child_weight',int32(3),'subsample',0.8,'random_state',int32(42)));
M.XGB.fit(py.numpy.array(X_lab_s),py.numpy.array(y_lab));
M.LGB = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
    'learning_rate',0.1,'min_child_samples',int32(5),'subsample',0.8,...
    'random_state',int32(42),'verbose',int32(-1)));
M.LGB.fit(py.numpy.array(X_lab_s),py.numpy.array(y_lab));
fprintf('  Done.\n');

%% Field Data Processing
dataFiles = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021', 'old', datetime(2021,6,3),  [600 300];
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023', 'new', datetime(2023,10,16),[300 150];
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024', 'new', datetime(2024,9,9),  [300 300];
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'), 'Y2025', 'new', datetime(2025,7,11), [300 150];
};
years = {'Y2021','Y2023','Y2024','Y2025'};
year_x = [2021, 2023, 2024, 2025];

FR = struct(); % FieldResults

for k = 1:size(dataFiles,1)
    fpath=dataFiles{k,1}; yr=dataFiles{k,2}; dtype=dataFiles{k,3};
    base_date=dataFiles{k,4}; min_dur=dataFiles{k,5};
    if ~exist(fpath,'file'), fprintf('  Missing: %s\n',fpath); continue; end
    fprintf('Processing %s...\n', yr);

    S=load(fpath);
    if strcmp(dtype,'old')
        D=S.Raw.Rack01;
        if isfield(D,'Time'), t=datetime(D.Time);
        elseif isduration(D.Date_Time), t=base_date+D.Date_Time;
        else, t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent_A(:); V_avg=D.AverageCV_V(:);
        raw_soh = getfield_safe(D,'SOHPct');
    else
        D=S.Raw;
        if isduration(D.Date_Time), t=base_date+D.Date_Time; else, t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent(:); V_avg=D.CVavg(:);
        raw_soh = getfield_safe(D,'SOH_BMS');
    end

    tsec=seconds(t-t(1)); I_cell=I_rack/Np;
    thr_A=Q_nom*0.05/Np;
    vs=raw_soh; vs(vs<=0)=NaN; soh_bms=median(vs,'omitnan');

    chgS=local_segs(I_cell>thr_A);  if ~isempty(chgS), dur=chgS(:,2)-chgS(:,1)+1; chgS=chgS(dur>=min_dur(1),:); end
    dchS=local_segs(I_cell<-thr_A); if ~isempty(dchS), dur=dchS(:,2)-dchS(:,1)+1; dchS=dchS(dur>=min_dur(2),:); end
    if isempty(chgS), fprintf('  No charge segs, skip.\n'); continue; end

    % Best charge segment (smoothed V, movmean 30)
    [~,bi]=max(chgS(:,2)-chgS(:,1));
    v_c=V_avg(chgS(bi,1):chgS(bi,2)); i_c=I_cell(chgS(bi,1):chgS(bi,2)); t_c=tsec(chgS(bi,1):chgS(bi,2));
    v_c=smoothdata(v_c,'movmean',30);
    q_c=cumtrapz(t_c,abs(i_c))/3600;
    vm=v_c(1); qm=q_c(1);
    for ii=2:length(v_c), if v_c(ii)>vm(end), vm(end+1)=v_c(ii); qm(end+1)=q_c(ii); end; end
    dQ_chg=nan(1,num_segs);
    if length(vm)>1, QR=interp1(vm,qm,V_chg,'linear',NaN); dQ_chg=abs(diff(QR)); end
    C_eff_chg_v=mean(abs(i_c))/Q_nom;

    dQ_dch=nan(1,num_segs); C_eff_dch_v=NaN;
    if ~isempty(dchS)
        [~,bi_d]=max(dchS(:,2)-dchS(:,1));
        v_d=V_avg(dchS(bi_d,1):dchS(bi_d,2)); i_d=I_cell(dchS(bi_d,1):dchS(bi_d,2)); t_d=tsec(dchS(bi_d,1):dchS(bi_d,2));
        v_d=smoothdata(v_d,'movmean',30);
        q_d=cumtrapz(t_d,abs(i_d))/3600;
        vd=v_d(1); qd=q_d(1);
        for ii=2:length(v_d), if v_d(ii)<vd(end), vd(end+1)=v_d(ii); qd(end+1)=q_d(ii); end; end
        va=flip(vd); qa=flip(qd);
        if length(va)>1, QR_d=interp1(va,qa,V_dch,'linear',NaN); dQ_dch=abs(diff(QR_d)); end
        C_eff_dch_v=mean(abs(i_d))/Q_nom;
    end

    dc=dQ_chg(3:12); dd=dQ_dch(1:11);
    dc_r=dc/sum(dc,'omitnan'); dd_r=dd/sum(dd,'omitnan');
    X_f=[dc_r, dd_r, C_eff_chg_v, C_eff_dch_v];
    X_fs=(X_f-mu_lab)./sigma_lab;
    X_imp=X_fs; X_imp(isnan(X_imp))=0;

    FR.(yr).BMS   = soh_bms;
    FR.(yr).LASSO = X_imp*M.LASSO.coef + M.LASSO.ic;
    FR.(yr).RF    = predict(M.RF, X_imp);
    FR.(yr).GBM   = predict(M.GBM, X_imp);
    X_py=py.numpy.array(X_fs).reshape(int32(1),int32(-1));
    FR.(yr).XGB   = double(M.XGB.predict(X_py));
    FR.(yr).LGB   = double(M.LGB.predict(X_py));
    fprintf('  BMS=%.1f%% | LASSO=%.1f%% RF=%.1f%% GBM=%.1f%% XGB=%.1f%% LGB=%.1f%%\n',...
        soh_bms, FR.(yr).LASSO, FR.(yr).RF, FR.(yr).GBM, FR.(yr).XGB, FR.(yr).LGB);
end
save(fullfile(rmDir,'RM01_FieldResults_NoMask.mat'),'FR');

%% ===== Visualization =====
modelNames = {'LASSO','RF','GBM','XGB','LGB'};
mColors = {[0.22 0.43 0.74],[0.29 0.62 0.35],[0.85 0.40 0.30],[0.59 0.44 0.72],[0.95 0.70 0.25]};
mMarkers= {'o','s','d','^','v'};
bms_color = [0.1 0.1 0.1];

% Collect data
yr_avail = {}; yr_xi = [];
bms_vals = []; pred = zeros(numel(years),5);
for k=1:numel(years)
    yr=years{k};
    if isfield(FR,yr)
        yr_avail{end+1}=yr; yr_xi(end+1)=year_x(k);
        bms_vals(end+1) = FR.(yr).BMS;
        pred(k,1)=FR.(yr).LASSO; pred(k,2)=FR.(yr).RF; pred(k,3)=FR.(yr).GBM;
        pred(k,4)=FR.(yr).XGB; pred(k,5)=FR.(yr).LGB;
    end
end
% Remove empty rows
valid = bms_vals>0;
yr_avail = yr_avail(valid); yr_xi = yr_xi(valid);
bms_vals = bms_vals(valid); pred = pred(valid,:);

n = length(yr_avail);

%% Fig 1: Combined (all 5 models + BMS)
fig1=figure('Position',[50 50 920 480],'Color','w');
ax=axes('Parent',fig1,'FontSize',11);
hold(ax,'on'); grid(ax,'on'); box(ax,'on');
set(ax,'GridColor',[0.85 0.85 0.85],'GridAlpha',1,'TickDir','out');

% BMS
plot(ax, yr_xi, bms_vals, 'k--o','LineWidth',2.5,'MarkerSize',9,'MarkerFaceColor','k','DisplayName','BMS SOH');

% Models
for m=1:5
    plot(ax, yr_xi, pred(:,m), '-','Color',mColors{m},'LineWidth',1.8,...
        'Marker',mMarkers{m},'MarkerSize',8,'MarkerFaceColor',mColors{m},...
        'DisplayName',modelNames{m});
end

% Labels on BMS
for i=1:n
    text(ax, yr_xi(i), bms_vals(i)+0.3, sprintf('%.1f%%',bms_vals(i)),...
        'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',bms_color);
end

set(ax,'XTick',yr_xi,'XTickLabel',yr_avail,'XLim',[yr_xi(1)-0.5, yr_xi(end)+0.5]);
ylabel(ax,'SOH (%)','FontSize',12);
title(ax,'Field SOH Estimation: RM-01 (No Masking Augmentation)','FontSize',13,'FontWeight','bold');
legend(ax,'Location','southwest','FontSize',10,'Box','off');
ylim(ax,[min([bms_vals pred(:)'])-3, max([bms_vals pred(:)'])+3]);

saveas(fig1, fullfile(visDir,'RM01_Field_Combined.png'));
fprintf('Saved: RM01_Field_Combined.png\n');

%% Fig 2: Per-model subplots (5 panels)
fig2=figure('Position',[50 50 1400 300],'Color','w');
sgtitle('Field SOH Estimation per Model (RM-01, No Masking Augmentation)','FontSize',12,'FontWeight','bold');

for m=1:5
    ax=subplot(1,5,m);
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    set(ax,'GridColor',[0.88 0.88 0.88],'GridAlpha',1,'TickDir','out','FontSize',9);

    fill([yr_xi fliplr(yr_xi)], [bms_vals-1 fliplr(bms_vals+1)],[0.8 0.8 0.8],'FaceAlpha',0.3,'EdgeColor','none');
    plot(ax, yr_xi, bms_vals,'k--o','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','k');
    plot(ax, yr_xi, pred(:,m),'-','Color',mColors{m},'LineWidth',2,...
        'Marker',mMarkers{m},'MarkerSize',8,'MarkerFaceColor',mColors{m});

    err = pred(:,m) - bms_vals';
    rmse_f = sqrt(mean(err.^2));
    mae_f  = mean(abs(err));
    maxe_f = max(abs(err));

    % Value labels — larger font
    for i=1:n
        text(ax,yr_xi(i),pred(i,m)+0.9,sprintf('%.1f%%',pred(i,m)),...
            'HorizontalAlignment','center','FontSize',11,'FontWeight','bold','Color',mColors{m});
    end

    % RMSE/MAE/MaxErr text box
    if m==1
        txt_pos = [0.45 0.90];
    else
        txt_pos = [0.05 0.16];
    end
    text(ax,txt_pos(1),txt_pos(2),sprintf('RMSE=%.2f%%\nMAE =%.2f%%\nMax =%.2f%%',rmse_f,mae_f,maxe_f),...
        'Units','normalized','FontSize',8,'Color',[0.2 0.2 0.2],...
        'BackgroundColor',[0.97 0.97 0.97],'EdgeColor',[0.8 0.8 0.8],...
        'HorizontalAlignment','left');

    title(ax, modelNames{m},'FontSize',11,'FontWeight','bold','Color',mColors{m});
    set(ax,'XTick',yr_xi,'XTickLabel',yr_avail);
    ylabel(ax,'SOH (%)','FontSize',9);

    if m == 1
        % LASSO: auto y-scale + warning note
        text(ax,0.5,0.5,{'※ 필드에서 발산','(NaN imputation 문제)'},'Units','normalized',...
            'HorizontalAlignment','center','FontSize',9,'Color',[0.7 0.1 0.1],...
            'FontWeight','bold','BackgroundColor',[1 0.95 0.95],'EdgeColor',[0.8 0.3 0.3]);
        legend(ax,{'±1% band','BMS SOH','LASSO'},'FontSize',7,'Location','southeast','Box','off');
    else
        ylim(ax,[90 105]);
        legend(ax,{'±1% band','BMS SOH',modelNames{m}},'FontSize',7,'Location','southeast','Box','off');
    end
end

saveas(fig2, fullfile(visDir,'RM01_Field_PerModel.png'));
saveas(fig2, fullfile(visDir,'RM01_Field_PerModel.fig'));

fprintf('Saved: RM01_Field_PerModel.png\n');

fprintf('\nDone.\n');

%% Helpers
function segs = local_segs(mask)
    d=[0;diff(mask(:))];
    starts=find(d==1); ends=find(d==-1)-1;
    if mask(1), starts=[1;starts]; end
    if mask(end), ends=[ends;length(mask)]; end
    if isempty(starts)||isempty(ends), segs=[]; return; end
    segs=[starts, ends];
end

function v=getfield_safe(S,fname)
    if isfield(S,fname), v=S.(fname)(:); else, v=nan; end
end
