% MA_01_Train_and_Field.m
% Masking Augmentation 실험
% Lab 학습 데이터에 랜덤 세그먼트 마스킹을 적용하여 필드 OOD 대응력 향상
% -------------------------------------------------------------------------
% 마스킹 증강 전략:
%   - 각 Lab 샘플에서 충전/방전 세그먼트를 무작위로 NaN 처리
%   - 마스킹 후 남은 세그먼트 비율로 재정규화
%   - 증강 데이터를 원본과 합산하여 학습
% -------------------------------------------------------------------------

clear; clc; close all;

%% Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir,'ExperimentalData','FeatureEngineering','Lab_RPT_Analysis','ver0317');
maDir   = fullfile(verDir,'MaskingAugmentation');
visDir  = fullfile(maDir,'Visualization');
rawDir  = fullfile(projDir,'Rack_raw2mat');
Q_nom=64; Np=2; num_segs=12;

pyenv('Version','C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

%% ===== Part 1: Lab Data Loading =====
d=load(fullfile(verDir,'FeatureMatrix_ver0317.mat')); FM=d.FM;
mr=load(fullfile(verDir,'MasterRulers_PseudoOCV.mat'));
V_chg=mr.MasterRuler_ver0317.V_bounds_chg; V_dch=mr.MasterRuler_ver0317.V_bounds_dch;

chg_cols=arrayfun(@(i)sprintf('dQ_c_%02d',i),3:12,'UniformOutput',false);
dch_cols=arrayfun(@(i)sprintf('dQ_d_%02d',i),1:11,'UniformOutput',false);
dQ_c_orig=FM{:,chg_cols}./sum(FM{:,chg_cols},2,'omitnan');
dQ_d_orig=FM{:,dch_cols}./sum(FM{:,dch_cols},2,'omitnan');
X_orig=[dQ_c_orig, dQ_d_orig, FM.C_eff_chg, FM.C_eff_dch];
y_orig=FM.Static_Capacity/Q_nom*100;

N_orig=size(X_orig,1);
n_chg=10; n_dch=11; % feature counts

fprintf('Lab data: %d samples, %d features\n',N_orig,size(X_orig,2));

%% ===== Part 2: Masking Augmentation =====
rng(42);
N_aug=20;  % augmented copies per original sample
X_aug_all=[]; y_aug_all=[];

for i=1:N_orig
    xi=X_orig(i,:); yi=y_orig(i);
    xc=xi(1:n_chg); xd=xi(n_chg+1:n_chg+n_dch);
    Ceff_c=xi(end-1); Ceff_d=xi(end);

    for a=1:N_aug
        % Randomly mask 0 to (n-1) charge segments
        n_mask_c=randi([0 n_chg-1]);
        n_mask_d=randi([0 n_dch-1]);

        xc_a=xc; xd_a=xd;
        if n_mask_c>0
            idx_c=randperm(n_chg,n_mask_c);
            xc_a(idx_c)=NaN;
        end
        if n_mask_d>0
            idx_d=randperm(n_dch,n_mask_d);
            xd_a(idx_d)=NaN;
        end

        % Re-normalize after masking
        sc=sum(xc_a,'omitnan'); if sc>0, xc_a=xc_a/sc; end
        sd=sum(xd_a,'omitnan'); if sd>0, xd_a=xd_a/sd; end

        X_aug_all(end+1,:)=[xc_a, xd_a, Ceff_c, Ceff_d];
        y_aug_all(end+1)=yi;
    end
end

% Combine original + augmented
X_train=[X_orig; X_aug_all];
y_train=[y_orig; y_aug_all'];

fprintf('Training set: %d samples (%d orig + %d augmented)\n',...
    size(X_train,1), N_orig, size(X_aug_all,1));

%% ===== Part 3: Normalization =====
mu_tr=mean(X_train,1,'omitnan');
sg_tr=std(X_train,0,1,'omitnan');
sg_tr(sg_tr==0)=1;
X_tr_s=(X_train-mu_tr)./sg_tr;
X_tr_imp=X_tr_s; X_tr_imp(isnan(X_tr_imp))=0;

%% ===== Part 4: Model Training =====
fprintf('\nTraining 5 models with Masking Augmentation...\n');

% LASSO
[B,FI]=lasso(X_tr_imp, y_train,'CV',4);
M.LASSO.coef=B(:,FI.IndexMinMSE); M.LASSO.ic=FI.Intercept(FI.IndexMinMSE);
fprintf('  LASSO done.\n');

% RF
M.RF=fitrensemble(X_tr_imp,y_train,'Method','Bag','NumLearningCycles',100,...
    'Learners',templateTree('MaxNumSplits',20,'MinLeafSize',5,'Surrogate','on'));
fprintf('  RF done.\n');

% GBM
M.GBM=fitrensemble(X_tr_imp,y_train,'Method','LSBoost','NumLearningCycles',100,...
    'LearnRate',0.1,'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',5,'Surrogate','on'));
fprintf('  GBM done.\n');

% XGBoost (pass NaN directly — now trained with NaN!)
X_tr_xgb=py.numpy.array(X_tr_s);  % NaN included
y_tr_py=py.numpy.array(y_train);
M.XGB=py.xgboost.XGBRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
    'learning_rate',0.1,'min_child_weight',int32(3),'subsample',0.8,'random_state',int32(42)));
M.XGB.fit(X_tr_xgb, y_tr_py);
fprintf('  XGBoost done.\n');

% LightGBM (pass NaN)
M.LGB=py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
    'learning_rate',0.1,'min_child_samples',int32(5),'subsample',0.8,...
    'random_state',int32(42),'verbose',int32(-1)));
M.LGB.fit(X_tr_xgb, y_tr_py);
fprintf('  LightGBM done.\n');

%% ===== Part 5: Field Data Processing =====
dataFiles={
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'),'Y2021','old',datetime(2021,6,3),[600 300];
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'),'Y2023','new',datetime(2023,10,16),[300 150];
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'),'Y2024','new',datetime(2024,9,9),[300 300];
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'),'Y2025','new',datetime(2025,7,11),[300 150];
};
years={'Y2021','Y2023','Y2024','Y2025'};
year_x=[2021,2023,2024,2025];

FR=struct();

for k=1:size(dataFiles,1)
    fpath=dataFiles{k,1}; yr=dataFiles{k,2}; dtype=dataFiles{k,3};
    bd=dataFiles{k,4}; mdur=dataFiles{k,5};
    if ~exist(fpath,'file'), fprintf('  Missing: %s\n',fpath); continue; end
    fprintf('Processing %s...\n',yr);

    S=load(fpath);
    if strcmp(dtype,'old')
        D=S.Raw.Rack01;
        if isfield(D,'Time'),t=datetime(D.Time);
        elseif isduration(D.Date_Time),t=bd+D.Date_Time; else,t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent_A(:); V_avg=D.AverageCV_V(:);
        raw_soh=getfield_safe(D,'SOHPct');
    else
        D=S.Raw;
        if isduration(D.Date_Time),t=bd+D.Date_Time; else,t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent(:); V_avg=D.CVavg(:);
        raw_soh=getfield_safe(D,'SOH_BMS');
    end
    tsec=seconds(t-t(1)); I_cell=I_rack/Np; thr=Q_nom*0.05/Np;
    vs=raw_soh; vs(vs<=0)=NaN; soh_bms=median(vs,'omitnan');

    chgS=local_segs(I_cell>thr); if ~isempty(chgS),dur=chgS(:,2)-chgS(:,1)+1;chgS=chgS(dur>=mdur(1),:);end
    dchS=local_segs(I_cell<-thr); if ~isempty(dchS),dur=dchS(:,2)-dchS(:,1)+1;dchS=dchS(dur>=mdur(2),:);end
    if isempty(chgS), continue; end

    [~,bi]=max(chgS(:,2)-chgS(:,1));
    v_c=V_avg(chgS(bi,1):chgS(bi,2)); i_c=I_cell(chgS(bi,1):chgS(bi,2)); t_c=tsec(chgS(bi,1):chgS(bi,2));
    v_c=smoothdata(v_c,'movmean',30);
    q_c=cumtrapz(t_c,abs(i_c))/3600;
    vm=v_c(1);qm=q_c(1);
    for ii=2:length(v_c),if v_c(ii)>vm(end),vm(end+1)=v_c(ii);qm(end+1)=q_c(ii);end;end
    dQ_chg=nan(1,num_segs);
    if length(vm)>1,QR=interp1(vm,qm,V_chg,'linear',NaN);dQ_chg=abs(diff(QR));end
    Ce_c=mean(abs(i_c))/Q_nom;

    dQ_dch=nan(1,num_segs); Ce_d=NaN;
    if ~isempty(dchS)
        [~,bi_d]=max(dchS(:,2)-dchS(:,1));
        v_d=V_avg(dchS(bi_d,1):dchS(bi_d,2)); i_d=I_cell(dchS(bi_d,1):dchS(bi_d,2)); t_d=tsec(dchS(bi_d,1):dchS(bi_d,2));
        v_d=smoothdata(v_d,'movmean',30);
        q_d=cumtrapz(t_d,abs(i_d))/3600;
        vd=v_d(1);qd=q_d(1);
        for ii=2:length(v_d),if v_d(ii)<vd(end),vd(end+1)=v_d(ii);qd(end+1)=q_d(ii);end;end
        va=flip(vd);qa=flip(qd);
        if length(va)>1,QR_d=interp1(va,qa,V_dch,'linear',NaN);dQ_dch=abs(diff(QR_d));end
        Ce_d=mean(abs(i_d))/Q_nom;
    end

    dc=dQ_chg(3:12); dd=dQ_dch(1:11);
    dc_r=dc/sum(dc,'omitnan'); dd_r=dd/sum(dd,'omitnan');
    X_f=[dc_r, dd_r, Ce_c, Ce_d];
    X_fs=(X_f-mu_tr)./sg_tr;
    X_imp=X_fs; X_imp(isnan(X_imp))=0;

    X_py=py.numpy.array(X_fs).reshape(int32(1),int32(-1));

    FR.(yr).BMS   = soh_bms;
    FR.(yr).LASSO = X_imp*M.LASSO.coef + M.LASSO.ic;
    FR.(yr).RF    = predict(M.RF, X_imp);
    FR.(yr).GBM   = predict(M.GBM, X_imp);
    FR.(yr).XGB   = double(M.XGB.predict(X_py));
    FR.(yr).LGB   = double(M.LGB.predict(X_py));
    fprintf('  BMS=%.1f%% | LASSO=%.1f%% RF=%.1f%% GBM=%.1f%% XGB=%.1f%% LGB=%.1f%%\n',...
        soh_bms,FR.(yr).LASSO,FR.(yr).RF,FR.(yr).GBM,FR.(yr).XGB,FR.(yr).LGB);
end

save(fullfile(maDir,'MA01_FieldResults.mat'),'FR');

%% ===== Part 6: Visualization =====
modelNames={'LASSO','RF','GBM','XGB','LGB'};
mColors={[0.22 0.43 0.74],[0.29 0.62 0.35],[0.85 0.40 0.30],[0.59 0.44 0.72],[0.95 0.70 0.25]};
mMarkers={'o','s','d','^','v'};

yr_avail={}; yr_xi=[]; bms_vals=[]; pred=zeros(4,5);
for k=1:4
    yr=years{k};
    if isfield(FR,yr)
        yr_avail{end+1}=yr; yr_xi(end+1)=year_x(k);
        bms_vals(end+1)=FR.(yr).BMS;
        m_idx=length(yr_avail);
        pred(m_idx,1)=FR.(yr).LASSO; pred(m_idx,2)=FR.(yr).RF;
        pred(m_idx,3)=FR.(yr).GBM; pred(m_idx,4)=FR.(yr).XGB; pred(m_idx,5)=FR.(yr).LGB;
    end
end
n=length(yr_avail);
pred=pred(1:n,:);

fig=figure('Position',[50 50 1400 300],'Color','w');
sgtitle('Field SOH Estimation per Model (MA-01, With Masking Augmentation)','FontSize',12,'FontWeight','bold');

for m=1:5
    ax=subplot(1,5,m);
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    set(ax,'GridColor',[0.88 0.88 0.88],'GridAlpha',1,'TickDir','out','FontSize',9);

    fill([yr_xi fliplr(yr_xi)],[bms_vals-1 fliplr(bms_vals+1)],...
        [0.8 0.8 0.8],'FaceAlpha',0.3,'EdgeColor','none');
    plot(ax,yr_xi,bms_vals,'k--o','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','k');
    plot(ax,yr_xi,pred(:,m),'-','Color',mColors{m},'LineWidth',2,...
        'Marker',mMarkers{m},'MarkerSize',8,'MarkerFaceColor',mColors{m});

    err=pred(:,m)-bms_vals';
    rmse_f=sqrt(mean(err.^2));
    mae_f=mean(abs(err));
    maxe_f=max(abs(err));

    for i=1:n
        text(ax,yr_xi(i),pred(i,m)+0.9,sprintf('%.1f%%',pred(i,m)),...
            'HorizontalAlignment','center','FontSize',11,'FontWeight','bold','Color',mColors{m});
    end

    if m==1
        txt_pos=[0.45 0.90];
    else
        txt_pos=[0.05 0.16];
    end
    text(ax,txt_pos(1),txt_pos(2),sprintf('RMSE=%.2f%%\nMAE =%.2f%%\nMax =%.2f%%',rmse_f,mae_f,maxe_f),...
        'Units','normalized','FontSize',8,'Color',[0.2 0.2 0.2],...
        'BackgroundColor',[0.97 0.97 0.97],'EdgeColor',[0.8 0.8 0.8],...
        'HorizontalAlignment','left');

    title(ax,modelNames{m},'FontSize',11,'FontWeight','bold','Color',mColors{m});
    set(ax,'XTick',yr_xi,'XTickLabel',yr_avail);
    ylabel(ax,'SOH (%)','FontSize',9);

    if m==1
        text(ax,0.5,0.5,{'※ LASSO 결과','확인 필요'},'Units','normalized',...
            'HorizontalAlignment','center','FontSize',9,'Color',[0.7 0.1 0.1],...
            'FontWeight','bold','BackgroundColor',[1 0.95 0.95],'EdgeColor',[0.8 0.3 0.3]);
        legend(ax,{'±1% band','BMS SOH','LASSO'},'FontSize',7,'Location','southeast','Box','off');
    else
        ylim(ax,[90 105]);
        legend(ax,{'±1% band','BMS SOH',modelNames{m}},'FontSize',7,'Location','southeast','Box','off');
    end
end

saveas(fig,fullfile(visDir,'MA01_Field_PerModel.png'));
fprintf('\nSaved: MA01_Field_PerModel.png\n');
fprintf('Done.\n');

%% Helpers
function segs=local_segs(mask)
    d=[0;diff(mask(:))]; starts=find(d==1); ends=find(d==-1)-1;
    if mask(1),starts=[1;starts];end; if mask(end),ends=[ends;length(mask)];end
    if isempty(starts)||isempty(ends),segs=[];return;end; segs=[starts,ends];
end
function v=getfield_safe(S,fname)
    if isfield(S,fname),v=S.(fname)(:); else,v=nan; end
end
