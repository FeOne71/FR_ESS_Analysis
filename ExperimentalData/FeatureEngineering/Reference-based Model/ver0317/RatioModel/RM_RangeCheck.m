% RM_RangeCheck.m — Check if field features are within Lab training range
clear; clc;
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir,'ExperimentalData','FeatureEngineering','Lab_RPT_Analysis','ver0317');
rawDir  = fullfile(projDir,'Rack_raw2mat');
Q_nom=64; Np=2; num_segs=12;

d=load(fullfile(verDir,'FeatureMatrix_ver0317.mat')); FM=d.FM;
mr=load(fullfile(verDir,'MasterRulers_PseudoOCV.mat'));
V_chg=mr.MasterRuler_ver0317.V_bounds_chg; V_dch=mr.MasterRuler_ver0317.V_bounds_dch;
chg_cols=arrayfun(@(i)sprintf('dQ_c_%02d',i),3:12,'UniformOutput',false);
dch_cols=arrayfun(@(i)sprintf('dQ_d_%02d',i),1:11,'UniformOutput',false);
dQ_c_r=FM{:,chg_cols}./sum(FM{:,chg_cols},2,'omitnan');
dQ_d_r=FM{:,dch_cols}./sum(FM{:,dch_cols},2,'omitnan');
X_lab=[dQ_c_r, dQ_d_r, FM.C_eff_chg, FM.C_eff_dch];
lab_min=min(X_lab,[],1,'omitnan'); lab_max=max(X_lab,[],1,'omitnan');

feat_names=[arrayfun(@(i)sprintf('dQr_c%02d',i),3:12,'UniformOutput',false),...
            arrayfun(@(i)sprintf('dQr_d%02d',i),1:11,'UniformOutput',false),...
            {'Ceff_chg','Ceff_dch'}];

dataFiles={
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'),'Y2021','old',datetime(2021,6,3),[600 300];
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'),'Y2023','new',datetime(2023,10,16),[300 150];
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'),'Y2024','new',datetime(2024,9,9),[300 300];
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'),'Y2025','new',datetime(2025,7,11),[300 150];
};

for k=1:size(dataFiles,1)
    fpath=dataFiles{k,1}; yr=dataFiles{k,2}; dtype=dataFiles{k,3};
    bd=dataFiles{k,4}; mdur=dataFiles{k,5};
    if ~exist(fpath,'file'), continue; end
    S=load(fpath);
    if strcmp(dtype,'old')
        D=S.Raw.Rack01;
        if isfield(D,'Time'),t=datetime(D.Time);
        elseif isduration(D.Date_Time),t=bd+D.Date_Time; else,t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent_A(:); V_avg=D.AverageCV_V(:);
    else
        D=S.Raw;
        if isduration(D.Date_Time),t=bd+D.Date_Time; else,t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent(:); V_avg=D.CVavg(:);
    end
    tsec=seconds(t-t(1)); I_cell=I_rack/Np; thr=Q_nom*0.05/Np;
    chgS=local_segs(I_cell>thr); if ~isempty(chgS),dur=chgS(:,2)-chgS(:,1)+1;chgS=chgS(dur>=mdur(1),:);end
    dchS=local_segs(I_cell<-thr); if ~isempty(dchS),dur=dchS(:,2)-dchS(:,1)+1;dchS=dchS(dur>=mdur(2),:);end
    if isempty(chgS), continue; end

    [~,bi]=max(chgS(:,2)-chgS(:,1));
    v_c=V_avg(chgS(bi,1):chgS(bi,2)); i_c=I_cell(chgS(bi,1):chgS(bi,2)); t_c=tsec(chgS(bi,1):chgS(bi,2));
    v_c=smoothdata(v_c,'movmean',30); q_c=cumtrapz(t_c,abs(i_c))/3600;
    vm=v_c(1);qm=q_c(1);
    for ii=2:length(v_c),if v_c(ii)>vm(end),vm(end+1)=v_c(ii);qm(end+1)=q_c(ii);end;end
    dQ_chg=nan(1,num_segs);
    if length(vm)>1,QR=interp1(vm,qm,V_chg,'linear',NaN);dQ_chg=abs(diff(QR));end
    Ce_c=mean(abs(i_c))/Q_nom;

    dQ_dch=nan(1,num_segs); Ce_d=NaN;
    if ~isempty(dchS)
        [~,bi_d]=max(dchS(:,2)-dchS(:,1));
        v_d=V_avg(dchS(bi_d,1):dchS(bi_d,2)); i_d=I_cell(dchS(bi_d,1):dchS(bi_d,2)); t_d=tsec(dchS(bi_d,1):dchS(bi_d,2));
        v_d=smoothdata(v_d,'movmean',30); q_d=cumtrapz(t_d,abs(i_d))/3600;
        vd=v_d(1);qd=q_d(1);
        for ii=2:length(v_d),if v_d(ii)<vd(end),vd(end+1)=v_d(ii);qd(end+1)=q_d(ii);end;end
        va=flip(vd);qa=flip(qd);
        if length(va)>1,QR_d=interp1(va,qa,V_dch,'linear',NaN);dQ_dch=abs(diff(QR_d));end
        Ce_d=mean(abs(i_d))/Q_nom;
    end
    dc=dQ_chg(3:12); dd=dQ_dch(1:11);
    X_f=[dc/sum(dc,'omitnan'), dd/sum(dd,'omitnan'), Ce_c, Ce_d];

    fprintf('\n=== %s ===\n',yr);
    n_out=0;
    for f=1:length(X_f)
        if isnan(X_f(f)), continue; end
        if X_f(f)<lab_min(f) || X_f(f)>lab_max(f)
            n_out=n_out+1;
            fprintf('  OUT  %-12s  field=%.4f  lab=[%.4f, %.4f]\n',...
                feat_names{f},X_f(f),lab_min(f),lab_max(f));
        end
    end
    fprintf('  Valid features: %d/%d | Out-of-range: %d\n',...
        sum(~isnan(X_f)),length(feat_names),n_out);
end

function segs=local_segs(mask)
    d=[0;diff(mask(:))]; starts=find(d==1); ends=find(d==-1)-1;
    if mask(1),starts=[1;starts];end; if mask(end),ends=[ends;length(mask)];end
    if isempty(starts)||isempty(ends),segs=[];return;end; segs=[starts,ends];
end
