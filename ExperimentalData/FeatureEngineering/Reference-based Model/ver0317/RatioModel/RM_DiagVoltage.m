% RM_DiagVoltage.m — 연도별 선택된 충전 세그먼트 전압 범위 확인
clear; clc;
projDir='C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir=fullfile(projDir,'ExperimentalData','FeatureEngineering','Lab_RPT_Analysis','ver0317');
rawDir=fullfile(projDir,'Rack_raw2mat');
mr=load(fullfile(verDir,'MasterRulers_PseudoOCV.mat'));
V_chg=mr.MasterRuler_ver0317.V_bounds_chg;
fprintf('V_chg segment boundaries:\n');
for s=1:12
    fprintf('  Seg%02d: %.3f ~ %.3f V\n',s,V_chg(s),V_chg(s+1));
end
fprintf('\n');

Q_nom=64; Np=2;
dataFiles={
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'),'Y2021','old',datetime(2021,6,3),[600 300];
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'),'Y2023','new',datetime(2023,10,16),[300 150];
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'),'Y2024','new',datetime(2024,9,9),[300 300];
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'),'Y2025','new',datetime(2025,7,11),[300 150];
};

for k=1:4
    fpath=dataFiles{k,1}; yr=dataFiles{k,2}; dtype=dataFiles{k,3};
    bd=dataFiles{k,4}; mdur=dataFiles{k,5};
    if ~exist(fpath,'file'), continue; end
    S=load(fpath);
    if strcmp(dtype,'old')
        D=S.Raw.Rack01;
        fns=fieldnames(D);
        t_fld=fns(contains(fns,'Time','IgnoreCase',true)|contains(fns,'Date','IgnoreCase',true));
        tf=D.(t_fld{1});
        if isduration(tf), t=bd+tf; else, t=datetime(tf); end
        I_rack=D.DCCurrent_A(:); V_avg=D.AverageCV_V(:);
    else
        D=S.Raw;
        if isduration(D.Date_Time), t=bd+D.Date_Time; else, t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent(:); V_avg=D.CVavg(:);
    end
    tsec=seconds(t-t(1)); I_cell=I_rack/Np; thr=Q_nom*0.05/Np;
    d2=[0;diff(double(I_cell>thr))];
    starts=find(d2==1); ends=find(d2==-1)-1;
    if I_cell(1)>thr, starts=[1;starts]; end
    if I_cell(end)>thr, ends=[ends;length(I_cell)]; end
    if isempty(starts)||isempty(ends), fprintf('%s: no chg segs\n',yr); continue; end
    chgS=[starts,ends]; dur=chgS(:,2)-chgS(:,1)+1;
    chgS=chgS(dur>=mdur(1),:);
    if isempty(chgS), fprintf('%s: no long enough chg seg\n',yr); continue; end
    [~,bi]=max(chgS(:,2)-chgS(:,1));
    v_c=V_avg(chgS(bi,1):chgS(bi,2));
    v_sm=smoothdata(v_c,'movmean',30);
    fprintf('%s: V_raw=[%.3f~%.3f]V  V_smooth=[%.3f~%.3f]V  dur=%ds\n',...
        yr,min(v_c),max(v_c),min(v_sm),max(v_sm),length(v_c));
    % Check which V_chg segments are covered
    covered=false(1,12);
    for s=1:12
        in_seg = v_sm >= V_chg(s) & v_sm <= V_chg(s+1);
        covered(s) = any(in_seg);
    end
    cov_str=''; for s=1:12, if covered(s), cov_str=[cov_str sprintf('Seg%02d ',s)]; end; end
    fprintf('  Covered segs: %s\n',cov_str);
    fprintf('  NaN segs: ');
    for s=1:12, if ~covered(s), fprintf('Seg%02d ',s); end; end
    fprintf('\n\n');
end
