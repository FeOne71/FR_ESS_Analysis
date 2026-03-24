% RPT_Field_Estimation_RF_v2.m  — dQ만 사용한 필드 SOH 추정
% [v2] 피처: dQ_chg_S1~5 + dQ_dch_S1~5 = 10개
%       라벨: SOH만
%       Peak, Energy, C_eff 피처 없음
%       C_eff는 정규화 그룹 결정에만 사용
clear; clc; close all; warning off;

%% Section 1: 모델 로드
fprintf('=== Section 1: Loading Trained RF v2 Model ===\n');
baseDir = fullfile('D:','JCW','Projects','KEPCO_ESS_Local','ExperimentalData','FeatureEngineering','Lab_RPT_Analysis','ver02');
rfDir   = fullfile(baseDir, 'ML_RandomForest_v2');
load(fullfile(rfDir, 'Result_RF_v2.mat'), 'Results_RF');
load(fullfile('D:','JCW','Projects','KEPCO_ESS_Local','ExperimentalData','FeatureEngineering','MasterRulers_v3.mat'), 'MasterRulers');

fns = fieldnames(MasterRulers);
VR_chg = MasterRulers.(fns{1}).V_bounds_chg;
VR_dch = MasterRulers.(fns{1}).V_bounds_dch;

NormStats      = Results_RF.NormStats;
crate_groups   = NormStats.crate_groups;
crate_vals_num = NormStats.crate_vals_num;
chg_idx = NormStats.chg_idx;
dch_idx = NormStats.dch_idx;
final_rf = Results_RF.SOH.Model;

%% Section 2: 필드 데이터 로드
fprintf('\n=== Section 2: Loading Field Data ===\n');
rawDir = fullfile('D:','JCW','Projects','KEPCO_ESS_Local','Rack_raw2mat');
dataFiles = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021','old', datetime(2021,6,3);
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023','new', datetime(2023,10,16);
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024_Auto','new', datetime(2024,9,9);
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024_Manual','new', datetime(2024,9,9);
    fullfuie(rawDir,'New','2025','202507','Raw_20250711.mat'), 'Y2025','new', datetime(2025,7,11);
};
min_charge_secs  = [600, 300, 300, 300];
min_discharge_secs = [300, 150, 300, 300];
Np = 2; C_cell_Ah = 64; thr_A = C_cell_Ah * 0.05;
ma_window = 60;

Results_Ev = struct();

for k = 1:size(dataFiles, 1)
    fpath = dataFiles{k,1}; yr = dataFiles{k,2}; dtype = dataFiles{k,3}; base_date = dataFiles{k,4};
    if ~exist(fpath,'file'), continue; end
    fprintf('  Processing %s...\n', yr);
    S = load(fpath);
    
    % --- 데이터 파싱 (old/new 포맷) ---
    if strcmp(dtype,'old')
        D = S.Raw.Rack01;
        if isfield(D,'Time'), t = datetime(D.Time);
        elseif isduration(D.Date_Time), t = base_date + D.Date_Time;
        else, t = datetime(D.Date_Time); end
        I_rack = D.DCCurrent_A(:); V_avg = D.AverageCV_V(:);
    else
        D = S.Raw;
        if isduration(D.Date_Time), t = base_date + D.Date_Time;
        else, t = datetime(D.Date_Time); end
        I_rack = D.DCCurrent(:); V_avg = D.CVavg(:);
    end
    tsec = seconds(t - t(1)); I_cell = I_rack / Np;
    
    % --- 세그먼트 탐색 ---
    isChg = I_cell > thr_A; isDch = I_cell < -thr_A;
    chgSegs = local_find_segments(isChg);
    dchSegs = local_find_segments(isDch);
    if ~isempty(chgSegs), dur=chgSegs(:,2)-chgSegs(:,1)+1; chgSegs=chgSegs(dur>=min_charge_secs(k),:); end
    if ~isempty(dchSegs), dur=dchSegs(:,2)-dchSegs(:,1)+1; dchSegs=dchSegs(dur>=min_discharge_secs(k),:); end
    
    % Y2024_Manual 수동 세그먼트
    if strcmp(yr, 'Y2024_Manual')
        tcs=datetime(2024,9,9,12,55,0); tce=datetime(2024,9,9,13,13,0);
        ci=find(abs(t(chgSegs(:,1))-tcs)<minutes(5)&abs(t(chgSegs(:,2))-tce)<minutes(5),1);
        if isempty(ci), chgSegs=[]; else, chgSegs=chgSegs(ci,:); end
        tds=datetime(2024,9,9,14,14,0); tde=datetime(2024,9,9,14,25,49);
        di=find(abs(t(dchSegs(:,1))-tds)<minutes(1),1);
        if ~isempty(di)
            if t(dchSegs(di,2))>tde, il=find(t<=tde,1,'last'); if ~isempty(il)&&il>=dchSegs(di,1), dchSegs(di,2)=il; end; end
            dchSegs=dchSegs(di,:);
        else, dchSegs=[]; end
    end
    if isempty(chgSegs)||isempty(dchSegs), fprintf('  %s: Skipped.\\n',yr); continue; end
    
    [~,bi]=max(chgSegs(:,2)-chgSegs(:,1)); cs=chgSegs(bi,1); ce=chgSegs(bi,2);
    [~,bi]=max(dchSegs(:,2)-dchSegs(:,1)); ds=dchSegs(bi,1); de=dchSegs(bi,2);
    
    % --- V-Q 곡선 전처리 ---
    t_c=tsec(cs:ce); v_c=V_avg(cs:ce); i_c=I_cell(cs:ce);
    t_d=tsec(ds:de); v_d=V_avg(ds:de); i_d=I_cell(ds:de);
    i_c_sm=movmean(abs(i_c),ma_window); q_c=cumtrapz(t_c,i_c_sm)/3600;
    i_d_sm=movmean(abs(i_d),ma_window); q_d=cumtrapz(t_d,i_d_sm)/3600;
    v_c_sm=movmean(v_c,ma_window); v_d_sm=movmean(v_d,ma_window);
    q_c_sm=movmean(q_c,ma_window); q_d_sm=movmean(q_d,ma_window);
    [v_c_u,ui]=unique(v_c_sm,'stable'); q_c_u=q_c_sm(ui);
    [v_d_u,ui]=unique(v_d_sm,'stable'); q_d_u=q_d_sm(ui);
    
    % monotonicity 보장
    mc=true(size(v_c_u)); for ii=2:length(v_c_u), if v_c_u(ii)<=v_c_u(ii-1), mc(ii)=false; end; end
    v_c_u=v_c_u(mc); q_c_u=q_c_u(mc);
    md=true(size(v_d_u)); for ii=2:length(v_d_u), if v_d_u(ii)>=v_d_u(ii-1), md(ii)=false; end; end
    v_d_u=v_d_u(md); q_d_u=q_d_u(md);
    
    % --- dQ 세그먼트 추출 (MasterRuler 경계 기준) ---
    num_segs = length(VR_chg) - 1;
    V_min_c=min(v_c); V_max_c=max(v_c); V_min_d=min(v_d); V_max_d=max(v_d);
    valid_chg=false(1,num_segs); valid_dch=false(1,num_segs);
    for s=1:num_segs
        valid_chg(s)=(VR_chg(s)>=V_min_c-0.02)&&(VR_chg(s+1)<=V_max_c+0.02);
        valid_dch(s)=(VR_dch(s)<=V_max_d+0.02)&&(VR_dch(s+1)>=V_min_d-0.02);
    end
    
    dQ_chg=nan(1,num_segs); dQ_dch=nan(1,num_segs);
    if length(v_c_u)>1
        vg=(VR_chg(1)-0.05):0.001:(VR_chg(end)+0.05);
        Qi=interp1(v_c_u,q_c_u,vg,'linear','extrap'); Qs=movmean(Qi,30);
        QR=nan(1,length(VR_chg));
        for b=1:length(VR_chg),[~,mi]=min(abs(vg-VR_chg(b)));QR(b)=Qs(mi);end
        dQ_chg_raw=abs(diff(QR)); dQ_chg_raw(~valid_chg)=NaN;
        dQ_chg=segment_fallback(dQ_chg_raw, valid_chg);
    end
    if length(v_d_u)>1
        vg=(VR_dch(1)+0.05):-0.001:(VR_dch(end)-0.05);
        Qi=interp1(v_d_u,q_d_u,vg,'linear','extrap'); Qs=movmean(Qi,30);
        QR=nan(1,length(VR_dch));
        for b=1:length(VR_dch),[~,mi]=min(abs(vg-VR_dch(b)));QR(b)=Qs(mi);end
        dQ_dch_raw=abs(diff(QR)); dQ_dch_raw(~valid_dch)=NaN;
        dQ_dch=segment_fallback(dQ_dch_raw, valid_dch);
    end
    
    C_eff_chg = mean(abs(I_cell(cs:ce))) / C_cell_Ah;
    C_eff_dch = mean(abs(I_cell(ds:de))) / C_cell_Ah;
    
    % --- Split Normalization (dQ만, C_eff로 그룹 결정) ---
    X_raw = [dQ_chg, dQ_dch];  % 10 features
    if any(isnan(X_raw)), X_raw(isnan(X_raw))=0; end
    
    % 충전: C_eff_chg 기준 보간
    [~,sc]=sort(abs(crate_vals_num-C_eff_chg));
    cr1c=crate_groups{sc(1)}; cr2c=crate_groups{sc(2)};
    d1c=abs(C_eff_chg-crate_vals_num(sc(1))); d2c=abs(C_eff_chg-crate_vals_num(sc(2)));
    if (d1c+d2c)<eps, w1c=1;w2c=0; else, w1c=d2c/(d1c+d2c); w2c=d1c/(d1c+d2c); end
    mu_c  = w1c*NormStats.(cr1c).mu_chg  + w2c*NormStats.(cr2c).mu_chg;
    sig_c = w1c*NormStats.(cr1c).sigma_chg + w2c*NormStats.(cr2c).sigma_chg;
    
    % 방전: C_eff_dch 기준 보간
    [~,sd]=sort(abs(crate_vals_num-C_eff_dch));
    cr1d=crate_groups{sd(1)}; cr2d=crate_groups{sd(2)};
    d1d=abs(C_eff_dch-crate_vals_num(sd(1))); d2d=abs(C_eff_dch-crate_vals_num(sd(2)));
    if (d1d+d2d)<eps, w1d=1;w2d=0; else, w1d=d2d/(d1d+d2d); w2d=d1d/(d1d+d2d); end
    mu_d  = w1d*NormStats.(cr1d).mu_dch  + w2d*NormStats.(cr2d).mu_dch;
    sig_d = w1d*NormStats.(cr1d).sigma_dch + w2d*NormStats.(cr2d).sigma_dch;
    
    X_norm = X_raw;
    X_norm(chg_idx) = (X_raw(chg_idx) - mu_c) ./ (sig_c + eps);
    X_norm(dch_idx) = (X_raw(dch_idx) - mu_d) ./ (sig_d + eps);
    
    soh_est = predict(final_rf, X_norm);
    Results_Ev.(yr).SOH  = soh_est;
    Results_Ev.(yr).Ceff = struct('chg', C_eff_chg, 'dch', C_eff_dch);
    fprintf('  %s | SOH=%.2f%% | C_eff_chg=%.3fC→%s(%.0f%%) C_eff_dch=%.3fC→%s(%.0f%%)\n', ...
        yr, soh_est, C_eff_chg, cr1c, w1c*100, C_eff_dch, cr1d, w1d*100);
end

%% Section 3: Trajectory 시각화
fprintf('\n=== Section 3: Trajectory ===\n');
yr_list = fieldnames(Results_Ev); n_yrs = length(yr_list);
yr_nums = zeros(n_yrs,1);
for i=1:n_yrs
    tmp=regexp(yr_list{i},'\d{4}','match'); yr_nums(i)=str2double(tmp{1});
    if contains(yr_list{i},'Manual'), yr_nums(i)=yr_nums(i)+0.4; end
end

fig = figure('Name','RF_v2 Field SOH','Position',[100 100 700 450],'Visible','off');
hold on; grid on; box on;
soh_vals = arrayfun(@(i) Results_Ev.(yr_list{i}).SOH, 1:n_yrs);
plot(yr_nums, soh_vals, '-o', 'Color',[0.2 0.4 0.85], 'LineWidth',2, 'MarkerSize',9, 'MarkerFaceColor',[0.2 0.4 0.85]);
for i=1:n_yrs
    text(yr_nums(i), soh_vals(i)+0.3, sprintf('%.1f%%',soh_vals(i)), 'HorizontalAlignment','center','FontWeight','bold','FontSize',10);
end
x_labels = cellfun(@(s) strrep(strrep(s,'Y',''),'_',' '), yr_list, 'UniformOutput', false);
set(gca,'XTick',yr_nums,'XTickLabel',x_labels,'FontSize',11); xtickangle(30);
xlabel('Year','FontSize',12); ylabel('SOH (%)','FontSize',12);
title('RF v2 (dQ only) — Field SOH Trajectory','FontSize',13,'FontWeight','bold');
xlim([min(yr_nums)-0.5 max(yr_nums)+0.5]);
saveas(fig, fullfile(rfDir, 'RF_v2_Field_Trajectory.fig'));
close(fig);
fprintf('Figure saved.\n');

%% Helper Functions
function dQ=segment_fallback(dQ_raw,valid)
    dQ=dQ_raw; vi=find(valid);
    if isempty(vi),return;end
    for s=1:length(dQ_raw)
        if valid(s),continue;end
        lv=vi(vi<s); rv=vi(vi>s);
        if ~isempty(lv)&&~isempty(rv), l=lv(end);r=rv(1);dQ(s)=dQ_raw(l)+(dQ_raw(r)-dQ_raw(l))*(s-l)/(r-l);
        elseif ~isempty(lv), dQ(s)=dQ_raw(lv(end)); elseif ~isempty(rv), dQ(s)=dQ_raw(rv(1)); end
    end
end

function segs=local_find_segments(mask)
    segs=[]; i=1; n=length(mask);
    while i<=n; if mask(i),j=i;while j<n&&mask(j+1),j=j+1;end;segs=[segs;i,j];i=j+1;else,i=i+1;end;end
end
