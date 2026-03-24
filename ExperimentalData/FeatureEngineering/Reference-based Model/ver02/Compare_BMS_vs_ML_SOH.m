% Compare_BMS_vs_ML_SOH.m
% BMS SOH와 RF/SVM 모델 추정 SOH 비교 시각화
% - BMS SOH: 필드 데이터의 SOHPct(old) / SOH_BMS(new) 필드에서 추출
% - ML SOH: Result_RandomForest.mat, Result_SVM.mat의 필드 추정 결과

clear; clc; close all;

baseDir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
    'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
rfDir  = fullfile(baseDir, 'ML_RandomForest');
svmDir = fullfile(baseDir, 'ML_SVM');
rawDir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'Rack_raw2mat');

%% 1. ML 추정 결과 로드
fprintf('=== Loading ML Results ===\n');
load(fullfile(rfDir,  'Result_RandomForest.mat'), 'Results_RF');
load(fullfile(svmDir, 'Result_SVM.mat'),          'Results_SVM');

% 필드 추정 결과(Results_Ev)는 스크립트 내부 변수이므로, 
% 각 모델의 FinalModel과 NormStats를 사용해 다시 추정 수행
% → RPT_Field_Estimation_RF 와 동일한 피처 추출 로직을 내장

%% 2. 필드 데이터 정의 (Y2025 제외)
dataFiles = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021', 'old',  datetime(2021,6,3);
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023', 'new',  datetime(2023,10,16);
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024_Auto',   'new',  datetime(2024,9,9);
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024_Manual', 'new',  datetime(2024,9,9);
};

min_charge_secs    = [600, 300, 300, 300];
min_discharge_secs = [300, 150, 300, 300];

Np = 2; C_cell_Ah = 64; thr_A = C_cell_Ah * 0.05;

% MasterRuler 로드
path_mr = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
    'FeatureEngineering', 'MasterRulers_v3.mat');
load(path_mr, 'MasterRulers');
fns = fieldnames(MasterRulers);
VR_chg = MasterRulers.(fns{1}).V_bounds_chg;
VR_dch = MasterRulers.(fns{1}).V_bounds_dch;

fallback_priority = [3, 4, 2, 5, 1];
ma_window = 60; Q_0 = 64;

%% 3. BMS SOH 추출 + ML SOH 추정
fprintf('\n=== Extracting BMS SOH & Computing ML SOH ===\n');

crate_groups   = Results_RF.NormStats.crate_groups;
crate_vals_num = Results_RF.NormStats.crate_vals_num;
chg_idx  = [1:5, 11, 13, 15];
dch_idx  = [6:10, 12, 14, 16, 17];
label_names = Results_RF.label_names;

ResultsComp = struct();

for k = 1:size(dataFiles, 1)
    fpath      = dataFiles{k,1};
    yr         = dataFiles{k,2};
    dtype      = dataFiles{k,3};
    base_date  = dataFiles{k,4};
    if ~exist(fpath, 'file'), continue; end
    fprintf('  Processing %s...\n', yr);

    S = load(fpath);

    %% -- BMS SOH 추출 --
    if strcmp(dtype, 'old')
        D = S.Raw.Rack01;
        if isfield(D,'Time'), t = datetime(D.Time);
        elseif isduration(D.Date_Time), t = base_date + D.Date_Time;
        else, t = datetime(D.Date_Time); end
        I_rack = D.DCCurrent_A(:); V_avg = D.AverageCV_V(:);
        if isfield(D,'AverageMT_degC'), T_avg = D.AverageMT_degC(:);
        else, T_avg = nan(size(I_rack)); end
        if isfield(D,'SOHPct'), raw_soh_bms = D.SOHPct(:);
        else, raw_soh_bms = nan(size(I_rack)); end
    else
        D = S.Raw;
        if isduration(D.Date_Time), t = base_date + D.Date_Time;
        else, t = datetime(D.Date_Time); end
        I_rack = D.DCCurrent(:); V_avg = D.CVavg(:);
        if isfield(D,'MTavg'), T_avg = D.MTavg(:);
        else, T_avg = nan(size(I_rack)); end
        if isfield(D,'SOH_BMS'), raw_soh_bms = D.SOH_BMS(:);
        else, raw_soh_bms = nan(size(I_rack)); end
    end

    tsec   = seconds(t - t(1));
    I_cell = I_rack / Np;

    % BMS SOH: 0은 센서 이상값 → NaN 처리 후 하루 대표값 = 유효값의 최빈값/중앙값
    valid_soh = raw_soh_bms;
    valid_soh(valid_soh <= 0) = NaN;
    soh_bms_day = round(median(valid_soh, 'omitnan'));  % 하루 대표값

    %% -- 충방전 세그먼트 탐색 --
    isChg = I_cell > thr_A; isDch = I_cell < -thr_A;
    chgSegs = local_find_segments(isChg);
    dchSegs = local_find_segments(isDch);

    if ~isempty(chgSegs)
        dur_c = chgSegs(:,2)-chgSegs(:,1)+1;
        chgSegs = chgSegs(dur_c >= ceil(min_charge_secs(k)/1), :);
    end
    if ~isempty(dchSegs)
        dur_d = dchSegs(:,2)-dchSegs(:,1)+1;
        dchSegs = dchSegs(dur_d >= ceil(min_discharge_secs(k)/1), :);
    end

    % Y2024_Manual: 특정 시간대 세그먼트 수동 지정
    if strcmp(yr, 'Y2024_Manual')
        tcs = datetime(2024,9,9,12,55,0); tce = datetime(2024,9,9,13,13,0);
        ci = find(abs(t(chgSegs(:,1))-tcs)<minutes(5) & abs(t(chgSegs(:,2))-tce)<minutes(5),1);
        if isempty(ci), chgSegs = []; else, chgSegs = chgSegs(ci,:); end
        tds = datetime(2024,9,9,14,14,0); tde = datetime(2024,9,9,14,25,49);
        di = find(abs(t(dchSegs(:,1))-tds)<minutes(1),1);
        if ~isempty(di)
            if t(dchSegs(di,2)) > tde
                il = find(t <= tde, 1, 'last');
                if ~isempty(il) && il >= dchSegs(di,1), dchSegs(di,2)=il; end
            end
            dchSegs = dchSegs(di,:);
        else, dchSegs = []; end
    end

    if isempty(chgSegs) || isempty(dchSegs), continue; end

    [~,bi] = max(chgSegs(:,2)-chgSegs(:,1)); chg_s=chgSegs(bi,1); chg_e=chgSegs(bi,2);
    [~,bi] = max(dchSegs(:,2)-dchSegs(:,1)); dch_s=dchSegs(bi,1); dch_e=dchSegs(bi,2);

    t_c=tsec(chg_s:chg_e); v_c=V_avg(chg_s:chg_e); i_c=I_cell(chg_s:chg_e);
    t_d=tsec(dch_s:dch_e); v_d=V_avg(dch_s:dch_e); i_d=I_cell(dch_s:dch_e);

    i_c_sm=movmean(abs(i_c),ma_window); i_d_sm=movmean(abs(i_d),ma_window);
    q_c=cumtrapz(t_c,i_c_sm)/3600; q_d=cumtrapz(t_d,i_d_sm)/3600;
    v_c_sm=movmean(v_c,ma_window); q_c_sm=movmean(q_c,ma_window);
    v_d_sm=movmean(v_d,ma_window); q_d_sm=movmean(q_d,ma_window);

    [v_c_u,ui]=unique(v_c_sm,'stable'); q_c_u=q_c_sm(ui);
    [v_d_u,ui]=unique(v_d_sm,'stable'); q_d_u=q_d_sm(ui);

    mono_c=true(size(v_c_u));
    for ii=2:length(v_c_u), if v_c_u(ii)<=v_c_u(ii-1), mono_c(ii)=false; end; end
    v_c_u=v_c_u(mono_c); q_c_u=q_c_u(mono_c);
    mono_d=true(size(v_d_u));
    for ii=2:length(v_d_u), if v_d_u(ii)>=v_d_u(ii-1), mono_d(ii)=false; end; end
    v_d_u=v_d_u(mono_d); q_d_u=q_d_u(mono_d);

    %% -- dQ 피처 추출 --
    num_segs=length(VR_chg)-1;
    V_min_c=min(v_c); V_max_c=max(v_c);
    V_min_d=min(v_d); V_max_d=max(v_d);
    valid_chg=false(1,num_segs); valid_dch=false(1,num_segs);
    for s=1:num_segs
        valid_chg(s)=(VR_chg(s)>=V_min_c-0.02)&&(VR_chg(s+1)<=V_max_c+0.02);
        valid_dch(s)=(VR_dch(s)<=V_max_d+0.02)&&(VR_dch(s+1)>=V_min_d-0.02);
    end

    dQ_chg_raw=nan(1,num_segs); dQ_dch_raw=nan(1,num_segs);
    if length(v_c_u)>1
        try
            vg=(VR_chg(1)-0.05):0.001:(VR_chg(end)+0.05);
            Qi=interp1(v_c_u,q_c_u,vg,'linear','extrap'); Qs=movmean(Qi,30);
            QR=nan(1,length(VR_chg));
            for b=1:length(VR_chg),[~,mi]=min(abs(vg-VR_chg(b)));QR(b)=Qs(mi);end
            dQ_chg_raw=abs(diff(QR)); dQ_chg_raw(~valid_chg)=NaN;
        catch; end
    end
    if length(v_d_u)>1
        try
            vg=(VR_dch(1)+0.05):-0.001:(VR_dch(end)-0.05);
            Qi=interp1(v_d_u,q_d_u,vg,'linear','extrap'); Qs=movmean(Qi,30);
            QR=nan(1,length(VR_dch));
            for b=1:length(VR_dch),[~,mi]=min(abs(vg-VR_dch(b)));QR(b)=Qs(mi);end
            dQ_dch_raw=abs(diff(QR)); dQ_dch_raw(~valid_dch)=NaN;
        catch; end
    end
    dQ_chg=segment_fallback(dQ_chg_raw,valid_chg);
    dQ_dch=segment_fallback(dQ_dch_raw,valid_dch);

    %% -- dQ/dV Peak 추출 --
    [PkH_chg,PkA_chg,PkPos_chg]=field_extract_peak(v_c_u,q_c_u);
    [PkH_dch,PkA_dch,PkPos_dch]=field_extract_peak(v_d_u,q_d_u);

    Qdch_seg = cumtrapz(abs(I_cell(dch_s:dch_e)))/3600;
    V_dch_seg = V_avg(dch_s:dch_e);
    idx_d = V_dch_seg >= min(VR_dch) & V_dch_seg <= max(VR_dch);
    if sum(idx_d) > 5
        Energy_dch = abs(trapz(Qdch_seg(idx_d), V_dch_seg(idx_d)));
    else, Energy_dch = NaN; end

    C_eff_chg = mean(abs(I_cell(chg_s:chg_e))) / Q_0;
    C_eff_dch = mean(abs(I_cell(dch_s:dch_e))) / Q_0;

    X_raw=[dQ_chg,dQ_dch,PkH_chg,PkH_dch,PkA_chg,PkA_dch,PkPos_chg,PkPos_dch,Energy_dch,C_eff_chg,C_eff_dch];
    if any(isnan(X_raw)), X_raw(isnan(X_raw))=0; end

    %% -- Split Normalization --
    NormS_rf = Results_RF.NormStats;
    [mu_c,si_c,mu_d,si_d] = interp_norm(C_eff_chg,C_eff_dch,NormS_rf,crate_groups,crate_vals_num,chg_idx,dch_idx);
    X_norm_rf=X_raw; X_norm_rf(chg_idx)=(X_raw(chg_idx)-mu_c)./(si_c+eps); X_norm_rf(dch_idx)=(X_raw(dch_idx)-mu_d)./(si_d+eps);

    NormS_sv = Results_SVM.NormStats;
    [mu_c2,si_c2,mu_d2,si_d2]=interp_norm(C_eff_chg,C_eff_dch,NormS_sv,crate_groups,crate_vals_num,chg_idx,dch_idx);
    X_norm_sv=X_raw; X_norm_sv(chg_idx)=(X_raw(chg_idx)-mu_c2)./(si_c2+eps); X_norm_sv(dch_idx)=(X_raw(dch_idx)-mu_d2)./(si_d2+eps);

    %% -- ML 추정 --
    soh_rf  = predict(Results_RF.Merged.SOH.Model,  X_norm_rf);
    soh_svm = predict(Results_SVM.SOH.FinalModel,   X_norm_sv);

    ResultsComp.(yr).SOH_BMS    = soh_bms_day;
    ResultsComp.(yr).SOH_RF     = soh_rf;
    ResultsComp.(yr).SOH_SVM    = soh_svm;
    ResultsComp.(yr).C_eff_chg  = C_eff_chg;
    ResultsComp.(yr).C_eff_dch  = C_eff_dch;

    fprintf('  %s | BMS=%.1f%%  RF=%.2f%%  SVM=%.2f%%\n', yr, soh_bms_day, soh_rf, soh_svm);
end

%% 4. 비교 시각화
fprintf('\n=== Generating Comparison Figure ===\n');

yr_list = fieldnames(ResultsComp);
n_yrs   = length(yr_list);

yr_nums = zeros(n_yrs,1);
for k=1:n_yrs
    tmp=regexp(yr_list{k},'\d{4}','match'); yr_nums(k)=str2double(tmp{1});
    if contains(yr_list{k},'Manual'), yr_nums(k)=yr_nums(k)+0.4; end
end

x_labels = cellfun(@(s) strrep(strrep(s,'Y',''),'_',' '), yr_list, 'UniformOutput', false);

soh_bms = arrayfun(@(k) ResultsComp.(yr_list{k}).SOH_BMS, 1:n_yrs);
soh_rf  = arrayfun(@(k) ResultsComp.(yr_list{k}).SOH_RF,  1:n_yrs);
soh_svm = arrayfun(@(k) ResultsComp.(yr_list{k}).SOH_SVM, 1:n_yrs);

fig = figure('Name', 'BMS vs ML SOH Comparison', 'Position', [100, 100, 1000, 500]);
sgtitle('KIMJ ESS – BMS SOH vs ML 추정 SOH 비교', 'FontSize', 15, 'FontWeight', 'bold');

ax = axes; hold on; grid on; box on;

% BMS
plot(yr_nums, soh_bms, 'k--s', 'LineWidth', 2.5, 'MarkerSize', 10, ...
    'MarkerFaceColor', 'k', 'DisplayName', 'BMS SOH');

% RF
plot(yr_nums, soh_rf, '-o', 'Color', [0.2 0.4 0.85], 'LineWidth', 2, ...
    'MarkerSize', 9, 'MarkerFaceColor', [0.2 0.4 0.85], 'DisplayName', 'RF 추정 SOH');

% SVM
plot(yr_nums, soh_svm, '-^', 'Color', [0.9 0.35 0.15], 'LineWidth', 2, ...
    'MarkerSize', 9, 'MarkerFaceColor', [0.9 0.35 0.15], 'DisplayName', 'SVM 추정 SOH');

% 데이터 라벨
for k=1:n_yrs
    text(yr_nums(k), soh_bms(k)+0.5,  sprintf('%.1f%%', soh_bms(k)),  'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color','k');
    text(yr_nums(k)-0.08, soh_rf(k)-0.8,  sprintf('%.1f%%',soh_rf(k)),  'HorizontalAlignment','center','FontSize',9,'Color',[0.2 0.4 0.85]);
    text(yr_nums(k)+0.08, soh_svm(k)-0.8, sprintf('%.1f%%',soh_svm(k)), 'HorizontalAlignment','center','FontSize',9,'Color',[0.9 0.35 0.15]);
end

set(gca, 'XTick', yr_nums, 'XTickLabel', x_labels, 'FontSize', 11);
xtickangle(30);
xlabel('측정 연도', 'FontSize', 12);
ylabel('SOH (%)', 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 11);
xlim([min(yr_nums)-0.5, max(yr_nums)+0.5]);

% 오차 출력
fprintf('\n=== BMS vs ML 오차 요약 ===\n');
fprintf('%-15s  %-12s  %-12s  %-12s\n', 'Year', 'BMS SOH', 'RF Err', 'SVM Err');
for k=1:n_yrs
    fprintf('%-15s  %-12.1f  %-+12.2f  %-+12.2f\n', yr_list{k}, soh_bms(k), soh_rf(k)-soh_bms(k), soh_svm(k)-soh_bms(k));
end

outDir = fullfile(baseDir, 'ML_RandomForest');
saveas(fig, fullfile(outDir, 'Compare_BMS_vs_ML_SOH.fig'));
fprintf('\nFigure saved: Compare_BMS_vs_ML_SOH.fig\n');

%% Helper Functions
function [mu_c, si_c, mu_d, si_d] = interp_norm(Cc, Cd, NormS, groups, vals, ci, di)
    [~,sc]=sort(abs(vals-Cc)); cr1=groups{sc(1)}; cr2=groups{sc(2)};
    d1=abs(Cc-vals(sc(1))); d2=abs(Cc-vals(sc(2)));
    if (d1+d2)<eps, w1=1;w2=0; else, w1=d2/(d1+d2);w2=d1/(d1+d2); end
    mu_c  = w1*NormS.(cr1).mu_chg    + w2*NormS.(cr2).mu_chg;
    si_c  = w1*NormS.(cr1).sigma_chg + w2*NormS.(cr2).sigma_chg;
    [~,sd]=sort(abs(vals-Cd)); dr1=groups{sd(1)}; dr2=groups{sd(2)};
    d1d=abs(Cd-vals(sd(1))); d2d=abs(Cd-vals(sd(2)));
    if (d1d+d2d)<eps, w1d=1;w2d=0; else, w1d=d2d/(d1d+d2d);w2d=d1d/(d1d+d2d); end
    mu_d  = w1d*NormS.(dr1).mu_dch    + w2d*NormS.(dr2).mu_dch;
    si_d  = w1d*NormS.(dr1).sigma_dch + w2d*NormS.(dr2).sigma_dch;
end

function [pk_h,pk_a,pk_p]=field_extract_peak(V_u,Q_u)
    pk_h=NaN;pk_a=NaN;pk_p=NaN;
    if length(V_u)<=5,return;end
    Vg=(ceil(min(V_u)*1000)/1000:0.001:floor(max(V_u)*1000)/1000)';
    Qg=interp1(V_u,Q_u,Vg,'linear','extrap');
    dQdV=abs(diff(Qg))/0.001; Vm=Vg(1:end-1)+0.0005;
    dQdV=movmean(dQdV,10); dQdV=max(dQdV,0);
    et=10; if length(dQdV)>2*et+1, dQdV=dQdV(et+1:end-et);Vm=Vm(et+1:end-et);end
    pk_a=abs(trapz(Vm,dQdV));
    [pks,locs]=findpeaks(dQdV,'SortStr','descend','NPeaks',1);
    if ~isempty(pks),pk_h=pks(1);pk_p=Vm(locs(1));end
end

function dQ=segment_fallback(dQ_raw,valid)
    dQ=dQ_raw; N=length(dQ_raw); vi=find(valid);
    if isempty(vi),return;end
    for s=1:N
        if valid(s),continue;end
        lv=vi(vi<s); rv=vi(vi>s);
        if ~isempty(lv)&&~isempty(rv)
            l=lv(end);r=rv(1);
            dQ(s)=dQ_raw(l)+(dQ_raw(r)-dQ_raw(l))*(s-l)/(r-l);
        elseif ~isempty(lv), dQ(s)=dQ_raw(lv(end));
        elseif ~isempty(rv), dQ(s)=dQ_raw(rv(1));
        end
    end
end

function segs=local_find_segments(mask)
    segs=[]; n=length(mask); i=1;
    while i<=n
        if mask(i), j=i;
            while j<n&&mask(j+1),j=j+1;end
            segs=[segs;i,j]; i=j+1;
        else, i=i+1;
        end
    end
end
