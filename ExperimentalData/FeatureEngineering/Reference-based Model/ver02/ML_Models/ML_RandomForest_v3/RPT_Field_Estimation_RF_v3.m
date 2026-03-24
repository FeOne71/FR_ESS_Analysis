% RPT_Field_Estimation_RF_v3.m — [3.4~4.0]V 11세그먼트, 부분 세그먼트 지원
%
% [v3 핵심] 필드 V-Q 곡선에서 11개 dQ 세그먼트 추출 시:
%   - 실제 전압 범위에 해당하는 세그먼트만 유효 (나머지 NaN)
%   - RF 모델이 Surrogate splits로 NaN 피처 처리 → 몇 개 세그먼트만 있어도 예측 가능
%   - Fallback 선형보간 사용하지 않음 (NaN 그대로 전달)
clear; clc; close all; warning off;

%% Section 1: 모델 로드
fprintf('=== Section 1: Loading Trained RF v3 Model ===\n');
baseDir = fullfile('D:','JCW','Projects','KEPCO_ESS_Local','ExperimentalData', ...
    'FeatureEngineering','Lab_RPT_Analysis','ver02');
rfDir = fullfile(baseDir, 'ML_RandomForest_v3');
load(fullfile(rfDir, 'Result_RF_v3.mat'), 'Results_RF', 'V_bounds');

NormStats      = Results_RF.NormStats;
crate_groups   = NormStats.crate_groups;
crate_vals_num = NormStats.crate_vals_num;
chg_idx        = NormStats.chg_idx;
dch_idx        = NormStats.dch_idx;
final_rf       = Results_RF.SOH.Model;
num_segments   = length(V_bounds) - 1;

fprintf('  V_bounds: [%.2f ~ %.2f] V, %d segments\n', V_bounds(1), V_bounds(end), num_segments);

%% Section 2: 필드 데이터 로드
fprintf('\n=== Section 2: Loading Field Data ===\n');
rawDir = fullfile('D:','JCW','Projects','KEPCO_ESS_Local','Rack_raw2mat');
dataFiles = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021','old', datetime(2021,6,3);
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023','new', datetime(2023,10,16);
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024_Auto','new', datetime(2024,9,9);
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024_Manual','new', datetime(2024,9,9);
};
min_charge_secs  = [600, 300, 300, 300];
min_discharge_secs = [300, 150, 300, 300];
Np = 2; C_cell_Ah = 64; thr_A = C_cell_Ah * 0.05;
ma_window = 60;

Results_Ev = struct();

for k = 1:size(dataFiles, 1)
    fpath=dataFiles{k,1}; yr=dataFiles{k,2}; dtype=dataFiles{k,3}; base_date=dataFiles{k,4};
    if ~exist(fpath,'file'), continue; end
    fprintf('  Processing %s...\n', yr);
    S = load(fpath);

    % 데이터 파싱
    if strcmp(dtype,'old')
        D=S.Raw.Rack01;
        if isfield(D,'Time'),t=datetime(D.Time);
        elseif isduration(D.Date_Time),t=base_date+D.Date_Time;
        else,t=datetime(D.Date_Time);end
        I_rack=D.DCCurrent_A(:); V_avg=D.AverageCV_V(:);
    else
        D=S.Raw;
        if isduration(D.Date_Time),t=base_date+D.Date_Time;
        else,t=datetime(D.Date_Time);end
        I_rack=D.DCCurrent(:); V_avg=D.CVavg(:);
    end
    tsec=seconds(t-t(1)); I_cell=I_rack/Np;

    % 세그먼트 탐색
    isChg=I_cell>thr_A; isDch=I_cell<-thr_A;
    chgSegs=local_find_segments(isChg);
    dchSegs=local_find_segments(isDch);
    if ~isempty(chgSegs),dur=chgSegs(:,2)-chgSegs(:,1)+1;chgSegs=chgSegs(dur>=min_charge_secs(k),:);end
    if ~isempty(dchSegs),dur=dchSegs(:,2)-dchSegs(:,1)+1;dchSegs=dchSegs(dur>=min_discharge_secs(k),:);end

    % Y2024_Manual 수동 세그먼트
    if strcmp(yr,'Y2024_Manual')
        tcs=datetime(2024,9,9,12,55,0);tce=datetime(2024,9,9,13,13,0);
        ci=find(abs(t(chgSegs(:,1))-tcs)<minutes(5)&abs(t(chgSegs(:,2))-tce)<minutes(5),1);
        if isempty(ci),chgSegs=[];else,chgSegs=chgSegs(ci,:);end
        tds=datetime(2024,9,9,14,14,0);tde=datetime(2024,9,9,14,25,49);
        di=find(abs(t(dchSegs(:,1))-tds)<minutes(1),1);
        if ~isempty(di)
            if t(dchSegs(di,2))>tde,il=find(t<=tde,1,'last');if~isempty(il)&&il>=dchSegs(di,1),dchSegs(di,2)=il;end;end
            dchSegs=dchSegs(di,:);
        else,dchSegs=[];end
    end
    if isempty(chgSegs)||isempty(dchSegs),fprintf('  %s: Skipped.\n',yr);continue;end

    [~,bi]=max(chgSegs(:,2)-chgSegs(:,1));cs=chgSegs(bi,1);ce=chgSegs(bi,2);
    [~,bi]=max(dchSegs(:,2)-dchSegs(:,1));ds=dchSegs(bi,1);de=dchSegs(bi,2);

    % V-Q 전처리
    t_c=tsec(cs:ce);v_c=V_avg(cs:ce);i_c=I_cell(cs:ce);
    t_d=tsec(ds:de);v_d=V_avg(ds:de);i_d=I_cell(ds:de);
    i_c_sm=movmean(abs(i_c),ma_window); q_c=cumtrapz(t_c,i_c_sm)/3600;
    i_d_sm=movmean(abs(i_d),ma_window); q_d=cumtrapz(t_d,i_d_sm)/3600;
    v_c_sm=movmean(v_c,ma_window); v_d_sm=movmean(v_d,ma_window);
    q_c_sm=movmean(q_c,ma_window); q_d_sm=movmean(q_d,ma_window);
    [v_c_u,ui]=unique(v_c_sm,'stable'); q_c_u=q_c_sm(ui);
    [v_d_u,ui]=unique(v_d_sm,'stable'); q_d_u=q_d_sm(ui);
    mc=true(size(v_c_u));for ii=2:length(v_c_u),if v_c_u(ii)<=v_c_u(ii-1),mc(ii)=false;end;end
    v_c_u=v_c_u(mc);q_c_u=q_c_u(mc);
    md=true(size(v_d_u));for ii=2:length(v_d_u),if v_d_u(ii)>=v_d_u(ii-1),md(ii)=false;end;end
    v_d_u=v_d_u(md);q_d_u=q_d_u(md);

    V_min_c=min(v_c);V_max_c=max(v_c); V_min_d=min(v_d);V_max_d=max(v_d);

    % === 11 dQ 세그먼트 추출 (NaN 허용, fallback 없음) ===
    dQ_chg = nan(1, num_segments);
    dQ_dch = nan(1, num_segments);

    if length(v_c_u) > 1
        vg = (V_bounds(1)-0.02):0.001:(V_bounds(end)+0.02);
        Qi = interp1(v_c_u, q_c_u, vg, 'linear', 'extrap');
        Qs = movmean(Qi, 30);
        for s = 1:num_segments
            vlo = V_bounds(s); vhi = V_bounds(s+1);
            % 양쪽 경계 모두 실제 V 범위 내에 있을 때만 유효
            if vlo >= V_min_c - 0.02 && vhi <= V_max_c + 0.02
                [~,i1] = min(abs(vg - vlo));
                [~,i2] = min(abs(vg - vhi));
                dQ_chg(s) = abs(Qs(i2) - Qs(i1));
            end
            % else: NaN 유지 → surrogate split에서 처리
        end
    end

    if length(v_d_u) > 1
        vg = (V_bounds(end)+0.02):-0.001:(V_bounds(1)-0.02);
        Qi = interp1(v_d_u, q_d_u, vg, 'linear', 'extrap');
        Qs = movmean(Qi, 30);
        for s = 1:num_segments
            vhi = V_bounds(s+1); vlo = V_bounds(s);
            if vlo >= V_min_d - 0.02 && vhi <= V_max_d + 0.02
                [~,i1] = min(abs(vg - vhi));
                [~,i2] = min(abs(vg - vlo));
                dQ_dch(s) = abs(Qs(i2) - Qs(i1));
            end
        end
    end

    valid_chg = ~isnan(dQ_chg);
    valid_dch = ~isnan(dQ_dch);
    fprintf('  %s | Valid Chg: %d/%d  Dch: %d/%d\n', yr, sum(valid_chg), num_segments, sum(valid_dch), num_segments);

    C_eff_chg = mean(abs(I_cell(cs:ce))) / C_cell_Ah;
    C_eff_dch = mean(abs(I_cell(ds:de))) / C_cell_Ah;

    % === Split Normalization (NaN은 건드리지 않음) ===
    X_raw = [dQ_chg, dQ_dch];   % 22 features (NaN 포함 가능)

    [~,sc]=sort(abs(crate_vals_num-C_eff_chg));
    cr1c=crate_groups{sc(1)};cr2c=crate_groups{sc(2)};
    d1c=abs(C_eff_chg-crate_vals_num(sc(1)));d2c=abs(C_eff_chg-crate_vals_num(sc(2)));
    if(d1c+d2c)<eps,w1c=1;w2c=0;else,w1c=d2c/(d1c+d2c);w2c=d1c/(d1c+d2c);end
    mu_c=w1c*NormStats.(cr1c).mu_chg+w2c*NormStats.(cr2c).mu_chg;
    sig_c=w1c*NormStats.(cr1c).sigma_chg+w2c*NormStats.(cr2c).sigma_chg;

    [~,sd]=sort(abs(crate_vals_num-C_eff_dch));
    cr1d=crate_groups{sd(1)};cr2d=crate_groups{sd(2)};
    d1d=abs(C_eff_dch-crate_vals_num(sd(1)));d2d=abs(C_eff_dch-crate_vals_num(sd(2)));
    if(d1d+d2d)<eps,w1d=1;w2d=0;else,w1d=d2d/(d1d+d2d);w2d=d1d/(d1d+d2d);end
    mu_d=w1d*NormStats.(cr1d).mu_dch+w2d*NormStats.(cr2d).mu_dch;
    sig_d=w1d*NormStats.(cr1d).sigma_dch+w2d*NormStats.(cr2d).sigma_dch;

    X_norm = X_raw;
    X_norm(chg_idx) = (X_raw(chg_idx) - mu_c) ./ (sig_c + eps);  % NaN - mu = NaN → OK
    X_norm(dch_idx) = (X_raw(dch_idx) - mu_d) ./ (sig_d + eps);

    soh_est = predict(final_rf, X_norm);
    Results_Ev.(yr).SOH = soh_est;
    Results_Ev.(yr).Ceff = struct('chg',C_eff_chg,'dch',C_eff_dch);
    Results_Ev.(yr).valid_chg = valid_chg;
    Results_Ev.(yr).valid_dch = valid_dch;
    fprintf('  %s | SOH=%.2f%% | C_eff chg=%.3fC dch=%.3fC\n', yr, soh_est, C_eff_chg, C_eff_dch);
end

%% Section 3: Trajectory 시각화
fprintf('\n=== Section 3: Trajectory ===\n');
yr_list=fieldnames(Results_Ev); n_yrs=length(yr_list);
yr_nums=zeros(n_yrs,1);
for i=1:n_yrs
    tmp=regexp(yr_list{i},'\d{4}','match');yr_nums(i)=str2double(tmp{1});
    if contains(yr_list{i},'Manual'),yr_nums(i)=yr_nums(i)+0.4;end
end

fig=figure('Name','RF_v3 Field SOH','Position',[100 100 700 450],'Visible','off');
hold on; grid on; box on;
soh_vals=arrayfun(@(i)Results_Ev.(yr_list{i}).SOH,1:n_yrs);
plot(yr_nums,soh_vals,'-o','Color',[0.2 0.5 0.3],'LineWidth',2,'MarkerSize',9,'MarkerFaceColor',[0.2 0.5 0.3]);
for i=1:n_yrs
    text(yr_nums(i),soh_vals(i)+0.3,sprintf('%.1f%%',soh_vals(i)),'HorizontalAlignment','center','FontWeight','bold','FontSize',10);
    n_valid = sum(Results_Ev.(yr_list{i}).valid_chg) + sum(Results_Ev.(yr_list{i}).valid_dch);
    text(yr_nums(i),soh_vals(i)-0.8,sprintf('%d/%d segs',n_valid,num_segments*2),'HorizontalAlignment','center','FontSize',8,'Color',[0.5 0.5 0.5]);
end
x_labels=cellfun(@(s)strrep(strrep(s,'Y',''),'_',' '),yr_list,'UniformOutput',false);
set(gca,'XTick',yr_nums,'XTickLabel',x_labels,'FontSize',11); xtickangle(30);
xlabel('Year','FontSize',12); ylabel('SOH (%)','FontSize',12);
title('RF v3 (11-seg, partial support) — Field SOH Trajectory','FontSize',13,'FontWeight','bold');
xlim([min(yr_nums)-0.5 max(yr_nums)+0.5]);
saveas(fig,fullfile(rfDir,'RF_v3_Field_Trajectory.fig')); close(fig);

%% Section 4: Segment Validity Heatmap
fig2=figure('Name','v3 Segment Validity','Position',[100 100 900 300],'Visible','off');
sgtitle('v3 Segment Validity per Year','FontSize',13,'FontWeight','bold');
data_valid = zeros(n_yrs, num_segments*2);
for i=1:n_yrs
    data_valid(i,:) = [Results_Ev.(yr_list{i}).valid_chg, Results_Ev.(yr_list{i}).valid_dch];
end
seg_names = {};
for s=1:num_segments, seg_names{end+1}=sprintf('chg_%d',s); end
for s=1:num_segments, seg_names{end+1}=sprintf('dch_%d',s); end
imagesc(data_valid); colormap([0.9 0.3 0.3; 0.3 0.8 0.3]); colorbar('Ticks',[0,1],'TickLabels',{'NaN','Valid'});
set(gca,'XTick',1:length(seg_names),'XTickLabel',seg_names,'FontSize',7,...
    'YTick',1:n_yrs,'YTickLabel',yr_list);
xtickangle(45);
title(sprintf('V range [%.2f~%.2f]V',V_bounds(1),V_bounds(end)),'FontSize',11);
saveas(fig2,fullfile(rfDir,'RF_v3_SegmentValidity.fig')); close(fig2);

fprintf('Figures saved.\n');

%% Helper Functions
function segs=local_find_segments(mask)
    segs=[];i=1;n=length(mask);
    while i<=n;if mask(i),j=i;while j<n&&mask(j+1),j=j+1;end;segs=[segs;i,j];i=j+1;else,i=i+1;end;end
end
