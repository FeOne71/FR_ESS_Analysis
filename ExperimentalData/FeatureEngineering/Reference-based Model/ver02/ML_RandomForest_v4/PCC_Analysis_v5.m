%% PCC_Analysis_v5.m — v5 세그먼트 ΔQ vs SOH PCC 분석
% - 11개 전기화학 기반 세그먼트 (3.4~4.1V)
% - ECM 보정 후 ΔQ 추출 (전채널 × 전사이클 × 전C-rate)
% - PCC 분석: 전체 풀링 + C-rate별
% - 히트맵: 충전 따로, 방전 따로
clear all; close all; clc;

%% ========================================
%  Load Data
% =========================================
d_vq  = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat');
d_ecm = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4\ECM_2RC_AllChannels_cyc0.mat');
d_cap = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Capacity_Trend_Figures\Capacity_Data_Static.mat');

%% ========================================
%  Config
% =========================================
channels = fieldnames(d_ecm.All_ECM);
crates = {'c01','c05','c1','c2','c3'};
crate_vals = [0.1, 0.5, 1, 2, 3];
clabels = {'0.1C','0.5C','1C','2C','3C'};
cyc_fields = fieldnames(d_vq.RPT_VQ_grid);

% v5 세그먼트 경계 (충전 기준: 오름차순)
VB = [3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.78, 3.84, 3.90, 4.00, 4.10];
N_seg = length(VB) - 1;
seg_labels = arrayfun(@(i) sprintf('S%d\n%.2f-%.2f', i, VB(i), VB(i+1)), 1:N_seg, 'UniformOutput', false);

saveDir = fileparts(mfilename('fullpath'));
sm_win  = 61;

%% ========================================
%  Extract ΔQ + SOH labels
% =========================================
data_CellID   = {};
data_Cycle    = [];
data_CrateNum = [];
dQ_chg_all    = [];
dQ_dch_all    = [];
data_SOH      = [];
cnt = 0;

for ch_idx = 1:length(channels)
    ch = channels{ch_idx};
    if ~isfield(d_cap.allChannelsCapacity, ch), continue; end
    cap_d = d_cap.allChannelsCapacity.(ch);
    idx0  = find(cap_d.cycles == 0, 1);
    if isempty(idx0), idx0 = 1; end
    Q_0 = max(cap_d.Q{1, idx0});
    if isnan(Q_0) || isempty(Q_0), continue; end

    ecm_chg = d_ecm.All_ECM.(ch).charge;
    ecm_dch = d_ecm.All_ECM.(ch).discharge;

    for ci = 1:length(cyc_fields)
        cyc_key = cyc_fields{ci};
        cyc_num = sscanf(cyc_key, 'cyc%d');
        if ~isfield(d_vq.RPT_VQ_grid.(cyc_key), ch), continue; end
        ch_data = d_vq.RPT_VQ_grid.(cyc_key).(ch);

        % SOH 라벨
        idx_c = find(cap_d.cycles == cyc_num, 1);
        if isempty(idx_c), continue; end
        Q_cyc = max(cap_d.Q{1, idx_c});
        SOH   = Q_cyc / Q_0 * 100;
        if isnan(SOH), continue; end

        for r = 1:length(crates)
            dQ_chg = nan(1, N_seg);
            dQ_dch = nan(1, N_seg);

            f_chg = [crates{r} '_charge'];
            if isfield(ch_data, f_chg)
                [Vc, Qp] = ecm_correct(ch_data.(f_chg), ecm_chg, Q_0, true,  sm_win);
                dQ_chg   = extract_dQ(Vc, Qp, VB, true);
            end

            f_dch = [crates{r} '_discharge'];
            if isfield(ch_data, f_dch)
                [Vc, Qp] = ecm_correct(ch_data.(f_dch), ecm_dch, Q_0, false, sm_win);
                dQ_dch   = extract_dQ(Vc, Qp, VB, false);
            end

            cnt = cnt + 1;
            data_CellID{cnt,1}  = ch;
            data_Cycle(cnt,1)   = cyc_num;
            data_CrateNum(cnt,1)= crate_vals(r);
            dQ_chg_all(cnt,:)   = dQ_chg;
            dQ_dch_all(cnt,:)   = dQ_dch;
            data_SOH(cnt,1)     = SOH;
        end
    end
    fprintf('Done: %s\n', ch);
end
fprintf('Total samples: %d\n\n', cnt);

%% ========================================
%  PCC 계산
% =========================================
% 1) 전체 풀링
pcc_chg = compute_pcc(dQ_chg_all, data_SOH);
pcc_dch = compute_pcc(dQ_dch_all, data_SOH);

% 2) C-rate별
pcc_chg_cr = nan(length(crates), N_seg);
pcc_dch_cr = nan(length(crates), N_seg);
for r = 1:length(crates)
    idx_r = data_CrateNum == crate_vals(r);
    pcc_chg_cr(r,:) = compute_pcc(dQ_chg_all(idx_r,:), data_SOH(idx_r));
    pcc_dch_cr(r,:) = compute_pcc(dQ_dch_all(idx_r,:), data_SOH(idx_r));
end

% 출력
fprintf('=== PCC: ΔQ_chg vs SOH (전체) ===\n');
for s=1:N_seg, fprintf('  Seg%2d: %+.3f\n', s, pcc_chg(s)); end
fprintf('\n=== PCC: ΔQ_dch vs SOH (전체) ===\n');
for s=1:N_seg, fprintf('  Seg%2d: %+.3f\n', s, pcc_dch(s)); end

%% ========================================
%  그림1: 전체 풀링 PCC 히트맵 (충전 / 방전)
% =========================================
fig1 = figure('Name','PCC Heatmap All C-rates','Position',[50,50,1100,320]);
modes = {'CHG','DCH'};
pcc_mat = {pcc_chg; pcc_dch};
for m = 1:2
    subplot(1,2,m);
    p = pcc_mat{m}(:)';
    imagesc(p);
    colormap(gca, redblue(256));
    clim([-1 1]); colorbar;
    xticks(1:N_seg);
    xticklabels(arrayfun(@(i) sprintf('S%d',i), 1:N_seg, 'UniformOutput',false));
    yticks(1); yticklabels({'SOH'});
    title(sprintf('PCC: \\DeltaQ_{%s} vs SOH (전체 C-rate 풀링)',modes{m}),'FontWeight','bold');
    for s=1:N_seg
        col = 'w'; if abs(p(s)) < 0.5, col='k'; end
        text(s, 1, sprintf('%.2f',p(s)), 'Color',col, ...
            'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
    end
end
sgtitle('v5 Segment PCC Heatmap','FontSize',13,'FontWeight','bold');
saveas(fig1, fullfile(saveDir,'PCC_Heatmap_v5.fig'));

%% ========================================
%  그림2: C-rate별 PCC 히트맵 (충전 / 방전)
% =========================================
fig2 = figure('Name','PCC Heatmap by C-rate','Position',[50,50,1100,480]);
for m = 1:2
    if m==1, pm = pcc_chg_cr; lbl='CHG';
    else,     pm = pcc_dch_cr; lbl='DCH'; end
    subplot(1,2,m);
    imagesc(pm); colormap(gca,redblue(256)); clim([-1 1]); colorbar;
    xticks(1:N_seg); yticks(1:5);
    xticklabels(arrayfun(@(i) sprintf('S%d',i), 1:N_seg,'UniformOutput',false));
    yticklabels(clabels);
    xlabel('Segment'); ylabel('C-rate');
    title(sprintf('PCC: \\DeltaQ_{%s} vs SOH (C-rate별)', lbl),'FontWeight','bold');
    for r=1:5, for s=1:N_seg
        if ~isnan(pm(r,s))
            col='w'; if abs(pm(r,s))<0.5, col='k'; end
            text(s,r,sprintf('%.2f',pm(r,s)),'Color',col,...
                'HorizontalAlignment','center','FontSize',7.5);
        end
    end; end
end
sgtitle('v5 Segment PCC (C-rate별)','FontSize',13,'FontWeight','bold');
saveas(fig2, fullfile(saveDir,'PCC_Heatmap_CRate_v5.fig'));

%% 데이터 저장
save(fullfile(saveDir,'PCC_v5_data.mat'), ...
    'dQ_chg_all','dQ_dch_all','data_SOH','data_CellID','data_Cycle','data_CrateNum', ...
    'pcc_chg','pcc_dch','pcc_chg_cr','pcc_dch_cr','VB','N_seg');
fprintf('\nSaved: PCC_Heatmap_v5.fig, PCC_Heatmap_CRate_v5.fig, PCC_v5_data.mat\n');

%% ========================================
%  Helper Functions
% =========================================
function [V_corr, Q_plot] = ecm_correct(s, ecm, Q_0, is_charge, sm_win)
    V_raw = double(s.V_raw(:)); I_raw = double(s.I_raw(:));
    t_s   = seconds(s.t_raw(:)); Q_raw = double(s.Q_raw(:));
    Q_cum = cumtrapz(t_s, abs(I_raw)) / 3600;
    if is_charge, SOC = Q_cum / Q_0 * 100;
    else,         SOC = 100 - Q_cum / Q_0 * 100; end
    SOC = max(0, min(100, SOC));
    [ss, si] = sort(ecm.SOC);
    Rt = (ecm.R0(si) + ecm.R1(si) + ecm.R2(si)) / 1e3;
    vld = ss >= 10 & ss <= 95; ss = ss(vld); Rt = Rt(vld);
    SC  = max(min(ss), min(max(ss), SOC));
    Ri  = interp1(ss, Rt, SC, 'linear');
    V_corr = movmean(V_raw - I_raw .* Ri, sm_win);
    Q_plot = Q_raw - min(Q_raw);
end

function dQ = extract_dQ(V_corr, Q_plot, VB, is_charge)
    N  = length(VB) - 1;
    dQ = nan(1, N);
    [Vu, uid] = unique(V_corr, 'stable');
    Qu = Q_plot(uid);
    mono = true(size(Vu));
    if is_charge
        for ii = 2:length(Vu), if Vu(ii) <= Vu(ii-1), mono(ii) = false; end; end
    else
        for ii = 2:length(Vu), if Vu(ii) >= Vu(ii-1), mono(ii) = false; end; end
    end
    Vu = Vu(mono); Qu = Qu(mono);
    if Vu(1) > Vu(end), Vu = flipud(Vu); Qu = flipud(Qu); end
    if length(Vu) < 10, return; end
    for s = 1:N
        if VB(s) >= min(Vu) && VB(s+1) <= max(Vu)
            q1 = interp1(Vu, Qu, VB(s),   'linear');
            q2 = interp1(Vu, Qu, VB(s+1), 'linear');
            if ~isnan(q1) && ~isnan(q2), dQ(s) = abs(q2 - q1); end
        end
    end
end

function pcc = compute_pcc(X, y)
    N_seg = size(X, 2);
    pcc   = nan(1, N_seg);
    for s = 1:N_seg
        valid = ~isnan(X(:,s)) & ~isnan(y);
        if sum(valid) > 10
            r = corrcoef(X(valid,s), y(valid));
            pcc(s) = r(1,2);
        end
    end
end

function c = redblue(n)
    if nargin < 1, n = 256; end
    half = floor(n/2);
    r = [linspace(0,1,half)'; ones(n-half,1)];
    b = [ones(half,1); linspace(1,0,n-half)'];
    g = [linspace(0,1,half)'; linspace(1,0,n-half)'];
    c = [r g b];
end
