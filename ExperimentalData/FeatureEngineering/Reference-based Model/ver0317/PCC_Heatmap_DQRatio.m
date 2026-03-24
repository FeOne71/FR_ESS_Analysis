% PCC_Heatmap_DQRatio.m
% dQ ratio 피처와 SOH 간의 PCC 히트맵 생성

clear; clc; close all;

%% Load
verDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver0317';
d = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat')); FM = d.FM;
mr = load(fullfile(verDir, 'MasterRulers_PseudoOCV.mat'));
Vb_c = mr.MasterRuler_ver0317.V_bounds_chg;
Vb_d = mr.MasterRuler_ver0317.V_bounds_dch;
Q_nom = 64;
visDir = fullfile(verDir, 'RatioModel', 'Visualization');

%% Compute dQ ratio
dQ_c_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i), 3:12, 'UniformOutput', false);
dQ_d_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false);
dQ_c_raw = FM{:, dQ_c_cols};
dQ_d_raw = FM{:, dQ_d_cols};
dQ_c_ratio = dQ_c_raw ./ sum(dQ_c_raw, 2, 'omitnan');
dQ_d_ratio = dQ_d_raw ./ sum(dQ_d_raw, 2, 'omitnan');
y = FM.Static_Capacity / Q_nom * 100;
conds = FM.Condition;

uConds    = {'c01','c05','c1','c2','c3'};
cr_labels = {'All','0.1C','0.5C','1.0C','2.0C','3.0C'};
nSeg_c = 10; nSeg_d = 11; nRows = 6;

%% Compute PCC
PCC_c = nan(nRows, nSeg_c); Pval_c = nan(nRows, nSeg_c);
PCC_d = nan(nRows, nSeg_d); Pval_d = nan(nRows, nSeg_d);

for s = 1:nSeg_c
    x = dQ_c_ratio(:,s); m = ~isnan(x);
    if sum(m)>3, [r,p]=corr(x(m),y(m)); PCC_c(1,s)=r; Pval_c(1,s)=p; end
    for ci=1:5
        m2 = m & strcmp(conds, uConds{ci});
        if sum(m2)>3, [r,p]=corr(x(m2),y(m2)); PCC_c(ci+1,s)=r; Pval_c(ci+1,s)=p; end
    end
end
for s = 1:nSeg_d
    x = dQ_d_ratio(:,s); m = ~isnan(x);
    if sum(m)>3, [r,p]=corr(x(m),y(m)); PCC_d(1,s)=r; Pval_d(1,s)=p; end
    for ci=1:5
        m2 = m & strcmp(conds, uConds{ci});
        if sum(m2)>3, [r,p]=corr(x(m2),y(m2)); PCC_d(ci+1,s)=r; Pval_d(ci+1,s)=p; end
    end
end

%% Segment labels
seg_lbl_c = arrayfun(@(i) sprintf('Seg%02d', i), 3:12, 'UniformOutput', false);
seg_lbl_d = arrayfun(@(i) sprintf('Seg%02d', i), 1:11, 'UniformOutput', false);

%% Helper: star
function s = star(p)
    if isnan(p), s='NaN';
    elseif p<0.001, s='***';
    elseif p<0.01,  s='**';
    elseif p<0.05,  s='*';
    else, s=''; end
end

%% Plot charge heatmap
fig1 = figure('Position',[50 50 1300 440]);
ax1 = axes('Parent',fig1);
imagesc(ax1, PCC_c, [-1 1]);
colormap(ax1, redblue_cmap(256));
cb = colorbar(ax1); cb.Label.String = 'Pearson Correlation (r)'; cb.FontSize=10;
for r=1:nRows
    for s=1:nSeg_c
        if ~isnan(PCC_c(r,s))
            txt = sprintf('%.2f%s', PCC_c(r,s), star(Pval_c(r,s)));
            clr = [1 1 1]; if abs(PCC_c(r,s))<0.5, clr=[0 0 0]; end
            text(ax1, s, r, txt, 'HorizontalAlignment','center','FontSize',9,'Color',clr,'FontWeight','bold');
        else
            text(ax1, s, r, 'NaN', 'HorizontalAlignment','center','FontSize',8,'Color',[0.9 0.9 0.9]);
        end
    end
end
set(ax1, 'XTick',1:nSeg_c, 'XTickLabel', arrayfun(@(i) sprintf('Seg%02d',i),3:12,'UniformOutput',false), ...
    'YTick',1:nRows, 'YTickLabel', cr_labels, 'FontSize',10, 'TickLength',[0 0]);
xtickangle(ax1, 30);
title(ax1, 'PCC(dQ_{charge} ratio, SOH) per Segment and Charging Rate', 'FontSize',12,'FontWeight','bold');
xlabel(ax1,'Charge Segment','FontSize',11);
ylabel(ax1,'Charging Rate (C-rate)','FontSize',11);
saveas(fig1, fullfile(visDir, 'PCC_Heatmap_Charge_Ratio.png'));

%% Plot discharge heatmap
fig2 = figure('Position',[50 50 1400 440]);
ax2 = axes('Parent',fig2);
imagesc(ax2, PCC_d, [-1 1]);
colormap(ax2, redblue_cmap(256));
cb2 = colorbar(ax2); cb2.Label.String = 'Pearson Correlation (r)'; cb2.FontSize=10;
for r=1:nRows
    for s=1:nSeg_d
        if ~isnan(PCC_d(r,s))
            txt = sprintf('%.2f%s', PCC_d(r,s), star(Pval_d(r,s)));
            clr = [1 1 1]; if abs(PCC_d(r,s))<0.5, clr=[0 0 0]; end
            text(ax2, s, r, txt, 'HorizontalAlignment','center','FontSize',9,'Color',clr,'FontWeight','bold');
        else
            text(ax2, s, r, 'NaN', 'HorizontalAlignment','center','FontSize',8,'Color',[0.9 0.9 0.9]);
        end
    end
end
set(ax2, 'XTick',1:nSeg_d, 'XTickLabel', arrayfun(@(i) sprintf('Seg%02d',i),1:11,'UniformOutput',false), ...
    'YTick',1:nRows, 'YTickLabel', cr_labels, 'FontSize',10, 'TickLength',[0 0]);
xtickangle(ax2, 30);
title(ax2, 'PCC(dQ_{discharge} ratio, SOH) per Segment and Discharging Rate', 'FontSize',12,'FontWeight','bold');
xlabel(ax2,'Discharge Segment','FontSize',11);
ylabel(ax2,'Discharging Rate (C-rate)','FontSize',11);
saveas(fig2, fullfile(visDir, 'PCC_Heatmap_Discharge_Ratio.png'));

fprintf('Saved both heatmaps.\n');

%% Colormap
function cmap = redblue_cmap(n)
    n1 = floor(n/2); n2 = ceil(n/2);
    c1 = [linspace(0,1,n1)', linspace(0,1,n1)', ones(n1,1)];
    c2 = [ones(n2,1), linspace(1,0,n2)', linspace(1,0,n2)'];
    cmap = [c1; c2];
end
