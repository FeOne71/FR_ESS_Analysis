% Plot_SOH_vs_DQRatio_Heatmap.m
% SOH (y축) × dQ ratio 피처 (x축) 히트맵
% C-rate별 분리, 색 = dQ ratio 값

clear; clc; close all;

%% Paths
verDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver0317';
d = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat')); FM = d.FM;
visDir = fullfile(verDir, 'RatioModel', 'Visualization');
Q_nom = 64;

%% Feature names
chg_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i), 3:12, 'UniformOutput', false);
dch_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false);
feat_labels_c = arrayfun(@(i) sprintf('Seg%02d',i), 3:12, 'UniformOutput', false);
feat_labels_d = arrayfun(@(i) sprintf('Seg%02d',i), 1:11, 'UniformOutput', false);

%% C-rate conditions
uConds    = {'c01','c05','c1','c2','c3'};
cr_labels = {'0.1C','0.5C','1.0C','2.0C','3.0C'};

%% ===== CHARGE heatmap =====
fig_c = figure('Position',[30 30 1600 800]);
sgtitle('SOH vs. dQ_{charge} ratio 피처 (Lab 데이터, C-rate별)', 'FontSize',13,'FontWeight','bold');

for ci = 1:5
    cond = uConds{ci};
    idx  = strcmp(FM.Condition, cond);
    FM_sub = FM(idx,:);

    % Compute dQ ratio
    dQ_raw = FM_sub{:, chg_cols};
    dQ_rat = dQ_raw ./ sum(dQ_raw, 2, 'omitnan');

    % SOH
    soh = FM_sub.Static_Capacity / Q_nom * 100;

    % Sort by SOH descending
    [soh_sorted, si] = sort(soh, 'descend');
    mat = dQ_rat(si, :);

    ax = subplot(1,5,ci);
    imagesc(ax, mat', [0 0.5]);
    colormap(ax, hot);
    if ci==5, colorbar(ax); end

    % y-tick: feature labels
    nS = length(soh_sorted);
    xtick_idx = round(linspace(1, nS, min(7,nS)));
    set(ax, 'YTick',1:10, 'YTickLabel', feat_labels_c, ...
        'XTick', xtick_idx, ...
        'XTickLabel', arrayfun(@(s) sprintf('%.1f%%',s), soh_sorted(xtick_idx), 'UniformOutput',false), ...
        'FontSize', 9, 'TickLength',[0 0]);
    xtickangle(ax, 45);
    title(ax, cr_labels{ci}, 'FontSize',11,'FontWeight','bold');
    if ci==1
        ylabel(ax,'Charge Segment','FontSize',10);
    end
    xlabel(ax,'SOH (%)','FontSize',9);
end

saveas(fig_c, fullfile(visDir, 'SOH_vs_DQRatio_Charge_Heatmap.png'));

%% ===== DISCHARGE heatmap =====
fig_d = figure('Position',[30 30 1600 800]);
sgtitle('SOH vs. dQ_{discharge} ratio 피처 (Lab 데이터, C-rate별)', 'FontSize',13,'FontWeight','bold');

for ci = 1:5
    cond = uConds{ci};
    idx  = strcmp(FM.Condition, cond);
    FM_sub = FM(idx,:);

    dQ_raw = FM_sub{:, dch_cols};
    dQ_rat = dQ_raw ./ sum(dQ_raw, 2, 'omitnan');

    soh = FM_sub.Static_Capacity / Q_nom * 100;
    [soh_sorted, si] = sort(soh, 'descend');
    mat = dQ_rat(si, :);

    ax = subplot(1,5,ci);
    imagesc(ax, mat', [0 0.5]);
    colormap(ax, hot);
    if ci==5, colorbar(ax); end

    nS = length(soh_sorted);
    xtick_idx = round(linspace(1, nS, min(7,nS)));
    set(ax, 'YTick',1:11, 'YTickLabel', feat_labels_d, ...
        'XTick', xtick_idx, ...
        'XTickLabel', arrayfun(@(s) sprintf('%.1f%%',s), soh_sorted(xtick_idx), 'UniformOutput',false), ...
        'FontSize', 9, 'TickLength',[0 0]);
    xtickangle(ax, 45);
    title(ax, cr_labels{ci}, 'FontSize',11,'FontWeight','bold');
    if ci==1
        ylabel(ax,'Discharge Segment','FontSize',10);
    end
    xlabel(ax,'SOH (%)','FontSize',9);
end

saveas(fig_d, fullfile(visDir, 'SOH_vs_DQRatio_Discharge_Heatmap.png'));

fprintf('Saved both SOH vs feature heatmaps.\n');
