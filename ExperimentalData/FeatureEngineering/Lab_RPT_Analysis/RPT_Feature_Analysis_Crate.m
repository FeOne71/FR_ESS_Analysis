% RPT Feature Analysis (Publication-Quality Correlation Plots)
% Charge / Discharge: separate datasets

clear; clc; close all;
warning off;

%% 0. Figure defaults
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultAxesLineWidth', 1);
set(0, 'DefaultAxesBox', 'on');
set(0, 'DefaultAxesGridAlpha', 0.3);
set(0, 'DefaultTextFontName', 'Arial');

% Diverging colormap for correlation (-1 to 1): blue-white-red
ncmap = 256;
r = [linspace(0, 1, ncmap/2), ones(1, ncmap/2)];
g = [linspace(0, 1, ncmap/2), linspace(1, 0, ncmap/2)];
b = [ones(1, ncmap/2), linspace(1, 0, ncmap/2)];
CORR_CMAP = [r(:), g(:), b(:)];

%% 1. Paths
resultRoot = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\Feature_Dataset';
figDir = fullfile(resultRoot, 'Correlation_Plots');
if ~exist(figDir, 'dir'); mkdir(figDir); end

% SOH = static capacity (same as RPT_Segment_Visualization) for consistent correlation
staticCapFile = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Capacity_Trend_Figures\Capacity_Data_Static.mat';

datasetNames = {'Charge', 'Discharge'};
datasetFiles = {fullfile(resultRoot, 'Dataset', 'Feature_Dataset_Charge.mat'), ...
                fullfile(resultRoot, 'Dataset', 'Feature_Dataset_Discharge.mat')};

fprintf('=== Starting Advanced Feature Analysis ===\n');

for d = 1:length(datasetNames)
    dName = datasetNames{d};
    datasetFile = datasetFiles{d};
    if ~exist(datasetFile, 'file')
        warning('Dataset not found: %s. Skip.', datasetFile);
        continue;
    end
    load(datasetFile);

    if d == 1
        X_data = X_data_chg; Y_data = Y_data_chg; Meta_data = Meta_data_chg; feature_names = feature_names_chg;
    else
        X_data = X_data_dchg; Y_data = Y_data_dchg; Meta_data = Meta_data_dchg; feature_names = feature_names_dchg;
    end

    % Optionally use static capacity (SOH) as Y for correlation (consistent with visualization)
    if isfile(staticCapFile)
        ld = load(staticCapFile, 'allChannelsCapacity');
        acc = ld.allChannelsCapacity;
        chKeys = fieldnames(acc);
        SOH_ref_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
        for k = 1:length(chKeys)
            ch_key = chKeys{k};
            cc = acc.(ch_key);
            if isfield(cc, 'cycles') && isfield(cc, 'capacity')
                for j = 1:length(cc.cycles)
                    key = sprintf('%s_cyc%d', ch_key, cc.cycles(j));
                    SOH_ref_map(key) = cc.capacity(j);
                end
            end
        end
        Y_data = nan(size(X_data, 1), 1);
        for r = 1:size(Meta_data, 1)
            ch = Meta_data.Channel(r); if iscell(ch), ch = ch{1}; end
            cy = Meta_data.Cycle(r);   if iscell(cy), cy = cy{1}; end
            key = sprintf('%s_%s', ch, cy);
            if SOH_ref_map.isKey(key)
                Y_data(r) = SOH_ref_map(key);
            end
        end
        fprintf('  Using static capacity (MAT) as SOH for correlation.\n');
    end

    unique_crates = unique(X_data(:, end));
    feature_labels = feature_names(1:end-1);
    n_feats = length(feature_labels);

    corr_matrix = [];
    pval_matrix = [];
    row_labels = {};
    case_stats = struct();
    idx = 1;

    %% 2. Compute Correlations and p-values
    for c = 1:length(unique_crates)
        target_c = unique_crates(c);
        mask = X_data(:, end) == target_c;
        X_sub = X_data(mask, 1:end-1);
        Y_sub = Y_data(mask);

        corr_row = nan(1, n_feats);
        pval_row = nan(1, n_feats);
        if length(Y_sub) >= 5
            for f = 1:n_feats
                valid = ~isnan(X_sub(:, f)) & ~isnan(Y_sub);
                if sum(valid) > 5
                    [corr_row(f), pval_row(f)] = corr(X_sub(valid, f), Y_sub(valid));
                end
            end
        end
        corr_matrix = [corr_matrix; corr_row];
        pval_matrix = [pval_matrix; pval_row];
        row_labels{end+1} = sprintf('%.1f C', target_c);
        case_stats(idx).Crate = target_c;
        case_stats(idx).Corrs = corr_row;
        case_stats(idx).Pvals = pval_row;
        case_stats(idx).Label = row_labels{end};
        idx = idx + 1;
    end
    n_rates = size(corr_matrix, 1);

    %% 3. Heatmap (publication-quality) with grid, R, p-value, |R|>=0.7 highlight
    f1 = figure('Visible', 'off', 'Position', [50 50 1100 420], 'Color', 'w');
    ax = axes('Parent', f1, 'Position', [0.12 0.18 0.75 0.72]);

    h = imagesc(ax, corr_matrix);
    colormap(ax, CORR_CMAP);
    caxis(ax, [-1 1]);
    ax.XLim = [0.5 n_feats+0.5];
    ax.YLim = [0.5 n_rates+0.5];
    ax.XAxis.TickValues = 1:n_feats;
    ax.XAxis.TickLabels = strrep(feature_labels, '_', ' ');
    ax.XAxis.TickLabelRotation = 45;
    ax.XAxis.FontSize = 9;
    ax.YAxis.TickValues = 1:n_rates;
    ax.YAxis.TickLabels = row_labels;
    ax.YAxis.FontSize = 10;
    xlabel(ax, 'Feature');
    ylabel(ax, 'C-rate');
    ax.LineWidth = 1;
    ax.Box = 'on';
    ax.YDir = 'normal';
    title(ax, sprintf('Feature–SOH correlation  |  %s', dName), 'FontSize', 12, 'FontWeight', 'bold');

    % Grid lines between cells
    hold(ax, 'on');
    for x = 0.5 : 1 : n_feats+0.5
        plot(ax, [x x], [0.5 n_rates+0.5], 'k-', 'LineWidth', 0.6);
    end
    for y = 0.5 : 1 : n_rates+0.5
        plot(ax, [0.5 n_feats+0.5], [y y], 'k-', 'LineWidth', 0.6);
    end
    hold(ax, 'off');

    % Highlight |R|>=0.7 with thick border
    hold(ax, 'on');
    for rr = 1:n_rates
        for cc = 1:n_feats
            val = corr_matrix(rr, cc);
            if ~isnan(val) && abs(val) >= 0.7
                rectangle(ax, 'Position', [cc-0.5 rr-0.5 1 1], 'EdgeColor', [0 0 0], 'LineWidth', 2);
            end
        end
    end
    hold(ax, 'off');

    cb = colorbar(ax, 'Location', 'eastoutside', 'AxisLocation', 'out');
    cb.Label.String = 'Correlation with SOH (R)';
    cb.Label.FontSize = 11;
    cb.LineWidth = 1;
    cb.FontSize = 10;

    % Value annotations: R and p-value (always visible with white background)
    for rr = 1:n_rates
        for cc = 1:n_feats
            val = corr_matrix(rr, cc);
            pv = pval_matrix(rr, cc);
            if ~isnan(val)
                txtR = sprintf('%.2f', val);
                txtOpt = {'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'FontSize', 8, 'FontWeight', 'bold', 'Color', [0 0 0], ...
                    'BackgroundColor', [1 1 1], 'EdgeColor', [0.6 0.6 0.6], 'Margin', 1};
                if ~isnan(pv)
                    if pv < 0.001
                        txtP = 'p<0.001';
                    else
                        txtP = sprintf('p=%.3f', pv);
                    end
                    text(ax, cc, rr+0.12, txtR, txtOpt{:});
                    text(ax, cc, rr-0.12, txtP, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                        'FontSize', 6, 'Color', [0 0 0], 'BackgroundColor', [1 1 1], 'EdgeColor', [0.6 0.6 0.6], 'Margin', 0.5);
                else
                    text(ax, cc, rr, txtR, txtOpt{:});
                end
            end
        end
    end

    saveas(f1, fullfile(figDir, sprintf('Correlation_Heatmap_%s.fig', dName)));
    try, exportgraphics(f1, fullfile(figDir, sprintf('Correlation_Heatmap_%s.png', dName)), 'Resolution', 150); catch, end

    %% 4. Bar charts (per C-rate) with p-value, |R|>=0.7 highlight, grid
    CAT_COLORS = struct('Time', [0.85 0.37 0.01], 'Cap', [0.00 0.45 0.70], ...
        'Peak', [0.47 0.67 0.19], 'Peak_Volt', [0.55 0.35 0.65], 'Other', [0.5 0.5 0.5]);

    for i = 1:length(case_stats)
        stats = case_stats(i);
        corrs = stats.Corrs;
        pvals = stats.Pvals;
        [~, sort_idx] = sort(abs(corrs), 'descend');
        top_R = corrs(sort_idx);
        top_P = pvals(sort_idx);
        top_names = feature_labels(sort_idx);
        n_plot = length(top_names);

        bar_colors = zeros(n_plot, 3);
        for k = 1:n_plot
            nm = top_names{k};
            if contains(nm, 'Time'), bar_colors(k,:) = CAT_COLORS.Time;
            elseif contains(nm, 'Cap'), bar_colors(k,:) = CAT_COLORS.Cap;
            elseif contains(nm, 'Peak_Volt'), bar_colors(k,:) = CAT_COLORS.Peak_Volt;
            elseif contains(nm, 'Peak'), bar_colors(k,:) = CAT_COLORS.Peak;
            else, bar_colors(k,:) = CAT_COLORS.Other; end
        end

        f2 = figure('Visible', 'off', 'Position', [50 50 720 540], 'Color', 'w');
        ax = axes('Parent', f2);
        b = barh(ax, top_R, 0.72);
        b.FaceColor = 'flat';
        b.CData = bar_colors;
        b.EdgeColor = 'none';
        ax.YTick = 1:n_plot;
        ax.YTickLabel = strrep(top_names, '_', ' ');
        ax.FontSize = 10;
        ax.XLim = [-1.2 1.2];
        ax.YLim = [0.4 n_plot+0.6];
        ax.YDir = 'reverse';
        ax.LineWidth = 1;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = [0.4 0.4 0.4];
        ax.GridAlpha = 0.6;
        ax.MinorGridLineStyle = '-';
        xlabel(ax, 'Correlation coefficient (R)');
        title(ax, sprintf('%s  |  %s  (* = |R|\\geq0.7)', stats.Label, dName), 'FontSize', 12, 'FontWeight', 'bold');
        grid(ax, 'on');
        hold(ax, 'on');
        plot(ax, [0 0], ax.YLim, 'k-', 'LineWidth', 1);
        hold(ax, 'off');

        % Value labels: R, p-value; bold + * when |R|>=0.7
        for k = 1:n_plot
            v = top_R(k);
            pv = top_P(k);
            xpos = v + 0.04*sign(v);
            if xpos < -1.15, xpos = -1.15 + 0.08; end
            if xpos > 1.15, xpos = 1.15 - 0.08; end
            if ~isnan(pv)
                if pv < 0.001
                    pstr = 'p<0.001';
                else
                    pstr = sprintf('p=%.3f', pv);
                end
                lbl = sprintf('%.2f  %s', v, pstr);
            else
                lbl = sprintf('%.2f', v);
            end
            if abs(v) >= 0.7
                text(ax, xpos, k, [lbl ' *'], 'FontSize', 9, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Color', [0 0 0]);
            else
                text(ax, xpos, k, lbl, 'FontSize', 8, ...
                    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
            end
        end

        % Legend (dummy markers)
        hold(ax, 'on');
        h1 = plot(ax, nan, nan, 's', 'MarkerFaceColor', CAT_COLORS.Time, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
        h2 = plot(ax, nan, nan, 's', 'MarkerFaceColor', CAT_COLORS.Cap, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
        h3 = plot(ax, nan, nan, 's', 'MarkerFaceColor', CAT_COLORS.Peak, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
        h4 = plot(ax, nan, nan, 's', 'MarkerFaceColor', CAT_COLORS.Peak_Volt, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
        lg = legend(ax, [h1 h2 h3 h4], {'Time segment', 'Capacity segment', 'Peak H/A', 'Peak voltage'}, ...
            'Location', 'southeast', 'FontSize', 9, 'Box', 'on');
        hold(ax, 'off');

        saveBase = sprintf('AllFeatures_%.1fC_%s', stats.Crate, dName);
        saveas(f2, fullfile(figDir, [saveBase '.fig']));
        try, exportgraphics(f2, fullfile(figDir, [saveBase '.png']), 'Resolution', 150); catch, end
    end

    fprintf('  %s: Heatmap + %d bar charts saved.\n', dName, length(case_stats));
end

fprintf('Analysis Complete. Results saved in %s\n', figDir);
