% RPT Segment Detailed Visualization (Publication-Quality)
% Charge: 5 segments (3.7-3.95V), Discharge: 4 segments (3.75-3.87V)

clear; clc; close all;
warning off;

%% 0. Figure defaults (professional styling)
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultAxesLineWidth', 1);
set(0, 'DefaultAxesBox', 'on');
set(0, 'DefaultAxesGridAlpha', 0.25);
set(0, 'DefaultAxesMinorGridAlpha', 0.1);
set(0, 'DefaultLineLineWidth', 1.8);
set(0, 'DefaultTextFontName', 'Arial');

% Segment colors (distinct, colorblind-friendly)
SEG_COLORS = [
    0.00 0.45 0.70   % blue
    0.85 0.37 0.01   % orange
    0.47 0.67 0.19   % green
    0.80 0.47 0.74   % magenta
    0.50 0.50 0.50   % gray
    ];

%% 1. Paths
resultRoot = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\Feature_Dataset';
figDir = fullfile(resultRoot, 'Segment_Plots');
if ~exist(figDir, 'dir'); mkdir(figDir); end

% Static capacity (SOH reference): from Reference/Cap_Trend_Summary_Dynamic.m output
staticCapFile = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Capacity_Trend_Figures\Capacity_Data_Static.mat';

datasetNames = {'Charge', 'Discharge'};
datasetFiles = {fullfile(resultRoot, 'Dataset', 'Feature_Dataset_Charge.mat'), ...
                fullfile(resultRoot, 'Dataset', 'Feature_Dataset_Discharge.mat')};

%% 2. Visualize per dataset
fprintf('=== Generating Segment Visualization Plots ===\n');

for d = 1:length(datasetNames)
    dName = datasetNames{d};
    datasetFile = datasetFiles{d};
    if ~exist(datasetFile, 'file')
        warning('Dataset not found: %s. Skip.', datasetFile);
        continue;
    end
    load(datasetFile);
    if d == 1
        X_data = X_data_chg; Y_data = Y_data_chg; Meta_data = Meta_data_chg;
        feature_names = feature_names_chg; segments = segments_chg;
    else
        X_data = X_data_dchg; Y_data = Y_data_dchg; Meta_data = Meta_data_dchg;
        feature_names = feature_names_dchg; segments = segments_dchg;
    end

    n_seg = length(segments) - 1;
    unique_crates = unique(X_data(:, end));
    % Feature cols: [n_seg Cap, n_seg Time, Peak_Height, Peak_Volt, Peak_Area, Crate]
    idx_height = size(X_data, 2) - 3;
    idx_volt   = size(X_data, 2) - 2;
    idx_area   = size(X_data, 2) - 1;

    % SOH = Static capacity (reference per Channel, Cycle) from Capacity_Data_Static.mat
    SOH_ref_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    if isfile(staticCapFile)
        ld = load(staticCapFile, 'allChannelsCapacity');
        acc = ld.allChannelsCapacity;
        chKeys = fieldnames(acc);
        for k = 1:length(chKeys)
            ch_key = chKeys{k};
            cc = acc.(ch_key);
            if isfield(cc, 'cycles') && isfield(cc, 'capacity')
                cyc = cc.cycles;
                cap = cc.capacity;
                for j = 1:length(cyc)
                    key = sprintf('%s_cyc%d', ch_key, cyc(j));
                    SOH_ref_map(key) = cap(j);
                end
            end
        end
        fprintf('  SOH (static capacity from MAT) lookup: %d (Channel,Cycle) pairs.\n', SOH_ref_map.Count);
    else
        % Fallback: 0.5C capacity from current dataset
        idx_05 = abs(Meta_data.Crate - 0.5) < 0.01;
        for k = find(idx_05)'
            ch = Meta_data.Channel(k); if iscell(ch), ch = ch{1}; end
            cy = Meta_data.Cycle(k);   if iscell(cy), cy = cy{1}; end
            key = sprintf('%s_%s', ch, cy);
            SOH_ref_map(key) = Meta_data.SOH(k);
        end
        fprintf('  SOH (0.5C capacity fallback) lookup: %d (Channel,Cycle) pairs.\n', SOH_ref_map.Count);
    end

    fprintf('\n--- %s (%d segments) ---\n', dName, n_seg);

    for i = 1:length(unique_crates)
        target_c = unique_crates(i);
        fprintf('  Plotting %.1f C %s...\n', target_c, dName);

        mask = X_data(:, end) == target_c;
        X_sub = X_data(mask, :);
        Meta_sub = Meta_data(mask, :);
        n_pts = size(Meta_sub, 1);
        soh_ref = nan(n_pts, 1);
        for j = 1:n_pts
            ch = Meta_sub.Channel(j); if iscell(ch), ch = ch{1}; end
            cy = Meta_sub.Cycle(j);   if iscell(cy), cy = cy{1}; end
            key = sprintf('%s_%s', ch, cy);
            if SOH_ref_map.isKey(key)
                soh_ref(j) = SOH_ref_map(key);
            end
        end
        [soh_sorted, sidx] = sort(soh_ref, 'descend');
        valid = ~isnan(soh_sorted);
        if sum(valid) < 2, continue; end
        soh_plot = soh_sorted(valid);
        X_sorted = X_sub(sidx, :);

        rows_sp = ceil(sqrt(n_seg));
        cols_sp = ceil(n_seg / rows_sp);

        % ----- Segment Time vs SOH -----
        f1 = figure('Visible', 'off', 'Position', [50 50 1000 720], 'Color', 'w');
        for s = 1:n_seg
            ax = subplot(rows_sp, cols_sp, s);
            time_col = n_seg + s;
            val = X_sorted(valid, time_col);
            c = SEG_COLORS(mod(s-1, size(SEG_COLORS,1))+1, :);
            plot(ax, soh_plot, val, '-', 'Color', c, 'LineWidth', 2);
            hold(ax, 'on');
            scatter(ax, soh_plot, val, 28, c, 'filled', 'MarkerEdgeColor', 'none');
            hold(ax, 'off');
            set(ax, 'XDir', 'reverse');
            title(ax, strrep(feature_names{n_seg+s}, '_', ' '), 'FontWeight', 'normal', 'FontSize', 11);
            xlabel(ax, 'SOH (static capacity, Ah) \rightarrow decrease');
            ylabel(ax, 'Time (s)');
            set(ax, 'FontSize', 10, 'LineWidth', 1, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.15);
            grid(ax, 'on');
            if length(val) > 1 && val(1) > 0
                deg = (val(end) - val(1)) / val(1) * 100;
                xlims = get(ax, 'XLim'); ylims = get(ax, 'YLim');
                text(ax, xlims(1) + 0.02*diff(xlims), ylims(2) - 0.08*diff(ylims), ...
                    sprintf('\\Delta = %+.1f%%', deg), 'FontSize', 9, 'FontWeight', 'bold', ...
                    'BackgroundColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7], 'VerticalAlignment', 'top');
            end
        end
        sgtitle(f1, sprintf('Segment time vs. SOH  |  %s  |  %.1f C', dName, target_c), 'FontSize', 13, 'FontWeight', 'bold');
        saveName1 = fullfile(figDir, sprintf('SegTime_vs_SOH_%.1fC_%s.fig', target_c, dName));
        saveas(f1, saveName1);
        try, exportgraphics(f1, strrep(saveName1, '.fig'), 'Resolution', 150); catch, end

        % ----- Segment Capacity vs SOH -----
        f2 = figure('Visible', 'off', 'Position', [50 50 1000 720], 'Color', 'w');
        for s = 1:n_seg
            ax = subplot(rows_sp, cols_sp, s);
            val = X_sorted(valid, s);
            c = SEG_COLORS(mod(s-1, size(SEG_COLORS,1))+1, :);
            plot(ax, soh_plot, val, '-', 'Color', c, 'LineWidth', 2);
            hold(ax, 'on');
            scatter(ax, soh_plot, val, 28, c, 'filled', 'MarkerEdgeColor', 'none');
            hold(ax, 'off');
            set(ax, 'XDir', 'reverse');
            title(ax, strrep(feature_names{s}, '_', ' '), 'FontWeight', 'normal', 'FontSize', 11);
            xlabel(ax, 'SOH (static capacity, Ah) \rightarrow decrease');
            ylabel(ax, 'Capacity (Ah)');
            set(ax, 'FontSize', 10, 'LineWidth', 1, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.15);
            grid(ax, 'on');
            if length(val) > 1 && val(1) > 0
                deg = (val(end) - val(1)) / val(1) * 100;
                xlims = get(ax, 'XLim'); ylims = get(ax, 'YLim');
                text(ax, xlims(1) + 0.02*diff(xlims), ylims(2) - 0.08*diff(ylims), ...
                    sprintf('\\Delta = %+.1f%%', deg), 'FontSize', 9, 'FontWeight', 'bold', ...
                    'BackgroundColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7], 'VerticalAlignment', 'top');
            end
        end
        sgtitle(f2, sprintf('Segment capacity vs. SOH  |  %s  |  %.1f C', dName, target_c), 'FontSize', 13, 'FontWeight', 'bold');
        saveName2 = fullfile(figDir, sprintf('SegCap_vs_SOH_%.1fC_%s.fig', target_c, dName));
        saveas(f2, saveName2);
        try, exportgraphics(f2, strrep(saveName2, '.fig'), 'Resolution', 150); catch, end
    end

    % ----- Peak Height, Peak Voltage, Peak Area vs SOH (static capacity from MAT) -----
    n_rates = length(unique_crates);
    cols_plot = ceil(sqrt(n_rates));
    rows_plot = ceil(n_rates / cols_plot);

    f3 = figure('Visible', 'off', 'Position', [50 50 1100 640], 'Color', 'w');
    f4 = figure('Visible', 'off', 'Position', [50 50 1100 640], 'Color', 'w');
    f5 = figure('Visible', 'off', 'Position', [50 50 1100 640], 'Color', 'w');

    for i = 1:n_rates
        target_c = unique_crates(i);
        mask = X_data(:, end) == target_c;
        X_sub = X_data(mask, :);
        Meta_sub = Meta_data(mask, :);
        n_pts = size(Meta_sub, 1);
        soh_ref = nan(n_pts, 1);
        for j = 1:n_pts
            ch = Meta_sub.Channel(j); if iscell(ch), ch = ch{1}; end
            cy = Meta_sub.Cycle(j);   if iscell(cy), cy = cy{1}; end
            key = sprintf('%s_%s', ch, cy);
            if SOH_ref_map.isKey(key)
                soh_ref(j) = SOH_ref_map(key);
            end
        end
        % x-axis: left = high SOH, right = low SOH (SOH decreases along x)
        [soh_sorted, sidx] = sort(soh_ref, 'descend');
        valid = ~isnan(soh_sorted);
        if sum(valid) < 2, continue; end
        soh_plot = soh_sorted(valid);
        X_sorted = X_sub(sidx, :);
        c = SEG_COLORS(mod(i-1, size(SEG_COLORS,1))+1, :);

        % Peak vs SOH: x-axis reversed so SOH decreases left → right
        figure(f3);
        ax = subplot(rows_plot, cols_plot, i);
        plot(ax, soh_plot, X_sorted(valid, idx_volt), '-', 'Color', c, 'LineWidth', 2);
        hold(ax, 'on');
        scatter(ax, soh_plot, X_sorted(valid, idx_volt), 32, c, 'filled', 'MarkerEdgeColor', 'none');
        hold(ax, 'off');
        set(ax, 'XDir', 'reverse');  % high SOH left, low SOH right
        title(ax, sprintf('%.1f C', target_c), 'FontWeight', 'normal');
        xlabel(ax, 'SOH (static capacity, Ah) \rightarrow decrease');
        ylabel(ax, 'Peak voltage (V)');
        set(ax, 'FontSize', 10, 'LineWidth', 1, 'GridAlpha', 0.3);
        grid(ax, 'on');

        figure(f4);
        ax = subplot(rows_plot, cols_plot, i);
        plot(ax, soh_plot, X_sorted(valid, idx_area), '-', 'Color', c, 'LineWidth', 2);
        hold(ax, 'on');
        scatter(ax, soh_plot, X_sorted(valid, idx_area), 32, c, 'filled', 'MarkerEdgeColor', 'none');
        hold(ax, 'off');
        set(ax, 'XDir', 'reverse');  % high SOH left, low SOH right
        title(ax, sprintf('%.1f C', target_c), 'FontWeight', 'normal');
        xlabel(ax, 'SOH (static capacity, Ah) \rightarrow decrease');
        ylabel(ax, 'Peak area (Ah/V)');
        set(ax, 'FontSize', 10, 'LineWidth', 1, 'GridAlpha', 0.3);
        grid(ax, 'on');

        figure(f5);
        ax = subplot(rows_plot, cols_plot, i);
        plot(ax, soh_plot, X_sorted(valid, idx_height), '-', 'Color', c, 'LineWidth', 2);
        hold(ax, 'on');
        scatter(ax, soh_plot, X_sorted(valid, idx_height), 32, c, 'filled', 'MarkerEdgeColor', 'none');
        hold(ax, 'off');
        set(ax, 'XDir', 'reverse');  % high SOH left, low SOH right
        title(ax, sprintf('%.1f C', target_c), 'FontWeight', 'normal');
        xlabel(ax, 'SOH (static capacity, Ah) \rightarrow decrease');
        ylabel(ax, 'Peak height (Ah/V)');
        set(ax, 'FontSize', 10, 'LineWidth', 1, 'GridAlpha', 0.3);
        grid(ax, 'on');
    end

    sgtitle(f3, sprintf('dQ/dV peak voltage vs SOH  |  %s', dName), 'FontSize', 13, 'FontWeight', 'bold');
    sgtitle(f4, sprintf('dQ/dV peak area vs SOH  |  %s', dName), 'FontSize', 13, 'FontWeight', 'bold');
    sgtitle(f5, sprintf('dQ/dV peak height vs SOH  |  %s', dName), 'FontSize', 13, 'FontWeight', 'bold');

    saveas(f3, fullfile(figDir, sprintf('Peak_Volt_vs_SOH_%s.fig', dName)));
    saveas(f4, fullfile(figDir, sprintf('Peak_Area_vs_SOH_%s.fig', dName)));
    saveas(f5, fullfile(figDir, sprintf('Peak_Height_vs_SOH_%s.fig', dName)));

end

fprintf('All plots saved in: %s\n', figDir);
