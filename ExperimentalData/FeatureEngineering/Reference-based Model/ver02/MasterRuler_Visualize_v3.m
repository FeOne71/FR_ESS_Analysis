%% MasterRuler_Visualize_v3.m
% Visualizes the 5-segment voltage boundaries from MasterRulers_v3.mat
% over 0.5C Charge & Discharge curves for ALL CHANNELS & ALL CYCLES.
% (Similar to ver01 RPT_Pipeline_Visualization Phase 1-4)
% Location: ver02 folder

clear; clc; close all;

%% ========================================================================
% Section 1: Configuration & File Paths
% ========================================================================
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_rpt_vq_mat = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');
path_master_ruler = fullfile(baseDir, 'FeatureEngineering', 'MasterRulers_v3.mat');

disp('Loading Data...');
load(path_rpt_vq_mat, 'RPT_VQ_grid');
load(path_master_ruler, 'MasterRulers');

channels = fieldnames(MasterRulers);
cyc_fields = fieldnames(RPT_VQ_grid);

% Master Rulers are generated to be universal across all channels.
% Thus, any channel's V_bounds_chg/dch can be used for drawing lines.
% We pick the first valid channel's ruler.
sample_ch = channels{1};
Global_V_bounds_chg = MasterRulers.(sample_ch).V_bounds_chg;
Global_V_bounds_dch = MasterRulers.(sample_ch).V_bounds_dch;
num_segments = length(Global_V_bounds_chg) - 1;

%% ========================================================================
% Section 2: Drawing 4-Subplot Visualization
% ========================================================================
fig = figure('Position', [50, 50, 1800, 900], 'Name', 'Master Ruler v3 Final Result (All Channels)');

% Define 5 Segment Colors (Pastel-like)
seg_colors = [
    0.6 0.8 1.0;  % Seg1 (Light Blue)
    1.0 0.7 0.7;  % Seg2 (Light Red)
    0.7 0.7 1.0;  % Seg3 (Periwinkle)
    1.0 0.8 0.6;  % Seg4 (Light Orange)
    0.8 0.6 0.9;  % Seg5 (Light Purple)
];

ch_colors = lines(length(channels));

%% --- Subplot 1: Charge V-Q (Full Range: 3.0 ~ 4.2V) ---
subplot(2,2,1); hold on;
for ch_idx = 1:length(channels)
    ch = channels{ch_idx};
    for cyc_idx = 1:length(cyc_fields)
        cyc = cyc_fields{cyc_idx};
        if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && ...
           isfield(RPT_VQ_grid.(cyc).(ch), 'c05_charge')
            data = RPT_VQ_grid.(cyc).(ch).c05_charge;
            if strcmp(cyc, 'cyc0')
                plot(data.Q, data.V_grid, 'Color', ch_colors(ch_idx,:), 'LineWidth', 1.5);
            else
                % High transparency for degrading cycles
                plot(data.Q, data.V_grid, 'Color', [ch_colors(ch_idx,:) 0.3], 'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
        end
    end
end
xl = xlim;
% Draw Master Ruler Patches (Charge: Ascending)
for i = 1:num_segments
    V_start = Global_V_bounds_chg(i);
    V_end = Global_V_bounds_chg(i+1);
    patch([xl(1) xl(2) xl(2) xl(1)], [V_start V_start V_end V_end], ...
          seg_colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
% Draw Master Ruler Lines
for i = 1:length(Global_V_bounds_chg)
    line(xl, [Global_V_bounds_chg(i) Global_V_bounds_chg(i)], ...
         'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
end
xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
title('Charge Voltage Segment (All Cycles 0.5C)', 'FontSize', 13, 'FontWeight', 'bold');
ylim([3.0 4.2]);
grid on; hold off;

%% --- Subplot 2: Discharge V-Q (Full Range: 3.0 ~ 4.2V) ---
subplot(2,2,2); hold on;
for ch_idx = 1:length(channels)
    ch = channels{ch_idx};
    for cyc_idx = 1:length(cyc_fields)
        cyc = cyc_fields{cyc_idx};
        if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && ...
           isfield(RPT_VQ_grid.(cyc).(ch), 'c05_discharge')
            data = RPT_VQ_grid.(cyc).(ch).c05_discharge;
            if strcmp(cyc, 'cyc0')
                plot(data.Q, data.V_grid, 'Color', ch_colors(ch_idx,:), 'LineWidth', 1.5);
            else
                plot(data.Q, data.V_grid, 'Color', [ch_colors(ch_idx,:) 0.3], 'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
        end
    end
end
xl = xlim;
% Draw Master Ruler Patches (Discharge: Descending)
for i = 1:num_segments
    V_max = Global_V_bounds_dch(i);    % Top
    V_min = Global_V_bounds_dch(i+1);  % Bottom
    patch([xl(1) xl(2) xl(2) xl(1)], [V_min V_min V_max V_max], ...
          seg_colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
for i = 1:length(Global_V_bounds_dch)
    line(xl, [Global_V_bounds_dch(i) Global_V_bounds_dch(i)], ...
         'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
end
xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
title('Discharge Voltage Segment (All Cycles 0.5C)', 'FontSize', 13, 'FontWeight', 'bold');
ylim([3.0 4.2]);
grid on; hold off;

%% --- Subplot 3: Charge V-Q (Zoom Range within Rulers) ---
subplot(2,2,3); hold on;
for ch_idx = 1:length(channels)
    ch = channels{ch_idx};
    for cyc_idx = 1:length(cyc_fields)
        cyc = cyc_fields{cyc_idx};
        if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && ...
           isfield(RPT_VQ_grid.(cyc).(ch), 'c05_charge')
            data = RPT_VQ_grid.(cyc).(ch).c05_charge;
            if strcmp(cyc, 'cyc0')
                plot(data.Q, data.V_grid, 'Color', ch_colors(ch_idx,:), 'LineWidth', 1.5);
            else
                plot(data.Q, data.V_grid, 'Color', [ch_colors(ch_idx,:) 0.3], 'LineWidth', 0.8);
            end
        end
    end
end
xl = xlim;
for i = 1:num_segments
    V_start = Global_V_bounds_chg(i);
    V_end = Global_V_bounds_chg(i+1);
    patch([xl(1) xl(2) xl(2) xl(1)], [V_start V_start V_end V_end], ...
          seg_colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
for i = 1:length(Global_V_bounds_chg)
    V_line = Global_V_bounds_chg(i);
    line(xl, [V_line V_line], ...
         'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
end
xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 11);
title('Charge Segment Zoom', 'FontSize', 12, 'FontWeight', 'bold');
% Determine precise zoom window based on active common bounds
ylim([Global_V_bounds_chg(1)-0.01, Global_V_bounds_chg(end)+0.01]);
xlim([20 60]); % Approximate zoom-in range for capacity
grid on; hold off;

%% --- Subplot 4: Discharge V-Q (Zoom Range within Rulers) ---
subplot(2,2,4); hold on;
for ch_idx = 1:length(channels)
    ch = channels{ch_idx};
    for cyc_idx = 1:length(cyc_fields)
        cyc = cyc_fields{cyc_idx};
        if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && ...
           isfield(RPT_VQ_grid.(cyc).(ch), 'c05_discharge')
            data = RPT_VQ_grid.(cyc).(ch).c05_discharge;
            if strcmp(cyc, 'cyc0')
                plot(data.Q, data.V_grid, 'Color', ch_colors(ch_idx,:), 'LineWidth', 1.5);
            else
                plot(data.Q, data.V_grid, 'Color', [ch_colors(ch_idx,:) 0.3], 'LineWidth', 0.8);
            end
        end
    end
end
xl = xlim;
for i = 1:num_segments
    V_max = Global_V_bounds_dch(i);
    V_min = Global_V_bounds_dch(i+1);
    patch([xl(1) xl(2) xl(2) xl(1)], [V_min V_min V_max V_max], ...
          seg_colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
for i = 1:length(Global_V_bounds_dch)
    V_line = Global_V_bounds_dch(i);
    line(xl, [V_line V_line], ...
         'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
end
xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 11);
title('Discharge Segment Zoom', 'FontSize', 12, 'FontWeight', 'bold');
% Determine precise zoom window
ylim([Global_V_bounds_dch(end)-0.01, Global_V_bounds_dch(1)+0.01]);
xlim([10 50]); % Approximate zoom-in range for capacity
grid on; hold off;

sgtitle('MasterRuler v3 Final Segments Rendered over All Channels & Cycles (0.5C)', 'FontSize', 16, 'FontWeight', 'bold');

% Save Section 2 figure
currentScriptPath = mfilename('fullpath');
[ver02Dir, ~, ~] = fileparts(currentScriptPath);
saveDir = fullfile(ver02Dir, 'MasterRulerVisual');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
saveas(fig, fullfile(saveDir, 'MasterRuler_v3_VQ_Segments.fig'));

%% ========================================================================
% Section 3: dQ/dV Curves by C-rate with Segment Coloring
% ========================================================================
disp('Generating dQ/dV visualizations by C-rate...');

% Peak window = segment boundary range (phase transition tracking)
% Charge: [V_bounds_chg(1), V_bounds_chg(end)]
% Discharge: [min(V_bounds_dch), max(V_bounds_dch)]
ma_window = 21;

target_crates = {'c01', 'c05', 'c1', 'c2', 'c3'};
crate_labels  = {'0.1C', '0.5C', '1C', '2C', '3C'};
n_crates = length(target_crates);

fig_dqdv = figure('Position', [30, 30, 2200, 800], 'Name', 'dQ/dV by C-rate with Segment Colors');

for r = 1:n_crates
    f_chg = [target_crates{r} '_charge'];
    f_dch = [target_crates{r} '_discharge'];

    % --- Row 1: Charge dQ/dV ---
    subplot(2, n_crates, r); hold on;
    % Draw segment patches first (background)
    yl_tmp = [0 500]; % temporary y-range, will adjust after plotting
    for s = 1:num_segments
        patch([Global_V_bounds_chg(s) Global_V_bounds_chg(s+1) Global_V_bounds_chg(s+1) Global_V_bounds_chg(s)], ...
              [yl_tmp(1) yl_tmp(1) yl_tmp(2) yl_tmp(2)], ...
              seg_colors(s,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
    end


    % Plot dQ/dV curves
    for ch_idx = 1:length(channels)
        ch = channels{ch_idx};
        for cyc_idx = 1:length(cyc_fields)
            cyc = cyc_fields{cyc_idx};
            if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && ...
               isfield(RPT_VQ_grid.(cyc).(ch), f_chg)
                data = RPT_VQ_grid.(cyc).(ch).(f_chg);
                [V_u, uid] = unique(data.V_grid);
                Q_u = data.Q(uid);
                if V_u(1) > V_u(end), V_u = flipud(V_u); Q_u = flipud(Q_u); end
                dV = gradient(V_u); dQ_g = gradient(Q_u);
                dV(dV == 0) = NaN;
                dQdV = dQ_g ./ dV; dQdV(isinf(dQdV)|isnan(dQdV)) = 0;
                dQdV = movmean(dQdV, ma_window);
                if strcmp(cyc, 'cyc1000')
                    plot(V_u, dQdV, 'Color', ch_colors(ch_idx,:), 'LineWidth', 1.2);
                    % findpeaks on cyc0 curves
                    [pks, locs] = findpeaks(dQdV, 'MinPeakProminence', 5);
                    if ~isempty(pks)
                        plot(V_u(locs), pks, 'r.', 'MarkerSize', 12, 'HandleVisibility','on');
                    end
                else
                    plot(V_u, dQdV, 'Color', [ch_colors(ch_idx,:) 0.2], 'LineWidth', 0.5, 'HandleVisibility','off');
                end
            end
        end
    end
    title(sprintf('Charge %s', crate_labels{r}), 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('V'); ylabel('dQ/dV (Ah/V)');
    xlim([3.0 4.2]); grid on;
    % Auto-adjust y to data range
    yl = ylim; ylim([0 yl(2)*1.05]);
    hold off;

    % --- Row 2: Discharge dQ/dV ---
    subplot(2, n_crates, n_crates + r); hold on;
    yl_tmp = [0 500];
    % Discharge segment patches (V_bounds_dch is descending, sort ascending for patch)
    V_dch_sorted = sort(Global_V_bounds_dch);
    for s = 1:(length(V_dch_sorted)-1)
        patch([V_dch_sorted(s) V_dch_sorted(s+1) V_dch_sorted(s+1) V_dch_sorted(s)], ...
              [yl_tmp(1) yl_tmp(1) yl_tmp(2) yl_tmp(2)], ...
              seg_colors(s,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
    end


    for ch_idx = 1:length(channels)
        ch = channels{ch_idx};
        for cyc_idx = 1:length(cyc_fields)
            cyc = cyc_fields{cyc_idx};
            if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && ...
               isfield(RPT_VQ_grid.(cyc).(ch), f_dch)
                data = RPT_VQ_grid.(cyc).(ch).(f_dch);
                [V_u, uid] = unique(data.V_grid);
                Q_u = data.Q(uid);
                if V_u(1) > V_u(end), V_u = flipud(V_u); Q_u = flipud(Q_u); end
                dV = gradient(V_u); dQ_g = gradient(Q_u);
                dV(dV == 0) = NaN;
                dQdV = dQ_g ./ dV; dQdV(isinf(dQdV)|isnan(dQdV)) = 0;
                dQdV = movmean(dQdV, ma_window);
                dQdV_abs = abs(dQdV);
                if strcmp(cyc, 'cyc0')
                    plot(V_u, dQdV_abs, 'Color', ch_colors(ch_idx,:), 'LineWidth', 1.2);
                    % findpeaks on cyc0 curves
                    [pks, locs] = findpeaks(dQdV_abs, 'MinPeakProminence', 5);
                    if ~isempty(pks)
                        plot(V_u(locs), pks, 'r.', 'MarkerSize', 12, 'HandleVisibility','off');
                    end
                else
                    plot(V_u, dQdV_abs, 'Color', [ch_colors(ch_idx,:) 0.2], 'LineWidth', 0.5, 'HandleVisibility','off');
                end
            end
        end
    end
    title(sprintf('Discharge %s', crate_labels{r}), 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('V'); ylabel('|dQ/dV| (Ah/V)');
    xlim([3.0 4.2]); grid on;
    yl = ylim; ylim([0 yl(2)*1.05]);
    hold off;
end

sgtitle({'dQ/dV Curves by C-rate with MasterRuler v3 Segments', ...
    sprintf('Peak Window = Segment Range: Chg[%.3f–%.3fV], Dch[%.3f–%.3fV]', ...
    Global_V_bounds_chg(1), Global_V_bounds_chg(end), min(Global_V_bounds_dch), max(Global_V_bounds_dch))}, ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig_dqdv, fullfile(saveDir, 'MasterRuler_v3_dQdV_byCrate.fig'));

%% ========================================================================
% Section 4-6: Feature Trend Visualization (19 Features)
% ========================================================================
disp('Loading Feature Matrix for trend visualization...');
featFile = fullfile(ver02Dir, 'RPT_FeatureLabelExtractor', 'Feature_Matrix_ver02.mat');
if ~exist(featFile, 'file')
    warning('Feature_Matrix_ver02.mat not found. Skipping feature trend plots.');
else

load(featFile, 'FeatureTable_ver02', 'feature_names');
cells_feat  = unique(FeatureTable_ver02.CellID);
crates      = {'c01', 'c05', 'c1', 'c2', 'c3'};
num_crates  = length(crates);
num_seg     = 5;
color_cells = lines(length(cells_feat));

%% --- Section 4: dQ Segment Capacity (Charge & Discharge) ---
fig_dQ_chg = figure('Name','Feature: dQ_chg Segments','Position',[25,25,1800,900],'Visible','off');
sgtitle('Charge dQ by Segment × C-rate (MasterRuler v3)','FontSize',14,'FontWeight','bold');
fig_dQ_dch = figure('Name','Feature: dQ_dch Segments','Position',[75,75,1800,900],'Visible','off');
sgtitle('Discharge dQ by Segment × C-rate (MasterRuler v3)','FontSize',14,'FontWeight','bold');

for seg = 1:num_seg
    for r = 1:num_crates
        plot_idx = (seg-1)*num_crates + r;
        idx_cr = strcmp(FeatureTable_ver02.CrateLabel, crates{r});
        F_sub  = FeatureTable_ver02(idx_cr, :);

        chg_col = sprintf('dQ_chg_S%d', seg);
        dch_col = sprintf('dQ_dch_S%d', seg);

        set(0,'CurrentFigure', fig_dQ_chg);
        subplot(num_seg, num_crates, plot_idx); hold on; grid on;
        for i = 1:length(cells_feat)
            ci = strcmp(F_sub.CellID, cells_feat{i});
            plot(F_sub.Cycle(ci), F_sub.(chg_col)(ci), '-o', 'Color',color_cells(i,:), 'MarkerSize',4, 'DisplayName',cells_feat{i});
        end
        if seg==1, title(crates{r},'FontSize',10); end
        if r==1,   ylabel(sprintf('Seg%d (Ah)',seg)); end

        set(0,'CurrentFigure', fig_dQ_dch);
        subplot(num_seg, num_crates, plot_idx); hold on; grid on;
        for i = 1:length(cells_feat)
            ci = strcmp(F_sub.CellID, cells_feat{i});
            plot(F_sub.Cycle(ci), F_sub.(dch_col)(ci), '-s', 'Color',color_cells(i,:), 'MarkerSize',4, 'DisplayName',cells_feat{i});
        end
        if seg==1, title(crates{r},'FontSize',10); end
        if r==1,   ylabel(sprintf('Seg%d (Ah)',seg)); end
    end
end
saveas(fig_dQ_chg, fullfile(saveDir, 'Feature_dQ_chg_Segments.fig'));
saveas(fig_dQ_dch, fullfile(saveDir, 'Feature_dQ_dch_Segments.fig'));

%% --- Section 5: Peak Features (Height, Area, Position) ---
peak_pairs = { ...
    'Peak_H_chg', 'Peak_H_dch', 'Peak Height (dQ/dV)', 'max(dQ/dV)'; ...
    'Peak_A_chg', 'Peak_A_dch', 'Peak Area (\DeltaQ)',  'Ah'; ...
    'Peak_Pos_chg','Peak_Pos_dch','Peak Position (V)',    'V' };

for p = 1:size(peak_pairs,1)
    col_chg   = peak_pairs{p,1};
    col_dch   = peak_pairs{p,2};
    fig_title = peak_pairs{p,3};
    y_label   = peak_pairs{p,4};

    fig_pk = figure('Name', sprintf('Feature: %s', fig_title), ...
        'Position', [50+p*30, 50+p*30, 1400, 500], 'Visible','off');
    sgtitle(fig_title, 'FontSize',14, 'FontWeight','bold');

    for r = 1:num_crates
        idx_cr = strcmp(FeatureTable_ver02.CrateLabel, crates{r});
        F_sub  = FeatureTable_ver02(idx_cr, :);

        subplot(2, num_crates, r); hold on; grid on;
        for i = 1:length(cells_feat)
            ci = strcmp(F_sub.CellID, cells_feat{i});
            plot(F_sub.Cycle(ci), F_sub.(col_chg)(ci), '-o', 'Color',color_cells(i,:), 'MarkerSize',4, 'DisplayName',cells_feat{i});
        end
        title(sprintf('Chg — %s', crates{r}),'FontSize',9);
        ylabel(y_label); xlabel('Cycle');

        subplot(2, num_crates, num_crates+r); hold on; grid on;
        for i = 1:length(cells_feat)
            ci = strcmp(F_sub.CellID, cells_feat{i});
            plot(F_sub.Cycle(ci), F_sub.(col_dch)(ci), '-s', 'Color',color_cells(i,:), 'MarkerSize',4, 'DisplayName',cells_feat{i});
        end
        title(sprintf('Dch — %s', crates{r}),'FontSize',9);
        ylabel(y_label); xlabel('Cycle');
    end
    tag = strrep(col_chg, '_chg', '');
    saveas(fig_pk, fullfile(saveDir, sprintf('Feature_%s.fig', tag)));
end

%% --- Section 6: Scalar Features (Energy_dch, C_eff, T_avg) ---
scalar_feats   = {'Energy_dch', 'C_eff_chg', 'C_eff_dch', 'T_chg_avg', 'T_dch_avg'};
scalar_ylabels = {'Energy_{dch} (Wh)', 'C_{eff,chg}', 'C_{eff,dch}', 'T_{chg} (°C)', 'T_{dch} (°C)'};

% Filter to only existing columns
valid_mask = ismember(scalar_feats, FeatureTable_ver02.Properties.VariableNames);
scalar_feats   = scalar_feats(valid_mask);
scalar_ylabels = scalar_ylabels(valid_mask);

for sf = 1:length(scalar_feats)
    col = scalar_feats{sf};
    fig_sf = figure('Name', sprintf('Feature: %s', col), ...
        'Position', [100+sf*30, 300, 1200, 400], 'Visible','off');
    sgtitle(scalar_ylabels{sf}, 'FontSize',14, 'FontWeight','bold');

    for r = 1:num_crates
        idx_cr = strcmp(FeatureTable_ver02.CrateLabel, crates{r});
        F_sub  = FeatureTable_ver02(idx_cr, :);

        subplot(1, num_crates, r); hold on; grid on;
        for i = 1:length(cells_feat)
            ci = strcmp(F_sub.CellID, cells_feat{i});
            plot(F_sub.Cycle(ci), F_sub.(col)(ci), '-o', 'Color',color_cells(i,:), 'MarkerSize',4, 'DisplayName',cells_feat{i});
        end
        title(crates{r},'FontSize',10); xlabel('Cycle'); ylabel(scalar_ylabels{sf});
    end
    saveas(fig_sf, fullfile(saveDir, sprintf('Feature_%s.fig', col)));
end

end % end of feature matrix existence check

fprintf('\n=== MasterRuler Visualization Complete ===\n');
fprintf('  All figures saved to: %s\n', saveDir);
