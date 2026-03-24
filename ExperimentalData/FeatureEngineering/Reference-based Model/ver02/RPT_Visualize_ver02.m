%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Feature & Label Visualizer (ver02 Final)
% - 21 Features: dQ_seg(10), Peak_H/A/Pos(6), Energy_dch, C_eff_chg/dch, T_avg(2)
% - 4  Labels:   SOH, LLI, LAM, SOP_dch_10s (max)
% - Saves all figures to 'ver02/FeatureLabelVisual'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
warning off;

%% ========================================================================
% Section 1: Configuration & Load
% ========================================================================
currentScriptPath = mfilename('fullpath');
[ver02Dir, ~, ~] = fileparts(currentScriptPath);
dataDir = fullfile(ver02Dir, 'RPT_FeatureLabelExtractor');
visDir  = fullfile(ver02Dir, 'FeatureLabelVisual');
if ~exist(visDir, 'dir'), mkdir(visDir); end

featFile  = fullfile(dataDir, 'Feature_Matrix_ver02.mat');
labelFile = fullfile(dataDir, 'Label_Matrix_ver02.mat');
if ~exist(featFile, 'file') || ~exist(labelFile, 'file')
    error('Feature or Label Matrix not found. Run extractors first.');
end

load(featFile,  'FeatureTable_ver02', 'feature_names');
load(labelFile, 'LabelTable_ver02', 'label_names');

cells     = unique(FeatureTable_ver02.CellID);
cycles    = unique(FeatureTable_ver02.Cycle);
crates    = {'c01', 'c05', 'c1', 'c2', 'c3'};
num_crates    = length(crates);
num_segments  = 5;

color_cells  = lines(length(cells));
color_cycles = jet(max(length(cycles), 2));
markers = {'o', 's', '^', 'd', 'v', 'p', 'h', '*', 'x', '+'};

set(0, 'DefaultFigureWindowStyle', 'docked');

%% ========================================================================
% Section 2: Feature — dQ Segmented Capacity (Charge & Discharge)
% ========================================================================
fig_dQ_chg = figure('Name', 'Feature: dQ_chg', 'Position', [25,25,1800,900], 'Visible','off');
sgtitle('Charge Capacity dQ\_chg by Segment × C-rate', 'FontSize',16, 'FontWeight','bold');

fig_dQ_dch = figure('Name', 'Feature: dQ_dch', 'Position', [75,75,1800,900], 'Visible','off');
sgtitle('Discharge Capacity dQ\_dch by Segment × C-rate', 'FontSize',16, 'FontWeight','bold');

for seg = 1:num_segments
    % Per-segment ylim: same segment across all C-rates
    chg_col = sprintf('dQ_chg_S%d', seg);
    dch_col = sprintf('dQ_dch_S%d', seg);
    v_c = FeatureTable_ver02.(chg_col); v_c = v_c(~isnan(v_c));
    v_d = FeatureTable_ver02.(dch_col); v_d = v_d(~isnan(v_d));
    m_c = (max(v_c)-min(v_c))*0.05; m_d = (max(v_d)-min(v_d))*0.05;
    ylim_seg_chg = [min(v_c)-m_c, max(v_c)+m_c];
    ylim_seg_dch = [min(v_d)-m_d, max(v_d)+m_d];

    for r = 1:num_crates
        plot_idx = (seg-1)*num_crates + r;
        idx_cr = strcmp(FeatureTable_ver02.CrateLabel, crates{r});
        F_sub  = FeatureTable_ver02(idx_cr, :);

        set(0,'CurrentFigure', fig_dQ_chg);
        subplot(num_segments, num_crates, plot_idx); hold on; grid on;
        for i = 1:length(cells)
            ci = strcmp(F_sub.CellID, cells{i});
            plot(F_sub.Cycle(ci), F_sub.(chg_col)(ci), '-o', 'Color', color_cells(i,:), 'DisplayName', cells{i});
        end
        ylim(ylim_seg_chg);
        if seg==1, title(crates{r},'FontSize',10); end
        if r==1,   ylabel(sprintf('Seg %d (Ah)',seg)); end

        set(0,'CurrentFigure', fig_dQ_dch);
        subplot(num_segments, num_crates, plot_idx); hold on; grid on;
        for i = 1:length(cells)
            ci = strcmp(F_sub.CellID, cells{i});
            plot(F_sub.Cycle(ci), F_sub.(dch_col)(ci), '-s', 'Color', color_cells(i,:), 'DisplayName', cells{i});
        end
        ylim(ylim_seg_dch);
        if seg==1, title(crates{r},'FontSize',10); end
        if r==1,   ylabel(sprintf('Seg %d (Ah)',seg)); end
    end
end
saveas(fig_dQ_chg, fullfile(visDir, 'Feature_dQ_chg.fig'));

saveas(fig_dQ_dch, fullfile(visDir, 'Feature_dQ_dch.fig'));


%% ========================================================================
% Section 3: Feature — Peak Height, Area, Position (Charge & Discharge)
% ========================================================================
peak_pairs = { ...
    'Peak_H_chg', 'Peak_H_dch', 'Peak Height (dQ/dV)', 'max(dQ/dV)'; ...
    'Peak_A_chg', 'Peak_A_dch', 'Peak Area (\DeltaQ)', 'Ah'; ...
    'Peak_Pos_chg','Peak_Pos_dch','Peak Position (V)', 'V' };

for p = 1:size(peak_pairs,1)
    col_chg  = peak_pairs{p,1};
    col_dch  = peak_pairs{p,2};
    fig_title = peak_pairs{p,3};
    y_label   = peak_pairs{p,4};

    % Separate ylim for charge row and discharge row
    v_chg = FeatureTable_ver02.(col_chg); v_chg = v_chg(~isnan(v_chg));
    v_dch = FeatureTable_ver02.(col_dch); v_dch = v_dch(~isnan(v_dch));
    m_chg = (max(v_chg)-min(v_chg))*0.05; m_dch = (max(v_dch)-min(v_dch))*0.05;
    ylim_pk_chg = [min(v_chg)-m_chg, max(v_chg)+m_chg];
    ylim_pk_dch = [min(v_dch)-m_dch, max(v_dch)+m_dch];

    fig_pk = figure('Name', sprintf('Feature: %s', fig_title), ...
        'Position', [50+p*30, 50+p*30, 1400, 500], 'Visible','off');
    sgtitle(fig_title, 'FontSize',14, 'FontWeight','bold');

    for r = 1:num_crates
        idx_cr = strcmp(FeatureTable_ver02.CrateLabel, crates{r});
        F_sub  = FeatureTable_ver02(idx_cr, :);

        subplot(2, num_crates, r); hold on; grid on;
        for i = 1:length(cells)
            ci = strcmp(F_sub.CellID, cells{i});
            plot(F_sub.Cycle(ci), F_sub.(col_chg)(ci), '-o', 'Color', color_cells(i,:), 'DisplayName', cells{i});
        end
        ylim(ylim_pk_chg);
        title(sprintf('Chg — %s', crates{r}),'FontSize',9);
        ylabel(y_label); xlabel('Cycle');

        subplot(2, num_crates, num_crates+r); hold on; grid on;
        for i = 1:length(cells)
            ci = strcmp(F_sub.CellID, cells{i});
            plot(F_sub.Cycle(ci), F_sub.(col_dch)(ci), '-s', 'Color', color_cells(i,:), 'DisplayName', cells{i});
        end
        ylim(ylim_pk_dch);
        title(sprintf('Dch — %s', crates{r}),'FontSize',9);
        ylabel(y_label); xlabel('Cycle');
    end
    tag = strrep(col_chg, '_chg', '');
    saveas(fig_pk, fullfile(visDir, sprintf('Feature_%s.fig', tag)));
    
end

%% ========================================================================
% Section 4: Feature — C_eff & T_avg
% ========================================================================
scalar_feats = {'C_eff_chg', 'C_eff_dch', 'T_chg_avg', 'T_dch_avg'};
scalar_ylabels = {'C_{eff,chg} = I_{chg}/Q_0', 'C_{eff,dch} = I_{dch}/Q_0', 'T_{chg} (°C)', 'T_{dch} (°C)'};

for sf = 1:length(scalar_feats)
    col = scalar_feats{sf};

    % Global ylim across all C-rates for this feature
    all_sf_vals = FeatureTable_ver02.(col);
    all_sf_vals = all_sf_vals(~isnan(all_sf_vals));
    sf_margin = (max(all_sf_vals) - min(all_sf_vals)) * 0.05;
    ylim_sf = [min(all_sf_vals) - sf_margin, max(all_sf_vals) + sf_margin];

    fig_sf = figure('Name', sprintf('Feature: %s', col), ...
        'Position', [100+sf*30, 300, 1200, 400], 'Visible','off');
    sgtitle(scalar_ylabels{sf}, 'FontSize',14, 'FontWeight','bold');

    for r = 1:num_crates
        idx_cr = strcmp(FeatureTable_ver02.CrateLabel, crates{r});
        F_sub  = FeatureTable_ver02(idx_cr, :);

        subplot(1, num_crates, r); hold on; grid on;
        for i = 1:length(cells)
            ci = strcmp(F_sub.CellID, cells{i});
            plot(F_sub.Cycle(ci), F_sub.(col)(ci), '-o', 'Color', color_cells(i,:), 'DisplayName', cells{i});
        end
        ylim(ylim_sf);
        title(crates{r},'FontSize',10); xlabel('Cycle'); ylabel(scalar_ylabels{sf});
    end
    saveas(fig_sf, fullfile(visDir, sprintf('Feature_%s.fig', col)));
    
end

%% ========================================================================
% Section 5-7: Label — SOH / LLI / LAM / SOP (통합 Subplot)
% ========================================================================
L_c1 = LabelTable_ver02(strcmp(LabelTable_ver02.CrateLabel, 'c1'), :);

fig_labels = figure('Name', 'Labels: SOH / LLI / LAM / SOP', ...
    'Position', [50, 200, 2000, 450]);
sgtitle('Battery Diagnostic Labels vs Cycle', 'FontSize', 16, 'FontWeight', 'bold');

label_cols   = {'SOH',   'LLI',   'LAM',   'SOP_dch_10s'};
label_titles = {'SOH (%)', 'LLI (%)', 'LAM (%)', 'SOP_{dch,10s} (W)'};
label_ylims  = {[70 101], [0 8], [0 50], []};

for d = 1:4
    subplot(1, 4, d); hold on; grid on; box on;
    col = label_cols{d};
    
    if ismember(col, L_c1.Properties.VariableNames)
        for i = 1:length(cells)
            ci = strcmp(L_c1.CellID, cells{i});
            mk = markers{mod(i-1, length(markers)) + 1};
            plot(L_c1.Cycle(ci), L_c1.(col)(ci), ['-' mk], ...
                'Color', color_cells(i,:), 'LineWidth', 1.5, ...
                'MarkerSize', 6, 'DisplayName', cells{i});
        end
        title(label_titles{d}, 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Cycle', 'FontSize', 10);
        ylabel(label_titles{d}, 'FontSize', 10);
        if ~isempty(label_ylims{d}), ylim(label_ylims{d}); end
        legend('Location', 'best', 'FontSize', 8);
    else
        title(sprintf('%s (N/A)', label_titles{d}), 'FontSize', 12);
        text(0.5, 0.5, 'Column not found', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
    end
end

saveas(fig_labels, fullfile(visDir, 'Labels_All.fig'));
saveas(fig_labels, fullfile(visDir, 'Labels_All.png'));


%% ========================================================================
% Section 8: Feature × Label — Correlation Heatmap
% ========================================================================
% Merge Feature + Label on CellID / Cycle / CrateLabel
X_all = table2array(FeatureTable_ver02(:, feature_names));
Y_all = table2array(LabelTable_ver02(:, label_names));

[~, ia, ib] = intersect( ...
    strcat(FeatureTable_ver02.CellID, num2str(FeatureTable_ver02.Cycle), FeatureTable_ver02.CrateLabel), ...
    strcat(LabelTable_ver02.CellID,   num2str(LabelTable_ver02.Cycle),   LabelTable_ver02.CrateLabel));

X_mat = X_all(ia,:);
Y_mat = Y_all(ib,:);
nan_rows = any(isnan([X_mat, Y_mat]), 2);
X_mat(nan_rows,:) = [];
Y_mat(nan_rows,:) = [];

% Feature groups
chg_feat_names = {'dQ_chg_S1','dQ_chg_S2','dQ_chg_S3','dQ_chg_S4','dQ_chg_S5', ...
                  'Peak_H_chg','Peak_A_chg','C_eff_chg'};
dch_feat_names = {'dQ_dch_S1','dQ_dch_S2','dQ_dch_S3','dQ_dch_S4','dQ_dch_S5', ...
                  'Peak_H_dch','Peak_A_dch','Energy_dch','C_eff_dch'};

chg_idx = cellfun(@(f) find(strcmp(feature_names,f), 1), chg_feat_names);
dch_idx = cellfun(@(f) find(strcmp(feature_names,f), 1), dch_feat_names);
% all_idx: union of chg+dch, excluding T_avg and C_eff
exclude_feats = {'T_chg_avg','T_dch_avg'};
all_feat_names = feature_names(~cellfun(@(f) ismember(f, exclude_feats), feature_names));
all_idx = cellfun(@(f) find(strcmp(feature_names,f), 1), all_feat_names);

% Colormap (Blue-White-Red)
n_cm = 64;
blue=[0.15 0.25 0.55]; wht=[1 1 1]; red=[0.7 0.15 0.15];
cmap_bwr = [ ...
    [linspace(blue(1),wht(1),n_cm)', linspace(blue(2),wht(2),n_cm)', linspace(blue(3),wht(3),n_cm)']; ...
    [linspace(wht(1),red(1),n_cm)',  linspace(wht(2),red(2),n_cm)',  linspace(wht(3),red(3),n_cm)']  ];

% Define 6 heatmap configurations: {feat_idx, feat_names, corr_type, tag, title_str}
hmap_configs = {
    chg_idx,  chg_feat_names,  'Pearson',  'PCC_chg',      'PCC — Charge Features x Labels';
    dch_idx,  dch_feat_names,  'Pearson',  'PCC_dch',      'PCC — Discharge Features x Labels';
    all_idx,  all_feat_names,  'Pearson',  'PCC_all',      'PCC — All Features x Labels';
    chg_idx,  chg_feat_names,  'Spearman', 'Spearman_chg', 'Spearman — Charge Features x Labels';
    dch_idx,  dch_feat_names,  'Spearman', 'Spearman_dch', 'Spearman — Discharge Features x Labels';
    all_idx,  all_feat_names,  'Spearman', 'Spearman_all', 'Spearman — All Features x Labels';
};

    % Get unique C-rates for individual heatmaps
    unique_crates = unique(LabelTable_ver02.CrateLabel);
    target_groups = [{'All'}; unique_crates];

    for g = 1:length(target_groups)
        grp = target_groups{g};
        
        if strcmp(grp, 'All')
            % Use all rows
            X_grp = X_mat;
            Y_grp = Y_mat;
            grp_title = '[All C-rates]';
            grp_tag = 'All';
        else
            % Filter rows by specific C-rate
            % ib is the index vector that mapped original LabelTable to Y_mat
            % We need to find which rows in Y_mat belong to this crate.
            % An easier way: LabelTable_ver02.CrateLabel(ib) corresponds to Y_mat
            grp_idx = strcmp(LabelTable_ver02.CrateLabel(ib), grp);
            % But we removed nan_rows, so we must align indices:
            ib_valid = ib(~nan_rows);
            grp_idx = strcmp(LabelTable_ver02.CrateLabel(ib_valid), grp);
            
            X_grp = X_mat(grp_idx, :);
            Y_grp = Y_mat(grp_idx, :);
            grp_title = sprintf('[%s]', grp);
            grp_tag = grp;
        end
        
        if isempty(X_grp) || size(X_grp, 1) < 3
            fprintf('Skipping heatmap for %s (not enough data)\n', grp);
            continue;
        end

        for k = 1:size(hmap_configs,1)
            f_idx    = hmap_configs{k,1};
            f_names  = hmap_configs{k,2};
            corr_type = hmap_configs{k,3};
            tag      = hmap_configs{k,4};
            ttl      = sprintf('%s %s', grp_title, hmap_configs{k,5});

            n_feat_k = length(f_idx);
            n_label  = length(label_names);

            if strcmp(grp_tag, 'All')
                vis_option = 'on';
            else
                vis_option = 'off';
            end

            fig_h = figure('Name', ttl, 'Position',[50,50, 200+n_label*80, 150+n_feat_k*40], ...
                'Visible', vis_option, 'Color','w');
            plot_corr_heatmap(X_grp(:,f_idx), Y_grp, f_names, label_names, cmap_bwr, ttl, corr_type);
            saveas(fig_h, fullfile(visDir, sprintf('CorrHeatmap_%s_%s.fig', grp_tag, tag)));
            if strcmp(vis_option, 'off')
               close(fig_h);
            end
        end
    end

fprintf('\n=== Visualization Complete ===\n');
fprintf('  All figures saved to: %s\n', visDir);

%% ========================================================================
% Local Function: plot_corr_heatmap
% ========================================================================
function plot_corr_heatmap(X_block, Y_block, x_names, y_names, cmap_bwr, ttl, corr_type)
    % Compute feature × label correlation + p-values
    n_x = size(X_block,2);
    n_y = size(Y_block,2);
    XY  = [X_block, Y_block];
    [R_full, P_full] = corr(XY, 'Rows','pairwise', 'Type', corr_type);
    R_block = R_full(1:n_x, n_x+1:end);  % [n_feat × n_label]
    P_block = P_full(1:n_x, n_x+1:end);

    ax = gca;
    imagesc(R_block);
    colormap(ax.Parent, cmap_bwr); caxis([-1 1]);
    cb = colorbar('eastoutside');
    cb.Label.String = sprintf('%s r', corr_type);
    cb.Label.FontSize = 11;

    title({ttl, '(*p<0.05, **p<0.01, ***p<0.001)'}, 'FontSize', 11, 'FontWeight', 'bold');

    set(ax, 'XTick', 1:n_y, 'XTickLabel', y_names, 'FontSize', 11, ...
        'FontWeight', 'bold', 'XAxisLocation', 'top');
    set(ax, 'YTick', 1:n_x, 'YTickLabel', x_names, 'FontSize', 10, ...
        'FontWeight', 'bold');
    set(ax, 'TickLabelInterpreter', 'none');
    xtickangle(30);

    for i = 1:n_x
        for j = 1:n_y
            val = R_block(i,j);
            p   = P_block(i,j);

            % Significance stars
            if     p < 0.001, stars = '***';
            elseif p < 0.01,  stars = '**';
            elseif p < 0.05,  stars = '*';
            else,             stars = '';
            end

            label_str = sprintf('%.2f%s', val, stars);
            fc = 'k'; if abs(val) > 0.55, fc = 'w'; end
            text(j, i, label_str, ...
                'HorizontalAlignment', 'center', 'FontSize', 9, ...
                'FontWeight', 'bold', 'Color', fc);
        end
    end
end
