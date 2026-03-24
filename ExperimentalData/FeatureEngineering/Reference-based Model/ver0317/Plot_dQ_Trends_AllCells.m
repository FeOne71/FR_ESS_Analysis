% Plot_dQ_Trends_AllCells.m
% 좌: Δ(ΔQ) vs Fresh Cell (세그먼트별, 사이클별)
% 우: 전 세그먼트 ΔQ 꺽은선 그래프 (x=사이클, y=ΔQ, 세그먼트별 색상)
% 채널별(Ch09~Ch16), C-rate별 저장

clear; clc; close all;

%% 1. Load Data
baseDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData';
verDir = fullfile(baseDir, 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
vqGridFile = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');

g = load(vqGridFile, 'RPT_VQ_grid');
grid = g.RPT_VQ_grid;

% Load Master Ruler
mr = load(fullfile(verDir, 'MasterRulers_PseudoOCV.mat'));
V_bounds = mr.MasterRuler_ver0317.V_bounds_chg;
num_segs = length(V_bounds) - 1;

%% 2. Output directory
visDir = fullfile(verDir, 'Visualization');
if ~exist(visDir, 'dir'), mkdir(visDir); end

%% 3. Configuration
conditions = {'c01', 'c05', 'c1', 'c2', 'c3'};
cond_labels = {'0.1C', '0.5C', '1.0C', '2.0C', '3.0C'};

% All cells
all_cells = {'Ch09','Ch10','Ch11','Ch12','Ch13','Ch14','Ch15','Ch16'};

% All cycles
cycles_all = fieldnames(grid);
cycles_all = sort(cycles_all(startsWith(cycles_all, 'cyc')));

% 12-color palette for segments (distinct colors)
seg_cmap = [
    0.90 0.10 0.10;  % Seg01 red
    0.00 0.45 0.74;  % Seg02 blue
    0.47 0.67 0.19;  % Seg03 green
    0.93 0.69 0.13;  % Seg04 orange
    0.49 0.18 0.56;  % Seg05 purple
    0.30 0.75 0.93;  % Seg06 cyan
    0.64 0.08 0.18;  % Seg07 dark red
    0.07 0.62 0.45;  % Seg08 teal
    0.85 0.33 0.10;  % Seg09 dark orange
    0.40 0.40 0.40;  % Seg10 gray
    0.72 0.27 1.00;  % Seg11 violet
    0.06 0.31 0.55;  % Seg12 navy
];

%% 4. Generate per Cell, per C-rate
for ci = 1:length(conditions)
    cond = conditions{ci};
    cond_label = cond_labels{ci};
    chg_field = [cond '_charge'];
    
    for ch_idx = 1:length(all_cells)
        target_cell = all_cells{ch_idx};
        
        % Collect dQ per cycle per segment
        cyc_nums = [];
        dQ_mat = [];
        
        for k = 1:length(cycles_all)
            cycStr = cycles_all{k};
            if ~isfield(grid.(cycStr), target_cell), continue; end
            cellData = grid.(cycStr).(target_cell);
            if ~isfield(cellData, chg_field), continue; end
            
            data = cellData.(chg_field);
            V = data.V_raw;
            Q = data.Q_raw;
            
            cyc_num = str2double(replace(cycStr, 'cyc', ''));
            cyc_nums(end+1) = cyc_num;
            
            row = zeros(1, num_segs);
            for s = 1:num_segs
                idx = V >= V_bounds(s) & V <= V_bounds(s+1);
                if any(idx)
                    row(s) = abs(max(Q(idx)) - min(Q(idx)));
                else
                    row(s) = NaN;
                end
            end
            dQ_mat(end+1, :) = row;
        end
        
        if isempty(cyc_nums) || length(cyc_nums) < 2
            continue;
        end
        
        % Sort by cycle
        [cyc_nums, sortIdx] = sort(cyc_nums);
        dQ_mat = dQ_mat(sortIdx, :);
        
        % Reference: first cycle
        dQ_ref = dQ_mat(1, :);
        delta_dQ = dQ_mat - dQ_ref;
        
        %% FIGURE
        fig = figure('Name', sprintf('%s_%s_%s', target_cell, cond, 'dQ_Trend'), ...
                     'Position', [30, 80, 1500, 550], 'Visible', 'off');
        
        %% LEFT: Δ(ΔQ) vs Segments per cycle
        subplot(1,2,1); hold on;
        
        % Pink background for positive (capacity INCREASE)
        y_max = max(delta_dQ(:), [], 'omitnan');
        y_min = min(delta_dQ(:), [], 'omitnan');
        if y_max > 0
            patch([0.5, num_segs+0.5, num_segs+0.5, 0.5], ...
                  [0, 0, y_max*1.3+0.3, y_max*1.3+0.3], ...
                  [1, 0.85, 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
            text(num_segs/2, y_max*1.15+0.2, 'Capacity INCREASE (Peak Shift Anomaly)', ...
                 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.8,0,0], 'FontWeight', 'bold');
        end
        
        yline(0, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        
        cyc_colors = jet(length(cyc_nums));
        for r = 1:length(cyc_nums)
            plot(1:num_segs, delta_dQ(r,:), '-o', 'Color', cyc_colors(r,:), ...
                 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', cyc_colors(r,:), ...
                 'DisplayName', sprintf('Cycle %d', cyc_nums(r)));
        end
        
        title(sprintf('NCM %s (%s) \\Delta(\\DeltaQ) vs Fresh Cell', target_cell, cond_label), ...
              'FontSize', 13, 'FontWeight', 'bold');
        xlabel('Voltage Segments (1 to 12)', 'FontSize', 11);
        ylabel('\Delta(\DeltaQ) vs Fresh Cell (Ah)', 'FontSize', 11);
        set(gca, 'XTick', 1:num_segs);
        xlim([0.5, num_segs+0.5]);
        grid on;
        legend('Location', 'southwest', 'FontSize', 7);
        hold off;
        
        %% RIGHT: All segments ΔQ line graph (x=cycle, y=ΔQ)
        subplot(1,2,2); hold on;
        
        markers = {'o','s','d','^','v','>','<','p','h','+','x','*'};
        for s = 1:num_segs
            mk = markers{mod(s-1, length(markers))+1};
            plot(cyc_nums, dQ_mat(:, s), ['-' mk], ...
                 'Color', seg_cmap(s,:), 'LineWidth', 1.5, ...
                 'MarkerSize', 5, 'MarkerFaceColor', seg_cmap(s,:), ...
                 'DisplayName', sprintf('Seg%d [%.2f~%.2fV]', s, V_bounds(s), V_bounds(s+1)));
        end
        
        title(sprintf('NCM %s (%s) \\DeltaQ per Segment over Cycles', target_cell, cond_label), ...
              'FontSize', 13, 'FontWeight', 'bold');
        xlabel('Cycle Number', 'FontSize', 11);
        ylabel('\DeltaQ (Ah)', 'FontSize', 11);
        grid on;
        legend('Location', 'eastoutside', 'FontSize', 7);
        hold off;
        
        %% Save per channel
        cellDir = fullfile(visDir, target_cell);
        if ~exist(cellDir, 'dir'), mkdir(cellDir); end
        
        fname = sprintf('dQ_Trend_%s_%s', target_cell, cond);
        saveas(fig, fullfile(cellDir, [fname '.png']), 'png');
        savefig(fig, fullfile(cellDir, [fname '.fig']));
        close(fig);
        fprintf('Saved: %s/%s.png\n', target_cell, fname);
    end
end

fprintf('\n=== All dQ Trend plots saved to: %s ===\n', visDir);
