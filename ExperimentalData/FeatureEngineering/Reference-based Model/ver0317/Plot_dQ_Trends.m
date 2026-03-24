% Plot_dQ_Trends.m
% C-rate별 세그먼트/사이클 dQ 변화 시각화
% 좌: Δ(ΔQ) vs Fresh Cell (세그먼트별, 사이클별)
% 우: Diverging Degradation Paths (이상 vs 단조 세그먼트)

clear; clc; close all;

%% 1. Load Feature Matrix
baseDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData';
verDir = fullfile(baseDir, 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
d = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat'));
FM = d.FM;

%% 2. Configuration
conditions = {'c01', 'c05', 'c1', 'c2', 'c3'};
cond_labels = {'0.1C', '0.5C', '1.0C', '2.0C', '3.0C'};
num_segs = 12;

% dQ feature column names (charge only)
feat_c = cell(1, num_segs);
for i = 1:num_segs
    feat_c{i} = sprintf('dQ_c_%02d', i);
end

% Representative cell for left plot
target_cell = 'Ch09';

% Cycle list (as stored in FM)
all_cycles = unique(FM.CycleNum);
cycle_colors = jet(length(all_cycles));

%% 3. Generate plots per C-rate
for ci = 1:length(conditions)
    cond = conditions{ci};
    cond_label = cond_labels{ci};
    
    % Filter for target cell & this condition
    idx = FM.CellID == target_cell & FM.Condition == cond;
    subFM = FM(idx, :);
    
    if height(subFM) < 2
        fprintf('Skipping %s: insufficient data\n', cond);
        continue;
    end
    
    % Sort by cycle
    subFM = sortrows(subFM, 'CycleNum');
    cyc_nums = double(subFM.CycleNum);
    
    % Extract dQ matrix: [num_cycles x num_segs]
    dQ_mat = zeros(height(subFM), num_segs);
    for s = 1:num_segs
        dQ_mat(:, s) = subFM.(feat_c{s});
    end
    
    % Reference: Cycle 0 (first row)
    dQ_ref = dQ_mat(1, :);
    
    % Δ(ΔQ) = dQ(cycle) - dQ(cycle0)
    delta_dQ = dQ_mat - dQ_ref;
    
    %% LEFT PLOT: Δ(ΔQ) across segments per cycle
    fig = figure('Name', sprintf('dQ Trend - %s - %s', target_cell, cond_label), ...
                 'Position', [50, 100, 1400, 550]);
    
    subplot(1,2,1); hold on;
    
    % Pink background for positive region (capacity INCREASE)
    patch([0.5, num_segs+0.5, num_segs+0.5, 0.5], [0, 0, max(delta_dQ(:))*1.3+0.5, max(delta_dQ(:))*1.3+0.5], ...
          [1, 0.85, 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
    
    % Zero reference line
    yline(0, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % Plot each cycle
    cyc_colors_local = jet(length(cyc_nums));
    for r = 1:length(cyc_nums)
        plot(1:num_segs, delta_dQ(r,:), '-o', 'Color', cyc_colors_local(r,:), ...
             'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', cyc_colors_local(r,:), ...
             'DisplayName', sprintf('Cycle %d', cyc_nums(r)));
    end
    
    % Annotate
    text(num_segs/2, max(delta_dQ(:))*1.1+0.3, 'Capacity INCREASE (Peak Shift Anomaly)', ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', [0.8, 0, 0], 'FontWeight', 'bold');
    
    title(sprintf('NCM Cell %s (%s) Segment Capacity Change over Cycles', target_cell, cond_label), 'FontSize', 13, 'FontWeight', 'bold');
    xlabel('Voltage Segments (1 to 12)', 'FontSize', 12);
    ylabel('\Delta(\DeltaQ) vs Fresh Cell (Ah)', 'FontSize', 12);
    set(gca, 'XTick', 1:num_segs);
    xlim([0.5, num_segs+0.5]);
    grid on;
    legend('Location', 'southwest', 'FontSize', 8);
    hold off;
    
    %% RIGHT PLOT: Diverging Degradation Paths
    subplot(1,2,2); hold on;
    
    % Find segments with strongest positive delta (anomaly), strongest negative, and mixed
    mean_delta_last = delta_dQ(end, :);
    [~, idx_inc] = max(mean_delta_last);   % Most increasing segment
    [~, idx_dec] = min(mean_delta_last);   % Most decreasing segment
    
    % Mixed: segment with middle behavior
    abs_delta = abs(mean_delta_last);
    abs_delta(idx_inc) = Inf;
    abs_delta(idx_dec) = Inf;
    [~, idx_mix] = min(abs_delta);
    
    % Plot absolute dQ over cycles
    plot(cyc_nums, dQ_mat(:, idx_dec), '-o', 'Color', [0.9, 0.1, 0.1], 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', [0.9,0.1,0.1], ...
         'DisplayName', sprintf('Seg %d (\\DeltaQ Decreases steadily)', idx_dec));
    plot(cyc_nums, dQ_mat(:, idx_inc), '-s', 'Color', [0.1, 0.1, 0.9], 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', [0.1,0.1,0.9], ...
         'DisplayName', sprintf('Seg %d (Anomaly: \\DeltaQ Increases!!)', idx_inc));
    plot(cyc_nums, dQ_mat(:, idx_mix), '-^', 'Color', [0.1, 0.8, 0.8], 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', [0.1,0.8,0.8], ...
         'DisplayName', sprintf('Seg %d (Mixed behavior)', idx_mix));
    
    title('Diverging Degradation Paths (Monotonic vs Anomaly)', 'FontSize', 13, 'FontWeight', 'bold');
    xlabel('Cycle', 'FontSize', 12);
    ylabel('Absolute \DeltaQ (Ah)', 'FontSize', 12);
    grid on;
    legend('Location', 'best', 'FontSize', 9);
    hold off;
    
    % Save
    saveas(fig, fullfile(verDir, sprintf('dQ_Trend_%s_%s.png', target_cell, cond)), 'png');
    savefig(fig, fullfile(verDir, sprintf('dQ_Trend_%s_%s.fig', target_cell, cond)));
    fprintf('Saved: dQ_Trend_%s_%s.png/fig\n', target_cell, cond);
end

fprintf('\n=== All dQ Trend plots generated! ===\n');
