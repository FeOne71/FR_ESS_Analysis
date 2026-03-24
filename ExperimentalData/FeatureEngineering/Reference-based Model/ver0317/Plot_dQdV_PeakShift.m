% Plot_dQdV_PeakShift.m
% dQ/dV 곡선을 사이클별로 오버레이하여 피크 쉬프트 및 ΔQ 증가 현상 직접 확인
% 상단: dQ/dV vs V (사이클별 오버레이, Master Ruler 경계 표시)
% 하단: 세그먼트별 ΔQ 사이클 트렌드 (증가/감소 확인)

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

%% 2. Configuration
target_cell = 'Ch09';
cycles = fieldnames(grid);
cycles = cycles(startsWith(cycles, 'cyc'));
cycles = sort(cycles);

conditions = {'c01', 'c05', 'c1', 'c2', 'c3'};
cond_labels = {'0.1C', '0.5C', '1.0C', '2.0C', '3.0C'};

% Color map for cycles
cyc_colors = jet(length(cycles));

%% 3. Generate per-C-rate plots
for ci = 1:length(conditions)
    cond = conditions{ci};
    cond_label = cond_labels{ci};
    chg_field = [cond '_charge'];
    
    fig = figure('Name', sprintf('dQ/dV Peak Shift - %s (%s)', target_cell, cond_label), ...
                 'Position', [50, 50, 1500, 800]);
    
    %% TOP: dQ/dV vs V overlay per cycle
    subplot(2,1,1); hold on;
    
    % Draw Master Ruler boundaries (vertical lines)
    for b = 1:length(V_bounds)
        xline(V_bounds(b), '--', sprintf('%.3f', V_bounds(b)), ...
               'Color', [0.6 0.6 0.6], 'Alpha', 0.5, 'LabelOrientation', 'aligned', ...
               'LabelVerticalAlignment', 'bottom', 'FontSize', 7, 'HandleVisibility', 'off');
    end
    
    % Shade segments where ΔQ increases (to be filled after analysis)
    dQ_per_cycle = zeros(length(cycles), num_segs);
    valid_cycles = false(length(cycles), 1);
    cyc_labels = {};
    
    for k = 1:length(cycles)
        cycStr = cycles{k};
        if ~isfield(grid.(cycStr), target_cell)
            continue;
        end
        cellData = grid.(cycStr).(target_cell);
        if ~isfield(cellData, chg_field)
            continue;
        end
        
        data = cellData.(chg_field);
        V = data.V_raw;
        Q = data.Q_raw;
        
        % Compute dQ/dV
        dV = diff(V);
        dQ = diff(Q);
        dQdV = dQ ./ dV;
        V_mid = (V(1:end-1) + V(2:end)) / 2;
        
        % Smooth dQ/dV for clean visualization
        dQdV_smooth = smoothdata(dQdV, 'sgolay', 51);
        
        % Plot
        plot(V_mid, dQdV_smooth, '-', 'Color', cyc_colors(k,:), 'LineWidth', 1.2, ...
             'DisplayName', strrep(cycStr, 'cyc', 'Cycle '));
        
        % Compute ΔQ per segment
        valid_cycles(k) = true;
        cyc_labels{end+1} = strrep(cycStr, 'cyc', '');
        for s = 1:num_segs
            idx = V >= V_bounds(s) & V <= V_bounds(s+1);
            if any(idx)
                dQ_per_cycle(k, s) = abs(max(Q(idx)) - min(Q(idx)));
            else
                dQ_per_cycle(k, s) = NaN;
            end
        end
    end
    
    title(sprintf('dQ/dV vs Voltage - %s (%s Charge) with Master Ruler Boundaries', target_cell, cond_label), ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Voltage (V)', 'FontSize', 12);
    ylabel('dQ/dV (Ah/V)', 'FontSize', 12);
    xlim([3.0 4.2]);
    % Limit y-axis to reasonable range
    yl = ylim;
    ylim([max(yl(1), -300), min(yl(2), 300)]);
    grid on;
    legend('Location', 'northwest', 'FontSize', 8, 'NumColumns', 2);
    hold off;
    
    %% BOTTOM: ΔQ trend per segment across cycles
    subplot(2,1,2); hold on;
    
    dQ_valid = dQ_per_cycle(valid_cycles, :);
    cyc_nums = cellfun(@str2double, cyc_labels);
    
    % Determine increasing vs decreasing: compare last to first
    delta_total = dQ_valid(end,:) - dQ_valid(1,:);
    
    seg_colors = zeros(num_segs, 3);
    for s = 1:num_segs
        if delta_total(s) > 0.05  % ΔQ INCREASES → anomaly (red, thick)
            seg_colors(s,:) = [0.9 0.1 0.1];
            lw = 2.5;
            marker = 's';
            label = sprintf('Seg%d [%.3f~%.3fV] ΔQ↑ ANOMALY', s, V_bounds(s), V_bounds(s+1));
        elseif delta_total(s) < -0.05  % ΔQ DECREASES → normal (blue)
            seg_colors(s,:) = [0.3 0.3 0.8];
            lw = 1.2;
            marker = 'o';
            label = sprintf('Seg%d [%.3f~%.3fV] ΔQ↓', s, V_bounds(s), V_bounds(s+1));
        else  % Mixed/flat
            seg_colors(s,:) = [0.5 0.5 0.5];
            lw = 1.0;
            marker = '.';
            label = sprintf('Seg%d [%.3f~%.3fV] flat', s, V_bounds(s), V_bounds(s+1));
        end
        
        plot(cyc_nums, dQ_valid(:, s), ['-' marker], 'Color', seg_colors(s,:), ...
             'LineWidth', lw, 'MarkerSize', 5, 'MarkerFaceColor', seg_colors(s,:), ...
             'DisplayName', label);
    end
    
    title(sprintf('Fragmented ΔQ per Segment over Cycles (%s, %s)', target_cell, cond_label), ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Cycle Number', 'FontSize', 12);
    ylabel('ΔQ (Ah)', 'FontSize', 12);
    grid on;
    legend('Location', 'eastoutside', 'FontSize', 8);
    hold off;
    
    % Print summary
    fprintf('\n=== %s (%s) Peak Shift Analysis ===\n', target_cell, cond_label);
    for s = 1:num_segs
        if delta_total(s) > 0.05
            fprintf('  ** Seg%02d [%.3f~%.3fV]: ΔQ INCREASES by +%.3f Ah (ANOMALY!)\n', ...
                    s, V_bounds(s), V_bounds(s+1), delta_total(s));
        elseif delta_total(s) < -0.05
            fprintf('     Seg%02d [%.3f~%.3fV]: ΔQ decreases by %.3f Ah\n', ...
                    s, V_bounds(s), V_bounds(s+1), delta_total(s));
        end
    end
    
    % Save
    saveas(fig, fullfile(verDir, sprintf('dQdV_PeakShift_%s_%s.png', target_cell, cond)), 'png');
    savefig(fig, fullfile(verDir, sprintf('dQdV_PeakShift_%s_%s.fig', target_cell, cond)));
    fprintf('  Saved: dQdV_PeakShift_%s_%s.png/fig\n', target_cell, cond);
end

fprintf('\n=== All dQ/dV Peak Shift plots generated! ===\n');
