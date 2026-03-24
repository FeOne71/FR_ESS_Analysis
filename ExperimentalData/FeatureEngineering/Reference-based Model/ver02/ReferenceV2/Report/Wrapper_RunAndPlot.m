% Wrapper to run field estimation and then plot Q-V
clear; clc; close all;

% 1. Run the main field estimation script to populate Results_Ev in base workspace
run('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest\RPT_Field_Estimation_RF.m');

fprintf('\n=== Generating Business-Style Q-V Overlay ===\n');

% 2. Extract and Plot from Results_Ev
base_dir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
save_dir = fullfile(base_dir, 'Report', 'Figures');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

yr_strs = fieldnames(Results_Ev);
n_yrs = length(yr_strs);

colors = lines(n_yrs); 
colors = colors * 0.85;

fig = figure('Name', 'Q-V Curves Overlay', 'Position', [100, 100, 1200, 600], 'Color', 'w');

for j = 1:2
    subplot(1, 2, j);
    hold on; grid on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'GridAlpha', 0.15, 'Box', 'on');
    
    if j == 1
        title('Charge Q-V Curves', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);
        xlabel('Capacity (Ah)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);
        ylabel('Voltage (V)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);
        xlim([0, 50]);
        ylim([3.6, 4.0]);
    else
        title('Discharge Q-V Curves', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);
        xlabel('Capacity (Ah)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);
        ylabel('Voltage (V)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);
        xlim([0, 50]);
        ylim([3.6, 3.9]);
    end
    
    h_lines = zeros(1, n_yrs);
    
    for k = 1:n_yrs
        yr = yr_strs{k};
        if j == 1
            V = Results_Ev.(yr).Curves.V_chg;
            dQ = Results_Ev.(yr).Curves.dQdV_chg;
            if ~isempty(V)
                dV = 0.001; 
                Q = cumsum(dQ * dV);
                h_lines(k) = plot(Q, V, 'Color', colors(k,:), 'LineWidth', 2.5, 'DisplayName', yr);
            end
        else
            V = Results_Ev.(yr).Curves.V_dch;
            dQ = Results_Ev.(yr).Curves.dQdV_dch;
            if ~isempty(V)
                dV = 0.001; 
                Q = cumsum(dQ * dV);
                h_lines(k) = plot(Q, V, 'Color', colors(k,:), 'LineWidth', 2.5, 'DisplayName', yr);
            end
        end
    end
    
    ax = gca;
    ax.XColor = [0.3 0.3 0.3];
    ax.YColor = [0.3 0.3 0.3];
    
    valid_lines = h_lines(h_lines~=0);
    if ~isempty(valid_lines)
        leg = legend(valid_lines, 'Location', 'best', 'FontSize', 11, 'TextColor', [0.2 0.2 0.2]);
        leg.ItemTokenSize = [30, 18];
        leg.Box = 'off';
    end
end

sgtitle('Field Evaluation Data: Overlaid Q-V Curves (2021-2025)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);

save_path_p = fullfile(save_dir, 'Report_Vis_QV_Curves.png');
saveas(fig, save_path_p);
fprintf('\nSaved %s\n', save_path_p);
