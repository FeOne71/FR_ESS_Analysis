% Report_Vis_QV_Curves.m
% Generates business-style overlaid Q-V curves for 2021-2025 field data
clear; clc; close all;

%% 1. Configuration
base_dir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
rf_script = fullfile(base_dir, 'ML_RandomForest', 'RPT_Field_Estimation_RF.m');

fprintf('Running Field Estimation to extract data in memory...\n');
run(rf_script);

% Re-define paths because rf_script clears the workspace
base_dir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
save_dir = fullfile(base_dir, 'Report', 'Figures');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

yr_strs = fieldnames(Results_Ev);
n_yrs = length(yr_strs);

colors = lines(n_yrs); 
colors = colors * 0.85; % Enhance for business style

fig = figure('Name', 'Q-V Curves Overlay', 'Position', [100, 100, 1200, 600], 'Color', 'w');

%% 2. Plotting loop
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
                % Reconstruct Q from dQ/dV: Q = int(dQ/dV * dV)
                dV = 0.001; % fixed grid size used in extraction
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

%% 3. Master Title and Saving
sgtitle('Field Evaluation Data: Overlaid Q-V Curves (2021-2025)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);

save_path_f = fullfile(save_dir, 'Report_Vis_QV_Curves.fig');
save_path_p = fullfile(save_dir, 'Report_Vis_QV_Curves.png');
saveas(fig, save_path_f);
saveas(fig, save_path_p);

fprintf('Saved business-style Q-V overalys to:\n  %s\n', save_path_p);
