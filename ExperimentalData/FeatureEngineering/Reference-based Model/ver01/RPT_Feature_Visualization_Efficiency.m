%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Feature Visualization (Efficiency)
% - Visualization of dQ_seg Distribution (Stacked Bar)
% - Correlation between Efficiency (eta_U, eta_Wh) and SOH
% - C-rate sensitivity Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% ========================================================================
% Section 1: Configuration & Data Loading
% ========================================================================
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis';
path_physics = fullfile(baseDir, 'Dataset', 'Feature_Matrix_Physics.mat');
path_labels = fullfile(baseDir, 'Dataset', 'Feature_Matrix_Final.mat');
saveDir = fullfile(baseDir, 'Pipeline_Visualizations');

if ~exist(saveDir, 'dir'), mkdir(saveDir); end

fprintf('Loading Physics Features and Labels...\n');
load(path_physics, 'FeatureTable_Physics');
load(path_labels, 'FeatureTable'); % For SOH labels

% Match Physics Features with Labels (assuming same row order)
% Or join tables
FeatureTable_All = [FeatureTable_Physics, FeatureTable(:, 'Y_Labels')];

%% ========================================================================
% Section 2: Efficiency vs SOH Correlation
% ========================================================================
fprintf('Generating Efficiency vs SOH Correlation Plots...\n');
fig1 = figure('Position', [100, 100, 1400, 600]);

target_crates = unique(FeatureTable_All.CrateNum);
colors = lines(length(target_crates));

% Subplot 1: Voltage Efficiency (eta_U) vs SOH
subplot(1,2,1); hold on;
for r = 1:length(target_crates)
    c_num = target_crates(r);
    idx = (FeatureTable_All.CrateNum == c_num) & strcmp(FeatureTable_All.CellID, 'Ch09');
    
    soh = FeatureTable_All.Y_Labels(idx, 1);
    eta_u = FeatureTable_All.eta_U(idx);
    
    scatter(soh, eta_u, 60, colors(r,:), 'filled', 'DisplayName', sprintf('%.1fC', c_num));
    
    % Trend line
    p = polyfit(soh, eta_u, 1);
    f = polyval(p, soh);
    plot(soh, f, 'Color', colors(r,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
end
xlabel('SoH (Ah)', 'FontSize', 12); ylabel('\eta_U (V_{cha,avg} / V_{dis,avg})', 'FontSize', 12);
title('Voltage Efficiency vs. SOH (Ch09)', 'FontSize', 14); grid on; legend('Location', 'best');

% Subplot 2: Energy Efficiency (eta_Wh) vs SOH
subplot(1,2,2); hold on;
for r = 1:length(target_crates)
    c_num = target_crates(r);
    idx = (FeatureTable_All.CrateNum == c_num) & strcmp(FeatureTable_All.CellID, 'Ch09');
    
    soh = FeatureTable_All.Y_Labels(idx, 1);
    eta_wh = FeatureTable_All.eta_Wh(idx);
    
    scatter(soh, eta_wh, 40, colors(r,:), 'filled', 'DisplayName', sprintf('%.1fC', c_num));
    
    % Trend line
    p = polyfit(soh, eta_wh, 1);
    f = polyval(p, soh);
    plot(soh, f, 'Color', colors(r,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
end
xlabel('SoH (Ah)', 'FontSize', 12); ylabel('\eta_{Wh} (Energy Out / Energy In)', 'FontSize', 12);
title('Energy Efficiency vs. SOH (Ch09)', 'FontSize', 14); grid on;
saveas(fig1, fullfile(saveDir, 'Efficiency_vs_SOH.fig'));

%% ========================================================================
% Section 3: Segmented Capacity Distribution (Stacked Bar)
% ========================================================================
fprintf('Generating Segmented Capacity Distribution Plots...\n');
fig2 = figure('Position', [100, 100, 1400, 800]);

% Analyze 1C Charge
idx_1c = (FeatureTable_All.CrateNum == 1.0) & strcmp(FeatureTable_All.CellID, 'Ch09');
tbl_1c = FeatureTable_All(idx_1c, :);
cycles = tbl_1c.Cycle;
dq_seg = tbl_1c.dQ_chg_seg; % [N x 5]

subplot(2,1,1);
b = bar(cycles, dq_seg, 'stacked');
xlabel('Cycle'); ylabel('dQ_{seg} (Ah)'); title('Charge dQ Segment Distribution (1C, Ch09)');
legend({'Seg1', 'Seg2', 'Seg3', 'Seg4', 'Seg5'}, 'Location', 'eastoutside');
grid on;

% Analyze 1C Discharge
dq_seg_dch = tbl_1c.dQ_dch_seg; % [N x 5]
subplot(2,1,2);
bar(cycles, dq_seg_dch, 'stacked');
xlabel('Cycle'); ylabel('dQ_{seg} (Ah)'); title('Discharge dQ Segment Distribution (1C, Ch09)');
grid on;

saveas(fig2, fullfile(saveDir, 'Segmented_Capacity_Stacked.fig'));
fprintf('All Visualizations Generated Successfully.\n');
