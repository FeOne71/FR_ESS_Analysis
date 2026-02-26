
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Cell Trace Comparison: Actual vs MLR vs No-PCA RF
% Purpose: Visualize how each model tracks SOH/LLI/LAM per cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% 1. Load Data
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis';
path_main  = fullfile(baseDir, 'Modeling_Results', 'Modeling_Results.mat');
path_nopca = fullfile(baseDir, 'Modeling_Results', 'Modeling_Results_NoPCA.mat');
saveDir    = fullfile(baseDir, 'Modeling_Results', 'Paper_Figures_Comparison');

if ~exist(path_main, 'file') || ~exist(path_nopca, 'file')
    error('Modeling result files missing. Run both RPT_Modeling.m and RPT_Modeling_NoPCA.m first.');
end
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

% Load MLR results
main = load(path_main);
% Load No-PCA RF results
nopca = load(path_nopca);

label_names = main.label_names;
n_labels = length(label_names);
unique_cells = unique(main.cellIDs);
n_cells = length(unique_cells);

%% 2. Reconstruct Full Prediction Vectors from K-Fold
% Since MLR and No-PCA RF might have different cv_indices, 
% we must reconstruct based on their respective indices.

% 2.1 Reconstruct MLR Predictions
Y_pred_MLR_full = zeros(size(main.Y));
for fold = 1:5
    idx = (main.cv_indices == fold);
    Y_pred_MLR_full(idx, :) = main.kfold_Y_pred_MLR{fold};
end

% 2.2 Reconstruct No-PCA RF Predictions
Y_pred_RF_full = zeros(size(nopca.Y));
for fold = 1:5
    idx = (nopca.cv_indices == fold);
    Y_pred_RF_full(idx, :) = nopca.kfold_Y_pred_RF{fold};
end

% Note: main.Y and nopca.Y should be identical (ordered same as FeatureTable)
Y_actual = main.Y;
cellIDs  = main.cellIDs;

%% 3. Generate Trace Plots for each Label
for lbl = 1:n_labels
    fprintf('Generating Trace Plot for %s...\n', label_names{lbl});
    
    fig = figure('Position', [50, 50, 1600, 900], 'Name', ['Cell Trace Comparison: ' label_names{lbl}]);
    
    for c = 1:n_cells
        subplot(2, 4, c);
        hold on;
        
        cell_mask = strcmp(cellIDs, unique_cells{c});
        
        y_act = Y_actual(cell_mask, lbl);
        y_mlr = Y_pred_MLR_full(cell_mask, lbl);
        y_rf  = Y_pred_RF_full(cell_mask, lbl);
        
        % Plotting
        xticks = 1:length(y_act);
        plot(xticks, y_act, '-ko', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Actual');
        plot(xticks, y_mlr, '--bv', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', 'MLR');
        plot(xticks, y_rf,  '--r^', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', 'RF');
        
        % Calculate cell-specific R2 (for title)
        r2_cell = 1 - sum((y_act - y_rf).^2)/sum((y_act - mean(y_act)).^2);
        
        title(sprintf('%s (R²=%.2f)', unique_cells{c}, r2_cell), 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Sample Index', 'FontSize', 10);
        ylabel([label_names{lbl} ' (%)'], 'FontSize', 10);
        grid on;
        
        if c == 1, legend('Location', 'best', 'FontSize', 8); end
    end
    
    sgtitle(sprintf('5-Fold CV: %s Trace per Cell (Actual vs MLR vs RF)', label_names{lbl}), ...
        'FontSize', 18, 'FontWeight', 'bold');
    
    saveas(fig, fullfile(saveDir, ['Cell_Trace_Comparison_' label_names{lbl} '.fig']));
end

fprintf('\nCell Trace Visualization Completed.\n');
fprintf('Figures saved in: %s\n', saveDir);
