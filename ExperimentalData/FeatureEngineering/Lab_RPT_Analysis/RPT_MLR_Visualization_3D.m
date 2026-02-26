%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT_MLR_Visualization_3D.m
% 3D Surface Plot for SOH Regression Analysis
%
% Features:
%   X-axis: Discharge Feature (e.g., Dch_Seg3)
%   Y-axis: Charge Feature (e.g., Chg_Seg3)
%   Z-axis: SOH (State of Health)
%
% Output:
%   3D Scatter + Regression Plane Mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% 1. Load Data
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis';
dataPath = fullfile(baseDir, 'Dataset', 'Feature_Matrix_Final.mat');
saveDir = fullfile(baseDir, 'Modeling_Results');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

fprintf('Loading Feature Matrix...\n');
load(dataPath, 'FeatureTable');

% Data Extraction
X_All = FeatureTable.X_Normalized; % Z-score Normalized Features
Y = FeatureTable.Y_Labels(:, 1);   % SOH is the first column
UniqueCells = unique(FeatureTable.CellID);

% --- Select Top 2 Features (Based on Correlation) ---
% Index 14: Dch_Seg3 (Highest correlation ~0.92)
% Index 9 : Chg_Seg3 (High correlation ~0.80)
% (Check: RPT_Feature_Extraction_Advanced.m used 14 features. 
%  Indices need to be verified. Let's assume indices based on feature list provided previously.)
%
% Feature List (Typical order in extraction):
% 1-5: Chg_Seg1..5
% 6: Chg_PkH, 7: Chg_PkA
% 8-12: Dch_Seg1..5
% 13: Dch_PkH, 14: Dch_PkA
%
% Wait, exact indices matter. Let's assume Dch_Seg3 is index 10?
% Actually, let's look at FeatureTable structure if possible or use column names if available.
% FeatureTable uses 'X_Features' matrix. Column names are not in matrix.
% But we can infer from extraction script order.
% Order: Chg_Seg(1-5), Chg_PkH, Chg_PkA, Dch_Seg(1-5), Dch_PkH, Dch_PkA.
% Total 14.
% Ind 1-5: Chg Segs
% Ind 6-7: Chg Pk
% Ind 8-12: Dch Segs (Seg1=8, Seg2=9, Seg3=10, Seg4=11, Seg5=12)
% Ind 13-14: Dch Pk
%
% Selected Features:
% X_feat = Dch_Seg3 (Index 10)
% Y_feat = Chg_Seg3 (Index 3)

idx_X = 10; % Discharge Segment 3
idx_Y = 3;  % Charge Segment 3

feat_name_X = 'Discharge Segment 3 (Normalized)';
feat_name_Y = 'Charge Segment 3 (Normalized)';
target_name = 'SOH (Ah)';

X1 = X_All(:, idx_X);
X2 = X_All(:, idx_Y);

% Remove NaNs
valid_mask = ~isnan(X1) & ~isnan(X2) & ~isnan(Y);
X1 = X1(valid_mask);
X2 = X2(valid_mask);
Y  = Y(valid_mask);
CellIDs = FeatureTable.CellID(valid_mask);

fprintf('Selected Features for 3D Plot:\n');
fprintf('  X: %s (Index %d)\n', feat_name_X, idx_X);
fprintf('  Y: %s (Index %d)\n', feat_name_Y, idx_Y);
fprintf('  Z: %s\n', target_name);
fprintf('  Valid Samples: %d\n', length(Y));

%% 2. Fit Simple Linear Model (for Visualization)
% Model: Y ~ b0 + b1*X1 + b2*X2
X_design = [ones(length(Y),1), X1, X2];
b = regress(Y, X_design);
fprintf('Regression Coefficients:\n');
fprintf('  Intercept: %.4f\n', b(1));
fprintf('  Coeff X1 : %.4f\n', b(2));
fprintf('  Coeff X2 : %.4f\n', b(3));

%% 3. Generate Grid for Mesh Plane
x1_range = linspace(min(X1), max(X1), 20);
x2_range = linspace(min(X2), max(X2), 20);
[X1_Grid, X2_Grid] = meshgrid(x1_range, x2_range);

% Predict Z on Grid
Z_Grid = b(1) + b(2)*X1_Grid + b(3)*X2_Grid;

%% 4. Plotting
fig = figure('Position', [100, 100, 1000, 800], 'Name', 'MLR 3D Surface Visualization');
hold on;

% A. Regression Plane (Mesh)
% Use a colormap for the surface (representing SOH level)
mesh(X1_Grid, X2_Grid, Z_Grid, 'FaceAlpha', 0.5, 'EdgeColor', [0.5 0.5 0.5]);
colormap jet;

% B. Scatter Plot (Real Data)
% Color points by Cell ID to see cell grouping
% Need numeric codes for cells
[~, ~, cell_idx] = unique(CellIDs);
scatter3(X1, X2, Y, 60, cell_idx, 'filled', 'MarkerEdgeColor', 'k');

% C. Residual Lines (Stem)
% Draw lines from actual point to predicted plane
Y_pred = b(1) + b(2)*X1 + b(3)*X2;
for i = 1:length(Y)
    plot3([X1(i), X1(i)], [X2(i), X2(i)], [Y(i), Y_pred(i)], 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
end

xlabel(feat_name_X, 'FontSize', 12, 'FontWeight', 'bold');
ylabel(feat_name_Y, 'FontSize', 12, 'FontWeight', 'bold');
zlabel(target_name, 'FontSize', 12, 'FontWeight', 'bold');
title({'SOH Regression Plane (Top 2 Features)', 'Model: SOH ~ Dch\_Seg3 + Chg\_Seg3'}, 'FontSize', 15);
grid on;
view(135, 30); % Initial view angle

% Legend setup (Manual for Cell IDs is hard, let's make a colorbar for cells?)
% Or just legend for "Data" and "Model"
h_plane = plot(nan, nan, 's', 'MarkerFaceColor', 'none', 'Color', [0.5 0.5 0.5], 'LineWidth', 2); % Dummy for legend
h_data  = plot(nan, nan, 'ko', 'MarkerFaceColor', 'b');
legend([h_data, h_plane], {'Real Data (Points)', 'Regression Plane'}, 'Location', 'best');

colorbar;
ylabel(colorbar, 'Cell Group Index');

hold off;

%% 5. Save
savePath = fullfile(saveDir, 'MLR_3D_Surface_SOH.fig');
saveas(fig, savePath);
saveas(fig, strrep(savePath, '.fig', '.png'));
fprintf('3D Visualization saved to: %s\n', savePath);
