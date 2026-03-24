% ML_00_Correlation_Analysis.m
% Analyzes the Pearson Correlation Coefficient (PCC) between 
% Fragmented Capacity Features (dQ) and Static Capacity (SOH) across C-rates.
% Goal: Prove that the optimal feature segment shifts dynamically with C-rate.

clear; clc; close all; warning off;

%% 1. Directories & Load Data
baseDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData';
verDir = fullfile(baseDir, 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
fmFilePath = fullfile(verDir, 'FeatureMatrix_ver0317.mat');

if ~exist(fmFilePath, 'file')
    error('FeatureMatrix_ver0317.mat not found! Run RPT_FeatureExtractor first.');
end

load(fmFilePath, 'FM');

%% 2. Define Features and Target Label
targetLabel = 'Static_Capacity';
allVars = FM.Properties.VariableNames;
feat_c = allVars(startsWith(allVars, 'dQ_c_'));
feat_d = allVars(startsWith(allVars, 'dQ_d_'));
num_seg_chg = length(feat_c);
num_seg_dch = length(feat_d);
features = [feat_c, feat_d];

%% 3. Split by Condition (C-rate) and Calculate PCC Matrix
conditions = unique(FM.Condition); 
num_conds = length(conditions);
num_feats = length(features);

% Matrix to store PCC results: [num_conds x num_feats]
PCC_matrix = zeros(num_conds, num_feats);

for c = 1:num_conds
    cond = conditions(c);
    % Filter data by exact C-rate group
    idx = FM.Condition == cond;
    subFM = FM(idx, :);
    
    Y = subFM.(targetLabel);
    
    for f = 1:num_feats
        X = subFM.(features{f});
        
        % NaN-aware: exclude NaN rows, then check std
        valid = ~isnan(X) & ~isnan(Y);
        if sum(valid) >= 3 && std(X(valid)) > 0 && std(Y(valid)) > 0
            R = corrcoef(X(valid), Y(valid));
            PCC_matrix(c, f) = R(1,2);
        else
            PCC_matrix(c, f) = NaN;
        end
    end
end

%% 4. Visualization 1: PCC Grouped Bar Chart
fig_bar = figure('Name', 'PCC Correlation: dQ vs SOH', 'Position', [100, 100, 1400, 800]);

% -- Charge dQ Subplot --
subplot(2,1,1);
bar(PCC_matrix(:, 1:num_seg_chg)');
title('PCC of Charge Segments (dQ\_c) vs SOH across C-rates', 'FontSize', 14, 'FontWeight', 'bold');
xlabel(sprintf('Segment Number (1~%d)', num_seg_chg), 'FontSize', 12);
ylabel('Pearson Correlation (r)', 'FontSize', 12);
ylim([-1 1]); grid on;
set(gca, 'XTick', 1:num_seg_chg);
legend(string(conditions), 'Location', 'southoutside', 'Orientation', 'horizontal');

% -- Discharge dQ Subplot --
subplot(2,1,2);
bar(PCC_matrix(:, num_seg_chg+1:end)');
title('PCC of Discharge Segments (dQ\_d) vs SOH across C-rates', 'FontSize', 14, 'FontWeight', 'bold');
xlabel(sprintf('Segment Number (1~%d)', num_seg_dch), 'FontSize', 12);
ylabel('Pearson Correlation (r)', 'FontSize', 12);
ylim([-1 1]); grid on;
set(gca, 'XTick', 1:num_seg_dch);
legend(string(conditions), 'Location', 'southoutside', 'Orientation', 'horizontal');

saveas(fig_bar, fullfile(verDir, 'ML_00_PCC_BarChart.fig'));

%% 5. Visualization 2: Correlation Heatmap Matrix
% Color-coded matrix showing global linear correlation patterns
fig_heat = figure('Name', 'Correlation Heatmap', 'Position', [100, 100, 1200, 400]);
h = heatmap(features, string(conditions), PCC_matrix);
h.Title = 'Pearson Correlation: Fragmented Features vs SOH';
h.XLabel = 'Fragmented Segment Features (Charge & Discharge)';
h.YLabel = 'Tested C-rate Condition';
h.Colormap = jet;
h.ColorLimits = [-1 1];

    saveas(fig_heat, fullfile(verDir, 'ML_00_PCC_Heatmap.fig'));
    
    %% 6. Export PCC Results for Analysis
    pcc_table = array2table(PCC_matrix, 'RowNames', string(conditions), 'VariableNames', string(features));
    writetable(pcc_table, fullfile(verDir, 'ML_00_PCC_Results.csv'), 'WriteRowNames', true);
    fprintf('Correlation Analysis Complete. Visualizations saved as .fig in ver0317 directory.\n');
    fprintf('Exact PCC matrix saved to ML_00_PCC_Results.csv\n');
