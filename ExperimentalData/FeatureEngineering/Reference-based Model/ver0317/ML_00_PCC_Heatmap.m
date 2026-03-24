% ML_00_PCC_Heatmap.m
% PCC 히트맵: ΔQ vs SOH(Static_Capacity) 상관계수 + p-value 유의성 표시
% X축: 세그먼트별 전압 범위
% Y축: Charging Rate (All, 0.1C ~ 3.0C)
% 유의성: * p<0.05, ** p<0.01, *** p<0.001

clear; clc; close all;

%% 1. Load
baseDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData';
verDir = fullfile(baseDir, 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
d = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat'));
FM = d.FM;

mr = load(fullfile(verDir, 'MasterRulers_PseudoOCV.mat'));
V_chg = mr.MasterRuler_ver0317.V_bounds_chg;
V_dch = mr.MasterRuler_ver0317.V_bounds_dch;

visDir = fullfile(verDir, 'Visualization');
if ~exist(visDir, 'dir'), mkdir(visDir); end

%% 2. Feature & Label
num_segs = 12;
feat_c = cell(1, num_segs);
feat_d = cell(1, num_segs);
for i = 1:num_segs
    feat_c{i} = sprintf('dQ_c_%02d', i);
    feat_d{i} = sprintf('dQ_d_%02d', i);
end
features = [feat_c, feat_d];
num_feats = length(features);
targetLabel = 'Static_Capacity';

%% 3. X-axis labels (voltage ranges)
xlabels_chg = cell(1, num_segs);
xlabels_dch = cell(1, num_segs);
for i = 1:num_segs
    xlabels_chg{i} = sprintf('%.2f~%.2fV', V_chg(i), V_chg(i+1));
    xlabels_dch{i} = sprintf('%.2f~%.2fV', V_dch(i), V_dch(i+1));
end

%% 4. Y-axis: All + 5 C-rates
conditions = {'All', 'c01', 'c05', 'c1', 'c2', 'c3'};
ylabels = {'All', '0.1C', '0.5C', '1.0C', '2.0C', '3.0C'};
num_rows = length(conditions);

%% 5. Compute PCC + p-value matrix
PCC = NaN(num_rows, num_feats);
PVAL = NaN(num_rows, num_feats);

for r = 1:num_rows
    if strcmp(conditions{r}, 'All')
        subFM = FM;
    else
        subFM = FM(FM.Condition == conditions{r}, :);
    end
    Y = subFM.(targetLabel);
    
    for f = 1:num_feats
        X = subFM.(features{f});
        valid = ~isnan(X) & ~isnan(Y);
        n = sum(valid);
        if n >= 3 && std(X(valid)) > 0 && std(Y(valid)) > 0
            [R, P] = corrcoef(X(valid), Y(valid));
            PCC(r, f) = R(1,2);
            PVAL(r, f) = P(1,2);
        end
    end
end

%% 6. Field coverage
field_chg = [0 0 0 0 0 0 0 1 1 0 0 0];
field_dch = [0 0 0 0 0 0 0 0 0 0 0 0];

%% 7. Helper: significance string
sig_str = @(p) get_sig(p);

%% 8. Plot Heatmap - Charge
fig1 = figure('Position', [50, 100, 1200, 500], 'Name', 'PCC Heatmap - Charge');

PCC_chg = PCC(:, 1:num_segs);
PVAL_chg = PVAL(:, 1:num_segs);
imagesc(PCC_chg);
colormap(redblue(256));
caxis([-1, 1]);
cb = colorbar;
cb.Label.String = 'Pearson Correlation Coefficient (r)';
cb.Label.FontSize = 11;

set(gca, 'XTick', 1:num_segs, 'XTickLabel', xlabels_chg, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:num_rows, 'YTickLabel', ylabels);
xlabel('Charge Segment Voltage Range', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Charging Rate (C-rate)', 'FontSize', 12, 'FontWeight', 'bold');
title('PCC(\DeltaQ_{charge}, SOH) per Segment and Charging Rate', 'FontSize', 14, 'FontWeight', 'bold');

% Annotate PCC values + significance
for r = 1:num_rows
    for c = 1:num_segs
        val = PCC_chg(r, c);
        pv = PVAL_chg(r, c);
        if isnan(val)
            text(c, r, 'NaN', 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.5 0.5 0.5]);
        else
            if abs(val) > 0.5
                txtColor = 'w';
            else
                txtColor = 'k';
            end
            stars = get_sig(pv);
            label = sprintf('%.2f%s', val, stars);
            text(c, r, label, 'HorizontalAlignment', 'center', 'FontSize', 7, ...
                 'FontWeight', 'bold', 'Color', txtColor);
        end
    end
end

% Field coverage markers
for s = 1:num_segs
    if field_chg(s)
        text(s, 0.55, char(9679), 'HorizontalAlignment', 'center', 'FontSize', 12, ...
             'Color', [0 0.7 0], 'FontWeight', 'bold', 'Clipping', 'off');
    end
end

set(gca, 'FontSize', 10);
saveas(fig1, fullfile(visDir, 'PCC_Heatmap_Charge.png'), 'png');
savefig(fig1, fullfile(visDir, 'PCC_Heatmap_Charge.fig'));
fprintf('Saved: PCC_Heatmap_Charge.png\n');

%% 9. Plot Heatmap - Discharge
fig2 = figure('Position', [50, 100, 1200, 500], 'Name', 'PCC Heatmap - Discharge');

PCC_dch = PCC(:, num_segs+1:end);
PVAL_dch = PVAL(:, num_segs+1:end);
imagesc(PCC_dch);
colormap(redblue(256));
caxis([-1, 1]);
cb = colorbar;
cb.Label.String = 'Pearson Correlation Coefficient (r)';
cb.Label.FontSize = 11;

set(gca, 'XTick', 1:num_segs, 'XTickLabel', xlabels_dch, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:num_rows, 'YTickLabel', ylabels);
xlabel('Discharge Segment Voltage Range', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Discharging Rate (C-rate)', 'FontSize', 12, 'FontWeight', 'bold');
title('PCC(\DeltaQ_{discharge}, SOH) per Segment and Discharging Rate', 'FontSize', 14, 'FontWeight', 'bold');

for r = 1:num_rows
    for c = 1:num_segs
        val = PCC_dch(r, c);
        pv = PVAL_dch(r, c);
        if isnan(val)
            text(c, r, 'NaN', 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.5 0.5 0.5]);
        else
            if abs(val) > 0.5
                txtColor = 'w';
            else
                txtColor = 'k';
            end
            stars = get_sig(pv);
            label = sprintf('%.2f%s', val, stars);
            text(c, r, label, 'HorizontalAlignment', 'center', 'FontSize', 7, ...
                 'FontWeight', 'bold', 'Color', txtColor);
        end
    end
end

set(gca, 'FontSize', 10);
saveas(fig2, fullfile(visDir, 'PCC_Heatmap_Discharge.png'), 'png');
savefig(fig2, fullfile(visDir, 'PCC_Heatmap_Discharge.fig'));
fprintf('Saved: PCC_Heatmap_Discharge.png\n');

%% 10. Save PCC + p-value to CSV
rowNames = ylabels';
T_pcc = array2table(PCC, 'VariableNames', features, 'RowNames', rowNames);
writetable(T_pcc, fullfile(verDir, 'ML_00_PCC_Results.csv'), 'WriteRowNames', true);

T_pval = array2table(PVAL, 'VariableNames', features, 'RowNames', rowNames);
writetable(T_pval, fullfile(verDir, 'ML_00_PCC_Pvalues.csv'), 'WriteRowNames', true);

fprintf('Saved: ML_00_PCC_Results.csv, ML_00_PCC_Pvalues.csv\n');
fprintf('\n=== PCC Heatmap generation complete! ===\n');

%% Helper Functions
function cmap = redblue(n)
    if nargin < 1, n = 256; end
    half = floor(n/2);
    r1 = linspace(0.2, 1, half)';
    g1 = linspace(0.2, 1, half)';
    b1 = linspace(0.8, 1, half)';
    r2 = linspace(1, 0.8, n-half)';
    g2 = linspace(1, 0.2, n-half)';
    b2 = linspace(1, 0.2, n-half)';
    cmap = [r1 g1 b1; r2 g2 b2];
end

function s = get_sig(p)
    if isnan(p)
        s = '';
    elseif p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = '';
    end
end
