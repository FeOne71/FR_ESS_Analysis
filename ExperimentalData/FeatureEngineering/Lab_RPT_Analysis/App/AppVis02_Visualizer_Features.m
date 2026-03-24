function App_Visualizer_Features(App_VQ_grid, FeatureTable, checkedNodes, saveDir)
% APP_VISUALIZER_FEATURES Visualizes extracted features and labels
%
% Reference: RPT_Label_Visualization.m
%
% Inputs:
%   App_VQ_grid: struct from App_DataLoader
%   FeatureTable: table from App_FeatureExtractor (CellID, Cycle, X_Features, Y_Labels, X_Normalized)
%   checkedNodes: cell array of checked Tree_3 node names (e.g. {'Equilibrium', 'Label'})
%   saveDir: path to save figures
%
% Figures:
%   [Always] Feature Heatmap, Label Summary (3-subplot)
%   [Equilibrium checked] OCV dQ/dV Peak Evolution (8-channel subplot)
%   [Label sub-items checked] Individual SOH/LLI/LAM trend plots

fprintf('--- App Visualizer: Features & Labels ---\n');

if nargin < 3 || isempty(checkedNodes), checkedNodes = {}; end
if nargin < 4 || isempty(saveDir)
    saveDir = fullfile(fileparts(mfilename('fullpath')), 'Results');
end
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

channels = unique(FeatureTable.CellID);
ch_colors = lines(length(channels));
marker_list = {'o', 's', 'd', '^', 'v', '>', '<', 'p'};

label_names = {'SOH (Capacity)', 'LLI', 'LAM'};
label_units = {'Capacity (Ah)', 'LLI (%)', 'LAM (%)'};

%% ===== Always Shown: Feature Distribution (Boxplot, Chg/Dch separated) =====
% Reference: RPT_Pipeline_Visualization.m Phase 4
fig1 = figure('Position', [50, 50, 1600, 900], 'Name', 'Feature Distribution');

X_raw = FeatureTable.X_Features;
X_norm = FeatureTable.X_Normalized;

feature_labels = {'Chg Seg1','Chg Seg2','Chg Seg3','Chg Seg4','Chg Seg5', ...
                  'Dch Seg1','Dch Seg2','Dch Seg3','Dch Seg4','Dch Seg5', ...
                  'Chg PkH','Chg PkA','Dch PkH','Dch PkA'};

% Before normalization
subplot(2,1,1);
boxplot(X_raw, 'Labels', feature_labels);
ylabel('Feature Value (Raw)', 'FontWeight', 'bold', 'FontSize', 12);
title('Feature Distribution (Raw)', 'FontSize', 13, 'FontWeight', 'bold');
grid on;

% After normalization (Z-score)
subplot(2,1,2);
boxplot(X_norm, 'Labels', feature_labels);
ylabel('Feature Value (Normalized)', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Features (Chg dQ 1-5, Dch dQ 1-5, Chg dQdV 2, Dch dQdV 2)', 'FontWeight', 'bold', 'FontSize', 11);
title('Feature Distribution (Z-score Normalized)', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
ylim([-4 4]);

sgtitle('Feature Distribution: Raw vs Normalized', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig1, fullfile(saveDir, 'Features_Distribution.fig'));

%% ===== Always Shown: Label Summary (3-subplot) =====
fig2 = figure('Position', [100, 100, 1500, 500], 'Name', 'Label Trends Summary');

for l = 1:3
    subplot(1, 3, l); hold on;
    
    for i = 1:length(channels)
        ch = channels{i};
        % Filter by C-rate (labels are C-rate independent, use 0.5C to avoid duplicates)
        if ismember('CrateNum', FeatureTable.Properties.VariableNames)
            idx = strcmp(FeatureTable.CellID, ch) & FeatureTable.CrateNum == 0.5;
        else
            idx = strcmp(FeatureTable.CellID, ch);
        end
        subT = FeatureTable(idx, :);
        if isempty(subT), continue; end
        [~, si] = sort(subT.Cycle);
        
        mk = marker_list{mod(i-1, length(marker_list))+1};
        plot(subT.Cycle(si), subT.Y_Labels(si, l), ['-', mk], ...
            'Color', ch_colors(i,:), 'LineWidth', 1.5, 'MarkerSize', 6, ...
            'MarkerFaceColor', ch_colors(i,:), 'DisplayName', ch);
    end
    
    title(label_names{l}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Cycle', 'FontWeight', 'bold');
    ylabel(label_units{l}, 'FontWeight', 'bold');
    grid on;
    if l == 1, legend('Location', 'best', 'FontSize', 8); end
end

sgtitle('Label Trends (SOH / LLI / LAM)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig2, fullfile(saveDir, 'Features_Label_Trends.fig'));

%% ===== Equilibrium Checked: OCV dQ/dV Peak Evolution =====
if any(strcmpi(checkedNodes, 'Equilibrium'))
    fprintf('  [Equilibrium] Generating OCV dQ/dV Peak Evolution...\n');
    
    cyc_fields = fieldnames(App_VQ_grid);
    cyc_colors = jet(length(cyc_fields));
    
    nCh = length(channels);
    nRows = ceil(nCh / 4);
    fig3 = figure('Position', [50, 50, 1600, 450*nRows], 'Name', 'OCV dQ/dV Peak Evolution');
    
    for i = 1:nCh
        ch = channels{i};
        subplot(nRows, min(nCh,4), i); hold on;
        
        for c = 1:length(cyc_fields)
            cyc = cyc_fields{c};
            if isfield(App_VQ_grid.(cyc), ch) && isfield(App_VQ_grid.(cyc).(ch), 'OCV_charge') && ...
               isfield(App_VQ_grid.(cyc).(ch).OCV_charge, 'V_grid')
                
                data = App_VQ_grid.(cyc).(ch).OCV_charge;
                [V_u, uid] = unique(data.V_grid);
                Q_u = data.Q(uid);
                
                dV = gradient(V_u); dQ = gradient(Q_u);
                dV(dV==0) = NaN;
                dQdV = dQ ./ dV;
                dQdV_filt = movmean(dQdV, 21);
                
                mask = V_u >= 3.6 & V_u <= 4.0;
                plot(V_u(mask), dQdV_filt(mask), 'Color', cyc_colors(c,:), ...
                    'LineWidth', 1.5, 'DisplayName', cyc);
                
                % Peak marker
                [pk_V, pk_H] = local_get_peak(data, 3.40, 4.00);
                if ~isnan(pk_V)
                    plot(pk_V, pk_H, 'o', 'MarkerEdgeColor', cyc_colors(c,:), ...
                        'MarkerFaceColor', cyc_colors(c,:), 'MarkerSize', 5, ...
                        'HandleVisibility', 'off');
                end
            end
        end
        
        title(sprintf('%s', ch), 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Voltage (V)'); ylabel('dQ/dV (Ah/V)');
        grid on; xlim([3.6 4.0]);
        if i == 1, legend('Location', 'northeast', 'FontSize', 7); end
    end
    
    sgtitle('OCV dQ/dV Peak Shift & Shrinkage', 'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig3, fullfile(saveDir, 'Features_OCV_Peak_Evolution.fig'));
end

%% ===== Individual Label Plots (if specific Label sub-items checked) =====
label_map = {'Available Capacity (Ah)', 1; 'LLI', 2; 'LAMp', 3};

for lm = 1:size(label_map, 1)
    node_name = label_map{lm, 1};
    lbl_idx = label_map{lm, 2};
    
    if any(strcmpi(checkedNodes, node_name))
        fprintf('  [%s] Generating individual trend plot...\n', node_name);
        fig_lbl = figure('Position', [100, 100, 800, 500], 'Name', sprintf('Label: %s', node_name));
        hold on;
        
        for i = 1:length(channels)
            ch = channels{i};
            if ismember('CrateNum', FeatureTable.Properties.VariableNames)
                idx = strcmp(FeatureTable.CellID, ch) & FeatureTable.CrateNum == 0.5;
            else
                idx = strcmp(FeatureTable.CellID, ch);
            end
            subT = FeatureTable(idx, :);
            if isempty(subT), continue; end
            [~, si] = sort(subT.Cycle);
            
            mk = marker_list{mod(i-1, length(marker_list))+1};
            plot(subT.Cycle(si), subT.Y_Labels(si, lbl_idx), ['-', mk], ...
                'Color', ch_colors(i,:), 'LineWidth', 1.5, 'MarkerSize', 6, ...
                'MarkerFaceColor', ch_colors(i,:), 'DisplayName', ch);
        end
        
        title(sprintf('%s Trend', label_names{lbl_idx}), 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('Cycle', 'FontWeight', 'bold', 'FontSize', 12);
        ylabel(label_units{lbl_idx}, 'FontWeight', 'bold', 'FontSize', 12);
        legend('Location', 'best', 'FontSize', 9); grid on; hold off;
        
        saveas(fig_lbl, fullfile(saveDir, sprintf('Features_%s_Trend.fig', label_names{lbl_idx})));
    end
end

fprintf('--- Features & Labels Visualization Complete ---\n');
end

%% =========================================================================
% Local Helper
% =========================================================================
function [pk_V, pk_H] = local_get_peak(data_struct, minV, maxV)
    V = data_struct.V_grid;
    Q = data_struct.Q;
    if numel(V) < 10, pk_V = NaN; pk_H = NaN; return; end
    
    [V_u, uid] = unique(V);
    Q_u = Q(uid);
    dV = gradient(V_u); dQ = gradient(Q_u);
    dV(dV==0) = NaN;
    dQdV = dQ ./ dV;
    dQdV(isinf(dQdV)) = NaN;
    dQdV = fillmissing(dQdV, 'linear');
    dQdV_filt = movmean(dQdV, 21);
    
    mask = V_u >= minV & V_u <= maxV;
    V_roi = V_u(mask);
    dH_roi = dQdV_filt(mask);
    if isempty(V_roi), pk_V = NaN; pk_H = NaN; return; end
    
    [pks, locs] = findpeaks(dH_roi);
    if isempty(pks)
        [pk_H, idx] = max(dH_roi);
        pk_V = V_roi(idx);
    else
        [pk_H, idx] = max(pks);
        pk_V = V_roi(locs(idx));
    end
end
