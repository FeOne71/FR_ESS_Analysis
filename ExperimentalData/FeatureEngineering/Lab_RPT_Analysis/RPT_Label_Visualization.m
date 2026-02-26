%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Label Trend Visualization
% Plots SOH, LLI, and LAM trends from Feature_Matrix_Final.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% 1. Load Data
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis';
featureMatPath = fullfile(baseDir, 'Dataset', 'Feature_Matrix_Final.mat');
saveDir = fullfile(baseDir, 'Pipeline_Visualizations', 'Label_Trends');

if ~exist(saveDir, 'dir'), mkdir(saveDir); end

fprintf('Loading Feature Matrix...\n');
load(featureMatPath, 'FeatureTable');
fprintf('Data Loaded.\n');

channels = unique(FeatureTable.CellID);
labels = {'SOH_Capacity', 'LLI', 'LAM'};
label_titles = {'Available Capacity)', 'LLI (Loss of Lithium Inventory)', 'LAM (Loss of Active Material)'};
label_units = {'Capacity (Ah)', 'LLI (%)', 'LAM (%)'};

% 색상 설정 (8개 채널)
colors = lines(length(channels));

%% 2. Visualization
fig = figure('Position', [100, 100, 1500, 500], 'Name', 'Battery Aging Label Trends');

for l = 1:3
    subplot(1, 3, l); hold on;
    label_name = labels{l};
    
    for i = 1:length(channels)
        ch = channels{i};
        
        % 특정 채널 데이터 필터링 
        % (라벨은 C-rate에 무관하므로 첫 번째 C-rate 데이터만 사용)
        unique_crates = unique(FeatureTable.CrateLabel);
        idx = strcmp(FeatureTable.CellID, ch) & strcmp(FeatureTable.CrateLabel, unique_crates{1});
        
        subTable = FeatureTable(idx, :);
        
        % 사이클 순으로 정렬
        [~, sort_idx] = sort(subTable.Cycle);
        
        % Channel-specific marker styles in addition to color
        marker_list = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
        mk_style = marker_list{mod(i-1, length(marker_list))+1};
        
        plot(subTable.Cycle(sort_idx), subTable.Y_Labels(sort_idx, l), ...
            ['-', mk_style], 'Color', colors(i,:), 'LineWidth', 1.5, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', ch);
    end
    
    title(label_titles{l}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Cycle', 'FontWeight', 'bold');
    ylabel(label_units{l}, 'FontWeight', 'bold');
    grid on;
    
    if l == 1
        legend('Location', 'best', 'FontSize', 8);
    end
end

sgtitle('Aging Label Trends across All Channels (SOH / LLI / LAM)', 'FontSize', 14, 'FontWeight', 'bold');

% 파일 저장
saveas(fig, fullfile(saveDir, 'Label_Trends_All_Channels.fig'));

%% 3. OCV dQ/dV Peak Visualization (ICA Analysis)
fprintf('\nGenerating OCV dQ/dV Peak Visualizations...\n');

% RPT_VQ_grid.mat 로드 (OCV 데이터 필요)
path_rpt_vq_mat = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat';
load(path_rpt_vq_mat, 'RPT_VQ_grid');

fig2 = figure('Position', [50, 50, 1600, 900], 'Name', 'OCV dQ/dV Peak Evolution');
cyc_fields = fieldnames(RPT_VQ_grid);
% Use Sequential Colormap (jet/parula) for Aging Trends (Cycle 0 -> 1000)
cyc_colors = jet(length(cyc_fields));

for i = 1:length(channels)
    ch = channels{i};
    subplot(2, 4, i); hold on;
    
    for c = 1:length(cyc_fields)
        cyc = cyc_fields{c};
        
        % OCV를 대표하는 0.05C 데이터 사용 (또는 ocv 필드)
        if isfield(RPT_VQ_grid, cyc) && isfield(RPT_VQ_grid.(cyc), ch) && ...
           isfield(RPT_VQ_grid.(cyc).(ch), 'OCV_charge')
            
            data = RPT_VQ_grid.(cyc).(ch).OCV_charge;
            [V_u, uid] = unique(data.V_grid);
            Q_u = data.Q(uid);
            
            % dQ/dV 계산 및 스무딩 (Window=21)
            dV = gradient(V_u);
            dQ = gradient(Q_u);
            dV(dV==0) = NaN;
            dQdV = dQ ./ dV;
            dQdV_filt = movmean(dQdV, 21);
            
            % 관심 구간 (3.6~4.0V)
            mask = V_u >= 3.6 & V_u <= 4.0;
            plot(V_u(mask), dQdV_filt(mask), 'Color', cyc_colors(c,:), 'LineWidth', 1.5, 'DisplayName', cyc);
            
            % 피크 검출 및 표시 (LLI/LAM 시각적 근거)
            [pk_V, pk_H] = get_main_peak(data, 3.40, 4.00);
            if ~isnan(pk_V)
                plot(pk_V, pk_H, 'o', 'MarkerEdgeColor', cyc_colors(c,:), 'MarkerFaceColor', cyc_colors(c,:), ...
                    'MarkerSize', 5, 'HandleVisibility', 'off');
            end
        end
    end
    
    title(sprintf('Channel: %s', ch), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Voltage (V)');
    ylabel('dQ/dV (Ah/V)');
    grid on;
    xlim([3.6 4.0]); ylim([40 240]);
    
    if i == 1
        legend('Location', 'northeast', 'FontSize', 7);
    end
end

sgtitle('OCV dQ/dV Peak Shift & Shrinkage', 'FontSize', 14, 'FontWeight', 'bold');

% 파일 저장
saveas(fig2, fullfile(saveDir, 'OCV_ICA_Peak_Evolution.fig'));

%% Section 4: Data Summary Tables (Vpeak, Hpeak, SOH)
fprintf('\n==========================================================\n');
fprintf('Section 4: Data Summary Tables (0 ~ 1000 cyc)\n');
fprintf('==========================================================\n');

% Initialize data storage for tables
cyc_vals = [0, 200, 400, 600, 800, 1000];
Vpeak_mat = nan(length(channels), length(cyc_vals));
Hpeak_mat = nan(length(channels), length(cyc_vals));
SOH_mat = nan(length(channels), length(cyc_vals));
LLI_mat = nan(length(channels), length(cyc_vals));
LAM_mat = nan(length(channels), length(cyc_vals));
% Initialize detailed metrics
dV_mat = nan(length(channels), length(cyc_vals)); % V_aged - V_fresh
Hfresh_vec = nan(length(channels), 1);
Qrated_vec = nan(length(channels), 1);

for i = 1:length(channels)
    ch = channels{i};
    ch_idx = strcmp(FeatureTable.CellID, ch) & strcmp(FeatureTable.CrateLabel, 'c01');
    ch_data_table = FeatureTable(ch_idx, :);
    
    % Get Q_rated & H_fresh from Cycle 0 data within FeatureTable or Calculated Matrix
    % Since we fill Vpeak_mat, Hpeak_mat, SOH_mat inside the loop, we can extract 0cyc values after the loop
    
    for j = 1:length(cyc_vals)
        c_val = cyc_vals(j);
        c_field = sprintf('cyc%d', c_val);
        
        if isfield(RPT_VQ_grid, c_field) && isfield(RPT_VQ_grid.(c_field), ch) && ...
           isfield(RPT_VQ_grid.(c_field).(ch), 'OCV_charge')
            
            ocv_struct = RPT_VQ_grid.(c_field).(ch).OCV_charge;
            [pV, pH] = get_main_peak(ocv_struct, 3.40, 4.00);
            Vpeak_mat(i, j) = pV;
            Hpeak_mat(i, j) = pH;
            % dV calculation deferred to after loop when V_fresh becomes known
        end
        
        r_idx = ch_data_table.Cycle == c_val;
        if any(r_idx)
            SOH_mat(i, j) = ch_data_table.Y_Labels(r_idx, 1);
            LLI_mat(i, j) = ch_data_table.Y_Labels(r_idx, 2);
            LAM_mat(i, j) = ch_data_table.Y_Labels(r_idx, 3);
        end
    end
    
    % Extract Fresh Parameters (Cycle 0)
    idx_0cyc = find(cyc_vals == 0, 1);
    if ~isempty(idx_0cyc)
        if ~isnan(SOH_mat(i, idx_0cyc))
             Qrated_vec(i) = SOH_mat(i, idx_0cyc);
        else
             Qrated_vec(i) = 55.6; % Fallback
        end
        
        if ~isnan(Hpeak_mat(i, idx_0cyc))
             Hfresh_vec(i) = Hpeak_mat(i, idx_0cyc);
             pV_fresh = Vpeak_mat(i, idx_0cyc);
        else
             Hfresh_vec(i) = NaN;
             pV_fresh = NaN;
        end
        
        % Calculate dV based on V_fresh extracted from loop
        if ~isnan(pV_fresh)
            dV_mat(i, :) = Vpeak_mat(i, :) - pV_fresh;
        end
    end
end

% Print Tables
display_table('Vpeak (V)', channels, cyc_vals, Vpeak_mat);
display_table('Hpeak (Ah/V)', channels, cyc_vals, Hpeak_mat);
display_table('SOH (Available Capacity)', channels, cyc_vals, SOH_mat);
display_table('LLI (%)', channels, cyc_vals, LLI_mat);
display_table('LAM (%)', channels, cyc_vals, LAM_mat);

% Detailed LLI Components
fprintf('\n[ LLI Detailed Components ]\n');
fprintf('%-8s | %10s | %12s\n', 'Channel', 'Q_rated(Ah)', 'H_fresh(Ah/V)');
fprintf('%s\n', repmat('-', 1, 36));
for i=1:length(channels)
    fprintf('%-8s | %10.4f | %12.4f\n', channels{i}, Qrated_vec(i), Hfresh_vec(i));
end
display_table('dV_peak (V_aged - V_fresh)', channels, cyc_vals, dV_mat);

fprintf('\nAll visualizations and data summaries complete.\n');
fprintf('Label Trends: %s\n', fullfile(saveDir, 'Label_Trends_All_Channels.fig'));
fprintf('OCV Peak Evolution: %s\n', fullfile(saveDir, 'OCV_ICA_Peak_Evolution.fig'));

%% Helper Functions
function display_table(title_str, channels, cyc_vals, data)
    fprintf('\n[ %s ]\n', title_str);
    fprintf('%-8s', 'Channel');
    for c = 1:length(cyc_vals)
        fprintf('| %7dcyc ', cyc_vals(c));
    end
    fprintf('\n%s\n', repmat('-', 1, 8 + 12*length(cyc_vals)));
    for i = 1:length(channels)
        fprintf('%-8s', channels{i});
        for j = 1:length(cyc_vals)
            if isnan(data(i,j)), fprintf('| %10s ', 'N/A');
            else, fprintf('| %10.4f ', data(i,j)); end
        end
        fprintf('\n');
    end
end

function [pk_V, pk_H] = get_main_peak(data_struct, minV, maxV)
    % Extracts main peak from OCV dQ/dV curve for LLI/LAM labeling
    V = data_struct.V_grid;
    Q = data_struct.Q;
    [V_u, uid] = unique(V);
    Q_u = Q(uid);
    dV = gradient(V_u);
    dQ = gradient(Q_u);
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
