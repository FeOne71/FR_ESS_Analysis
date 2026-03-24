%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Feature Extractor
% - Extraction of Physical Features (Table 2)
% - Features: dQ_seg, Voltage Efficiency (eta_U), 
%   Energy Efficiency (eta_Wh), Effective C-rate (C_eff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
warning off;

%% ========================================================================
% Section 1: Configuration & Data Loading
% ========================================================================
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_rpt_vq_mat = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');
path_static_mat = fullfile(baseDir, 'Capacity_Trend_Figures', 'Capacity_Data_Static.mat');
path_master_ruler = fullfile(baseDir, 'FeatureEngineering\Lab_RPT_Analysis\Dataset', 'MasterRulers.mat');
saveDir = fullfile(baseDir, 'FeatureEngineering\Lab_RPT_Analysis\Dataset');

if ~exist(saveDir, 'dir'), mkdir(saveDir); end

fprintf('Loading Data...\n');
load(path_rpt_vq_mat, 'RPT_VQ_grid');
load(path_static_mat, 'allChannelsCapacity');
load(path_master_ruler, 'MasterRulers');

channels = fieldnames(allChannelsCapacity);
target_crates = {'c01', 'c05', 'c1', 'c2', 'c3'};
target_crates_val = [0.1, 0.5, 1, 2, 3];
num_segments = 5;

% Constants
V_min_sop = 3.0; % Per cell
V_max_sop = 4.2; 
C_F = 1.0; % Correction factor (default)

%% ========================================================================
% Section 2: Feature Extraction Loop
% ========================================================================
fprintf('\nExtracting Physics-based Features...\n');
data_CellID = {};
data_Cycle = [];
data_CrateLabel = {};
data_CrateNum = [];
data_X = []; % [dQ_seg(10), eta_U, eta_Wh, C_eff] -> Total 13 features

cnt = 0;
for i = 1:length(channels)
    ch = channels{i};
    if ~isfield(MasterRulers, ch), continue; end
    
    Q_nom = 63.5; % Default nominal Ah for 64Ah cell
    
    % Get Ruler
    VR_chg = MasterRulers.(ch).V_bounds_chg;
    VR_dch = MasterRulers.(ch).V_bounds_dch;
    
    cyc_fields = fieldnames(RPT_VQ_grid);
    for c = 1:length(cyc_fields)
        cyc_key = cyc_fields{c};
        cyc_num = sscanf(cyc_key, 'cyc%d');
        if ~isfield(RPT_VQ_grid.(cyc_key), ch), continue; end
        
        ch_data = RPT_VQ_grid.(cyc_key).(ch);
        
        for r = 1:length(target_crates)
            crate_label = target_crates{r};
            crate_val = target_crates_val(r);
            f_chg = [crate_label '_charge'];
            f_dch = [crate_label '_discharge'];
            
            if ~isfield(ch_data, f_chg) || ~isfield(ch_data, f_dch)
                continue;
            end
            
            % --- [1] dQ_seg (Segmented Capacity) ---
            % Formula: dQ_seg,k = integral(I dt) within Master Ruler Voltage Boundaries
            s_chg = ch_data.(f_chg);
            s_dch = ch_data.(f_dch);
            
            dQ_chg = nan(1, num_segments);
            dQ_dch = nan(1, num_segments);
            
            % Interpolate Q at Voltage Boundaries
            [V_u_c, uid_c] = unique(s_chg.V_raw);
            Q_u_c = s_chg.Q_raw(uid_c);
            [V_u_d, uid_d] = unique(s_dch.V_raw);
            Q_u_d = s_dch.Q_raw(uid_d);
            
            if length(V_u_c) > 1 && length(V_u_d) > 1
                QR_chg = interp1(V_u_c, Q_u_c, VR_chg, 'linear');
                dQ_chg = abs(diff(QR_chg));
                
                QR_dch = interp1(V_u_d, Q_u_d, VR_dch, 'linear');
                dQ_dch = abs(diff(QR_dch));
            end
            
            % --- [2] η_U (Voltage Efficiency) ---
            % η_U = (∫V_cha·I dt / ∫I dt) / (∫V_dis·I dt / ∫I dt) × C_F
            % This is Avg_V_cha / Avg_V_dis
            V_cha_avg = trapz(s_chg.t_raw, s_chg.V_raw .* abs(s_chg.I_raw)) / trapz(s_chg.t_raw, abs(s_chg.I_raw));
            V_dis_avg = trapz(s_dch.t_raw, s_dch.V_raw .* abs(s_dch.I_raw)) / trapz(s_dch.t_raw, abs(s_dch.I_raw));
            eta_U = (V_cha_avg / V_dis_avg) * C_F;
            
            % --- [3] η_Wh (Energy Efficiency) ---
            % η_Wh = (∫V_dis·I_dis dt) / (∫V_cha·I_cha dt) 
            E_dis = trapz(s_dch.t_raw, s_dch.V_raw .* abs(s_dch.I_raw));
            E_cha = trapz(s_chg.t_raw, s_chg.V_raw .* abs(s_chg.I_raw));
            eta_Wh = (E_dis / E_cha) * C_F;
            
            % --- [4] C_eff (Effective C-rate) ---
            % C_eff = I_test / Q_nominal
            I_avg = (mean(abs(s_chg.I_raw)) + mean(abs(s_dch.I_raw))) / 2;
            C_eff = I_avg / Q_nom;
            
            % Collect row
            cnt = cnt + 1;
            data_CellID{cnt,1} = ch;
            data_Cycle(cnt,1) = cyc_num;
            data_CrateLabel{cnt,1} = crate_label;
            data_CrateNum(cnt,1) = crate_val;
            data_X(cnt,:) = [dQ_chg, dQ_dch, eta_U, eta_Wh, C_eff];
        end
    end
end

% Finalize Table
FeatureTable_Physics = table(data_CellID, data_Cycle, data_CrateLabel, data_CrateNum, ...
    data_X(:, 1:5), data_X(:, 6:10), data_X(:, 11), data_X(:, 12), data_X(:, 13), ...
    'VariableNames', {'CellID', 'Cycle', 'CrateLabel', 'CrateNum', 'dQ_chg_seg', 'dQ_dch_seg', 'eta_U', 'eta_Wh', 'C_eff'});

save(fullfile(saveDir, 'Feature_Matrix_Physics.mat'), 'FeatureTable_Physics');
fprintf('Physics-based Features Extracted and Saved (%d rows).\n', cnt);

%% ========================================================================
% Section 3: Visualization (Efficiency Trends)
% ========================================================================
fprintf('\nGenerating Efficiency Visualization...\n');
fig = figure('Position', [100, 100, 1200, 500]);

% Plot eta_U vs Cycle
subplot(1,2,1); hold on;
crates = unique(data_CrateNum);
colors = lines(length(crates));
for r = 1:length(crates)
    idx = (data_CrateNum == crates(r)) & strcmp(data_CellID, 'Ch09');
    plot(data_Cycle(idx), data_X(idx, 11), '-o', 'Color', colors(r,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%.1fC', crates(r)));
end
xlabel('Cycle'); ylabel('\eta_U (V_{cha,avg} / V_{dis,avg})');
title('Voltage Efficiency Transition (Ch09)'); grid on; legend('Location', 'best');

% Plot eta_Wh vs Cycle
subplot(1,2,2); hold on;
for r = 1:length(crates)
    idx = (data_CrateNum == crates(r)) & strcmp(data_CellID, 'Ch09');
    plot(data_Cycle(idx), data_X(idx, 12), '-s', 'Color', colors(r,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%.1fC', crates(r)));
end
xlabel('Cycle'); ylabel('\eta_{Wh} (Energy Efficiency)');
title('Energy Efficiency Transition (Ch09)'); grid on;
saveas(fig, fullfile(saveDir, 'Efficiency_Trend.fig'));
fprintf('Visualization Saved.\n');
