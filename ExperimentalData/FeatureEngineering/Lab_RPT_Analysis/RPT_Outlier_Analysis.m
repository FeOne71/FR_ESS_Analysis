%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Outlier Analysis: Non-monotonic Behavior Verification (dQ/dV Version)
% Comparison: Ch14 (Normal) vs Ch16 (Outlier) - All 5 Segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% ========================================================================
% Section 1: Data Loading & Configuration
% ========================================================================
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_rpt_vq_mat = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\Pipeline_Visualizations\Outlier';

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

fprintf('Loading Data...\n');
load(path_rpt_vq_mat, 'RPT_VQ_grid');
fprintf('Data Loaded.\n');

% Constants (Sync with Master Ruler logic)
win_chg_min = 3.7;
win_chg_max = 3.95;
num_segments = 5;
target_ch = {'Ch14', 'Ch16'}; 
crate = 'c01_charge'; 

%% ========================================================================
% Section 2: Master Ruler Reconstruction (Global)
% ========================================================================
channels = {'Ch09','Ch10','Ch11','Ch12','Ch13','Ch14','Ch15','Ch16'};

% Collect V-t data for Master Ruler (using cyc0 average)
standard_V_grid_chg = RPT_VQ_grid.cyc0.Ch09.c05_charge.V_grid;
mask_std_chg = standard_V_grid_chg >= win_chg_min & standard_V_grid_chg <= win_chg_max;
V_standard_chg = standard_V_grid_chg(mask_std_chg);

all_T_grids_chg = zeros(length(channels), length(V_standard_chg));

for i = 1:length(channels)
    ch = channels{i};
    data_chg = RPT_VQ_grid.cyc0.(ch).c05_charge;
    V_grid_chg = data_chg.V_grid;
    t_grid_chg = data_chg.t;
    mask_chg = V_grid_chg >= win_chg_min & V_grid_chg <= win_chg_max;
    t_chg_interp = interp1(V_grid_chg(mask_chg), t_grid_chg(mask_chg), V_standard_chg, 'linear', 'extrap');
    all_T_grids_chg(i, :) = t_chg_interp';
end

avg_T_grid_chg = mean(all_T_grids_chg, 1, 'omitnan');
T_start_chg = min(avg_T_grid_chg);
T_end_chg = max(avg_T_grid_chg);
target_Ts_chg = linspace(T_start_chg, T_end_chg, num_segments + 1);
Global_V_bounds_chg = interp1(avg_T_grid_chg, V_standard_chg, target_Ts_chg, 'linear');

MasterRulers = struct();
for i = 1:length(channels)
    MasterRulers.(channels{i}).V_bounds_chg = Global_V_bounds_chg;
end

%% ========================================================================
% Section 3: Visualization for All 5 Segments (dQ/dV Version)
% ========================================================================
cycs = fieldnames(RPT_VQ_grid);
[~, sort_idx] = sort(cellfun(@(x) str2double(regexprep(x, '\D', '')), cycs));
cycs = cycs(sort_idx);
colors = lines(length(cycs)); % jet 대신 구분이 뚜렷한 lines colormap 사용

for seg_idx = 1:num_segments
    fprintf('Generating Plot for Segment %d...\n', seg_idx);
    
    fig = figure('Position', [100, 100, 1400, 600], 'Name', sprintf('Peak Shift Proof - Seg %d', seg_idx));
    
    for i = 1:2
        subplot(1, 2, i); hold on;
        ch = target_ch{i};
        
        % Segment 전압 경계
        V_start = MasterRulers.(ch).V_bounds_chg(seg_idx);
        V_end = MasterRulers.(ch).V_bounds_chg(seg_idx+1);
        
        y_max = 0;
        
        for j = 1:length(cycs)
            cyc_key = cycs{j};
            if isfield(RPT_VQ_grid, cyc_key) && isfield(RPT_VQ_grid.(cyc_key), ch) && ...
               isfield(RPT_VQ_grid.(cyc_key).(ch), crate)
                
                data = RPT_VQ_grid.(cyc_key).(ch).(crate);
                [V_u, uid] = unique(data.V_grid);
                Q_u = data.Q(uid);
                
                % dQ/dV 계산 (Smoothing 적용)
                dV = gradient(V_u);
                dQ = gradient(Q_u);
                dV(dV==0) = NaN;
                dQdV = dQ ./ dV;
                dQdV_filt = movmean(dQdV, 51);
                
                mask_win = V_u >= 3.7 & V_u <= 3.95;
                plot(V_u(mask_win), dQdV_filt(mask_win), 'Color', colors(j,:), 'LineWidth', 1.5, 'DisplayName', cyc_key);
                
                % Y축 범위 결정을 위한 최대값 추적
                curr_max = max(dQdV_filt(mask_win));
                if curr_max > y_max, y_max = curr_max; end
            end
        end
        
        % Y축 상한 설정 (조금 여유 있게)
        ylim_upper = y_max * 1.2;
        if isnan(ylim_upper) || ylim_upper <= 0, ylim_upper = 25; end
        ylim([0 ylim_upper]);
        
        % 고정된 세그먼트 경계 표시 (점선 및 투명 박스)
        patch([V_start V_end V_end V_start], [0 0 ylim_upper ylim_upper], 'k', ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        line([V_start V_start], [0 ylim_upper], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2, 'HandleVisibility', 'off');
        line([V_end V_end], [0 ylim_upper], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2, 'HandleVisibility', 'off');
        
        title(sprintf('%s: dQ/dV Peak Shift (Seg%d)', ch, seg_idx), 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Voltage (V)', 'FontWeight', 'bold');
        ylabel('dQ/dV (Ah/V)', 'FontWeight', 'bold');
        grid on;
        
        % 레전드 표시 (모든 subplot)
        legend('Location', 'northeast', 'FontSize', 8, 'Interpreter', 'none');
    end
    
    sgtitle(sprintf('Proof of Non-monotonicity (C-rate: %s): Peak Moving across Seg%d Boundary', crate, seg_idx), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % 파일 저장    
    saveas(fig, fullfile(saveDir, sprintf('PeakShift_Proof_Seg%d.fig', seg_idx)));
    close(fig);
end

fprintf('\nAll segment (1~5) dQ/dV visualizations complete.\n');

%% ========================================================================
% Section 4: V-Q 시각화 (용량 너비 변화 확인 - 전 구간)
% ========================================================================
fprintf('\n=== Section 4: V-Q Visualization (Voltage Segment Bounds) ===\n');

for seg_idx = 1:num_segments
    fprintf('Generating V-Q Plot with Boundaries for Segment %d...\n', seg_idx);
    
    fig = figure('Position', [100, 100, 1400, 600], 'Name', sprintf('V-Q Boundary Check - Seg %d', seg_idx));
    
    for i = 1:2
        subplot(1, 2, i); hold on;
        ch = target_ch{i};
        
        % Segment 전압 경계
        V_start = MasterRulers.(ch).V_bounds_chg(seg_idx);
        V_end = MasterRulers.(ch).V_bounds_chg(seg_idx+1);
        
        for j = 1:length(cycs)
            cyc_key = cycs{j};
            if isfield(RPT_VQ_grid, cyc_key) && isfield(RPT_VQ_grid.(cyc_key), ch) && ...
               isfield(RPT_VQ_grid.(cyc_key).(ch), crate)
                
                data = RPT_VQ_grid.(cyc_key).(ch).(crate);
                
                % 분석 윈도우(3.7~3.95V) 전 구간 출력
                mask_win = data.V_grid >= 3.7 & data.V_grid <= 3.95;
                
                % 가독성을 위해 시작 용량을 0으로 맞춤 (Relative Q)
                Q_base = min(data.Q(mask_win));
                plot(data.Q(mask_win) - Q_base, data.V_grid(mask_win), ...
                    'Color', colors(j,:), 'LineWidth', 1.5, 'DisplayName', cyc_key);
            end
        end
        
        % 세그먼트 구간 표시 (가로 음영 구역)
        xl = xlim;
        patch([xl(1) xl(2) xl(2) xl(1)], [V_start V_start V_end V_end], 'k', ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        yline(V_start, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        yline(V_end, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        
        title(sprintf('%s: Full V-Q with Seg%d Highlight', ch, seg_idx), 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Relative Capacity (Ah)', 'FontWeight', 'bold');
        ylabel('Voltage (V)', 'FontWeight', 'bold');
        grid on;
        
        % 레전드 표시 (모두 표시)
        legend('Location', 'northwest', 'FontSize', 8, 'Interpreter', 'none');
    end
    
    sgtitle(sprintf('V-Q Segment Analysis (C-rate: %s): How Peak Shift affects dQ in Seg%d', crate, seg_idx), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % 파일 저장 (.fig만 저장)
    saveas(fig, fullfile(saveDir, sprintf('VQ_Boundary_Proof_Seg%d.fig', seg_idx)));
    close(fig);
end

fprintf('\nAll visual analysis complete.\nSaved at: %s\n', saveDir);
