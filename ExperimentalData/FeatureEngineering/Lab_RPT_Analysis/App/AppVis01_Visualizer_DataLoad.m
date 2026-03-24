function App_Visualizer_DataLoad(App_VQ_grid, saveDir)
% APP_VISUALIZER_DATALOAD Visualizes loaded V-Q data from App_DataLoader
%
% Reference: RPT_Pipeline_Visualization.m Phase 0
%
% Figures:
%   1. V-Q Overlay (OCV Charge/Discharge, all channels × all cycles)
%   2. Voltage Window (analysis window shading on representative channel)
%   3. SOH Trend (Static Capacity vs Cycle)

fprintf('--- App Visualizer: Data Load ---\n');

if nargin < 2 || isempty(saveDir)
    saveDir = fullfile(fileparts(mfilename('fullpath')), 'Results');
end
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

cyc_fields = fieldnames(App_VQ_grid);
% Collect all channels across all cycles
all_channels = {};
for c = 1:length(cyc_fields)
    ch_list = fieldnames(App_VQ_grid.(cyc_fields{c}));
    all_channels = union(all_channels, ch_list);
end
channels = sort(all_channels);

cyc_colors = lines(length(cyc_fields));
ch_colors = lines(length(channels));

%% Figure 1: V-Q Overlay (OCV Charge / Discharge)
fig1 = figure('Position', [50, 50, 1600, 700], 'Name', 'Data Load: V-Q Overlay');

% --- Charge ---
subplot(1,2,1); hold on;
for c = 1:length(cyc_fields)
    cyc = cyc_fields{c};
    ch_list = fieldnames(App_VQ_grid.(cyc));
    for i = 1:length(ch_list)
        ch = ch_list{i};
        if isfield(App_VQ_grid.(cyc).(ch), 'OCV_charge') && ...
           isfield(App_VQ_grid.(cyc).(ch).OCV_charge, 'V_grid')
            data = App_VQ_grid.(cyc).(ch).OCV_charge;
            if i == 1
                plot(data.Q, data.V_grid, 'Color', cyc_colors(c,:), ...
                    'LineWidth', 1, 'DisplayName', cyc);
            else
                plot(data.Q, data.V_grid, 'Color', cyc_colors(c,:), ...
                    'LineWidth', 1, 'HandleVisibility', 'off');
            end
        end
    end
end
xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
title('OCV Charge V-Q (All Channels)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9); grid on; hold off;

% --- Discharge ---
subplot(1,2,2); hold on;
for c = 1:length(cyc_fields)
    cyc = cyc_fields{c};
    ch_list = fieldnames(App_VQ_grid.(cyc));
    for i = 1:length(ch_list)
        ch = ch_list{i};
        if isfield(App_VQ_grid.(cyc).(ch), 'OCV_discharge') && ...
           isfield(App_VQ_grid.(cyc).(ch).OCV_discharge, 'V_grid')
            data = App_VQ_grid.(cyc).(ch).OCV_discharge;
            if i == 1
                plot(data.Q, data.V_grid, 'Color', cyc_colors(c,:), ...
                    'LineWidth', 1, 'DisplayName', cyc);
            else
                plot(data.Q, data.V_grid, 'Color', cyc_colors(c,:), ...
                    'LineWidth', 1, 'HandleVisibility', 'off');
            end
        end
    end
end
xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
title('OCV Discharge V-Q (All Channels)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9); grid on; hold off;

sgtitle(sprintf('V-Q Overlay (%d Channels x %d Cycles)', length(channels), length(cyc_fields)), ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig1, fullfile(saveDir, 'DataLoad_VQ_Overlay.fig'));

%% Figure 2: Voltage Window
fig2 = figure('Position', [50, 50, 1600, 600], 'Name', 'Data Load: Voltage Windows');

win_chg_min = 3.7; win_chg_max = 3.95;
win_dch_min = 3.75; win_dch_max = 3.88;

% Use first available channel
ref_ch = channels{1};

% Charge Window
subplot(1,2,1); hold on;
for c = 1:length(cyc_fields)
    cyc = cyc_fields{c};
    if isfield(App_VQ_grid.(cyc), ref_ch) && isfield(App_VQ_grid.(cyc).(ref_ch), 'OCV_charge') && ...
       isfield(App_VQ_grid.(cyc).(ref_ch).OCV_charge, 'V_grid')
        data = App_VQ_grid.(cyc).(ref_ch).OCV_charge;
        plot(data.Q, data.V_grid, 'Color', cyc_colors(c,:), 'LineWidth', 1.5, 'DisplayName', cyc);
    end
end
xl = xlim;
fill([xl(1) xl(2) xl(2) xl(1)], [win_chg_min win_chg_min win_chg_max win_chg_max], ...
    'g', 'FaceAlpha', 0.1, 'EdgeColor', 'g', 'LineWidth', 2, 'LineStyle', '--', ...
    'DisplayName', sprintf('Window [%.2f-%.2fV]', win_chg_min, win_chg_max));
xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
title('Charge Analysis Window', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9); grid on; hold off;

% Discharge Window
subplot(1,2,2); hold on;
for c = 1:length(cyc_fields)
    cyc = cyc_fields{c};
    if isfield(App_VQ_grid.(cyc), ref_ch) && isfield(App_VQ_grid.(cyc).(ref_ch), 'OCV_discharge') && ...
       isfield(App_VQ_grid.(cyc).(ref_ch).OCV_discharge, 'V_grid')
        data = App_VQ_grid.(cyc).(ref_ch).OCV_discharge;
        plot(data.Q, data.V_grid, 'Color', cyc_colors(c,:), 'LineWidth', 1.5, 'DisplayName', cyc);
    end
end
xl = xlim;
fill([xl(1) xl(2) xl(2) xl(1)], [win_dch_min win_dch_min win_dch_max win_dch_max], ...
    'r', 'FaceAlpha', 0.1, 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--', ...
    'DisplayName', sprintf('Window [%.2f-%.2fV]', win_dch_min, win_dch_max));
xlabel('Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 12);
title('Discharge Analysis Window', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9); grid on; hold off;

sgtitle(sprintf('Analysis Voltage Windows (%s)', ref_ch), 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig2, fullfile(saveDir, 'DataLoad_Voltage_Windows.fig'));

%% Figure 3: SOH Trend (Static Capacity)
fig3 = figure('Position', [50, 50, 1400, 600], 'Name', 'Data Load: SOH Trend');

subplot(1,2,1); hold on;
for i = 1:length(channels)
    ch = channels{i};
    cyc_nums = [];
    cap_vals = [];
    for c = 1:length(cyc_fields)
        cyc = cyc_fields{c};
        cyc_num = sscanf(cyc, 'cyc%d');
        if isfield(App_VQ_grid.(cyc), ch) && isfield(App_VQ_grid.(cyc).(ch), 'Static') && ...
           isfield(App_VQ_grid.(cyc).(ch).Static, 'Q') && ~isempty(App_VQ_grid.(cyc).(ch).Static.Q)
            cyc_nums(end+1) = cyc_num;
            cap_vals(end+1) = max(abs(App_VQ_grid.(cyc).(ch).Static.Q));
        end
    end
    if ~isempty(cyc_nums)
        [cyc_nums, si] = sort(cyc_nums);
        cap_vals = cap_vals(si);
        plot(cyc_nums, cap_vals, 'o-', 'Color', ch_colors(i,:), ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', ch);
    end
end
xlabel('Cycle', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Maximum Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 12);
title('SOH: Static Capacity', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9); grid on; hold off;

subplot(1,2,2); hold on;
for i = 1:length(channels)
    ch = channels{i};
    cyc_nums = [];
    cap_vals = [];
    for c = 1:length(cyc_fields)
        cyc = cyc_fields{c};
        cyc_num = sscanf(cyc, 'cyc%d');
        if isfield(App_VQ_grid.(cyc), ch) && isfield(App_VQ_grid.(cyc).(ch), 'Static') && ...
           isfield(App_VQ_grid.(cyc).(ch).Static, 'Q') && ~isempty(App_VQ_grid.(cyc).(ch).Static.Q)
            cyc_nums(end+1) = cyc_num;
            cap_vals(end+1) = max(abs(App_VQ_grid.(cyc).(ch).Static.Q));
        end
    end
    if length(cyc_nums) > 1
        [cyc_nums, si] = sort(cyc_nums);
        cap_vals = cap_vals(si);
        cap_pct = (cap_vals / cap_vals(1)) * 100;
        plot(cyc_nums, cap_pct, 'o-', 'Color', ch_colors(i,:), ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', ch);
    end
end
xlabel('Cycle', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Capacity Retention (%)', 'FontWeight', 'bold', 'FontSize', 12);
title('SOH: Capacity Retention', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9); grid on;
ylim([85 105]); hold off;

sgtitle('SOH Trend (Static Capacity, All Channels)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig3, fullfile(saveDir, 'DataLoad_SOH_Trend.fig'));

fprintf('--- Data Load Visualization Complete (3 figures saved) ---\n');
end
