% Generate_MasterRulers_PseudoOCV.m
% Description: Generates capacity-balanced voltage segment boundaries (Master Rulers)
% using the Pseudo-OCV (0.05C) cyc0 charge and discharge curves.
% The logic enforces a strict 3.0V - 4.2V range divided into 10 capacity-equal segments.

clear; clc; close all;

%% 1. Configuration
% Target Voltage Range defined by user 
target_V_min = 3.0;
target_V_max = 4.2;
num_segments = 12; % Changed from 10 to 12 based on dQ/dV valley strategy
num_boundaries = num_segments + 1;

% Paths
input_file = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\OCV_integrated.mat';
save_path = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver0317\MasterRulers_PseudoOCV.mat';


%"H:\내 드라이브\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat"
% Ensure directory exists
[save_dir, ~, ~] = fileparts(save_path);
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%% 2. Load Data
fprintf('Loading Pseudo-OCV Data...\n');
% .mat 내부 함수 핸들이 RPT_Postprocessing.m을 참조하므로 해당 경로를 path에 추가
addpath(fullfile(fileparts(input_file)));
data = load(input_file);
ocv0 = data.OCV_data.OCV_integrated_0;

% OCV curves are mapped to 201 SOC points (0 to 100)
soc_array = linspace(0, 100, 201)';
V_chg = ocv0.V_avg_SOC_charge(:);
V_dch = ocv0.V_avg_SOC_discharge(:);

%% 3. Master Ruler: 워크스루 3.3절 최종 확정 경계 (dQ/dV Valley 기반, 12 Segments)
% 8개 셀 앙상블 평균 OCV에서 findpeaks(-dQ/dV)로 검출된 Valley 기반 경계
fprintf('Applying confirmed 12-segment Master Ruler boundaries...\n');
V_bounds_chg = [3.000, 3.140, 3.270, 3.410, 3.535, 3.567, 3.600, 3.690, 3.774, 3.880, 3.990, 4.100, 4.200];
V_bounds_dch = [3.000, 3.130, 3.250, 3.380, 3.498, 3.531, 3.570, 3.660, 3.746, 3.860, 3.970, 4.090, 4.200];

%% 4. Construction and Saving
num_seg_chg = length(V_bounds_chg) - 1;
num_seg_dch = length(V_bounds_dch) - 1;
fprintf('Charge: %d segments, Discharge: %d segments\n', num_seg_chg, num_seg_dch);

MasterRuler_ver0317 = struct();
MasterRuler_ver0317.V_bounds_chg = V_bounds_chg;
MasterRuler_ver0317.V_bounds_dch = V_bounds_dch;
MasterRuler_ver0317.num_seg_chg = num_seg_chg;
MasterRuler_ver0317.num_seg_dch = num_seg_dch;
MasterRuler_ver0317.description = sprintf('dQ/dV Valley based: %g~%gV, Chg %d segs, Dch %d segs', target_V_min, target_V_max, num_seg_chg, num_seg_dch);

save(save_path, 'MasterRuler_ver0317');
fprintf('Successfully generated and saved to:\n %s\n', save_path);

%% 5. Visualization (Multi-Cycle Gradient & Segment Overlay)
fig = figure('Name', 'Pseudo-OCV Master Ruler (Multi-Cycle)', 'Position', [100, 100, 1000, 500]);

% Cycles to plot
cyc_list = [0, 200, 400, 600, 800, 1000];
colors = jet(length(cyc_list)); % Gradient from blue (cyc0) to red (cyc1000)

subplot(1,2,1); hold on;
% Draw alternating background colors and horizontal lines
seg_colors = [0.9 0.9 0.9; 0.95 0.95 1]; % Alternating light gray and light blue
for i = 1:num_seg_chg
    y1 = V_bounds_chg(i);
    y2 = V_bounds_chg(i+1);
    c_idx = mod(i-1, 2) + 1;
    patch([-10 100 100 -10], [y1 y1 y2 y2], seg_colors(c_idx,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    yline(y1, 'k--', 'HandleVisibility', 'off');
end
yline(V_bounds_chg(end), 'k--', 'HandleVisibility', 'off');

for i = 1:length(cyc_list)
    cyc = cyc_list(i);
    ocv_data = data.OCV_data.(sprintf('OCV_integrated_%d', cyc));
    cap_array_chg = soc_array / 100 * ocv_data.mean_capacity;
    plot(cap_array_chg, ocv_data.V_avg_SOC_charge(:), 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('cyc%d', cyc));
end
title('Charge Voltage Segment');
xlabel('Charge Capacity (Ah)'); ylabel('Voltage (V)');
grid on; ylim([target_V_min-0.05, target_V_max+0.05]); xlim([0, 70]);
legend('Location', 'southeast');

subplot(1,2,2); hold on;
% Draw alternating background colors and horizontal lines
for i = 1:num_seg_dch
    y1 = V_bounds_dch(i);
    y2 = V_bounds_dch(i+1);
    c_idx = mod(i-1, 2) + 1;
    patch([-10 100 100 -10], [y1 y1 y2 y2], seg_colors(c_idx,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    yline(y1, 'k--', 'HandleVisibility', 'off');
end
yline(V_bounds_dch(end), 'k--', 'HandleVisibility', 'off');

for i = 1:length(cyc_list)
    cyc = cyc_list(i);
    ocv_data = data.OCV_data.(sprintf('OCV_integrated_%d', cyc));
    cap_array_dch = soc_array / 100 * ocv_data.mean_capacity;
    plot(cap_array_dch, ocv_data.V_avg_SOC_discharge(:), 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('cyc%d', cyc));
end
title('Discharge Voltage Segment');
xlabel('Discharge Capacity (Ah)'); ylabel('Voltage (V)');
grid on; ylim([target_V_min-0.05, target_V_max+0.05]); xlim([0, 70]);
legend('Location', 'southwest');

% Save figure
savefig(fig, fullfile(save_dir, 'MasterRuler_Visualization.fig'));
fprintf('Figure saved as MasterRuler_Visualization.fig\n');

% Dynamic dQ/dV valley extraction (ver0317) - no hardcoded boundaries
