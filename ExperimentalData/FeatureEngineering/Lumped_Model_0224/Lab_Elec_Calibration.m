% Lab_Elec_Calibration.m
% Implements Phase 1-A of the lumped thermal model workflow.
% Extracts intrinsic cell electrical parameters (R0, R1, C1) 
% from laboratory DCIR (n1C discharge pulse) data.

clear; clc; close all;

fprintf('--- Starting Phase 1-A: Electrical Model Calibration (DCIR) ---\n');

%% 1. Configuration & Data Loading
% Load parsed RPT data for Ch09, 0cyc
parsed_file = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Parsed\RPT0_ch09_parsed.mat';
if ~exist(parsed_file, 'file')
    error('Parsed DCIR file not found: %s', parsed_file);
end

fprintf('Loading parsed DCIR data...\n');
load(parsed_file); % Loads 'pdata', 'data', 'intrim', etc.

% Convert pdata to the structured format used in RPT_DCIR_Final
data_struct = [];
for i = 1:length(pdata)
    data_struct(i).t = pdata(i).t;
    data_struct(i).V = pdata(i).V;
    data_struct(i).I = pdata(i).I;
    data_struct(i).type = char(pdata(i).type);
end

I_1C = 64;              % 1C Current [A]
tolC_main = 0.05;       % 5% tolerance
dt_thresh_n1C = 0.2;    % 0.1s sampling condition

% Identify discharge 1C (n1C) pulses
n1C_low = -(1+tolC_main) * I_1C;
n1C_high = -(1-tolC_main) * I_1C;

for k = 1:numel(data_struct)
    data_struct(k).avg_I = mean(data_struct(k).I);
    if numel(data_struct(k).t) >= 2
        data_struct(k).dt_mean = mean(diff(data_struct(k).t));
    else
        data_struct(k).dt_mean = Inf;
    end
end

is_n1C = arrayfun(@(x) ...
    (x.avg_I >= n1C_low && x.avg_I <= n1C_high) && ...
    (x.dt_mean <= dt_thresh_n1C) && numel(x.t) >= 2, data_struct);

n1C_idx = find(is_n1C);

if isempty(n1C_idx)
    error('No discharge 1C DCIR pulse found.');
end

% We will use the very first discharge 1C pulse for calibration
target_idx = n1C_idx(1);
pulse_data = data_struct(target_idx);

t_exp = pulse_data.t(:);
I_exp = pulse_data.I(:);
V_exp = pulse_data.V(:);

% Normalize time to start exactly at 0
if ~isempty(t_exp)
    t_exp = t_exp - t_exp(1);
end

% Assuming Voc is the voltage right before the pulse
% We can approximate it as the first point V_exp(1) 
% since the pulse just started, or we could take it from the previous rest step.
% RPT_DCIR_Final.m uses V_start = V_vec(1) for R calculations.
Voc_fit = V_exp(1) * ones(size(t_exp)); 

%% 2. Electrical Model Parameter Extraction (R0, R1, C1)
fprintf('\nStarting Electrical Optimization (R0, R1, C1) on DCIR Pulse...\n');
fprintf('Pulse index: %d, Duration: %.1fs, Avg Current: %.2fA\n', target_idx, max(t_exp), mean(I_exp));

% Variables: x(1) = R0, x(2) = R1, x(3) = C1
lb_elec = [1e-5, 1e-4, 10]; 
ub_elec = [1e-1, 1e-1, 50000];

% Objective Function: Sum of Squared Errors between V_sim and V_exp
obj_elec = @(x) calc_electrical_error(x, t_exp, I_exp, V_exp, Voc_fit);

% Optimize using Particle Swarm
options_ps = optimoptions('particleswarm', 'SwarmSize', 30, 'Display', 'iter', 'MaxIterations', 50);
[x_elec, fval_elec] = particleswarm(obj_elec, 3, lb_elec, ub_elec, options_ps);

R0_opt = x_elec(1);
R1_opt = x_elec(2);
C1_opt = x_elec(3);

fprintf('--> Electrical Calibration Results:\n');
fprintf('    R0 = %.5f Ohm\n', R0_opt);
fprintf('    R1 = %.5f Ohm\n', R1_opt);
fprintf('    C1 = %.1f F\n', C1_opt);
fprintf('    Final SSE = %.4f\n', fval_elec);

%% 3. Save Extracted Parameters
save_path = fullfile(pwd, 'Cell_Parameters_Elec_Lab.mat');
save(save_path, 'R0_opt', 'R1_opt', 'C1_opt', 'pulse_data');
fprintf('\nSuccessfully saved Electrical Parameters to: %s\n', save_path);

%% 4. Visualization (Saved strictly as .fig)
fig_handle = figure('Name', 'Lab Electrical Calibration', 'Position', [100, 100, 800, 500], 'Color', 'w');

[~, V_sim] = calc_electrical_error(x_elec, t_exp, I_exp, V_exp, Voc_fit);
plot(t_exp, V_exp, 'k', 'LineWidth', 2); hold on;
plot(t_exp, V_sim, 'r--', 'LineWidth', 2);
title('Electrical Validation (1C DCIR Discharge Pulse)');
xlabel('Time [s]'); ylabel('Voltage [V]');
legend('Experimental', 'Simulated', 'Location', 'best');
grid on;

fig_save_path = fullfile(pwd, 'Lab_Elec_Calibration_Results.fig');
savefig(fig_handle, fig_save_path);
close(fig_handle);
fprintf('Saved figure strictly as .fig to: %s\n', fig_save_path);

%% ========================================================================
% HELPER FUNCTIONS
function [sse, V_sim] = calc_electrical_error(x, t, I, V_exp, Voc)
    R0 = x(1);
    R1 = x(2);
    C1 = x(3);
    
    dt = [0; diff(t)];
    N = length(t);
    V_rc = zeros(N, 1);
    V_sim = zeros(N, 1);
    
    for k = 2:N
        % dV_rc/dt = -V_rc/(R1*C1) + I/C1
        dV_rc = (-V_rc(k-1)/(R1*C1) + I(k-1)/C1) * dt(k);
        V_rc(k) = V_rc(k-1) + dV_rc;
        
        % V_cell(t) = V_oc + I_cell(t)*R0 + V_rc(t)
        % Note: I is negative during discharge
        V_sim(k) = Voc(k) + I(k)*R0 + V_rc(k);
    end
    
    V_sim(1) = V_exp(1);
    
    % Sum of Squared Errors
    sse = sum((V_exp - V_sim).^2);
end
